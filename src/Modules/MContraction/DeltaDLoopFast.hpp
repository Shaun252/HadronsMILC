/*
 * DeltaDLoopFast.hpp, part of HadronsMILC
 *
 * StagDeltaDLoop with M computed by dense matrix products on the device
 * (GridBLAS) instead of meson-field runs. Parameters, smearing, traces and
 * output are inherited unchanged from StagDeltaDLoop, which stays as the
 * reference implementation.
 *
 * M_ij = < w_i , D[U~] v_j - w_j > = < w_i , deltaD v_j >. Each vector is a
 * contiguous array of complex numbers on the device (3V numbers, 3V/2 for one
 * parity), and a whole-lattice inner product is sum conj(a) b over that array
 * whatever the site layout, so a block of M is one product A^dag B of
 * column-stacked vectors:
 *   rows from eigenvectors: <(e +- o)/sqrt2, deltaDv>
 *       = |e|^-1/sqrt2 [ <e, deltaDv> +- (i/lambda) <e, MeooeDag deltaDv> ]
 *     with e and both right vectors on e's parity only (half the length);
 *   rows from noise: <w, deltaDv> on the full lattice.
 */
#ifndef HadronsMILC_MContraction_DeltaDLoopFast_hpp_
#define HadronsMILC_MContraction_DeltaDLoopFast_hpp_

#include <Modules/MContraction/DeltaDLoop.hpp>

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MContraction)

// out[i] = sum_c part[c * n + i], c < nChunk (a free function: nvcc does not
// allow device lambdas inside protected member functions or other lambdas)
inline void sumChunks(ComplexD *out, const ComplexD *part, size_t n,
                      unsigned int nChunk) {
  accelerator_for(idx, n, 1, {
    ComplexD s = part[idx];
    for (unsigned int c = 1; c < nChunk; ++c) {
      s += part[c * n + idx];
    }
    out[idx] = s;
  });
}

template <typename FImpl, typename Pack>
class TDeltaDLoopFastMILC : public TDeltaDLoopMILC<FImpl, Pack> {
public:
  FERM_TYPE_ALIASES(FImpl, );
  typedef TDeltaDLoopMILC<FImpl, Pack> Base;

public:
  TDeltaDLoopFastMILC(const std::string name) : Base(name) {}
  virtual ~TDeltaDLoopFastMILC(void) {};

protected:
  virtual void computeM(std::vector<Complex> &M, FMat &Dtilde);
};

MODULE_REGISTER(StagDeltaDLoopFast,
                ARG(TDeltaDLoopFastMILC<STAGIMPL, MassShiftEigenPack<STAGIMPL>>),
                MContraction);

/******************************************************************************
 *                    TDeltaDLoopFastMILC implementation                       *
 ******************************************************************************/
// M stored row-major, M[i * N + j], N = nLow_ + nHigh_
template <typename FImpl, typename Pack>
void TDeltaDLoopFastMILC<FImpl, Pack>::computeM(std::vector<Complex> &M,
                                                FMat &Dtilde) {
  auto *grid = envGetGrid(FermionField);
  auto *rbGrid = envGetRbGrid(FermionField);
  auto &mat = envGet(FMat, this->par().action);
  std::vector<FermionField> none;
  const unsigned int nLow = this->nLow_, nHigh = this->nHigh_;
  auto &wHigh =
      nHigh ? envGet(std::vector<FermionField>, this->par().left) : none;
  auto &vHigh =
      nHigh ? envGet(std::vector<FermionField>, this->par().right) : none;
  const Pack *epack = nullptr;
  if (!this->par().lowModes.empty()) {
    epack = &envGet(Pack, this->par().lowModes);
  }
  const unsigned int N = nLow + nHigh, B = this->par().block;
  // v_j per step of the outer loop: half of B, as each v_j gives two right
  // vectors (at least one eigenvector pair)
  const unsigned int BR = std::max(2u, B / 2);
  const unsigned int nEvec = nLow / 2;
  M.assign((size_t)N * N, 0.);

  // complex numbers per vector, full lattice and one parity
  typedef typename FermionField::vector_object SiteVec;
  const size_t Kfull = grid->oSites() * sizeof(SiteVec) / sizeof(ComplexD);
  const size_t Khalf = rbGrid->oSites() * sizeof(SiteVec) / sizeof(ComplexD);

  // column-stacked device buffers (about 2 B full vectors in all):
  //   left:      B noise w_i, or 2 B eigenvectors e_k (half length)
  //   rightFull: deltaDv_j                                  (BR columns)
  //   rightHalf: [deltaDv_j on e's parity | MeooeDag_deltaDv_j] (2 BR columns)
  deviceVector<ComplexD> left((size_t)B * Kfull);
  deviceVector<ComplexD> rightFull((size_t)BR * Kfull);
  deviceVector<ComplexD> rightHalf(nLow ? (size_t)2 * BR * Khalf : 1);
  const unsigned int maxRows = std::max(B, 2 * B), maxCols = 2 * BR;
  deviceVector<ComplexD> prod((size_t)maxRows * maxCols);
  std::vector<ComplexD> prodHost((size_t)maxRows * maxCols);
  deviceVector<ComplexD *> pA(1), pB(1), pC(1);
  GridBLAS blas;

  auto copyIn = [&](const auto &f, ComplexD *dst, size_t K) {
    autoView(f_v, f, AcceleratorRead);
    acceleratorCopyDeviceToDevice((void *)&f_v[0], (void *)dst,
                                  K * sizeof(ComplexD));
  };
  // the summed dimension K (millions) is cut into nChunk pieces done as one
  // strided batch, then the partial m x n blocks are added: a single product
  // has so few m x n output tiles that most of the device sits idle
  const unsigned int maxChunk = 128;
  deviceVector<ComplexD> partial((size_t)maxChunk * maxRows * maxCols);
  // prodHost[r + m * c] = sum_x conj(A[x, r]) Bm[x, c], summed over ranks
  auto product = [&](ComplexD *A, ComplexD *Bm, unsigned int m,
                     unsigned int n, size_t K) {
#ifdef GRID_CUDA
    unsigned int nChunk = maxChunk;
    while (K % nChunk) {
      nChunk /= 2;
    }
    const size_t Kc = K / nChunk, mn = (size_t)m * n;
    cuDoubleComplex one = make_cuDoubleComplex(1., 0.);
    cuDoubleComplex zero = make_cuDoubleComplex(0., 0.);
    auto err = cublasZgemmStridedBatched(
        GridBLAS::gridblasHandle, CUBLAS_OP_C, CUBLAS_OP_N, (int)m, (int)n,
        (int)Kc, &one, (cuDoubleComplex *)A, (int)K, (long long)Kc,
        (cuDoubleComplex *)Bm, (int)K, (long long)Kc, &zero,
        (cuDoubleComplex *)&partial[0], (int)m, (long long)mn, (int)nChunk);
    GRID_ASSERT(err == CUBLAS_STATUS_SUCCESS);
    blas.synchronise();
    sumChunks(&prod[0], &partial[0], mn, nChunk);
#else
    acceleratorPut(pA[0], A);
    acceleratorPut(pB[0], Bm);
    acceleratorPut(pC[0], &prod[0]);
    blas.gemmBatched(GridBLAS_OP_C, GridBLAS_OP_N, (int)m, (int)n, (int)K,
                     ComplexD(1.), pA, pB, ComplexD(0.), pC);
    blas.synchronise();
#endif
    acceleratorCopyFromDevice(&prod[0], prodHost.data(),
                              (size_t)m * n * sizeof(ComplexD));
    grid->GlobalSumVector(prodHost.data(), (int)(m * n));
  };

  // a set of BR v_j starting at index j0 (low pairs come in twos)
  auto blockEnd = [&](unsigned int i0) {
    unsigned int i1 = std::min(i0 + BR, N);
    if (i0 < nLow) {
      i1 = std::min(i1, nLow);
      if ((i1 - i0) % 2) {
        i1 -= 1; // keep pairs together
      }
      if (i1 == i0) {
        i1 = std::min(i0 + 2, nLow);
      }
    }
    return i1;
  };

  int cbE = Even;
  std::vector<RealD> lowScale;
  std::vector<ComplexD> lowCoef;
  if (nLow) {
    cbE = epack->evec[0].Checkerboard();
    for (unsigned int k = 0; k < nEvec; ++k) {
      if (epack->evec[k].Checkerboard() != cbE) {
        HADRONS_ERROR(Argument, "lowModes eigenvectors not all on one parity");
      }
      lowScale.push_back(1. / std::sqrt(2. * norm2(epack->evec[k])));
      lowCoef.push_back(ComplexD(0., 1. / epack->eval[k].imag()));
    }
  }
  int cbNeg = (cbE == Even) ? Odd : Even;
  std::vector<FermionField> pair(2, grid);
  FermionField deltaDv(grid);
  FermionField halfDeltaDv(rbGrid), halfMeooeDag(rbGrid);

  for (unsigned int j0 = 0; j0 < N; j0 = blockEnd(j0)) {
    unsigned int j1 = blockEnd(j0), nj = j1 - j0;
    // right vectors: deltaDv_j = D[U~] v_j - w_j, and its two halves
    for (unsigned int jj = 0; jj < nj; ++jj) {
      unsigned int j = j0 + jj;
      this->startTimer("M: right D[U~] v - w");
      if (j < nLow) {
        unsigned int k = j / 2;
        if (j % 2 == 0 || jj == 0) {
          this->stopTimer("M: right D[U~] v - w");
          this->startTimer("M: right low pairs");
          this->makeLowPairs(pair.data(), *epack, mat, k, 1);
          this->stopTimer("M: right low pairs");
          this->startTimer("M: right D[U~] v - w");
        }
        ComplexD ev = epack->eval[k]; // 2m + i lambda
        ComplexD d = (j % 2) ? ComplexD(ev.real(), -ev.imag()) : ev;
        Dtilde.M(pair[j % 2], deltaDv);
        deltaDv = (1. / d) * deltaDv - pair[j % 2];
      } else {
        unsigned int a = j - nLow;
        Dtilde.M(vHigh[a], deltaDv);
        deltaDv = deltaDv - wHigh[a];
      }
      this->stopTimer("M: right D[U~] v - w");
      this->startTimer("M: copy to buffers");
      copyIn(deltaDv, &rightFull[(size_t)jj * Kfull], Kfull);
      this->stopTimer("M: copy to buffers");
      if (nLow) {
        this->startTimer("M: right MeooeDag");
        pickCheckerboard(cbE, halfDeltaDv, deltaDv);
        this->stopTimer("M: right MeooeDag");
        this->startTimer("M: copy to buffers");
        copyIn(halfDeltaDv, &rightHalf[(size_t)jj * Khalf], Khalf);
        this->stopTimer("M: copy to buffers");
        this->startTimer("M: right MeooeDag");
        pickCheckerboard(cbNeg, halfDeltaDv, deltaDv);
        mat.MeooeDag(halfDeltaDv, halfMeooeDag);
        halfMeooeDag.Checkerboard() = cbE;
        this->stopTimer("M: right MeooeDag");
        this->startTimer("M: copy to buffers");
        copyIn(halfMeooeDag, &rightHalf[(size_t)(nj + jj) * Khalf], Khalf);
        this->stopTimer("M: copy to buffers");
      }
    }
    // rows from eigenvectors: 2 B e_k at a time (rows 2k and 2k+1)
    for (unsigned int k0 = 0; k0 < nEvec; k0 += 2 * B) {
      unsigned int nk = std::min(2 * B, nEvec - k0);
      this->startTimer("M: copy to buffers");
      for (unsigned int kk = 0; kk < nk; ++kk) {
        copyIn(epack->evec[k0 + kk], &left[(size_t)kk * Khalf], Khalf);
      }
      this->stopTimer("M: copy to buffers");
      this->startTimer("M: matrix products");
      product(&left[0], &rightHalf[0], nk, 2 * nj, Khalf);
      this->stopTimer("M: matrix products");
      for (unsigned int kk = 0; kk < nk; ++kk)
        for (unsigned int jj = 0; jj < nj; ++jj) {
          ComplexD eDeltaDv = prodHost[kk + (size_t)nk * jj];
          ComplexD eMeooeDag = prodHost[kk + (size_t)nk * (nj + jj)];
          ComplexD odd = lowCoef[k0 + kk] * eMeooeDag;
          RealD sc = lowScale[k0 + kk];
          M[(size_t)(2 * (k0 + kk)) * N + j0 + jj] = sc * (eDeltaDv + odd);
          M[(size_t)(2 * (k0 + kk) + 1) * N + j0 + jj] = sc * (eDeltaDv - odd);
        }
    }
    // rows from noise: B w_i at a time, against deltaDv on the full lattice
    for (unsigned int i0 = nLow; i0 < N; i0 += B) {
      unsigned int ni = std::min(B, N - i0);
      this->startTimer("M: copy to buffers");
      for (unsigned int ii = 0; ii < ni; ++ii) {
        copyIn(wHigh[i0 + ii - nLow], &left[(size_t)ii * Kfull], Kfull);
      }
      this->stopTimer("M: copy to buffers");
      this->startTimer("M: matrix products");
      product(&left[0], &rightFull[0], ni, nj, Kfull);
      this->stopTimer("M: matrix products");
      for (unsigned int ii = 0; ii < ni; ++ii)
        for (unsigned int jj = 0; jj < nj; ++jj) {
          M[(size_t)(i0 + ii) * N + j0 + jj] = prodHost[ii + (size_t)ni * jj];
        }
    }
  }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MContraction_DeltaDLoopFast_hpp_
