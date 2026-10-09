/*
 * DeltaDLoop.hpp, part of HadronsMILC
 *
 * Sea-quark loops with one and two insertions of
 *     deltaD = D[U~] - D[U],
 * the difference of the HISQ operator on phased-and-resmeared links (U~:
 * thin links times exp(i q A), then HISQ smearing) and on the plain links.
 *
 * With the complete (exact low + diluted stochastic high) vectors of the
 * pure-QCD propagator S = D[U]^-1 = sum_i v_i w_i^dag, every right vector
 * satisfies D[U] v_j = w_j:
 *   low modes   w = (e +- o)/sqrt(2), o = Meooe(e)/(i lambda), v = w/(2m +- i lambda)
 *   high modes  w = P_H eta, v = D[U]^-1 w           (A2AVectors, projection route)
 * so only D[U~] is applied:
 *   M_ij = < w_i , D[U~] v_j - w_j >      (whole-lattice inner product)
 * and the loops are
 *   T_1 = Tr[S deltaD]          = sum_i  M_ii
 *   T_2 = Tr[S deltaD S deltaD] = sum_ij M_ij M_ji
 * written split by block (low / high), noise source and dilution component so
 * that the analysis can pair noise vectors correctly (different sources, or
 * the same source with different dilution components).
 *
 * U(1) fields: externalEmFields (fields passed in by name, e.g. the photon
 * shared with the valence observables) and/or internalEmFields (count fields
 * generated inside the module from an MGauge::StochEmFunc, seeded from
 * seedString).
 */
#ifndef HadronsMILC_MContraction_DeltaDLoop_hpp_
#define HadronsMILC_MContraction_DeltaDLoop_hpp_

#include <A2AMatrix.hpp>
#include <EigenPack.hpp>
#include <GridMilc/GridMilc.h>
#include <Grid/qcd/utils/HighlyImprovedStaggeredFermionImpl.h>
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MContraction)

class DeltaDInternalEmPar : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(DeltaDInternalEmPar, unsigned int, count,
                                  std::string, emFunc, std::string,
                                  seedString);
  DeltaDInternalEmPar() : count(0) {}
};

class DeltaDLoopPar : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(DeltaDLoopPar, std::string, gauge,
                                  std::string, boundary, double, mass,
                                  std::string, action, std::string, lowModes,
                                  std::string, left, std::string, right,
                                  unsigned int, nSources, std::string, charges,
                                  std::vector<std::string>, externalEmFields,
                                  DeltaDInternalEmPar, internalEmFields,
                                  unsigned int, block, unsigned int,
                                  checkShortcut, std::string, output);
};
// gauge:       thin, unphased links (as fed to the HISQSmear of 'action')
// boundary:    boundary phases for the smearing, as HISQSmear ("" = 1 1 1 -1)
// mass:        quark mass m (the operator carries 2m, MILC convention)
// action:      the unphased D[U] (only used for the low-mode partners and for
//              checkShortcut)
// lowModes:    mass-shifted eigenpack ("" = no exact low modes)
// left/right:  A2AVectors outputs <name>_w / <name>_v (projection route)
// nSources:    number of noise sources in left/right; each is diluted into
//              right.size()/nSources components, index source*nComp + comp
// charges:     whitespace-separated list of q e values, both signs listed
// block:       left vectors per run, and right vectors per run (half
//              deltaD v, half MeooeDag deltaD v); memory about 2 x block
//              full fields
// checkShortcut: number of low and of high vectors on which D[U] v - w is
//              also computed explicitly and logged (0 = off)

class DeltaDLoopResult : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(
      DeltaDLoopResult, std::string, field, double, charge, unsigned int,
      nLow, unsigned int, nSources, unsigned int, nComp, Complex, T1_low,
      std::vector<std::vector<Complex>>, T1_high, Complex, T2_LL,
      std::vector<Complex>, T2_LH, std::vector<std::vector<Complex>>, T2_HH,
      std::vector<std::vector<Complex>>, T2_HH_samecomp);
};
// T1_low            sum_{i low} M_ii
// T1_high[h][d]     M_aa, a = (source h, component d)
// T2_LL             sum_{i,j low} M_ij M_ji
// T2_LH[h]          sum_{i low, a in h} (M_ia M_ai + M_ai M_ia)
// T2_HH[h][h']      sum_{a in h, b in h'} M_ab M_ba
// T2_HH_samecomp[h][h'] the part of T2_HH with comp(a) == comp(b)
//                   (h == h': the a == b terms)

template <typename FImpl, typename Pack>
class TDeltaDLoopMILC : public Module<DeltaDLoopPar> {
public:
  FERM_TYPE_ALIASES(FImpl, );
  typedef typename FImpl::GaugeLinkField GaugeLinkField;
  typedef PhotonR::GaugeField EmField;
  typedef PhotonR::GaugeLinkField EmComp;

public:
  TDeltaDLoopMILC(const std::string name) : Module<DeltaDLoopPar>(name) {}
  virtual ~TDeltaDLoopMILC(void) {};
  virtual std::vector<std::string> getInput(void);
  virtual std::vector<std::string> getOutput(void);

protected:
  virtual void setup(void);
  virtual void execute(void);

protected: // StagDeltaDLoopFast derives from this module and replaces computeM
  // build the full-volume low-mode pair vectors k0..k0+n-1 (2n vectors,
  // slot 2k = (e+o)/sqrt2, 2k+1 = (e-o)/sqrt2, unit norm) into out[0..2n)
  void makeLowPairs(FermionField *out, const Pack &epack,
                    FMat &mat, unsigned int k0, unsigned int n);
  // phased + HISQ-smeared links for one field and charge
  void phasedLinks(GaugeField &Ufat, GaugeField &Ulong, const GaugeField &U,
                   const EmField &A, double charge);
  // M stored row-major, M[i * N + j], N = nLow_ + nHigh_
  virtual void computeM(std::vector<Complex> &M, FMat &Dtilde);
  DeltaDLoopResult traces(const std::vector<Complex> &M,
                          const std::string &field, double charge);

  unsigned int nLow_{0}, nHigh_{0}, nComp_{0};
};

MODULE_REGISTER(StagDeltaDLoop,
                ARG(TDeltaDLoopMILC<STAGIMPL, MassShiftEigenPack<STAGIMPL>>),
                MContraction);

/******************************************************************************
 *                    TDeltaDLoopMILC implementation                           *
 ******************************************************************************/
template <typename FImpl, typename Pack>
std::vector<std::string> TDeltaDLoopMILC<FImpl, Pack>::getInput(void) {
  std::vector<std::string> in = {par().gauge, par().action};
  if (!par().left.empty() || !par().right.empty()) {
    in.push_back(par().left);
    in.push_back(par().right);
  }
  if (!par().lowModes.empty()) {
    in.push_back(par().lowModes);
  }
  for (auto &f : par().externalEmFields) {
    in.push_back(f);
  }
  if (par().internalEmFields.count > 0) {
    in.push_back(par().internalEmFields.emFunc);
  }
  return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TDeltaDLoopMILC<FImpl, Pack>::getOutput(void) {
  return {};
}

template <typename FImpl, typename Pack>
void TDeltaDLoopMILC<FImpl, Pack>::setup(void) {
  if (par().block == 0) {
    HADRONS_ERROR(Argument, "block must be > 0");
  }
  if (par().nSources == 0) {
    HADRONS_ERROR(Argument, "nSources must be > 0");
  }
  if (par().externalEmFields.empty() && par().internalEmFields.count == 0) {
    HADRONS_ERROR(Argument, "no U(1) fields: set externalEmFields and/or "
                            "internalEmFields.count");
  }
  if (par().internalEmFields.count > 0 && par().internalEmFields.emFunc.empty()) {
    HADRONS_ERROR(Argument, "internalEmFields.count > 0 needs emFunc");
  }
  envTmpLat(EmField, "AField");
}

template <typename FImpl, typename Pack>
void TDeltaDLoopMILC<FImpl, Pack>::makeLowPairs(
    FermionField *out, const Pack &epack, FMat &mat,
    unsigned int k0, unsigned int n) {
  auto *rbGrid = envGetRbGrid(FermionField);
  FermionField half(rbGrid), partner(rbGrid);
  RealD invSqrt2 = 1. / std::sqrt(2.);

  for (unsigned int k = 0; k < n; ++k) {
    const FermionField &e = epack.evec[k0 + k];
    int cb = e.Checkerboard();
    int cbNeg = (cb == Even) ? Odd : Even;
    ComplexD eval_D = ComplexD(0., epack.eval[k0 + k].imag());
    RealD enorm = 1. / std::sqrt(norm2(e));

    // o = Meooe(e)/(i lambda): the opposite-parity partner, same norm as e
    partner = Zero();
    partner.Checkerboard() = cbNeg;
    mat.Meooe(e, partner);
    partner = (enorm * invSqrt2 / eval_D) * partner;
    half = (enorm * invSqrt2) * e;
    // expression templates do not carry the checkerboard label: set it
    half.Checkerboard() = cb;
    partner.Checkerboard() = cbNeg;

    out[2 * k] = Zero();
    setCheckerboard(out[2 * k], half);
    setCheckerboard(out[2 * k], partner);
    partner = -partner;
    partner.Checkerboard() = cbNeg;
    out[2 * k + 1] = Zero();
    setCheckerboard(out[2 * k + 1], half);
    setCheckerboard(out[2 * k + 1], partner);
  }
}

template <typename FImpl, typename Pack>
void TDeltaDLoopMILC<FImpl, Pack>::phasedLinks(GaugeField &Ufat,
                                               GaugeField &Ulong,
                                               const GaugeField &U,
                                               const EmField &A,
                                               double charge) {
  // U1Phase: thin links times exp(i q A_mu)
  startTimer("smear: U(1) phase");
  GaugeField Uq(U.Grid());
  EmComp Amu(A.Grid()), phase(U.Grid());
  GaugeLinkField Umu(U.Grid());
  for (unsigned int mu = 0; mu < Nd; ++mu) {
    Amu = peekLorentz(A, mu);
    Umu = peekLorentz(U, mu);
    phase = exp(timesI(charge * Amu));
    Umu = phase * Umu;
    pokeLorentz(Uq, Umu, mu);
  }
  stopTimer("smear: U(1) phase");

  // HISQSmear: rephase, fat7, U(3) projection, asqtad level 2
  StaggeredImplParams stagParams;
  if (!par().boundary.empty()) {
    stagParams.boundary_phases = strToVec<Complex>(par().boundary);
  } else {
    stagParams.boundary_phases = std::vector<Complex>{
        Complex(1.), Complex(1.), Complex(1.), Complex(-1.)};
  }
  auto *grid = envGetGrid(GaugeField);
  HighlyImprovedStaggeredFermionImpl<FImpl> hisq(grid, stagParams);
  GaugeField R(grid), V(grid), W(grid);
  startTimer("smear: rephase");
  hisq.rephase(R, Uq);
  stopTimer("smear: rephase");
  startTimer("smear: fat7");
  hisq.smear(V, R);
  stopTimer("smear: fat7");
  startTimer("smear: U(3) projection");
  hisq.project(W, V);
  stopTimer("smear: U(3) projection");
  startTimer("smear: asqtad level 2");
  hisq.smear(Ufat, Ulong, W);
  stopTimer("smear: asqtad level 2");
}

// M_ij = < w_i , D[U~] v_j - w_j >, i, j over [low pairs | high modes]
template <typename FImpl, typename Pack>
void TDeltaDLoopMILC<FImpl, Pack>::computeM(std::vector<Complex> &M,
                                            FMat &Dtilde) {
  auto *grid = envGetGrid(FermionField);
  auto &mat = envGet(FMat, par().action);
  std::vector<FermionField> none;
  auto &wHigh = nHigh_ ? envGet(std::vector<FermionField>, par().left) : none;
  auto &vHigh = nHigh_ ? envGet(std::vector<FermionField>, par().right) : none;
  const Pack *epack = nullptr;
  if (!par().lowModes.empty()) {
    epack = &envGet(Pack, par().lowModes);
  }
  const unsigned int N = nLow_ + nHigh_, B = par().block;
  // v_j per step of the outer loop: half of B, as each v_j gives two right
  // vectors (at least one eigenvector pair)
  const unsigned int BR = std::max(2u, B / 2);
  const int nt = env().getDim(Tp);

  M.assign((size_t)N * N, 0.);
  // right vectors for the current set of v_j, in ONE array so that a single
  // meson-field run covers both kinds:
  //   [0, nj)     deltaDv_j          = D[U~] v_j - w_j = deltaD v_j
  //   [nj, 2 nj)  MeooeDag_deltaDv_j = MeooeDag (deltaDv_j on the parity
  //               opposite to e), placed on the parity of e
  // when the v_j come from eigenvectors, [nj, 2 nj) first holds their pairs w_j
  std::vector<FermionField> rightVecs(2 * BR, grid);
  // left: up to B plain eigenvectors e_k (zero on the other parity)
  std::vector<FermionField> leftEvecs(B, grid);
  std::vector<StagGamma> gammas = {StagGamma()}; // (G1 G1): Gamma = 1
  std::vector<ComplexField> noMom;
  A2AWorkerSpinTasteStencil<FImpl> worker(grid, noMom, gammas, nullptr, Tp);
  std::vector<Complex> buf;

  // a set of BR v_j starting at index j0 (low pairs come in twos)
  auto blockEnd = [&](unsigned int i0) {
    unsigned int i1 = std::min(i0 + BR, N);
    if (i0 < nLow_) {
      i1 = std::min(i1, nLow_);
      if ((i1 - i0) % 2) {
        i1 -= 1; // keep pairs together
      }
      if (i1 == i0) {
        i1 = std::min(i0 + 2, nLow_);
      }
    }
    return i1;
  };
  // left low rows without building the pairs: with o = Meooe e/(i lambda),
  //   <(e +- o)/sqrt2, deltaDv> = |e|^-1/sqrt2 [ <e, deltaDv>
  //                                   +- (i/lambda) <e, MeooeDag_deltaDv> ]
  // (e normalised by |e|^-1 as in makeLowPairs; all evecs on one parity)
  auto *rbGrid = envGetRbGrid(FermionField);
  FermionField halfDeltaDv(rbGrid), halfMeooeDag(rbGrid);
  const unsigned int nEvec = nLow_ / 2;
  int cbE = Even;
  std::vector<RealD> lowScale;
  std::vector<ComplexD> lowCoef;
  if (nLow_) {
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

  for (unsigned int j0 = 0; j0 < N; j0 = blockEnd(j0)) {
    unsigned int j1 = blockEnd(j0), nj = j1 - j0;
    FermionField *deltaDv = rightVecs.data();
    FermionField *MeooeDag_deltaDv = rightVecs.data() + nj;
    // deltaDv_j = D[U~] v_j - w_j
    if (j0 < nLow_) {
      startTimer("M: right low pairs");
      makeLowPairs(MeooeDag_deltaDv, *epack, mat, j0 / 2, nj / 2); // w_j
      stopTimer("M: right low pairs");
      startTimer("M: right D[U~] v - w");
      for (unsigned int jj = 0; jj < nj; ++jj) {
        unsigned int k = (j0 + jj) / 2;
        ComplexD ev = epack->eval[k];          // 2m + i lambda
        ComplexD d = ((j0 + jj) % 2) ? ComplexD(ev.real(), -ev.imag()) : ev;
        Dtilde.M(MeooeDag_deltaDv[jj], deltaDv[jj]);
        deltaDv[jj] = (1. / d) * deltaDv[jj] - MeooeDag_deltaDv[jj];
      }
      stopTimer("M: right D[U~] v - w");
    } else {
      startTimer("M: right D[U~] v - w");
      for (unsigned int jj = 0; jj < nj; ++jj) {
        unsigned int a = j0 + jj - nLow_;
        Dtilde.M(vHigh[a], deltaDv[jj]);
        deltaDv[jj] = deltaDv[jj] - wHigh[a];
      }
      stopTimer("M: right D[U~] v - w");
    }
    // MeooeDag_deltaDv_j, for the left low rows (overwrites the pairs w_j)
    if (nLow_) {
      startTimer("M: right MeooeDag");
      for (unsigned int jj = 0; jj < nj; ++jj) {
        pickCheckerboard(cbNeg, halfDeltaDv, deltaDv[jj]);
        mat.MeooeDag(halfDeltaDv, halfMeooeDag);
        halfMeooeDag.Checkerboard() = cbE;
        MeooeDag_deltaDv[jj] = Zero();
        setCheckerboard(MeooeDag_deltaDv[jj], halfMeooeDag);
      }
      stopTimer("M: right MeooeDag");
    }
    // left low rows: B eigenvectors e_k at a time (rows 2k and 2k+1), one
    // run against all 2 nj right vectors
    for (unsigned int k0 = 0; k0 < nEvec; k0 += B) {
      unsigned int nk = std::min(B, nEvec - k0);
      startTimer("M: left evec copy");
      for (unsigned int kk = 0; kk < nk; ++kk) {
        leftEvecs[kk] = Zero();
        setCheckerboard(leftEvecs[kk], epack->evec[k0 + kk]);
      }
      stopTimer("M: left evec copy");
      startTimer("M: overlaps low");
      buf.assign((size_t)nt * nk * 2 * nj, 0.);
      A2AMatrixSet<Complex> m(buf.data(), 1, 1, nt, nk, 2 * nj);
      worker.resetCache(); // buffers are reused with new contents
      worker.StagMesonField(m, leftEvecs.data(), rightVecs.data(), (int)nk,
                            (int)(2 * nj));
      for (unsigned int kk = 0; kk < nk; ++kk)
        for (unsigned int jj = 0; jj < nj; ++jj) {
          Complex eDeltaDv = 0., eMeooeDag = 0.;
          for (int t = 0; t < nt; ++t) {
            eDeltaDv += m(0, 0, t, kk, jj);
            eMeooeDag += m(0, 0, t, kk, nj + jj);
          }
          Complex odd = lowCoef[k0 + kk] * eMeooeDag;
          RealD sc = lowScale[k0 + kk];
          M[(size_t)(2 * (k0 + kk)) * N + j0 + jj] = sc * (eDeltaDv + odd);
          M[(size_t)(2 * (k0 + kk) + 1) * N + j0 + jj] = sc * (eDeltaDv - odd);
        }
      stopTimer("M: overlaps low");
    }
    // left high rows: B noise vectors w_i at a time, against deltaDv only
    for (unsigned int i0 = nLow_; i0 < N; i0 += B) {
      unsigned int ni = std::min(B, N - i0);
      startTimer("M: overlaps high");
      buf.assign((size_t)nt * ni * nj, 0.);
      A2AMatrixSet<Complex> m(buf.data(), 1, 1, nt, ni, nj);
      worker.resetCache(); // buffers are reused with new contents
      worker.StagMesonField(m, &wHigh[i0 - nLow_], deltaDv, (int)ni, (int)nj);
      for (unsigned int ii = 0; ii < ni; ++ii)
        for (unsigned int jj = 0; jj < nj; ++jj) {
          Complex sum = 0.;
          for (int t = 0; t < nt; ++t) {
            sum += m(0, 0, t, ii, jj);
          }
          M[(size_t)(i0 + ii) * N + j0 + jj] = sum;
        }
      stopTimer("M: overlaps high");
    }
  }
}

template <typename FImpl, typename Pack>
DeltaDLoopResult TDeltaDLoopMILC<FImpl, Pack>::traces(
    const std::vector<Complex> &Mv, const std::string &field, double charge) {
  const unsigned int nS = par().nSources, N = nLow_ + nHigh_;
  auto M = [&](unsigned int i, unsigned int j) { return Mv[(size_t)i * N + j]; };
  DeltaDLoopResult r;
  r.field = field;
  r.charge = charge;
  r.nLow = nLow_;
  r.nSources = nS;
  r.nComp = nComp_;
  r.T1_low = 0.;
  r.T2_LL = 0.;
  r.T1_high.assign(nS, std::vector<Complex>(nComp_, 0.));
  r.T2_LH.assign(nS, 0.);
  r.T2_HH.assign(nS, std::vector<Complex>(nS, 0.));
  r.T2_HH_samecomp.assign(nS, std::vector<Complex>(nS, 0.));

  for (unsigned int i = 0; i < nLow_; ++i) {
    r.T1_low += M(i, i);
    for (unsigned int j = 0; j < nLow_; ++j) {
      r.T2_LL += M(i, j) * M(j, i);
    }
  }
  for (unsigned int a = 0; a < nHigh_; ++a) {
    unsigned int h = a / nComp_, d = a % nComp_, ia = nLow_ + a; // nHigh_ > 0 here
    r.T1_high[h][d] = M(ia, ia);
    for (unsigned int i = 0; i < nLow_; ++i) {
      r.T2_LH[h] += M(i, ia) * M(ia, i) + M(ia, i) * M(i, ia);
    }
    for (unsigned int b = 0; b < nHigh_; ++b) {
      unsigned int hb = b / nComp_, db = b % nComp_, ib = nLow_ + b;
      Complex x = M(ia, ib) * M(ib, ia);
      r.T2_HH[h][hb] += x;
      if (d == db) {
        r.T2_HH_samecomp[h][hb] += x;
      }
    }
  }
  return r;
}

template <typename FImpl, typename Pack>
void TDeltaDLoopMILC<FImpl, Pack>::execute(void) {
  auto &U = envGet(GaugeField, par().gauge);
  auto &mat = envGet(FMat, par().action);
  bool hasHigh = !par().left.empty();
  std::vector<FermionField> none;
  auto &wHigh = hasHigh ? envGet(std::vector<FermionField>, par().left) : none;
  auto &vHigh = hasHigh ? envGet(std::vector<FermionField>, par().right) : none;
  auto *grid = envGetGrid(FermionField);
  auto *rbGrid = envGetRbGrid(FermionField);
  auto *ggrid = envGetGrid(GaugeField);

  nLow_ = 0;
  if (!par().lowModes.empty()) {
    auto &epack = envGet(Pack, par().lowModes);
    nLow_ = 2 * epack.evec.size();
    if (std::abs(epack.eval[0].real() - 2. * par().mass) > 1e-12) {
      HADRONS_ERROR(Argument, "lowModes mass " +
                                  std::to_string(epack.eval[0].real() / 2.) +
                                  " differs from mass");
    }
  }
  nHigh_ = vHigh.size();
  if (wHigh.size() != nHigh_ || nHigh_ % par().nSources) {
    HADRONS_ERROR(Size, "left/right sizes do not match nSources");
  }
  nComp_ = nHigh_ / par().nSources; // 0 without high modes
  LOG(Message) << "deltaD loops: " << nLow_ << " low-mode vectors, " << nHigh_
               << " high-mode vectors (" << par().nSources << " sources x "
               << nComp_ << " components), block " << par().block << std::endl;

  // optional check of the shortcut D[U] v_j = w_j
  if (par().checkShortcut > 0) {
    FermionField Dv(grid), v(grid);
    std::vector<FermionField> pairs(2, grid);
    unsigned int nl = std::min(par().checkShortcut, nLow_ / 2);
    for (unsigned int k = 0; k < nl; ++k) {
      auto &ep = envGet(Pack, par().lowModes);
      makeLowPairs(pairs.data(), ep, mat, k, 1);
      ComplexD ev = ep.eval[k];
      for (int s = 0; s < 2; ++s) {
        ComplexD d = s ? ComplexD(ev.real(), -ev.imag()) : ev;
        v = (1. / d) * pairs[s];
        mat.M(v, Dv);
        Dv = Dv - pairs[s];
        LOG(Message) << "shortcut check low " << 2 * k + s << ": |D[U]v - w|/|w| = "
                     << std::sqrt(norm2(Dv) / norm2(pairs[s])) << std::endl;
      }
    }
    unsigned int nh = std::min(par().checkShortcut, nHigh_);
    for (unsigned int a = 0; a < nh; ++a) {
      mat.M(vHigh[a], Dv);
      Dv = Dv - wHigh[a];
      LOG(Message) << "shortcut check high " << a << ": |D[U]v - w|/|w| = "
                   << std::sqrt(norm2(Dv) / norm2(wHigh[a])) << std::endl;
    }
  }

  std::vector<double> charges = strToVec<double>(par().charges);
  std::vector<DeltaDLoopResult> results;
  GaugeField Ufat(ggrid), Ulong(ggrid);
  std::vector<Complex> M;

  auto doField = [&](const EmField &A, const std::string &label) {
    for (double q : charges) {
      LOG(Message) << "field '" << label << "', charge " << q << std::endl;
      phasedLinks(Ufat, Ulong, U, A, q);
      startTimer("D[U~] build");
      ImprovedStaggeredFermion<FImpl> Dtilde(*grid, *rbGrid, 2. * par().mass);
      Dtilde.ImportGaugeSimple(Ulong, Ufat);
      stopTimer("D[U~] build");
      startTimer("M: total");
      computeM(M, Dtilde);
      stopTimer("M: total");
      startTimer("traces");
      results.push_back(traces(M, label, q));
      stopTimer("traces");
      RealD f2 = 0.;
      for (auto &x : M) {
        f2 += real(conjugate(x) * x);
      }
      LOG(Message) << "  T1_low = " << results.back().T1_low
                   << " ; |M|_F = " << std::sqrt(f2) << std::endl;
    }
  };

  for (auto &name : par().externalEmFields) {
    doField(envGet(EmField, name), "external/" + name);
  }
  if (par().internalEmFields.count > 0) {
    envGetTmp(EmField, AField);
    auto &photon = envGet(PhotonR, par().internalEmFields.emFunc);
    auto &weight =
        envGet(EmComp, "_" + par().internalEmFields.emFunc + "_weight");
    auto &rng = rng4d();
    if (!par().internalEmFields.seedString.empty()) {
      rng.SeedUniqueString(par().internalEmFields.seedString);
    }
    for (unsigned int f = 0; f < par().internalEmFields.count; ++f) {
      photon.StochasticField(AField, rng, weight);
      doField(AField, "internal/" + std::to_string(f));
    }
  }

  saveResult(par().output, "deltaD", results);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MContraction_DeltaDLoop_hpp_
