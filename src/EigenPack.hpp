/*
 * EigenPack.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: Ryan Hill <rchrys.hill@gmail.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution
 * directory.
 */

/*  END LEGAL */
#ifndef HadronsMILC_EigenPack_hpp_
#define HadronsMILC_EigenPack_hpp_

#include <Hadrons/Global.hpp>
#include <Grid/parallelIO/QudaPartfileIO.h>

BEGIN_HADRONS_NAMESPACE

// ---------------------------------------------------------------------------
// AdaptorEigenPackMILC
//
// Stores Nvec eigenvectors of the staggered Dslash^dag Dslash operator
// (or the full Dslash for inflated files), together with the MILC-shifted
// eigenvalues used by QUDA:
//
//   eval[i] = ComplexD(2*mass, sqrt(lambda_i))
//
// where lambda_i is the i-th eigenvalue of M^dag M.
// ---------------------------------------------------------------------------
template <typename Field> class AdaptorEigenPackMILC {
public:
  AdaptorEigenPackMILC(std::vector<Field> &_evec,
                       const std::vector<RealD> &_eval, RealD _mass = 0.0)
      : evec(_evec), mass(_mass) {
    // Store shifted M eigenvalues instead of massless M^dagM eigenvalues
    eval.resize(_eval.size(), 0.0);

    Real m = 2 * _mass;

    for (int i = 0; i < (int)eval.size(); i++) {
      eval[i] = ComplexD(m, sqrt(_eval[i]));
    }
  }

public:
  std::vector<Field> &evec;
  std::vector<ComplexD> eval;
  RealD mass;
};

template <typename FImpl>
using MassShiftEigenPack = AdaptorEigenPackMILC<typename FImpl::FermionField>;

// ---------------------------------------------------------------------------
// QUDA partfile I/O for staggered eigenpacks
//
// writeQudaStaggeredEigenpack:
//   Write a QUDA/MILC eigenvector partfile set from Grid staggered fields.
//   evecs: full-lattice colour-vector fields (one per eigenvector).
//   evals: MILC-shifted eigenvalues as ComplexD (from AdaptorEigenPackMILC).
//   stem:  base filename; each rank writes <stem>.vol<NNNN>.
//   singlePrec: write in single (float) precision (default: true).
//
// readQudaStaggeredEigenpack:
//   Read a QUDA/MILC eigenvector partfile set into Grid staggered fields.
//   evecs must be pre-allocated to the expected number of eigenvectors.
//   evals is resized and filled from file (or zeroed if no eval record).
// ---------------------------------------------------------------------------

template <typename FImpl>
void writeQudaStaggeredEigenpack(
    const std::string                                     &stem,
    std::vector<typename FImpl::FermionField>             &evecs,
    const std::vector<ComplexD>                           &evals,
    GridBase                                              *grid,
    bool singlePrec = true)
{
  Grid::writeQudaEigenpackPartfile(stem, evecs, evals, grid,
                                   /*nColor=*/3, /*nSpin=*/1, singlePrec);
}

template <typename FImpl>
void readQudaStaggeredEigenpack(
    const std::string                           &stem,
    std::vector<typename FImpl::FermionField>   &evecs,
    std::vector<ComplexD>                       &evals,
    GridBase                                    *grid)
{
  Grid::readQudaEigenpackPartfile(stem, evecs, evals, grid,
                                  /*nColor=*/3, /*nSpin=*/1);
}

END_HADRONS_NAMESPACE

#endif // HadronsMILC_EigenPack_hpp_
