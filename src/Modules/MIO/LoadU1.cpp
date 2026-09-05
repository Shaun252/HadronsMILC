/*
 * LoadU1.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Michael Lynch <michaellynch628@gmail.com>
 * Author: Shaun Lahert <shaun.lahert@gmail.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
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
#include "LoadU1.hpp"

using namespace Grid;
using namespace Hadrons;
using namespace MIO;

/******************************************************************************
*                       TLoadU1 implementation                                *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
TLoadU1::TLoadU1(const std::string name)
: Module<LoadU1Par>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TLoadU1::getInput(void)
{
    std::vector<std::string> in;

    return in;
}

std::vector<std::string> TLoadU1::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void TLoadU1::setup(void)
{
    envCreateLat(EmField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
void TLoadU1::execute(void)
{
    LOG(Message) << "Loading MILC U(1) photon field from file '" << par().file
                 << "'" << std::endl;

    auto &A = envGet(EmField, getName());

    // Build the single-precision read grid directly rather than through
    // envGetGrid: Hadrons pairs every grid it caches with a red-black one, and
    // GridRedBlackCartesian asserts that all reduced dimensions are even. On a
    // small lattice the vRealF SIMD layout leaves a reduced dimension of 1 and
    // that assertion aborts (it bites at 4^4, not at 64^3x96). No checkerboard
    // is needed to read a field.
    auto        *g4 = envGetGrid(EmField);
    GridCartesian rgrid(g4->GlobalDimensions(),
                        GridDefaultSimd(Nd, vRealF::Nsimd()),
                        g4->ProcessorGrid());
    U1AField      Araw(&rgrid);

    // readLimeLatticeBinaryObject verifies the SciDAC checksum itself and
    // asserts that the record size matches the lattice, so a wrong file or a
    // wrong volume is caught here rather than downstream.
    GridLimeReader reader;

    reader.open(par().file);
    reader.readLimeLatticeBinaryObject(Araw, std::string(SCIDAC_BINARY_DATA));
    reader.close();

    // A_mu is real: widen to the default-precision complex EM field, leaving
    // the imaginary part zero.
    //
    // Grid's toComplex() CANNOT be used for this. It is the strict inverse of
    // toReal(), which duplicates a complex field's real part into both halves
    // of every complex SIMD slot; toComplex() therefore asserts rr==ri and
    // aborts on a real field holding independent values in adjacent lanes
    // (Grid_vector_types.h:1047). Nor can the two fields be walked with a
    // shared site index: vRealF and vComplex have different SIMD widths, so
    // envGetGrid gives them different SIMD layouts.
    //
    // unvectorizeToLexOrdArray/vectorizeFromLexOrdArray are the layout-safe
    // route: both grids share the MPI decomposition and hence lSites() and the
    // local lexicographic ordering, so the scalar arrays correspond site by
    // site.
    typedef typename U1AField::scalar_object U1ASite;
    typedef typename EmField::scalar_object  EmSite;

    std::vector<U1ASite> abuf;
    std::vector<EmSite>  ebuf;

    unvectorizeToLexOrdArray(abuf, Araw);
    ebuf.resize(abuf.size());
    thread_for(i, abuf.size(),
    {
        for (unsigned int mu = 0; mu < Nd; ++mu)
        {
            ebuf[i](mu)()() = Complex(abuf[i](mu)()(), 0.);
        }
    });
    vectorizeFromLexOrdArray(ebuf, A);

    LOG(Message) << "Photon field norm2 = " << norm2(A) << std::endl;
}
