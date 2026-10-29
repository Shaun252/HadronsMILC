/*
 * PureGaugeEm.cpp, part of HadronsMILC
 *
 * Copyright (C) 2015 - 2026
 *
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
#include "PureGaugeEm.hpp"

using namespace Grid;
using namespace Hadrons;
using namespace MGauge;

TPureGaugeEm::TPureGaugeEm(const std::string name)
    : Module<PureGaugeEmPar>(name) {}

std::vector<std::string> TPureGaugeEm::getInput(void) { return {}; }

std::vector<std::string> TPureGaugeEm::getOutput(void) { return {getName()}; }

void TPureGaugeEm::setup(void) { envCreateLat(EmField, getName()); }

void TPureGaugeEm::execute(void) {
    LOG(Message) << "Pure-gauge U(1) field, A_mu = " << par().scale
                 << " [lambda(x+mu) - lambda(x)]" << std::endl;

    auto &A   = envGet(EmField, getName());
    auto &rng = rng4d();
    if (!par().seedString.empty()) {
        rng.SeedUniqueString(par().seedString);
    }

    // real Gaussian lambda: the complex Gaussian's real part
    EmComp lambda(A.Grid()), Amu(A.Grid());
    gaussian(rng, lambda);
    lambda = 0.5*(lambda + conjugate(lambda));

    for (unsigned int mu = 0; mu < Nd; ++mu) {
        Amu = par().scale*(Cshift(lambda, mu, 1) - lambda);
        pokeLorentz(A, Amu, mu);
    }
    LOG(Message) << "norm2(lambda) = " << norm2(lambda)
                 << ", norm2(A) = " << norm2(A) << std::endl;
}
