/*
 * PureGaugeEm.hpp, part of HadronsMILC
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
#ifndef HadronsMILC_MGauge_PureGaugeEm_hpp_
#define HadronsMILC_MGauge_PureGaugeEm_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 Pure-gauge U(1) field, for the Ward-identity test of the sea-quark loops

 scale       s in A_mu(x) = s [lambda(x + mu) - lambda(x)]
 seedString  seeds the 4d RNG for lambda ("" = leave the RNG as it is)

 lambda(x) is a real unit Gaussian field on the (periodic) lattice. The output
 has zero field strength, so phasing the thin links with charge q is the U(1)
 gauge transformation
     U_mu(x) -> G(x) U_mu(x) G(x + mu)^*,   G = exp(i q s lambda),
 and, since HISQ smearing and the U(3) projection are gauge covariant,
 D[U~] = G D[U] G^dag: the sea-loop determinant ratio is exactly 1. The
 output is a PhotonR::GaugeField, so it is passed to the sea loops (or
 U1Phase) like any other photon. Its zero mode vanishes (periodic lambda).
 ******************************************************************************/

BEGIN_MODULE_NAMESPACE(MGauge)

class PureGaugeEmPar : Serializable {
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PureGaugeEmPar,
                                    double,      scale,
                                    std::string, seedString);
};

class TPureGaugeEm : public Module<PureGaugeEmPar> {
public:
    typedef PhotonR::GaugeField     EmField;
    typedef PhotonR::GaugeLinkField EmComp;

public:
    TPureGaugeEm(const std::string name);
    virtual ~TPureGaugeEm(void) {};
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual void setup(void);
    virtual void execute(void);
};

MODULE_REGISTER(PureGaugeEm, TPureGaugeEm, MGauge);

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MGauge_PureGaugeEm_hpp_
