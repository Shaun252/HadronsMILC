/*
 * U1Phase.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef HadronsMILC_MGauge_U1Phase_hpp_
#define HadronsMILC_MGauge_U1Phase_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 Multiply thin gauge links by a U(1) photon phase

 gauge    Name of the (unphased, thin) SU(3) gauge field.
 emField  Name of the photon field A_mu(x), from MIO::LoadU1 (a MILC u1 field)
          or MGauge::StochEmFunc (generated in Grid) -- both produce
          PhotonR::GaugeField, so either may be used.
 charge   The full Qe product, following MILC's convention: A is stored bare and
          the physical charge is applied here, at insertion. Production values
          are ( 0.202 0.101 0.0 -0.101 -0.202 ), i.e. (2/3)e and (1/3)e with
          e = sqrt(4*pi*alpha) ~ 0.3028.

 Applies, per direction,

     U_mu(x) -> exp(i * charge * A_mu(x)) * U_mu(x)

 i.e. a colour-diagonal overall phase, matching MILC's u1phase_on()
 (generic_u1/u1link.c:32, c_scalar_mult_su3mat by ce_itheta(charge*A)). The
 phased link is unitary but NOT special unitary: det = exp(3*i*theta).

 The phase MUST be applied to the THIN links, before smearing: the fat link is
 a sum over paths carrying different products of phases, and the U(3)
 projection acts on that sum, so phasing fat links is not equivalent. Feed the
 output of this module to MGauge::HISQSmear, which redoes the fattening for
 each (photon, charge). Ordering among the KS phases and the boundary twist is
 irrelevant, since a scalar phase commutes with all of them; only the smearing
 must come after.

 charge = 0 is deliberately NOT special-cased: it must traverse the whole code
 path and come out as an exact no-op, which is the cheapest regression test of
 this module (production already runs ncharges=5 including 0.0).
 ******************************************************************************/

BEGIN_MODULE_NAMESPACE(MGauge)

/******************************************************************************
 *                         U1Phase module parameters                          *
 ******************************************************************************/
class U1PhasePar : Serializable {
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(U1PhasePar,
                                    std::string, gauge,
                                    std::string, emField,
                                    double,      charge);
};

/******************************************************************************
 *                            TU1Phase module                                 *
 ******************************************************************************/
template <typename GImpl> class TU1Phase : public Module<U1PhasePar> {
public:
    GAUGE_TYPE_ALIASES(GImpl, );
    typedef PhotonR::GaugeField     EmField;
    typedef PhotonR::GaugeLinkField EmComp;

public:
    // constructor
    TU1Phase(const std::string name);
    // destructor
    virtual ~TU1Phase(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(U1Phase, TU1Phase<GIMPL>, MGauge);

/******************************************************************************
 *                       TU1Phase implementation                              *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl>
TU1Phase<GImpl>::TU1Phase(const std::string name)
    : Module<U1PhasePar>(name) {}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl>
std::vector<std::string> TU1Phase<GImpl>::getInput(void) {
    return {par().gauge, par().emField};
}

template <typename GImpl>
std::vector<std::string> TU1Phase<GImpl>::getOutput(void) {
    return {getName()};
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl> void TU1Phase<GImpl>::setup(void) {
    envCreateLat(GaugeField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl> void TU1Phase<GImpl>::execute(void) {
    LOG(Message) << "Phasing gauge field '" << par().gauge << "' with photon '"
                 << par().emField << "', charge Qe = " << par().charge
                 << std::endl;

    auto &U    = envGet(GaugeField, par().gauge);
    auto &A    = envGet(EmField, par().emField);
    auto &Uout = envGet(GaugeField, getName());

    // Records which photon was actually applied, and is the statistic that
    // compares a MILC-generated photon against a Grid-generated one: <A^2> per
    // component = norm2/(Nd*V) is gauge-dependent and so is sensitive to the
    // normalisation convention of A.
    LOG(Message) << "Photon field norm2 = " << norm2(A) << std::endl;

    EmComp         Amu(A.Grid());
    EmComp         phase(U.Grid());
    GaugeLinkField Umu(U.Grid());

    for (unsigned int mu = 0; mu < Nd; ++mu) {
        Amu   = peekLorentz(A, mu);
        Umu   = peekLorentz(U, mu);
        phase = exp(timesI(par().charge*Amu));
        Umu   = phase*Umu;
        pokeLorentz(Uout, Umu, mu);
    }

    // Mirrors MILC's "CHECK U(1) PLAQ" diagnostic: at charge 0 the two
    // plaquettes must agree to the last digit.
    LOG(Message) << "Plaquette: unphased = "
                 << WilsonLoops<GImpl>::avgPlaquette(U) << ", phased = "
                 << WilsonLoops<GImpl>::avgPlaquette(Uout) << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MGauge_U1Phase_hpp_
