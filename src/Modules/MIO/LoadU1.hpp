/*
 * LoadU1.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef HadronsMILC_MIO_LoadU1_hpp_
#define HadronsMILC_MIO_LoadU1_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 Load a MILC U(1) photon field A_mu(x) in SciDAC/QIO format

 file   Name of the MILC u1 file to read.

 Written by MILC's u1_gauge application through
   save_real_scidac_from_field(filename, "", "", QIO_SINGLEFILE, u1_A, 4)
 (generic_u1/io_u1lat.c:330), which produces a QIO record of Nd=4
 single-precision reals per site:

   scidac-private-record-xml : datatype QLA_F_Real, precision F,
                               typesize 4, datacount 4
   scidac-binary-data        : 4*4 bytes per site, lexicographic site order,
                               network (big-endian) byte order
   scidac-checksum           : suma/sumb, verified on read

 Note the MILC input file must use "save_u1_serial_scidac"; the plain
 "save_u1_serial" option writes MILC's own binary format, which this does not
 read.

 The stored field is bare A_mu. The physical charge Qe multiplies it when the
 links are phased (MGauge::U1Phase), not here, so one A field serves every
 charge.

 The output type is PhotonR::GaugeField, the same type MGauge::StochEmFunc
 produces, so a MILC-generated photon and a Grid-generated one are
 interchangeable downstream. A_mu is real and is stored in the real part with
 zero imaginary part, mirroring Grid's Photon::StochasticField, which likewise
 finishes with "out = real(out)".
 ******************************************************************************/

BEGIN_MODULE_NAMESPACE(MIO)

/******************************************************************************
 *                          LoadU1 module parameters                          *
 ******************************************************************************/
class LoadU1Par : Serializable {
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadU1Par,
                                    std::string, file);
};

/******************************************************************************
 *                             TLoadU1 module                                 *
 ******************************************************************************/
class TLoadU1 : public Module<LoadU1Par> {
public:
    typedef PhotonR::GaugeField EmField;
    // The on-disk record type. This is not a free choice: GridLimeReader's
    // readLimeLatticeBinaryObject asserts sizeof(sobj)*gsites == record bytes
    // and derives the word format from the lattice type itself
    // (getFormatString -> "IEEE32BIG" for a single-precision scalar type), so
    // only Nd single-precision reals per site matches the MILC record.
    typedef Lattice<iVector<iScalar<iScalar<vRealF>>, Nd>> U1AField;
public:
    // constructor
    TLoadU1(const std::string name);
    // destructor
    virtual ~TLoadU1(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER(LoadU1, TLoadU1, MIO);

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MIO_LoadU1_hpp_
