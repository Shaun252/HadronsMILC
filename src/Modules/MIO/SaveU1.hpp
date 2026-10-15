/*
 * SaveU1.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef HadronsMILC_MIO_SaveU1_hpp_
#define HadronsMILC_MIO_SaveU1_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 Save a U(1) photon field A_mu(x) in MILC's SciDAC/QIO format (inverse of
 MIO::LoadU1)

 emField   Name of the photon field object (PhotonR::GaugeField, e.g. from
           MGauge::StochasticQedL or MIO::LoadU1).
 fileStem  Namestem of the file to write. The current trajectory number is
           appended as "<fileStem>.<traj>".

 Writes the same LIME layout as MILC's
   save_real_scidac_from_field(filename, "", "", QIO_SINGLEFILE, u1_A, 4)
 (generic_u1/io_u1lat.c), which MILC's reload_u1_serial reads back with
 restore_real_scidac_to_field:

   message 1: scidac-private-file-xml, scidac-file-xml (empty)
   message 2: scidac-private-record-xml : datatype QLA_F_Real, precision F,
                                          typesize 4, datacount 4
              scidac-record-xml (empty)
              scidac-binary-data        : 4*4 bytes per site, lexicographic
                                          site order, big-endian
              scidac-checksum           : suma/sumb

 Grid's ScidacWriter is not used because it adds Grid-only records
 (FieldMetaData, field norm) and names the binary record ildg-binary-data.

 The stored field is bare A_mu (the real part of the photon field; the
 imaginary part is zero, see MIO::LoadU1). The charge enters only when the
 links are phased (MGauge::U1Phase in Grid, u1phase_on in MILC).
 ******************************************************************************/

BEGIN_MODULE_NAMESPACE(MIO)

/******************************************************************************
 *                          SaveU1 module parameters                          *
 ******************************************************************************/
class SaveU1Par : Serializable {
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SaveU1Par,
                                    std::string, emField,
                                    std::string, fileStem);
};

/******************************************************************************
 *                             TSaveU1 module                                 *
 ******************************************************************************/
class TSaveU1 : public Module<SaveU1Par> {
public:
    typedef PhotonR::GaugeField EmField;
    // On-disk record type: Nd single-precision reals per site, as in MIO::LoadU1
    typedef Lattice<iVector<iScalar<iScalar<vRealF>>, Nd>> U1AField;
public:
    // constructor
    TSaveU1(const std::string name);
    // destructor
    virtual ~TSaveU1(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER(SaveU1, TSaveU1, MIO);

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MIO_SaveU1_hpp_
