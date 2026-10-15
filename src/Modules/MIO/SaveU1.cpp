/*
 * SaveU1.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#include "SaveU1.hpp"

using namespace Grid;
using namespace Hadrons;
using namespace MIO;

/******************************************************************************
*                       TSaveU1 implementation                                *
******************************************************************************/
// raw LIME record with a given string payload (boss node only)
static void writeRawLimeRecord(GridLimeWriter &w, int MB, int ME,
                               const std::string &type, const std::string &data)
{
    if (!w.boss_node) return;
    uint64_t         nbytes = data.size();
    LimeRecordHeader *h     = limeCreateHeader(MB, ME,
                                  const_cast<char *>(type.c_str()), nbytes);
    GRID_ASSERT(h != NULL);
    GRID_ASSERT(limeWriteRecordHeader(h, w.LimeW) >= 0);
    GRID_ASSERT(limeWriteRecordData(const_cast<char *>(data.data()), &nbytes,
                                    w.LimeW) >= 0);
    GRID_ASSERT(limeWriterCloseRecord(w.LimeW) >= 0);
    limeDestroyHeader(h);
}

// constructor /////////////////////////////////////////////////////////////////
TSaveU1::TSaveU1(const std::string name)
: Module<SaveU1Par>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TSaveU1::getInput(void)
{
    std::vector<std::string> in = {par().emField};

    return in;
}

std::vector<std::string> TSaveU1::getOutput(void)
{
    std::vector<std::string> out;

    return out;
}

std::vector<std::string> TSaveU1::getOutputFiles(void)
{
    return {par().fileStem + "." + std::to_string(vm().getTrajectory())};
}

// setup ///////////////////////////////////////////////////////////////////////
void TSaveU1::setup(void)
{}

// execution ///////////////////////////////////////////////////////////////////
void TSaveU1::execute(void)
{
    std::string fileName = par().fileStem + "." + std::to_string(vm().getTrajectory());

    LOG(Message) << "Saving U(1) photon field '" << par().emField
                 << "' in MILC SciDAC format to file '" << fileName << "'"
                 << std::endl;

    auto &A = envGet(EmField, par().emField);

    // Single-precision real grid, built directly for the reason given in
    // MIO::LoadU1 (no red-black partner needed, small-lattice SIMD layout).
    auto          *g4 = A.Grid();
    GridCartesian  rgrid(g4->GlobalDimensions(),
                         GridDefaultSimd(Nd, vRealF::Nsimd()),
                         g4->ProcessorGrid());
    U1AField       Araw(&rgrid);

    // Narrow the complex field to its real part, site by site through the
    // layout-safe lexicographic arrays (the inverse of MIO::LoadU1).
    typedef typename U1AField::scalar_object U1ASite;
    typedef typename EmField::scalar_object  EmSite;

    std::vector<EmSite>  ebuf;
    std::vector<U1ASite> abuf;
    RealD                maxIm = 0.;

    unvectorizeToLexOrdArray(ebuf, A);
    abuf.resize(ebuf.size());
    for (size_t i = 0; i < ebuf.size(); ++i)
    {
        for (unsigned int mu = 0; mu < Nd; ++mu)
        {
            abuf[i](mu)()() = (RealF)real(ebuf[i](mu)()());
            maxIm = std::max(maxIm, std::abs(imag(ebuf[i](mu)()())));
        }
    }
    g4->GlobalMax(maxIm);
    if (maxIm != 0.)
    {
        HADRONS_ERROR(Argument, "photon field '" + par().emField
                      + "' has a non-zero imaginary part (max |Im A| = "
                      + std::to_string(maxIm) + ")");
    }
    vectorizeFromLexOrdArray(abuf, Araw);

    // MILC/QIO record headers
    std::stringstream dims;
    for (int d = 0; d < Nd; ++d) dims << g4->GlobalDimensions()[d] << " ";
    std::string fileXml =
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?><scidacFile><version>1.1"
        "</version><spacetime>" + std::to_string(Nd) + "</spacetime><dims>"
        + dims.str() + "</dims><volfmt>0</volfmt></scidacFile>";
    std::time_t now = std::time(nullptr);
    char        date[64];
    std::strftime(date, sizeof(date), "%a %b %e %H:%M:%S %Y UTC", std::gmtime(&now));
    std::string recordXml =
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?><scidacRecord><version>1.1"
        "</version><date>" + std::string(date) + "</date><recordtype>0"
        "</recordtype><datatype>QLA_F_Real</datatype><precision>F</precision>"
        "<typesize>4</typesize><datacount>" + std::to_string(Nd)
        + "</datacount></scidacRecord>";

    makeFileDir(fileName, g4);

    GridLimeWriter writer(g4->IsBoss());
    writer.open(fileName);
    writeRawLimeRecord(writer, 1, 0, SCIDAC_PRIVATE_FILE_XML, fileXml);
    writeRawLimeRecord(writer, 0, 1, SCIDAC_FILE_XML, std::string(1, '\0'));
    writeRawLimeRecord(writer, 1, 0, SCIDAC_PRIVATE_RECORD_XML, recordXml);
    writeRawLimeRecord(writer, 0, 0, SCIDAC_RECORD_XML, std::string(1, '\0'));

    // Binary record + checksum, as GridLimeWriter::writeLimeLatticeBinaryObject
    // but without its Grid-only field-norm record.
    uint64_t PayloadSize = sizeof(U1ASite) * rgrid._gsites;
    uint64_t offset1     = 0;
    uint32_t nersc_csum, scidac_csuma, scidac_csumb;

    if (writer.boss_node)
    {
        writer.createLimeRecordHeader(SCIDAC_BINARY_DATA, 0, 0, PayloadSize);
        fflush(writer.File);
        offset1 = ftello(writer.File);
    }
    rgrid.Broadcast(0, (void *)&offset1, sizeof(offset1));

    std::string format = getFormatString<typename U1AField::vector_object>();
    BinarySimpleMunger<U1ASite, U1ASite> munge;
    BinaryIO::writeLatticeObject<typename U1AField::vector_object, U1ASite>(
        Araw, fileName, munge, offset1, format, nersc_csum, scidac_csuma,
        scidac_csumb, BinaryIO::BINARYIO_LEXICOGRAPHIC | BinaryIO::BINARYIO_AGGREGATE);

    if (writer.boss_node)
    {
        fseek(writer.File, 0, SEEK_END);
        GRID_ASSERT((uint64_t)ftello(writer.File) - offset1 == PayloadSize);
        GRID_ASSERT(limeWriterCloseRecord(writer.LimeW) >= 0);
    }

    std::stringstream sa, sb;
    sa << std::hex << scidac_csuma;
    sb << std::hex << scidac_csumb;
    std::string checksumXml =
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?><scidacChecksum><version>1.0"
        "</version><suma>" + sa.str() + "</suma><sumb>" + sb.str()
        + "</sumb></scidacChecksum>";
    writeRawLimeRecord(writer, 0, 1, SCIDAC_CHECKSUM, checksumXml);
    writer.close();

    LOG(Message) << "Photon field norm2 = " << norm2(A) << ", format " << format
                 << ", suma " << sa.str() << " sumb " << sb.str() << std::endl;
}
