/*
 * SaveQudaEigenPackMILC.hpp, part of HadronsMILC
 *
 * Write a staggered eigenpack in QUDA/MILC partfile format.
 * The eigenvectors are taken from a BaseFermionEigenPack produced by
 * StagFermionIRL (red-black / checkerboard fields).  The eigenvalues
 * are converted to the MILC convention
 *
 *   eval_MILC[i] = ComplexD(2*mass,  sqrt(lambda_i))
 *
 * where lambda_i is the i-th eigenvalue of M^dag M stored in the pack.
 *
 * XML parameters
 * --------------
 *   eigenPack  : name of the BaseFermionEigenPack object in the environment
 *   mass       : bare quark mass (used for the MILC eigenvalue convention)
 *   output     : file-stem; each MPI rank writes <output>.vol<NNNN>
 *   singlePrec : write in single (float) precision (default true)
 */

#ifndef HadronsMILC_MUtilities_SaveQudaEigenPack_hpp_
#define HadronsMILC_MUtilities_SaveQudaEigenPack_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>
#include <EigenPack.hpp>   // AdaptorEigenPackMILC, writeQudaStaggeredEigenpack

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MUtilities)

class SaveQudaEigenPackPar : Serializable {
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SaveQudaEigenPackPar,
                                    std::string, eigenPack,
                                    double,      mass,
                                    std::string, output,
                                    bool,        singlePrec);
};

template <typename FImpl, typename Pack>
class TSaveQudaEigenPack : public Module<SaveQudaEigenPackPar> {
public:
    typedef typename Pack::Field  Field;
    typedef BaseEigenPack<Field>  BasePack;
    FERM_TYPE_ALIASES(FImpl, );

public:
    TSaveQudaEigenPack(const std::string name);
    virtual ~TSaveQudaEigenPack(void) {}

    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual void setup(void);
    virtual void execute(void);
};

MODULE_REGISTER_TMP(SaveQudaEigenPack,
    ARG(TSaveQudaEigenPack<STAGIMPL, BaseFermionEigenPack<STAGIMPL>>),
    MUtilities);

/******************************************************************************
 *                TSaveQudaEigenPack implementation                           *
 ******************************************************************************/

template <typename FImpl, typename Pack>
TSaveQudaEigenPack<FImpl, Pack>::TSaveQudaEigenPack(const std::string name)
    : Module<SaveQudaEigenPackPar>(name)
{}

template <typename FImpl, typename Pack>
std::vector<std::string> TSaveQudaEigenPack<FImpl, Pack>::getInput(void)
{
    return {par().eigenPack};
}

template <typename FImpl, typename Pack>
std::vector<std::string> TSaveQudaEigenPack<FImpl, Pack>::getOutput(void)
{
    return {};
}

template <typename FImpl, typename Pack>
void TSaveQudaEigenPack<FImpl, Pack>::setup(void)
{
    if (par().output.empty())
    {
        HADRONS_ERROR(Parsing, "SaveQudaEigenPack: 'output' must not be empty.");
    }
}

template <typename FImpl, typename Pack>
void TSaveQudaEigenPack<FImpl, Pack>::execute(void)
{
    auto &epack = envGet(BasePack, par().eigenPack);
    int   Nvec  = epack.evec.size();

    LOG(Message) << "Writing " << Nvec << " eigenvectors in QUDA/MILC partfile"
                 << " format to stem '" << par().output << "'"
                 << " (mass=" << par().mass
                 << ", " << (par().singlePrec ? "single" : "double")
                 << " precision)" << std::endl;

    // Build MILC-shifted complex eigenvalues: eval[i] = (2m, sqrt(lambda_i))
    AdaptorEigenPackMILC<Field> adaptor(epack.evec, epack.eval, par().mass);

    // Grid lives on the eigenvectors themselves (red-black grid from IRL)
    GridBase *grid = epack.evec[0].Grid();

    writeQudaStaggeredEigenpack<FImpl>(par().output, epack.evec,
                                       adaptor.eval, grid,
                                       par().singlePrec);

    LOG(Message) << "Done writing QUDA eigenpack." << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MUtilities_SaveQudaEigenPack_hpp_
