/*
 * SparseSpinColorDiagonal.hpp, part of HadronsMILC
 *
 * Spin-colour diagonal noise diluted in space and time (SparseNoiseMILC from
 * DilutedNoise.hpp): each source is split into nSparseL^3 * nSparseT patterns,
 * pattern = one value of every coordinate (x_d == r_d, r_d < nSparseL in
 * space, < nSparseT in time), times colour. With nSparseL = L and
 * nSparseT = T every pattern is a single site, i.e. complete dilution, which
 * makes the high-mode propagator exact (used to validate StagDeltaDLoop).
 *
 * Output <name>_vec: index ((source * nPattern) + pattern) * Nc + colour.
 */
#ifndef HadronsMILC_MNoise_SparseSpinColorDiagonal_hpp_
#define HadronsMILC_MNoise_SparseSpinColorDiagonal_hpp_

#include <DilutedNoise.hpp>
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MNoise)

class SparseSpinColorDiagonalMILCPar : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(SparseSpinColorDiagonalMILCPar, unsigned int,
                                  nsrc, unsigned int, nSparseL, unsigned int,
                                  nSparseT);
};

template <typename FImpl>
class TSparseSpinColorDiagonalMILC
    : public Module<SparseSpinColorDiagonalMILCPar> {
public:
  FERM_TYPE_ALIASES(FImpl, );

public:
  TSparseSpinColorDiagonalMILC(const std::string name)
      : Module<SparseSpinColorDiagonalMILCPar>(name) {}
  virtual ~TSparseSpinColorDiagonalMILC(void) {};
  virtual std::vector<std::string> getInput(void) { return {}; }
  virtual std::vector<std::string> getOutput(void) {
    return {getName(), getName() + "_vec"};
  }
  virtual void setup(void);
  virtual void execute(void);
};

MODULE_REGISTER_TMP(StagSparseSpinColorDiagonal,
                    TSparseSpinColorDiagonalMILC<STAGIMPL>, MNoise);

template <typename FImpl>
void TSparseSpinColorDiagonalMILC<FImpl>::setup(void) {
  if (par().nSparseL == 0 || par().nSparseT == 0) {
    HADRONS_ERROR(Argument, "nSparseL and nSparseT must be > 0");
  }
  envCreateDerived(SpinColorDiagonalNoiseMILC<FImpl>, SparseNoiseMILC<FImpl>,
                   getName(), 1, envGetGrid(FermionField), par().nsrc,
                   par().nSparseL, par().nSparseT);
  envCreate(std::vector<FermionField>, getName() + "_vec", 1, 0,
            envGetGrid(FermionField));
}

template <typename FImpl>
void TSparseSpinColorDiagonalMILC<FImpl>::execute(void) {
  auto &noise = envGet(SpinColorDiagonalNoiseMILC<FImpl>, getName());
  LOG(Message) << "Generating sparse-diluted spin-color diagonal noise: "
               << par().nsrc << " sources, nSparseL " << par().nSparseL
               << ", nSparseT " << par().nSparseT << std::endl;
  noise.generateNoise(rng4d());

  auto &noisevec = envGet(std::vector<FermionField>, getName() + "_vec");
  int nferm = noise.fermSize();
  noisevec.resize(nferm, envGetGrid(FermionField));
  for (int i = 0; i < nferm; i++) {
    noisevec[i] = noise.getFerm(i);
  }
  LOG(Message) << nferm << " noise vectors" << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HadronsMILC_MNoise_SparseSpinColorDiagonal_hpp_
