/*
 * SaveQudaEigenPackMILC.cpp, part of HadronsMILC
 */

/*  END LEGAL */
#include "SaveQudaEigenPack.hpp"

using namespace Grid;
using namespace Hadrons;
using namespace MUtilities;

template class Grid::Hadrons::MUtilities::TSaveQudaEigenPack<
    STAGIMPL, BaseFermionEigenPack<STAGIMPL>>;
