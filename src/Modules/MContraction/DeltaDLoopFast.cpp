#include "DeltaDLoopFast.hpp"

using namespace Grid;
using namespace Hadrons;
using namespace MContraction;

template class Grid::Hadrons::MContraction::TDeltaDLoopFastMILC<
    STAGIMPL, MassShiftEigenPack<STAGIMPL>>;
