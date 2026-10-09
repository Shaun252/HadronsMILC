#include "DeltaDLoop.hpp"

using namespace Grid;
using namespace Hadrons;
using namespace MContraction;

template class Grid::Hadrons::MContraction::TDeltaDLoopMILC<
    STAGIMPL, MassShiftEigenPack<STAGIMPL>>;
