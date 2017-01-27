
#include "InteractionManager.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Simulation.hpp"
#include "Model.hpp"

#include "debug.h"

#include <algorithm>

void InteractionManager::insertNonSmoothLaw(SP::NonSmoothLaw nslaw,
                                            long unsigned int group1,
                                            long unsigned int group2)
{
  // ublas::matrix size type is not the same on 32 bits and 64 bits
  NSLawMatrix::size_type maxgroup = std::max((NSLawMatrix::size_type) group1,
                                             (NSLawMatrix::size_type) group2);
  _nslaws.resize( std::max(_nslaws.size1(), maxgroup+1) );
  _nslaws(group1, group2) = nslaw;
}

SP::NonSmoothLaw InteractionManager::nonSmoothLaw(long unsigned int group1,
                                                  long unsigned int group2)
{
  if (group1 < _nslaws.size1() && group2 < _nslaws.size2())
    return _nslaws(group1, group2);
  else
    return SP::NonSmoothLaw();
}
