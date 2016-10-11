
#include "InteractionManager.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Simulation.hpp"
#include "Model.hpp"

#include "debug.h"

void InteractionManager::link(SP::Interaction inter,
                              SP::DynamicalSystem ds1,
                              SP::DynamicalSystem ds2)
{
  DEBUG_PRINTF("link interaction : %d\n", inter->number());

  _simulation->nonSmoothDynamicalSystem()->link(inter, ds1, ds2);

  _simulation->initializeInteraction(_simulation->nextTime(), inter);

  // Note FP : ds init should probably be done once and only once for
  // all ds (like in simulation->initialize()) but where/when?
  // Note SS : in InteractionManager::buildGraph()?
  unsigned int levelMinForInput = inter->lowerLevelForInput();
  unsigned int levelMaxForInput = inter->upperLevelForInput();
  bool has2DS = inter->has2Bodies();
  for (unsigned int k = levelMinForInput ; k < levelMaxForInput + 1; k++)
  {
    ds1->initializeNonSmoothInput(k);
    if(has2DS)
      ds2->initializeNonSmoothInput(k);
  }
}

void InteractionManager::unlink(SP::Interaction inter)
{
  _simulation->nonSmoothDynamicalSystem()->removeInteraction(inter);
}
