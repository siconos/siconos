
#include "SiconosBroadphase.hpp"
#include "BodyDS.hpp"

#include "Model.hpp"
#include "Interaction.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "TimeStepping.hpp"

#include <debug.h>

void SiconosBroadphase::link(SP::Interaction inter,
                             SP::DynamicalSystem ds1,
                             SP::DynamicalSystem ds2)
{
  DEBUG_PRINTF("link interaction : %d\n", inter->number());
  model()->nonSmoothDynamicalSystem()->link(inter, ds1, ds2);
  model()->simulation()->computeLevelsForInputAndOutput(inter);

  // Note FP : ds init should probably be done once and only once for
  // all ds (like in simulation->initialize()) but where/when?
  // Note SS : in SiconosBroadphase::buildGraph()?
  unsigned int levelMinForInput = inter->lowerLevelForInput();
  unsigned int levelMaxForInput = inter->upperLevelForInput();
  bool has2DS = inter->has2Bodies();
  for (unsigned int k = levelMinForInput ; k < levelMaxForInput + 1; k++)
  {
    ds1->initializeNonSmoothInput(k);
    if(has2DS)
      ds2->initializeNonSmoothInput(k);
  }

  SP::InteractionsGraph indexSet0 = model()->nonSmoothDynamicalSystem()->topology()->indexSet0();
  InteractionsGraph::VDescriptor ui = indexSet0->descriptor(inter);

  inter->initialize(model()->simulation()->nextTime(), indexSet0->properties(ui));
}

void SiconosBroadphase::unlink(SP::Interaction inter)
{
  model()->nonSmoothDynamicalSystem()->removeInteraction(inter);
}

void SiconosBroadphase::buildGraph(std::vector<SP::SiconosContactor> contactors)
{
  std::vector<SP::SiconosContactor>::iterator it;
  for (it=contactors.begin(); it!=contactors.end(); ++it)
  {
    buildGraph(*it);
  }
}
