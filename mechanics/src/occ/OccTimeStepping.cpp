class OccBody;

#include <SiconosVisitables.hpp>
#undef SICONOS_VISITABLES

#define SICONOS_VISITABLES()                    \
  KERNEL_CLASSES()                              \
  REGISTER(OccBody)



#include "OccTimeStepping.hpp"
#include "OccBody.hpp"

#include <Model.hpp>
#include <NonSmoothDynamicalSystem.hpp>

#include <SiconosVisitor.hpp>

#define VISITOR_CLASSES()                       \
  REGISTER(OccBody)

#include <VisitorMaker.hpp>



using namespace Experimental;

struct UpdateShapes : public SiconosVisitor
{
  using SiconosVisitor::visit;

  template<typename T>
  void operator() (const T& ds)
  {
    const_cast<T&>(ds).updateShapes();
    const_cast<T&>(ds).updateContactShapes();
  }
};


void OccTimeStepping::updateWorldFromDS()
{
  DynamicalSystemsGraph& dsg = *_nsds->dynamicalSystems();
  DynamicalSystemsGraph::VIterator dsi, dsiend;
  std11::tie(dsi, dsiend) = dsg.vertices();

  Visitor< Classes < OccBody >, UpdateShapes >::Make up;

  for (; dsi != dsiend; ++dsi)
  {
    dsg.bundle(*dsi)->accept(up);
  }

}
