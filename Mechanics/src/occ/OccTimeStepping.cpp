#include "OccTimeStepping.hpp"

#include "OccBody.hpp"

#include <Model.hpp>
#include <NonSmoothDynamicalSystem.hpp>

struct UpdateContactShapes : public SiconosVisitor
{
  using SiconosVisitor::visit;

  void visit(const OccBody& ods)
  {
    ods.updateContactShapes();
  }
};


void OccTimeStepping::updateWorldFromDS()
{
  DynamicalSystemsGraph& dsg = *model()->nonSmoothDynamicalSystem()->dynamicalSystems();
  DynamicalSystemsGraph::VIterator dsi, dsiend;
  std11::tie(dsi, dsiend) = dsg.vertices();

  static UpdateContactShapes up;

  for (; dsi != dsiend; ++dsi)
  {
    dsg.bundle(*dsi)->accept(up);
  }

}
