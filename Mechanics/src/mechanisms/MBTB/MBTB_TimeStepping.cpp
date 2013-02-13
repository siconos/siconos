#include "MBTB_TimeStepping.hpp"
#include "MBTB_PYTHON_API.hpp"
#include "MBTB_internalTool.hpp"
#include <boost/math/quaternion.hpp>
//#define TS_DEBUG
MBTB_TimeStepping::MBTB_TimeStepping(SP::TimeDiscretisation td,
                                     SP::OneStepIntegrator osi,
                                     SP::OneStepNSProblem osnspb_velo):TimeStepping(td,osi,osnspb_velo)
{
}


void MBTB_TimeStepping::updateWorldFromDS()
{
#ifdef TS_DEBUG
  printf("MBTB_TimeStepping::updateWordFromDS \n");
#endif
  MBTB_updateDSFromSiconos();
  _MBTB_updateContactFromDS();
}
