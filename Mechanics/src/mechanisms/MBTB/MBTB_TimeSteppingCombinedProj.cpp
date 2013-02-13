#include "MBTB_TimeSteppingCombinedProj.hpp"
#include "MBTB_PYTHON_API.hpp"
#include "MBTB_internalTool.hpp"
#include <boost/math/quaternion.hpp>
//#define TSPROJ_DEBUG

void MBTB_TimeSteppingCombinedProj::updateWorldFromDS()
{
#ifdef TSPROJ_DEBUG
  printf("MBTB_TimeSteppingCombinedProj::updateWordFromDS \n");
#endif
  MBTB_updateDSFromSiconos();
  _MBTB_updateContactFromDS();
}
