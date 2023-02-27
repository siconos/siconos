#include "CableTools.h"

#include "Point.h"
#include "SiconosVector.hpp"

void siconos::mechanics::cables::tools::pointsToSiconosVector(
    const std::vector<Point> vecin, std::shared_ptr<SiconosVector> vecout)
{
  assert(vecout);
  assert(vecin.size() * 3 == vecout->size());
  size_t i = 0;
  for (auto &point : vecin) {
    vecout->setValue(i++, point.x);
    vecout->setValue(i++, point.y);
    vecout->setValue(i++, point.z);
  }
}
