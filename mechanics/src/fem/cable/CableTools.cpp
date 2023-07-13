#include "CableTools.h"

#include "Point.h"
#include "SiconosVector.hpp"

void siconos::fem::cable::tools::pointsToSiconosVector(
    const std::vector<siconos::fem::cable::Point> vecin,
    std::shared_ptr<siconos::algebra::SiconosVector> vecout) {
  assert(vecout);
  assert(vecin.size() * 3 == vecout->size());
  size_t i = 0;
  for (auto &point : vecin) {
    vecout->setValue(i++, point.x);
    vecout->setValue(i++, point.y);
    vecout->setValue(i++, point.z);
  }
}
