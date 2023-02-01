#ifndef CABLETOOLS
#define CABLETOOLS
#include <memory>
#include <vector>

class SiconosVector;
class Point;

namespace siconos::mechanics::cables::tools {

void pointsToSiconosVector(const std::vector<Point> vecin,
                           std::shared_ptr<SiconosVector> vecout);

} // namespace siconos::mechanics::cables::tools

#endif
