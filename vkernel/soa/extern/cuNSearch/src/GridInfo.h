#pragma once
#include "Types.h"

namespace cuNSearch {
struct GridInfo {
  // The grid information of a point set is computed in computeCellInformation,
  // which is skipped for an empty point set. Default initialize the members so
  // that the grid information of an empty point set is well defined (and not
  // some random garbage) even if it is read.
  cuNSearch::Real3 GridMin = {};
  uint ParticleCount = 0;
  cuNSearch::Real3 GridDelta = {};
  UInt3 GridDimension = {};
  UInt3 MetaGridDimension = {};
  float SquaredSearchRadius = 0.f;
};
}  // namespace cuNSearch
