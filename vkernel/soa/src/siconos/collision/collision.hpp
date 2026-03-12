#pragma once

#include "collision_head.hpp"

namespace siconos::collision {

/* 2D */
decltype(auto) distance_sq(match::vector auto& q1, match::vector auto& q2)
{
  const auto dx = q2(0) - q1(0);
  const auto dy = q2(1) - q1(1);
  return dx * dx + dy * dy;
}

decltype(auto) distance(match::vector auto& q1, match::vector auto& q2)
{
  return sqrt(distance_sq(q1, q2));
}

}  // namespace siconos::collision
