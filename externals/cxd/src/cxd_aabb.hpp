/* CXD is a program dedicated to convex collision detection.
 * It is part of the Siconos Software.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef CXD_AABB_HPP
#define CXD_AABB_HPP

#include <array>
#include <cstddef>
#include <limits>
#include <algorithm>

namespace cxd
{
  template<typename Scalar, std::size_t Dim>
  struct aabb
  {
    std::array<Scalar, Dim> minimum;
    std::array<Scalar, Dim> maximum;

    constexpr aabb()
    {
      minimum.fill(std::numeric_limits<Scalar>::infinity());
      maximum.fill(-std::numeric_limits<Scalar>::infinity());
    }

    constexpr explicit aabb(const std::array<Scalar, Dim>& min,
                            const std::array<Scalar, Dim>& max)
      : minimum(min), maximum(max) {}

    constexpr bool overlaps(const aabb& other) const
    {
      for (std::size_t i = 0; i < Dim; ++i)
        if (maximum[i] < other.minimum[i] || other.maximum[i] < minimum[i])
          return false;
      return true;
    }

    constexpr bool contains(const aabb& other) const
    {
      for (std::size_t i = 0; i < Dim; ++i)
        if (minimum[i] > other.minimum[i] || maximum[i] < other.maximum[i])
          return false;
      return true;
    }

    template<std::size_t Axis>
    constexpr Scalar center() const
    {
      return (minimum[Axis] + maximum[Axis]) * Scalar{0.5};
    }

    constexpr aabb& expand(const aabb& other)
    {
      for (std::size_t i = 0; i < Dim; ++i)
      {
        minimum[i] = std::min(minimum[i], other.minimum[i]);
        maximum[i] = std::max(maximum[i], other.maximum[i]);
      }
      return *this;
    }
  };

  // No implicit default-dimension alias: callers must write aabb<Scalar, 2> or aabb<Scalar, 3>.
}

#endif
