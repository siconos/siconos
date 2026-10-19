/* CXD is a program dedicated to convex collision detection
 * (currently part of the Siconos software).
 *
 * Copyright 2026 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
