/* CXD is a program dedicated to convex collision detection.
 * It is part of the Siconos Software.
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
#ifndef CXD_BROAD_HPP
#define CXD_BROAD_HPP

#include <algorithm>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "cxd.hpp"
#include "cxd_aabb.hpp"

namespace cxd
{
  template<std::size_t Dim, typename Shape>
  constexpr auto bounding_box(const Shape& shape)
  {
    using env = decltype(env(shape));
    using scalar = typename env::scalar;

    if (shape.empty())
      throw std::invalid_argument("shape must not be empty");

    aabb<scalar, Dim> result;
    for (const auto& point : shape)
    {
      if (point.size() != Dim)
        throw std::invalid_argument("shape points must match requested dimension");
      for (std::size_t i = 0; i < Dim; ++i)
      {
        result.minimum[i] = std::min(result.minimum[i], point[i]);
        result.maximum[i] = std::max(result.maximum[i], point[i]);
      }
    }
    return result;
  }
}

#endif
