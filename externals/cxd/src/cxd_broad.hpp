/* CXD is a program dedicated to convex collision detection.
 * It is part of the Siconos Software.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
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
