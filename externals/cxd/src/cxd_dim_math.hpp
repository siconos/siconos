/* CXD is a program dedicated to convex collision detection.
 * It is part of the Siconos Software.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef CXD_DIM_MATH_HPP
#define CXD_DIM_MATH_HPP

#include <cstddef>
#include <type_traits>
#include "cxd_traits.hpp"
#include "cxd_utility.hpp"
#include "cxd_concepts.hpp"

namespace cxd
{
  using namespace cxd::concepts;

  // 3D cross product: returns vector orthogonal to a and b
  template<any_vector T>
  static inline constexpr auto cross(const T& a, const T& b)
  {
    return T{
      a[1] * b[2] - a[2] * b[1],
      a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0]
    };
  }

  // Dimension-aware triple product:
  // - 2D: (a x b) x c, returns vector
  // - 3D: c . (a x b), returns scalar
  template<any_vector T>
  static inline constexpr auto triple_product(const T& a,
                                              const T& b,
                                              const T& c)
  {
    if constexpr (decltype(env(a))::dimension == 2)
    {
      return (-num::inner_prod(c, b) * a + num::inner_prod(c, a) * b);
    }
    else
    {
      return num::inner_prod(c, cross(a, b));
    }
  }
}

#endif
