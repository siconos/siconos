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
