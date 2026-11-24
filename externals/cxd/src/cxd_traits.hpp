/* CXD is a program dedicated to convex collision detection.
 * It is part of the Siconos Software.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef CXD_TRAITS_HPP
#define CXD_TRAITS_HPP
#include <type_traits>
#include <tuple>
#include <stdexcept>
#include <cstdlib>


namespace cxd
{
  namespace traits
  {
    template<typename T>
    struct ENVIRONMENT_IS_NOT_DEFINED_FOR_THAT_TYPE
    {
    };

    // for specialization
    template<typename T, typename R = bool>
    struct env_t
    {
      static constexpr std::size_t dimension = 0; // unknown/default
    };

    auto env(auto x) -> ENVIRONMENT_IS_NOT_DEFINED_FOR_THAT_TYPE<typename std::decay_t<decltype(x)>>
    {
    };
  }

  namespace concepts
  {
  }
}

#endif
