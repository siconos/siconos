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
