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
#ifndef CXD_EIGEN_HPP
#define CXD_EIGEN_HPP

#include <concepts>
#include <tuple>
#include <type_traits>

#include <Eigen/Core>

#include "cxd_traits.hpp"

namespace cxd
{
  namespace concepts
  {
    template<typename T>
    concept eigen_vector = requires
    {
      typename std::decay_t<T>::Scalar;
      requires std::floating_point<typename std::decay_t<T>::Scalar>;
      requires (std::decay_t<T>::RowsAtCompileTime == 1 ||
                std::decay_t<T>::ColsAtCompileTime == 1);
    };
  }

  namespace traits
  {
    template<concepts::eigen_vector T>
    struct env_t<T>
    {
      using vector = T;
      using scalar = typename vector::Scalar;
      using indice = Eigen::Index;
      using support = std::tuple<indice, indice, vector>;
      static constexpr std::size_t dimension =
        (vector::RowsAtCompileTime == 1) ? vector::ColsAtCompileTime : vector::RowsAtCompileTime;
    };

    constexpr auto env(concepts::eigen_vector auto v)
      -> env_t<std::decay_t<decltype(v)>>{};
  }

  namespace num
  {
    constexpr const auto inner_prod(const concepts::eigen_vector auto& a,
                                    const concepts::eigen_vector auto& b)
    {
      return a.dot(b);
    }
  }
}

#endif
