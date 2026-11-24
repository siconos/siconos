/* CXD is a program dedicated to convex collision detection.
 * It is part of the Siconos Software.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
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
