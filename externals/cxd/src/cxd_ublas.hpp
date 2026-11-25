/* CXD is a program dedicated to convex collision detection.
 * It is part of the Siconos Software.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef CXD_UBLAS_HPP
#define CXD_UBLAS_HPP

#include <tuple> // for tuple_size_v
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/functional.hpp>

namespace cxd
{
  namespace concepts
  {
    template<typename T>
    concept ublas_vector = requires (T v)
    {
      { boost::numeric::ublas::inner_prod(v, v) };
    };
  }
  namespace traits
  {
    template<concepts::ublas_vector T>
    struct env_t<T>
    {
      using vector = T;
      using scalar = typename vector::value_type;
      using indice = typename vector::size_type;
      using support = std::tuple<indice, indice, vector>;
      static constexpr std::size_t dimension = 0;
    };

    template<typename T, std::size_t N>
    struct env_t<boost::numeric::ublas::fixed_vector<T, N>>
    {
      using vector = boost::numeric::ublas::fixed_vector<T, N>;
      using scalar = typename vector::value_type;
      using indice = typename vector::size_type;
      using support = std::tuple<indice, indice, vector>;
      static constexpr std::size_t dimension = N;
    };

    constexpr auto env(concepts::ublas_vector auto v) -> env_t<std::decay_t<decltype(v)>>{};
  };
  namespace num
  {
    constexpr const auto inner_prod(const concepts::ublas_vector auto& a,
                                    const concepts::ublas_vector auto& b)
    {
      return boost::numeric::ublas::inner_prod(a, b);
    };
  }
}
#endif
