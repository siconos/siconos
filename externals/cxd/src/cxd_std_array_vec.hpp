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
#ifndef CXD_STD_ARRAY_VEC_HPP
#define CXD_STD_ARRAY_VEC_HPP

#include <ranges>

namespace cxd
{
  namespace array_vec
  {
    template<typename T, size_t N>
    struct vector
    {
      vector() : data{} {};

      template<typename ...Ts>
      vector(Ts... values) : data{{(T)values...}} {};
      vector(const std::array<T, N>& v) : data(v) {};
      using type = std::array<T, N>;
      using size_type = type::size_type;
      using value_type = T; /* for container of containers concept */

      value_type operator [](size_type i) const
      {
        return data[i];
      }

      constexpr size_type size() const
      {
        return data.size();
      }
      std::array<T, N> data;
    };
  }
  namespace num
  {
    namespace operators
    {
      template<typename T, size_t N>
      static constexpr const array_vec::vector<T, N>
      operator -
      (const array_vec::vector<T, N>& a, const array_vec::vector<T, N>& b)
      {
        return a.data - b.data;
      }

      template<typename T, size_t N>
      static constexpr const array_vec::vector<T, N>
      operator -
      (const array_vec::vector<T, N>& a)
      {
        return -a.data;
      }

      template<typename T, size_t N>
      static constexpr const array_vec::vector<T, N>
      operator +
      (const array_vec::vector<T, N>&a, const array_vec::vector<T, N>& b)
      {
        return a.data + b.data;
      }

      template<typename T, size_t N>
      static constexpr const array_vec::vector<T, N>
      operator *
      (const auto& v, const array_vec::vector<T, N>& a)
      {
        return v*a.data;
      }

    }
  }

  namespace traits
  {
    struct array_vec_t {};
    template<typename T, std::size_t N>
    struct env_t<array_vec::vector<T, N>>
    {
      using library = array_vec_t;
      using vector = array_vec::vector<T, N>;
      using scalar = typename vector::type::value_type;
      using indice = typename vector::type::size_type;
      using support = std::tuple<indice, indice, vector>;
    };
  }

  namespace concepts
  {
    template<typename T>
    concept array_vec_vector = requires (T a)
    {
      typename T::type::value_type;
      { typename traits::env_t<T>::library() } -> std::convertible_to<traits::array_vec_t>;
    };

    template<typename T>
    concept container_of_array_vec_vector = requires(T a)
    {
      { a[0] } -> array_vec_vector;
    };
  }
  namespace traits
  {
    constexpr auto env(concepts::array_vec_vector auto x)
      -> env_t<typename std::decay_t<decltype(x)>>{};

    constexpr auto env(concepts::container_of_array_vec_vector auto v)
      -> env_t<typename std::decay_t<decltype(v[0])>::value_type>{};
  }
  namespace num
  {
    template<typename T, size_t N>
    constexpr const auto inner_prod(const array_vec::vector<T, N>& v1,
                                    const array_vec::vector<T, N>& v2)
      {
        using env = decltype(cxd::traits::env(v1));
        using scalar = typename env::scalar;

        // https://stackoverflow.com/questions/64408217/dot-product-of-2-vectorors-c
        auto tmp = cxd::util::transform(v1.data,
                                        [&v2](const auto i, const auto& v1i)
                                        { return v1i*v2.data[i]; });
        std::ranges::sort(tmp,
                          [](const auto& a, const auto&b)
                          { return fabs(a) < fabs(b); });
        return cxd::util::ifold_left(tmp, static_cast<scalar>(0.0), [&]
                                     (const auto& state, const auto i, const auto& tmpi)
        { return state + tmpi; });
      }
  }
}

#endif
