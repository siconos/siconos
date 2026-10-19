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
#ifndef CXD_STD_HPP
#define CXD_STD_HPP

#include <vector>
#include <array>
#include <tuple>

#include "cxd_traits.hpp"
#include "cxd_utility.hpp"

namespace cxd
{
  namespace traits
  {
    struct stl {};

    template<typename T, std::size_t N>
    struct env_t<std::array<T, N>>
    {
      using library = stl;
      using vector = std::array<T, N>;
      using scalar = typename vector::value_type;
      using indice = typename vector::size_type;
      using support = std::tuple<indice, indice, vector>;
      static constexpr std::size_t dimension = N;
    };

    template<typename T>
    struct env_t<std::vector<T>>
    {
      using library = stl;
      using vector = std::vector<T>;
      using scalar = typename vector::value_type;
      using indice = typename vector::size_type;
      using support = std::tuple<indice, indice, vector>;
      static constexpr std::size_t dimension = 2; // historical default
    };

  }
  namespace concepts
  {
    template<typename T>
    concept std_container = requires (T a)
    {
      typename T::value_type;
      { typename traits::env_t<T>::library() } -> std::convertible_to<traits::stl>;
    };

    // std::vector
    // cf https://en.cppreference.com/w/cpp/named_req/AllocatorAwareContainer
    template<typename T>
    concept std_container_with_allocator = std_container<T> && requires
    {
      typename T::allocator_type;
    };

    // std::array
    template<typename T>
    concept std_container_with_fixed_size = requires
    {
      { std::tuple_size<T>::value };
    };

    template<typename T>
    concept dynamic_floating_point_container =
      std_container_with_allocator<T> &&
      std::floating_point<typename T::value_type>;

    template<typename T>
    concept static_floating_point_container =
      std_container_with_fixed_size<T> && std::floating_point<typename T::value_type>;

    // vector of vectors
    template<typename T>
    concept std_container_of_container = requires
    {
      typename T::value_type;
      typename T::value_type::value_type;
    };

    template<typename T>
    concept container_of_floating_point_container =
      std_container_of_container<T> && std::floating_point<typename T::value_type::value_type>;

  }

  namespace num
  {

    namespace operators
    {
      static constexpr const auto operator -
      (const concepts::std_container auto& a, const concepts::std_container auto& b)
      {
        return cxd::util::transform(a,
                                    [&b](const auto i, const auto& ai)
                                    { return ai-b[i];});
      }

      static constexpr const auto operator -
      (const concepts::std_container auto&a)
      {
        return cxd::util::transform(a,
                                    [](const auto i, const auto& ai)
                                    { return -ai;});
      }


      static constexpr const auto operator +
      (const concepts::std_container auto &a, const concepts::std_container auto& b)
      {
        return cxd::util::transform(a,
                                    [&b](const auto i, const auto& ai)
                                    { return ai+b[i];});
      }

      static constexpr const auto operator *
      (const auto& v, const concepts::std_container auto& a)
      {
        return cxd::util::transform(a,
                                    [&v](const auto i, const auto& ai)
                                    { return v*ai;});
      }
    }
  }

  namespace traits
  {

    auto env(concepts::container_of_floating_point_container auto x)
      -> env_t<typename std::decay_t<decltype(x)>::value_type>{};

    auto env(concepts::dynamic_floating_point_container auto x)
      -> env_t<std::decay_t<decltype(x)>> {};
    auto env(concepts::static_floating_point_container auto x)
      -> env_t<std::decay_t<decltype(x)>> {};
  }

  namespace num
  {
    using namespace cxd::concepts;

    template<std_container_with_allocator T>
    constexpr const auto inner_prod(const T& v1,
                                    const T& v2)
    {
      using env = decltype(cxd::traits::env(v1));
      using scalar = typename env::scalar;
      return std::inner_product(v1.begin(), v1.end(), v2.begin(),
                                static_cast<scalar>(0.0));
    }

    template<std_container_with_fixed_size T>
    constexpr const auto inner_prod(const T& a,
                                    const T& b)
    {
      return util::sum_ct(a, b,
                          [](const auto& ai, const auto& bi)
                          { return ai*bi; });
    }
  }

}



#endif
