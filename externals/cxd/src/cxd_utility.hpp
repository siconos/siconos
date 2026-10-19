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
#ifndef CXD_UTILITY_HPP
#define CXD_UTILITY_HPP

#include <utility>
#include <algorithm>
#include <initializer_list>
#include <type_traits>
#include <numeric>

#include <boost/hana.hpp>
#include <boost/hana/fwd/fold_left.hpp>
#include <boost/hana/fwd/transform.hpp>
#include <boost/hana/fwd/zip.hpp>
#include <boost/hana/ext/std/array.hpp>
#include <boost/hana/ext/std/integer_sequence.hpp>
#include "boost/hana/ext/boost/ublas.hpp"

namespace cxd
{
  namespace util
  {
    // compile-time iteration on one array
    // f(a[i]), i= 1..N
    static constexpr void iter_ct(const auto& a, auto&& f)
    {
      using array_type = std::decay_t<decltype(a)>;
      using size_type = typename array_type::size_type;

      array_type ta; // 'a' passed as const ref is not a constant expression
      [&f, &a]<size_type... I>(std::index_sequence<I...>)
      {
        (f(a[I]), ... );
      }(std::make_integer_sequence<size_type, std::size(ta)>{});
    }
    // compile-time iteration on 2 array
    // f(a[i],b[i]), i = 1..N
    static constexpr void iter_ct(const auto& a, const auto& b, auto&& f)
    {
      using array_type = std::decay_t<decltype(a)>;
      using size_type = typename array_type::size_type;

      array_type ta; // 'a' passed as const ref is not a constant expression
      [&f, &a, &b]<size_type... I>(std::index_sequence<I...>)
      {
        ((f(a[I], b[I])), ... );
      }(std::make_integer_sequence<size_type, std::size(ta)>{});
    }

    // compile-time sum
    // sum(f(a[i],b[i])), i = 1..N
    static constexpr auto sum_ct(const auto& a, const auto& b, auto&& f)
    {
      using array_type = std::decay_t<decltype(a)>;
      using size_type = typename array_type::size_type;

      array_type ta; // 'a' passed as const ref is not a constant expression
      return
        [&f,&a,&b]<size_type... I>(std::index_sequence<I...>)
      {
        return ((f(a[I], b[I])) + ... );
      }(std::make_integer_sequence<size_type, std::size(ta)>{});
    }

    // fold with indices
    static constexpr decltype(auto) ifold_left(const auto& array,
                                               const auto& initial_state,
                                               const auto&& func)
    {
      using array_type = std::decay_t<decltype(array)>;
      using size_type = typename array_type::size_type;

      if constexpr (boost::hana::Foldable<array_type>::value)
      // length known at compile-time
      {
        array_type ta; // 'array' passed as const ref is not a constant expression
        using indices = std::make_integer_sequence<size_type,
                                                   std::size(ta)>;
        return
          boost::hana::fold_left(indices{},
                                 initial_state,
                                 [&array,&func](decltype(initial_state) state,
                                                size_type i)
                                 {
                                   return func(state, i, array[i]);
                                 });
      }
      else if constexpr (
        std::copy_constructible<array_type> &&
        std::equality_comparable<decltype(array.begin())> &&
        std::input_iterator<decltype(array.begin())>)
      {
        return std::accumulate(array.begin(), array.end(), initial_state,
                               [&array,&func](decltype(initial_state) state,
                                              const auto& x)
                               {
                                 size_type i = &x - &array[0];
                                 return func(state, i, std::move(x));
                               });
      }
      else
      {
        // cf https://stackoverflow.com/questions/38304847/constexpr-if-and-static-assert
        []<bool flag = false>()
          {
            static_assert(flag, "cannot fold_left with these parameters");
          }();
      }
    }

    static constexpr const auto imin(const auto& array,
                                     const auto&& func)
    {
      using array_type = std::decay_t<decltype(array)>;
      using value_type = typename array_type::value_type;
      using size_type = typename array_type::size_type;
      return std::get<1>(ifold_left(array,
                                    std::tuple{std::numeric_limits<value_type>::infinity(),
                                      size_type{0}},
                                    [&func](const auto& state,
                                            const size_type i,
                                            const auto& x)
                                    {
                                      auto v = func(std::get<0>(state), x);
                                      return  v ? std::tuple{v, i} : state;
                                    }));
    }

    static constexpr const auto imax(const auto& array,
                                     const auto&& func)
    {
      using array_type = std::decay_t<decltype(array)>;
      using value_type = typename array_type::value_type;
      using size_type = typename array_type::size_type;
      return std::get<1>(ifold_left(array,
                                    -std::tuple{std::numeric_limits<value_type>::infinity(),
                                      size_type{0}},
                                    [&func](const auto& state,
                                            const size_type i,
                                            const auto& x)
                                    {
                                      auto v = func(std::get<0>(state), x);
                                      return  v ? std::tuple{v, i} : state;
                                    }));
    }


    // compile-time transform
    static constexpr const auto transform_ct(const auto& a, auto&& f)
    {
      using array_type = std::decay_t<decltype(a)>;
      using size_type = typename array_type::size_type;
      array_type ta; // 'a' passed as const ref is not a constant expression
      return
        [&f, &a]<size_type... I>(std::index_sequence<I...>)
        {
          return (array_type{f(I, a[I])... });
        }(std::make_integer_sequence<size_type, std::size(ta)>{});
    }

    static constexpr const auto transform(const auto& array,
                                          auto&& func)
    {
      using array_type = std::decay_t<decltype(array)>;
      if constexpr (boost::hana::Foldable<array_type>::value)
      {
        return transform_ct(array, func);
      }
      else if constexpr (
        std::equality_comparable<decltype(array.begin())> &&
        std::input_iterator<decltype(array.begin())>)
      {
        array_type res;
        std::transform(array.begin(), array.end(),
                       std::back_inserter(res),
                       [&array,&func](const auto& x)
                       {
                         auto i = &x - &array[0];
                         return func(i, x);
                       });
        return res;
      }
      else
      {
        // cf https://stackoverflow.com/questions/38304847/constexpr-if-and-static-assert
        []<bool flag = false>()
          {
            static_assert(flag, "cannot transform with these parameters");
          }();
      }
    }
  }
}
#endif
