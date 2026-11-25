/* CXD is a program dedicated to convex collision detection.
 * It is part of the Siconos Software.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef CXD_CONCEPTS_HPP
#define CXD_CONCEPTS_HPP
#include <concepts>
#include <type_traits>

#include "cxd_traits.hpp"

namespace cxd
{
  namespace num { namespace operators {}}

  namespace concepts
  {
    using namespace cxd::num::operators;

    // after all env specializations
    template<typename T>
    concept known = requires (T a)
    {
      typename decltype(traits::env(a))::indice;
      typename decltype(traits::env(a))::scalar;
      typename decltype(traits::env(a))::vector;
      typename decltype(traits::env(a))::support;
    };

    template<typename T>
    concept any_vector = known<T> && std::is_floating_point_v<typename decltype(traits::env(std::declval<std::remove_cvref_t<T>>()))::scalar> && requires (T v)
    {
      { v.size() } ->
      std::convertible_to<typename decltype(traits::env(v))::indice>;
      { v[0] } ->
      std::convertible_to<typename decltype(traits::env(v))::scalar>;
      { v + v } -> known;
      { v - v } -> known;
      { - v } -> known;
      { typename decltype(traits::env(v))::scalar{1.} * v } -> known;
    };

    template<typename T>
    concept any_vectors = requires (T a)
    {
      { typename decltype(traits::env(a))::vector() } -> any_vector;
      { a[0] } -> any_vector; // [] operator
    };

    // Compile-time dimension-constrained concepts
    template<typename T>
    concept vector_2d = any_vector<T> &&
      requires { requires decltype(traits::env(std::declval<T>()))::dimension == 2; };

    template<typename T>
    concept vector_3d = any_vector<T> &&
      requires { requires decltype(traits::env(std::declval<T>()))::dimension == 3; };

    // For shapes (containers of vectors), check the INNER vector's dimension
    template<typename T>
    concept shape_2d = any_vectors<T> &&
      vector_2d<std::remove_cvref_t<decltype(std::declval<T>()[0])>>;

    template<typename T>
    concept shape_3d = any_vectors<T> &&
      vector_3d<std::remove_cvref_t<decltype(std::declval<T>()[0])>>;

    template<typename T>
    concept any_scalar = requires (T a)
    {
      { a } -> std::convertible_to<float>;
    };

    template<typename T>
    concept any_indice = requires (T i)
    {
      { i } -> std::convertible_to<std::size_t>;
    };

    template<typename T>
    concept any_support = requires (T a)
    {
      { std::get<0>(a) } -> any_indice;
      { std::get<1>(a) } -> any_indice;
      { std::get<2>(a) } -> any_vector;
    };
  }
}


#endif
