/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cassert>

#include "siconos/vkernel/algebra/linear_algebra.hpp"

namespace siconos::vkernel::algebra {

template <typename T, typename M>
decltype(auto) cast(T, M&& m) {
  return m.template cast<typename T::type>();
}

template <typename T, size_t M, size_t N>
using matrix = Eigen::Matrix<T, M, N>;  // column storage

template <typename T, size_t N>
using unbounded_col_matrix = Eigen::Matrix<T, Eigen::Dynamic, N, Eigen::RowMajor>;

template <typename T>
using unbounded_matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template <typename T>
using unbounded_vector = Eigen::Vector<T, Eigen::Dynamic>;

template <typename T, size_t M>
using vector = Eigen::Vector<T, M>;  // column vector

template <match::any_matrix T>
using matrix_view = Eigen::Map<T>;

template <typename T>
using matrix_ref = Eigen::Ref<T>;

template <typename T, size_t M>
using diagonal_matrix = Eigen::DiagonalMatrix<T, M>;

template <typename T>
using sparse_matrix = Eigen::SparseMatrix<T>;

template <typename T>
struct value_type {
  using type = decltype([]<bool flag = false>() {
    if constexpr (requires(T m) { typename T::value_type; }) {
      return typename T::value_type{};  // ok with gcc!
    } else {
      static_assert(flag, "no value_type");
    }
  }());
};

template <size_t M>
static constexpr decltype(auto) head(match::vector auto v) {
  return v.template head<M>(v);
}
// template specialization ok with clang, fails with gcc:
//
// template <typename T, size_t M, size_t N>
// struct value_type<Eigen::Matrix<T, M, N>> {
//   using type = typename Eigen::template Matrix<T, M, N>::value_type;
// };

// template <typename T, size_t M>
// struct value_type<Eigen::DiagonalMatrix<T, M>> {
//   using type = typename Eigen::template DiagonalMatrix<T, M>::Scalar;
// };

template <typename A>
using trans_t =
    matrix<typename value_type<A>::type, A::ColsAtCompileTime, A::RowsAtCompileTime>;

template <typename A, typename B>
using prod_t =
    matrix<typename value_type<B>::type, A::RowsAtCompileTime, B::ColsAtCompileTime>;

template <typename T>
static constexpr decltype(auto) nrows(T& value) {
  namespace match = siconos::vkernel::storage::pattern::match;

  if constexpr (match::any_matrix<T>) {
    return value.rows();  // specific to Eigen
  } else {
    []<typename Attr = T, bool flag = false>() { static_assert(flag, "no value type"); }();
  }
}

template <typename T>
static constexpr decltype(auto) ncols(T& value) {
  namespace match = siconos::vkernel::storage::pattern::match;

  if constexpr (match::any_matrix<T>) {
    return value.cols();  // specific to Eigen
  } else {
    []<typename Attr = T, bool flag = false>() { static_assert(flag, "no value type"); }();
  }
}

void set_zero(match::matrix auto& m) { m.setZero(); }

decltype(auto) dot(match::vector auto& v, match::vector auto& w) { return v.dot(w); }

void solve_in_place(match::diagonal_matrix auto& m, auto& b) { b = m.inverse() * b; }

void solve_in_place(match::matrix auto& m, auto& b) { assert(false); }

void solve_linear_system(match::matrix auto& m, auto& b, auto& c) {
  auto lu = m.partialPivLu();
  c = lu.solve(b);
}

// General sparse solver: c = m^{-1} * b
template <match::sparse_matrix M>
void solve_linear_system(M& m, auto& b, auto& c) {
  Eigen::SparseLU<std::remove_cvref_t<M>> solver;
  solver.compute(m);
  c = solver.solve(b);
}

// Sparse in-place solver: b = m^{-1} * b
template <match::sparse_matrix M>
void solve_in_place(M& m, auto& b) {
  Eigen::SparseLU<std::remove_cvref_t<M>> solver;
  solver.compute(m);
  b = solver.solve(b);
}

auto rows(auto& m) { return m.rows(); };
auto cols(auto& m) { return m.cols(); };

auto norm(auto& m) { return m.norm(); };

void clear(auto& m) { m.resize(0); };
}  // namespace siconos::vkernel::algebra
