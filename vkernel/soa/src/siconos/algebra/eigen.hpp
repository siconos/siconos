#pragma once

#include <Eigen/Dense>

#include "siconos/algebra/linear_algebra.hpp"

// concepts
namespace siconos::storage::pattern::match {

template <typename T>
concept matrix = requires(T m) { m(0, 0); };

template <typename T>
concept vector = matrix<T> && T::ColsAtCompileTime == 1;

// Is there a better way to check if a matrix from the Eigen C++ library  is
// diagonal :
template <typename T>
concept diagonal_matrix = !matrix<T> && requires(T m) { m.diagonal()[0]; };

// template <typename T>
// concept unbounded_matrix = requires(T m) {
//   { m.rows() } -> std::convertible_to<int>;
//   { m.cols() } -> std::convertible_to<int>;
// };

// template <typename T>
// concept unbounded_vector = requires(T m) {
//   { m.rows() } -> std::convertible_to<int>;
// } && T::ColsAtCompileTime == 1;

template <typename T>
concept any_matrix = (diagonal_matrix<T> || matrix<T>);

// template <typename T>
// concept fixed_size_matrix = any_matrix<T> && requires {
//   T::RowsAtCompileTime != Eigen::Dynamic;
//   T::ColsAtCompileTime != Eigen::Dynamic;
// };

template <typename T>
concept fixed_size_matrix =
    requires(T matrix) {
      // Check if the matrix has fixed rows and columns
      { matrix.rows() } -> std::convertible_to<std::size_t>;
      { matrix.cols() } -> std::convertible_to<std::size_t>;
    } && (T::RowsAtCompileTime != Eigen::Dynamic) &&
    (T::ColsAtCompileTime != Eigen::Dynamic);

template <typename T>
concept variable_size_matrix = !fixed_size_matrix<T>;

template <typename T>
concept fixed_size_vector = fixed_size_matrix<T> && vector<T>;

template <typename T>
concept variable_size_vector = vector<T> && !fixed_size_matrix<T>;

}  // namespace siconos::storage::pattern::match

namespace siconos::algebra {
namespace match = siconos::storage::pattern::match;

template <typename T, size_t M, size_t N>
using matrix = Eigen::Matrix<T, M, N>;  // column storage

template <typename T, size_t N>
using unbounded_col_matrix =
    Eigen::Matrix<T, Eigen::Dynamic, N, Eigen::RowMajor>;

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
struct value_type {
  using type = decltype([]<bool flag = false>() {
    if constexpr (requires(T m) { typename T::value_type; }) {
      return typename T::value_type{};  // ok with gcc!
    }
    else {
      static_assert(flag, "no value_type");
    }
  }());
};

template <size_t M>
static constexpr decltype(auto) head(match::vector auto v)
{
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
using trans_t = matrix<typename value_type<A>::type, A::ColsAtCompileTime,
                       A::RowsAtCompileTime>;

template <typename A, typename B>
using prod_t = matrix<typename value_type<B>::type, A::RowsAtCompileTime,
                      B::ColsAtCompileTime>;


template <typename T>
static constexpr decltype(auto) nrows(T)
{
  namespace match = siconos::storage::pattern::match;

  if constexpr (match::matrix<T> || match::diagonal_matrix<T>) {
    return T::RowsAtCompileTime;  // specific to Eigen
  }
  else {
    []<typename Attr = T, bool flag = false>() {
      static_assert(flag, "no value type");
    }();
  }
}

template <typename T>
static constexpr decltype(auto) ncols(T)
{
  namespace match = siconos::storage::pattern::match;

  if constexpr (match::matrix<T> || match::diagonal_matrix<T>) {
    return T::ColsAtCompileTime;  // specific to Eigen
  }
  else {
    []<typename Attr = T, bool flag = false>() {
      static_assert(flag, "no value type");
    }();
  }
}

void set_zero(match::matrix auto& m) { m.setZero(); };

decltype(auto) dot(match::vector auto& v, match::vector auto& w)
{
  return v.dot(w);
};

void solve_in_place(match::diagonal_matrix auto& m, match::vector auto& b)
{
  b = m.inverse() * b;
}

void solve_in_place(match::matrix auto& m, match::vector auto& b)
{
  auto lu = m.partialPivLu();
  b = lu.solve(b);
}
}  // namespace siconos::algebra
