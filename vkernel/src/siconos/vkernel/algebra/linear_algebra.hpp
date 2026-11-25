// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#pragma once

#include <concepts>
#include <type_traits>

#include "siconos/vkernel/storage/pattern/pattern.hpp"
#include "siconos/vkernel/storage/traits/traits.hpp"

// The algebra contract.
//
// Two layers share one namespace, siconos::vkernel::algebra, and two different
// kinds of value:
//
//   the per-system layer (algebra/eigen_backend.hpp)
//     small, fixed size, stack value types: the 3x3 mass block of a rigid
//     system, a 3-vector of velocities, the vertices of a polygon.  Eigen
//     types are used directly (algebra::matrix, algebra::vector,
//     algebra::diagonal_matrix, ...), so the 2x2 and 3x3 products are
//     compile-time sized.  The storage layer asks a type for its extents
//     (vnrows/vncols) and the concepts below are what it matches.
//
//   the assembled layer (algebra/numerics_backend.hpp)
//     one large sparse matrix for the whole simulation, indexed by dof, plus
//     its vectors.  The value types (algebra::mat<T>, algebra::vec<T>,
//     algebra::diag_mat<T>) wrap a NumericsMatrix* and an offset, so that a
//     system's block is a view into the assembled matrix (mat_view/vec_view)
//     and assembly is a copy into it.  These are the match::mat / match::vec /
//     match::any_mat / match::diag_mat concepts.
//
// The free functions are the contract; a backend is a header that provides
// them for its types.  For a matrix type M and a vector type V:
//
//   resize(M&, rows, cols)            setup(M&)               size0(M&)
//   size1(M&)                         raw_size0/raw_size1     nnz(M&)
//   set_value(M&, i, j, value)        insert(M&, const M&, r0, c0)
//   transpose(M&)                     prod(a, b, c)           prodt1(a, b, c)
//   prodt2(a, b, c)                   add(alpha, A, beta, B)
//   add(alpha, A, beta, B, C)         solve_linear_system(a, b, c)
//   solve_in_place(a, b)              solve_linear_system_with_transpose(d, b, c)
//   mat_view<M>(m, r0, c0)
//   resize(V&, n)                     size0(V&)               raw_size0(V&)
//   set_value(V&, i, value)           insert(V&, const V&, i0)
//   get_vector(V&, i)                 copy(a, b)              add(a, b)
//   scal(h, v)                        dot(v, w)               norm(v)
//   display(m)                        display(v)              rows/cols/clear
//
// set_value's value is a scalar, a diagonal matrix, a dense matrix or a sparse
// matrix: the block structure of the vkernel, which is why the per-system and
// the assembled layers meet there.
//
// What is NOT part of the contract: block_diagonal_structure,
// block_diagonal_solver, compute_kkt_matrix and solve_block_diagonal_system.
// These read the pattern of the assembled matrix (connected components, the
// factorization of each FEM system, the compliance slots) and are written for
// the Numerics/CSparse layout.  A backend that cannot give that structure must
// provide its own W = H A^-1 H^T, not reuse these.
namespace siconos::vkernel::storage::pattern::match {

// A compile-time extent equal to this value means "sized at run time".  It is
// the value of Eigen::Dynamic, and of any backend that sizes at run time.
inline constexpr auto dynamic_extent = -1;

// ---- the per-system layer: matrices of a system ----

template <typename T>
concept matrix = requires(T m) { m(0, 0); };

template <typename T>
concept vector_raw = matrix<T> && T::ColsAtCompileTime == 1;

template <typename T>
concept vector = vector_raw<std::decay_t<T>>;

// general sparse case
template <typename T>
concept sparse_matrix_raw = requires(T m) {
  typename T::Scalar;
  { m.rows() } -> std::convertible_to<std::size_t>;
  { m.cols() } -> std::convertible_to<std::size_t>;
  { m.nonZeros() } -> std::convertible_to<std::size_t>;
  { m.makeCompressed() } -> std::same_as<void>;
  // data pointers
  { m.innerIndexPtr() } -> std::convertible_to<void*>;
  { m.valuePtr() } -> std::convertible_to<typename T::Scalar*>;
} && !requires {
  // exclude fixed-size dense matrices
  requires T::RowsAtCompileTime != dynamic_extent&& T::ColsAtCompileTime != dynamic_extent;
};

template <typename T>
concept sparse_matrix = sparse_matrix_raw<std::decay_t<T>>;

// A diagonal matrix: it has a diagonal() and it is neither a dense nor a sparse
// matrix.  A sparse matrix also has .diagonal() (of type Eigen::Diagonal), which
// is a sparse-backed diagonal expression, not a diagonal storage: it must stay a
// sparse matrix (the dense storage wrappers bind it as one, and Python can cast
// Eigen::SparseMatrix but not Eigen::Diagonal).
template <typename T>
concept diagonal_matrix_raw =
    !matrix<T> && !sparse_matrix<T> && requires(T m) { m.diagonal()[0]; };

template <typename T>
concept diagonal_matrix = diagonal_matrix_raw<std::decay_t<T>>;

template <typename T>
concept any_matrix_raw = (diagonal_matrix<T> || matrix<T> || sparse_matrix<T>);

template <typename T>
concept any_matrix = any_matrix_raw<std::decay_t<T>>;

template <typename T>
concept fixed_size_matrix_raw = requires(T matrix) {
  // Check if the matrix has fixed rows and columns
  { matrix.rows() } -> std::convertible_to<std::size_t>;
  { matrix.cols() } -> std::convertible_to<std::size_t>;
} && (T::RowsAtCompileTime != dynamic_extent) && (T::ColsAtCompileTime != dynamic_extent);

template <typename T>
concept fixed_size_matrix = fixed_size_matrix_raw<std::decay_t<T>>;

template <typename T>
concept matrix_1x1_raw = matrix<T> && T::RowsAtCompileTime == 1 && T::ColsAtCompileTime == 1;

template <typename T>
concept matrix_1x1 = matrix_1x1_raw<std::decay_t<T>>;

template <typename T>
concept diagonal_or_1x1_raw = diagonal_matrix<T> || matrix_1x1<T>;

template <typename T>
concept diagonal_or_1x1 = diagonal_or_1x1_raw<std::decay_t<T>>;

template <typename T>
concept variable_size_matrix = !fixed_size_matrix<T>;

template <typename T>
concept fixed_size_vector = fixed_size_matrix<T> && vector<T>;

template <typename T>
concept variable_size_vector = vector<T> && !fixed_size_matrix<T>;

template <typename T>
concept unbounded = variable_size_vector<T> || variable_size_matrix<T>;

// ---- the assembled layer: the wrapper types of the whole simulation ----

template <typename T>
concept vec = requires { typename std::decay_t<T>::vec_t; };

template <typename T>
concept any_mat = !vec<T> && requires { typename std::decay_t<T>::any_mat_t; };

template <typename T>
concept mat = any_mat<T> && requires { typename std::decay_t<T>::mat_t; };

template <typename T>
concept diag_mat = any_mat<T> && requires { typename std::decay_t<T>::diag_mat_t; };

}  // namespace siconos::vkernel::storage::pattern::match

namespace siconos::vkernel::algebra {
namespace match = siconos::vkernel::storage::pattern::match;
}
