#pragma once

#include <siconos/numerics/CSparseMatrix_internal.h>  // for CSparseMatrix
#include <siconos/numerics/NumericsMatrix.h>
#include <siconos/numerics/NumericsSparseMatrix.h>

#include "siconos/algebra/algebra.hpp"
#include "siconos/algebra/linear_algebra.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/traits/traits.hpp"

namespace siconos::storage::pattern::match {
template <typename T>
concept vec = requires { typename T::vec_t; };

template <typename T>
concept any_mat = !vec<T> && requires { typename T::any_mat_t; };

template <typename T>
concept mat = any_mat<T> && requires { typename T::mat_t; };

template <typename T>
concept diag_mat = any_mat<T> && requires { typename T::diag_mat_t; };

}  // namespace siconos::storage::pattern::match

namespace siconos::algebra {

namespace match = siconos::storage::pattern::match;

static constexpr auto zero_threshold = 1e-30;

struct any_mat {
  using any_mat_t = void;
};

template <typename T>
struct mat : any_mat {
  using mat_t = void;
  static constexpr auto vncols = []() {
    if constexpr (match::matrix<T>) {
      // only eigen
      return T::ColsAtCompileTime;
    }
    else {
      return 1;
    }
  }();
  static constexpr auto vnrows = []() {
    if constexpr (match::matrix<T>) {
      // only eigen
      return T::RowsAtCompileTime;
    }
    else {
      return 1;
    }
  }();

  NumericsMatrix* _m = nullptr;
  bool _inversed = false;         // for diagonal format
  NumericsMatrix* _mt = nullptr;  // transposed is allocated with Numerics
  constexpr mat() {}

  ~mat()
  {
    if (_m) {
      _m = NM_free(_m);
    }
    if (_mt) {
      _m = NM_free(_mt);
    }
  }
};

template <typename T>
struct diag_mat : mat<T> {
  using any_mat_t = typename mat<T>::any_mat_t;
  using diag_mat_t = void;

  static constexpr auto vncols = T::ColsAtCompileTime;
  static constexpr auto vnrows = T::RowsAtCompileTime;
};

template <typename T>
struct vec {
  using vec_t = void;
  static constexpr auto vncols = 1;
  static constexpr auto vnrows = T::RowsAtCompileTime;

  NumericsMatrix* _v = nullptr;

  constexpr vec() {}

  ~vec()
  {
    if (_v) {
      _v = NM_free(_v);
    }
  }
};

template <>
struct vec<double> {
  using vec_t = void;
  static constexpr auto vncols = 1;
  static constexpr auto vnrows = 1;

  NumericsMatrix* _v = nullptr;

  constexpr vec() {}

  ~vec()
  {
    if (_v) {
      _v = NM_free(_v);
    }
  }
};

const auto size0(match::any_mat auto& m) { return m._m->size0 / m.vnrows; };
const auto size0(match::vec auto& v) { return v._v->size0 / v.vnrows; };
const auto size1(match::any_mat auto& m) { return m._m->size1 / m.vncols; };

static_assert(vec<vector<int, 2>>::vncols == 1);

void resize(match::any_mat auto& m, match::indice auto nrows,
            match::indice auto ncols)
{
  if (m._m) m._m = NM_free(m._m);
  if (m._mt) m._mt = NM_free(m._mt);

  m._inversed = false;

  m._m = NM_create(NM_SPARSE, nrows * m.vnrows, ncols * m.vncols);
  NM_triplet_alloc(m._m, 1);
}

void resize(match::vec auto& v, match::indice auto nrows)
{
  if (v._v) v._v = NM_free(v._v);

  //  static_assert(m.vncols == 1);  // only vector of vectors
  // dense vector
  assert(nrows * v.vnrows >= 1);
  assert(v.vncols == 1);
  v._v = NM_create(NM_DENSE, nrows * v.vnrows, v.vncols);
}

void transpose(match::any_mat auto& m)
{
  if (!m._mt) {
    m._mt = NM_transpose(m._m);
  }
}

void setup(match::any_mat auto& m) { resize(m, 1, 1); }

template <typename T>
void set_value(match::any_mat auto& m, match::indice auto i,
               match::indice auto j, const T& value)
{
  if constexpr (match::scalar<T>) {
    NM_zentry(m._m, i * m.vnrows, j * m.vncols, value, zero_threshold);
  }
  // diagonal block
  else if constexpr (match::diagonal_matrix<T>) {
    for (decltype(i) k = 0; k < ncols(T{}); ++k) {
      NM_zentry(m._m, i * m.vnrows + k, j * m.vncols + k, value.diagonal()(k),
                zero_threshold);
    }
  }
  // full block
  else if constexpr (match::matrix<T>) {
    for (decltype(i) k = 0; k < nrows(T{}); ++k) {
      for (decltype(j) l = 0; l < ncols(T{}); ++l) {
        NM_zentry(m._m, i * m.vnrows + k, j * m.vncols + l, value(k, l),
                  zero_threshold);
      }
    }
  }
  else {
    []<bool flag = false>() {
      static_assert(flag, "set_value: cannot insert this value");
    }();
  }
}

template <typename T>
void set_value(match::vec auto& m, match::indice auto i, const T& value)
{
  if constexpr (match::scalar<T>) {
    NM_zentry(m._m, i * m.vnrows, 0, value, zero_threshold);
  }
  // vector block
  else if constexpr (match::vector<T>) {
    for (decltype(i) k = 0; k < nrows(T{}); ++k) {
      NM_zentry(m._v, i * m.vnrows + k, 0, value(k), zero_threshold);
    }
  }
  // compile time error
  else {
    []<bool flag = false>() {
      static_assert(flag, "set_value: cannot insert this value");
    }();
  }
}

template <match::diagonal_matrix A>
void inverse(diag_mat<A>& a)
{
  if (!a._inversed) {
    for (auto i = 0; i < NM_triplet(a._m)->nz; ++i)
      a._m->matrix2->triplet->x[i] = 1.0 / a._m->matrix2->triplet->x[i];
  }
  a._inversed = true;
}

// b += a
template <typename T>
void add(const vec<T>& a, vec<T>& b)
{
  // improve Numerics  cblas_daxpy(size0(a), 1, a._v->matrix0, 1,
  // b._v->matrix0);

  for (size_t i = 0; i < size0(a) * a.vnrows; ++i) {
    b._v->matrix0[i] += a._v->matrix0[i];
  }
}

template <typename T>
void scal(match::scalar auto h, vec<T>& v)
{
  NM_scal(h, v._v);
}

template <typename T>
decltype(auto) get_vector(vec<T>& v, match::indice auto i)
{
  return matrix_view<T>(v._v->matrix0 + i * v.vnrows);
}

//
// c <- a b
// Matrix Matrix
template <typename A, typename B>
void prod(mat<A>& a, mat<B>& b, mat<prod_t<A, B>>& c)
{
  NM_gemm(1, a._m, b._m, 1, c._m);
}

// Matrix Vector
template <typename A, typename B>
void prod(mat<A>& a, vec<B>& b, vec<prod_t<A, B>>& c)
{
  NM_gemv(1, a._m, b._v->matrix0, 1, c._v->matrix0);
}

// c <- a^t b
template <typename A, typename B>
void prodt1(const mat<A>& a, const vec<B>& b, vec<prod_t<trans_t<A>, B>>& c)
{
  assert(a._mt);
  assert(a._mt->size1 == b._v->size0);  // transpose mult
  assert(c._v->size0 == a._mt->size0);
  NM_gemv(1, a._mt, b._v->matrix0, 1, c._v->matrix0);
}
// c <- a b^t
template <typename A, typename B>
void prodt2(const diag_mat<A>& a, const mat<B>& b,
            mat<prod_t<A, trans_t<B>>>& c)
{
  assert(b._mt);
  NM_gemm(1, a._m, b._mt, 1, c._m);
}

// c <- a^-1 b
template <match::diagonal_matrix A, typename B>
void solve(diag_mat<A>& a, vec<B>& b, vec<B>& c)
{
  inverse(a);
  prod(a, b, c);
}

template <match::diagonal_matrix A, match::matrix B>
void solvet(diag_mat<A>& a, mat<B>& b, mat<trans_t<B>>& c)
{
  inverse(a);
  transpose(b);
  prodt2(a, b, c);
}

void display(match::any_mat auto& a) { NM_display(a._m); }

void display(match::vec auto& v) { NM_display(v._v); }

}  // namespace siconos::algebra
