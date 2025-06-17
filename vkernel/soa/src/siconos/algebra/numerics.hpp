#pragma once

#include "CSparseMatrix.h"  // for CSparseMatrix
#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"
#include "siconos/algebra/algebra.hpp"
#include "siconos/algebra/eigen.hpp"
#include "siconos/algebra/linear_algebra.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/traits/traits.hpp"

namespace siconos::storage::pattern::match {
template <typename T>
concept vec = requires { typename std::decay_t<T>::vec_t; };

template <typename T>
concept any_mat =
    !vec<T> && requires { typename std::decay_t<T>::any_mat_t; };

template <typename T>
concept mat = any_mat<T> && requires { typename std::decay_t<T>::mat_t; };

template <typename T>
concept diag_mat =
    any_mat<T> && requires { typename std::decay_t<T>::diag_mat_t; };

}  // namespace siconos::storage::pattern::match

namespace siconos::algebra {

namespace match = siconos::storage::pattern::match;

using indice_t = std::size_t;

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

  indice_t _offsets[2] = {0, 0};

  NumericsMatrix* _m = nullptr;
  NumericsMatrix* _mt = nullptr;  // transposed is allocated with Numerics

  bool _inversed = false;  // for diagonal format
  bool _view = false;      // for destruction

  constexpr mat() {}

  ~mat()
  {
    if (_m) {
      if (_view) {
        _m = nullptr;
      }
      else {
        _m = NM_free(_m);
      }
    }
    if (_mt) {
      if (_view) {
        _mt = nullptr;
      }
      else {
        _mt = NM_free(_mt);
      }
    }
    _inversed = false;
    _view = false;
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

  bool _view = false;

  constexpr vec() {}

  ~vec()
  {
    if (_v) {
      if (_view) {
        _v = nullptr;
      }
      else {
        _v = NM_free(_v);
      }
    }

    _view = false;
    _offset = 0;
  }

  indice_t _offset = 0;
};

// template <>
// struct vec<double> {
//   using vec_t = void;
//   static constexpr auto vncols = 1;
//   static constexpr auto vnrows = 1;

//   NumericsMatrix* _v = nullptr;
//   _view = false;

//   constexpr vec() {}

//   ~vec()
//   {
//     if (_v) {
//       _v = NM_free(_v);
//     }
//   }
// };

template <typename T>
mat<T> mat_view(match::any_mat auto& m, indice_t row_offset,
                indice_t col_offset)
{
  mat<T> vm;
  /* pointers copy */
  vm._m = m._m;
  vm._mt = m._mt;
  vm._inversed = m._inversed;
  vm._offsets[0] = row_offset;
  vm._offsets[1] = col_offset;
  vm._view = true;

  return vm; /* RVO */
}

template <typename T>
mat<T> mat_view(T, match::any_mat auto& m, indice_t row_offset,
                indice_t col_offset)
{
  return mat_view<T>(m, row_offset, col_offset);
}

template <typename T>
vec<T> vec_view(match::vec auto& v, indice_t offset)
{
  vec<T> vv;
  /* pointers copy */
  vv._v = v._v;
  vv._offset = offset;
  vv._view = true;
  return vv; /* RVO */
}

template <typename T>
decltype(auto) vec_view(T, match::vec auto& v, indice_t offset)
{
  return vec_view<T>(v, offset);
}

const auto size0(match::any_mat auto& m) { return m._m->size0 / m.vnrows; };
const auto size0(match::vec auto& v) { return v._v->size0 / v.vnrows; };
const auto size1(match::any_mat auto& m) { return m._m->size1 / m.vncols; };

const auto raw_size0(match::any_mat auto& m) { return m._m->size0; };
const auto raw_size0(match::vec auto& v) { return v._v->size0; };
const auto raw_size1(match::any_mat auto& m) { return m._m->size1; };

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

void insert(match::any_mat auto& am, match::any_mat auto& bm,
            match::indice auto offset0, match::indice auto offset1)
{
  NM_insert(am, bm, offset0, offset1);
}

void insert(match::vec auto& va, match::vec auto& vb,
            match::indice auto offset0)
{
  NM_insert(va, vb, offset0, 1);
}

void transpose(match::any_mat auto& m)
{
  if (!m._mt) {
    m._mt = NM_transpose(m._m);
  }
}

void setup(match::any_mat auto& m) { resize(m, 1, 1); }

template <match::any_mat M, typename T>
void set_value(M&& m, match::indice auto i, match::indice auto j,
               const T& value)
{
  if constexpr (match::scalar<T>) {
    NM_zentry(m._m, i * m.vnrows + m._offsets[0],
              j * m.vncols + m._offsets[1], value, zero_threshold);
  }
  // diagonal block
  else if constexpr (match::diagonal_matrix<T>) {
    for (decltype(i) k = 0; k < ncols(T{}); ++k) {
      NM_zentry(m._m, i * m.vnrows + k + m._offsets[0],
                j * m.vncols + k + m._offsets[1], value.diagonal()(k),
                zero_threshold);
    }
  }
  // full block
  else if constexpr (match::matrix<T>) {
    for (decltype(i) k = 0; k < nrows(T{}); ++k) {
      for (decltype(j) l = 0; l < ncols(T{}); ++l) {
        NM_zentry(m._m, i * m.vnrows + k + m._offsets[0],
                  j * m.vncols + l + m._offsets[1], value(k, l),
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
void set_value(match::vec auto&& m, match::indice auto i, const T& value)
{
  if constexpr (match::scalar<T>) {
    NM_zentry(m._m, i * m.vnrows + m._offset, 0, value, zero_threshold);
  }
  // vector block
  else if constexpr (match::vector<T>) {
    for (decltype(i) k = 0; k < nrows(T{}); ++k) {
      NM_zentry(m._v, i * m.vnrows + k + m._offset, 0, value(k),
                zero_threshold);
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

// b <- a
template <typename V>
void copy(const vec<V>& a, vec<V>& b)
{
  for (auto i = 0; i < size0(a) * a.vnrows; ++i) {
    b._v->matrix0[i] = a._v->matrix0[i];
  }
}

// v <- h*v
template <typename T>
void scal(match::scalar auto h, vec<T>& v)
{
  NM_scal(h, v._v);
}

template <typename T>
decltype(auto) get_vector(vec<T>& v, match::indice auto i)
{
  return matrix_view<T>(v._v->matrix0 + i * v.vnrows + v._offset);
}

template <typename T>
decltype(auto) get_vector(const vec<T>& v, match::indice auto i)
{
  return matrix_view<T>(v._v->matrix0 + i * v.vnrows + v._offset);
}


template <typename T>
decltype(auto) get_vector(vec<T>&& v, match::indice auto i)
{
  return get_vector(v, i);
}

template <typename T>
decltype(auto) get_vector(vec<T>& v, match::indice auto i,
                          match::indice auto vector_size)
{
  return matrix_view<T>(v._v->matrix0 + i * v.vnrows + v._offset,
                        vector_size);
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

template <match::any_matrix A, typename B>
void solve(const mat<A>& a, const vec<B>& b, vec<B>& c)
{
  copy(b, c);
  if (NM_gesv_expert(a._m, c._v->matrix0, 1) != 0) {
    throw std::runtime_error(
        "NumericsMatrix solve failed in NM_gesv_expert.");
  }
}

template <match::any_matrix A, typename B>
void solve_in_place(const mat<A>& a, const vec<B>& b)
{
  if (NM_gesv_expert(a._m, b._v->matrix0, 1) != 0) {
    throw std::runtime_error(
        "NumericsMatrix solve_in_place failed in NM_gesv_expert.");
  }
}

template <match::diagonal_matrix A, match::matrix B>
void solvet(diag_mat<A>& a, mat<B>& b, mat<trans_t<B>>& c)
{
  inverse(a);
  transpose(b);
  prodt2(a, b, c);
}

template <match::any_matrix H, match::any_matrix M, typename W>
void compute_kkt_matrix(mat<H>& h_matrix, mat<M>& mass_matrix,
                        mat<W>& w_matrix)
{
  // Ensure mass_matrix is square
  assert(mass_matrix._m->size0 == mass_matrix._m->size1);

  // Transpose H if not already done
  transpose(h_matrix);

  // Resize w_matrix (output): w = H M^{-1} H^T
  // H: (n x d), M: (d x d), H^T: (d x n), W: (n x n)
  int n = h_matrix._m->size0;
  int d = h_matrix._m->size1;

  if (!w_matrix._m || w_matrix._m->size0 != n || w_matrix._m->size1 != n) {
    if (w_matrix._m) w_matrix._m = NM_free(w_matrix._m);
    w_matrix._m = NM_create(NM_SPARSE, n, n);
    NM_triplet_alloc(w_matrix._m, 1);
  }
  else if (w_matrix._m->storageType != NM_SPARSE) {
    w_matrix._m = NM_free(w_matrix._m);
    w_matrix._m = NM_create(NM_SPARSE, n, n);
    NM_triplet_alloc(w_matrix._m, 1);
  }

  // temp NumericsMatrix for storing M^{-1} H^T
  NumericsMatrix* Minv_ht = NM_create(NM_SPARSE, d, n);
  NM_triplet_alloc(Minv_ht, 1);

  // Dense buffer for rhs (right hand side), reused for each solve
  std::vector<double> rhs_vec(d);

  // Support for Siconos NumericsSparseMatrix either triplet or compressed row
  NumericsSparseMatrix* H_sparse = h_matrix._m->matrix2;
  bool is_sparse = (h_matrix._m->storageType == NM_SPARSE && H_sparse);

  for (auto k = 0; k < n; ++k) {
    std::fill(rhs_vec.begin(), rhs_vec.end(), 0.0);

    if (is_sparse) {
      if (H_sparse->triplet) {
        // Triplet format: scan all triplets to collect row k
        for (auto l = 0; l < H_sparse->triplet->nz; ++l) {
          auto row_idx = H_sparse->triplet->i[l];
          auto col_idx = H_sparse->triplet->p[l];
          if (row_idx == k && col_idx < d) {
            rhs_vec[col_idx] = H_sparse->triplet->x[l];
          }
        }
      }
      else if (H_sparse->csr) {
        // Compressed Sparse Row format: iterate shared row structure
        CS_INT* row_ptr = H_sparse->csr->p;
        CS_INT* col_ind = H_sparse->csr->i;
        double* values = H_sparse->csr->x;
        for (auto idx = row_ptr[k]; idx < row_ptr[k + 1]; ++idx) {
          auto col_idx = col_ind[idx];
          if (col_idx < d) {
            rhs_vec[col_idx] = values[idx];
          }
        }
      }
      else {
        std::runtime_error("solve_kkt: unknown sparse format");
      }
    }
    else {
      // fallback to dense (should not happen if h_matrix is actually sparse,
      // but for completeness)
      if (h_matrix._mt && h_matrix._mt->matrix0) {
        for (int j = 0; j < d; ++j) {
          rhs_vec[j] = h_matrix._mt->matrix0[k * d + j];
        }
      }
      else {
        for (int j = 0; j < d; ++j) {
          rhs_vec[j] = 0.0;
        }
      }
    }

    if (NM_gesv_expert(mass_matrix._m, rhs_vec.data(), 0) != 0) {
      NM_free(Minv_ht);
      throw std::runtime_error(
          "NumericsMatrix solve failed in solve_kkt: M^{-1} * col(H^T)");
    }

    for (int i = 0; i < d; ++i) {
      double val = rhs_vec[i];
      if (std::abs(val) > zero_threshold) {
        NM_zentry(Minv_ht, i, k, val, zero_threshold);
      }
    }
  }

  NM_gemm(1.0, h_matrix._m, Minv_ht, 0.0, w_matrix._m);

  NM_free(Minv_ht);
}

void display(match::any_mat auto& a) { NM_display(a._m); }

void display(match::vec auto& v) { NM_display(v._v); }

}  // namespace siconos::algebra
