#pragma once

#include <algorithm>
#include <cstdio>
#include <map>
#include <memory>
#include <span>

#include "CSparseMatrix.h"  // for CSparseMatrix
#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"
#include "siconos/algebra/algebra.hpp"
#include "siconos/algebra/eigen.hpp"
#include "siconos/algebra/linear_algebra.hpp"
#include "siconos/algebra/sparse_factor.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/traits/traits.hpp"

namespace siconos::storage::pattern::match {
template <typename T>
concept vec = requires { typename std::decay_t<T>::vec_t; };

template <typename T>
concept any_mat = !vec<T> && requires { typename std::decay_t<T>::any_mat_t; };

template <typename T>
concept mat = any_mat<T> && requires { typename std::decay_t<T>::mat_t; };

template <typename T>
concept diag_mat = any_mat<T> && requires { typename std::decay_t<T>::diag_mat_t; };

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
    if constexpr (match::fixed_size_matrix<T>) {
      // only eigen
      return T::ColsAtCompileTime;
    } else {
      // runtime
      return 1;
    }
  }();
  static constexpr auto vnrows = []() {
    if constexpr (match::fixed_size_matrix<T>) {
      // only eigen
      return T::RowsAtCompileTime;
    } else {
      // runtime
      return 1;
    }
  }();

  indice_t _offsets[2] = {0, 0};

  NumericsMatrix* _m = nullptr;
  NumericsMatrix* _mt = nullptr;  // transposed is allocated with Numerics

  bool _view = false;  // for destruction

  constexpr mat() {}

  virtual ~mat() {
    if (_m) {
      if (_view) {
        _m = nullptr;
      } else {
        _m = NM_free(_m);
      }
    }
    if (_mt) {
      _mt = NM_free(_mt);
    }
    _view = false;
  }
};

template <typename T>
struct diag_mat : mat<T> {
  using any_mat_t = typename mat<T>::any_mat_t;
  using diag_mat_t = void;
};

template <typename T>
diag_mat<T> cast_to_diag_mat(mat<T>& m) {
  diag_mat<T> result;
  result._m = m._m;
  result._offsets[0] = m._offsets[0];
  result._offsets[1] = m._offsets[1];
  result._view = m._view;
  return result;
}

template <typename T>
struct vec {
  using vec_t = void;

  static constexpr auto vncols = 1;

  static constexpr auto vnrows = []() {
    if constexpr (match::fixed_size_matrix<T>) {
      // only eigen
      return T::RowsAtCompileTime;
    } else {
      // runtime
      return 1;
    }
  }();

  NumericsMatrix* _v = nullptr;

  bool _view = false;

  constexpr vec() {}

  ~vec() {
    if (_v) {
      if (_view) {
        _v = nullptr;
      } else {
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
diag_mat<T> mat_view(match::diag_mat auto& m, indice_t row_offset, indice_t col_offset) {
  diag_mat<T> vm;
  /* pointers copy */
  vm._m = m._m;
  vm._mt = m._mt;
  vm._offsets[0] = row_offset;
  vm._offsets[1] = col_offset;
  vm._view = true;

  return vm; /* RVO */
}

template <typename T>
diag_mat<T> mat_view(T, match::diag_mat auto& m, indice_t row_offset, indice_t col_offset) {
  return mat_view<T>(m, row_offset, col_offset);
}

template <typename T>
mat<T> mat_view(match::any_mat auto& m, indice_t row_offset, indice_t col_offset) {
  mat<T> vm;
  /* pointers copy */
  vm._m = m._m;
  vm._mt = m._mt;
  vm._offsets[0] = row_offset;
  vm._offsets[1] = col_offset;
  vm._view = true;

  return vm; /* RVO */
}

template <typename T>
mat<T> mat_view(T, match::any_mat auto& m, indice_t row_offset, indice_t col_offset) {
  return mat_view<T>(m, row_offset, col_offset);
}

template <typename T>
vec<T> vec_view(match::vec auto& v, indice_t offset) {
  vec<T> vv;
  /* pointers copy */
  vv._v = v._v;
  vv._offset = offset;
  vv._view = true;
  return vv; /* RVO */
}

template <typename T>
decltype(auto) vec_view(T, match::vec auto& v, indice_t offset) {
  return vec_view<T>(v, offset);
}

const auto size0(match::any_mat auto& m) { return m._m->size0 / m.vnrows; };
const auto size0(match::vec auto& v) { return v._v->size0 / v.vnrows; };
const auto size1(match::any_mat auto& m) { return m._m->size1 / m.vncols; };

const auto raw_size0(match::any_mat auto& m) { return m._m->size0; };
const auto raw_size0(match::vec auto& v) { return v._v->size0; };
const auto raw_size1(match::any_mat auto& m) { return m._m->size1; };

void resize(match::any_mat auto& m, match::indice auto nrows, match::indice auto ncols) {
  if (m._m) m._m = NM_free(m._m);
  if (m._mt) m._mt = NM_free(m._mt);

  m._m = NM_create(NM_SPARSE, nrows * m.vnrows, ncols * m.vncols);

  NM_triplet_alloc(m._m, 1);
}

void resize(match::vec auto& v, match::indice auto nrows) {
  if (v._v) v._v = NM_free(v._v);

  //  static_assert(m.vncols == 1);  // only vector of vectors
  // dense vector
  assert(nrows * v.vnrows >= 1);
  assert(v.vncols == 1);
  v._v = NM_create(NM_DENSE, nrows * v.vnrows, v.vncols);
}

auto nnz(match::any_mat auto& m) { return NM_triplet(m._m)->nz; }

template <typename T>
vec<T> copy(vec<T>& v) {
  vec<T> new_vec;

  auto* new_v = NM_new();
  NM_copy(v._v, new_v);

  // new_vec destructor will free new_v
  new_vec._v = new_v;

  return new_vec;
}

void insert(match::any_mat auto& am, match::any_mat auto& bm, match::indice auto offset0,
            match::indice auto offset1) {
  NM_insert(am, bm, offset0, offset1);
}

void insert(match::vec auto& va, match::vec auto& vb, match::indice auto offset0) {
  NM_insert(va, vb, offset0, 1);
}

void transpose(match::any_mat auto& m) {
  if (!m._mt) {
    m._mt = NM_transpose(m._m);
  }
}

void setup(match::any_mat auto& m) { resize(m, 1, 1); }

template <match::any_mat M, typename T>
void set_value(M&& m, match::indice auto i, match::indice auto j, const T& value) {
  assert(i >= 0 && j >= 0);
  assert(i * m.vnrows + m._offsets[0] < m._m->size0);
  assert(j * m.vncols + m._offsets[1] < m._m->size1);

  if constexpr (match::scalar<T>) {
    NM_zentry(m._m, i * m.vnrows + m._offsets[0], j * m.vncols + m._offsets[1], value,
              zero_threshold);
  }
  // diagonal block
  else if constexpr (match::diagonal_matrix<T>) {
    for (decltype(i) k = 0; k < ncols(value); ++k) {
      NM_zentry(m._m, i * m.vnrows + k + m._offsets[0], j * m.vncols + k + m._offsets[1],
                value.diagonal()(k), zero_threshold);
    }
  }
  // full block
  else if constexpr (match::matrix<T>) {
    for (decltype(i) k = 0; k < nrows(value); ++k) {
      for (decltype(j) l = 0; l < ncols(value); ++l) {
        NM_zentry(m._m, i * m.vnrows + k + m._offsets[0], j * m.vncols + l + m._offsets[1],
                  value(k, l), zero_threshold);
      }
    }
  } else if constexpr (match::sparse_matrix<T>) {
    for (typename T::Index k = 0; k < value.outerSize(); ++k) {
      // /!\ specific to eigen
      for (typename T::InnerIterator it(value, k); it; ++it) {
        NM_zentry(m._m, i * m.vnrows + it.row() + m._offsets[0],
                  j * m.vncols + it.col() + m._offsets[1], it.value(), zero_threshold);
      }
    }
  } else {
    []<bool flag = false>() { static_assert(flag, "set_value: cannot insert this value"); }();
  }
}

template <typename T>
void set_value(match::vec auto&& m, match::indice auto i, const T& value) {
  assert(i >= 0);
  assert(i * m.vnrows + m._offset < m._v->size0);

  if constexpr (match::scalar<T>) {
    NM_zentry(m._v, i * m.vnrows + m._offset, 0, value, zero_threshold);
  }
  // vector block
  else if constexpr (match::vector<T>) {
    for (decltype(i) k = 0; k < nrows(value); ++k) {
      NM_zentry(m._v, i * m.vnrows + k + m._offset, 0, value(k), zero_threshold);
    }
  }
  // compile time error
  else {
    []<bool flag = false>() { static_assert(flag, "set_value: cannot insert this value"); }();
  }
}

template <match::diagonal_or_1x1 A>
void invert_diagonal_matrix(diag_mat<A>& a) {
  if (!NM_internalData(a._m)->isInversed) {
    for (auto i = 0; i < NM_triplet(a._m)->nz; ++i)
      a._m->matrix2->triplet->x[i] = 1.0 / a._m->matrix2->triplet->x[i];
  }
  NM_internalData(a._m)->isInversed = true;
}

// b += a
template <typename T>
void add(const vec<T>& a, vec<T>& b) {
  // improve Numerics  cblas_daxpy(size0(a), 1, a._v->matrix0, 1,
  // b._v->matrix0);

  for (size_t i = 0; i < size0(a) * a.vnrows; ++i) {
    b._v->matrix0[i] += a._v->matrix0[i];
  }
}

// alpha * A + beta * B
template <match::scalar scalar, match::any_matrix M>
mat<M> add(scalar alpha, const mat<M>& A, scalar beta, const mat<M>& B) {
  NumericsMatrix* nm;

  if (nnz(B) > 0) {
    nm = NM_add(alpha, A._m, beta, B._m);
  } else {
    nm = NM_new();
    NM_copy(A._m, nm);
  }
  mat<M> rm;
  rm._m = nm;
  return rm;
}

// c = alpha * A + beta * B
template <match::scalar scalar, match::any_matrix M>
void add(scalar alpha, const mat<M>& A, scalar beta, const mat<M>& B, mat<M>& C) {
  if (nnz(B) > 0) {
    C._m = NM_add(alpha, A._m, beta, B._m);
  } else {
    NumericsMatrix* nm = NM_new();
    NM_copy(A._m, nm);
    C._m = nm;
  }
  C._mt = nullptr;
  C._view = false;
}

// b <- a
template <typename V>
void copy(const vec<V>& a, vec<V>& b) {
  for (auto i = 0; i < size0(a) * a.vnrows; ++i) {
    b._v->matrix0[i] = a._v->matrix0[i];
  }
}

// v <- h*v
template <typename T>
void scal(match::scalar auto h, vec<T>& v) {
  NM_scal(h, v._v);
}

template <typename T>
decltype(auto) get_vector(vec<T>& v, match::indice auto i) {
  return matrix_view<T>(v._v->matrix0 + i * v.vnrows + v._offset);
}

template <typename T>
decltype(auto) get_vector(const vec<T>& v, match::indice auto i) {
  return matrix_view<T>(v._v->matrix0 + i * v.vnrows + v._offset);
}

template <typename T>
decltype(auto) get_vector(vec<T>&& v, match::indice auto i) {
  return get_vector(v, i);
}

template <typename T>
decltype(auto) get_vector(vec<T>& v, match::indice auto i, match::indice auto vector_size) {
  return matrix_view<T>(v._v->matrix0 + i * v.vnrows + v._offset, vector_size);
}
//
// c <- a b
// Matrix Matrix
template <typename A, typename B>
void prod(const mat<A>& a, const mat<B>& b, mat<prod_t<A, B>>& c) {
  NM_gemm(1, a._m, b._m, 1, c._m);
}

// Matrix Vector
template <typename A, typename B>
void prod(const mat<A>& a, const vec<B>& b, vec<prod_t<A, B>>&& c) {
  NM_gemv(1, a._m, b._v->matrix0, 1, c._v->matrix0);
}

// c <- a^t b
template <typename A, typename B>
void prodt1(const mat<A>& a, const vec<B>& b, vec<prod_t<trans_t<A>, B>>& c) {
  assert(a._mt);
  assert(a._mt->size1 == b._v->size0);  // transpose mult
  assert(c._v->size0 == a._mt->size0);
  NM_gemv(1, a._mt, b._v->matrix0, 1, c._v->matrix0);
}
// c <- a b^t
template <typename A, typename B>
void prodt2(const diag_mat<A>& a, const mat<B>& b, mat<prod_t<A, trans_t<B>>>& c) {
  assert(b._mt);
  NM_gemm(1, a._m, b._mt, 1, c._m);
}

// c <- a^-1 b
template <match::diagonal_matrix A, typename B>
void solve_linear_system(diag_mat<A>& a, vec<B>& b, vec<B>& c) {
  invert_diagonal_matrix(a);
  prod(a, b, c);
}

template <match::any_matrix A, typename B>
void solve_linear_system(const mat<A>& a, const vec<B>& b, vec<B>& c) {
  if (!NM_internalData(a._m)->isInversed) {
    copy(b, c);
    if (NM_gesv_expert(a._m, c._v->matrix0, 1) != 0) {
      throw std::runtime_error("NumericsMatrix solve failed in NM_gesv_expert.");
    }
  } else {
    // this case occurs if a solve has been done with a diagonal view.
    prod(a, b, vec_view<prod_t<A, B>>(c, 0));
  }
}

template <match::any_matrix A, typename B>
void solve_in_place(const mat<A>& a, const vec<B>& b) {
  if (NM_gesv_expert(a._m, b._v->matrix0, 1) != 0) {
    throw std::runtime_error("NumericsMatrix solve_in_place failed in NM_gesv_expert.");
  }
}

template <match::diagonal_or_1x1 A, match::matrix B>
void solve_linear_system_with_transpose(diag_mat<A>& a, mat<B>& b, mat<trans_t<B>>& c) {
  invert_diagonal_matrix(a);
  transpose(b);
  prodt2(a, b, c);
}

template <match::diagonal_or_1x1 A, match::matrix B>
void solve_linear_system_with_transpose(diag_mat<A>&& a, mat<B>& b, mat<trans_t<B>>& c) {
  solve_linear_system_with_transpose(a, b, c);
}

// Diagonal blocks of a block diagonal sparse matrix, found as the connected
// components of its pattern.  For the assembled iteration matrix of the
// vkernel, rigid systems with a diagonal mass matrix give 1x1 blocks and
// each involved FEM system gives one large block.
struct block_diagonal_structure {
  std::vector<CS_INT> comp;      // dof -> component
  std::vector<CS_INT> comp_ptr;  // component -> range in dofs
  std::vector<CS_INT> dofs;      // dofs grouped by component, increasing
  std::vector<CS_INT> local;     // dof -> local index in its component

  explicit block_diagonal_structure(const CSparseMatrix* a) {
    CS_INT d = a->n;
    std::vector<CS_INT> parent(d);
    for (CS_INT j = 0; j < d; ++j) parent[j] = j;
    auto find = [&parent](CS_INT x) {
      while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
      }
      return x;
    };
    for (CS_INT j = 0; j < d; ++j) {
      for (CS_INT p = a->p[j]; p < a->p[j + 1]; ++p) {
        auto ri = find(a->i[p]);
        auto rj = find(j);
        if (ri != rj) parent[std::max(ri, rj)] = std::min(ri, rj);
      }
    }
    // components numbered by their smallest dof
    comp.assign(d, -1);
    std::vector<CS_INT> root_comp(d, -1);
    std::vector<CS_INT> count;
    for (CS_INT j = 0; j < d; ++j) {
      auto r = find(j);
      if (root_comp[r] < 0) {
        root_comp[r] = static_cast<CS_INT>(count.size());
        count.push_back(0);
      }
      comp[j] = root_comp[r];
      count[comp[j]]++;
    }
    comp_ptr.assign(count.size() + 1, 0);
    for (std::size_t c = 0; c < count.size(); ++c) comp_ptr[c + 1] = comp_ptr[c] + count[c];
    dofs.resize(d);
    local.resize(d);
    std::vector<CS_INT> fill(comp_ptr.begin(), comp_ptr.end() - 1);
    for (CS_INT j = 0; j < d; ++j) {
      auto pos = fill[comp[j]]++;
      dofs[pos] = j;
      local[j] = pos - comp_ptr[comp[j]];
    }
  }

  CS_INT number_of_components() const { return static_cast<CS_INT>(comp_ptr.size()) - 1; }
  CS_INT size(CS_INT c) const { return comp_ptr[c + 1] - comp_ptr[c]; }
  CS_INT dof(CS_INT c, CS_INT l) const { return dofs[comp_ptr[c] + l]; }

  // value of a 1x1 block
  static double diagonal(const CSparseMatrix* a, CS_INT j) {
    double v = 0.;
    for (CS_INT p = a->p[j]; p < a->p[j + 1]; ++p) v += a->x[p];
    return v;
  }

  Eigen::MatrixXd dense_block(const CSparseMatrix* a, CS_INT c) const {
    Eigen::MatrixXd m = Eigen::MatrixXd::Zero(size(c), size(c));
    for (CS_INT l = 0; l < size(c); ++l) {
      auto j = dof(c, l);
      for (CS_INT p = a->p[j]; p < a->p[j + 1]; ++p) m(local[a->i[p]], l) += a->x[p];
    }
    return m;
  }

  sparse_factor::matrix_t sparse_block(const CSparseMatrix* a, CS_INT c) const {
    std::vector<Eigen::Triplet<double, sparse_factor::index_t>> triplets;
    for (CS_INT l = 0; l < size(c); ++l) {
      auto j = dof(c, l);
      for (CS_INT p = a->p[j]; p < a->p[j + 1]; ++p) {
        triplets.emplace_back(static_cast<sparse_factor::index_t>(local[a->i[p]]),
                              static_cast<sparse_factor::index_t>(l), a->x[p]);
      }
    }
    sparse_factor::matrix_t m(size(c), size(c));
    m.setFromTriplets(triplets.begin(), triplets.end());
    m.makeCompressed();
    return m;
  }
};

// blocks up to this size are solved densely at each call, larger ones with
// a sparse factorization.
static constexpr CS_INT small_block_size = 32;

// The diagonal blocks of a block diagonal matrix A, and how to solve with
// each of them.  Blocks with a known factorization (FEM systems, see
// fem_operator) are given by their offset in A.  The others are found as
// connected components of the pattern of A: 1x1 blocks (rigid systems with
// a diagonal mass matrix) are a division, small blocks a dense LU, large
// blocks a sparse factorization local to this call.
class block_diagonal_solver {
 public:
  block_diagonal_solver(const CSparseMatrix* a, std::span<const factored_block> factored)
      : _a(a), _blocks(a), _known(a->n, -1) {
    for (const auto& fb : factored) {
      if (check(fb)) {
        auto k = static_cast<CS_INT>(_factored.size());
        _factored.push_back(fb);
        for (CS_INT j = fb.offset; j < fb.offset + fb.size; ++j) _known[j] = k;
      } else {
        std::fprintf(stderr,
                     "block_diagonal_solver: factored block at offset %td (size %td) does "
                     "not match the assembled matrix, it is refactorized\n",
                     fb.offset, fb.size);
      }
    }
  }

  // group of a dof: a factored block (0 .. nf-1) or a component (nf + c)
  CS_INT group(CS_INT j) const {
    return _known[j] >= 0 ? _known[j] : number_of_factored() + _blocks.comp[j];
  }
  CS_INT number_of_factored() const { return static_cast<CS_INT>(_factored.size()); }
  CS_INT number_of_groups() const {
    return number_of_factored() + _blocks.number_of_components();
  }

  CS_INT size(CS_INT g) const {
    return g < number_of_factored() ? _factored[g].size
                                    : _blocks.size(g - number_of_factored());
  }
  CS_INT local(CS_INT j) const {
    return _known[j] >= 0 ? j - _factored[_known[j]].offset : _blocks.local[j];
  }
  CS_INT dof(CS_INT g, CS_INT l) const {
    return g < number_of_factored() ? _factored[g].offset + l
                                    : _blocks.dof(g - number_of_factored(), l);
  }
  // a component whose dofs belong to a factored block: the blocks are
  // checked to be diagonal blocks of A, so a component is either fully
  // inside one of them or outside all of them
  bool covered(CS_INT g) const { return g >= number_of_factored() && _known[dof(g, 0)] >= 0; }
  // a 1x1 component that is not part of a factored block
  bool diagonal(CS_INT g) const {
    return g >= number_of_factored() && size(g) == 1 && !covered(g);
  }
  double diagonal_value(CS_INT g) const {
    return block_diagonal_structure::diagonal(_a, dof(g, 0));
  }

  // sparse factorization of a group, nullptr for a small (dense) one
  sparse_factor* factor(CS_INT g) {
    if (g < number_of_factored()) return _factored[g].factor;
    if (size(g) <= small_block_size) return nullptr;
    auto& f = _local_factors[g];
    if (!f) {
      f = std::make_unique<sparse_factor>(_blocks.sparse_block(_a, g - number_of_factored()));
    }
    return f.get();
  }

  Eigen::MatrixXd solve_dense(CS_INT g, const Eigen::MatrixXd& rhs) const {
    return _blocks.dense_block(_a, g - number_of_factored()).partialPivLu().solve(rhs);
  }

 private:
  // the factored block must be a diagonal block of A (no coupling outside
  // its range) with the same diagonal: guards against a wrong offset.
  bool check(const factored_block& fb) const {
    if (fb.offset < 0 || fb.offset + fb.size > _a->n || !fb.factor ||
        fb.factor->matrix().rows() != fb.size) {
      return false;
    }
    const auto& w = fb.factor->matrix();
    for (CS_INT l = 0; l < fb.size; ++l) {
      auto j = fb.offset + l;
      double ajj = 0.;
      for (CS_INT p = _a->p[j]; p < _a->p[j + 1]; ++p) {
        auto i = _a->i[p];
        if (i < fb.offset || i >= fb.offset + fb.size) return false;
        if (i == j) ajj += _a->x[p];
      }
      double wjj = w.coeff(l, l);
      if (std::abs(ajj - wjj) > 1e-10 * std::max(std::abs(wjj), 1.)) return false;
    }
    return true;
  }

  const CSparseMatrix* _a;
  block_diagonal_structure _blocks;
  std::vector<CS_INT> _known;  // dof -> factored block or -1
  std::vector<factored_block> _factored;
  std::map<CS_INT, std::unique_ptr<sparse_factor>> _local_factors;
};

// w <- H A^{-1} H^T for a block diagonal A (the assembled iteration matrix
// M + h^2 theta^2 K).
//
// Y = A^{-1} H^T is built block by block, then W = H Y.  1x1 blocks (rigid
// systems) are a division.  Large blocks (FEM systems) only involve the rows
// of H that touch them, and use the compliance A^{-1} restricted to the
// dofs present in H, kept with the factorization: once every contact dof
// has been seen, no linear solve is done here.  Y rows are only needed at
// the dofs present in H, since W = H Y never reads the others.
template <match::any_matrix H, match::any_matrix M, typename W>
void compute_kkt_matrix(mat<H>& h_matrix, mat<M>& iteration_matrix, mat<W>& w_matrix,
                        std::span<const factored_block> factored = {}) {
  assert(iteration_matrix._m->size0 == iteration_matrix._m->size1);

  // Transpose H if not already done
  transpose(h_matrix);

  // H: (n x d), A: (d x d), W: (n x n)
  int n = h_matrix._m->size0;
  int d = h_matrix._m->size1;
  assert(iteration_matrix._m->size0 == d);

  if (!w_matrix._m || w_matrix._m->size0 != n || w_matrix._m->size1 != n ||
      w_matrix._m->storageType != NM_SPARSE) {
    if (w_matrix._m) w_matrix._m = NM_free(w_matrix._m);
    w_matrix._m = NM_create(NM_SPARSE, n, n);
    NM_triplet_alloc(w_matrix._m, 1);
  }

  if (h_matrix._m->storageType != NM_SPARSE) {
    throw std::runtime_error("compute_kkt_matrix: h_matrix must be sparse");
  }

  const CSparseMatrix* a = NM_csc(iteration_matrix._m);
  const CSparseMatrix* hc = NM_csc(h_matrix._m);  // duplicates summed

  block_diagonal_solver blocks(a, factored);

  // Y = A^{-1} H^T (d x n)
  NumericsMatrix* minv_ht = NM_create(NM_SPARSE, d, n);
  NM_triplet_alloc(minv_ht, 1);

  struct h_entry {
    CS_INT row;
    CS_INT local;
    double value;
  };
  std::map<CS_INT, std::vector<h_entry>> block_entries;

  for (CS_INT j = 0; j < d; ++j) {
    if (hc->p[j] == hc->p[j + 1]) continue;
    auto g = blocks.group(j);
    if (blocks.diagonal(g)) {
      double ajj = blocks.diagonal_value(g);
      if (ajj == 0.) {
        NM_free(minv_ht);
        throw std::runtime_error("compute_kkt_matrix: zero diagonal block");
      }
      for (CS_INT p = hc->p[j]; p < hc->p[j + 1]; ++p) {
        NM_zentry(minv_ht, j, hc->i[p], hc->x[p] / ajj, zero_threshold);
      }
    } else {
      for (CS_INT p = hc->p[j]; p < hc->p[j + 1]; ++p) {
        block_entries[g].push_back({hc->i[p], blocks.local(j), hc->x[p]});
      }
    }
  }

  for (auto& [g, entries] : block_entries) {
    // rows of H touching block g, and the block dofs present in H
    std::vector<CS_INT> rows;
    std::vector<sparse_factor::index_t> locals;
    for (auto& e : entries) {
      rows.push_back(e.row);
      locals.push_back(static_cast<sparse_factor::index_t>(e.local));
    }
    std::ranges::sort(rows);
    rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
    std::ranges::sort(locals);
    locals.erase(std::unique(locals.begin(), locals.end()), locals.end());

    auto row_slot = [&rows](CS_INT r) {
      return std::lower_bound(rows.begin(), rows.end(), r) - rows.begin();
    };

    auto bsize = blocks.size(g);
    sparse_factor* factor = blocks.factor(g);

    if (factor && factor->symmetric() && factor->can_track(locals)) {
      factor->ensure_compliance(locals);

      std::vector<sparse_factor::index_t> local_slots(locals.size());
      for (std::size_t l = 0; l < locals.size(); ++l) local_slots[l] = factor->slot(locals[l]);
      for (auto& e : entries)
        e.local = factor->slot(static_cast<sparse_factor::index_t>(e.local));

      // group entries by row
      std::ranges::stable_sort(entries, {}, &h_entry::row);
      for (std::size_t b = 0; b < entries.size();) {
        auto e = b;
        while (e < entries.size() && entries[e].row == entries[b].row) ++e;
        auto row = entries[b].row;
        for (std::size_t l = 0; l < locals.size(); ++l) {
          double y = 0.;
          for (auto q = b; q < e; ++q) {
            // entries[q].local now holds the compliance slot
            y += factor->compliance(local_slots[l], entries[q].local) * entries[q].value;
          }
          NM_zentry(minv_ht, blocks.dof(g, locals[l]), row, y, zero_threshold);
        }
        b = e;
      }
    } else {
      Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(bsize, rows.size());
      for (auto& e : entries) rhs(e.local, row_slot(e.row)) += e.value;

      Eigen::MatrixXd x = factor ? factor->solve(rhs) : blocks.solve_dense(g, rhs);
      for (std::size_t k = 0; k < rows.size(); ++k) {
        for (auto lj : locals) {
          NM_zentry(minv_ht, blocks.dof(g, lj), rows[k], x(lj, k), zero_threshold);
        }
      }
    }
  }

  NM_gemm(1.0, h_matrix._m, minv_ht, 0.0, w_matrix._m);

  NM_free(minv_ht);
}

// c <- A^{-1} b for a block diagonal A, see compute_kkt_matrix.
template <match::any_matrix M, typename B>
void solve_block_diagonal_system(mat<M>& matrix, vec<B>& b, vec<B>& c,
                                 std::span<const factored_block> factored = {}) {
  const CSparseMatrix* a = NM_csc(matrix._m);
  assert(b._v->size0 == a->n && c._v->size0 == a->n);

  const double* bx = b._v->matrix0;
  double* cx = c._v->matrix0;

  block_diagonal_solver blocks(a, factored);

  for (CS_INT g = 0; g < blocks.number_of_groups(); ++g) {
    if (blocks.covered(g)) continue;  // solved with its factored block
    auto bsize = blocks.size(g);
    if (blocks.diagonal(g)) {
      auto j = blocks.dof(g, 0);
      double ajj = blocks.diagonal_value(g);
      if (ajj == 0.) {
        throw std::runtime_error("solve_block_diagonal_system: zero diagonal block");
      }
      cx[j] = bx[j] / ajj;
    } else {
      Eigen::MatrixXd rhs(bsize, 1);
      for (CS_INT l = 0; l < bsize; ++l) rhs(l, 0) = bx[blocks.dof(g, l)];

      auto* factor = blocks.factor(g);
      Eigen::MatrixXd x = factor ? factor->solve(rhs) : blocks.solve_dense(g, rhs);
      for (CS_INT l = 0; l < bsize; ++l) cx[blocks.dof(g, l)] = x(l, 0);
    }
  }
}

void display(match::any_mat auto& a) { NM_display(a._m); }

void display(match::vec auto& v) { NM_display(v._v); }

}  // namespace siconos::algebra
