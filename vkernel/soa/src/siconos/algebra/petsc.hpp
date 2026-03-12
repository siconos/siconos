#pragma once

#include <petscmat.h>
#include <petscvec.h>
#include <cassert>
#include <stdexcept>
#include <vector>

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
  static constexpr auto vncols = 1;
  static constexpr auto vnrows = 1;

  Mat _m = nullptr;
  bool _inversed = false;
  Mat _mt = nullptr;

  constexpr mat() {}

  ~mat()
  {
    if (_m) MatDestroy(&_m);
    if (_mt) MatDestroy(&_mt);
  }
};

template <typename T>
struct diag_mat : mat<T> {
  using any_mat_t = typename mat<T>::any_mat_t;
  using diag_mat_t = void;

  static constexpr auto vncols = 1;
  static constexpr auto vnrows = 1;
};

template <typename T>
struct vec {
  using vec_t = void;
  static constexpr auto vncols = 1;
  static constexpr auto vnrows = 1;

  Vec _v = nullptr;

  constexpr vec() {}

  ~vec()
  {
    if (_v) VecDestroy(&_v);
  }
};

template <>
struct vec<double> {
  using vec_t = void;
  static constexpr auto vncols = 1;
  static constexpr auto vnrows = 1;

  Vec _v = nullptr;

  constexpr vec() {}

  ~vec()
  {
    if (_v) VecDestroy(&_v);
  }
};

const auto size0(match::any_mat auto& m) {
  PetscInt rows;
  MatGetSize(m._m, &rows, nullptr);
  return rows / m.vnrows;
};
const auto size0(match::vec auto& v) {
  PetscInt size;
  VecGetSize(v._v, &size);
  return size / v.vnrows;
};
const auto size1(match::any_mat auto& m) {
  PetscInt cols;
  MatGetSize(m._m, nullptr, &cols);
  return cols / m.vncols;
};

const auto raw_size0(match::any_mat auto& m) {
  PetscInt rows;
  MatGetSize(m._m, &rows, nullptr);
  return rows;
};
const auto raw_size0(match::vec auto& v) {
  PetscInt size;
  VecGetSize(v._v, &size);
  return size;
};
const auto raw_size1(match::any_mat auto& m) {
  PetscInt cols;
  MatGetSize(m._m, nullptr, &cols);
  return cols;
};

static_assert(vec<std::vector<int>>::vncols == 1);

void resize(match::any_mat auto& m, PetscInt nrows, PetscInt ncols)
{
  if (m._m) MatDestroy(&(m._m));
  if (m._mt) MatDestroy(&(m._mt));

  m._inversed = false;

  MatCreateSeqAIJ(PETSC_COMM_SELF, nrows * m.vnrows, ncols * m.vncols, PETSC_DECIDE, nullptr, &m._m);
  // No explicit allocation for transposed
}

void resize(match::vec auto& v, PetscInt nrows)
{
  if (v._v) VecDestroy(&(v._v));
  assert(nrows * v.vnrows >= 1);
  assert(v.vncols == 1);
  VecCreateSeq(PETSC_COMM_SELF, nrows * v.vnrows, &v._v);
}

void insert(match::any_mat auto& am, match::any_mat auto& bm,
            PetscInt offset0, PetscInt offset1)
{
  PetscInt brow, bcol;
  MatGetSize(bm._m, &brow, &bcol);
  std::vector<PetscInt> rows(brow);
  std::vector<PetscInt> cols(bcol);

  for(PetscInt i=0; i<brow; ++i) rows[i]=i;
  for(PetscInt i=0; i<bcol; ++i) cols[i]=i;

  Mat mat_bm = bm._m;
  for(PetscInt i=0; i<brow; ++i)
    for(PetscInt j=0; j<bcol; ++j)
    {
      PetscScalar v;
      MatGetValues(mat_bm, 1, &i, 1, &j, &v);
      if (std::abs(v) > zero_threshold)
        MatSetValue(am._m, offset0 + i, offset1 + j, v, INSERT_VALUES);
    }
  MatAssemblyBegin(am._m, MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(am._m, MAT_FLUSH_ASSEMBLY);
}

void insert(match::vec auto& va, match::vec auto& vb,
            PetscInt offset0)
{
  PetscInt sz;
  VecGetSize(vb._v, &sz);
  std::vector<PetscScalar> buf(sz);
  VecGetArray(vb._v, buf.data());
  for (PetscInt i = 0; i < sz; ++i)
  {
    if (std::abs(buf[i]) > zero_threshold)
      VecSetValue(va._v, offset0 + i, buf[i], INSERT_VALUES);
  }
  VecRestoreArray(vb._v, buf.data());
  VecAssemblyBegin(va._v);
  VecAssemblyEnd(va._v);
}

void transpose(match::any_mat auto& m)
{
  if (!m._mt)
  {
    MatTranspose(m._m, MAT_INITIAL_MATRIX, &m._mt);
  }
}

void setup(match::any_mat auto& m) { resize(m, 1, 1); }

template <typename T>
void set_value(match::any_mat auto& m, PetscInt i,
               PetscInt j, const T& value)
{
  MatSetValue(m._m, i * m.vnrows, j * m.vncols, (PetscScalar)value, INSERT_VALUES);
}

template <typename T>
void set_value(match::vec auto& v, PetscInt i, const T& value)
{
  VecSetValue(v._v, i * v.vnrows, (PetscScalar)value, INSERT_VALUES);
}

template <typename A>
void inverse(diag_mat<A>& a)
{
  // In PETSc invert a diagonal matrix for later multiplication
  // Not recommended in practice: use MatDiagonalScale or a PC
  PetscInt nrows, ncols;
  MatGetSize(a._m, &nrows, &ncols);
  for(PetscInt i=0; i<nrows; ++i)
  {
    PetscScalar v;
    MatGetValues(a._m, 1, &i, 1, &i, &v);
    v = (std::abs(v) > zero_threshold) ? 1.0/v : 0.0;
    MatSetValue(a._m, i, i, v, INSERT_VALUES);
  }
  MatAssemblyBegin(a._m, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(a._m, MAT_FINAL_ASSEMBLY);
  a._inversed = true;
}

// b += a
template <typename T>
void add(const vec<T>& a, vec<T>& b)
{
  VecAXPY(b._v, 1.0, a._v);
}

template <typename T>
void scal(PetscScalar h, vec<T>& v)
{
  VecScale(v._v, h);
}

template <typename T>
decltype(auto) get_vector(vec<T>& v, PetscInt i)
{
  PetscScalar* arr;
  VecGetArray(v._v, &arr);
  return arr + i * v.vnrows; // user responsible for VecRestoreArray
}

// c <- a b
// Matrix Matrix
template <typename A, typename B>
void prod(mat<A>& a, mat<B>& b, mat<void>& c)
{
  if (!c._m) {
    PetscInt rows, cols;
    MatGetSize(a._m, &rows, nullptr);
    MatGetSize(b._m, nullptr, &cols);
    MatCreateSeqAIJ(PETSC_COMM_SELF, rows, cols, PETSC_DECIDE, 0, &c._m);
  }
  MatMatMult(a._m, b._m, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &c._m);
}

// Matrix Vector
template <typename A, typename B>
void prod(mat<A>& a, vec<B>& b, vec<void>& c)
{
  if (!c._v)
  {
    PetscInt sz;
    MatGetSize(a._m, nullptr, &sz);
    VecCreateSeq(PETSC_COMM_SELF, sz, &c._v);
  }
  MatMult(a._m, b._v, c._v);
}

// c <- a^t b
template <typename A, typename B>
void prodt1(const mat<A>& a, const vec<B>& b, vec<void>& c)
{
  if (!a._mt)
    throw std::runtime_error("Transpose matrix was not assembled");
  MatMult(a._mt, b._v, c._v);
}

// c <- a b^t
template <typename A, typename B>
void prodt2(const diag_mat<A>& a, const mat<B>& b, mat<void>& c)
{
  if (!b._mt)
    throw std::runtime_error("Transpose matrix was not assembled");
  MatMatMult(a._m, b._mt, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &c._m);
}

// c <- a^-1 b
template <typename A, typename B>
void solve(diag_mat<A>& a, vec<B>& b, vec<B>& c)
{
  inverse(a);
  prod(a, b, c);
}

template <typename A, typename B>
void solve(mat<A>& a, vec<B>& b, vec<B>& c)
{
  KSP ksp;
  KSPCreate(PETSC_COMM_SELF, &ksp);
  KSPSetOperators(ksp, a._m, a._m);
  KSPSetFromOptions(ksp);
  KSPSolve(ksp, b._v, c._v);
  KSPDestroy(&ksp);
}

template <typename A, typename B>
void solvet(diag_mat<A>& a, mat<B>& b, mat<void>& c)
{
  inverse(a);
  transpose(b);
  prodt2(a, b, c);
}

void display(match::any_mat auto& a)
{
  MatView(a._m, PETSC_VIEWER_STDOUT_SELF);
}

void display(match::vec auto& v)
{
  VecView(v._v, PETSC_VIEWER_STDOUT_SELF);
}

}  // namespace siconos::algebra
