#include "SimpleMatrix.h"
#include "SiconosMatrixException.h"
#include "SimpleVector.h"
#include "ioMatrix.h"
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "boost/numeric/bindings/atlas/clapack.hpp"
#include "boost/numeric/bindings/traits/std_vector.hpp"
#include "boost/numeric/bindings/traits/ublas_matrix.hpp"
#include "boost/numeric/bindings/traits/ublas_vector2.hpp" // vector considered as matrix => necessary for bindings atlas-lapack
#include "boost/numeric/bindings/atlas/cblas2.hpp"
#include "boost/numeric/bindings/atlas/cblas3.hpp"

namespace atlas = boost::numeric::bindings::atlas;

// Build a Simple Matrix from its type (ie DENSE, TRIANGULAR, BANDED, SPARSE or SYMMETRIC)
// => default (private) constructor
SimpleMatrix::SimpleMatrix(TYP typ): SiconosMatrix(false), num(1), isPLUFactorized(false), isPLUInversed(false)
{
  if (typ == DENSE)
  {
    mat.Dense = new DenseMat();
    num = 1;
  }
  else if (typ == TRIANGULAR)
  {
    mat.Triang = new TriangMat();
    num = 2;
  }
  else if (typ == SYMMETRIC)
  {
    mat.Sym = new SymMat();
    num = 3;
  }
  else if (typ == SPARSE)
  {
    mat.Sparse = new SparseMat();
    num = 4;
  }
  else if (typ == BANDED)
  {
    mat.Banded = new BandedMat();
    num = 5;
  }
  else if (typ == ZERO)
  {
    mat.Zero = new ZeroMat();
  }
  else if (typ == IDENTITY)
  {
    mat.Identity = new IdentityMat();
  }
  else
    SiconosMatrixException::selfThrow("constructor(TYP) : invalid type given");
  dim[0] = 0;
  dim[1] = 0;
}

// Copy constructors
SimpleMatrix::SimpleMatrix(const SimpleMatrix &smat): SiconosMatrix(false), num(smat.getNum()), isPLUFactorized(false), isPLUInversed(false)
{
  if (num == 1)
    mat.Dense = new DenseMat(*smat.getDensePtr());

  else if (num == 2)
    mat.Triang = new TriangMat(*smat.getTriangPtr());

  else if (num == 3)
    mat.Sym = new SymMat(*smat.getSymPtr());

  else if (smat.getNum() == 4)
    mat.Sparse = new SparseMat(*smat.getSparsePtr());

  else if (smat.getNum() == 5)
    mat.Banded = new BandedMat(*smat.getBandedPtr());

  else if (smat.getNum() == 6)
    mat.Zero = new ZeroMat(*smat.getZeroPtr());

  else if (smat.getNum() == 7)
    mat.Identity = new IdentityMat(*smat.getIdentityPtr());

  else
    SiconosMatrixException::selfThrow("constructor(const SimpleMatrix) : invalid parameter given");
  dim[0] = smat.size(0);
  dim[1] = smat.size(1);
}

SimpleMatrix::SimpleMatrix(const SiconosMatrix &smat): SiconosMatrix(false), num(smat.getNum()), isPLUFactorized(false), isPLUInversed(false)
{
  assert(smat.isBlock() == false);
  if (num == 1)
    mat.Dense = new DenseMat(*smat.getDensePtr());

  else if (num == 2)
    mat.Triang = new TriangMat(*smat.getTriangPtr());

  else if (num == 3)
    mat.Sym = new SymMat(*smat.getSymPtr());

  else if (num == 4)
    mat.Sparse = new SparseMat(*smat.getSparsePtr());

  else if (num == 5)
    mat.Banded = new BandedMat(*smat.getBandedPtr());

  else if (num == 6)
    mat.Zero = new ZeroMat(*smat.getZeroPtr());

  else if (num == 7)
    mat.Identity = new IdentityMat(*smat.getIdentityPtr());

  else
    SiconosMatrixException::selfThrow("constructor(const SiconosMatrix) : invalid parameter given");
  dim[0] = smat.size(0);
  dim[1] = smat.size(1);
}

SimpleMatrix::SimpleMatrix(const std::string& in, const SiconosMatrix &smat): SiconosMatrix(false), num(smat.getNum()), isPLUFactorized(false), isPLUInversed(false)
{
  if (in != "transpose")
    SiconosMatrixException::selfThrow("constructor(trans,const SiconosMatrix) : wrong value for trans.");
  assert(smat.isBlock() == false);
  if (num == 1)
    mat.Dense = new DenseMat(ublas::trans(*smat.getDensePtr()));

  else if (num == 3)
    mat.Sym = new SymMat(*smat.getSymPtr());

  else if (num == 4)
    mat.Sparse = new SparseMat(ublas::trans(*smat.getSparsePtr()));

  else if (num == 5)
    mat.Banded = new BandedMat(ublas::trans(*smat.getBandedPtr()));

  else if (num == 6)
    mat.Zero = new ZeroMat(*smat.getZeroPtr());

  else if (num == 7)
    mat.Identity = new IdentityMat(*smat.getIdentityPtr());

  else
    SiconosMatrixException::selfThrow("constructor(const SiconosMatrix) : invalid parameter given");

  dim[0] = smat.size(1);
  dim[1] = smat.size(0);
}


SimpleMatrix::SimpleMatrix(unsigned int row, unsigned int col, TYP typ, unsigned int upper, unsigned int lower):
  SiconosMatrix(false), num(1), isPLUFactorized(false), isPLUInversed(false)
{
  dim[0] = row;
  dim[1] = col;
  if (typ == DENSE)
  {
    mat.Dense = new DenseMat(row, col);
    num = 1;
    zero();
  }
  else if (typ == TRIANGULAR)
  {
    mat.Triang = new TriangMat(row, col);
    num = 2;
    zero();
  }
  else if (typ == SYMMETRIC)
  {
    mat.Sym = new SymMat(row, col);
    num = 3;
    zero();
  }
  else if (typ == SPARSE)
  {
    mat.Sparse = new SparseMat(row, col, upper);
    num = 4;
    zero();
  }
  else if (typ == BANDED)
  {
    mat.Banded = new BandedMat(row, col, upper, lower);
    num = 5;
    zero();
  }
  else if (typ == ZERO)
  {
    mat.Zero = new ZeroMat(row, col);
    num = 6;
  }
  else if (typ == IDENTITY)
  {
    mat.Identity = new IdentityMat(row, col);
    num = 7;
  }
  else
    SiconosMatrixException::selfThrow("constructor(TYP type, unsigned int row, unsigned int col) : invalid type or dimensions given");
}

SimpleMatrix::SimpleMatrix(const DenseMat& m): SiconosMatrix(false), num(1), isPLUFactorized(false), isPLUInversed(false)
{
  mat.Dense = new DenseMat(m);
  dim[0] = m.size1();
  dim[1] = m.size2();
}

SimpleMatrix::SimpleMatrix(const TriangMat& m): SiconosMatrix(false), num(2), isPLUFactorized(false), isPLUInversed(false)
{
  mat.Triang = new TriangMat(m);
  dim[0] = m.size1();
  dim[1] = m.size2();
}

SimpleMatrix::SimpleMatrix(const SymMat& m): SiconosMatrix(false), num(3), isPLUFactorized(false), isPLUInversed(false)
{
  mat.Sym = new SymMat(m);
  dim[0] = m.size1();
  dim[1] = m.size2();
}

SimpleMatrix::SimpleMatrix(const SparseMat& m): SiconosMatrix(false), num(4), isPLUFactorized(false), isPLUInversed(false)
{
  mat.Sparse = new SparseMat(m);
  dim[0] = m.size1();
  dim[1] = m.size2();
}

SimpleMatrix::SimpleMatrix(const BandedMat& m): SiconosMatrix(false), num(5), isPLUFactorized(false), isPLUInversed(false)
{
  mat.Banded = new BandedMat(m);
  dim[0] = m.size1();
  dim[1] = m.size2();
}

SimpleMatrix::SimpleMatrix(const ZeroMat& m): SiconosMatrix(false), num(6), isPLUFactorized(false), isPLUInversed(false)
{
  mat.Zero = new ZeroMat(m);
  dim[0] = m.size1();
  dim[1] = m.size2();
}

SimpleMatrix::SimpleMatrix(const IdentityMat& m): SiconosMatrix(false), num(7), isPLUFactorized(false), isPLUInversed(false)
{
  mat.Identity = new IdentityMat(m);
  dim[0] = m.size1();
  dim[1] = m.size2();
}

SimpleMatrix::SimpleMatrix(const std::vector<double> &v, unsigned int row, unsigned int col, TYP typ, unsigned int lower, unsigned int upper):
  SiconosMatrix(false), num(1), isPLUFactorized(false), isPLUInversed(false)
{
  if (((v.size()) != (unsigned int)row * col && (typ != SYMMETRIC && typ != BANDED)) || (v.size() != (unsigned int)row * row && typ == SYMMETRIC) || (typ == BANDED && ((v.size()) != (unsigned int)(std::max)(row, col) * (lower + 1 + upper))))
    SiconosMatrixException::selfThrow("constructor(TYP, const std::vector<double>, int, int) : invalid vector size");

  if (typ == DENSE)
  {
    mat.Dense = new DenseMat(row, col, v);
    num = 1;
    dim[0] = row;
    dim[1] = col;
  }
  else if (typ == TRIANGULAR)
  {
    mat.Triang = new TriangMat(row, col, v);
    num = 2;
    dim[0] = row;
    dim[1] = col;
  }
  else if (typ == SYMMETRIC)
  {
    mat.Sym = new SymMat(row, v);
    num = 3;
    dim[0] = row;
    dim[1] = row;
  }
  else if (typ == SPARSE)
  {
    SiconosMatrixException::selfThrow("constructor(TYP, const std::vector<double>, int row, int col, int lower, int upper) : warning -- use constructor(const SparseMat &m) or constructor(TYP, int row, int col) with TYP = SPARSE");

  }
  else if (typ == BANDED)
  {
    mat.Banded = new BandedMat(row, col, lower, upper, v);
    num = 5;
    dim[0] = row;
    dim[1] = col;
  }
  else
    SiconosMatrixException::selfThrow("constructor(TYP, const std::vector<double>, int, int) : invalid type of matrix given");
}

SimpleMatrix::SimpleMatrix(const std::string &file, bool ascii): SiconosMatrix(false), num(1), isPLUFactorized(false), isPLUInversed(false)
{
  mat.Dense = new DenseMat();
  if (ascii)
  {
    ioMatrix io(file, "ascii");
    io.read(*this);
  }
  else
  {
    ioMatrix io(file, "binary");
    io.read(*this);
  }
  dim[0] = (mat.Dense)->size1();
  dim[1] = (mat.Dense)->size2();

}

/****************************** DESTRUCTOR  ****************************/
SimpleMatrix::~SimpleMatrix()
{
  if (num == 1)
    delete(mat.Dense);
  else if (num == 2)
    delete(mat.Triang);
  else if (num == 3)
    delete(mat.Sym);
  else if (num == 4)
    delete(mat.Sparse);
  else if (num == 5)
    delete(mat.Banded);
  else if (num == 6)
    delete(mat.Zero);
  else if (num == 7)
    delete(mat.Identity);
}

void SimpleMatrix::computeDim()
{
  if (num == 1)
  {
    dim[0] = (*mat.Dense).size1();
    dim[1] = (*mat.Dense).size2();
  }
  else if (num == 2)
  {
    dim[0] = (*mat.Triang).size1();
    dim[1] = (*mat.Triang).size2();
  }
  else if (num == 3)
  {
    dim[0] = (*mat.Sym).size1();
    dim[1] = (*mat.Sym).size2();
  }
  else if (num == 4)
  {
    dim[0] = (*mat.Sparse).size1();
    dim[1] = (*mat.Sparse).size2();
  }
  else if (num == 5)
  {
    dim[0] = (*mat.Banded).size1();
    dim[1] = (*mat.Banded).size2();
  }
  else if (num == 6)
  {
    dim[0] = (*mat.Zero).size1();
    dim[1] = (*mat.Zero).size2();
  }
  else if (num == 7)
  {
    dim[0] = (*mat.Identity).size1();
    dim[1] = (*mat.Identity).size2();
  }
}

void SimpleMatrix::resize(unsigned int row, unsigned int col, unsigned int lower, unsigned int upper, bool preserve)
{

  if (num == 1)
  {
    (*mat.Dense).resize(row, col, preserve);
  }
  else if (num == 2)
  {
    (*mat.Triang).resize(row, col, preserve);
  }
  else if (num == 3)
  {
    (*mat.Sym).resize(row, col, preserve);
  }
  else if (num == 4)
  {
    (*mat.Sparse).resize(row, col, preserve);
  }
  else if (num == 5)
  {
    (*mat.Banded).resize(row, col, lower, upper, preserve);
  }
  else if (num == 6)
  {
    (*mat.Zero).resize(row, col, preserve);
  }
  else if (num == 7)
  {
    (*mat.Identity).resize(row, col, preserve);
  }
  resetLU();
  dim[0] = row;
  dim[1] = col;
}

const DenseMat SimpleMatrix::getDense(unsigned int row, unsigned int col)const
{

  if (num != 1)
    SiconosMatrixException::selfThrow("DenseMat getDense(unsigned int row, unsigned int col) : the current matrix is not a Dense matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("DenseMat getDense(unsigned int row, unsigned int col) : row or col not equal to 0.0");

  return *mat.Dense;
}

const TriangMat SimpleMatrix::getTriang(unsigned int row, unsigned int col)const
{

  if (num != 2)
    SiconosMatrixException::selfThrow("TriangMat getTriang(unsigned int row, unsigned int col) : the current matrix is not a Triangular matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("TriangMat getTriang(unsigned int row, unsigned int col) : row or col not equal to 0.0");

  return *mat.Triang;
}

const SymMat SimpleMatrix::getSym(unsigned int row, unsigned int col)const
{

  if (num != 3)
    SiconosMatrixException::selfThrow("SymMat getSym(unsigned int row, unsigned int col) : the current matrix is not a Symmetric matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("SymMat getSym(unsigned int row, unsigned int col) : row or col not equal to 0.0");

  return *mat.Sym;
}

const SparseMat SimpleMatrix::getSparse(unsigned int row, unsigned int col)const
{

  if (num != 4)
    SiconosMatrixException::selfThrow("SparseMat getSparse(unsigned int row, unsigned int col) : the current matrix is not a Sparse matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("SparseMat getSparse(unsigned int row, unsigned int col) : row or col not equal to 0.0");

  return *mat.Sparse;
}

const BandedMat SimpleMatrix::getBanded(unsigned int row, unsigned int col)const
{

  if (num != 5)
    SiconosMatrixException::selfThrow("BandedMat getBanded(unsigned int row, unsigned int col) : the current matrix is not a Banded matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("BandedMat getBanded(unsigned int row, unsigned int col) : row or col not equal to 0.0");

  return *mat.Banded;
}

const ZeroMat SimpleMatrix::getZero(unsigned int row, unsigned int col)const
{

  if (num != 6)
    SiconosMatrixException::selfThrow("ZeroMat getZero(unsigned int row, unsigned int col) : the current matrix is not a Zero matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("ZeroMat getZero(unsigned int row, unsigned int col) : row or col not equal to 0.0");

  return *mat.Zero;
}

const IdentityMat SimpleMatrix::getIdentity(unsigned int row, unsigned int col)const
{

  if (num != 7)
    SiconosMatrixException::selfThrow("IdentityMat getIdentity(unsigned int row, unsigned int col) : the current matrix is not a Identity matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("IdentityMat getIdentity(unsigned int row, unsigned int col) : row or col not equal to 0.0");

  return *mat.Identity;
}

DenseMat* SimpleMatrix::getDensePtr(unsigned int row, unsigned int col)const
{

  if (num != 1)
    SiconosMatrixException::selfThrow("DenseMat* getDensePtr(unsigned int row, unsigned int col) : the current matrix is not a Dense matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("DenseMat* getDensePtr(unsigned int row, unsigned int col) : row or col not equal to 0.0");

  return mat.Dense;
}

TriangMat* SimpleMatrix::getTriangPtr(unsigned int row, unsigned int col)const
{

  if (num != 2)
    SiconosMatrixException::selfThrow("TriangMat* getTriangPtr(unsigned int row, unsigned int col) : the current matrix is not a Triangular matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("TriangMat* getTriangPtr(unsigned int row, unsigned int col) : row or col not equal to 0.0");

  return mat.Triang;
}

SymMat* SimpleMatrix::getSymPtr(unsigned int row, unsigned int col)const
{

  if (num != 3)
    SiconosMatrixException::selfThrow("SymMat* getSymPtr(unsigned int row, unsigned int col) : the current matrix is not a Symmetric matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("SymMat* getSymPtr(unsigned int row, unsigned int col) : row or col not equal to 0.0");

  return mat.Sym;
}

SparseMat* SimpleMatrix::getSparsePtr(unsigned int row, unsigned int col)const
{

  if (num != 4)
    SiconosMatrixException::selfThrow("SparseMat* getSparsePtr(unsigned int row, unsigned int col) : the current matrix is not a Sparse matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("SparseMat* getSparsePtr(unsigned int row, unsigned int col) : row or col not equal to 0.0");

  return mat.Sparse;
}

BandedMat* SimpleMatrix::getBandedPtr(unsigned int row, unsigned int col)const
{

  if (num != 5)
    SiconosMatrixException::selfThrow("BandedMat* getBandedPtr(unsigned int row, unsigned int col) : the current matrix is not a Banded matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("BandedMat* getBandedPtr(unsigned int row, unsigned int col) : row or col not equal to 0.0");

  return mat.Banded;
}

ZeroMat* SimpleMatrix::getZeroPtr(unsigned int row, unsigned int col)const
{

  if (num != 6)
    SiconosMatrixException::selfThrow("ZeroMat* getZeroPtr(unsigned int row, unsigned int col) : the current matrix is not a Zero matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("ZeroMat* getZeroPtr(unsigned int row, unsigned int col) : row or col not equal to 0.0");

  return mat.Zero;
}

IdentityMat* SimpleMatrix::getIdentityPtr(unsigned int row, unsigned int col)const
{

  if (num != 7)
    SiconosMatrixException::selfThrow("IdentityMat* getIdentityPtr(unsigned int row, unsigned int col) : the current matrix is not a Identity matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("IdentityMat* getIdentityPtr(unsigned int row, unsigned int col) : row or col not equal to 0.0");

  return mat.Identity;
}

const BlocksMat SimpleMatrix::getAllBlocks(void)const
{
  SiconosMatrixException::selfThrow("BlocksMat getAllBlocks : getAllBlocks is forbidden for SimpleMatrix");
  const BlocksMat a;
  return a;
}

void SimpleMatrix::matrixCopy(const SiconosMatrix &m, unsigned int i, unsigned int j)
{

  if (num != 1)
    SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : the current matrix is not dense, a copy into its data may change its type.");

  if (i >= size(0) || i < 0)
    SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : row_min given is out of range");

  if (j >= size(1) || j < 0)
    SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : col_min given is out of range");

  unsigned int num2 = m.getNum();
  unsigned int row_max = i, col_max = j;

  unsigned int sizeM1, sizeM2;
  if (num2 == 1)
  {
    sizeM1 = (m.getDensePtr())->size1();
    sizeM2 = (m.getDensePtr())->size2();
  }
  else if (num2 == 2)
  {
    sizeM1 = (m.getTriangPtr())->size1();
    sizeM2 = (m.getTriangPtr())->size2();
  }
  else if (num2 == 3)
  {
    sizeM1 = (m.getSymPtr())->size1();
    sizeM2 = (m.getSymPtr())->size2();
  }
  else if (num2 == 4)
  {
    sizeM1 = (m.getSparsePtr())->size1();
    sizeM2 = (m.getSparsePtr())->size2();
  }
  else if (num2 == 5)
  {
    sizeM1 = (m.getBandedPtr())->size1();
    sizeM2 = (m.getBandedPtr())->size2();
  }
  else if (num2 == 6)
  {
    sizeM1 = (m.getZeroPtr())->size1();
    sizeM2 = (m.getZeroPtr())->size2();
  }
  else if (num2 == 7)
  {
    sizeM1 = (m.getIdentityPtr())->size1();
    sizeM2 = (m.getIdentityPtr())->size2();
  }

  row_max += sizeM1;
  col_max += sizeM2;

  if (row_max > size(0) || row_max < 0)
    SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");
  if (col_max > size(1) || col_max < 0)
    SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  if (num2 == 1)
    ublas::subrange(*mat.Dense, i, row_max, j, col_max) = *m.getDensePtr();
  else if (num2 == 2)
    ublas::subrange(*mat.Dense, i, row_max, j, col_max) = *m.getTriangPtr();
  else if (num2 == 3)
    ublas::subrange(*mat.Dense, i, row_max, j, col_max) = *m.getSymPtr();
  else if (num2 == 4)
    ublas::subrange(*mat.Dense, i, row_max, j, col_max) = *m.getSparsePtr();
  else if (num2 == 5)
    ublas::subrange(*mat.Dense, i, row_max, j, col_max) = *m.getBandedPtr();
  else if (num2 == 6)
    ublas::subrange(*mat.Dense, i, row_max, j, col_max) = *m.getZeroPtr();
  else if (num2 == 7)
    ublas::subrange(*mat.Dense, i, row_max, j, col_max) = *m.getIdentityPtr();
}

void SimpleMatrix::getBlock(unsigned int row_min, unsigned int col_min, SiconosMatrix &m)const
{
  // We only accept dense matrix for m.
  if (m.getNum() != 1)
    SiconosMatrixException::selfThrow("getBlock(i,j,m) : m must be a dense matrix.");

  if (row_min >= size(0) || row_min < 0)
    SiconosMatrixException::selfThrow("getBlock : row_min given is out of range");

  if (col_min >= size(1) || col_min < 0)
    SiconosMatrixException::selfThrow("getBlock : col_min given is out of range");
  unsigned int row_max, col_max;
  row_max = m.getDensePtr()->size1() + row_min;
  col_max = m.getDensePtr()->size2() + col_min;

  if (row_max > size(0) || row_max < 0)
    SiconosMatrixException::selfThrow("getBlock : inconsistent sizes");

  if (col_max > size(1) || col_max < 0)
    SiconosMatrixException::selfThrow("getBlock : inconsistent sizes");

  DenseMat * q = m.getDensePtr();
  if (num == 1)
    *q = ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max);
  else if (num == 2)
    *q = ublas::subrange(*mat.Triang, row_min, row_max, col_min, col_max);
  else if (num == 3)
    *q = ublas::subrange(*mat.Sym, row_min, row_max, col_min, col_max);
  else if (num == 4)
    *q = ublas::subrange(*mat.Sparse, row_min, row_max, col_min, col_max);
  else if (num == 5)
    *q = ublas::subrange(*mat.Banded, row_min, row_max, col_min, col_max);
  else if (num == 6)
    *q = ublas::subrange(*mat.Zero, row_min, row_max, col_min, col_max);
  else if (num == 7)
    *q = ublas::subrange(*mat.Identity, row_min, row_max, col_min, col_max);
}

const std::deque<bool> SimpleMatrix::getBlockAllocated(void)const
{
  SiconosMatrixException::selfThrow("std::deque<bool> getBlockAllocated : getBlockAllocated is forbidden for SimpleMatrix");
  std::deque<bool> tmp;
  return tmp; // to avoid warning
}

void SimpleMatrix::getRow(unsigned int r, SimpleVector &vect) const
{

  if (r >= size(0) || r < 0)
    SiconosMatrixException::selfThrow("getRow : row is out of range");

  if (vect.size() != size(1))
    SiconosMatrixException::selfThrow("getRow : inconsistent sizes");

  if (num == 1)
  {
    *(vect.getDensePtr()) = ublas::row(*mat.Dense, r);
  }
  else if (num == 2)
  {
    *(vect.getDensePtr()) = ublas::row(*mat.Triang, r);
  }
  else if (num == 3)
  {
    *(vect.getDensePtr()) = ublas::row(*mat.Sym, r);
  }
  else if (num == 4)
  {
    *(vect.getDensePtr()) = ublas::row(*mat.Sparse, r);
  }
  else if (num == 5)
  {
    *(vect.getDensePtr()) = ublas::row(*mat.Banded, r);
  }
  else if (num == 6)
  {
    *(vect.getDensePtr()) = ublas::row(*mat.Zero, r);
  }
  else if (num == 7)
  {
    *(vect.getDensePtr()) = ublas::row(*mat.Identity, r);
  }
}

void SimpleMatrix::setRow(unsigned int r, const SimpleVector &vect)
{
  unsigned int numV = vect.getNum();
  if (r >= size(0) || r < 0)
    SiconosMatrixException::selfThrow("setRow : row is out of range");

  if (vect.size() != size(1))
    SiconosMatrixException::selfThrow("setRow : inconsistent sizes");

  if (num == 1)
  {
    if (numV == 1)
    {
      ublas::row(*mat.Dense, r) = *vect.getDensePtr();
    }
    else if (numV == 2)
    {
      ublas::row(*mat.Dense, r) = *vect.getSparsePtr();
    }
  }
  else if (num == 2)
  {
    if (numV == 1)
    {
      ublas::row(*mat.Triang, r) = *vect.getDensePtr();
    }
    else if (numV == 2)
    {
      ublas::row(*mat.Triang, r) = *vect.getSparsePtr();
    }
  }
  else if (num == 3)
  {
    if (numV == 1)
    {
      ublas::row(*mat.Sym, r) = *vect.getDensePtr();
    }
    else if (numV == 2)
    {
      ublas::row(*mat.Sym, r) = *vect.getSparsePtr();
    }
  }
  else if (num == 4)
  {
    if (numV == 1)
    {
      ublas::row(*mat.Sparse, r) = *vect.getDensePtr();
    }
    else if (numV == 2)
    {
      ublas::row(*mat.Sparse, r) = *vect.getSparsePtr();
    }
  }
  else if (num == 5)
  {
    if (numV == 1)
    {
      ublas::row(*mat.Sparse, r) = *vect.getDensePtr();
    }
    else if (numV == 2)
    {
      ublas::row(*mat.Banded, r) = *vect.getSparsePtr();
    }
  }
  else // if(num==6 || num == 7)
    SiconosMatrixException::selfThrow("setRow : forbidden for this type of matrix, num = " + num);
  resetLU();
}

void SimpleMatrix::getCol(unsigned int r, SimpleVector &vect)const
{
  if (r >= size(1) || r < 0)
    SiconosMatrixException::selfThrow("getCol : col is out of range");

  if (vect.size() != size(0))
    SiconosMatrixException::selfThrow("getCol : inconsistent sizes");

  if (num == 1)
  {
    *(vect.getDensePtr()) = ublas::column(*mat.Dense, r);
  }
  else if (num == 2)
  {
    *(vect.getDensePtr()) = ublas::column(*mat.Triang, r);
  }
  else if (num == 3)
  {
    *(vect.getDensePtr()) = ublas::column(*mat.Sym, r);
  }
  else if (num == 4)
  {
    *(vect.getDensePtr()) = ublas::column(*mat.Sparse, r);
  }
  else if (num == 5)
  {
    *(vect.getDensePtr()) = ublas::column(*mat.Banded, r);
  }
  else if (num == 6)
  {
    *(vect.getDensePtr()) = ublas::column(*mat.Zero, r);
  }
  else if (num == 7)
  {
    *(vect.getDensePtr()) = ublas::column(*mat.Identity, r);
  }
}

void SimpleMatrix::setCol(unsigned int r, const SimpleVector &vect)
{

  unsigned int numV = vect.getNum();

  if (r >= size(1) || r < 0)
    SiconosMatrixException::selfThrow("setCol : col is out of range");

  if (vect.size() != size(0))
    SiconosMatrixException::selfThrow("setCol : inconsistent sizes");

  if (num == 1)
  {
    if (numV == 1)
    {
      ublas::column(*mat.Dense, r) = *vect.getDensePtr();
    }
    else if (numV == 2)
    {
      ublas::column(*mat.Dense, r) = *vect.getSparsePtr();
    }
  }
  else if (num == 2)
  {
    if (numV == 1)
    {
      ublas::column(*mat.Triang, r) = *vect.getDensePtr();
    }
    else if (numV == 2)
    {
      ublas::column(*mat.Triang, r) = *vect.getSparsePtr();
    }
  }
  else if (num == 3)
  {
    if (numV == 1)
    {
      ublas::column(*mat.Sym, r) = *vect.getDensePtr();
    }
    else if (numV == 2)
    {
      ublas::column(*mat.Sym, r) = *vect.getSparsePtr();
    }
  }
  else if (num == 4)
  {
    if (numV == 1)
    {
      ublas::column(*mat.Sparse, r) = *vect.getDensePtr();
    }
    else if (numV == 2)
    {
      ublas::column(*mat.Sparse, r) = *vect.getSparsePtr();
    }
  }
  else if (num == 5)
  {
    if (numV == 1)
    {
      ublas::column(*mat.Banded, r) = *vect.getDensePtr();
    }
    else if (numV == 2)
    {
      ublas::column(*mat.Banded, r) = *vect.getSparsePtr();
    }
  }
  else // if(num==6 || num == 7)
    SiconosMatrixException::selfThrow("setCol : forbidden for this type of matrix, num = " + num);
  resetLU();
}

const double SimpleMatrix::normInf(void)const
{
  double d = 0;
  if (num == 1)
    d = norm_inf(*mat.Dense);
  else if (num == 2)
    d = norm_inf(*mat.Triang);
  else if (num == 3)
    d = norm_inf(*mat.Sym);
  else if (num == 4)
    d = norm_inf(*mat.Sparse);
  else if (num == 5)
    d = norm_inf(*mat.Banded);
  else if (num == 6)
    d = 0;
  else if (num == 7)
    d = 1;
  return d;
}

void SimpleMatrix::trans()
{
  switch (num)
  {
  case 1:
    *mat.Dense = ublas::trans(*mat.Dense);
    break;
  case 2:
    SiconosMatrixException::selfThrow("SimpleMatrix::trans() failed, the matrix is triangular matrix and can not be transposed in place.");
    break;
  case 3:
    break;
  case 4:
    *mat.Sparse = ublas::trans(*mat.Sparse);
  case 5:
    *mat.Banded = ublas::trans(*mat.Banded);
    break;
  case 6:
    break;
  case 7:
    break;
  }
  unsigned int tmp = dim[0];
  dim[0] = dim[1];
  dim[1] = tmp;
}

void SimpleMatrix::trans(const SiconosMatrix &m)
{
  if (&m == this)
    SiconosMatrixException::selfThrow("SimpleMatrix::trans(m) failed, m = this, use this->trans().");

  unsigned int numM = m.getNum();
  switch (numM)
  {
  case 1:
    if (num != 1)
      SiconosMatrixException::selfThrow("SimpleMatrix::trans(m) failed, try to transpose a dense matrix into another type.");
    noalias(*mat.Dense) = ublas::trans(*m.getDensePtr());
    break;
  case 2:
    if (num != 1)
      SiconosMatrixException::selfThrow("SimpleMatrix::trans(m) failed, try to transpose a triangular matrix into a non-dense one.");
    noalias(*mat.Dense) = ublas::trans(*m.getTriangPtr());
    break;
  case 3:
    *this = m;
    break;
  case 4:
    if (num == 1)
      noalias(*mat.Dense) = ublas::trans(*m.getSparsePtr());
    else if (num == 4)
      noalias(*mat.Sparse) = ublas::trans(*m.getSparsePtr());
    else
      SiconosMatrixException::selfThrow("SimpleMatrix::trans(m) failed, try to transpose a sparse matrix into a forbidden type (not dense nor sparse).");
    break;
  case 5:
    if (num == 1)
      noalias(*mat.Dense) = ublas::trans(*m.getBandedPtr());
    else if (num == 5)
      noalias(*mat.Banded) = ublas::trans(*m.getBandedPtr());
    else
      SiconosMatrixException::selfThrow("SimpleMatrix::trans(m) failed, try to transpose a banded matrix into a forbidden type (not dense nor banded).");
    break;
  case 6:
    *this = m;
    break;
  case 7:
    *this = m;
  }
  unsigned int tmp = dim[0];
  dim[0] = dim[1];
  dim[1] = tmp;
}

void SimpleMatrix::display(void)const
{
  std::cout << "mat: " ;
  if (num == 1)
    std::cout << *mat.Dense << std::endl;
  else if (num == 2)
    std::cout << *mat.Triang << std::endl;
  else if (num == 3)
    std::cout << *mat.Sym << std::endl;
  else if (num == 4)
    std::cout << *mat.Sparse << std::endl;
  else if (num == 5)
    std::cout << *mat.Banded << std::endl;
  else if (num == 6)
    std::cout << *mat.Zero << std::endl;
  else if (num == 7)
    std::cout << *mat.Identity << std::endl;
}

double* SimpleMatrix::getArray(unsigned int, unsigned int) const
{
  double * d = 0;
  if (num == 1)
    d = &(((*mat.Dense).data())[0]);
  else if (num == 2)
    d = &(((*mat.Triang).data())[0]);
  else if (num == 3)
    d = &(((*mat.Sym).data())[0]);
  else if (num == 4)
    SiconosMatrixException::selfThrow("SimpleMatrix::getArray() : not yet implemented for sparse matrix.");
  else if (num == 6)
  {
    ZeroMat::iterator1 it = (*mat.Zero).begin1();
    d = const_cast<double*>(&(*it));
  }
  else if (num == 7)
  {
    IdentityMat::iterator1 it = (*mat.Identity).begin1();
    d = const_cast<double*>(&(*it));
  }
  else
    d = &(((*mat.Banded).data())[0]);

  return d;
}

void SimpleMatrix::zero()
{
  unsigned int size1 = dim[0];
  unsigned int size2 = dim[1];
  if (num == 1)
    *mat.Dense = ublas::zero_matrix<double>(size1, size2);

  else if (num == 2)
    *mat.Triang = ublas::zero_matrix<double>(size1, size2);

  else if (num == 3)
    *mat.Sym = ublas::zero_matrix<double>(size1, size2);

  else if (num == 4)
    *mat.Sparse = ublas::zero_matrix<double>(size1, size2);

  else if (num == 5)
    *mat.Banded = ublas::zero_matrix<double>(size1, size2);

  else if (num == 7)
    SiconosMatrixException::selfThrow("SimpleMatrix::zero() : you can not set to zero a matrix of type Identity!.");
  resetLU();
  // if num == 6: nothing
}

void SimpleMatrix::eye(void)
{
  unsigned int size1 = dim[0];
  unsigned int size2 = dim[1];
  if (num == 1)
    *mat.Dense = ublas::identity_matrix<double>(size1, size2);

  else if (num == 2)
    *mat.Triang = ublas::identity_matrix<double>(size1, size2);

  else if (num == 3)
    *mat.Sym = ublas::identity_matrix<double>(size1, size2);

  else if (num == 4)
    *mat.Sparse = ublas::identity_matrix<double>(size1, size2);

  else if (num == 5)
    *mat.Banded = ublas::identity_matrix<double>(size1, size2);

  else if (num == 6)
    SiconosMatrixException::selfThrow("SimpleMatrix::eye() : you can not set to identity a matrix of type Zero!.");
  resetLU();
}

/***************************** OPERATORS ******************************/

double SimpleMatrix::getValue(unsigned int row, unsigned int col)
{
  if (row >= size(0) || col >= size(1))
    SiconosMatrixException::selfThrow("SimpleMatrix:getValue(index) : Index out of range");

  if (num == 1)
    return (*mat.Dense)(row, col);
  else if (num == 2)
    return (*mat.Triang)(row, col);
  else if (num == 3)
    return (*mat.Sym)(row, col);
  else if (num == 4)
    return (*mat.Sparse)(row, col);
  else if (num == 5)
    return (*mat.Banded)(row, col);
  else if (num == 6)
    return 0;
  else //if (num==7)
    return(row == col);
}

void SimpleMatrix::setValue(unsigned int row, unsigned int col, double value)
{
  if (row >= size(0) || col >= size(1))
    SiconosMatrixException::selfThrow("SimpleMatrix:setValue : Index out of range");

  if (num == 1)
    (*mat.Dense)(row, col) = value;
  else if (num == 2)
    (*mat.Triang)(row, col) = value;
  else if (num == 3)
    (*mat.Sym)(row, col) = value ;
  else if (num == 4)
    (*mat.Sparse)(row, col) = value;
  else if (num == 5)
    (*mat.Banded)(row, col) = value;
  else if (num == 6 || num == 7)
    SiconosMatrixException::selfThrow("SimpleMatrix:setValue : forbidden for Identity or Zero type matrices.");
  resetLU();

}

double& SimpleMatrix::operator()(unsigned int row, unsigned int col)
{
  if (row >= dim[0] || col >= dim[1])
    SiconosMatrixException::selfThrow("SimpleMatrix:operator() : Index out of range");

  if (num == 1)
    return (*mat.Dense)(row, col);
  else if (num == 2)
    return (*mat.Triang)(row, col);
  else if (num == 3)
    return (*mat.Sym)(row, col);
  else if (num == 4)
  {
    double *d = (*mat.Sparse).find_element(row, col);
    double & ref = *d;
    return ref;
  }
  else if (num == 5)
    return (*mat.Banded)(row, col);
  else if (num == 6)
    return const_cast<double&>((*mat.Zero)(row, col));
  else // i(num==7)
    return const_cast<double&>((*mat.Identity)(row, col));
}

double SimpleMatrix::operator()(unsigned int row, unsigned int col) const
{
  if (row >= dim[0] || col >= dim[1])
    SiconosMatrixException::selfThrow("SimpleMatrix:operator() : Index out of range");
  double d = 0;
  switch (num)
  {
  case 1:
    d = (*mat.Dense)(row, col);
    break;
  case 2:
    d = (*mat.Triang)(row, col);
    break;
  case 3:
    d = (*mat.Sym)(row, col);
    break;
  case 4:
    d = (*mat.Sparse)(row, col);
    break;
  case 5:
    d = (*mat.Banded)(row, col);
    break;
  case 6:
    d = 0.0;
  case 7:
    d = (row == col);
  default:
    SiconosMatrixException::selfThrow("op() (unsigned int, unsigned int) : invalid type of matrix");
    break;
  }
  return d;
}

SimpleMatrix& SimpleMatrix::operator = (const SimpleMatrix& m)
{
  // Warning!!! If sizes are inconsistent between m and this, boost operator = results in resize of this to the dim of m !!!
  // Add an exception to prevent this???
  //  if(dim[0]!=m.size(0) || dim[1] !=m.size(1))
  //    SiconosMatrixException::selfThrow("SimpleMatrix operator = failed. Inconsistent sizes.");
  // if exception, strange behavio in things like A = A*B ...

  unsigned int numM = m.getNum();
  switch (num)
  {
  case 1:
    switch (numM)
    {
    case 1:
      *mat.Dense = *m.getDensePtr();
      break;
    case 2:
      *mat.Dense = *m.getTriangPtr();
      break;
    case 3:
      *mat.Dense = *m.getSymPtr();
      break;
    case 4:
      *mat.Dense = *m.getSparsePtr();
      break;
    case 5:
      *mat.Dense = *m.getBandedPtr();
      break;
    case 6:
      *mat.Dense = *m.getZeroPtr(); // warning: this is not equivalent to a zero() call, because boost = results in resizing if required.
      break;
    case 7:
      *mat.Dense = *m.getIdentityPtr();
      break;
    default:
      SiconosMatrixException::selfThrow("SimpleMatrix::op= (const SimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 2:
    switch (numM)
    {
    case 2:
      *mat.Triang = *m.getTriangPtr();
      break;
    case 6:
      *mat.Triang = *m.getZeroPtr();
      break;
    case 7:
      *mat.Triang = *m.getIdentityPtr();
      break;
    default:
      SiconosMatrixException::selfThrow("SimpleMatrix::assignment of a bad type of matrix into a triangular one.");
      break;
    }
    break;
  case 3:
    if (numM == 3)
      *mat.Sym = *m.getSymPtr();
    else if (numM == 6)
      *mat.Sym = *m.getZeroPtr();
    else if (numM == 7)
      *mat.Sym = *m.getIdentityPtr();
    else
      SiconosMatrixException::selfThrow("SimpleMatrix::bad assignment of matrix (symetric one = dense or ...)");
    break;
  case 4:
    switch (numM)
    {
    case 2:
      *mat.Sparse = *m.getTriangPtr();
      break;
    case 3:
      *mat.Sparse = *m.getSymPtr();
      break;
    case 4:
      *mat.Sparse = *m.getSparsePtr();
      break;
    case 5:
      *mat.Sparse = *m.getBandedPtr();
      break;
    case 6:
      *mat.Sparse = *m.getZeroPtr();
      break;
    case 7:
      *mat.Sparse = *m.getIdentityPtr();
      break;
    default:
      SiconosMatrixException::selfThrow("SimpleMatrix::op= (const SimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 5:
    switch (numM)
    {
    case 5:
      *mat.Banded = *m.getBandedPtr();
      break;
    case 6:
      *mat.Banded = *m.getZeroPtr();
      break;
    case 7:
      *mat.Banded = *m.getIdentityPtr();
      break;
    default:
      SiconosMatrixException::selfThrow("SimpleMatrix::op= (const SimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  default:
    SiconosMatrixException::selfThrow("SimpleMatrix::op= (const SimpleMatrix) : invalid type of matrix");
    break;
  }
  resetLU();
  return *this;
}

SimpleMatrix& SimpleMatrix::operator = (const SiconosMatrix& m)
{

  // Warning!!! If sizes are inconsistent between m and this, boost operator = results in resize of this to the dim of m !!!
  // Add an exception to prevent this???
  // if exception, strange behavio in things like A = A*B ...
  //if(dim[0]!=m.size(0) || dim[1] !=m.size(1))
  //    SiconosMatrixException::selfThrow("SimpleMatrix operator = failed. Inconsistent sizes.");
  unsigned int numM = m.getNum();
  switch (num)
  {
  case 1:
    switch (numM)
    {
    case 1:
      *mat.Dense = *m.getDensePtr();
      break;
    case 2:
      *mat.Dense = *m.getTriangPtr();
      break;
    case 3:
      *mat.Dense = *m.getSymPtr();
      break;
    case 4:
      *mat.Dense = *m.getSparsePtr();
      break;
    case 5:
      *mat.Dense = *m.getBandedPtr();
      break;
    case 6:
      *mat.Dense = *m.getZeroPtr(); // warning: this is not equivalent to a zero() call, because boost = results in resizing if required.
      break;
    case 7:
      *mat.Dense = *m.getIdentityPtr();
      break;
    default:
      SiconosMatrixException::selfThrow("SimpleMatrix::op= (const SiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 2:
    switch (numM)
    {
    case 2:
      *mat.Triang = *m.getTriangPtr();
      break;
    case 6:
      *mat.Triang = *m.getZeroPtr();
      break;
    case 7:
      *mat.Triang = *m.getIdentityPtr();
      break;
    default:
      SiconosMatrixException::selfThrow("SimpleMatrix::op= (const SiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 3:
    switch (numM)
    {
    case 3:
      *mat.Sym = *m.getSymPtr();
      break;
    case 6:
      *mat.Sym = *m.getZeroPtr();
      break;
    case 7:
      *mat.Sym = *m.getIdentityPtr();
      break;
    default:
      SiconosMatrixException::selfThrow("SimpleMatrix::op= (const SiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 4:
    switch (numM)
    {
    case 2:
      *mat.Sparse = *m.getTriangPtr();
      break;
    case 3:
      *mat.Sparse = *m.getSymPtr();
      break;
    case 4:
      *mat.Sparse = *m.getSparsePtr();
      break;
    case 5:
      *mat.Sparse = *m.getBandedPtr();
      break;
    case 6:
      *mat.Sparse = *m.getZeroPtr();
      break;
    case 7:
      *mat.Sparse = *m.getIdentityPtr();
      break;
    default:
      SiconosMatrixException::selfThrow("SimpleMatrix::op= (const SimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 5:
    switch (numM)
    {
    case 5:
      *mat.Banded = *m.getBandedPtr();
      break;
    case 6:
      *mat.Banded = *m.getZeroPtr();
      break;
    case 7:
      *mat.Banded = *m.getIdentityPtr();
      break;
    default:
      SiconosMatrixException::selfThrow("SimpleMatrix::op= (const SimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  default:
    SiconosMatrixException::selfThrow("SimpleMatrix::op= (const SiconosMatrix) : invalid type of matrix");
    break;
  }
  resetLU();
  return *this;
}

SimpleMatrix& SimpleMatrix::operator += (const SiconosMatrix& m)
{
  unsigned int numM = m.getNum();
  switch (num)
  {
  case 1:
    switch (numM)
    {
    case 1:
      *mat.Dense += *m.getDensePtr();
      break;
    case 2:
      *mat.Dense += *m.getTriangPtr();
      break;
    case 3:
      *mat.Dense += *m.getSymPtr();
      break;
    case 4:
      *mat.Dense += *m.getSparsePtr();
      break;
    case 5:
      *mat.Dense += *m.getBandedPtr();
      break;
    case 6:
      break;
    case 7:
      *mat.Dense += *m.getIdentityPtr();
      break;
    default:
      SiconosMatrixException::selfThrow("op+= (const SiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 2:
    switch (numM)
    {
    case 2:
      *mat.Triang += *m.getTriangPtr();
      break;
    case 6:
      break;
    case 7:
      *mat.Triang += *m.getIdentityPtr();
      break;
    default:
      SiconosMatrixException::selfThrow("op+= (const SiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 3:
    switch (numM)
    {
    case 3:
      *mat.Sym += *m.getSymPtr();
      break;
    case 6:
      break;
    case 7:
      *mat.Sym += *m.getIdentityPtr();
      break;
    default:
      SiconosMatrixException::selfThrow("op+= (const SiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 4:
    switch (numM)
    {
    case 2:
      *mat.Sparse += *m.getTriangPtr();
      break;
    case 3:
      *mat.Sparse += *m.getSymPtr();
      break;
    case 4:
      *mat.Sparse += *m.getSparsePtr();
      break;
    case 5:
      *mat.Sparse += *m.getBandedPtr();
      break;
    case 6:
      break;
    default:
      SiconosMatrixException::selfThrow("op= (const SimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 5:
    switch (numM)
    {
    case 5:
      *mat.Banded += *m.getBandedPtr();
      break;
    case 6:
      break;
    case 7:
      *mat.Banded += *m.getIdentityPtr();
      break;
    default:
      SiconosMatrixException::selfThrow("op= (const SimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  default:
    SiconosMatrixException::selfThrow("op+= (const SiconosMatrix) : invalid type of matrix");
    break;
  }
  resetLU();
  return *this;
}

SimpleMatrix& SimpleMatrix::operator -= (const SiconosMatrix& m)
{

  unsigned int numM = m.getNum();

  switch (num)
  {
  case 1:
    switch (numM)
    {
    case 1:
      *mat.Dense -= *m.getDensePtr();
      break;
    case 2:
      *mat.Dense -= *m.getTriangPtr();
      break;
    case 3:
      *mat.Dense -= *m.getSymPtr();
      break;
    case 4:
      *mat.Dense -= *m.getSparsePtr();
      break;
    case 5:
      *mat.Dense -= *m.getBandedPtr();
      break;
    case 6:
      break;
    case 7:
      *mat.Dense -= *m.getIdentityPtr();
      break;
    default:
      SiconosMatrixException::selfThrow("op-= (const SiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 2:
    switch (numM)
    {
    case 2:
      *mat.Triang -= *m.getTriangPtr();
      break;
    case 6:
      break;
    case 7:
      *mat.Triang -= *m.getIdentityPtr();
      break;
    default:
      SiconosMatrixException::selfThrow("op-= (const SiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 3:
    switch (numM)
    {
    case 3:
      *mat.Sym -= *m.getSymPtr();
      break;
    case 6:
      break;
    case 7:
      *mat.Sym -= *m.getIdentityPtr();
      break;
    default:
      SiconosMatrixException::selfThrow("op-= (const SiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 4:
    switch (numM)
    {
    case 2:
      *mat.Sparse -= *m.getTriangPtr();
      break;
    case 3:
      *mat.Sparse -= *m.getSymPtr();
      break;
    case 4:
      *mat.Sparse -= *m.getSparsePtr();
      break;
    case 5:
      *mat.Sparse -= *m.getBandedPtr();
      break;
    case 6:
      break;
    default:
      SiconosMatrixException::selfThrow("op-= (const SimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 5:
    switch (numM)
    {
    case 5:
      *mat.Banded -= *m.getBandedPtr();
      break;
    case 6:
      break;
    case 7:
      *mat.Banded -= *m.getIdentityPtr();
      break;
    default:
      SiconosMatrixException::selfThrow("op-= (const SimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  default:
    SiconosMatrixException::selfThrow("op-= (const SiconosMatrix) : invalid type of matrix");
    break;
  }
  resetLU();
  return *this;
}

SimpleMatrix& SimpleMatrix::operator *= (double m)
{
  switch (num)
  {
  case 1:
    *mat.Dense *= m;
    break;
  case 2:
    *mat.Triang *= m;
    break;
  case 3:
    *mat.Sym *= m;
    break;
  case 4:
    *mat.Sparse *= m;
    break;
  case 5:
    *mat.Banded *= m;
    break;
  case 6:
    break;
  default:
    SiconosMatrixException::selfThrow("op*= (double) : invalid type of matrix");
    break;
  }
  return *this;
}

SimpleMatrix& SimpleMatrix::operator *= (int m)
{
  switch (num)
  {
  case 1:
    *mat.Dense *= m;
    break;
  case 2:
    *mat.Triang *= m;
    break;
  case 3:
    *mat.Sym *= m;
    break;
  case 4:
    *mat.Sparse *= m;
    break;
  case 5:
    *mat.Banded *= m;
    break;
  case 6:
    break;
  default:
    SiconosMatrixException::selfThrow("op*= (int) : invalid type of matrix");
    break;
  }
  return *this;
}

SimpleMatrix& SimpleMatrix::operator /= (double m)
{
  if (m == 0)
    SiconosMatrixException::selfThrow("op/= (double) : division by zero.");

  switch (num)
  {
  case 1:
    *mat.Dense /= m;
    break;
  case 2:
    *mat.Triang /= m;
    break;
  case 3:
    *mat.Sym /= m;
    break;
  case 4:
    *mat.Sparse /= m;
    break;
  case 5:
    *mat.Banded /= m;
    break;
  case 6:
    break;
  default:
    SiconosMatrixException::selfThrow("op/= (double) : invalid type of matrix");
    break;
  }
  return *this;
}

SimpleMatrix& SimpleMatrix::operator /= (int m)
{
  if (m == 0)
    SiconosMatrixException::selfThrow("op/= (int) : division by zero.");
  switch (num)
  {
  case 1:
    *mat.Dense /= m;
    break;
  case 2:
    *mat.Triang /= m;
    break;
  case 3:
    *mat.Sym /= m;
    break;
  case 4:
    *mat.Sparse /= m;
    break;
  case 5:
    *mat.Banded /= m;
    break;
  case 6:
    break;
  default:
    SiconosMatrixException::selfThrow("op/= (int) : invalid type of matrix");
    break;
  }
  return *this;
}

void SimpleMatrix::PLUFactorizationInPlace()
{
  if (num != 1)
    SiconosMatrixException::selfThrow(" SimpleMatrix::PLUFactorizationInPlace : only implemented for dense matrices.");

  ipiv.resize(size(0));
  int info = boost::numeric::bindings::atlas::getrf(*mat.Dense, ipiv);
  if (info != 0)
    std::cout << "SimpleMatrix::PLUFactorizationInPlace warning: the matrix is singular." << std::endl;
  isPLUFactorized = true;
}

void SimpleMatrix::PLUInverseInPlace()
{
  if (!isPLUFactorized)
    PLUFactorizationInPlace();

  int info = boost::numeric::bindings::atlas::getri(*mat.Dense, ipiv);   // solve from factorization
  if (info != 0)
    SiconosMatrixException::selfThrow("SimpleMatrix::PLUInverseInPlace failed, the matrix is singular.");

  isPLUInversed = true;
}

void SimpleMatrix::PLUForwardBackwardInPlace(SiconosMatrix &B)
{
  int info;
  if (!isPLUFactorized) // call gesv => LU-factorize+solve
  {
    // solve system:
    ipiv.resize(size(0));
    info = boost::numeric::bindings::atlas::gesv(*mat.Dense, ipiv, *(B.getDensePtr()));
    isPLUFactorized = true;
    // B now contains solution:
  }
  else // call getrs: only solve using previous lu-factorization
    info = boost::numeric::bindings::atlas::getrs(*mat.Dense, ipiv, *(B.getDensePtr()));

  if (info != 0)
    SiconosMatrixException::selfThrow("SimpleMatrix::PLUForwardBackwardInPlace failed.");
}

void SimpleMatrix::PLUForwardBackwardInPlace(SiconosVector &B)
{
  DenseMat tmpB(B.size(), 1);
  ublas::column(tmpB, 0) = *(B.getDensePtr()); // Conversion of vector to matrix. Temporary solution.
  int info;
  if (!isPLUFactorized) // call gesv => LU-factorize+solve
  {
    // solve system:
    ipiv.resize(size(0));
    info = boost::numeric::bindings::atlas::gesv(*mat.Dense, ipiv, tmpB);
    isPLUFactorized = true;
    // B now contains solution:
  }
  else // call getrs: only solve using previous lu-factorization
    info = boost::numeric::bindings::atlas::getrs(*mat.Dense, ipiv, tmpB);

  if (info != 0)
    SiconosMatrixException::selfThrow("SimpleMatrix::PLUForwardBackwardInPlace failed.");
  *(B.getDensePtr()) = ublas::column(tmpB, 0);
}

void SimpleMatrix::resetLU()
{
  ipiv.clear();
  isPLUFactorized = false;
  isPLUInversed = false;
}

bool operator == (const SiconosMatrix &m, const SiconosMatrix &x)
{
  if (m.isBlock() ||  x.isBlock())
    SiconosMatrixException::selfThrow("op == (const SiconosMatrix, const SiconosMatrix) : incompatible type of matrix");
  double norm = (sub(m, x)).normInf();
  return (norm < tolerance);
}

SimpleMatrix operator + (const SiconosMatrix &x, const SiconosMatrix &m)
{
  if ((x.size(0) != m.size(0)) || (x.size(1) != m.size(1)))
    SiconosMatrixException::selfThrow("Matrix addition: inconsistent sizes");

  unsigned int numX = x.getNum();
  if (numX != m.getNum())
    SiconosMatrixException::selfThrow("SimpleMatrix:Matrix addition: use function add in order to add matrices of different type");

  if (numX == 7)
    SiconosMatrixException::selfThrow("SimpleMatrix:Matrix addition: illegal matrix type for operator -. Try sub()?");

  if (numX == 1)
  {
    DenseMat p = *x.getDensePtr() + *m.getDensePtr();
    return p;
  }
  else if (numX == 2)
  {
    TriangMat t = *x.getTriangPtr() + *m.getTriangPtr();
    return t;
  }
  else if (numX == 3)
  {
    SymMat s = *x.getSymPtr() + *m.getSymPtr();
    return s;
  }
  else if (numX == 4)
  {
    SparseMat sp = *x.getSparsePtr() + *m.getSparsePtr();
    return sp;
  }
  else if (numX == 5)
  {
    BandedMat b;
    b.resize(m.size(0), m.size(1), (*m.getBandedPtr()).lower(), (*m.getBandedPtr()).upper(), false);
    b = *x.getBandedPtr() + *m.getBandedPtr();
    return b;
  }
  else //if(numX == 6){
  {
    ZeroMat z;
    return z;
  }
}

SimpleMatrix operator - (const SiconosMatrix &x, const SiconosMatrix &m)
{
  if ((x.size(0) != m.size(0)) || (x.size(1) != m.size(1)))
    SiconosMatrixException::selfThrow("Matrix subtraction: inconsistent sizes");

  unsigned int numX = x.getNum();
  if (numX != m.getNum())
    SiconosMatrixException::selfThrow("SimpleMatrix:Matrix subtraction: use function sub in order to subtract matrices of different type");

  if (numX == 7)
    SiconosMatrixException::selfThrow("SimpleMatrix:Matrix subtraction: illegal matrix type for operator -. Try sub()?");

  if (numX == 1)
  {
    DenseMat p = *x.getDensePtr() - *m.getDensePtr();
    return p;
  }
  else if (numX == 2)
  {
    TriangMat t = *x.getTriangPtr() - *m.getTriangPtr();
    return t;
  }
  else if (numX == 3)
  {
    SymMat s = *x.getSymPtr() - *m.getSymPtr();
    return s;
  }
  else if (numX == 4)
  {
    SparseMat sp = *x.getSparsePtr() - *m.getSparsePtr();
    return sp;
  }
  else if (numX == 5)
  {
    BandedMat b;
    b.resize(m.size(0), m.size(1), (*m.getBandedPtr()).lower(), (*m.getBandedPtr()).upper(), false);
    b = *x.getBandedPtr() - *m.getBandedPtr();
    return b;
  }
  else  // if(numX == 6){
  {
    ZeroMat z;
    return z;
  }
}

SimpleMatrix operator * (const SiconosMatrix &x, const SiconosMatrix &m)
{

  if ((x.size(1) != m.size(0)))
    SiconosMatrixException::selfThrow("Matrix product : inconsistent sizes");

  unsigned int numX = x.getNum();
  if (numX != m.getNum())
    SiconosMatrixException::selfThrow("SimpleMatrix:Matrix product : use function prod in order to multiply matrices of different type");

  if (numX == 1)
  {
    DenseMat p = prod(*x.getDensePtr(), *m.getDensePtr());
    return p;
  }
  else if (numX == 2)
  {
    TriangMat t = prod(*x.getTriangPtr(), *m.getTriangPtr());
    return t;
  }
  else if (numX == 3)
  {
    SymMat s = prod(*x.getSymPtr(), *m.getSymPtr());
    return s;
  }
  else if (numX == 4)
  {
    SparseMat sp = prod(*x.getSparsePtr(), *m.getSparsePtr());
    return sp;
  }
  else if (numX == 5)
  {
    DenseMat p = prod(*x.getBandedPtr(), *m.getBandedPtr());
    return p;
  }
  else if (numX == 6)
  {
    ZeroMat z;
    return z;
  }
  else  // num = 7
  {
    IdentityMat I;
    return I;
  }
}

SimpleMatrix add(const SiconosMatrix &x, const SiconosMatrix& m)
{
  DenseMat p;

  if ((x.size(0) != m.size(0)) || (x.size(1) != m.size(1)))
    SiconosMatrixException::selfThrow("SimpleMatrix: function add, inconsistent sizes.");

  unsigned int numM = m.getNum();
  unsigned int numX = x.getNum();

  if (numM == 1)
  {
    if (numX == 1)
      p = *x.getDensePtr() + *m.getDensePtr();
    else if (numX == 2)
      p = *x.getTriangPtr() + *m.getDensePtr();
    else if (numX == 3)
      p = *x.getSymPtr() + *m.getDensePtr();
    else if (numX == 4)
      p = *x.getSparsePtr() + *m.getDensePtr();
    else if (numX == 5)
      p = *x.getBandedPtr() + *m.getDensePtr();
    else if (numX == 6)
      p = *m.getDensePtr();
    else if (numX == 7)
      p = *x.getIdentityPtr() + *m.getDensePtr();
  }
  else if (numM == 2)
  {
    if (numX == 1)
      p = *x.getDensePtr() + *m.getTriangPtr();
    else if (numX == 2)
    {
      TriangMat t = *x.getTriangPtr() + *m.getTriangPtr();
      return t;
    }
    else if (numX == 3)
      p = *x.getSymPtr() + *m.getTriangPtr();
    else if (numX == 4)
      p = *x.getSparsePtr() + *m.getTriangPtr();
    else if (numX == 5)
      p = *x.getBandedPtr() + *m.getTriangPtr();
    else if (numX == 6)
      p = *m.getTriangPtr();
    else if (numX == 7)
      p = *x.getIdentityPtr() + *m.getTriangPtr();
  }
  else if (numM == 3)
  {
    if (numX == 1)
      p = *x.getDensePtr() + *m.getSymPtr();
    else if (numX == 2)
      p = *x.getTriangPtr() + *m.getSymPtr();
    else if (numX == 3)
    {
      SymMat s = *x.getSymPtr() + *m.getSymPtr();
      return s;
    }
    else if (numX == 4)
      p = *x.getSparsePtr() + *m.getSymPtr();
    else if (numX == 5)
      p = *x.getBandedPtr() + *m.getSymPtr();
    else if (numX == 6)
      p = *m.getSymPtr();
    else if (numX == 7)
      p = *x.getIdentityPtr() + *m.getSymPtr();
  }
  else if (numM == 4)
  {
    if (numX == 1)
      p = *x.getDensePtr() + *m.getSparsePtr();
    else if (numX == 2)
      p = *x.getTriangPtr() + *m.getSparsePtr();
    else if (numX == 3)
      p = *x.getSymPtr() + *m.getSparsePtr();
    else if (numX == 4)
      p = *x.getSparsePtr() + *m.getSparsePtr();
    else if (numX == 5)
      p = *x.getBandedPtr() + *m.getSparsePtr();
    else if (numX == 6)
      p = *m.getSparsePtr();
    else if (numX == 7)
      p = *x.getIdentityPtr() + *m.getSparsePtr();
  }
  else if (numM == 5)
  {
    if (numX == 1)
      p = *x.getDensePtr() + *m.getBandedPtr();
    else if (numX == 2)
      p = *x.getTriangPtr() + *m.getBandedPtr();
    else if (numX == 3)
      p = *x.getSymPtr() + *m.getBandedPtr();
    else if (numX == 4)
      p = *x.getSparsePtr() + *m.getBandedPtr();
    else if (numX == 5)
    {
      BandedMat b = *x.getBandedPtr() + *m.getBandedPtr();
      return b;
    }
    else if (numX == 6)
      p = *m.getBandedPtr();
    else if (numX == 7)
      p = *x.getIdentityPtr() + *m.getBandedPtr();
  }
  else if (numM == 6)
  {
    if (numX == 1)
      p = *x.getDensePtr();
    else if (numX == 2)
    {
      TriangMat t  = *x.getTriangPtr();
      return t;
    }
    else if (numX == 3)
    {
      SymMat s = *x.getSymPtr();
      return s;
    }
    else if (numX == 4)
    {
      SparseMat s = *x.getSparsePtr();
      return s;
    }
    else if (numX == 5)
    {
      BandedMat b = *x.getBandedPtr();
      return b;
    }
    else if (numX == 6)
    {
      ZeroMat z(x.size(0), x.size(1));
      return z;
    }
    else if (numX == 7)
    {
      IdentityMat i = *x.getIdentityPtr();
      return i;
    }
  }
  else if (numM == 7)
  {
    if (numX == 1)
      p = *x.getDensePtr() + *m.getIdentityPtr();
    else if (numX == 2)
    {
      TriangMat t = *x.getTriangPtr() + *m.getIdentityPtr();
      return t;
    }
    else if (numX == 3)
    {
      SymMat s = *x.getSymPtr() + *m.getIdentityPtr();
      return s;
    }
    else if (numX == 4)
    {
      SparseMat s = *x.getSparsePtr() + *m.getIdentityPtr();
      return s;
    }
    else if (numX == 5)
      p = *x.getBandedPtr() + *m.getIdentityPtr();
    else if (numX == 6)
      p = *m.getIdentityPtr();
    else if (numX == 7)
      p = 2 * *x.getIdentityPtr();
  }
  else
  {
    SiconosMatrixException::selfThrow("Matrix function sub: invalid type of matrix");
  }
  return p;
}

SimpleMatrix sub(const SiconosMatrix &x, const SiconosMatrix& m)
{
  DenseMat p;

  if ((x.size(0) != m.size(0)) || (x.size(1) != m.size(1)))
    SiconosMatrixException::selfThrow("SimpleMatrix: function sub, inconsistent sizes.");

  unsigned int numM = m.getNum();
  unsigned int numX = x.getNum();

  if (numM == 1)
  {
    if (numX == 1)
      p = *x.getDensePtr() - *m.getDensePtr();
    else if (numX == 2)
      p = *x.getTriangPtr() - *m.getDensePtr();
    else if (numX == 3)
      p = *x.getSymPtr() - *m.getDensePtr();
    else if (numX == 4)
      p = *x.getSparsePtr() - *m.getDensePtr();
    else if (numX == 5)
      p = *x.getBandedPtr() - *m.getDensePtr();
    else if (numX == 6)
      p = -*m.getDensePtr();
    else if (numX == 7)
      p = *x.getIdentityPtr() - *m.getDensePtr();
  }
  else if (numM == 2)
  {
    if (numX == 1)
      p = *x.getDensePtr() - *m.getTriangPtr();
    else if (numX == 2)
    {
      TriangMat t = *x.getTriangPtr() - *m.getTriangPtr();
      return t;
    }
    else if (numX == 3)
      p = *x.getSymPtr() - *m.getTriangPtr();
    else if (numX == 4)
      p = *x.getSparsePtr() - *m.getTriangPtr();
    else if (numX == 5)
      p = *x.getBandedPtr() - *m.getTriangPtr();
    else if (numX == 6)
      p = -*m.getTriangPtr();
    else if (numX == 7)
      p = *x.getIdentityPtr() - *m.getTriangPtr();
  }
  else if (numM == 3)
  {
    if (numX == 1)
      p = *x.getDensePtr() - *m.getSymPtr();
    else if (numX == 2)
      p = *x.getTriangPtr() - *m.getSymPtr();
    else if (numX == 3)
    {
      SymMat s = *x.getSymPtr() - *m.getSymPtr();
      return s;
    }
    else if (numX == 4)
      p = *x.getSparsePtr() - *m.getSymPtr();
    else if (numX == 5)
      p = *x.getBandedPtr() - *m.getSymPtr();
    else if (numX == 6)
      p = -*m.getSymPtr();
    else if (numX == 7)
      p = *x.getIdentityPtr() - *m.getSymPtr();
  }
  else if (numM == 4)
  {
    if (numX == 1)
      p = *x.getDensePtr() - *m.getSparsePtr();
    else if (numX == 2)
      p = *x.getTriangPtr() - *m.getSparsePtr();
    else if (numX == 3)
      p = *x.getSymPtr() - *m.getSparsePtr();
    else if (numX == 4)
      p = *x.getSparsePtr() - *m.getSparsePtr();
    else if (numX == 5)
      p = *x.getBandedPtr() - *m.getSparsePtr();
    else if (numX == 6)
      p = -*m.getSparsePtr();
    else if (numX == 7)
      p = *x.getIdentityPtr() - *m.getSparsePtr();
  }
  else if (numM == 5)
  {
    if (numX == 1)
      p = *x.getDensePtr() - *m.getBandedPtr();
    else if (numX == 2)
      p = *x.getTriangPtr() - *m.getBandedPtr();
    else if (numX == 3)
      p = *x.getSymPtr() - *m.getBandedPtr();
    else if (numX == 4)
      p = *x.getSparsePtr() - *m.getBandedPtr();
    else if (numX == 5)
    {
      BandedMat b = *x.getBandedPtr() - *m.getBandedPtr();
      return b;
    }
    else if (numX == 6)
      p = -*m.getBandedPtr();
    else if (numX == 7)
      p = *x.getIdentityPtr() - *m.getBandedPtr();
  }
  else if (numM == 6)
  {
    if (numX == 1)
      p = *x.getDensePtr();
    else if (numX == 2)
    {
      TriangMat t  = *x.getTriangPtr();
      return t;
    }
    else if (numX == 3)
    {
      SymMat s = *x.getSymPtr();
      return s;
    }
    else if (numX == 4)
    {
      SparseMat s = *x.getSparsePtr();
      return s;
    }
    else if (numX == 5)
    {
      BandedMat b = *x.getBandedPtr();
      return b;
    }
    else if (numX == 6)
    {
      ZeroMat z(x.size(0), x.size(1));
      return z;
    }
    else if (numX == 7)
    {
      IdentityMat i = *x.getIdentityPtr();
      return i;
    }
  }
  else if (numM == 7)
  {
    if (numX == 1)
      p = *x.getDensePtr() - *m.getIdentityPtr();
    else if (numX == 2)
    {
      TriangMat t = *x.getTriangPtr() - *m.getIdentityPtr();
      return t;
    }
    else if (numX == 3)
    {
      SymMat s = *x.getSymPtr() - *m.getIdentityPtr();
      return s;
    }
    else if (numX == 4)
    {
      SparseMat s = *x.getSparsePtr() - *m.getIdentityPtr();
      return s;
    }
    else if (numX == 5)
      p = *x.getBandedPtr() - *m.getIdentityPtr();

    else if (numX == 6)
      p = -*m.getIdentityPtr();
    else if (numX == 7)
    {
      ZeroMat z(x.size(0), x.size(1));
      return z;
    }
  }
  else
  {
    SiconosMatrixException::selfThrow("Matrix function sub: invalid type of matrix");
  }
  return p;
}

SimpleMatrix prod(const SiconosMatrix &x, const SiconosMatrix& m)
{
  if ((x.size(1) != m.size(0)))
    SiconosMatrixException::selfThrow("Matrix function prod : inconsistent sizes");

  if (x.isBlock() || m.isBlock())
    SiconosMatrixException::selfThrow("Matrix function prod : not yet implemented for block matrices.");

  unsigned int numM = m.getNum();
  unsigned int numX = x.getNum();

  if (numM == 6 || numX == 6)
  {
    DenseMat p(ublas::scalar_matrix<double>(x.size(0), m.size(1), 0.0));
    return p;
  }
  else
  {
    DenseMat p(x.size(0), m.size(1));
    if (numM == 1)
    {
      if (numX == 1)
        p = prod(*x.getDensePtr(), *m.getDensePtr());

      else if (numX == 2)
        p = prod(*x.getTriangPtr(), *m.getDensePtr());

      else if (numX == 3)
        p = prod(*x.getSymPtr(), *m.getDensePtr());

      else if (numX == 4)
        p = prod(*x.getSparsePtr(), *m.getDensePtr());

      else if (numX == 5)
        p = prod(*x.getBandedPtr(), *m.getDensePtr());

      else //if(numX==7)
        p = prod(*x.getIdentityPtr(), *m.getDensePtr());

      return p;
    }
    else if (numM == 2)
    {
      if (numX == 1)
      {
        p = prod(*x.getDensePtr(), *m.getTriangPtr());
        return p;
      }
      else if (numX == 2)
      {
        TriangMat t = prod(*x.getTriangPtr(), *m.getTriangPtr());
        return t;
      }
      else if (numX == 3)
      {
        p = prod(*x.getSymPtr(), *m.getTriangPtr());
        return p;
      }
      else if (numX == 4)
      {
        p = prod(*x.getSparsePtr(), *m.getTriangPtr());
        return p;
      }
      else if (numX == 5)
      {
        p = prod(*x.getBandedPtr(), *m.getTriangPtr());
        return p;
      }
      else //if(numX==7)
      {
        p = prod(*x.getIdentityPtr(), *m.getTriangPtr());
        return p;
      }
    }
    else if (numM == 3)
    {
      if (numX == 1)
      {
        p = prod(*x.getDensePtr(), *m.getSymPtr());
        return p;
      }
      else if (numX == 2)
      {
        p = prod(*x.getTriangPtr(), *m.getSymPtr());
        return p;
      }
      else if (numX == 3)
      {
        SymMat s = prod(*x.getSymPtr(), *m.getSymPtr());
        return s;
      }
      else if (numX == 4)
      {
        p = prod(*x.getSparsePtr(), *m.getSymPtr());
        return p;
      }
      else if (numX == 5)
      {
        p = prod(*x.getBandedPtr(), *m.getSymPtr());
        return p;
      }
      else //if(numX==7)
      {
        p = prod(*x.getIdentityPtr(), *m.getSymPtr());
        return p;
      }
    }
    else if (numM == 4)
    {
      if (numX == 1)
      {
        p = prod(*x.getDensePtr(), *m.getSparsePtr());
        return p;
      }
      else if (numX == 2)
      {
        p = prod(*x.getTriangPtr(), *m.getSparsePtr());
        return p;
      }
      else if (numX == 3)
      {
        p = prod(*x.getSymPtr(), *m.getSparsePtr());
        return p;
      }
      else if (numX == 4)
      {
        SparseMat sp = prod(*x.getSparsePtr(), *m.getSparsePtr());
        return sp;
      }
      else if (numX == 5)
      {
        p = prod(*x.getBandedPtr(), *m.getSparsePtr());
        return p;
      }
      else //if(numX==7)
      {
        p = prod(*x.getIdentityPtr(), *m.getSparsePtr());
        return p;
      }
    }
    else if (numM == 5)
    {
      if (numX == 1)
        p = prod(*x.getDensePtr(), *m.getBandedPtr());

      else if (numX == 2)
        p = prod(*x.getTriangPtr(), *m.getBandedPtr());

      else if (numX == 3)
        p = prod(*x.getSymPtr(), *m.getBandedPtr());

      else if (numX == 4)
        p = prod(*x.getSparsePtr(), *m.getBandedPtr());

      else if (numX == 5)
        p = prod(*x.getBandedPtr(), *m.getBandedPtr());

      else //if(numX==7)
        p = prod(*x.getIdentityPtr(), *m.getBandedPtr());

      return p;
    }
    else //if(numM==7)
    {
      if (numX == 1)
        p = prod(*x.getDensePtr(), *m.getIdentityPtr());

      else if (numX == 2)
        p = prod(*x.getTriangPtr(), *m.getIdentityPtr());

      else if (numX == 3)
        p = prod(*x.getSymPtr(), *m.getIdentityPtr());

      else if (numX == 4)
        p = prod(*x.getSparsePtr(), *m.getIdentityPtr());

      else if (numX == 5)
        p = prod(*x.getBandedPtr(), *m.getIdentityPtr());

      else //if(numX==7)
        p = prod(*x.getIdentityPtr(), *m.getIdentityPtr());
      return p;
    }
  }
}

SimpleMatrix operator * (const SiconosMatrix &m, double d)
{
  unsigned int num = m.getNum();

  if (num == 7)
    SiconosMatrixException::selfThrow("SimpleMatrix:operator * (m*double), forbidden for this type of matrix.");

  if (num == 1)
  {
    DenseMat p;
    p = *m.getDensePtr() * d;
    return p;
  }
  else if (num == 2)
  {
    TriangMat t;
    t = *m.getTriangPtr() * d;
    return t;
  }
  else if (num == 3)
  {
    SymMat s;
    s = *m.getSymPtr() * d;
    return s;
  }
  else if (num == 4)
  {
    SparseMat sp;
    sp = *m.getSparsePtr() * d;
    return sp;
  }
  else if (num == 5)
  {
    BandedMat b;
    b.resize(m.size(0), m.size(1), (*m.getBandedPtr()).lower(), (*m.getBandedPtr()).upper(), false);
    b = *m.getBandedPtr() * d;
    return b;
  }
  else //if(num==6)
  {
    ZeroMat z;
    return z;
  }
}

SimpleMatrix operator * (const SiconosMatrix &m, int d)
{

  unsigned int num = m.getNum();
  if (num == 7)
    SiconosMatrixException::selfThrow("SimpleMatrix:operator * (m*int), forbidden for this type of matrix.");

  if (num == 1)
  {
    DenseMat p = *m.getDensePtr() * d;
    return p;
  }
  else if (num == 2)
  {
    TriangMat t = *m.getTriangPtr() * d;
    return t;
  }
  else if (num == 3)
  {
    SymMat s = *m.getSymPtr() * d;
    return s;
  }
  else if (num == 4)
  {
    SparseMat sp = *m.getSparsePtr() * d;
    return sp;
  }
  else if (num == 5)
  {
    BandedMat b;
    b.resize(m.size(0), m.size(1), (*m.getBandedPtr()).lower(), (*m.getBandedPtr()).upper(), false);
    b = *m.getBandedPtr() * d;
    return b;
  }
  else // if(num==6)
  {
    ZeroMat z;
    return z;
  }
}

SimpleMatrix operator * (double d, const SiconosMatrix &m)
{
  unsigned int num = m.getNum();
  if (num == 7)
    SiconosMatrixException::selfThrow("SimpleMatrix:operator * (double*m), forbidden for this type of matrix.");

  if (num == 1)
  {
    DenseMat p = d * *m.getDensePtr();
    return p;
  }
  else if (num == 2)
  {
    TriangMat t = d * *m.getTriangPtr();
    return t;
  }
  else if (num == 3)
  {
    SymMat s = d * *m.getSymPtr();
    return s;
  }
  else if (num == 4)
  {
    SparseMat sp = d * *m.getSparsePtr();
    return sp;
  }
  else if (num == 5)
  {
    BandedMat b;
    b.resize(m.size(0), m.size(1), (*m.getBandedPtr()).lower(), (*m.getBandedPtr()).upper(), false);
    b = d * *m.getBandedPtr();
    return b;
  }
  else //if(num==6)
  {
    ZeroMat z;
    return z;
  }
}

SimpleMatrix operator * (int d, const SiconosMatrix &m)
{
  unsigned int num = m.getNum();
  if (num == 7)
    SiconosMatrixException::selfThrow("SimpleMatrix:operator * (int*m), forbidden for this type of matrix.");

  if (num == 1)
  {
    DenseMat p = d * *m.getDensePtr();
    return p;
  }
  else if (num == 2)
  {
    TriangMat t = d * *m.getTriangPtr();
    return t;
  }
  else if (num == 3)
  {
    SymMat s = d * *m.getSymPtr();
    return s;
  }
  else if (num == 4)
  {
    SparseMat sp = d * *m.getSparsePtr();
    return sp;
  }
  else if (num == 5)
  {
    BandedMat b;
    b.resize(m.size(0), m.size(1), (*m.getBandedPtr()).lower(), (*m.getBandedPtr()).upper(), false);
    b = d * *m.getBandedPtr();
    return b;
  }
  else //if(num==6)
  {
    ZeroMat z;
    return z;
  }
}

SimpleMatrix operator / (const SiconosMatrix &m, double d)
{
  if (d == 0)
    SiconosMatrixException::selfThrow("SimpleMatrix:operator /, division by zero.");

  unsigned int num = m.getNum();
  if (num == 7)
    SiconosMatrixException::selfThrow("SimpleMatrix:operator / (m/double), forbidden for this type of matrix.");

  if (num == 1)
  {
    DenseMat p = *m.getDensePtr() / d;
    return p;
  }
  else if (num == 2)
  {
    TriangMat t = *m.getTriangPtr() / d;
    return t;
  }
  else if (num == 3)
  {
    SymMat s = *m.getSymPtr() / d;
    return s;
  }
  else if (num == 4)
  {
    SparseMat sp = *m.getSparsePtr() / d;
    return sp;
  }
  else if (num == 5)
  {
    BandedMat b;
    b.resize(m.size(0), m.size(1), (*m.getBandedPtr()).lower(), (*m.getBandedPtr()).upper(), false);
    b = *m.getBandedPtr() / d;
    return b;
  }
  else //if(num==6)
  {
    ZeroMat z;
    return z;
  }
}

SimpleMatrix operator / (const SiconosMatrix &m, int d)
{
  if (d == 0)
    SiconosMatrixException::selfThrow("SimpleMatrix:operator /, division by zero.");

  unsigned int num = m.getNum();
  if (num == 7)
    SiconosMatrixException::selfThrow("SimpleMatrix:operator / (m/int), forbidden for this type of matrix.");

  if (num == 1)
  {
    DenseMat p = *m.getDensePtr() / d;
    return p;
  }
  else if (num == 2)
  {
    TriangMat t = *m.getTriangPtr() / d;
    return t;
  }
  else if (num == 3)
  {
    SymMat s = *m.getSymPtr() / d;
    return s;
  }
  else if (num == 4)
  {
    SparseMat sp = *m.getSparsePtr() / d;
    return sp;
  }
  else if (num == 5)
  {
    BandedMat b;
    b.resize(m.size(0), m.size(1), (*m.getBandedPtr()).lower(), (*m.getBandedPtr()).upper(), false);
    b = *m.getBandedPtr() / d;
    return b;
  }
  else //if(num==6)
  {
    ZeroMat z;
    return z;
  }
}

SimpleMatrix pow(const SimpleMatrix& m, unsigned int power)
{
  if (!m.isSquare())
    SiconosMatrixException::selfThrow("pow(SimpleMatrix), matrix is not square.");

  if (power < 0)
    SiconosMatrixException::selfThrow("pow(SimpleMatrix,n) with negative value is not supported");

  if (power > 0)
  {
    unsigned int num = m.getNum();
    if (num == 1)
    {
      DenseMat p ;
      p = *m.getDensePtr();
      for (unsigned int i = 1; i < power; i++)
        p = prod(p, *m.getDensePtr());
      return p;
    }
    else if (num == 2)
    {
      TriangMat t = *m.getTriangPtr();
      for (unsigned int i = 1; i < power; i++)
        t = prod(t, *m.getTriangPtr());
      return t;
    }
    else if (num == 3)
    {
      SymMat s = *m.getSymPtr();
      for (unsigned int i = 1; i < power; i++)
        s = prod(s, *m.getSymPtr());
      return s;
    }
    else if (num == 4)
    {
      SparseMat sp = *m.getSparsePtr();
      for (unsigned int i = 1; i < power; i++)
        sp = prod(sp, *m.getSparsePtr());
      return sp;
    }
    else if (num == 5)
    {
      DenseMat b = *m.getBandedPtr();
      for (unsigned int i = 1; i < power; i++)
        b = prod(b, *m.getBandedPtr());
      return b;
    }
    else if (num == 6)
    {
      ZeroMat z;
      return z;
    }
    else // if (num==7)
    {
      IdentityMat I;
      return I;
    }
  }
  else// if(power == 0)
  {
    SimpleMatrix p = m;
    p.eye();
    return p;
  }
}

SimpleVector prod(const SiconosMatrix& m, const SiconosVector& v)
{
  if (m.isBlock())
    SiconosMatrixException::selfThrow("prod(matrix,vector) error: not yet implemented for block matrix.");

  if (m.size(1) != v.size())
    SiconosMatrixException::selfThrow("prod(matrix,vector) error: inconsistent sizes.");

  DenseVect res;

  unsigned int numM = m.getNum();
  if (v.isBlock())
  {
    SiconosVector * tmp;
    // TMP: copy into a SimpleVector
    tmp = new SimpleVector(v);
    if (numM == 1)
      res = prod(*m.getDensePtr(), *tmp->getDensePtr());
    else if (numM == 2)
      res = prod(*m.getTriangPtr(), *tmp->getDensePtr());
    else if (numM == 3)
      res = prod(*m.getSymPtr(), *tmp->getDensePtr());
    else if (numM == 4)
      res = prod(*m.getSparsePtr(), *tmp->getDensePtr());
    else if (numM == 5)
      res = prod(*m.getBandedPtr(), *tmp->getDensePtr());
    else if (numM == 6)
      res = 0.0 * *tmp->getDensePtr();
    else if (numM == 7)
      res = *tmp->getDensePtr();
    else
      SiconosMatrixException::selfThrow("prod(matrix,vector) error: unknown matrix type.");

    delete tmp;
  }
  else
  {
    unsigned int numV = v.getNum();
    if (numV == 1)
    {
      if (numM == 1)
        res = prod(*m.getDensePtr(), *v.getDensePtr());
      else if (numM == 2)
        res = prod(*m.getTriangPtr(), *v.getDensePtr());
      else if (numM == 3)
        res = prod(*m.getSymPtr(), *v.getDensePtr());
      else if (numM == 4)
        res = prod(*m.getSparsePtr(), *v.getDensePtr());
      else if (numM == 5)
        res = prod(*m.getBandedPtr(), *v.getDensePtr());
      else if (numM == 6)
        res = 0.0 * *v.getDensePtr();
      else if (numM == 7)
        res = *v.getDensePtr();
      else
        SiconosMatrixException::selfThrow("prod(matrix,vector) error: unknown matrix type.");
    }
    else if (numV == 4)
    {
      if (numM == 1)
        res = prod(*m.getDensePtr(), *v.getSparsePtr());
      else if (numM == 2)
        res = prod(*m.getTriangPtr(), *v.getSparsePtr());
      else if (numM == 3)
        res = prod(*m.getSymPtr(), *v.getSparsePtr());
      else if (numM == 4)
        res = prod(*m.getSparsePtr(), *v.getSparsePtr());
      else if (numM == 5)
        res = prod(*m.getBandedPtr(), *v.getSparsePtr());
      else if (numM == 6)
        res = 0.0 * *v.getSparsePtr();
      else if (numM == 7)
        res = *v.getSparsePtr();
      else
        SiconosMatrixException::selfThrow("prod(matrix,vector) error: unknown matrix type.");
    }
    else
      SiconosMatrixException::selfThrow("prod(matrix,vector) error: unknown vector type.");
  }
  return res;
}

void gemv(const CBLAS_TRANSPOSE transA, double a, const SiconosMatrix& A, const SiconosVector& x, double b, SiconosVector& y)
{
  if (x.isBlock() || y.isBlock())
    SiconosMatrixException::selfThrow("gemv(...) not yet implemented for block vectors.");

  unsigned int numA = A.getNum();
  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();
  std::cout << "num...." << numA << " " << numX << " " << numY << std::endl;
  if (numA == numX && numX == numY && numX == 1) // if all are dense ...
    atlas::gemv(transA, a, *A.getDensePtr(), *x.getDensePtr(), b, *y.getDensePtr());
  else
    SiconosMatrixException::selfThrow("gemv(...) not yet implemented for matrix or vector that are not dense.");
}

void gemv(double a, const SiconosMatrix& A, const SiconosVector& x, double b, SiconosVector& y)
{
  if (x.isBlock() || y.isBlock())
    SiconosMatrixException::selfThrow("gemv(...) not yet implemented for block vectors.");
  unsigned int numA = A.getNum();
  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();
  std::cout << "num...." << numA << " " << numX << " " << numY << std::endl;
  if (numA == numX && numX == numY && numX == 1) // if all are dense ...
    atlas::gemv(a, *A.getDensePtr(), *x.getDensePtr(), b, *y.getDensePtr());
  else
    SiconosMatrixException::selfThrow("gemv(...) not yet implemented for matrix or vector that are not dense.");
}

void prod(const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y)
{
  if (x.isBlock() || y.isBlock())
    SiconosMatrixException::selfThrow("prod(...) not yet implemented for block vectors.");
  unsigned int numA = A.getNum();
  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();
  if (numA == numX && numX == numY && numX == 1) // if all are dense ...
    atlas::gemv(*A.getDensePtr(), *x.getDensePtr(), *y.getDensePtr());
  else
    y = prod(A, x);
}

void gemm(const CBLAS_TRANSPOSE transA, const CBLAS_TRANSPOSE transB, double a, const SiconosMatrix& A, const SiconosMatrix& B, double b, SiconosMatrix& C)
{
  unsigned int numA = A.getNum();
  unsigned int numB = B.getNum();
  unsigned int numC = C.getNum();
  if (numA == numB && numA == numC && numA == 1) // if all are dense ...
    atlas::gemm(transA, transB, a, *A.getDensePtr(), *B.getDensePtr(), b, *C.getDensePtr());
  else
    SiconosMatrixException::selfThrow("gemm(...) not yet implemented for matrices that are not dense.");
}

void gemm(double a, const SiconosMatrix& A, const SiconosMatrix& B, double b, SiconosMatrix& C)
{
  unsigned int numA = A.getNum();
  unsigned int numB = B.getNum();
  unsigned int numC = C.getNum();
  if (numA == numB && numA == numC && numA == 1) // if all are dense ...
    atlas::gemm(a, *A.getDensePtr(), *B.getDensePtr(), b, *C.getDensePtr());
  else
    SiconosMatrixException::selfThrow("gemm(...) not yet implemented for matrices that are not dense.");
}
void prod(const SiconosMatrix& A, const SiconosMatrix& B, SiconosMatrix& C)
{
  unsigned int numA = A.getNum();
  unsigned int numB = B.getNum();
  unsigned int numC = C.getNum();
  if (numA == numB && numA == numC && numA == 1) // if all are dense ...
    atlas::gemm(*A.getDensePtr(), *B.getDensePtr(), *C.getDensePtr());
  else
    C = prod(A, B);
}
