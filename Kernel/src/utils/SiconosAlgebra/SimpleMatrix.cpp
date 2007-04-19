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

// \todo: treat BlockMatrix case in all friend operators (at the time: exceptions).


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
  if (m.isBlock())
    SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy(M) failed: not yet implemented for M being a block matrix. You can try with BlockMatrix::getBlockPtr function?");

  if (num != 1)
    SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : the current matrix is not dense, a copy into its data may change its type.");

  if (i >= dim[0])
    SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : row_min given is out of range");

  if (j >= dim[1])
    SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : col_min given is out of range");

  unsigned int num2 = m.getNum();
  unsigned int row_max = i, col_max = j;

  unsigned int sizeM1 = 0, sizeM2 = 0;
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

  if (row_max > dim[0])
    SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");
  if (col_max > dim[1])
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
  if (m.isBlock())
    SiconosMatrixException::selfThrow("SimpleMatrix::getBlock(M) failed: not yet implemented for M being a block matrix. You can try with BlockMatrix::getBlockPtr function?");


  // We only accept dense matrix for m.
  if (m.getNum() != 1)
    SiconosMatrixException::selfThrow("getBlock(i,j,m) : m must be a dense matrix.");

  if (row_min >= dim[0] || row_min < 0)
    SiconosMatrixException::selfThrow("getBlock : row_min given is out of range");

  if (col_min >= dim[1] || col_min < 0)
    SiconosMatrixException::selfThrow("getBlock : col_min given is out of range");
  unsigned int row_max, col_max;
  row_max = m.getDensePtr()->size1() + row_min;
  col_max = m.getDensePtr()->size2() + col_min;

  if (row_max > dim[0] || row_max < 0)
    SiconosMatrixException::selfThrow("getBlock : inconsistent sizes");

  if (col_max > dim[1] || col_max < 0)
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

  if (r >= dim[0] || r < 0)
    SiconosMatrixException::selfThrow("getRow : row is out of range");

  if (vect.size() != dim[1])
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
  if (r >= dim[0] || r < 0)
    SiconosMatrixException::selfThrow("setRow : row is out of range");

  if (vect.size() != dim[1])
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
  if (r >= dim[1] || r < 0)
    SiconosMatrixException::selfThrow("getCol : col is out of range");

  if (vect.size() != dim[0])
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

  if (r >= dim[1] || r < 0)
    SiconosMatrixException::selfThrow("setCol : col is out of range");

  if (vect.size() != dim[0])
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

  if (m.isBlock())
    SiconosMatrixException::selfThrow("SimpleMatrix::trans(m) failed, not yet implemented for m being a BlockMatrix.");


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
  if (row >= dim[0] || col >= dim[1])
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
  if (row >= dim[0] || col >= dim[1])
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
  if (dim[0] != m.size(0) || dim[1] != m.size(1))
    SiconosMatrixException::selfThrow("SimpleMatrix operator = failed. Inconsistent sizes.");

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
  if (m.isBlock())
    SiconosMatrixException::selfThrow("SimpleMatrix operator = M failed. Not yet implemented for M being a BlockMatrix.");

  if (dim[0] != m.size(0) || dim[1] != m.size(1))
    SiconosMatrixException::selfThrow("SimpleMatrix operator = failed. Inconsistent sizes.");

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
  if (m.isBlock())
    SiconosMatrixException::selfThrow("SimpleMatrix operator += M failed. Not yet implemented for M being a BlockMatrix.");

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
  if (m.isBlock())
    SiconosMatrixException::selfThrow("SimpleMatrix operator -= M failed. Not yet implemented for M being a BlockMatrix.");

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

  ipiv.resize(dim[0]);
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
  if (B.isBlock())
    SiconosMatrixException::selfThrow("SimpleMatrix PLUForwardBackwardInPlace(M) failed. Not yet implemented for M being a BlockMatrix.");
  int info;
  if (!isPLUFactorized) // call gesv => LU-factorize+solve
  {
    // solve system:
    ipiv.resize(dim[0]);
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
  if (B.isBlock())
    SiconosMatrixException::selfThrow("SimpleMatrix PLUForwardBackwardInPlace(V) failed. Not yet implemented for V being a BlockVector.");
  DenseMat tmpB(B.size(), 1);
  ublas::column(tmpB, 0) = *(B.getDensePtr()); // Conversion of vector to matrix. Temporary solution.
  int info;
  if (!isPLUFactorized) // call gesv => LU-factorize+solve
  {
    // solve system:
    ipiv.resize(dim[0]);
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

void SimpleMatrix::scal(double a, const SiconosMatrix &A)
{
  if (this == &A)
    *this *= a;
  else if (A.isBlock())
  {
    SiconosMatrixException::selfThrow("SimpleMatrix::scal(a,A) failed. Not yet implemented for A a BlockMatrix.");
  }
  else
  {
    if (num == 6 || num == 7)
      SiconosMatrixException::selfThrow("SimpleMatrix::scal(a,A) failed. Try to assign value to a null or identity matrix.");

    unsigned int numA = A.getNum();

    if (numA == 6)
      zero();
    else if (numA == 7)
    {
      eye();
      *this *= a;
    }
    else
    {
      if (numA == num)
      {
        switch (num)
        {
        case 1:
          noalias(*mat.Dense) = a ** A.getDensePtr();
          break;
        case 2:
          noalias(*mat.Triang) = a ** A.getTriangPtr();
          break;
        case 3:
          noalias(*mat.Sym) = a ** A.getSymPtr();
          break;
        case 4:
          noalias(*mat.Sparse) = a ** A.getSparsePtr();
          break;
        case 5:
          noalias(*mat.Banded) = a ** A.getBandedPtr();
          break;
        }
      }
      else
      {
        if (num != 1)
          SiconosMatrixException::selfThrow("SimpleMatrix::B = a*A failed. A and B types do not fit together.");
        switch (numA)
        {
        case 1:
          noalias(*mat.Dense) = a ** A.getDensePtr();
          break;
        case 2:
          noalias(*mat.Dense) = a ** A.getTriangPtr();
          break;
        case 3:
          noalias(*mat.Dense) = a ** A.getSymPtr();
          break;
        case 4:
          noalias(*mat.Dense) = a ** A.getSparsePtr();
          break;
        case 5:
          noalias(*mat.Dense) = a ** A.getBandedPtr();
          break;
        }
      }
    }
  }
}

void SimpleMatrix::scal(int a, const SiconosMatrix &A)
{
  if (this == &A)
    *this *= a;
  else if (A.isBlock())
  {
    SiconosMatrixException::selfThrow("SimpleMatrix::scal(a,A) failed. Not yet implemented for A a BlockMatrix.");
  }
  else
  {
    if (num == 6 || num == 7)
      SiconosMatrixException::selfThrow("SimpleMatrix::scal(a,A) failed. Try to assign value to a null or identity matrix.");

    unsigned int numA = A.getNum();

    if (numA == 6)
      zero();
    else if (numA == 7)
    {
      eye();
      *this *= a;
    }
    else
    {
      if (numA == num)
      {
        switch (num)
        {
        case 1:
          noalias(*mat.Dense) = a ** A.getDensePtr();
          break;
        case 2:
          noalias(*mat.Triang) = a ** A.getTriangPtr();
          break;
        case 3:
          noalias(*mat.Sym) = a ** A.getSymPtr();
          break;
        case 4:
          noalias(*mat.Sparse) = a ** A.getSparsePtr();
          break;
        case 5:
          noalias(*mat.Banded) = a ** A.getBandedPtr();
          break;
        }
      }
      else
      {
        if (num != 1)
          SiconosMatrixException::selfThrow("SimpleMatrix::scal(a,A) failed. A and B types do not fit together.");
        switch (numA)
        {
        case 1:
          noalias(*mat.Dense) = a ** A.getDensePtr();
          break;
        case 2:
          noalias(*mat.Dense) = a ** A.getTriangPtr();
          break;
        case 3:
          noalias(*mat.Dense) = a ** A.getSymPtr();
          break;
        case 4:
          noalias(*mat.Dense) = a ** A.getSparsePtr();
          break;
        case 5:
          noalias(*mat.Dense) = a ** A.getBandedPtr();
          break;
        }
      }
    }
  }
}

bool operator == (const SiconosMatrix &m, const SiconosMatrix &x)
{
  if (m.isBlock() ||  x.isBlock())
    SiconosMatrixException::selfThrow("op == (const SiconosMatrix, const SiconosMatrix) : incompatible type of matrix");
  double norm = (sub(m, x)).normInf();
  return (norm < tolerance);
}

const SimpleMatrix operator + (const SiconosMatrix &x, const SiconosMatrix &m)
{
  if ((x.size(0) != m.size(0)) || x.size(1) != m.size(1))
    SiconosMatrixException::selfThrow("Matrix addition: inconsistent sizes");

  if (x.isBlock() || m.isBlock())
    SiconosMatrixException::selfThrow("Matrix addition: not yet implemented for BlockMatrix.");

  unsigned int numX = x.getNum();
  if (numX != m.getNum())
    SiconosMatrixException::selfThrow("SimpleMatrix:Matrix operator +: addition of matrices of different types. Use 'add' function.");

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
    BandedMat b ;
    b.resize(m.size(0), m.size(1), (*m.getBandedPtr()).lower(), (*m.getBandedPtr()).upper(), false);
    b = *x.getBandedPtr() + *m.getBandedPtr();
    return b;
  }
  else if (numX == 6)
  {
    ZeroMat z(x.size(0), x.size(1));
    return z;
  }
  else // if(numX ==7)
  {
    DenseMat p = *x.getIdentityPtr();
    p *= 2.0;
    return p;
  }
}

const SimpleMatrix add(const SiconosMatrix &x, const SiconosMatrix& m)
{
  if ((x.size(0) != m.size(0)) || (x.size(1) != m.size(1)))
    SiconosMatrixException::selfThrow("SimpleMatrix: function add, inconsistent sizes.");

  if (x.isBlock() || m.isBlock())
    SiconosMatrixException::selfThrow("Matrix, add function: not yet implemented for BlockMatrix.");

  unsigned int numM = m.getNum();
  unsigned int numX = x.getNum();

  if (numX == 6) // if x = zero
    return m;
  else if (numM == 6) // if m = 0
    return x;
  else
  {
    if (numX == numM) // if x and m are of the same type.
    {
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
      else if (numX == 3)
      {
        SparseMat s = *x.getSparsePtr() + *m.getSparsePtr();
        return s;
      }
      else if (numX == 5)
      {
        BandedMat b = *x.getBandedPtr() + *m.getBandedPtr();
        return b;
      }
      else // if(numX==7)
      {
        DenseMat p = *x.getIdentityPtr();
        p *= 2.0;
        return p;
      }
    }
    else
    {
      // if x and m are different, result is a DenseMat.
      DenseMat p(x.size(0), x.size(1));
      if (numX == 1) // if x dense.
        switch (numM)
        {
        case 2:
          noalias(p) = *x.getDensePtr() + *m.getTriangPtr();
          break;
        case 3:
          noalias(p) = *x.getDensePtr() + *m.getSymPtr();
          break;
        case 4:
          noalias(p) = *x.getDensePtr() + *m.getSparsePtr();
          break;
        case 5:
          noalias(p) = *x.getDensePtr() + *m.getBandedPtr();
          break;
        case 7:
          noalias(p) = *x.getDensePtr() + *m.getIdentityPtr();
          break;
        default:
          SiconosMatrixException::selfThrow("Matrix function add(A,B): invalid type of matrix");
        }
      else if (numX == 2)
        switch (numM)
        {
        case 1:
          noalias(p) = *x.getTriangPtr() + *m.getDensePtr();
          break;
        case 3:
          noalias(p) = *x.getTriangPtr() + *m.getSymPtr();
          break;
        case 4:
          noalias(p) = *x.getTriangPtr() + *m.getSparsePtr();
          break;
        case 5:
          noalias(p) = *x.getTriangPtr() + *m.getBandedPtr();
          break;
        case 7:
          noalias(p) = *x.getTriangPtr() + *m.getIdentityPtr();
          break;
        default:
          SiconosMatrixException::selfThrow("Matrix function add(A,B): invalid type of matrix");
        }
      else if (numX == 3)
        switch (numM)
        {
        case 1:
          noalias(p) = *x.getSymPtr() + *m.getDensePtr();
          break;
        case 2:
          noalias(p) = *x.getSymPtr() + *m.getTriangPtr();
          break;
        case 4:
          noalias(p) = *x.getSymPtr() + *m.getSparsePtr();
          break;
        case 5:
          noalias(p) = *x.getSymPtr() + *m.getBandedPtr();
          break;
        case 7:
          noalias(p) = *x.getSymPtr() + *m.getIdentityPtr();
          break;
        default:
          SiconosMatrixException::selfThrow("Matrix function add(A,B): invalid type of matrix");
        }
      else if (numX == 4)
        switch (numM)
        {
        case 1:
          noalias(p) = *x.getSparsePtr() + *m.getDensePtr();
          break;
        case 2:
          noalias(p) = *x.getSparsePtr() + *m.getTriangPtr();
          break;
        case 3:
          noalias(p) = *x.getSparsePtr() + *m.getSymPtr();
          break;
        case 5:
          noalias(p) = *x.getSparsePtr() + *m.getBandedPtr();
          break;
        case 7:
          noalias(p) = *x.getSparsePtr() + *m.getIdentityPtr();
          break;
        default:
          SiconosMatrixException::selfThrow("Matrix function add(A,B): invalid type of matrix");
        }
      else if (numX == 5)
        switch (numM)
        {
        case 1:
          noalias(p) = *x.getBandedPtr() + *m.getDensePtr();
          break;
        case 2:
          noalias(p) = *x.getBandedPtr() + *m.getTriangPtr();
          break;
        case 3:
          noalias(p) = *x.getBandedPtr() + *m.getSymPtr();
          break;
        case 4:
          noalias(p) = *x.getBandedPtr() + *m.getSparsePtr();
          break;
        case 7:
          noalias(p) = *x.getBandedPtr() + *m.getIdentityPtr();
          break;
        default:
          SiconosMatrixException::selfThrow("Matrix function add(A,B): invalid type of matrix");
        }
      else if (numX == 7)
        switch (numM)
        {
        case 1:
          noalias(p) = *x.getIdentityPtr() + *m.getDensePtr();
          break;
        case 2:
          noalias(p) = *x.getIdentityPtr() + *m.getTriangPtr();
          break;
        case 3:
          noalias(p) = *x.getIdentityPtr() + *m.getSymPtr();
          break;
        case 4:
          noalias(p) = *x.getIdentityPtr() + *m.getSparsePtr();
          break;
        case 7:
          noalias(p) = *x.getIdentityPtr() + *m.getIdentityPtr();
          break;
        default:
          SiconosMatrixException::selfThrow("Matrix function add(A,B): invalid type of matrix");
        }
      else
        SiconosMatrixException::selfThrow("Matrix function add(A,B): invalid type of matrix");
      return p;
    }
  }
}

void add(const SiconosMatrix &x, const SiconosMatrix& m, SiconosMatrix& res)
{
  if ((x.size(0) != m.size(0)) || (x.size(1) != m.size(1)))
    SiconosMatrixException::selfThrow("Matrix function add: inconsistent sizes");
  if ((x.size(0) != res.size(0)) || (x.size(1) != res.size(1)))
    SiconosMatrixException::selfThrow("Matrix function add: inconsistent sizes");
  if (x.isBlock() || m.isBlock() || res.isBlock())
    SiconosMatrixException::selfThrow("Matrix, add function: not yet implemented for BlockMatrix.");

  unsigned int numX = x.getNum();
  unsigned int numM = m.getNum();
  unsigned int numRes = res.getNum();

  if (&res == &m)
    res += x;
  else if (&res == &x)
    res += m;
  else
  {
    if (numX == 6) // if x = zero
      res = m;
    else if (numM == 6) // if m = 0
      res = x;
    else
    {
      if (numX == numM) // if x and m are of the same type.
      {
        if (numRes != numX && numX != 7)
          SiconosMatrixException::selfThrow("Matrix function add(A,B,C): wrong type for C.");
        if (numX == 1)
          noalias(*res.getDensePtr()) = *x.getDensePtr() + *m.getDensePtr();
        else if (numX == 2)
          noalias(*res.getTriangPtr()) = *x.getTriangPtr() + *m.getTriangPtr();
        else if (numX == 3)
          noalias(*res.getSymPtr()) = *x.getSymPtr() + *m.getSymPtr();
        else if (numX == 4)
          noalias(*res.getSparsePtr()) = *x.getSparsePtr() + *m.getSparsePtr();
        else if (numX == 5)
          noalias(*res.getBandedPtr()) = *x.getBandedPtr() + *m.getBandedPtr();
        else // if(numX==7)
        {
          res = x;
          res *= 2.0;
        }
      }
      else
      {
        // if m and x are of different type, res must be dense.
        if (numRes != 1)
          SiconosMatrixException::selfThrow("Matrix function add(A,B,C): wrong type for C.");

        if (numX == 1) // if x dense.
          switch (numM)
          {
          case 2:
            noalias(*res.getDensePtr()) = *x.getDensePtr() + *m.getTriangPtr();
            break;
          case 3:
            noalias(*res.getDensePtr()) = *x.getDensePtr() + *m.getSymPtr();
            break;
          case 4:
            noalias(*res.getDensePtr()) = *x.getDensePtr() + *m.getSparsePtr();
            break;
          case 5:
            noalias(*res.getDensePtr()) = *x.getDensePtr() + *m.getBandedPtr();
            break;
          case 7:
            noalias(*res.getDensePtr()) = *x.getDensePtr() + *m.getIdentityPtr();
            break;
          default:
            SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
          }
        else if (numX == 2)
          switch (numM)
          {
          case 1:
            noalias(*res.getDensePtr()) = *x.getTriangPtr() + *m.getDensePtr();
            break;
          case 3:
            noalias(*res.getDensePtr()) = *x.getTriangPtr() + *m.getSymPtr();
            break;
          case 4:
            noalias(*res.getDensePtr()) = *x.getTriangPtr() + *m.getSparsePtr();
            break;
          case 5:
            noalias(*res.getDensePtr()) = *x.getTriangPtr() + *m.getBandedPtr();
            break;
          case 7:
            noalias(*res.getDensePtr()) = *x.getTriangPtr() + *m.getIdentityPtr();
            break;
          default:
            SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
          }
        else if (numX == 3)
          switch (numM)
          {
          case 1:
            noalias(*res.getDensePtr()) = *x.getSymPtr() + *m.getDensePtr();
            break;
          case 2:
            noalias(*res.getDensePtr()) = *x.getSymPtr() + *m.getTriangPtr();
            break;
          case 4:
            noalias(*res.getDensePtr()) = *x.getSymPtr() + *m.getSparsePtr();
            break;
          case 5:
            noalias(*res.getDensePtr()) = *x.getSymPtr() + *m.getBandedPtr();
            break;
          case 7:
            noalias(*res.getDensePtr()) = *x.getSymPtr() + *m.getIdentityPtr();
            break;
          default:
            SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
          }
        else if (numX == 4)
          switch (numM)
          {
          case 1:
            noalias(*res.getDensePtr()) = *x.getSparsePtr() + *m.getDensePtr();
            break;
          case 2:
            noalias(*res.getDensePtr()) = *x.getSparsePtr() + *m.getTriangPtr();
            break;
          case 3:
            noalias(*res.getDensePtr()) = *x.getSparsePtr() + *m.getSymPtr();
            break;
          case 5:
            noalias(*res.getDensePtr()) = *x.getSparsePtr() + *m.getBandedPtr();
            break;
          case 7:
            noalias(*res.getDensePtr()) = *x.getSparsePtr() + *m.getIdentityPtr();
            break;
          default:
            SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
          }
        else if (numX == 5)
          switch (numM)
          {
          case 1:
            noalias(*res.getDensePtr()) = *x.getBandedPtr() + *m.getDensePtr();
            break;
          case 2:
            noalias(*res.getDensePtr()) = *x.getBandedPtr() + *m.getTriangPtr();
            break;
          case 3:
            noalias(*res.getDensePtr()) = *x.getBandedPtr() + *m.getSymPtr();
            break;
          case 4:
            noalias(*res.getDensePtr()) = *x.getBandedPtr() + *m.getSparsePtr();
            break;
          case 7:
            noalias(*res.getDensePtr()) = *x.getBandedPtr() + *m.getIdentityPtr();
            break;
          default:
            SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
          }
        else if (numX == 7)
          switch (numM)
          {
          case 1:
            noalias(*res.getDensePtr()) = *x.getIdentityPtr() + *m.getDensePtr();
            break;
          case 2:
            noalias(*res.getDensePtr()) = *x.getIdentityPtr() + *m.getTriangPtr();
            break;
          case 3:
            noalias(*res.getDensePtr()) = *x.getIdentityPtr() + *m.getSymPtr();
            break;
          case 4:
            noalias(*res.getDensePtr()) = *x.getIdentityPtr() + *m.getSparsePtr();
            break;
          case 7:
            noalias(*res.getDensePtr()) = *x.getIdentityPtr() + *m.getIdentityPtr();
            break;
          default:
            SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");
          }
        else
          SiconosMatrixException::selfThrow("Matrix function add(A,B,C): invalid type of matrix");

      }
    }
  }
}

const SimpleMatrix operator - (const SiconosMatrix &x, const SiconosMatrix &m)
{
  if ((x.size(0) != m.size(0)) || (x.size(1) != m.size(1)))
    SiconosMatrixException::selfThrow("Matrix subtraction: inconsistent sizes");
  if (x.isBlock() || m.isBlock())
    SiconosMatrixException::selfThrow("Matrix, operator -: not yet implemented for BlockMatrix.");

  unsigned int numX = x.getNum();
  if (numX != m.getNum())
    SiconosMatrixException::selfThrow("SimpleMatrix:Matrix subtraction: use function sub in order to subtract matrices of different type");

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
  else
  {
    //if(numX == 6 || numX ==7){
    ZeroMat z(x.size(0), x.size(1));
    return z;
  }
}

const SimpleMatrix sub(const SiconosMatrix &x, const SiconosMatrix& m)
{
  if ((x.size(0) != m.size(0)) || (x.size(1) != m.size(1)))
    SiconosMatrixException::selfThrow("Function sub(A,B), inconsistent sizes.");
  if (x.isBlock() || m.isBlock())
    SiconosMatrixException::selfThrow("Matrix, sub function: not yet implemented for BlockMatrix.");

  unsigned int numM = m.getNum();
  unsigned int numX = x.getNum();

  if (numM == 6) // if m = zero
    return x;
  else
  {
    if (numX == numM) // if x and m are of the same type.
    {
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
      else if (numX == 3)
      {
        SparseMat s = *x.getSparsePtr() - *m.getSparsePtr();
        return s;
      }
      else if (numX == 5)
      {
        BandedMat b = *x.getBandedPtr() - *m.getBandedPtr();
        return b;
      }
      else // if(numX==7)
      {
        ZeroMat z(x.size(0), x.size(1));
        return z;
      }
    }
    else
    {
      // if x and m are different, result is a DenseMat.
      DenseMat p(x.size(0), x.size(1));
      if (numX == 1) // if x dense.
        switch (numM)
        {
        case 2:
          noalias(p) = *x.getDensePtr() - *m.getTriangPtr();
          break;
        case 3:
          noalias(p) = *x.getDensePtr() - *m.getSymPtr();
          break;
        case 4:
          noalias(p) = *x.getDensePtr() - *m.getSparsePtr();
          break;
        case 5:
          noalias(p) = *x.getDensePtr() - *m.getBandedPtr();
          break;
        case 7:
          noalias(p) = *x.getDensePtr() - *m.getIdentityPtr();
          break;
        default:
          SiconosMatrixException::selfThrow("Matrix function sub(A,B): invalid type of matrix");
        }
      else if (numX == 2)
        switch (numM)
        {
        case 1:
          noalias(p) = *x.getTriangPtr() - *m.getDensePtr();
          break;
        case 3:
          noalias(p) = *x.getTriangPtr() - *m.getSymPtr();
          break;
        case 4:
          noalias(p) = *x.getTriangPtr() - *m.getSparsePtr();
          break;
        case 5:
          noalias(p) = *x.getTriangPtr() - *m.getBandedPtr();
          break;
        case 7:
          noalias(p) = *x.getTriangPtr() - *m.getIdentityPtr();
          break;
        default:
          SiconosMatrixException::selfThrow("Matrix function sub(A,B): invalid type of matrix");
        }
      else if (numX == 3)
        switch (numM)
        {
        case 1:
          noalias(p) = *x.getSymPtr() - *m.getDensePtr();
          break;
        case 2:
          noalias(p) = *x.getSymPtr() - *m.getTriangPtr();
          break;
        case 4:
          noalias(p) = *x.getSymPtr() - *m.getSparsePtr();
          break;
        case 5:
          noalias(p) = *x.getSymPtr() - *m.getBandedPtr();
          break;
        case 7:
          noalias(p) = *x.getSymPtr() - *m.getIdentityPtr();
          break;
        default:
          SiconosMatrixException::selfThrow("Matrix function sub(A,B): invalid type of matrix");
        }
      else if (numX == 4)
        switch (numM)
        {
        case 1:
          noalias(p) = *x.getSparsePtr() - *m.getDensePtr();
          break;
        case 2:
          noalias(p) = *x.getSparsePtr() - *m.getTriangPtr();
          break;
        case 3:
          noalias(p) = *x.getSparsePtr() - *m.getSymPtr();
          break;
        case 5:
          noalias(p) = *x.getSparsePtr() - *m.getBandedPtr();
          break;
        case 7:
          noalias(p) = *x.getSparsePtr() - *m.getIdentityPtr();
          break;
        default:
          SiconosMatrixException::selfThrow("Matrix function sub(A,B): invalid type of matrix");
        }
      else if (numX == 5)
        switch (numM)
        {
        case 1:
          noalias(p) = *x.getBandedPtr() - *m.getDensePtr();
          break;
        case 2:
          noalias(p) = *x.getBandedPtr() - *m.getTriangPtr();
          break;
        case 3:
          noalias(p) = *x.getBandedPtr() - *m.getSymPtr();
          break;
        case 4:
          noalias(p) = *x.getBandedPtr() - *m.getSparsePtr();
          break;
        case 7:
          noalias(p) = *x.getBandedPtr() - *m.getIdentityPtr();
          break;
        default:
          SiconosMatrixException::selfThrow("Matrix function sub(A,B): invalid type of matrix");
        }
      else if (numX == 6)
      {
        if (numM == 1)
        {
          p = *m.getDensePtr();
          p *= -1.0;
        }
        else if (numM == 2)
        {
          TriangMat t = *m.getTriangPtr();
          t *= -1.0;
          return t;
        }
        else if (numM == 3)
        {
          SymMat s = *m.getSymPtr();
          s *= -1.0;
          return s;
        }
        else if (numM == 4)
        {
          SparseMat s =  *m.getSparsePtr();
          s *= -1.0;
          return s;
        }
        else if (numM == 5)
        {
          BandedMat b =  *m.getBandedPtr();
          b *= -1.0;
          return b;
        }
        else if (numM == 7)
        {
          p =  *m.getIdentityPtr();
          p *= -1.0;
        }
      }
      else if (numX == 7)
        switch (numM)
        {
        case 1:
          noalias(p) = *x.getIdentityPtr() - *m.getDensePtr();
          break;
        case 2:
          noalias(p) = *x.getIdentityPtr() - *m.getTriangPtr();
          break;
        case 3:
          noalias(p) = *x.getIdentityPtr() - *m.getSymPtr();
          break;
        case 4:
          noalias(p) = *x.getIdentityPtr() - *m.getSparsePtr();
          break;
        case 5:
          noalias(p) = *x.getIdentityPtr() - *m.getBandedPtr();
          break;
        default:
          SiconosMatrixException::selfThrow("Matrix function sub(A,B): invalid type of matrix");
        }
      else
        SiconosMatrixException::selfThrow("Matrix function sub(A,B): invalid type of matrix");
      return p;
    }
  }
}

void sub(const SiconosMatrix &x, const SiconosMatrix& m, SiconosMatrix& res)
{
  if ((x.size(0) != m.size(0)) || (x.size(1) != m.size(1)))
    SiconosMatrixException::selfThrow("Matrix function sub(A,B,C): inconsistent sizes");
  if ((x.size(0) != res.size(0)) || (x.size(1) != res.size(1)))
    SiconosMatrixException::selfThrow("Matrix function sub(A,B,C): inconsistent sizes");
  if (x.isBlock() || m.isBlock() || res.isBlock())
    SiconosMatrixException::selfThrow("Matrix, sub function: not yet implemented for BlockMatrix.");

  if (&res == &m)
  {
    res -= x;
    res *= -1.0;
  }
  else if (&res == &x)
    res -= m;
  else
  {

    unsigned int numX = x.getNum();
    unsigned int numM = m.getNum();
    unsigned int numRes = res.getNum();

    if (numM == 6) // if m = zero
      res = x;
    else if (numX == 6) // if x = 0
    {
      res = m;
      res *= -1.0;
    }
    else
    {
      if (numX == numM) // if x and m are of the same type.
      {
        if (numRes != numX && numX != 7)
          SiconosMatrixException::selfThrow("Matrix function sub(A,B,C): wrong type for C.");
        if (numX == 1)
          noalias(*res.getDensePtr()) = *x.getDensePtr() - *m.getDensePtr();
        else if (numX == 2)
          noalias(*res.getTriangPtr()) = *x.getTriangPtr() - *m.getTriangPtr();
        else if (numX == 3)
          noalias(*res.getSymPtr()) = *x.getSymPtr() - *m.getSymPtr();
        else if (numX == 4)
          noalias(*res.getSparsePtr()) = *x.getSparsePtr() - *m.getSparsePtr();
        else if (numX == 5)
          noalias(*res.getBandedPtr()) = *x.getBandedPtr() - *m.getBandedPtr();
        else // if(numX==7)
          res.zero();
      }
      else
      {
        // if m and x are of different type, res must be dense.
        if (numRes != 1)
          SiconosMatrixException::selfThrow("Matrix function sub(A,B,C): wrong type for C.");

        if (numX == 1) // if x dense.
          switch (numM)
          {
          case 2:
            noalias(*res.getDensePtr()) = *x.getDensePtr() - *m.getTriangPtr();
            break;
          case 3:
            noalias(*res.getDensePtr()) = *x.getDensePtr() - *m.getSymPtr();
            break;
          case 4:
            noalias(*res.getDensePtr()) = *x.getDensePtr() - *m.getSparsePtr();
            break;
          case 5:
            noalias(*res.getDensePtr()) = *x.getDensePtr() - *m.getBandedPtr();
            break;
          case 7:
            noalias(*res.getDensePtr()) = *x.getDensePtr() - *m.getIdentityPtr();
            break;
          default:
            SiconosMatrixException::selfThrow("Matrix function sub(A,B,C): invalid type of matrix");
          }
        else if (numX == 2)
          switch (numM)
          {
          case 1:
            noalias(*res.getDensePtr()) = *x.getTriangPtr() - *m.getDensePtr();
            break;
          case 3:
            noalias(*res.getDensePtr()) = *x.getTriangPtr() - *m.getSymPtr();
            break;
          case 4:
            noalias(*res.getDensePtr()) = *x.getTriangPtr() - *m.getSparsePtr();
            break;
          case 5:
            noalias(*res.getDensePtr()) = *x.getTriangPtr() - *m.getBandedPtr();
            break;
          case 7:
            noalias(*res.getDensePtr()) = *x.getTriangPtr() - *m.getIdentityPtr();
            break;
          default:
            SiconosMatrixException::selfThrow("Matrix function sub(A,B,C): invalid type of matrix");
          }
        else if (numX == 3)
          switch (numM)
          {
          case 1:
            noalias(*res.getDensePtr()) = *x.getSymPtr() - *m.getDensePtr();
            break;
          case 2:
            noalias(*res.getDensePtr()) = *x.getSymPtr() - *m.getTriangPtr();
            break;
          case 4:
            noalias(*res.getDensePtr()) = *x.getSymPtr() - *m.getSparsePtr();
            break;
          case 5:
            noalias(*res.getDensePtr()) = *x.getSymPtr() - *m.getBandedPtr();
            break;
          case 7:
            noalias(*res.getDensePtr()) = *x.getSymPtr() - *m.getIdentityPtr();
            break;
          default:
            SiconosMatrixException::selfThrow("Matrix function sub(A,B,C): invalid type of matrix");
          }
        else if (numX == 4)
          switch (numM)
          {
          case 1:
            noalias(*res.getDensePtr()) = *x.getSparsePtr() - *m.getDensePtr();
            break;
          case 2:
            noalias(*res.getDensePtr()) = *x.getSparsePtr() - *m.getTriangPtr();
            break;
          case 3:
            noalias(*res.getDensePtr()) = *x.getSparsePtr() - *m.getSymPtr();
            break;
          case 5:
            noalias(*res.getDensePtr()) = *x.getSparsePtr() - *m.getBandedPtr();
            break;
          case 7:
            noalias(*res.getDensePtr()) = *x.getSparsePtr() - *m.getIdentityPtr();
            break;
          default:
            SiconosMatrixException::selfThrow("Matrix function sub(A,B,C): invalid type of matrix");
          }
        else if (numX == 5)
          switch (numM)
          {
          case 1:
            noalias(*res.getDensePtr()) = *x.getBandedPtr() - *m.getDensePtr();
            break;
          case 2:
            noalias(*res.getDensePtr()) = *x.getBandedPtr() - *m.getTriangPtr();
            break;
          case 3:
            noalias(*res.getDensePtr()) = *x.getBandedPtr() - *m.getSymPtr();
            break;
          case 4:
            noalias(*res.getDensePtr()) = *x.getBandedPtr() - *m.getSparsePtr();
            break;
          case 7:
            noalias(*res.getDensePtr()) = *x.getBandedPtr() - *m.getIdentityPtr();
            break;
          default:
            SiconosMatrixException::selfThrow("Matrix function sub(A,B,C): invalid type of matrix");
          }
        else if (numX == 7)
          switch (numM)
          {
          case 1:
            noalias(*res.getDensePtr()) = *x.getIdentityPtr() - *m.getDensePtr();
            break;
          case 2:
            noalias(*res.getDensePtr()) = *x.getIdentityPtr() - *m.getTriangPtr();
            break;
          case 3:
            noalias(*res.getDensePtr()) = *x.getIdentityPtr() - *m.getSymPtr();
            break;
          case 4:
            noalias(*res.getDensePtr()) = *x.getIdentityPtr() - *m.getSparsePtr();
            break;
          case 5:
            noalias(*res.getDensePtr()) = *x.getIdentityPtr() - *m.getBandedPtr();
            break;
          default:
            SiconosMatrixException::selfThrow("Matrix function sub(A,B,C): invalid type of matrix");
          }
        else
          SiconosMatrixException::selfThrow("Matrix function sub(A,B,C): invalid type of matrix");

      }
    }
  }
}

const SimpleMatrix operator * (const SiconosMatrix &x, const SiconosMatrix &m)
{

  if ((x.size(1) != m.size(0)))
    SiconosMatrixException::selfThrow("Matrix product : inconsistent sizes");

  if (x.isBlock() || m.isBlock())
    SiconosMatrixException::selfThrow("Matrix function prod : not yet implemented for block matrices.");

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
  else //if(numX == 6 || numX == 7){
    return x;
}

const SimpleMatrix prod(const SiconosMatrix &x, const SiconosMatrix& m)
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
  else if (numM == 7)
    return x;
  else if (numX == 7)
    return m;
  else
  {
    if (numM == 1)
    {
      if (numX == 1)
      {
        DenseMat p = prod(*x.getDensePtr(), *m.getDensePtr());
        return p;
      }
      else if (numX == 2)
      {
        DenseMat p = prod(*x.getTriangPtr(), *m.getDensePtr());
        return p;
      }

      else if (numX == 3)
      {
        DenseMat p  = prod(*x.getSymPtr(), *m.getDensePtr());
        return p;
      }

      else if (numX == 4)
      {
        DenseMat p = prod(*x.getSparsePtr(), *m.getDensePtr());
        return p;
      }

      else// if(numX==5)
      {
        DenseMat p  = prod(*x.getBandedPtr(), *m.getDensePtr());
        return p;
      }
    }
    else if (numM == 2)
    {
      if (numX == 1)
      {
        DenseMat p  = prod(*x.getDensePtr(), *m.getTriangPtr());
        return p;
      }
      else if (numX == 2)
      {
        TriangMat t = prod(*x.getTriangPtr(), *m.getTriangPtr());
        return t;
      }
      else if (numX == 3)
      {
        DenseMat p  = prod(*x.getSymPtr(), *m.getTriangPtr());
        return p;
      }
      else if (numX == 4)
      {
        DenseMat p  = prod(*x.getSparsePtr(), *m.getTriangPtr());
        return p;
      }
      else //if(numX==5)
      {
        DenseMat p  = prod(*x.getBandedPtr(), *m.getTriangPtr());
        return p;
      }
    }
    else if (numM == 3)
    {
      if (numX == 1)
      {
        DenseMat p  = prod(*x.getDensePtr(), *m.getSymPtr());
        return p;
      }
      else if (numX == 2)
      {
        DenseMat p  = prod(*x.getTriangPtr(), *m.getSymPtr());
        return p;
      }
      else if (numX == 3)
      {
        SymMat s = prod(*x.getSymPtr(), *m.getSymPtr());
        return s;
      }
      else if (numX == 4)
      {
        DenseMat p  = prod(*x.getSparsePtr(), *m.getSymPtr());
        return p;
      }
      else
      {
        //if(numX==5){
        DenseMat p  = prod(*x.getBandedPtr(), *m.getSymPtr());
        return p;
      }
    }
    else if (numM == 4)
    {
      if (numX == 1)
      {
        DenseMat p = prod(*x.getDensePtr(), *m.getSparsePtr());
        return p;
      }
      else if (numX == 2)
      {
        DenseMat p = prod(*x.getTriangPtr(), *m.getSparsePtr());
        return p;
      }
      else if (numX == 3)
      {
        DenseMat p = prod(*x.getSymPtr(), *m.getSparsePtr());
        return p;
      }
      else if (numX == 4)
      {
        SparseMat sp = prod(*x.getSparsePtr(), *m.getSparsePtr());
        return sp;
      }
      else //if(numX==5){
      {
        DenseMat p = prod(*x.getBandedPtr(), *m.getSparsePtr());
        return p;
      }
    }
    else //if(numM==5)
    {
      if (numX == 1)
      {
        DenseMat p = prod(*x.getDensePtr(), *m.getBandedPtr());
        return p;
      }

      else if (numX == 2)
      {
        DenseMat p = prod(*x.getTriangPtr(), *m.getBandedPtr());
        return p;
      }

      else if (numX == 3)
      {
        DenseMat p = prod(*x.getSymPtr(), *m.getBandedPtr());
        return p;
      }

      else if (numX == 4)
      {
        DenseMat p = prod(*x.getSparsePtr(), *m.getBandedPtr());
        return p;
      }

      else //if(numX==5)
      {
        DenseMat p = prod(*x.getBandedPtr(), *m.getBandedPtr());
        return p;
      }
    }
  }
}

void prod(const SiconosMatrix& A, const SiconosMatrix& B, SiconosMatrix& C)
{
  if ((A.size(1) != B.size(0)))
    SiconosMatrixException::selfThrow("Matrix function prod(A,B,C) : inconsistent sizes");

  if (A.size(0) != C.size(0) || B.size(1) != C.size(1))
    SiconosMatrixException::selfThrow("Matrix function prod(A,B,C) : inconsistent sizes");

  if (A.isBlock() || B.isBlock() || C.isBlock())
    SiconosMatrixException::selfThrow("Matrix function prod(A,B,C) : not yet implemented for block matrices.");

  unsigned int numA = A.getNum();
  unsigned int numB = B.getNum();
  unsigned int numC = C.getNum();

  if (&C == &A || &C == &B) // if common memory ...
  {
    if (numA == 7) // if A is identity
    {
      if (&C != &B) C = B;
    }
    else if (numB == 7) // if B is identity
    {
      if (&C != &A) C = A;
    }
    else if (numA == 6 || numB == 6) // if A or B is null
      C.zero();
    else
    {
      switch (numC)
      {
      case 1:
        if (numB == 1)
        {
          if (numA == 1)
            *C.getDensePtr() = prod(*A.getDensePtr(), *B.getDensePtr());
          else if (numA == 2)
            *C.getDensePtr() = prod(*A.getTriangPtr(), *B.getDensePtr());
          else if (numA == 3)
            *C.getDensePtr()  = prod(*A.getSymPtr(), *B.getDensePtr());
          else if (numA == 4)
            *C.getDensePtr() = prod(*A.getSparsePtr(), *B.getDensePtr());
          else// if(numA==5)
            *C.getDensePtr()  = prod(*A.getBandedPtr(), *B.getDensePtr());
        }
        else if (numB == 2)
        {
          if (numA == 1)
            *C.getDensePtr()  = prod(*A.getDensePtr(), *B.getTriangPtr());
          else if (numA == 2)
            *C.getDensePtr()  = prod(*A.getTriangPtr(), *B.getTriangPtr());
          else if (numA == 3)
            *C.getDensePtr()  = prod(*A.getSymPtr(), *B.getTriangPtr());
          else if (numA == 4)
            *C.getDensePtr()  = prod(*A.getSparsePtr(), *B.getTriangPtr());
          else //if(numA==5)
            *C.getDensePtr()  = prod(*A.getBandedPtr(), *B.getTriangPtr());
        }
        else if (numB == 3)
        {
          if (numA == 1)
            *C.getDensePtr()  = prod(*A.getDensePtr(), *B.getSymPtr());
          else if (numA == 2)
            *C.getDensePtr()  = prod(*A.getTriangPtr(), *B.getSymPtr());
          else if (numA == 3)
            *C.getDensePtr()  = prod(*A.getSymPtr(), *B.getSymPtr());
          else if (numA == 4)
            *C.getDensePtr()  = prod(*A.getSparsePtr(), *B.getSymPtr());
          else // if (numA == 5)
            *C.getDensePtr()  = prod(*A.getBandedPtr(), *B.getSymPtr());
        }
        else if (numB == 4)
        {
          if (numA == 1)
            *C.getDensePtr() = prod(*A.getDensePtr(), *B.getSparsePtr());
          else if (numA == 2)
            *C.getDensePtr() = prod(*A.getTriangPtr(), *B.getSparsePtr());
          else if (numA == 3)
            *C.getDensePtr() = prod(*A.getSymPtr(), *B.getSparsePtr());
          else if (numA == 4)
            *C.getDensePtr() = prod(*A.getSparsePtr(), *B.getSparsePtr());
          else //if(numA==5){
            *C.getDensePtr() = prod(*A.getBandedPtr(), *B.getSparsePtr());
        }
        else //if(numB==5)
        {
          if (numA == 1)
            *C.getDensePtr() = prod(*A.getDensePtr(), *B.getBandedPtr());
          else if (numA == 2)
            *C.getDensePtr() = prod(*A.getTriangPtr(), *B.getBandedPtr());
          else if (numA == 3)
            *C.getDensePtr() = prod(*A.getSymPtr(), *B.getBandedPtr());
          else if (numA == 4)
            *C.getDensePtr() = prod(*A.getSparsePtr(), *B.getBandedPtr());
          else //if(numA==5)
            *C.getDensePtr() = prod(*A.getBandedPtr(), *B.getBandedPtr());
        }
        break;
      case 2:
        if (numA != 2 || numB != 2)
          SiconosMatrixException::selfThrow("Matrix function prod(A,B,C) : wrong type for C (according to A and B types).");
        *C.getTriangPtr() = prod(*A.getTriangPtr(), *B.getTriangPtr());
        break;
      case 3:
        if (numA != 3 || numB != 3)
          SiconosMatrixException::selfThrow("Matrix function prod(A,B,C) : wrong type for C (according to A and B types).");
        *C.getSymPtr() = prod(*A.getSymPtr(), *B.getSymPtr());
        break;
      case 4:
        if (numA != 4 || numB != 4)
          SiconosMatrixException::selfThrow("Matrix function prod(A,B,C) : wrong type for C (according to A and B types).");
        *C.getSparsePtr() = prod(*A.getSparsePtr(), *B.getSparsePtr());
        break;
      default:
        SiconosMatrixException::selfThrow("Matrix function prod(A,B,C) : wrong type for C (according to A and B types).");
      }
    }
  }
  else // if no alias between C and A or B.
  {
    if (numA == 7) // if A is identity
      C = B;
    else if (numB == 7) // if B is identity
      C = A;
    else if (numA == 6 || numB == 6) // if A or B is null
      C.zero();
    else
    {
      switch (numC)
      {
      case 1:
        if (numB == 1)
        {
          if (numA == 1)
            noalias(*C.getDensePtr()) = prod(*A.getDensePtr(), *B.getDensePtr());
          else if (numA == 2)
            noalias(*C.getDensePtr()) = prod(*A.getTriangPtr(), *B.getDensePtr());
          else if (numA == 3)
            noalias(*C.getDensePtr())  = prod(*A.getSymPtr(), *B.getDensePtr());
          else if (numA == 4)
            noalias(*C.getDensePtr()) = prod(*A.getSparsePtr(), *B.getDensePtr());
          else// if(numA==5)
            noalias(*C.getDensePtr())  = prod(*A.getBandedPtr(), *B.getDensePtr());
        }
        else if (numB == 2)
        {
          if (numA == 1)
            noalias(*C.getDensePtr())  = prod(*A.getDensePtr(), *B.getTriangPtr());
          else if (numA == 2)
            noalias(*C.getDensePtr())  = prod(*A.getTriangPtr(), *B.getTriangPtr());
          else if (numA == 3)
            noalias(*C.getDensePtr())  = prod(*A.getSymPtr(), *B.getTriangPtr());
          else if (numA == 4)
            noalias(*C.getDensePtr())  = prod(*A.getSparsePtr(), *B.getTriangPtr());
          else //if(numA==5)
            noalias(*C.getDensePtr())  = prod(*A.getBandedPtr(), *B.getTriangPtr());
        }
        else if (numB == 3)
        {
          if (numA == 1)
            noalias(*C.getDensePtr())  = prod(*A.getDensePtr(), *B.getSymPtr());
          else if (numA == 2)
            noalias(*C.getDensePtr())  = prod(*A.getTriangPtr(), *B.getSymPtr());
          else if (numA == 3)
            noalias(*C.getDensePtr())  = prod(*A.getSymPtr(), *B.getSymPtr());
          else if (numA == 4)
            noalias(*C.getDensePtr())  = prod(*A.getSparsePtr(), *B.getSymPtr());
          else // if (numA == 5)
            noalias(*C.getDensePtr())  = prod(*A.getBandedPtr(), *B.getSymPtr());
        }
        else if (numB == 4)
        {
          if (numA == 1)
            noalias(*C.getDensePtr()) = prod(*A.getDensePtr(), *B.getSparsePtr());
          else if (numA == 2)
            noalias(*C.getDensePtr()) = prod(*A.getTriangPtr(), *B.getSparsePtr());
          else if (numA == 3)
            noalias(*C.getDensePtr()) = prod(*A.getSymPtr(), *B.getSparsePtr());
          else if (numA == 4)
            noalias(*C.getDensePtr()) = prod(*A.getSparsePtr(), *B.getSparsePtr());
          else //if(numA==5){
            noalias(*C.getDensePtr()) = prod(*A.getBandedPtr(), *B.getSparsePtr());
        }
        else //if(numB==5)
        {
          if (numA == 1)
            noalias(*C.getDensePtr()) = prod(*A.getDensePtr(), *B.getBandedPtr());
          else if (numA == 2)
            noalias(*C.getDensePtr()) = prod(*A.getTriangPtr(), *B.getBandedPtr());
          else if (numA == 3)
            noalias(*C.getDensePtr()) = prod(*A.getSymPtr(), *B.getBandedPtr());
          else if (numA == 4)
            noalias(*C.getDensePtr()) = prod(*A.getSparsePtr(), *B.getBandedPtr());
          else //if(numA==5)
            noalias(*C.getDensePtr()) = prod(*A.getBandedPtr(), *B.getBandedPtr());
        }
        break;
      case 2:
        if (numA != 2 || numB != 2)
          SiconosMatrixException::selfThrow("Matrix function prod(A,B,C) : wrong type for C (according to A and B types).");
        noalias(*C.getTriangPtr()) = prod(*A.getTriangPtr(), *B.getTriangPtr());
        break;
      case 3:
        if (numA != 3 || numB != 3)
          SiconosMatrixException::selfThrow("Matrix function prod(A,B,C) : wrong type for C (according to A and B types).");
        noalias(*C.getSymPtr()) = prod(*A.getSymPtr(), *B.getSymPtr());
        break;
      case 4:
        if (numA != 4 || numB != 4)
          SiconosMatrixException::selfThrow("Matrix function prod(A,B,C) : wrong type for C (according to A and B types).");
        noalias(*C.getSparsePtr()) = prod(*A.getSparsePtr(), *B.getSparsePtr());
        break;
      default:
        SiconosMatrixException::selfThrow("Matrix function prod(A,B,C) : wrong type for C (according to A and B types).");
      }
    }
  }
}

const SimpleMatrix operator * (const SiconosMatrix &m, double d)
{
  if (m.isBlock())
    SiconosMatrixException::selfThrow("Matrix, operator * (m*scalar): not yet implemented for BlockMatrix.");
  unsigned int num = m.getNum();

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
  else if (num == 6)
    return m; // ie zero
  else //if(num==7)
  {
    DenseMat p(*m.getIdentityPtr());
    p *= d;
    return p;
  }
}

const SimpleMatrix operator * (const SiconosMatrix &m, int d)
{
  if (m.isBlock())
    SiconosMatrixException::selfThrow("Matrix, operator * (m*scalar): not yet implemented for BlockMatrix.");

  unsigned int num = m.getNum();

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
  else if (num == 6)
    return m; // ie zero
  else //if(num==7)
  {
    DenseMat p(*m.getIdentityPtr());
    p *= d;
    return p;
  }
}

const SimpleMatrix operator * (double d, const SiconosMatrix &m)
{
  if (m.isBlock())
    SiconosMatrixException::selfThrow("Matrix, operator * (m*scalar): not yet implemented for BlockMatrix.");
  unsigned int num = m.getNum();

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
  else if (num == 6)
    return m; // ie zero
  else //if(num==7)
  {
    DenseMat p(*m.getIdentityPtr());
    p *= d;
    return p;
  }
}

const SimpleMatrix operator * (int d, const SiconosMatrix &m)
{
  if (m.isBlock())
    SiconosMatrixException::selfThrow("Matrix, operator * (m*scalar): not yet implemented for BlockMatrix.");
  unsigned int num = m.getNum();

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
  else if (num == 6)
    return m; // ie zero
  else //if(num==7)
  {
    DenseMat p(*m.getIdentityPtr());
    p *= d;
    return p;
  }
}

const SimpleMatrix operator / (const SiconosMatrix &m, double d)
{
  if (m.isBlock())
    SiconosMatrixException::selfThrow("Matrix, operator / (m/scalar): not yet implemented for BlockMatrix.");
  if (d == 0)
    SiconosMatrixException::selfThrow("SimpleMatrix:operator /, division by zero.");

  unsigned int num = m.getNum();
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
  else if (num == 6)
    return m; // ie zero
  else //if(num==7)
  {
    DenseMat p(*m.getIdentityPtr());
    p /= d;
    return p;
  }
}

const SimpleMatrix operator / (const SiconosMatrix &m, int d)
{
  if (m.isBlock())
    SiconosMatrixException::selfThrow("Matrix, operator / (m/scalar): not yet implemented for BlockMatrix.");
  if (d == 0)
    SiconosMatrixException::selfThrow("SimpleMatrix:operator /, division by zero.");

  unsigned int num = m.getNum();

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
  else if (num == 6)
    return m; // ie zero
  else //if(num==7)
  {
    DenseMat p(*m.getIdentityPtr());
    p /= d;
    return p;
  }
}

const SimpleMatrix pow(const SimpleMatrix& m, unsigned int power)
{
  if (m.isBlock())
    SiconosMatrixException::selfThrow("Matrix, pow function: not yet implemented for BlockMatrix.");
  if (!m.isSquare())
    SiconosMatrixException::selfThrow("pow(SimpleMatrix), matrix is not square.");

  if (power < 0)
    SiconosMatrixException::selfThrow("pow(SimpleMatrix,n) with negative value is not supported");

  if (power > 0)
  {
    unsigned int num = m.getNum();
    if (num == 1)
    {
      DenseMat p = *m.getDensePtr();
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
      ZeroMat z(m.size(0), m.size(1));
      return z;
    }
    else // if (num==7)
    {
      IdentityMat I(m.size(0), m.size(1));;
      return I;
    }
  }
  else// if(power == 0)
  {
    IdentityMat I = ublas::identity_matrix<double>(m.size(0), m.size(1));
    return I;
  }
}

// ========== Products matrix - vector ==========
//
// Computation of y = A*x
//
// Two specific functions are used to handle all the cases where x or y are blocks.
// All of their blocks can also be blocks. Then we use:
// - private_prod to "slice" A when y is block, ie according to its rows.
// - private_addprod to "slice" A when x is block, and to sum over the columns of blocks to compute y = sum subA x[i].

// The following function is private and used inside prod(...) public functions.
// It is required to deal with block vectors of blocks ( of blocks ...).
// It computes res = subA*x +res, subA being a submatrix of A (rows from startRow to startRow+sizeY and columns between startCol and startCol+sizeX).
// If x is a block vector, it call the present function for all blocks.
void private_addprod(const SiconosMatrix *A, unsigned int startRow, unsigned int startCol, const SiconosVector *x, DenseVect& res)
{
  if (A->isBlock())
    SiconosMatrixException::selfThrow("private_prod(A,start,x,y) error: not yet implemented for block matrix.");

  if (!x->isBlock()) // if input vector is not block
  {
    // we take a submatrix subA of A, starting from row startRow to row (startRow+sizeY) and between columns startCol and (startCol+sizeX).
    // Then computation of res = subA*x + res.
    unsigned int numA = A->getNum();
    unsigned int sizeX = x->size();
    unsigned int sizeY = res.size();
    if (numA == 1)
      res += prod(ublas::subrange(*A->getDensePtr(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->getDensePtr());
    else if (numA == 2)
      res += prod(ublas::subrange(*A->getTriangPtr(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->getDensePtr());
    else if (numA == 3)
      res += prod(ublas::subrange(*A->getSymPtr(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->getDensePtr());
    else if (numA == 4)
      res += prod(ublas::subrange(*A->getSparsePtr(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->getDensePtr());
    else //if(numA==5)
      res += prod(ublas::subrange(*A->getBandedPtr(), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->getDensePtr());
  }
  else // if block
  {
    ConstBlockVectIterator it;
    unsigned int startColBis = startCol;
    for (it = x->begin(); it != x->end(); ++it)
    {
      private_addprod(A, startRow, startColBis, (*it), res);
      startColBis += (*it)->size();
    }
  }
}

// x and y blocks
void private_prod(const SiconosMatrix *A, unsigned int startRow, const SiconosVector *x, SiconosVector * y)
{

  if (!y->isBlock()) // if y is not a block vector, private_addProd, to consider cases where x is block.
  {
    private_addprod(A, startRow, 0, x, *(y->getDensePtr()));
  }
  else // if y is a block: call private_prod on each block and so on until the considered block is a simple vector.
  {
    unsigned int row = startRow;
    ConstBlockVectIterator it;
    for (it = y->begin(); it != y->end(); ++it)
    {
      private_prod(A, row, x, *it);
      row += (*it)->size();
    }
  }
}

const SimpleVector prod(const SiconosMatrix& A, const SiconosVector& x)
{
  if (A.isBlock())
    SiconosMatrixException::selfThrow("prod(matrix,vector) error: not yet implemented for block matrix.");

  if (A.size(1) != x.size())
    SiconosMatrixException::selfThrow("prod(matrix,vector) error: inconsistent sizes.");

  unsigned int numA = A.getNum();

  if (numA == 6)
  {
    DenseVect res = ublas::zero_vector<double>(x.size());
    return res;
  }
  else if (numA == 7)
    return x;
  else
  {
    if (!x.isBlock()) // if x is not a block vector.
    {
      unsigned int numX = x.getNum();
      if (numX == 1)
      {
        if (numA == 1)
        {
          DenseVect res = prod(*A.getDensePtr(), *x.getDensePtr());
          return res;
        }
        else if (numA == 2)
        {
          DenseVect res = prod(*A.getTriangPtr(), *x.getDensePtr());
          return res;
        }
        else if (numA == 3)
        {
          DenseVect res = prod(*A.getSymPtr(), *x.getDensePtr());
          return res;
        }
        else if (numA == 4)
        {
          DenseVect res = prod(*A.getSparsePtr(), *x.getDensePtr());
          return res;
        }
        else if (numA == 5)
        {
          DenseVect res = prod(*A.getBandedPtr(), *x.getDensePtr());
          return res;
        }
        else if (numA == 6)
        {
          DenseVect res = ublas::zero_vector<double>(x.size());
          return res;
        }
        else //if(numA==7)
        {
          DenseVect  res = *x.getDensePtr();
          return res;
        }
      }
      else //if(numX == 4)
      {
        if (numA == 1)
        {
          DenseVect res = prod(*A.getDensePtr(), *x.getSparsePtr());
          return res;
        }
        else if (numA == 2)
        {
          DenseVect res = prod(*A.getTriangPtr(), *x.getSparsePtr());
          return res;
        }
        else if (numA == 3)
        {
          DenseVect res = prod(*A.getSymPtr(), *x.getSparsePtr());
          return res;
        }
        else if (numA == 4)
        {
          DenseVect res = prod(*A.getSparsePtr(), *x.getSparsePtr());
          return res;
        }
        else if (numA == 5)
        {
          DenseVect res = prod(*A.getBandedPtr(), *x.getSparsePtr());
          return res;
        }
        else if (numA == 6)
        {
          DenseVect res = ublas::zero_vector<double>(x.size());
          return res;
        }
        else //if(numA==7)
        {
          DenseVect res = *x.getSparsePtr();
          return res;
        }
      }
    }
    else // if (x.isBlock())
    {
      ConstBlockVectIterator it;
      DenseVect res = ublas::zero_vector<double>(A.size(0));
      unsigned int start = 0;
      for (it = x.begin(); it != x.end(); ++it)
      {
        private_addprod(&A, 0, start, *it, res);
        start += (*it)->size();
      }
      return res;
    }
  }
}

void prod(const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y)
{
  if (A.isBlock())
    SiconosMatrixException::selfThrow("prod(A,x,y) error: not yet implemented for block matrices.");

  if (A.size(1) != x.size())
    SiconosMatrixException::selfThrow("prod(A,x,y) error: inconsistent sizes between A and x.");

  if (A.size(0) != y.size())
    SiconosMatrixException::selfThrow("prod(A,x,y) error: inconsistent sizes between A and y.");
  unsigned int numA = A.getNum();

  if (numA == 6) // A = 0
    y.zero();
  else if (numA == 7) // A = I
  {
    if (&x != &y) y = x;
  }
  else
  {
    // === First case: y is not a block vector ===
    if (!y.isBlock())
    {
      if (y.getNum() != 1)
        SiconosMatrixException::selfThrow("prod(A,x,y) error: y (output) must be a dense vector.");

      // if x is block: call of a specific function to treat each block
      if (x.isBlock())
      {
        y.zero();
        unsigned int startRow = 0;
        unsigned int startCol = 0;
        // In private_addprod, the sum of all blocks of x, x[i], is computed : y = Sum_i (subA x[i]), with subA a submatrix of A,
        // starting from position startRow in rows and startCol in columns.
        // private_prod takes also into account the fact that each block of x can also be a block.
        ConstBlockVectIterator it;
        for (it = x.begin(); it != x.end(); ++it)
        {
          private_addprod(&A, startRow, startCol, *it, *y.getDensePtr());
          startCol += (*it)->size();
        }
      }
      else // If neither x nor y are block: direct call to ublas::prod.
      {
        unsigned int numX = x.getNum();
        if (&x != &y) // if no common memory between x and y.
        {
          if (numX == 1)
          {
            if (numA == 1)
              noalias(*y.getDensePtr()) = ublas::prod(*A.getDensePtr(), *x.getDensePtr());
            else if (numA == 2)
              noalias(*y.getDensePtr()) = ublas::prod(*A.getTriangPtr(), *x.getDensePtr());
            else if (numA == 3)
              noalias(*y.getDensePtr()) = ublas::prod(*A.getSymPtr(), *x.getDensePtr());
            else if (numA == 4)
              noalias(*y.getDensePtr()) = ublas::prod(*A.getSparsePtr(), *x.getDensePtr());
            else //if(numA==5)
              noalias(*y.getDensePtr()) = ublas::prod(*A.getBandedPtr(), *x.getDensePtr());
          }
          else //if(numX == 4)
          {
            if (numA == 1)
              noalias(*y.getDensePtr()) = ublas::prod(*A.getDensePtr(), *x.getSparsePtr());
            else if (numA == 2)
              noalias(*y.getDensePtr()) = ublas::prod(*A.getTriangPtr(), *x.getSparsePtr());
            else if (numA == 3)
              noalias(*y.getDensePtr()) = ublas::prod(*A.getSymPtr(), *x.getSparsePtr());
            else if (numA == 4)
              noalias(*y.getDensePtr()) = ublas::prod(*A.getSparsePtr(), *x.getSparsePtr());
            else //if(numA==5)
              noalias(*y.getDensePtr()) = ublas::prod(*A.getBandedPtr(), *x.getSparsePtr());
          }
        }
        else
        {
          if (numX == 1)
          {
            if (numA == 1)
              *y.getDensePtr() = ublas::prod(*A.getDensePtr(), *x.getDensePtr());
            else if (numA == 2)
              *y.getDensePtr() = ublas::prod(*A.getTriangPtr(), *x.getDensePtr());
            else if (numA == 3)
              *y.getDensePtr() = ublas::prod(*A.getSymPtr(), *x.getDensePtr());
            else if (numA == 4)
              *y.getDensePtr() = ublas::prod(*A.getSparsePtr(), *x.getDensePtr());
            else //if(numA==5)
              *y.getDensePtr() = ublas::prod(*A.getBandedPtr(), *x.getDensePtr());
          }
          else //if(numX == 4)
          {
            if (numA == 1)
              *y.getDensePtr() = ublas::prod(*A.getDensePtr(), *x.getSparsePtr());
            else if (numA == 2)
              *y.getDensePtr() = ublas::prod(*A.getTriangPtr(), *x.getSparsePtr());
            else if (numA == 3)
              *y.getDensePtr() = ublas::prod(*A.getSymPtr(), *x.getSparsePtr());
            else if (numA == 4)
              *y.getDensePtr() = ublas::prod(*A.getSparsePtr(), *x.getSparsePtr());
            else //if(numA==5)
              *y.getDensePtr() = ublas::prod(*A.getBandedPtr(), *x.getSparsePtr());
          }
        }
      }
    }
    else // === Second case: y is a block vector ===
    {
      unsigned int startRow = 0;
      ConstBlockVectIterator it;
      // For Each subvector of y, y[i], private_prod computes y[i] = subA x, subA being a submatrix of A corresponding to y[i] position.
      // private_prod takes into account the fact that x and y[i] may be block vectors.
      for (it = y.begin(); it != y.end(); ++it)
      {
        private_prod(&A, startRow, &x, *it);
        startRow += (*it)->size();
      }
    }
  }
}

void gemv(const CBLAS_TRANSPOSE transA, double a, const SiconosMatrix& A, const SiconosVector& x, double b, SiconosVector& y)
{
  if (x.isBlock() || y.isBlock() || A.isBlock())
    SiconosMatrixException::selfThrow("gemv(...) not yet implemented for block vectors or matrices.");

  unsigned int numA = A.getNum();
  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();
  if (numA != 1 || numX != 1 || numY != 1)
    SiconosMatrixException::selfThrow("gemv(...) failed: reserved to dense matrices or vectors.");

  atlas::gemv(transA, a, *A.getDensePtr(), *x.getDensePtr(), b, *y.getDensePtr());
}

void gemv(double a, const SiconosMatrix& A, const SiconosVector& x, double b, SiconosVector& y)
{
  if (x.isBlock() || y.isBlock() || A.isBlock())
    SiconosMatrixException::selfThrow("gemv(...) not yet implemented for block vectors or matrices.");
  unsigned int numA = A.getNum();
  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();
  if (numA != 1 || numX != 1 || numY != 1)
    SiconosMatrixException::selfThrow("gemv(...) failed: reserved to dense matrices or vectors.");

  atlas::gemv(a, *A.getDensePtr(), *x.getDensePtr(), b, *y.getDensePtr());
}

void gemv(const SiconosMatrix& A, const SiconosVector& x, SiconosVector& y)
{
  if (x.isBlock() || y.isBlock() || A.isBlock())
    SiconosMatrixException::selfThrow("gemv(...) not yet implemented for block vectors or matrices.");
  unsigned int numA = A.getNum();
  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();
  if (numA != 1 || numX != 1 || numY != 1)
    SiconosMatrixException::selfThrow("gemv(...) failed: reserved to dense matrices or vectors.");

  atlas::gemv(*A.getDensePtr(), *x.getDensePtr(), *y.getDensePtr());
}

void gemm(const CBLAS_TRANSPOSE transA, const CBLAS_TRANSPOSE transB, double a, const SiconosMatrix& A, const SiconosMatrix& B, double b, SiconosMatrix& C)
{
  if (A.isBlock() || B.isBlock() || C.isBlock())
    SiconosMatrixException::selfThrow("gemm(...) not yet implemented for block matrices.");
  unsigned int numA = A.getNum();
  unsigned int numB = B.getNum();
  unsigned int numC = C.getNum();
  if (numA != 1 || numB != 1 || numC != 1)
    SiconosMatrixException::selfThrow("gemm(...) failed: reserved to dense matrices.");

  atlas::gemm(transA, transB, a, *A.getDensePtr(), *B.getDensePtr(), b, *C.getDensePtr());
}

void gemm(double a, const SiconosMatrix& A, const SiconosMatrix& B, double b, SiconosMatrix& C)
{
  if (A.isBlock() || B.isBlock() || C.isBlock())
    SiconosMatrixException::selfThrow("gemm(...) not yet implemented for block matrices.");
  unsigned int numA = A.getNum();
  unsigned int numB = B.getNum();
  unsigned int numC = C.getNum();
  if (numA != 1 || numB != 1 || numC != 1)
    SiconosMatrixException::selfThrow("gemm(...) failed: reserved to dense matrices.");

  atlas::gemm(a, *A.getDensePtr(), *B.getDensePtr(), b, *C.getDensePtr());
}

void gemm(const SiconosMatrix& A, const SiconosMatrix& B, SiconosMatrix& C)
{
  if (A.isBlock() || B.isBlock() || C.isBlock())
    SiconosMatrixException::selfThrow("gemm(...) not yet implemented for block matrices.");
  unsigned int numA = A.getNum();
  unsigned int numB = B.getNum();
  unsigned int numC = C.getNum();
  if (numA != 1 || numB != 1 || numC != 1)
    SiconosMatrixException::selfThrow("gemm(...) failed: reserved to dense matrices.");

  atlas::gemm(*A.getDensePtr(), *B.getDensePtr(), *C.getDensePtr());
}
