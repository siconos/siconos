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

void SimpleMatrix::checkNum()
{
  if (num != 1 || num != 2 || num != 3 || num != 4 || num != 5)
    SiconosMatrixException::selfThrow("SimpleMatrix error: bad value for num (invalid matrix type)");
  // Remark: this function is usefull to avoid "warning: control reaches end of non-void function" in many functions of this class.
}

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
  else
    SiconosMatrixException::selfThrow("constructor(TYP) : invalid type given");
  zero();
}

// Copy constructors
SimpleMatrix::SimpleMatrix(const SimpleMatrix &smat): SiconosMatrix(false), num(smat.getNum()), isPLUFactorized(false), isPLUInversed(false)
{
  if (num == 1)
    mat.Dense = new DenseMat(smat.getDense());

  else if (num == 2)
    mat.Triang = new TriangMat(smat.getTriang());

  else if (num == 3)
    mat.Sym = new SymMat(smat.getSym());

  else if (smat.getNum() == 4)
    mat.Sparse = new SparseMat(smat.getSparse());

  else if (smat.getNum() == 5)
    mat.Banded = new BandedMat(smat.getBanded());

  else
    SiconosMatrixException::selfThrow("constructor(const SimpleMatrix) : invalid parameter given");
}

SimpleMatrix::SimpleMatrix(const SiconosMatrix &smat): SiconosMatrix(false), num(smat.getNum()), isPLUFactorized(false), isPLUInversed(false)
{
  assert(smat.isBlock() == false);
  if (num == 1)
    mat.Dense = new DenseMat(smat.getDense());

  else if (num == 2)
    mat.Triang = new TriangMat(smat.getTriang());

  else if (num == 3)
    mat.Sym = new SymMat(smat.getSym());

  else if (num == 4)
    mat.Sparse = new SparseMat(smat.getSparse());

  else if (num == 5)
    mat.Banded = new BandedMat(smat.getBanded());

  else
    SiconosMatrixException::selfThrow("constructor(const SiconosMatrix) : invalid parameter given");
}

SimpleMatrix::SimpleMatrix(unsigned int row, unsigned int col, TYP typ, unsigned int upper, unsigned int lower):
  SiconosMatrix(false), num(1), isPLUFactorized(false), isPLUInversed(false)
{
  if (typ == DENSE)
  {
    mat.Dense = new DenseMat(row, col);
    num = 1;
  }
  else if (typ == TRIANGULAR)
  {
    mat.Triang = new TriangMat(row, col);
    num = 2;
  }
  else if (typ == SYMMETRIC)
  {
    mat.Sym = new SymMat(row, col);
    num = 3;
  }
  else if (typ == SPARSE)
  {
    mat.Sparse = new SparseMat(row, col, upper);
    num = 4;
  }
  else if (typ == BANDED)
  {
    mat.Banded = new BandedMat(row, col, upper, lower);
    num = 5;
  }
  else
    SiconosMatrixException::selfThrow("constructor(TYP type, unsigned int row, unsigned int col) : invalid type or dimensions given");
  zero();
}

SimpleMatrix::SimpleMatrix(const DenseMat& m): SiconosMatrix(false), num(1), isPLUFactorized(false), isPLUInversed(false)
{
  mat.Dense = new DenseMat(m);
}

SimpleMatrix::SimpleMatrix(const TriangMat& m): SiconosMatrix(false), num(2), isPLUFactorized(false), isPLUInversed(false)
{
  mat.Triang = new TriangMat(m);
}

SimpleMatrix::SimpleMatrix(const SymMat& m): SiconosMatrix(false), num(3), isPLUFactorized(false), isPLUInversed(false)
{
  mat.Sym = new SymMat(m);
}

SimpleMatrix::SimpleMatrix(const SparseMat& m): SiconosMatrix(false), num(4), isPLUFactorized(false), isPLUInversed(false)
{
  mat.Sparse = new SparseMat(m);
}

SimpleMatrix::SimpleMatrix(const BandedMat& m): SiconosMatrix(false), num(5), isPLUFactorized(false), isPLUInversed(false)
{
  mat.Banded = new BandedMat(m);
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
  }
  else if (typ == TRIANGULAR)
  {
    mat.Triang = new TriangMat(row, col, v);
    num = 2;
  }
  else if (typ == SYMMETRIC)
  {
    mat.Sym = new SymMat(row, v);
    num = 3;
  }
  else if (typ == SPARSE)
  {
    SiconosMatrixException::selfThrow("constructor(TYP, const std::vector<double>, int row, int col, int lower, int upper) : warning -- use constructor(const SparseMat &m) or constructor(TYP, int row, int col) with TYP = SPARSE");

  }
  else if (typ == BANDED)
  {
    mat.Banded = new BandedMat(row, col, lower, upper, v);
    num = 5;
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
}

/******************************** METHODS ******************************/
unsigned int  SimpleMatrix::getNum(void)const
{
  return num;
}

unsigned int  SimpleMatrix::size(unsigned int dim) const
{
  unsigned int sizeOut = 0;
  if (dim == 0)
  {
    if (num == 1)
    {
      sizeOut = (*mat.Dense).size1();
    }
    else if (num == 2)
    {
      sizeOut = (*mat.Triang).size1();
    }
    else if (num == 3)
    {
      sizeOut = (*mat.Sym).size1();
    }
    else if (num == 4)
    {
      sizeOut = (*mat.Sparse).size1();
    }
    else if (num == 5)
    {
      sizeOut = (*mat.Banded).size1();
    }
  }
  else
  {
    if (num == 1)
    {
      sizeOut = (*mat.Dense).size2();
    }
    else if (num == 2)
    {
      sizeOut = (*mat.Triang).size2();
    }
    else if (num == 3)
    {
      sizeOut = (*mat.Sym).size2();
    }
    else if (num == 4)
    {
      sizeOut = (*mat.Sparse).size2();
    }
    else if (num == 5)
    {
      sizeOut = (*mat.Banded).size2();
    }
  }
  return sizeOut;
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

  //if(num==1){
  if (num2 == 1)
  {
    unsigned int row_max = i + (m.getDense()).size1();
    unsigned int col_max = j + (m.getDense()).size2();
    if (row_max > size(0) || row_max < 0)
      SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

    if (col_max > size(1) || col_max < 0)
      SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

    ublas::subrange(*mat.Dense, i, i + (m.getDense()).size1(), j, j + (m.getDense()).size2()) = m.getDense();
  }
  else if (num2 == 2)
  {
    unsigned int row_max = i + (m.getTriang()).size1();
    unsigned int col_max = j + (m.getTriang()).size2();
    if (row_max > size(0) || row_max < 0)
      SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

    if (col_max > size(1) || col_max < 0)
      SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

    ublas::subrange(*mat.Dense, i, i + (m.getTriang()).size1(), j, j + (m.getTriang()).size2()) = m.getTriang();
  }
  else if (num2 == 3)
  {
    unsigned int row_max = i + (m.getSym()).size1();
    unsigned int col_max = j + (m.getSym()).size2();
    if (row_max > size(0) || row_max < 0)
      SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

    if (col_max > size(1) || col_max < 0)
      SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

    ublas::subrange(*mat.Dense, i, i + (m.getSym()).size1(), j, j + (m.getSym()).size2()) = m.getSym();
  }
  else if (num2 == 4)
  {
    unsigned int row_max = i + (m.getSparse()).size1();
    unsigned int col_max = j + (m.getSparse()).size2();
    if (row_max > size(0) || row_max < 0)
      SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

    if (col_max > size(1) || col_max < 0)
      SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

    ublas::subrange(*mat.Dense, i, i + (m.getSparse()).size1(), j, j + (m.getSparse()).size2()) = m.getSparse();
  }
  else if (num2 == 5)
  {
    unsigned int row_max = i + (m.getBanded()).size1();
    unsigned int col_max = j + (m.getBanded()).size2();
    if (row_max > size(0) || row_max < 0)
      SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

    if (col_max > size(1) || col_max < 0)
      SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

    ublas::subrange(*mat.Dense, i, i + (m.getBanded()).size1(), j, j + (m.getBanded()).size2()) = m.getBanded();
  }
  //  }
  //   else if(num==2){
  //     if(num2==1){
  //       int row_max = i + (m.getDense()).size1();
  //       int col_max = j + (m.getDense()).size2();
  //       if(row_max > size(0) || row_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       if(col_max > size(1) || col_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       ublas::subrange(*mat.Triang, i, i+(m.getDense ()).size1 (), j, j+(m.getDense ()).size2 ()) = m.getDense ();
  //     }
  //     else if(num2==2){
  //       int row_max = i + (m.getTriang()).size1();
  //       int col_max = j + (m.getTriang()).size2();
  //       if(row_max > size(0) || row_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       if(col_max > size(1) || col_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       ublas::subrange(*mat.Triang, i, i+(m.getTriang ()).size1 (), j, j+(m.getTriang ()).size2 ()) = m.getTriang();
  //     }
  //     else if(num2==3){
  //       int row_max = i + (m.getSym()).size1();
  //       int col_max = j + (m.getSym()).size2();
  //       if(row_max > size(0) || row_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       if(col_max > size(1) || col_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       ublas::subrange(*mat.Triang, i, i+(m.getSym ()).size1 (), j, j+(m.getSym ()).size2 ()) = m.getSym ();
  //     }
  //     else if(num2==4){
  //       int row_max = i + (m.getSparse()).size1();
  //       int col_max = j + (m.getSparse()).size2();
  //       if(row_max > size(0) || row_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       if(col_max > size(1) || col_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       ublas::subrange(*mat.Triang, i, i+(m.getSparse ()).size1 (), j, j+(m.getSparse ()).size2 ()) = m.getSparse ();
  //     }
  //     else if(num2==5){
  //       int row_max = i + (m.getBanded()).size1();
  //       int col_max = j + (m.getBanded()).size2();
  //       if(row_max > size(0) || row_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       if(col_max > size(1) || col_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       ublas::subrange(*mat.Triang, i, i+(m.getBanded ()).size1 (), j, j+(m.getBanded ()).size2 ()) = m.getBanded ();
  //     }
  //   }
  //   else if(num==3){
  //     if(num2==1){
  //       int row_max = i + (m.getDense()).size1();
  //       int col_max = j + (m.getDense()).size2();
  //       if(row_max > size(0) || row_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       if(col_max > size(1) || col_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       ublas::subrange(*mat.Sym, i, i+(m.getDense ()).size1 (), j, j+(m.getDense ()).size2 ()) = m.getDense ();
  //     }
  //     else if(num2==2){
  //       int row_max = i + (m.getTriang()).size1();
  //       int col_max = j + (m.getTriang()).size2();
  //       if(row_max > size(0) || row_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       if(col_max > size(1) || col_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       ublas::subrange(*mat.Sym, i, i+(m.getTriang ()).size1 (), j, j+(m.getTriang ()).size2 ()) = m.getTriang ();
  //     }
  //     else if(num2==3){
  //       int row_max = i + (m.getSym()).size1();
  //       int col_max = j + (m.getSym()).size2();
  //       if(row_max > size(0) || row_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       if(col_max > size(1) || col_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       ublas::subrange(*mat.Sym, i, i+(m.getSym ()).size1 (), j, j+(m.getSym ()).size2 ()) = m.getSym ();
  //     }
  //     else if(num2==4){
  //       int row_max = i + (m.getSparse()).size1();
  //       int col_max = j + (m.getSparse()).size2();
  //       if(row_max > size(0) || row_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       if(col_max > size(1) || col_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       ublas::subrange(*mat.Sym, i, i+(m.getSparse ()).size1 (), j, j+(m.getSparse ()).size2 ()) = m.getSparse ();
  //     }
  //     else if(num2==5){
  //       int row_max = i + (m.getBanded()).size1();
  //       int col_max = j + (m.getBanded()).size2();
  //       if(row_max > size(0) || row_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       if(col_max > size(1) || col_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       ublas::subrange(*mat.Sym, i, i+(m.getBanded ()).size1 (), j, j+(m.getBanded ()).size2 ()) = m.getBanded ();
  //     }
  //   }
  //   if(num==4){
  //     if(num2==1){
  //       int row_max = i + (m.getDense()).size1();
  //       int col_max = j + (m.getDense()).size2();
  //       if(row_max > size(0) || row_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       if(col_max > size(1) || col_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       ublas::subrange(*mat.Sparse, i, i+(m.getDense ()).size1 (), j, j+(m.getDense ()).size2 ()) = m.getDense ();
  //     }
  //     else if(num2==2){
  //       int row_max = i + (m.getTriang()).size1();
  //       int col_max = j + (m.getTriang()).size2();
  //       if(row_max > size(0) || row_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       if(col_max > size(1) || col_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       ublas::subrange(*mat.Sparse, i, i+(m.getTriang ()).size1 (), j, j+(m.getTriang ()).size2 ()) = m.getTriang ();
  //     }
  //     else if(num2==3){
  //       int row_max = i + (m.getSym()).size1();
  //       int col_max = j + (m.getSym()).size2();
  //       if(row_max > size(0) || row_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       if(col_max > size(1) || col_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       ublas::subrange(*mat.Sparse, i, i+(m.getSym ()).size1 (), j, j+(m.getSym ()).size2 ()) = m.getSym ();
  //     }
  //     else if(num2==4){
  //       int row_max = i + (m.getSparse()).size1();
  //       int col_max = j + (m.getSparse()).size2();
  //       if(row_max > size(0) || row_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       if(col_max > size(1) || col_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       ublas::subrange(*mat.Sparse, i, i+(m.getSparse ()).size1 (), j, j+(m.getSparse ()).size2 ()) = m.getSparse ();
  //     }
  //     else if(num2==5){
  //       int row_max = i + (m.getBanded()).size1();
  //       int col_max = j + (m.getBanded()).size2();
  //       if(row_max > size(0) || row_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       if(col_max > size(1) || col_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       ublas::subrange(*mat.Sparse, i, i+(m.getBanded ()).size1 (), j, j+(m.getBanded ()).size2 ()) = m.getBanded ();
  //     }
  //   }
  //   if(num==5){
  //     if(num2==1){
  //       int row_max = i + (m.getDense()).size1();
  //       int col_max = j + (m.getDense()).size2();
  //       if(row_max > size(0) || row_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       if(col_max > size(1) || col_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       ublas::subrange(*mat.Banded, i, i+(m.getDense ()).size1 (), j, j+(m.getDense ()).size2 ()) = m.getDense ();
  //     }
  //     else if(num2==2){
  //       int row_max = i + (m.getTriang()).size1();
  //       int col_max = j + (m.getTriang()).size2();
  //       if(row_max > size(0) || row_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       if(col_max > size(1) || col_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       ublas::subrange(*mat.Banded, i, i+(m.getTriang ()).size1 (), j, j+(m.getTriang ()).size2 ()) = m.getTriang ();
  //     }
  //     else if(num2==3){
  //       int row_max = i + (m.getSym()).size1();
  //       int col_max = j + (m.getSym()).size2();
  //       if(row_max > size(0) || row_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       if(col_max > size(1) || col_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       ublas::subrange(*mat.Banded, i, i+(m.getSym ()).size1 (), j, j+(m.getSym ()).size2 ()) = m.getSym ();
  //     }
  //     else if(num2==4){
  //       int row_max = i + (m.getSparse()).size1();
  //       int col_max = j + (m.getSparse()).size2();
  //       if(row_max > size(0) || row_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       if(col_max > size(1) || col_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       ublas::subrange(*mat.Banded, i, i+(m.getSparse ()).size1 (), j, j+(m.getSparse ()).size2 ()) = m.getSparse ();
  //     }
  //     else if(num2==5){
  //       int row_max = i + (m.getBanded()).size1();
  //       int col_max = j + (m.getBanded()).size2();
  //       if(row_max > size(0) || row_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       if(col_max > size(1) || col_max<0)
  //  SiconosMatrixException::selfThrow("SimpleMatrix::matrixCopy : inconsistent sizes");

  //       ublas::subrange(*mat.Banded, i, i+(m.getBanded ()).size1 (), j, j+(m.getBanded ()).size2 ()) = m.getBanded ();
  //     }
  //   }
}

void SimpleMatrix::getBlock(unsigned int row_min, unsigned int col_min, SiconosMatrix &m)const
{

  // We only accept dense matrix for m.
  if (m.getNum() != 1)
    SiconosMatrixException::selfThrow("getBlock(i,j,m) : m must be a dense matrix.");

  if (row_min >= size(0) || row_min < 0)
    SiconosMatrixException::selfThrow("getBlock : row_min given is out of range");

  if (col_min >= size(1) || row_min < 0)
    SiconosMatrixException::selfThrow("getBlock : col_min given is out of range");
  DenseMat q;
  unsigned int row_max, col_max;
  if (num == 1)
  {
    row_max = m.getDense().size1() + row_min;
    col_max = m.getDense().size2() + col_min;

    if (row_max > size(0) || row_max < 0)
      SiconosMatrixException::selfThrow("getBlock : inconsistent sizes");

    if (col_max > size(1) || col_max < 0)
      SiconosMatrixException::selfThrow("getBlock : inconsistent sizes");

    q = ublas::subrange(*mat.Dense, row_min, row_max, col_min, col_max);
  }
  else if (num == 2)
  {
    row_max = m.getDense().size1() + row_min;
    col_max = m.getDense().size2() + col_min;

    if (row_max > size(0) || row_max < 0)
      SiconosMatrixException::selfThrow("getBlock : inconsistent sizes");

    if (col_max > size(1) || col_max < 0)
      SiconosMatrixException::selfThrow("getBlock : inconsistent sizes");
    q = ublas::subrange(*mat.Triang, row_min, row_max, col_min, col_max);
  }
  else if (num == 3)
  {
    row_max = m.getDense().size1() + row_min;
    col_max = m.getDense().size2() + col_min;

    if (row_max > size(0) || row_max < 0)
      SiconosMatrixException::selfThrow("getBlock : inconsistent sizes");

    if (col_max > size(1) || col_max < 0)
      SiconosMatrixException::selfThrow("getBlock : inconsistent sizes");
    q = ublas::subrange(*mat.Sym, row_min, row_max, col_min, col_max);
  }
  else if (num == 4)
  {
    row_max = m.getDense().size1() + row_min;
    col_max = m.getDense().size2() + col_min;

    if (row_max > size(0) || row_max < 0)
      SiconosMatrixException::selfThrow("getBlock : inconsistent sizes");

    if (col_max > size(1) || col_max < 0)
      SiconosMatrixException::selfThrow("getBlock : inconsistent sizes");
    q = ublas::subrange(*mat.Sparse, row_min, row_max, col_min, col_max);
  }
  else if (num == 5)
  {
    row_max = m.getDense().size1() + row_min;
    col_max = m.getDense().size2() + col_min;

    if (row_max > size(0) || row_max < 0)
      SiconosMatrixException::selfThrow("getBlock : inconsistent sizes");

    if (col_max > size(1) || col_max < 0)
      SiconosMatrixException::selfThrow("getBlock : inconsistent sizes");
    q = ublas::subrange(*mat.Banded, row_min, row_max, col_min, col_max);
  }
  SimpleMatrix p(q);
  m = p;
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
}

void SimpleMatrix::setRow(unsigned int r, const SimpleVector &vect)
{

  if (r >= size(0) || r < 0)
    SiconosMatrixException::selfThrow("setRow : row is out of range");

  if (vect.size() != size(1))
    SiconosMatrixException::selfThrow("setRow : inconsistent sizes");

  if (num == 1)
  {
    if (vect.getNum() == 1)
    {
      ublas::row(*mat.Dense, r) = vect.getDense();
    }
    else if (vect.getNum() == 2)
    {
      ublas::row(*mat.Dense, r) = vect.getSparse();
    }
  }
  else if (num == 2)
  {
    if (vect.getNum() == 1)
    {
      ublas::row(*mat.Triang, r) = vect.getDense();
    }
    else if (vect.getNum() == 2)
    {
      ublas::row(*mat.Triang, r) = vect.getSparse();
    }
  }
  else if (num == 3)
  {
    if (vect.getNum() == 1)
    {
      ublas::row(*mat.Sym, r) = vect.getDense();
    }
    else if (vect.getNum() == 2)
    {
      ublas::row(*mat.Sym, r) = vect.getSparse();
    }
  }
  if (num == 4)
  {
    if (vect.getNum() == 1)
    {
      ublas::row(*mat.Sparse, r) = vect.getDense();
    }
    else if (vect.getNum() == 2)
    {
      ublas::row(*mat.Sparse, r) = vect.getSparse();
    }
  }
  if (num == 5)
  {
    if (vect.getNum() == 1)
    {
      ublas::row(*mat.Sparse, r) = vect.getDense();
    }
    else if (vect.getNum() == 2)
    {
      ublas::row(*mat.Banded, r) = vect.getSparse();
    }
  }
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
}

void SimpleMatrix::setCol(unsigned int r, const SimpleVector &vect)
{

  if (r >= size(1) || r < 0)
    SiconosMatrixException::selfThrow("setCol : col is out of range");

  if (vect.size() != size(0))
    SiconosMatrixException::selfThrow("setCol : inconsistent sizes");

  if (num == 1)
  {
    if (vect.getNum() == 1)
    {
      ublas::column(*mat.Dense, r) = vect.getDense();
    }
    else if (vect.getNum() == 2)
    {
      ublas::column(*mat.Dense, r) = vect.getSparse();
    }
  }
  else if (num == 2)
  {
    if (vect.getNum() == 1)
    {
      ublas::column(*mat.Triang, r) = vect.getDense();
    }
    else if (vect.getNum() == 2)
    {
      ublas::column(*mat.Triang, r) = vect.getSparse();
    }
  }
  else if (num == 3)
  {
    if (vect.getNum() == 1)
    {
      ublas::column(*mat.Sym, r) = vect.getDense();
    }
    else if (vect.getNum() == 2)
    {
      ublas::column(*mat.Sym, r) = vect.getSparse();
    }
  }
  else if (num == 4)
  {
    if (vect.getNum() == 1)
    {
      ublas::column(*mat.Sparse, r) = vect.getDense();
    }
    else if (vect.getNum() == 2)
    {
      ublas::column(*mat.Sparse, r) = vect.getSparse();
    }
  }
  else if (num == 5)
  {
    if (vect.getNum() == 1)
    {
      ublas::column(*mat.Banded, r) = vect.getDense();
    }
    else if (vect.getNum() == 2)
    {
      ublas::column(*mat.Banded, r) = vect.getSparse();
    }
  }
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
  return d;
}

SimpleMatrix trans(const SiconosMatrix &m)
{
  if (m.isBlock())
    SiconosMatrixException::selfThrow("transpose not yet implemented for block matrices.");

  if (m.getNum() == 5) // Boost error
    SiconosMatrixException::selfThrow("transpose of a banded matrix not allowed.");

  if (m.getNum() == 1)
  {
    DenseMat p = trans(m.getDense());
    return p;
  }
  else if (m.getNum() == 2)
  {
    SiconosMatrixException::selfThrow("transpose of a triangular matrix not allowed since lower triangular matrix are not implemented.");
    TriangMat p = trans(m.getTriang());
    return p;
  }
  else if (m.getNum() == 3)
  {
    return m;
  }
  else  //if(m.getNum()==4){
  {
    SparseMat p = trans(m.getSparse());
    return p;
  }
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
  else
    d = &(((*mat.Banded).data())[0]);

  return d;
}

void SimpleMatrix::zero(void)
{
  unsigned int size1 = (*this).size(0);
  unsigned int size2 = (*this).size(1);
  if (num == 1)
  {
    ublas::zero_matrix<double> p(size1, size2);
    *mat.Dense = p;
  }
  else if (num == 2)
  {
    ublas::zero_matrix<double> p(size1, size2);
    *mat.Triang = p;
  }
  else if (num == 3)
  {
    ublas::zero_matrix<double> p(size1, size2);
    *mat.Sym = p;
  }
  else if (num == 4)
  {
    ublas::zero_matrix<double> p(size1, size2);
    *mat.Sparse = p;
  }
  else if (num == 5)
  {
    ublas::zero_matrix<double> p(size1, size2);
    *mat.Banded = p;
  }
}

void SimpleMatrix::eye(void)
{
  unsigned int size1 = (*this).size(0);
  unsigned int size2 = (*this).size(1);
  if (num == 1)
  {
    ublas::identity_matrix<double> p(size1, size2);
    *mat.Dense = p;
  }
  else if (num == 2)
  {
    ublas::identity_matrix<double> p(size1, size2);
    *mat.Triang = p;
  }
  else if (num == 3)
  {
    ublas::identity_matrix<double> p(size1, size2);
    *mat.Sym = p;
  }
  else if (num == 4)
  {
    ublas::identity_matrix<double> p(size1, size2);
    *mat.Sparse = p;
  }
  else if (num == 5)
  {
    ublas::identity_matrix<double> p(size1, size2);
    *mat.Banded = p;
  }
}

/***************************** OPERATORS ******************************/

double& SimpleMatrix::operator()(unsigned int row, unsigned int col)
{
  if (row >= size(0) || col >= size(1))
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
  else // if(num==5)
    return (*mat.Banded)(row, col);
}

double SimpleMatrix::operator()(unsigned int row, unsigned int col) const
{
  if (row >= size(0) || col >= size(1))
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
  default:
    SiconosMatrixException::selfThrow("op() (unsigned int, unsigned int) : invalid type of matrix");
    break;
  }
  return d;
}

SimpleMatrix& SimpleMatrix::operator = (const SimpleMatrix& m)
{
  switch (num)
  {
  case 1:
    switch (m.getNum())
    {
    case 1:
      *mat.Dense = m.getDense();
      break;
    case 2:
      *mat.Dense = m.getTriang();
      break;
    case 3:
      *mat.Dense = m.getSym();
      break;
    case 4:
      *mat.Dense = m.getSparse();
      break;
    case 5:
      *mat.Dense = m.getBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("SimpleMatrix::op= (const SimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 2:
    switch (m.getNum())
    {
    case 2:
      *mat.Triang = m.getTriang();
      break;
    default:
      SiconosMatrixException::selfThrow("SimpleMatrix::assignment of a bad type of matrix into a triangular one.");
      break;
    }
    break;
  case 3:
    if (m.getNum() == 3)
      *mat.Sym = m.getSym();
    else
      SiconosMatrixException::selfThrow("SimpleMatrix::bad assignment of matrix (symetric one = dense or ...)");
    break;
  case 4:
    switch (m.getNum())
    {
    case 2:
      *mat.Sparse = m.getTriang();
      break;
    case 3:
      *mat.Sparse = m.getSym();
      break;
    case 4:
      *mat.Sparse = m.getSparse();
      break;
    case 5:
      *mat.Sparse = m.getBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("SimpleMatrix::op= (const SimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 5:
    switch (m.getNum())
    {
    case 5:
      *mat.Banded = m.getBanded();
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
  return *this;
}

SimpleMatrix& SimpleMatrix::operator = (const SiconosMatrix& m)
{
  switch (num)
  {
  case 1:
    switch (m.getNum())
    {
    case 1:
      *mat.Dense = m.getDense();
      break;
    case 2:
      *mat.Dense = m.getTriang();
      break;
    case 3:
      *mat.Dense = m.getSym();
      break;
    case 4:
      *mat.Dense = m.getSparse();
      break;
    case 5:
      *mat.Dense = m.getBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("SimpleMatrix::op= (const SiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 2:
    switch (m.getNum())
    {
    case 2:
      *mat.Triang = m.getTriang();
      break;
    default:
      SiconosMatrixException::selfThrow("SimpleMatrix::op= (const SiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 3:
    switch (m.getNum())
    {
    case 3:
      *mat.Sym = m.getSym();
      break;
    default:
      SiconosMatrixException::selfThrow("SimpleMatrix::op= (const SiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 4:
    switch (m.getNum())
    {
    case 2:
      *mat.Sparse = m.getTriang();
      break;
    case 3:
      *mat.Sparse = m.getSym();
      break;
    case 4:
      *mat.Sparse = m.getSparse();
      break;
    case 5:
      *mat.Sparse = m.getBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("SimpleMatrix::op= (const SimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 5:
    switch (m.getNum())
    {
    case 5:
      *mat.Banded = m.getBanded();
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
  return *this;
}

SimpleMatrix& SimpleMatrix::operator += (const SiconosMatrix& m)
{
  switch (num)
  {
  case 1:
    switch (m.getNum())
    {
    case 1:
      *mat.Dense += m.getDense();
      break;
    case 2:
      *mat.Dense += m.getTriang();
      break;
    case 3:
      *mat.Dense += m.getSym();
      break;
    case 4:
      *mat.Dense += m.getSparse();
      break;
    case 5:
      *mat.Dense += m.getBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op+= (const SiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 2:
    switch (m.getNum())
    {
    case 2:
      *mat.Triang += m.getTriang();
      break;
    default:
      SiconosMatrixException::selfThrow("op+= (const SiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 3:
    switch (m.getNum())
    {
    case 3:
      *mat.Sym += m.getSym();
      break;
    default:
      SiconosMatrixException::selfThrow("op+= (const SiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 4:
    switch (m.getNum())
    {
    case 2:
      *mat.Sparse += m.getTriang();
      break;
    case 3:
      *mat.Sparse += m.getSym();
      break;
    case 4:
      *mat.Sparse += m.getSparse();
      break;
    case 5:
      *mat.Sparse += m.getBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op= (const SimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 5:
    switch (m.getNum())
    {
    case 5:
      *mat.Banded += m.getBanded();
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
  return *this;
}

SimpleMatrix& SimpleMatrix::operator -= (const SiconosMatrix& m)
{
  switch (num)
  {
  case 1:
    switch (m.getNum())
    {
    case 1:
      *mat.Dense -= m.getDense();
      break;
    case 2:
      *mat.Dense -= m.getTriang();
      break;
    case 3:
      *mat.Dense -= m.getSym();
      break;
    case 4:
      *mat.Dense -= m.getSparse();
      break;
    case 5:
      *mat.Dense -= m.getBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op-= (const SiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 2:
    switch (m.getNum())
    {
    case 2:
      *mat.Triang -= m.getTriang();
      break;
    default:
      SiconosMatrixException::selfThrow("op-= (const SiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 3:
    switch (m.getNum())
    {
    case 3:
      *mat.Sym -= m.getSym();
      break;
    default:
      SiconosMatrixException::selfThrow("op-= (const SiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 4:
    switch (m.getNum())
    {
    case 2:
      *mat.Sparse -= m.getTriang();
      break;
    case 3:
      *mat.Sparse -= m.getSym();
      break;
    case 4:
      *mat.Sparse -= m.getSparse();
      break;
    case 5:
      *mat.Sparse -= m.getBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op-= (const SimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 5:
    switch (m.getNum())
    {
    case 5:
      *mat.Banded -= m.getBanded();
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
  default:
    SiconosMatrixException::selfThrow("op*= (int) : invalid type of matrix");
    break;
  }
  return *this;
}

SimpleMatrix& SimpleMatrix::operator /= (double m)
{
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
  default:
    SiconosMatrixException::selfThrow("op/= (double) : invalid type of matrix");
    break;
  }
  return *this;
}

SimpleMatrix& SimpleMatrix::operator /= (int m)
{
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
  default:
    SiconosMatrixException::selfThrow("op/= (int) : invalid type of matrix");
    break;
  }
  return *this;
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

  if (x.getNum() != m.getNum())
    SiconosMatrixException::selfThrow("SimpleMatrix:Matrix addition: use function add in order to add matrices of different type");

  if (x.getNum() == 1)
  {
    DenseMat p = x.getDense() + m.getDense();
    return p;
  }
  else if (x.getNum() == 2)
  {
    TriangMat t = x.getTriang() + m.getTriang();
    return t;
  }
  else if (x.getNum() == 3)
  {
    SymMat s = x.getSym() + m.getSym();
    return s;
  }
  else if (x.getNum() == 4)
  {
    SparseMat sp = x.getSparse() + m.getSparse();
    return sp;
  }
  else // if(x.getNum() == 5){
  {
    BandedMat b;
    b.resize(m.size(0), m.size(1), (m.getBanded()).lower(), (m.getBanded()).upper(), false);
    b = x.getBanded() + m.getBanded();
    return b;
  }
}

SimpleMatrix operator - (const SiconosMatrix &x, const SiconosMatrix &m)
{
  if ((x.size(0) != m.size(0)) || (x.size(1) != m.size(1)))
    SiconosMatrixException::selfThrow("Matrix subtraction: inconsistent sizes");

  if (x.getNum() != m.getNum())
    SiconosMatrixException::selfThrow("SimpleMatrix:Matrix subtraction: use function sub in order to subtract matrices of different type");

  if (x.getNum() == 1)
  {
    DenseMat p = x.getDense() - m.getDense();
    return p;
  }
  else if (x.getNum() == 2)
  {
    TriangMat t = x.getTriang() - m.getTriang();
    return t;
  }
  else if (x.getNum() == 3)
  {
    SymMat s = x.getSym() - m.getSym();
    return s;
  }
  else if (x.getNum() == 4)
  {
    SparseMat sp = x.getSparse() - m.getSparse();
    return sp;
  }
  else // if(x.getNum() == 5){
  {
    BandedMat b;
    b.resize(m.size(0), m.size(1), (m.getBanded()).lower(), (m.getBanded()).upper(), false);
    b = x.getBanded() - m.getBanded();
    return b;
  }
}

SimpleMatrix operator * (const SiconosMatrix &x, const SiconosMatrix &m)
{

  if ((x.size(1) != m.size(0)))
    SiconosMatrixException::selfThrow("Matrix product : inconsistent sizes");

  if (x.getNum() != m.getNum())
    SiconosMatrixException::selfThrow("SimpleMatrix:Matrix product : use function prod in order to multiply matrices of different type");

  if (x.getNum() == 1)
  {
    DenseMat p = prod(x.getDense(), m.getDense());
    return p;
  }
  else if (x.getNum() == 2)
  {
    TriangMat t = prod(x.getTriang(), m.getTriang());
    return t;
  }
  else if (x.getNum() == 3)
  {
    SymMat s = prod(x.getSym(), m.getSym());
    return s;
  }
  else if (x.getNum() == 4)
  {
    SparseMat sp = prod(x.getSparse(), m.getSparse());
    return sp;
  }
  else // if(x.getNum() == 5){
  {
    DenseMat p = prod(x.getBanded(), m.getBanded());
    return p;
  }
}

SimpleMatrix add(const SiconosMatrix &x, const SiconosMatrix& m)
{
  DenseMat p;

  if ((x.size(0) != m.size(0)) || (x.size(1) != m.size(1)))
    SiconosMatrixException::selfThrow("Matrix function add: inconsistent sizes");

  if (m.getNum() == 1)
  {
    DenseMat q;
    q = m.getDense();
    if (x.getNum() == 1)
    {
      p = x.getDense() + q;
    }
    else if (x.getNum() == 2)
    {
      p = x.getTriang() + q;
    }
    else if (x.getNum() == 3)
    {
      p = x.getSym() + q;
    }
    else if (x.getNum() == 4)
    {
      p = x.getSparse() + q;
    }
    else if (x.getNum() == 5)
    {
      p = x.getBanded() + q;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function add: invalid type of matrix");

  }
  else if (m.getNum() == 2)
  {
    TriangMat q;
    q = m.getTriang();
    if (x.getNum() == 1)
    {
      p = x.getDense() + q;
    }
    else if (x.getNum() == 2)
    {
      p = x.getTriang() + q;
    }
    else if (x.getNum() == 3)
    {
      p = x.getSym() + q;
    }
    else if (x.getNum() == 4)
    {
      p = x.getSparse() + q;
    }
    else if (x.getNum() == 5)
    {
      p = x.getBanded() + q;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function add: invalid type of matrix");

  }
  else if (m.getNum() == 3)
  {
    SymMat q;
    q = m.getSym();
    if (x.getNum() == 1)
    {
      p = x.getDense() + q;
    }
    else if (x.getNum() == 2)
    {
      p = x.getTriang() + q;
    }
    else if (x.getNum() == 3)
    {
      p = x.getSym() + q;
    }
    else if (x.getNum() == 4)
    {
      p = x.getSparse() + q;
    }
    else if (x.getNum() == 5)
    {
      p = x.getBanded() + q;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function add: invalid type of matrix");

  }
  else if (m.getNum() == 4)
  {
    SparseMat q;
    q = m.getSparse();
    if (x.getNum() == 1)
    {
      p = x.getDense() + q;
    }
    else if (x.getNum() == 2)
    {
      p = x.getTriang() + q;
    }
    else if (x.getNum() == 3)
    {
      p = x.getSym() + q;
    }
    else if (x.getNum() == 4)
    {
      p = x.getSparse() + q;
    }
    else if (x.getNum() == 5)
    {
      p = x.getBanded() + q;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function add: invalid type of matrix");

  }
  else if (m.getNum() == 5)
  {
    BandedMat q;
    q = m.getBanded();
    if (x.getNum() == 1)
    {
      p = x.getDense() + q;
    }
    else if (x.getNum() == 2)
    {
      p = x.getTriang() + q;
    }
    else if (x.getNum() == 3)
    {
      p = x.getSym() + q;
    }
    else if (x.getNum() == 4)
    {
      p = x.getSparse() + q;
    }
    else if (x.getNum() == 5)
    {
      p = x.getBanded() + q;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function add: invalid type of matrix");

  }
  else
  {
    SiconosMatrixException::selfThrow("Matrix function add: invalid type of matrix");
  }
  return p;

}

SimpleMatrix sub(const SiconosMatrix &x, const SiconosMatrix& m)
{
  DenseMat p;

  if ((x.size(0) != m.size(0)) || (x.size(1) != m.size(1)))
    SiconosMatrixException::selfThrow("Matrix function sub: inconsistent sizes");

  if (m.getNum() == 1)
  {
    DenseMat q;
    q = m.getDense();
    if (x.getNum() == 1)
    {
      p = x.getDense() - q;
    }
    else if (x.getNum() == 2)
    {
      p = x.getTriang() - q;
    }
    else if (x.getNum() == 3)
    {
      p = x.getSym() - q;
    }
    else if (x.getNum() == 4)
    {
      p = x.getSparse() - q;
    }
    else if (x.getNum() == 5)
    {
      p = x.getBanded() - q;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function sub: invalid type of matrix");

  }
  else if (m.getNum() == 2)
  {
    TriangMat q;
    q = m.getTriang();
    if (x.getNum() == 1)
    {
      p = x.getDense() - q;
    }
    else if (x.getNum() == 2)
    {
      p = x.getTriang() - q;
    }
    else if (x.getNum() == 3)
    {
      p = x.getSym() - q;
    }
    else if (x.getNum() == 4)
    {
      p = x.getSparse() - q;
    }
    else if (x.getNum() == 5)
    {
      p = x.getBanded() - q;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function sub: invalid type of matrix");

  }
  else if (m.getNum() == 3)
  {
    SymMat q;
    q = m.getSym();
    if (x.getNum() == 1)
    {
      p = x.getDense() - q;
    }
    else if (x.getNum() == 2)
    {
      p = x.getTriang() - q;
    }
    else if (x.getNum() == 3)
    {
      p = x.getSym() - q;
    }
    else if (x.getNum() == 4)
    {
      p = x.getSparse() - q;
    }
    else if (x.getNum() == 5)
    {
      p = x.getBanded() - q;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function sub: invalid type of matrix");

  }
  else if (m.getNum() == 4)
  {
    SparseMat q;
    q = m.getSparse();
    if (x.getNum() == 1)
    {
      p = x.getDense() - q;
    }
    else if (x.getNum() == 2)
    {
      p = x.getTriang() - q;
    }
    else if (x.getNum() == 3)
    {
      p = x.getSym() - q;
    }
    else if (x.getNum() == 4)
    {
      p = x.getSparse() - q;
    }
    else if (x.getNum() == 5)
    {
      p = x.getBanded() - q;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function add: invalid type of matrix");

  }
  else if (m.getNum() == 5)
  {
    BandedMat q;
    q = m.getBanded();
    if (x.getNum() == 1)
    {
      p = x.getDense() - q;
    }
    else if (x.getNum() == 2)
    {
      p = x.getTriang() - q;
    }
    else if (x.getNum() == 3)
    {
      p = x.getSym() - q;
    }
    else if (x.getNum() == 4)
    {
      p = x.getSparse() - q;
    }
    else if (x.getNum() == 5)
    {
      p = x.getBanded() - q;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function add: invalid type of matrix");

  }
  else
  {
    SiconosMatrixException::selfThrow("Matrix function sub: invalid type of matrix");
  }
  return p;

}

SimpleMatrix prod(const SiconosMatrix &x, const SiconosMatrix& m)
{
  DenseMat p;

  if ((x.size(1) != m.size(0)))
    SiconosMatrixException::selfThrow("Matrix function prod : inconsistent sizes");

  if (x.isBlock() || m.isBlock())
    SiconosMatrixException::selfThrow("Matrix function prod : not yet implemented for block matrices.");

  unsigned int num0 = m.getNum();
  if (num0 == 1)
  {
    if (x.getNum() == 1)
      p = prod(x.getDense(), m.getDense());

    else if (x.getNum() == 2)
      p = prod(x.getTriang(), m.getDense());

    else if (x.getNum() == 3)
      p = prod(x.getSym(), m.getDense());

    else if (x.getNum() == 4)
      p = prod(x.getSparse(), m.getDense());

    else //if(x.getNum()==5)
      p = prod(x.getBanded(), m.getDense());
    return p;
  }
  else if (num0 == 2)
  {
    if (x.getNum() == 1)
    {
      p = prod(x.getDense(), m.getTriang());
      return p;
    }
    else if (x.getNum() == 2)
    {
      TriangMat t = prod(x.getTriang(), m.getTriang());
      return t;
    }
    else if (x.getNum() == 3)
    {
      p = prod(x.getSym(), m.getTriang());
      return p;
    }
    else if (x.getNum() == 4)
    {
      p = prod(x.getSparse(), m.getTriang());
      return p;
    }
    else //if(x.getNum()==5)
    {
      p = prod(x.getBanded(), m.getTriang());
      return p;
    }
  }
  else if (num0 == 3)
  {
    if (x.getNum() == 1)
    {
      p = prod(x.getDense(), m.getSym());
      return p;
    }
    else if (x.getNum() == 2)
    {
      p = prod(x.getTriang(), m.getSym());
      return p;
    }
    else if (x.getNum() == 3)
    {
      SymMat s = prod(x.getSym(), m.getSym());
      return s;
    }
    else if (x.getNum() == 4)
    {
      p = prod(x.getSparse(), m.getSym());
      return p;
    }
    else  //if(x.getNum()==5){
    {
      p = prod(x.getBanded(), m.getSym());
      return p;
    }
  }
  else if (num0 == 4)
  {
    if (x.getNum() == 1)
    {
      p = prod(x.getDense(), m.getSparse());
      return p;
    }
    else if (x.getNum() == 2)
    {
      p = prod(x.getTriang(), m.getSparse());
      return p;
    }
    else if (x.getNum() == 3)
    {
      p = prod(x.getSym(), m.getSparse());
      return p;
    }
    else if (x.getNum() == 4)
    {
      SparseMat sp = prod(x.getSparse(), m.getSparse());
      return sp;
    }
    else  //if(x.getNum()==5){
    {
      p = prod(x.getBanded(), m.getSparse());
      return p;
    }
  }
  else // if(num0==5)
  {
    if (x.getNum() == 1)
    {
      p = prod(x.getDense(), m.getBanded());
      return p;
    }
    else if (x.getNum() == 2)
    {
      p = prod(x.getTriang(), m.getBanded());
      return p;
    }
    else if (x.getNum() == 3)
    {
      p = prod(x.getSym(), m.getBanded());
      return p;
    }
    else if (x.getNum() == 4)
    {
      p = prod(x.getSparse(), m.getBanded());
      return p;
    }
    else  //if(x.getNum()==5){
    {
      p = prod(x.getBanded(), m.getBanded());
      return p;
    }
  }
}

SimpleMatrix multTranspose(const SiconosMatrix &x, const SiconosMatrix &m)
{
  return prod(x, trans(m));
}

SimpleMatrix operator * (const SiconosMatrix &m, double d)
{
  DenseMat p;
  TriangMat t;
  SymMat s;
  SparseMat sp;
  BandedMat b;
  if (m.getNum() == 1)
  {
    p = m.getDense() * d;
    return p;
  }
  else if (m.getNum() == 2)
  {
    t = m.getTriang() * d;
    return t;
  }
  else if (m.getNum() == 3)
  {
    s = m.getSym() * d;
    return s;
  }
  else if (m.getNum() == 4)
  {
    sp = m.getSparse() * d;
    return sp;
  }
  else // if(m.getNum()==5){
  {
    b.resize(m.size(0), m.size(1), (m.getBanded()).lower(), (m.getBanded()).upper(), false);
    b = m.getBanded() * d;
    return b;
  }
}

SimpleMatrix operator * (const SiconosMatrix &m, int d)
{

  if (m.getNum() == 1)
  {
    DenseMat p = m.getDense() * d;
    return p;
  }
  else if (m.getNum() == 2)
  {
    TriangMat t = m.getTriang() * d;
    return t;
  }
  else if (m.getNum() == 3)
  {
    SymMat s = m.getSym() * d;
    return s;
  }
  else if (m.getNum() == 4)
  {
    SparseMat sp = m.getSparse() * d;
    return sp;
  }
  else // if(m.getNum()==5){
  {
    BandedMat b;
    b.resize(m.size(0), m.size(1), (m.getBanded()).lower(), (m.getBanded()).upper(), false);
    b = m.getBanded() * d;
    return b;
  }
}

SimpleMatrix operator * (double d, const SiconosMatrix &m)
{
  if (m.getNum() == 1)
  {
    DenseMat p = d * m.getDense();
    return p;
  }
  else if (m.getNum() == 2)
  {
    TriangMat t = d * m.getTriang();
    return t;
  }
  else if (m.getNum() == 3)
  {
    SymMat s = d * m.getSym();
    return s;
  }
  else if (m.getNum() == 4)
  {
    SparseMat sp = d * m.getSparse();
    return sp;
  }
  else  //if(m.getNum()==5){
  {
    BandedMat b;
    b.resize(m.size(0), m.size(1), (m.getBanded()).lower(), (m.getBanded()).upper(), false);
    b = d * m.getBanded();
    return b;
  }
}

SimpleMatrix operator * (int d, const SiconosMatrix &m)
{
  if (m.getNum() == 1)
  {
    DenseMat p = d * m.getDense();
    return p;
  }
  else if (m.getNum() == 2)
  {
    TriangMat t = d * m.getTriang();
    return t;
  }
  else if (m.getNum() == 3)
  {
    SymMat s = d * m.getSym();
    return s;
  }
  else if (m.getNum() == 4)
  {
    SparseMat sp = d * m.getSparse();
    return sp;
  }
  else  //if(m.getNum()==5){
  {
    BandedMat b;
    b.resize(m.size(0), m.size(1), (m.getBanded()).lower(), (m.getBanded()).upper(), false);
    b = d * m.getBanded();
    return b;
  }
}

SimpleMatrix operator / (const SiconosMatrix &m, double d)
{
  if (m.getNum() == 1)
  {
    DenseMat p = m.getDense() / d;
    return p;
  }
  else if (m.getNum() == 2)
  {
    TriangMat t = m.getTriang() / d;
    return t;
  }
  else if (m.getNum() == 3)
  {
    SymMat s = m.getSym() / d;
    return s;
  }
  else if (m.getNum() == 4)
  {
    SparseMat sp = m.getSparse() / d;
    return sp;
  }
  else// if(m.getNum()==5){
  {
    BandedMat b;
    b.resize(m.size(0), m.size(1), (m.getBanded()).lower(), (m.getBanded()).upper(), false);
    b = m.getBanded() / d;
    return b;
  }
}

SimpleMatrix operator / (const SiconosMatrix &m, int d)
{
  if (m.getNum() == 1)
  {
    DenseMat p = m.getDense() / d;
    return p;
  }
  else if (m.getNum() == 2)
  {
    TriangMat t = m.getTriang() / d;
    return t;
  }
  else if (m.getNum() == 3)
  {
    SymMat s = m.getSym() / d;
    return s;
  }
  else if (m.getNum() == 4)
  {
    SparseMat sp = m.getSparse() / d;
    return sp;
  }
  else // if(m.getNum()==5){
  {
    BandedMat b;
    b.resize(m.size(0), m.size(1), (m.getBanded()).lower(), (m.getBanded()).upper(), false);
    b = m.getBanded() / d;
    return b;
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
    if (m.getNum() == 1)
    {
      DenseMat p ;
      p = m.getDense();
      for (unsigned int i = 1; i < power; i++)
        p = prod(p, m.getDense());
      return p;
    }
    else if (m.getNum() == 2)
    {
      TriangMat t = m.getTriang();
      for (unsigned int i = 1; i < power; i++)
        t = prod(t, m.getTriang());
      return t;
    }
    else if (m.getNum() == 3)
    {
      SymMat s = m.getSym();
      for (unsigned int i = 1; i < power; i++)
        s = prod(s, m.getSym());
      return s;
    }
    else if (m.getNum() == 4)
    {
      SparseMat sp = m.getSparse();
      for (unsigned int i = 1; i < power; i++)
        sp = prod(sp, m.getSparse());
      return sp;
    }
    else// if(m.getNum() == 5)
    {
      DenseMat b = m.getBanded();
      for (unsigned int i = 1; i < power; i++)
        b = prod(b, m.getBanded());
      return b;
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
      res = prod(m.getDense(), tmp->getDense());
    else if (numM == 2)
      res = prod(m.getTriang(), tmp->getDense());
    else if (numM == 3)
      res = prod(m.getSym(), tmp->getDense());
    else if (numM == 4)
      res = prod(m.getSparse(), tmp->getDense());
    else if (numM == 5)
      res = prod(m.getBanded(), tmp->getDense());
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
        res = prod(m.getDense(), v.getDense());
      else if (numM == 2)
        res = prod(m.getTriang(), v.getDense());
      else if (numM == 3)
        res = prod(m.getSym(), v.getDense());
      else if (numM == 4)
        res = prod(m.getSparse(), v.getDense());
      else if (numM == 5)
        res = prod(m.getBanded(), v.getDense());
      else
        SiconosMatrixException::selfThrow("prod(matrix,vector) error: unknown matrix type.");
    }
    else if (numV == 4)
    {
      if (numM == 1)
        res = prod(m.getDense(), v.getSparse());
      else if (numM == 2)
        res = prod(m.getTriang(), v.getSparse());
      else if (numM == 3)
        res = prod(m.getSym(), v.getSparse());
      else if (numM == 4)
        res = prod(m.getSparse(), v.getSparse());
      else if (numM == 5)
        res = prod(m.getBanded(), v.getSparse());
      else
        SiconosMatrixException::selfThrow("prod(matrix,vector) error: unknown matrix type.");
    }
    else
      SiconosMatrixException::selfThrow("prod(matrix,vector) error: unknown vector type.");
  }
  return res;
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
    info = boost::numeric::bindings::atlas::gesv(*mat.Dense, *(B.getDensePtr()));
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
    info = boost::numeric::bindings::atlas::gesv(*mat.Dense, tmpB);
    // B now contains solution:
  }
  else // call getrs: only solve using previous lu-factorization
    info = boost::numeric::bindings::atlas::getrs(*mat.Dense, ipiv, tmpB);

  if (info != 0)
    SiconosMatrixException::selfThrow("SimpleMatrix::PLUForwardBackwardInPlace failed.");
  *(B.getDensePtr()) = ublas::column(tmpB, 0);
}
