#include "MySimpleMatrix.h"
#include <cassert>
#include "ioMatrix.h"

//Default private
MySimpleMatrix::MySimpleMatrix(void): num(0)
{
  SetIsBlock(false);
}
/***************************** CONSTRUCTORS ****************************/


// Build a Simple Matrix from its type (ie DENSE, TRIANGULAR, BANDED, SPARSE or SYMMETRIC)
MySimpleMatrix::MySimpleMatrix(TYP typ)
{
  SetIsBlock(false);
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

}

MySimpleMatrix::MySimpleMatrix(TYP typ, int row, int col)
{
  SetIsBlock(false);
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
    mat.Sparse = new SparseMat();
    num = 4;
  }
  else if (typ == BANDED)
  {
    mat.Banded = new BandedMat();
    num = 5;
  }
  else
    SiconosMatrixException::selfThrow("constructor(TYP type, int row, int col) : invalid type or dimensions given");
}


MySimpleMatrix::MySimpleMatrix(const MySimpleMatrix &smat)
{
  SetIsBlock(false);
  if (smat.GetNum() == 1)
  {
    mat.Dense = new DenseMat(smat.GetDense());
    num = 1;
  }
  else if (smat.GetNum() == 2)
  {
    mat.Triang = new TriangMat(smat.GetTriang());
    num = 2;
  }
  else if (smat.GetNum() == 3)
  {
    mat.Sym = new SymMat(smat.GetSym());
    num = 3;
  }
  else if (smat.GetNum() == 4)
  {
    mat.Sparse = new SparseMat(smat.GetSparse());
    num = 4;
  }
  else if (smat.GetNum() == 5)
  {
    mat.Banded = new BandedMat(smat.GetBanded());
    num = 5;
  }
  else
    SiconosMatrixException::selfThrow("constructor(const MySimpleMatrix) : invalid parameter given");
}

MySimpleMatrix::MySimpleMatrix(const MySiconosMatrix &smat)
{
  SetIsBlock(false);
  assert(smat.isBlock() == false);
  if (smat.GetNum() == 1)
  {
    mat.Dense = new DenseMat(smat.GetDense());
    num = 1;
  }
  else if (smat.GetNum() == 2)
  {
    mat.Triang = new TriangMat(smat.GetTriang());
    num = 2;
  }
  else if (smat.GetNum() == 3)
  {
    mat.Sym = new SymMat(smat.GetSym());
    num = 3;
  }
  else if (smat.GetNum() == 4)
  {
    mat.Sparse = new SparseMat(smat.GetSparse());
    num = 4;
  }
  else if (smat.GetNum() == 5)
  {
    mat.Banded = new BandedMat(smat.GetBanded());
    num = 5;
  }
  else
    SiconosMatrixException::selfThrow("constructor(const MySiconosMatrix) : invalid parameter given");
}

MySimpleMatrix::MySimpleMatrix(const DenseMat& m)
{
  SetIsBlock(false);
  mat.Dense = new DenseMat(m);
  num = 1;
}

MySimpleMatrix::MySimpleMatrix(const TriangMat& m)
{
  SetIsBlock(false);
  mat.Triang = new TriangMat(m);
  num = 2;
}

MySimpleMatrix::MySimpleMatrix(const SymMat& m)
{
  SetIsBlock(false);
  mat.Sym = new SymMat(m);
  num = 3;
}

MySimpleMatrix::MySimpleMatrix(const SparseMat& m)
{
  SetIsBlock(false);
  mat.Sparse = new SparseMat(m);
  num = 4;
}

MySimpleMatrix::MySimpleMatrix(const BandedMat& m)
{
  SetIsBlock(false);
  mat.Banded = new BandedMat(m);
  num = 5;
}

MySimpleMatrix::MySimpleMatrix(TYP typ, const std::vector<double> &v, int row, int col, int lower, int upper)
{

  if (((v.size()) != (unsigned int)row * col && (typ != SYMMETRIC && typ != BANDED)) || (v.size() != (unsigned int)row * row && typ == SYMMETRIC) || (typ == BANDED && ((v.size()) != (unsigned int)(std::max)(row, col) * (lower + 1 + upper))))
    SiconosMatrixException::selfThrow("constructor(TYP, const std::vector<double>, int, int) : invalid vector size");

  SetIsBlock(false);
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

MySimpleMatrix::MySimpleMatrix(const std::string &file, bool ascii)
{
  SetIsBlock(false);
  num = 1;
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
MySimpleMatrix::~MySimpleMatrix(void)
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
int  MySimpleMatrix::GetNum(void)const
{
  return num;
}
void MySimpleMatrix::SetNum(int n)
{
  num = n;
}
int  MySimpleMatrix::size1(void)const
{
  int size1;
  if (num == 1)
  {
    size1 = (*mat.Dense).size1();
  }
  else if (num == 2)
  {
    size1 = (*mat.Triang).size1();
  }
  else if (num == 3)
  {
    size1 = (*mat.Sym).size1();
  }
  else if (num == 4)
  {
    size1 = (*mat.Sparse).size1();
  }
  else if (num == 5)
  {
    size1 = (*mat.Banded).size1();
  }
  return size1;
}
int  MySimpleMatrix::size2(void)const
{
  int size2;
  if (num == 1)
  {
    size2 = (*mat.Dense).size2();
  }
  else if (num == 2)
  {
    size2 = (*mat.Triang).size2();
  }
  else if (num == 3)
  {
    size2 = (*mat.Sym).size2();
  }
  else if (num == 4)
  {
    size2 = (*mat.Sparse).size2();
  }
  else if (num == 5)
  {
    size2 = (*mat.Banded).size2();
  }
  return size2;
}

void MySimpleMatrix::resize(int row, int col, int lower, int upper, bool preserve)
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

const DenseMat MySimpleMatrix::GetDense(int row, int col)const
{

  if (num != 1)
    SiconosMatrixException::selfThrow("DenseMat GetDense(int row, int col) : the current matrix is not a Dense matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("DenseMat GetDense(int row, int col) : row or col not equal to 0.0");

  return *mat.Dense;
}

const TriangMat MySimpleMatrix::GetTriang(int row, int col)const
{

  if (num != 2)
    SiconosMatrixException::selfThrow("TriangMat GetTriang(int row, int col) : the current matrix is not a Triangular matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("TriangMat GetTriang(int row, int col) : row or col not equal to 0.0");

  return *mat.Triang;
}

const SymMat MySimpleMatrix::GetSym(int row, int col)const
{

  if (num != 3)
    SiconosMatrixException::selfThrow("SymMat GetSym(int row, int col) : the current matrix is not a Symmetric matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("SymMat GetSym(int row, int col) : row or col not equal to 0.0");

  return *mat.Sym;
}

const SparseMat MySimpleMatrix::GetSparse(int row, int col)const
{

  if (num != 4)
    SiconosMatrixException::selfThrow("SparseMat GetSparse(int row, int col) : the current matrix is not a Sparse matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("SparseMat GetSparse(int row, int col) : row or col not equal to 0.0");

  return *mat.Sparse;
}

const BandedMat MySimpleMatrix::GetBanded(int row, int col)const
{

  if (num != 5)
    SiconosMatrixException::selfThrow("BandedMat GetBanded(int row, int col) : the current matrix is not a Banded matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("BandedMat GetBanded(int row, int col) : row or col not equal to 0.0");

  return *mat.Banded;
}

const DenseMat* MySimpleMatrix::GetDensePtr(int row, int col)const
{

  if (num != 1)
    SiconosMatrixException::selfThrow("DenseMat* GetDensePtr(int row, int col) : the current matrix is not a Dense matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("DenseMat* GetDensePtr(int row, int col) : row or col not equal to 0.0");

  return mat.Dense;
}

const TriangMat* MySimpleMatrix::GetTriangPtr(int row, int col)const
{

  if (num != 2)
    SiconosMatrixException::selfThrow("TriangMat* GetTriangPtr(int row, int col) : the current matrix is not a Triangular matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("TriangMat* GetTriangPtr(int row, int col) : row or col not equal to 0.0");

  return mat.Triang;
}

const SymMat* MySimpleMatrix::GetSymPtr(int row, int col)const
{

  if (num != 3)
    SiconosMatrixException::selfThrow("SymMat* GetSymPtr(int row, int col) : the current matrix is not a Symmetric matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("SymMat* GetSymPtr(int row, int col) : row or col not equal to 0.0");

  return mat.Sym;
}

const SparseMat* MySimpleMatrix::GetSparsePtr(int row, int col)const
{

  if (num != 4)
    SiconosMatrixException::selfThrow("SparseMat* GetSparsePtr(int row, int col) : the current matrix is not a Sparse matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("SparseMat* GetSparsePtr(int row, int col) : row or col not equal to 0.0");

  return mat.Sparse;
}

const BandedMat* MySimpleMatrix::GetBandedPtr(int row, int col)const
{

  if (num != 5)
    SiconosMatrixException::selfThrow("BandedMat* GetBandedPtr(int row, int col) : the current matrix is not a Banded matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("BandedMat* GetBandedPtr(int row, int col) : row or col not equal to 0.0");

  return mat.Banded;
}

const mapped MySimpleMatrix::GetMap(void)const
{
  SiconosMatrixException::selfThrow("mapped GetMap : GetMap is forbidden for MySimpleMatrix");
}

void MySimpleMatrix::BlockMatrixCopy(const MySiconosMatrix &m, int i, int j)
{

  if (i >= size1() || i < 0)
    SiconosMatrixException::selfThrow("BlockMatriCopy : row_min given is out of range");

  if (j >= size2() || j < 0)
    SiconosMatrixException::selfThrow("BlockMatriCopy : col_min given is out of range");
  if (num == 1)
  {
    if (num == 1)
    {
      int row_max = i + (m.GetDense()).size1();
      int col_max = j + (m.GetDense()).size2();
      if (row_max >= size1() || row_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      if (col_max >= size2() || col_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      subrange(*mat.Dense, i, i + (m.GetDense()).size1(), j, j + (m.GetDense()).size2()) = m.GetDense();
    }
    else if (num == 2)
    {
      int row_max = i + (m.GetTriang()).size1();
      int col_max = j + (m.GetTriang()).size2();
      if (row_max >= size1() || row_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      if (col_max >= size2() || col_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      subrange(*mat.Dense, i, i + (m.GetTriang()).size1(), j, j + (m.GetTriang()).size2()) = m.GetTriang();
    }
    else if (num == 3)
    {
      int row_max = i + (m.GetSym()).size1();
      int col_max = j + (m.GetSym()).size2();
      if (row_max >= size1() || row_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      if (col_max >= size2() || col_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      subrange(*mat.Dense, i, i + (m.GetSym()).size1(), j, j + (m.GetSym()).size2()) = m.GetSym();
    }
    else if (num == 4)
    {
      int row_max = i + (m.GetSparse()).size1();
      int col_max = j + (m.GetSparse()).size2();
      if (row_max >= size1() || row_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      if (col_max >= size2() || col_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      subrange(*mat.Dense, i, i + (m.GetSparse()).size1(), j, j + (m.GetSparse()).size2()) = m.GetSparse();
    }
    else if (num == 5)
    {
      int row_max = i + (m.GetBanded()).size1();
      int col_max = j + (m.GetBanded()).size2();
      if (row_max >= size1() || row_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      if (col_max >= size2() || col_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      subrange(*mat.Dense, i, i + (m.GetBanded()).size1(), j, j + (m.GetBanded()).size2()) = m.GetBanded();
    }
  }
  else if (num == 2)
  {
    if (num == 1)
    {
      int row_max = i + (m.GetDense()).size1();
      int col_max = j + (m.GetDense()).size2();
      if (row_max >= size1() || row_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      if (col_max >= size2() || col_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      subrange(*mat.Triang, i, i + (m.GetDense()).size1(), j, j + (m.GetDense()).size2()) = m.GetDense();
    }
    else if (num == 2)
    {
      int row_max = i + (m.GetTriang()).size1();
      int col_max = j + (m.GetTriang()).size2();
      if (row_max >= size1() || row_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      if (col_max >= size2() || col_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      subrange(*mat.Triang, i, i + (m.GetTriang()).size1(), j, j + (m.GetTriang()).size2()) = m.GetTriang();
    }
    else if (num == 3)
    {
      int row_max = i + (m.GetSym()).size1();
      int col_max = j + (m.GetSym()).size2();
      if (row_max >= size1() || row_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      if (col_max >= size2() || col_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      subrange(*mat.Triang, i, i + (m.GetSym()).size1(), j, j + (m.GetSym()).size2()) = m.GetSym();
    }
    else if (num == 4)
    {
      int row_max = i + (m.GetSparse()).size1();
      int col_max = j + (m.GetSparse()).size2();
      if (row_max >= size1() || row_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      if (col_max >= size2() || col_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      subrange(*mat.Triang, i, i + (m.GetSparse()).size1(), j, j + (m.GetSparse()).size2()) = m.GetSparse();
    }
    else if (num == 5)
    {
      int row_max = i + (m.GetBanded()).size1();
      int col_max = j + (m.GetBanded()).size2();
      if (row_max >= size1() || row_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      if (col_max >= size2() || col_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      subrange(*mat.Triang, i, i + (m.GetBanded()).size1(), j, j + (m.GetBanded()).size2()) = m.GetBanded();
    }
  }
  else if (num == 3)
  {
    if (num == 1)
    {
      int row_max = i + (m.GetDense()).size1();
      int col_max = j + (m.GetDense()).size2();
      if (row_max >= size1() || row_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      if (col_max >= size2() || col_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      subrange(*mat.Sym, i, i + (m.GetDense()).size1(), j, j + (m.GetDense()).size2()) = m.GetDense();
    }
    else if (num == 2)
    {
      int row_max = i + (m.GetTriang()).size1();
      int col_max = j + (m.GetTriang()).size2();
      if (row_max >= size1() || row_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      if (col_max >= size2() || col_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      subrange(*mat.Sym, i, i + (m.GetTriang()).size1(), j, j + (m.GetTriang()).size2()) = m.GetTriang();
    }
    else if (num == 3)
    {
      int row_max = i + (m.GetSym()).size1();
      int col_max = j + (m.GetSym()).size2();
      if (row_max >= size1() || row_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      if (col_max >= size2() || col_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      subrange(*mat.Sym, i, i + (m.GetSym()).size1(), j, j + (m.GetSym()).size2()) = m.GetSym();
    }
    else if (num == 4)
    {
      int row_max = i + (m.GetSparse()).size1();
      int col_max = j + (m.GetSparse()).size2();
      if (row_max >= size1() || row_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      if (col_max >= size2() || col_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      subrange(*mat.Sym, i, i + (m.GetSparse()).size1(), j, j + (m.GetSparse()).size2()) = m.GetSparse();
    }
    else if (num == 5)
    {
      int row_max = i + (m.GetBanded()).size1();
      int col_max = j + (m.GetBanded()).size2();
      if (row_max >= size1() || row_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      if (col_max >= size2() || col_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      subrange(*mat.Sym, i, i + (m.GetBanded()).size1(), j, j + (m.GetBanded()).size2()) = m.GetBanded();
    }
  }
  if (num == 4)
  {
    if (num == 1)
    {
      int row_max = i + (m.GetDense()).size1();
      int col_max = j + (m.GetDense()).size2();
      if (row_max >= size1() || row_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      if (col_max >= size2() || col_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      subrange(*mat.Sparse, i, i + (m.GetDense()).size1(), j, j + (m.GetDense()).size2()) = m.GetDense();
    }
    else if (num == 2)
    {
      int row_max = i + (m.GetTriang()).size1();
      int col_max = j + (m.GetTriang()).size2();
      if (row_max >= size1() || row_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      if (col_max >= size2() || col_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      subrange(*mat.Sparse, i, i + (m.GetTriang()).size1(), j, j + (m.GetTriang()).size2()) = m.GetTriang();
    }
    else if (num == 3)
    {
      int row_max = i + (m.GetSym()).size1();
      int col_max = j + (m.GetSym()).size2();
      if (row_max >= size1() || row_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      if (col_max >= size2() || col_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      subrange(*mat.Sparse, i, i + (m.GetSym()).size1(), j, j + (m.GetSym()).size2()) = m.GetSym();
    }
    else if (num == 4)
    {
      int row_max = i + (m.GetSparse()).size1();
      int col_max = j + (m.GetSparse()).size2();
      if (row_max >= size1() || row_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      if (col_max >= size2() || col_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      subrange(*mat.Sparse, i, i + (m.GetSparse()).size1(), j, j + (m.GetSparse()).size2()) = m.GetSparse();
    }
    else if (num == 5)
    {
      int row_max = i + (m.GetBanded()).size1();
      int col_max = j + (m.GetBanded()).size2();
      if (row_max >= size1() || row_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      if (col_max >= size2() || col_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      subrange(*mat.Sparse, i, i + (m.GetBanded()).size1(), j, j + (m.GetBanded()).size2()) = m.GetBanded();
    }
  }
  if (num == 5)
  {
    if (num == 1)
    {
      int row_max = i + (m.GetDense()).size1();
      int col_max = j + (m.GetDense()).size2();
      if (row_max >= size1() || row_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      if (col_max >= size2() || col_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      subrange(*mat.Banded, i, i + (m.GetDense()).size1(), j, j + (m.GetDense()).size2()) = m.GetDense();
    }
    else if (num == 2)
    {
      int row_max = i + (m.GetTriang()).size1();
      int col_max = j + (m.GetTriang()).size2();
      if (row_max >= size1() || row_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      if (col_max >= size2() || col_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      subrange(*mat.Banded, i, i + (m.GetTriang()).size1(), j, j + (m.GetTriang()).size2()) = m.GetTriang();
    }
    else if (num == 3)
    {
      int row_max = i + (m.GetSym()).size1();
      int col_max = j + (m.GetSym()).size2();
      if (row_max >= size1() || row_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      if (col_max >= size2() || col_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      subrange(*mat.Banded, i, i + (m.GetSym()).size1(), j, j + (m.GetSym()).size2()) = m.GetSym();
    }
    else if (num == 4)
    {
      int row_max = i + (m.GetSparse()).size1();
      int col_max = j + (m.GetSparse()).size2();
      if (row_max >= size1() || row_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      if (col_max >= size2() || col_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      subrange(*mat.Banded, i, i + (m.GetSparse()).size1(), j, j + (m.GetSparse()).size2()) = m.GetSparse();
    }
    else if (num == 5)
    {
      int row_max = i + (m.GetBanded()).size1();
      int col_max = j + (m.GetBanded()).size2();
      if (row_max >= size1() || row_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      if (col_max >= size2() || col_max < 0)
        SiconosMatrixException::selfThrow("BlockMatriCopy : inconsistent sizes");

      subrange(*mat.Banded, i, i + (m.GetBanded()).size1(), j, j + (m.GetBanded()).size2()) = m.GetBanded();
    }
  }
}
void MySimpleMatrix::GetBlock(int row_min, int col_min, MySiconosMatrix &m)const
{

  if (row_min >= size1() || row_min < 0)
    SiconosMatrixException::selfThrow("GetBlock : row_min given is out of range");

  if (col_min >= size2() || row_min < 0)
    SiconosMatrixException::selfThrow("GetBlock : col_min given is out of range");
  DenseMat q;
  int row_max, col_max;
  if (num == 1)
  {
    row_max = m.GetDense().size1() + row_min;
    col_max = m.GetDense().size2() + col_min;

    if (row_max >= size1() || row_max < 0)
      SiconosMatrixException::selfThrow("GetBlock : inconsistent sizes");

    if (col_max >= size2() || col_max < 0)
      SiconosMatrixException::selfThrow("GetBlock : inconsistent sizes");

    q = subrange(*mat.Dense, row_min, row_max, col_min, col_max);
  }
  else if (num == 2)
  {
    row_max = m.GetTriang().size1() + row_min;
    col_max = m.GetTriang().size2() + col_min;

    if (row_max >= size1() || row_max < 0)
      SiconosMatrixException::selfThrow("GetBlock : inconsistent sizes");

    if (col_max >= size2() || col_max < 0)
      SiconosMatrixException::selfThrow("GetBlock : inconsistent sizes");
    q = subrange(*mat.Triang, row_min, row_max, col_min, col_max);
  }
  else if (num == 3)
  {
    row_max = m.GetSym().size1() + row_min;
    col_max = m.GetSym().size2() + col_min;

    if (row_max >= size1() || row_max < 0)
      SiconosMatrixException::selfThrow("GetBlock : inconsistent sizes");

    if (col_max >= size2() || col_max < 0)
      SiconosMatrixException::selfThrow("GetBlock : inconsistent sizes");
    q = subrange(*mat.Sym, row_min, row_max, col_min, col_max);
  }
  else if (num == 4)
  {
    row_max = m.GetSparse().size1() + row_min;
    col_max = m.GetSparse().size2() + col_min;

    if (row_max >= size1() || row_max < 0)
      SiconosMatrixException::selfThrow("GetBlock : inconsistent sizes");

    if (col_max >= size2() || col_max < 0)
      SiconosMatrixException::selfThrow("GetBlock : inconsistent sizes");
    q = subrange(*mat.Sparse, row_min, row_max, col_min, col_max);
  }
  else if (num == 5)
  {
    row_max = m.GetBanded().size1() + row_min;
    col_max = m.GetBanded().size2() + col_min;

    if (row_max >= size1() || row_max < 0)
      SiconosMatrixException::selfThrow("GetBlock : inconsistent sizes");

    if (col_max >= size2() || col_max < 0)
      SiconosMatrixException::selfThrow("GetBlock : inconsistent sizes");
    q = subrange(*mat.Banded, row_min, row_max, col_min, col_max);
  }
  MySimpleMatrix p(q);
  m = p;
}

const std::deque<bool> MySimpleMatrix::GetBlockAllocated(void)const
{
  SiconosMatrixException::selfThrow("std::deque<bool> GetBlockAllocated : GetBlockAllocated is forbidden for MySimpleMatrix");
}

void MySimpleMatrix::GetRow(int r, MySimpleVector &vect)const
{

  if (r >= size1() || r < 0)
    SiconosMatrixException::selfThrow("GetRow : row is out of range");

  if (vect.size() != size2())
    SiconosMatrixException::selfThrow("GetRow : inconsistent sizes");

  DenseVect v1;
  if (num == 1)
  {
    v1 = row(*mat.Dense, r);
  }
  else if (num == 2)
  {
    v1 = row(*mat.Triang, r);
  }
  else if (num == 3)
  {
    v1 = row(*mat.Sym, r);
  }
  else if (num == 4)
  {
    v1 = row(*mat.Sparse, r);
  }
  else if (num == 5)
  {
    v1 = row(*mat.Banded, r);
  }
  MySimpleVector p(v1);
  vect = p;
}

void MySimpleMatrix::SetRow(int r, const MySimpleVector &vect)
{

  if (r >= size1() || r < 0)
    SiconosMatrixException::selfThrow("SetRow : row is out of range");

  if (vect.size() != size2())
    SiconosMatrixException::selfThrow("SetRow : inconsistent sizes");

  if (num == 1)
  {
    if (vect.GetNum() == 1)
    {
      row(*mat.Dense, r) = vect.GetDense();
    }
    else if (vect.GetNum() == 2)
    {
      row(*mat.Dense, r) = vect.GetSparse();
    }
  }
  else if (num == 2)
  {
    if (vect.GetNum() == 1)
    {
      row(*mat.Triang, r) = vect.GetDense();
    }
    else if (vect.GetNum() == 2)
    {
      row(*mat.Triang, r) = vect.GetSparse();
    }
  }
  else if (num == 3)
  {
    if (vect.GetNum() == 1)
    {
      row(*mat.Sym, r) = vect.GetDense();
    }
    else if (vect.GetNum() == 2)
    {
      row(*mat.Sym, r) = vect.GetSparse();
    }
  }
  if (num == 4)
  {
    if (vect.GetNum() == 1)
    {
      row(*mat.Sparse, r) = vect.GetDense();
    }
    else if (vect.GetNum() == 2)
    {
      row(*mat.Sparse, r) = vect.GetSparse();
    }
  }
  if (num == 5)
  {
    if (vect.GetNum() == 1)
    {
      row(*mat.Sparse, r) = vect.GetDense();
    }
    else if (vect.GetNum() == 2)
    {
      row(*mat.Banded, r) = vect.GetSparse();
    }
  }
}

void MySimpleMatrix::GetCol(int r, MySimpleVector &vect)const
{

  if (r >= size2() || r < 0)
    SiconosMatrixException::selfThrow("GetCol : col is out of range");

  if (vect.size() != size1())
    SiconosMatrixException::selfThrow("GetCol : inconsistent sizes");

  DenseVect v1;
  if (num == 1)
  {
    v1 = column(*mat.Dense, r);
  }
  else if (num == 2)
  {
    v1 = column(*mat.Triang, r);
  }
  else if (num == 3)
  {
    v1 = column(*mat.Sym, r);
  }
  else if (num == 4)
  {
    v1 = column(*mat.Sparse, r);
  }
  else if (num == 5)
  {
    v1 = column(*mat.Banded, r);
  }
  MySimpleVector p(v1);
  vect = p;
}

void MySimpleMatrix::SetCol(int r, const MySimpleVector &vect)
{

  if (r >= size2() || r < 0)
    SiconosMatrixException::selfThrow("SetCol : col is out of range");

  if (vect.size() != size1())
    SiconosMatrixException::selfThrow("SetCol : inconsistent sizes");

  if (num == 1)
  {
    if (vect.GetNum() == 1)
    {
      column(*mat.Dense, r) = vect.GetDense();
    }
    else if (vect.GetNum() == 2)
    {
      column(*mat.Dense, r) = vect.GetSparse();
    }
  }
  else if (num == 2)
  {
    if (vect.GetNum() == 1)
    {
      column(*mat.Triang, r) = vect.GetDense();
    }
    else if (vect.GetNum() == 2)
    {
      column(*mat.Triang, r) = vect.GetSparse();
    }
  }
  else if (num == 3)
  {
    if (vect.GetNum() == 1)
    {
      column(*mat.Sym, r) = vect.GetDense();
    }
    else if (vect.GetNum() == 2)
    {
      column(*mat.Sym, r) = vect.GetSparse();
    }
  }
  else if (num == 4)
  {
    if (vect.GetNum() == 1)
    {
      column(*mat.Sparse, r) = vect.GetDense();
    }
    else if (vect.GetNum() == 2)
    {
      column(*mat.Sparse, r) = vect.GetSparse();
    }
  }
  else if (num == 5)
  {
    if (vect.GetNum() == 1)
    {
      column(*mat.Banded, r) = vect.GetDense();
    }
    else if (vect.GetNum() == 2)
    {
      column(*mat.Banded, r) = vect.GetSparse();
    }
  }
}

const double MySimpleMatrix::normInf(void)const
{
  double d;
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

MySimpleMatrix trans(const MySiconosMatrix &m)
{
  if (m.GetNum() == 1)
  {
    DenseMat p = trans(m.GetDense());
    return p;
  }
  else if (m.GetNum() == 2)
  {
    TriangMat p = trans(m.GetTriang());
    return p;
  }
  else if (m.GetNum() == 3)
  {
    SymMat p = trans(m.GetSym());
    return p;
  }
  else if (m.GetNum() == 4)
  {
    SparseMat p = trans(m.GetSparse());
    return p;
  }
  else if (m.GetNum() == 5)
  {
    BandedMat p = trans(m.GetBanded());
    return p;
  }
}

void MySimpleMatrix::display(void)const
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

void MySimpleMatrix::zero(void)
{
  int size1 = (*this).size1();
  int size2 = (*this).size2();
  if (num == 1)
  {
    zero_matrix<double> p(size1, size2);
    *mat.Dense = p;
  }
  else if (num == 2)
  {
    zero_matrix<double> p(size1, size2);
    *mat.Triang = p;
  }
  else if (num == 3)
  {
    zero_matrix<double> p(size1, size2);
    *mat.Sym = p;
  }
  else if (num == 4)
  {
    zero_matrix<double> p(size1, size2);
    *mat.Sparse = p;
  }
  else if (num == 5)
  {
    zero_matrix<double> p(size1, size2);
    *mat.Banded = p;
  }
}

void MySimpleMatrix::eye(void)
{
  int size1 = (*this).size1();
  int size2 = (*this).size2();
  if (num == 1)
  {
    identity_matrix<double> p(size1, size2);
    *mat.Dense = p;
  }
  else if (num == 2)
  {
    identity_matrix<double> p(size1, size2);
    *mat.Triang = p;
  }
  else if (num == 3)
  {
    identity_matrix<double> p(size1, size2);
    *mat.Sym = p;
  }
  else if (num == 4)
  {
    identity_matrix<double> p(size1, size2);
    *mat.Sparse = p;
  }
  else if (num == 5)
  {
    identity_matrix<double> p(size1, size2);
    *mat.Banded = p;
  }
}

/***************************** OPERATORS ******************************/

double& MySimpleMatrix::operator()(int row, int col)
{
  double d;
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
    SiconosMatrixException::selfThrow("op() (int, int) : invalid type of matrix");
    break;
  }
  return d;
}

double MySimpleMatrix::operator()(int row, int col)const
{
  double d;
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
    SiconosMatrixException::selfThrow("op() (int, int) : invalid type of matrix");
    break;
  }
  return d;
}

const MySimpleMatrix& MySimpleMatrix::operator = (const MySimpleMatrix& m)
{
  switch (num)
  {
  case 1:
    switch (m.GetNum())
    {
    case 1:
      *mat.Dense = m.GetDense();
      break;
    case 2:
      *mat.Dense = m.GetTriang();
      break;
    case 3:
      *mat.Dense = m.GetSym();
      break;
    case 4:
      *mat.Dense = m.GetSparse();
      break;
    case 5:
      *mat.Dense = m.GetBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op= (const MySimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 2:
    switch (m.GetNum())
    {
    case 1:
      *mat.Triang = m.GetDense();
      break;
    case 2:
      *mat.Triang = m.GetTriang();
      break;
    case 3:
      *mat.Triang = m.GetSym();
      break;
    case 4:
      *mat.Triang = m.GetSparse();
      break;
    case 5:
      *mat.Triang = m.GetBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op= (const MySimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 3:
    switch (m.GetNum())
    {
    case 1:
      *mat.Sym = m.GetDense();
      break;
    case 2:
      *mat.Sym = m.GetTriang();
      break;
    case 3:
      *mat.Sym = m.GetSym();
      break;
    case 4:
      *mat.Sym = m.GetSparse();
      break;
    case 5:
      *mat.Sym = m.GetBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op= (const MySimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 4:
    switch (m.GetNum())
    {
    case 1:
      *mat.Sparse = m.GetDense();
      break;
    case 2:
      *mat.Sparse = m.GetTriang();
      break;
    case 3:
      *mat.Sparse = m.GetSym();
      break;
    case 4:
      *mat.Sparse = m.GetSparse();
      break;
    case 5:
      *mat.Sparse = m.GetBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op= (const MySimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 5:
    switch (m.GetNum())
    {
    case 1:
      *mat.Banded = m.GetDense();
      break;
    case 2:
      *mat.Banded = m.GetTriang();
      break;
    case 3:
      *mat.Banded = m.GetSym();
      break;
    case 4:
      *mat.Banded = m.GetSparse();
      break;
    case 5:
      *mat.Banded = m.GetBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op= (const MySimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  default:
    SiconosMatrixException::selfThrow("op= (const MySimpleMatrix) : invalid type of matrix");
    break;
  }
  return *this;
}

const MySimpleMatrix& MySimpleMatrix::operator = (const MySiconosMatrix& m)
{
  switch (num)
  {
  case 1:
    switch (m.GetNum())
    {
    case 1:
      *mat.Dense = m.GetDense();
      break;
    case 2:
      *mat.Dense = m.GetTriang();
      break;
    case 3:
      *mat.Dense = m.GetSym();
      break;
    case 4:
      *mat.Dense = m.GetSparse();
      break;
    case 5:
      *mat.Dense = m.GetBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op= (const MySiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 2:
    switch (m.GetNum())
    {
    case 1:
      *mat.Triang = m.GetDense();
      break;
    case 2:
      *mat.Triang = m.GetTriang();
      break;
    case 3:
      *mat.Triang = m.GetSym();
      break;
    case 4:
      *mat.Triang = m.GetSparse();
      break;
    case 5:
      *mat.Triang = m.GetBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op= (const MySiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 3:
    switch (m.GetNum())
    {
    case 1:
      *mat.Sym = m.GetDense();
      break;
    case 2:
      *mat.Sym = m.GetTriang();
      break;
    case 3:
      *mat.Sym = m.GetSym();
      break;
    case 4:
      *mat.Sym = m.GetSparse();
      break;
    case 5:
      *mat.Sym = m.GetBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op= (const MySiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 4:
    switch (m.GetNum())
    {
    case 1:
      *mat.Sparse = m.GetDense();
      break;
    case 2:
      *mat.Sparse = m.GetTriang();
      break;
    case 3:
      *mat.Sparse = m.GetSym();
      break;
    case 4:
      *mat.Sparse = m.GetSparse();
      break;
    case 5:
      *mat.Sparse = m.GetBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op= (const MySimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 5:
    switch (m.GetNum())
    {
    case 1:
      *mat.Banded = m.GetDense();
      break;
    case 2:
      *mat.Banded = m.GetTriang();
      break;
    case 3:
      *mat.Banded = m.GetSym();
      break;
    case 4:
      *mat.Banded = m.GetSparse();
      break;
    case 5:
      *mat.Banded = m.GetBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op= (const MySimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  default:
    SiconosMatrixException::selfThrow("op= (const MySiconosMatrix) : invalid type of matrix");
    break;
  }
  return *this;
}

const MySimpleMatrix& MySimpleMatrix::operator += (const MySiconosMatrix& m)
{
  switch (num)
  {
  case 1:
    switch (m.GetNum())
    {
    case 1:
      *mat.Dense += m.GetDense();
      break;
    case 2:
      *mat.Dense += m.GetTriang();
      break;
    case 3:
      *mat.Dense += m.GetSym();
      break;
    case 4:
      *mat.Dense += m.GetSparse();
      break;
    case 5:
      *mat.Dense += m.GetBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op+= (const MySiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 2:
    switch (m.GetNum())
    {
    case 1:
      *mat.Triang += m.GetDense();
      break;
    case 2:
      *mat.Triang += m.GetTriang();
      break;
    case 3:
      *mat.Triang += m.GetSym();
      break;
    case 4:
      *mat.Triang += m.GetSparse();
      break;
    case 5:
      *mat.Triang += m.GetBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op+= (const MySiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 3:
    switch (m.GetNum())
    {
    case 1:
      *mat.Sym += m.GetDense();
      break;
    case 2:
      *mat.Sym += m.GetTriang();
      break;
    case 3:
      *mat.Sym += m.GetSym();
      break;
    case 4:
      *mat.Sym += m.GetSparse();
      break;
    case 5:
      *mat.Sym += m.GetBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op+= (const MySiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 4:
    switch (m.GetNum())
    {
    case 1:
      *mat.Sparse += m.GetDense();
      break;
    case 2:
      *mat.Sparse += m.GetTriang();
      break;
    case 3:
      *mat.Sparse += m.GetSym();
      break;
    case 4:
      *mat.Sparse += m.GetSparse();
      break;
    case 5:
      *mat.Sparse += m.GetBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op= (const MySimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 5:
    switch (m.GetNum())
    {
    case 1:
      *mat.Banded += m.GetDense();
      break;
    case 2:
      *mat.Banded += m.GetTriang();
      break;
    case 3:
      *mat.Banded += m.GetSym();
      break;
    case 4:
      *mat.Banded += m.GetSparse();
      break;
    case 5:
      *mat.Banded += m.GetBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op= (const MySimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  default:
    SiconosMatrixException::selfThrow("op+= (const MySiconosMatrix) : invalid type of matrix");
    break;
  }
  return *this;
}

const MySimpleMatrix& MySimpleMatrix::operator -= (const MySiconosMatrix& m)
{
  switch (num)
  {
  case 1:
    switch (m.GetNum())
    {
    case 1:
      *mat.Dense -= m.GetDense();
      break;
    case 2:
      *mat.Dense -= m.GetTriang();
      break;
    case 3:
      *mat.Dense -= m.GetSym();
      break;
    case 4:
      *mat.Dense -= m.GetSparse();
      break;
    case 5:
      *mat.Dense -= m.GetBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op-= (const MySiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 2:
    switch (m.GetNum())
    {
    case 1:
      *mat.Triang -= m.GetDense();
      break;
    case 2:
      *mat.Triang -= m.GetTriang();
      break;
    case 3:
      *mat.Triang -= m.GetSym();
      break;
    case 4:
      *mat.Triang -= m.GetSparse();
      break;
    case 5:
      *mat.Triang -= m.GetBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op-= (const MySiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 3:
    switch (m.GetNum())
    {
    case 1:
      *mat.Sym -= m.GetDense();
      break;
    case 2:
      *mat.Sym -= m.GetTriang();
      break;
    case 3:
      *mat.Sym -= m.GetSym();
      break;
    case 4:
      *mat.Sym -= m.GetSparse();
      break;
    case 5:
      *mat.Sym -= m.GetBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op-= (const MySiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 4:
    switch (m.GetNum())
    {
    case 1:
      *mat.Sparse -= m.GetDense();
      break;
    case 2:
      *mat.Sparse -= m.GetTriang();
      break;
    case 3:
      *mat.Sparse -= m.GetSym();
      break;
    case 4:
      *mat.Sparse -= m.GetSparse();
      break;
    case 5:
      *mat.Sparse -= m.GetBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op= (const MySimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 5:
    switch (m.GetNum())
    {
    case 1:
      *mat.Banded -= m.GetDense();
      break;
    case 2:
      *mat.Banded -= m.GetTriang();
      break;
    case 3:
      *mat.Banded -= m.GetSym();
      break;
    case 4:
      *mat.Banded -= m.GetSparse();
      break;
    case 5:
      *mat.Banded -= m.GetBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op= (const MySimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  default:
    SiconosMatrixException::selfThrow("op-= (const MySiconosMatrix) : invalid type of matrix");
    break;
  }
  return *this;
}

const MySimpleMatrix& MySimpleMatrix::operator *= (double m)
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

const MySimpleMatrix& MySimpleMatrix::operator *= (int m)
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

const MySimpleMatrix& MySimpleMatrix::operator /= (double m)
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

const MySimpleMatrix& MySimpleMatrix::operator /= (int m)
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


bool operator == (const MySiconosMatrix &m, const MySiconosMatrix &x)
{
  if (m.isBlock() ||  x.isBlock())
    SiconosMatrixException::selfThrow("op == (const MySiconosMatrix, const MySiconosMatrix) : incompatible type of matrix");

  double norm = (m - x).normInf();
  return (norm < tolerance);
}

MySimpleMatrix operator + (const MySiconosMatrix &x, const MySiconosMatrix &m)
{
  DenseMat p;
  TriangMat t;
  SymMat s;
  SparseMat sp;
  BandedMat b;

  if ((x.size1() != m.size1()) || (x.size2() != m.size2()))
    SiconosMatrixException::selfThrow("Matrix addition: inconsistent sizes");

  if (x.GetNum() == m.GetNum())
  {
    if (x.GetNum() == 1)
    {
      p = x.GetDense() + m.GetDense();
      return p;
    }
    else if (x.GetNum() == 2)
    {
      t = x.GetTriang() + m.GetTriang();
      return t;
    }
    else if (x.GetNum() == 3)
    {
      s = x.GetSym() + m.GetSym();
      return s;
    }
    else if (x.GetNum() == 4)
    {
      sp = x.GetSparse() + m.GetSparse();
      return sp;
    }
    else if (x.GetNum() == 5)
    {
      b.resize(m.size1(), m.size2(), (m.GetBanded()).lower(), (m.GetBanded()).upper(), false);
      b = x.GetBanded() + m.GetBanded();
      return b;
    }
  }
  else
    SiconosMatrixException::selfThrow("Matrix addition: use function add in order to add matrices of different type");
}

MySimpleMatrix operator - (const MySiconosMatrix &x, const MySiconosMatrix &m)
{
  DenseMat p;
  TriangMat t;
  SymMat s;
  SparseMat sp;
  BandedMat b;

  if ((x.size1() != m.size1()) || (x.size2() != m.size2()))
    SiconosMatrixException::selfThrow("Matrix subtraction: inconsistent sizes");

  if (x.GetNum() == m.GetNum())
  {
    if (x.GetNum() == 1)
    {
      p = x.GetDense() - m.GetDense();
      return p;
    }
    else if (x.GetNum() == 2)
    {
      t = x.GetTriang() - m.GetTriang();
      return t;
    }
    else if (x.GetNum() == 3)
    {
      s = x.GetSym() - m.GetSym();
      return s;
    }
    else if (x.GetNum() == 4)
    {
      sp = x.GetSparse() - m.GetSparse();
      return sp;
    }
    else if (x.GetNum() == 5)
    {
      b.resize(m.size1(), m.size2(), (m.GetBanded()).lower(), (m.GetBanded()).upper(), false);
      b = x.GetBanded() - m.GetBanded();
      return b;
    }
  }
  else
    SiconosMatrixException::selfThrow("Matrix subtraction: use function sub in order to subtract matrices of different type");
}

MySimpleMatrix operator * (const MySiconosMatrix &x, const MySiconosMatrix &m)
{
  DenseMat p;
  TriangMat t;
  SymMat s;
  SparseMat sp;
  BandedMat b;

  if ((x.size2() != m.size1()))
    SiconosMatrixException::selfThrow("Matrix product : inconsistent sizes");

  if (x.GetNum() == m.GetNum())
  {
    if (x.GetNum() == 1)
    {
      p = prod(x.GetDense(), m.GetDense());
      return p;
    }
    else if (x.GetNum() == 2)
    {
      t = prod(x.GetTriang(), m.GetTriang());
      return t;
    }
    else if (x.GetNum() == 3)
    {
      s = prod(x.GetSym(), m.GetSym());
      return s;
    }
    else if (x.GetNum() == 4)
    {
      sp = prod(x.GetSparse(), m.GetSparse());
      return sp;
    }
    else if (x.GetNum() == 5)
    {
      b.resize(x.size1(), m.size2(), (x.GetBanded()).lower(), (m.GetBanded()).upper(), false);
      b = prod(x.GetBanded(), m.GetBanded());
      return b;
    }
  }
  else
    SiconosMatrixException::selfThrow("Matrix product : use function prod in order to multiply matrices of different type");
}

MySimpleMatrix add(const MySiconosMatrix &x, const MySiconosMatrix& m)
{
  DenseMat p;

  if ((x.size1() != m.size1()) || (x.size2() != m.size2()))
    SiconosMatrixException::selfThrow("Matrix function add: inconsistent sizes");

  if (m.GetNum() == 1)
  {
    DenseMat q;
    q = m.GetDense();
    if (x.GetNum() == 1)
    {
      p = x.GetDense() + q;
    }
    else if (x.GetNum() == 2)
    {
      p = x.GetTriang() + q;
    }
    else if (x.GetNum() == 3)
    {
      p = x.GetSym() + q;
    }
    else if (x.GetNum() == 4)
    {
      p = x.GetSparse() + q;
    }
    else if (x.GetNum() == 5)
    {
      p = x.GetBanded() + q;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function add: invalid type of matrix");

  }
  else if (m.GetNum() == 2)
  {
    TriangMat q;
    q = m.GetTriang();
    if (x.GetNum() == 1)
    {
      p = x.GetDense() + q;
    }
    else if (x.GetNum() == 2)
    {
      p = x.GetTriang() + q;
    }
    else if (x.GetNum() == 3)
    {
      p = x.GetSym() + q;
    }
    else if (x.GetNum() == 4)
    {
      p = x.GetSparse() + q;
    }
    else if (x.GetNum() == 5)
    {
      p = x.GetBanded() + q;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function add: invalid type of matrix");

  }
  else if (m.GetNum() == 3)
  {
    SymMat q;
    q = m.GetSym();
    if (x.GetNum() == 1)
    {
      p = x.GetDense() + q;
    }
    else if (x.GetNum() == 2)
    {
      p = x.GetTriang() + q;
    }
    else if (x.GetNum() == 3)
    {
      p = x.GetSym() + q;
    }
    else if (x.GetNum() == 4)
    {
      p = x.GetSparse() + q;
    }
    else if (x.GetNum() == 5)
    {
      p = x.GetBanded() + q;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function add: invalid type of matrix");

  }
  else if (m.GetNum() == 4)
  {
    SparseMat q;
    q = m.GetSparse();
    if (x.GetNum() == 1)
    {
      p = x.GetDense() + q;
    }
    else if (x.GetNum() == 2)
    {
      p = x.GetTriang() + q;
    }
    else if (x.GetNum() == 3)
    {
      p = x.GetSym() + q;
    }
    else if (x.GetNum() == 4)
    {
      p = x.GetSparse() + q;
    }
    else if (x.GetNum() == 5)
    {
      p = x.GetBanded() + q;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function add: invalid type of matrix");

  }
  else if (m.GetNum() == 5)
  {
    BandedMat q;
    q = m.GetBanded();
    if (x.GetNum() == 1)
    {
      p = x.GetDense() + q;
    }
    else if (x.GetNum() == 2)
    {
      p = x.GetTriang() + q;
    }
    else if (x.GetNum() == 3)
    {
      p = x.GetSym() + q;
    }
    else if (x.GetNum() == 4)
    {
      p = x.GetSparse() + q;
    }
    else if (x.GetNum() == 5)
    {
      p = x.GetBanded() + q;
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

MySimpleMatrix sub(const MySiconosMatrix &x, const MySiconosMatrix& m)
{
  DenseMat p;

  if ((x.size1() != m.size1()) || (x.size2() != m.size2()))
    SiconosMatrixException::selfThrow("Matrix function sub: inconsistent sizes");

  if (m.GetNum() == 1)
  {
    DenseMat q;
    q = m.GetDense();
    if (x.GetNum() == 1)
    {
      p = x.GetDense() - q;
    }
    else if (x.GetNum() == 2)
    {
      p = x.GetTriang() - q;
    }
    else if (x.GetNum() == 3)
    {
      p = x.GetSym() - q;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function sub: invalid type of matrix");

  }
  else if (m.GetNum() == 2)
  {
    TriangMat q;
    q = m.GetTriang();
    if (x.GetNum() == 1)
    {
      p = x.GetDense() - q;
    }
    else if (x.GetNum() == 2)
    {
      p = x.GetTriang() - q;
    }
    else if (x.GetNum() == 3)
    {
      p = x.GetSym() - q;
    }
    else if (x.GetNum() == 4)
    {
      p = x.GetSparse() - q;
    }
    else if (x.GetNum() == 5)
    {
      p = x.GetBanded() - q;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function sub: invalid type of matrix");

  }
  else if (m.GetNum() == 3)
  {
    SymMat q;
    q = m.GetSym();
    if (x.GetNum() == 1)
    {
      p = x.GetDense() - q;
    }
    else if (x.GetNum() == 2)
    {
      p = x.GetTriang() - q;
    }
    else if (x.GetNum() == 3)
    {
      p = x.GetSym() - q;
    }
    else if (x.GetNum() == 4)
    {
      p = x.GetSparse() - q;
    }
    else if (x.GetNum() == 5)
    {
      p = x.GetBanded() - q;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function sub: invalid type of matrix");

  }
  else if (m.GetNum() == 4)
  {
    SparseMat q;
    q = m.GetSparse();
    if (x.GetNum() == 1)
    {
      p = x.GetDense() - q;
    }
    else if (x.GetNum() == 2)
    {
      p = x.GetTriang() - q;
    }
    else if (x.GetNum() == 3)
    {
      p = x.GetSym() - q;
    }
    else if (x.GetNum() == 4)
    {
      p = x.GetSparse() - q;
    }
    else if (x.GetNum() == 5)
    {
      p = x.GetBanded() - q;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function add: invalid type of matrix");

  }
  else if (m.GetNum() == 5)
  {
    BandedMat q;
    q = m.GetBanded();
    if (x.GetNum() == 1)
    {
      p = x.GetDense() - q;
    }
    else if (x.GetNum() == 2)
    {
      p = x.GetTriang() - q;
    }
    else if (x.GetNum() == 3)
    {
      p = x.GetSym() - q;
    }
    else if (x.GetNum() == 4)
    {
      p = x.GetSparse() - q;
    }
    else if (x.GetNum() == 5)
    {
      p = x.GetBanded() - q;
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

MySimpleMatrix prod(const MySiconosMatrix &x, const MySiconosMatrix& m)
{
  DenseMat p;
  TriangMat t;
  SymMat s;
  SparseMat sp;
  BandedMat b;

  if ((x.size2() != m.size1()))
    SiconosMatrixException::selfThrow("Matrix function prod : inconsistent sizes");

  if (m.GetNum() == 1)
  {
    DenseMat q;
    q = m.GetDense();
    if (x.GetNum() == 1)
    {
      p = prod(x.GetDense(), q);
    }
    else if (x.GetNum() == 2)
    {
      p = prod(x.GetTriang(), q);
    }
    else if (x.GetNum() == 3)
    {
      p = prod(x.GetSym(), q);
    }
    else if (x.GetNum() == 4)
    {
      p = prod(x.GetSparse(), q);
    }
    else if (x.GetNum() == 5)
    {
      p = prod(x.GetBanded(), q);
    }
    else
      SiconosMatrixException::selfThrow("Matrix function prod: invalid type of matrix");

    return p;
  }
  else if (m.GetNum() == 2)
  {
    TriangMat q;
    q = m.GetTriang();
    if (x.GetNum() == 1)
    {
      p = prod(x.GetDense(), q);
      return p;
    }
    else if (x.GetNum() == 2)
    {
      t = prod(x.GetTriang(), q);
      return t;
    }
    else if (x.GetNum() == 3)
    {
      p = prod(x.GetSym(), q);
      return p;
    }
    else if (x.GetNum() == 4)
    {
      p = prod(x.GetSparse(), q);
      return p;
    }
    else if (x.GetNum() == 5)
    {
      p = prod(x.GetBanded(), q);
      return p;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function prod: invalid type of matrix");

  }
  else if (m.GetNum() == 3)
  {
    SymMat q;
    q = m.GetSym();
    if (x.GetNum() == 1)
    {
      p = prod(x.GetDense(), q);
      return p;
    }
    else if (x.GetNum() == 2)
    {
      p = prod(x.GetTriang(), q);
      return p;
    }
    else if (x.GetNum() == 3)
    {
      s = prod(x.GetSym(), q);
      return s;
    }
    else if (x.GetNum() == 4)
    {
      p = prod(x.GetSparse(), q);
      return p;
    }
    else if (x.GetNum() == 5)
    {
      p = prod(x.GetBanded(), q);
      return p;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function prod: invalid type of matrix");

  }
  else if (m.GetNum() == 4)
  {
    SparseMat q;
    q = m.GetSparse();
    if (x.GetNum() == 1)
    {
      p = prod(x.GetDense(), q);
      return p;
    }
    else if (x.GetNum() == 2)
    {
      t = prod(x.GetTriang(), q);
      return t;
    }
    else if (x.GetNum() == 3)
    {
      p = prod(x.GetSym(), q);
      return p;
    }
    else if (x.GetNum() == 4)
    {
      p = prod(x.GetSparse(), q);
      return p;
    }
    else if (x.GetNum() == 5)
    {
      p = prod(x.GetBanded(), q);
      return p;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function prod: invalid type of matrix");

  }
  else if (m.GetNum() == 5)
  {
    BandedMat q;
    q = m.GetBanded();
    if (x.GetNum() == 1)
    {
      p = prod(x.GetDense(), q);
      return p;
    }
    else if (x.GetNum() == 2)
    {
      t = prod(x.GetTriang(), q);
      return t;
    }
    else if (x.GetNum() == 3)
    {
      p = prod(x.GetSym(), q);
      return p;
    }
    else if (x.GetNum() == 4)
    {
      p = prod(x.GetSparse(), q);
      return p;
    }
    else if (x.GetNum() == 5)
    {
      p = prod(x.GetBanded(), q);
      return p;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function prod: invalid type of matrix");

  }
  else
  {
    SiconosMatrixException::selfThrow("Matrix function prod: invalid type of matrix");
  }

}

MySimpleMatrix multranspose(const MySiconosMatrix &x, const MySiconosMatrix &m)
{
  return prod(x, trans(m));
}

MySimpleMatrix operator * (const MySiconosMatrix &m, double d)
{
  DenseMat p;
  TriangMat t;
  SymMat s;
  SparseMat sp;
  BandedMat b;
  if (m.GetNum() == 1)
  {
    p = m.GetDense() * d;
    return p;
  }
  else if (m.GetNum() == 2)
  {
    t = m.GetTriang() * d;
    return t;
  }
  else if (m.GetNum() == 3)
  {
    s = m.GetSym() * d;
    return s;
  }
  else if (m.GetNum() == 4)
  {
    sp = m.GetSparse() * d;
    return sp;
  }
  else if (m.GetNum() == 5)
  {
    b.resize(m.size1(), m.size2(), (m.GetBanded()).lower(), (m.GetBanded()).upper(), false);
    b = m.GetBanded() * d;
    return b;
  }
  else
    SiconosMatrixException::selfThrow("Matrix op * (const MySiconosMatrix&, double): invalid type of matrix");
}

MySimpleMatrix operator * (const MySiconosMatrix &m, int d)
{
  DenseMat p;
  TriangMat t;
  SymMat s;
  SparseMat sp;
  BandedMat b;
  if (m.GetNum() == 1)
  {
    p = m.GetDense() * d;
    return p;
  }
  else if (m.GetNum() == 2)
  {
    t = m.GetTriang() * d;
    return t;
  }
  else if (m.GetNum() == 3)
  {
    s = m.GetSym() * d;
    return s;
  }
  else if (m.GetNum() == 4)
  {
    sp = m.GetSparse() * d;
    return sp;
  }
  else if (m.GetNum() == 5)
  {
    b = m.GetBanded() * d;
    return b;
  }
  else
    SiconosMatrixException::selfThrow("Matrix op * (const MySiconosMatrix&, int): invalid type of matrix");
}

MySimpleMatrix operator * (double d, const MySiconosMatrix &m)
{
  DenseMat p;
  TriangMat t;
  SymMat s;
  SparseMat sp;
  BandedMat b;
  if (m.GetNum() == 1)
  {
    p = d * m.GetDense();
    return p;
  }
  else if (m.GetNum() == 2)
  {
    t = d * m.GetTriang();
    return t;
  }
  else if (m.GetNum() == 3)
  {
    s = d * m.GetSym();
    return s;
  }
  else if (m.GetNum() == 4)
  {
    sp = d * m.GetSparse();
    return sp;
  }
  else if (m.GetNum() == 5)
  {
    b = d * m.GetBanded();
    return b;
  }
  else
    SiconosMatrixException::selfThrow("Matrix op * (double, const MySiconosMatrix&): invalid type of matrix");
}

MySimpleMatrix operator * (int d, const MySiconosMatrix &m)
{
  DenseMat p;
  TriangMat t;
  SymMat s;
  SparseMat sp;
  BandedMat b;
  if (m.GetNum() == 1)
  {
    p = d * m.GetDense();
    return p;
  }
  else if (m.GetNum() == 2)
  {
    t = d * m.GetTriang();
    return t;
  }
  else if (m.GetNum() == 3)
  {
    s = d * m.GetSym();
    return s;
  }
  else if (m.GetNum() == 4)
  {
    sp = d * m.GetSparse();
    return sp;
  }
  else if (m.GetNum() == 5)
  {
    b = d * m.GetBanded();
    return b;
  }
  else
    SiconosMatrixException::selfThrow("Matrix op * (int, const MySiconosMatrix&): invalid type of matrix");
}

MySimpleMatrix operator / (const MySiconosMatrix &m, double d)
{
  DenseMat p;
  TriangMat t;
  SymMat s;
  SparseMat sp;
  BandedMat b;
  if (m.GetNum() == 1)
  {
    p = m.GetDense() / d;
    return p;
  }
  else if (m.GetNum() == 2)
  {
    t = m.GetTriang() / d;
    return t;
  }
  else if (m.GetNum() == 3)
  {
    s = m.GetSym() / d;
    return s;
  }
  else if (m.GetNum() == 4)
  {
    sp = m.GetSparse() / d;
    return sp;
  }
  else if (m.GetNum() == 5)
  {
    b = m.GetBanded() / d;
    return b;
  }
  else
    SiconosMatrixException::selfThrow("Matrix op / (const MySiconosMatrix&, double): invalid type of matrix");
}

MySimpleMatrix operator / (const MySiconosMatrix &m, int d)
{
  DenseMat p;
  TriangMat t;
  SymMat s;
  SparseMat sp;
  BandedMat b;
  if (m.GetNum() == 1)
  {
    p = m.GetDense() / d;
    return p;
  }
  else if (m.GetNum() == 2)
  {
    t = m.GetTriang() / d;
    return t;
  }
  else if (m.GetNum() == 3)
  {
    s = m.GetSym() / d;
    return s;
  }
  else if (m.GetNum() == 4)
  {
    sp = m.GetSparse() / d;
    return sp;
  }
  else if (m.GetNum() == 5)
  {
    b = m.GetBanded() / d;
    return b;
  }
  else
    SiconosMatrixException::selfThrow("Matrix op / (const MySiconosMatrix&, int): invalid type of matrix");
}



