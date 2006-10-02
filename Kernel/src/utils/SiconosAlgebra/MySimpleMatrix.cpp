#include "MySimpleMatrix.h"
#include "ioMatrix.h"
#include <cassert>

// Build a Simple Matrix from its type (ie DENSE, TRIANGULAR, BANDED, SPARSE or SYMMETRIC)
// => default (private) constructor
MySimpleMatrix::MySimpleMatrix(TYP typ)
{
  setIsBlock(false);
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
/***************************** CONSTRUCTORS ****************************/

// Copy constructors
MySimpleMatrix::MySimpleMatrix(const MySimpleMatrix &smat)
{
  setIsBlock(false);
  if (smat.getNum() == 1)
  {
    mat.Dense = new DenseMat(smat.getDense());
    num = 1;
  }
  else if (smat.getNum() == 2)
  {
    mat.Triang = new TriangMat(smat.getTriang());
    num = 2;
  }
  else if (smat.getNum() == 3)
  {
    mat.Sym = new SymMat(smat.getSym());
    num = 3;
  }
  else if (smat.getNum() == 4)
  {
    mat.Sparse = new SparseMat(smat.getSparse());
    num = 4;
  }
  else if (smat.getNum() == 5)
  {
    mat.Banded = new BandedMat(smat.getBanded());
    num = 5;
  }
  else
    SiconosMatrixException::selfThrow("constructor(const MySimpleMatrix) : invalid parameter given");
}

MySimpleMatrix::MySimpleMatrix(const MySiconosMatrix &smat)
{
  setIsBlock(false);
  assert(smat.isBlock() == false);
  if (smat.getNum() == 1)
  {
    mat.Dense = new DenseMat(smat.getDense());
    num = 1;
  }
  else if (smat.getNum() == 2)
  {
    mat.Triang = new TriangMat(smat.getTriang());
    num = 2;
  }
  else if (smat.getNum() == 3)
  {
    mat.Sym = new SymMat(smat.getSym());
    num = 3;
  }
  else if (smat.getNum() == 4)
  {
    mat.Sparse = new SparseMat(smat.getSparse());
    num = 4;
  }
  else if (smat.getNum() == 5)
  {
    mat.Banded = new BandedMat(smat.getBanded());
    num = 5;
  }
  else
    SiconosMatrixException::selfThrow("constructor(const MySiconosMatrix) : invalid parameter given");
}

MySimpleMatrix::MySimpleMatrix(unsigned int row, unsigned int col, TYP typ, unsigned int upper, unsigned int lower)
{
  setIsBlock(false);
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

MySimpleMatrix::MySimpleMatrix(const DenseMat& m)
{
  setIsBlock(false);
  mat.Dense = new DenseMat(m);
  num = 1;
}

MySimpleMatrix::MySimpleMatrix(const TriangMat& m)
{
  setIsBlock(false);
  mat.Triang = new TriangMat(m);
  num = 2;
}

MySimpleMatrix::MySimpleMatrix(const SymMat& m)
{
  setIsBlock(false);
  mat.Sym = new SymMat(m);
  num = 3;
}

MySimpleMatrix::MySimpleMatrix(const SparseMat& m)
{
  setIsBlock(false);
  mat.Sparse = new SparseMat(m);
  num = 4;
}

MySimpleMatrix::MySimpleMatrix(const BandedMat& m)
{
  setIsBlock(false);
  mat.Banded = new BandedMat(m);
  num = 5;
}

MySimpleMatrix::MySimpleMatrix(const std::vector<double> &v, unsigned int row, unsigned int col, TYP typ, unsigned int lower, unsigned int upper)
{

  if (((v.size()) != (unsigned int)row * col && (typ != SYMMETRIC && typ != BANDED)) || (v.size() != (unsigned int)row * row && typ == SYMMETRIC) || (typ == BANDED && ((v.size()) != (unsigned int)(std::max)(row, col) * (lower + 1 + upper))))
    SiconosMatrixException::selfThrow("constructor(TYP, const std::vector<double>, int, int) : invalid vector size");

  setIsBlock(false);
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
  setIsBlock(false);
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
unsigned int  MySimpleMatrix::getNum(void)const
{
  return num;
}
void MySimpleMatrix::setNum(unsigned int n)
{
  num = n;
}
unsigned int  MySimpleMatrix::size1(void)const
{
  unsigned int size1;
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
unsigned int  MySimpleMatrix::size2(void)const
{
  unsigned int size2;
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

void MySimpleMatrix::resize(unsigned int row, unsigned int col, unsigned int lower, unsigned int upper, bool preserve)
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

const DenseMat MySimpleMatrix::getDense(unsigned int row, unsigned int col)const
{

  if (num != 1)
    SiconosMatrixException::selfThrow("DenseMat getDense(unsigned int row, unsigned int col) : the current matrix is not a Dense matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("DenseMat getDense(unsigned int row, unsigned int col) : row or col not equal to 0.0");

  return *mat.Dense;
}

const TriangMat MySimpleMatrix::getTriang(unsigned int row, unsigned int col)const
{

  if (num != 2)
    SiconosMatrixException::selfThrow("TriangMat getTriang(unsigned int row, unsigned int col) : the current matrix is not a Triangular matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("TriangMat getTriang(unsigned int row, unsigned int col) : row or col not equal to 0.0");

  return *mat.Triang;
}

const SymMat MySimpleMatrix::getSym(unsigned int row, unsigned int col)const
{

  if (num != 3)
    SiconosMatrixException::selfThrow("SymMat getSym(unsigned int row, unsigned int col) : the current matrix is not a Symmetric matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("SymMat getSym(unsigned int row, unsigned int col) : row or col not equal to 0.0");

  return *mat.Sym;
}

const SparseMat MySimpleMatrix::getSparse(unsigned int row, unsigned int col)const
{

  if (num != 4)
    SiconosMatrixException::selfThrow("SparseMat getSparse(unsigned int row, unsigned int col) : the current matrix is not a Sparse matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("SparseMat getSparse(unsigned int row, unsigned int col) : row or col not equal to 0.0");

  return *mat.Sparse;
}

const BandedMat MySimpleMatrix::getBanded(unsigned int row, unsigned int col)const
{

  if (num != 5)
    SiconosMatrixException::selfThrow("BandedMat getBanded(unsigned int row, unsigned int col) : the current matrix is not a Banded matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("BandedMat getBanded(unsigned int row, unsigned int col) : row or col not equal to 0.0");

  return *mat.Banded;
}

const DenseMat* MySimpleMatrix::getDensePtr(unsigned int row, unsigned int col)const
{

  if (num != 1)
    SiconosMatrixException::selfThrow("DenseMat* getDensePtr(unsigned int row, unsigned int col) : the current matrix is not a Dense matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("DenseMat* getDensePtr(unsigned int row, unsigned int col) : row or col not equal to 0.0");

  return mat.Dense;
}

const TriangMat* MySimpleMatrix::getTriangPtr(unsigned int row, unsigned int col)const
{

  if (num != 2)
    SiconosMatrixException::selfThrow("TriangMat* getTriangPtr(unsigned int row, unsigned int col) : the current matrix is not a Triangular matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("TriangMat* getTriangPtr(unsigned int row, unsigned int col) : row or col not equal to 0.0");

  return mat.Triang;
}

const SymMat* MySimpleMatrix::getSymPtr(unsigned int row, unsigned int col)const
{

  if (num != 3)
    SiconosMatrixException::selfThrow("SymMat* getSymPtr(unsigned int row, unsigned int col) : the current matrix is not a Symmetric matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("SymMat* getSymPtr(unsigned int row, unsigned int col) : row or col not equal to 0.0");

  return mat.Sym;
}

const SparseMat* MySimpleMatrix::getSparsePtr(unsigned int row, unsigned int col)const
{

  if (num != 4)
    SiconosMatrixException::selfThrow("SparseMat* getSparsePtr(unsigned int row, unsigned int col) : the current matrix is not a Sparse matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("SparseMat* getSparsePtr(unsigned int row, unsigned int col) : row or col not equal to 0.0");

  return mat.Sparse;
}

const BandedMat* MySimpleMatrix::getBandedPtr(unsigned int row, unsigned int col)const
{

  if (num != 5)
    SiconosMatrixException::selfThrow("BandedMat* getBandedPtr(unsigned int row, unsigned int col) : the current matrix is not a Banded matrix");

  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("BandedMat* getBandedPtr(unsigned int row, unsigned int col) : row or col not equal to 0.0");

  return mat.Banded;
}

const mapped MySimpleMatrix::getMap(void)const
{
  SiconosMatrixException::selfThrow("mapped getMap : getMap is forbidden for MySimpleMatrix");
}

void MySimpleMatrix::blockMatrixCopy(const MySiconosMatrix &m, unsigned int i, unsigned int j)
{

  if (num != 1)
    SiconosMatrixException::selfThrow("blockMatriCopy : the current matrix is not dense, a copy into its data may change its type.");

  if (i >= size1() || i < 0)
    SiconosMatrixException::selfThrow("blockMatriCopy : row_min given is out of range");

  if (j >= size2() || j < 0)
    SiconosMatrixException::selfThrow("blockMatriCopy : col_min given is out of range");

  unsigned int num2 = m.getNum();

  //if(num==1){
  if (num2 == 1)
  {
    unsigned int row_max = i + (m.getDense()).size1();
    unsigned int col_max = j + (m.getDense()).size2();
    if (row_max > size1() || row_max < 0)
      SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

    if (col_max > size2() || col_max < 0)
      SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

    subrange(*mat.Dense, i, i + (m.getDense()).size1(), j, j + (m.getDense()).size2()) = m.getDense();
  }
  else if (num2 == 2)
  {
    unsigned int row_max = i + (m.getTriang()).size1();
    unsigned int col_max = j + (m.getTriang()).size2();
    if (row_max > size1() || row_max < 0)
      SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

    if (col_max > size2() || col_max < 0)
      SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

    subrange(*mat.Dense, i, i + (m.getTriang()).size1(), j, j + (m.getTriang()).size2()) = m.getTriang();
  }
  else if (num2 == 3)
  {
    unsigned int row_max = i + (m.getSym()).size1();
    unsigned int col_max = j + (m.getSym()).size2();
    if (row_max > size1() || row_max < 0)
      SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

    if (col_max > size2() || col_max < 0)
      SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

    subrange(*mat.Dense, i, i + (m.getSym()).size1(), j, j + (m.getSym()).size2()) = m.getSym();
  }
  else if (num2 == 4)
  {
    unsigned int row_max = i + (m.getSparse()).size1();
    unsigned int col_max = j + (m.getSparse()).size2();
    if (row_max > size1() || row_max < 0)
      SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

    if (col_max > size2() || col_max < 0)
      SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

    subrange(*mat.Dense, i, i + (m.getSparse()).size1(), j, j + (m.getSparse()).size2()) = m.getSparse();
  }
  else if (num2 == 5)
  {
    unsigned int row_max = i + (m.getBanded()).size1();
    unsigned int col_max = j + (m.getBanded()).size2();
    if (row_max > size1() || row_max < 0)
      SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

    if (col_max > size2() || col_max < 0)
      SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

    subrange(*mat.Dense, i, i + (m.getBanded()).size1(), j, j + (m.getBanded()).size2()) = m.getBanded();
  }
  //  }
  //   else if(num==2){
  //     if(num2==1){
  //       int row_max = i + (m.getDense()).size1();
  //       int col_max = j + (m.getDense()).size2();
  //       if(row_max > size1() || row_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       if(col_max > size2() || col_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       subrange(*mat.Triang, i, i+(m.getDense ()).size1 (), j, j+(m.getDense ()).size2 ()) = m.getDense ();
  //     }
  //     else if(num2==2){
  //       int row_max = i + (m.getTriang()).size1();
  //       int col_max = j + (m.getTriang()).size2();
  //       if(row_max > size1() || row_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       if(col_max > size2() || col_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       subrange(*mat.Triang, i, i+(m.getTriang ()).size1 (), j, j+(m.getTriang ()).size2 ()) = m.getTriang();
  //     }
  //     else if(num2==3){
  //       int row_max = i + (m.getSym()).size1();
  //       int col_max = j + (m.getSym()).size2();
  //       if(row_max > size1() || row_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       if(col_max > size2() || col_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       subrange(*mat.Triang, i, i+(m.getSym ()).size1 (), j, j+(m.getSym ()).size2 ()) = m.getSym ();
  //     }
  //     else if(num2==4){
  //       int row_max = i + (m.getSparse()).size1();
  //       int col_max = j + (m.getSparse()).size2();
  //       if(row_max > size1() || row_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       if(col_max > size2() || col_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       subrange(*mat.Triang, i, i+(m.getSparse ()).size1 (), j, j+(m.getSparse ()).size2 ()) = m.getSparse ();
  //     }
  //     else if(num2==5){
  //       int row_max = i + (m.getBanded()).size1();
  //       int col_max = j + (m.getBanded()).size2();
  //       if(row_max > size1() || row_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       if(col_max > size2() || col_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       subrange(*mat.Triang, i, i+(m.getBanded ()).size1 (), j, j+(m.getBanded ()).size2 ()) = m.getBanded ();
  //     }
  //   }
  //   else if(num==3){
  //     if(num2==1){
  //       int row_max = i + (m.getDense()).size1();
  //       int col_max = j + (m.getDense()).size2();
  //       if(row_max > size1() || row_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       if(col_max > size2() || col_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       subrange(*mat.Sym, i, i+(m.getDense ()).size1 (), j, j+(m.getDense ()).size2 ()) = m.getDense ();
  //     }
  //     else if(num2==2){
  //       int row_max = i + (m.getTriang()).size1();
  //       int col_max = j + (m.getTriang()).size2();
  //       if(row_max > size1() || row_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       if(col_max > size2() || col_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       subrange(*mat.Sym, i, i+(m.getTriang ()).size1 (), j, j+(m.getTriang ()).size2 ()) = m.getTriang ();
  //     }
  //     else if(num2==3){
  //       int row_max = i + (m.getSym()).size1();
  //       int col_max = j + (m.getSym()).size2();
  //       if(row_max > size1() || row_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       if(col_max > size2() || col_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       subrange(*mat.Sym, i, i+(m.getSym ()).size1 (), j, j+(m.getSym ()).size2 ()) = m.getSym ();
  //     }
  //     else if(num2==4){
  //       int row_max = i + (m.getSparse()).size1();
  //       int col_max = j + (m.getSparse()).size2();
  //       if(row_max > size1() || row_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       if(col_max > size2() || col_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       subrange(*mat.Sym, i, i+(m.getSparse ()).size1 (), j, j+(m.getSparse ()).size2 ()) = m.getSparse ();
  //     }
  //     else if(num2==5){
  //       int row_max = i + (m.getBanded()).size1();
  //       int col_max = j + (m.getBanded()).size2();
  //       if(row_max > size1() || row_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       if(col_max > size2() || col_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       subrange(*mat.Sym, i, i+(m.getBanded ()).size1 (), j, j+(m.getBanded ()).size2 ()) = m.getBanded ();
  //     }
  //   }
  //   if(num==4){
  //     if(num2==1){
  //       int row_max = i + (m.getDense()).size1();
  //       int col_max = j + (m.getDense()).size2();
  //       if(row_max > size1() || row_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       if(col_max > size2() || col_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       subrange(*mat.Sparse, i, i+(m.getDense ()).size1 (), j, j+(m.getDense ()).size2 ()) = m.getDense ();
  //     }
  //     else if(num2==2){
  //       int row_max = i + (m.getTriang()).size1();
  //       int col_max = j + (m.getTriang()).size2();
  //       if(row_max > size1() || row_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       if(col_max > size2() || col_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       subrange(*mat.Sparse, i, i+(m.getTriang ()).size1 (), j, j+(m.getTriang ()).size2 ()) = m.getTriang ();
  //     }
  //     else if(num2==3){
  //       int row_max = i + (m.getSym()).size1();
  //       int col_max = j + (m.getSym()).size2();
  //       if(row_max > size1() || row_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       if(col_max > size2() || col_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       subrange(*mat.Sparse, i, i+(m.getSym ()).size1 (), j, j+(m.getSym ()).size2 ()) = m.getSym ();
  //     }
  //     else if(num2==4){
  //       int row_max = i + (m.getSparse()).size1();
  //       int col_max = j + (m.getSparse()).size2();
  //       if(row_max > size1() || row_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       if(col_max > size2() || col_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       subrange(*mat.Sparse, i, i+(m.getSparse ()).size1 (), j, j+(m.getSparse ()).size2 ()) = m.getSparse ();
  //     }
  //     else if(num2==5){
  //       int row_max = i + (m.getBanded()).size1();
  //       int col_max = j + (m.getBanded()).size2();
  //       if(row_max > size1() || row_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       if(col_max > size2() || col_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       subrange(*mat.Sparse, i, i+(m.getBanded ()).size1 (), j, j+(m.getBanded ()).size2 ()) = m.getBanded ();
  //     }
  //   }
  //   if(num==5){
  //     if(num2==1){
  //       int row_max = i + (m.getDense()).size1();
  //       int col_max = j + (m.getDense()).size2();
  //       if(row_max > size1() || row_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       if(col_max > size2() || col_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       subrange(*mat.Banded, i, i+(m.getDense ()).size1 (), j, j+(m.getDense ()).size2 ()) = m.getDense ();
  //     }
  //     else if(num2==2){
  //       int row_max = i + (m.getTriang()).size1();
  //       int col_max = j + (m.getTriang()).size2();
  //       if(row_max > size1() || row_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       if(col_max > size2() || col_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       subrange(*mat.Banded, i, i+(m.getTriang ()).size1 (), j, j+(m.getTriang ()).size2 ()) = m.getTriang ();
  //     }
  //     else if(num2==3){
  //       int row_max = i + (m.getSym()).size1();
  //       int col_max = j + (m.getSym()).size2();
  //       if(row_max > size1() || row_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       if(col_max > size2() || col_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       subrange(*mat.Banded, i, i+(m.getSym ()).size1 (), j, j+(m.getSym ()).size2 ()) = m.getSym ();
  //     }
  //     else if(num2==4){
  //       int row_max = i + (m.getSparse()).size1();
  //       int col_max = j + (m.getSparse()).size2();
  //       if(row_max > size1() || row_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       if(col_max > size2() || col_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       subrange(*mat.Banded, i, i+(m.getSparse ()).size1 (), j, j+(m.getSparse ()).size2 ()) = m.getSparse ();
  //     }
  //     else if(num2==5){
  //       int row_max = i + (m.getBanded()).size1();
  //       int col_max = j + (m.getBanded()).size2();
  //       if(row_max > size1() || row_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       if(col_max > size2() || col_max<0)
  //  SiconosMatrixException::selfThrow("blockMatriCopy : inconsistent sizes");

  //       subrange(*mat.Banded, i, i+(m.getBanded ()).size1 (), j, j+(m.getBanded ()).size2 ()) = m.getBanded ();
  //     }
  //   }
}

void MySimpleMatrix::getBlock(unsigned int row_min, unsigned int col_min, MySiconosMatrix &m)const
{

  // We only accept dense matrix for m.
  if (m.getNum() != 1)
    SiconosMatrixException::selfThrow("getBlock(i,j,m) : m must be a dense matrix.");

  if (row_min >= size1() || row_min < 0)
    SiconosMatrixException::selfThrow("getBlock : row_min given is out of range");

  if (col_min >= size2() || row_min < 0)
    SiconosMatrixException::selfThrow("getBlock : col_min given is out of range");
  DenseMat q;
  unsigned int row_max, col_max;
  if (num == 1)
  {
    row_max = m.getDense().size1() + row_min;
    col_max = m.getDense().size2() + col_min;

    if (row_max > size1() || row_max < 0)
      SiconosMatrixException::selfThrow("getBlock : inconsistent sizes");

    if (col_max > size2() || col_max < 0)
      SiconosMatrixException::selfThrow("getBlock : inconsistent sizes");

    q = subrange(*mat.Dense, row_min, row_max, col_min, col_max);
  }
  else if (num == 2)
  {
    row_max = m.getDense().size1() + row_min;
    col_max = m.getDense().size2() + col_min;

    if (row_max > size1() || row_max < 0)
      SiconosMatrixException::selfThrow("getBlock : inconsistent sizes");

    if (col_max > size2() || col_max < 0)
      SiconosMatrixException::selfThrow("getBlock : inconsistent sizes");
    q = subrange(*mat.Triang, row_min, row_max, col_min, col_max);
  }
  else if (num == 3)
  {
    row_max = m.getDense().size1() + row_min;
    col_max = m.getDense().size2() + col_min;

    if (row_max > size1() || row_max < 0)
      SiconosMatrixException::selfThrow("getBlock : inconsistent sizes");

    if (col_max > size2() || col_max < 0)
      SiconosMatrixException::selfThrow("getBlock : inconsistent sizes");
    q = subrange(*mat.Sym, row_min, row_max, col_min, col_max);
  }
  else if (num == 4)
  {
    row_max = m.getDense().size1() + row_min;
    col_max = m.getDense().size2() + col_min;

    if (row_max > size1() || row_max < 0)
      SiconosMatrixException::selfThrow("getBlock : inconsistent sizes");

    if (col_max > size2() || col_max < 0)
      SiconosMatrixException::selfThrow("getBlock : inconsistent sizes");
    q = subrange(*mat.Sparse, row_min, row_max, col_min, col_max);
  }
  else if (num == 5)
  {
    row_max = m.getDense().size1() + row_min;
    col_max = m.getDense().size2() + col_min;

    if (row_max > size1() || row_max < 0)
      SiconosMatrixException::selfThrow("getBlock : inconsistent sizes");

    if (col_max > size2() || col_max < 0)
      SiconosMatrixException::selfThrow("getBlock : inconsistent sizes");
    q = subrange(*mat.Banded, row_min, row_max, col_min, col_max);
  }
  MySimpleMatrix p(q);
  m = p;
}

const std::deque<bool> MySimpleMatrix::getBlockAllocated(void)const
{
  SiconosMatrixException::selfThrow("std::deque<bool> getBlockAllocated : getBlockAllocated is forbidden for MySimpleMatrix");
}

void MySimpleMatrix::getRow(unsigned int r, MySimpleVector &vect)const
{

  if (r >= size1() || r < 0)
    SiconosMatrixException::selfThrow("getRow : row is out of range");

  if (vect.size() != size2())
    SiconosMatrixException::selfThrow("getRow : inconsistent sizes");

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

void MySimpleMatrix::setRow(unsigned int r, const MySimpleVector &vect)
{

  if (r >= size1() || r < 0)
    SiconosMatrixException::selfThrow("setRow : row is out of range");

  if (vect.size() != size2())
    SiconosMatrixException::selfThrow("setRow : inconsistent sizes");

  if (num == 1)
  {
    if (vect.getNum() == 1)
    {
      row(*mat.Dense, r) = vect.getDense();
    }
    else if (vect.getNum() == 2)
    {
      row(*mat.Dense, r) = vect.getSparse();
    }
  }
  else if (num == 2)
  {
    if (vect.getNum() == 1)
    {
      row(*mat.Triang, r) = vect.getDense();
    }
    else if (vect.getNum() == 2)
    {
      row(*mat.Triang, r) = vect.getSparse();
    }
  }
  else if (num == 3)
  {
    if (vect.getNum() == 1)
    {
      row(*mat.Sym, r) = vect.getDense();
    }
    else if (vect.getNum() == 2)
    {
      row(*mat.Sym, r) = vect.getSparse();
    }
  }
  if (num == 4)
  {
    if (vect.getNum() == 1)
    {
      row(*mat.Sparse, r) = vect.getDense();
    }
    else if (vect.getNum() == 2)
    {
      row(*mat.Sparse, r) = vect.getSparse();
    }
  }
  if (num == 5)
  {
    if (vect.getNum() == 1)
    {
      row(*mat.Sparse, r) = vect.getDense();
    }
    else if (vect.getNum() == 2)
    {
      row(*mat.Banded, r) = vect.getSparse();
    }
  }
}

void MySimpleMatrix::getCol(unsigned int r, MySimpleVector &vect)const
{

  if (r >= size2() || r < 0)
    SiconosMatrixException::selfThrow("getCol : col is out of range");

  if (vect.size() != size1())
    SiconosMatrixException::selfThrow("getCol : inconsistent sizes");

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

void MySimpleMatrix::setCol(unsigned int r, const MySimpleVector &vect)
{

  if (r >= size2() || r < 0)
    SiconosMatrixException::selfThrow("setCol : col is out of range");

  if (vect.size() != size1())
    SiconosMatrixException::selfThrow("setCol : inconsistent sizes");

  if (num == 1)
  {
    if (vect.getNum() == 1)
    {
      column(*mat.Dense, r) = vect.getDense();
    }
    else if (vect.getNum() == 2)
    {
      column(*mat.Dense, r) = vect.getSparse();
    }
  }
  else if (num == 2)
  {
    if (vect.getNum() == 1)
    {
      column(*mat.Triang, r) = vect.getDense();
    }
    else if (vect.getNum() == 2)
    {
      column(*mat.Triang, r) = vect.getSparse();
    }
  }
  else if (num == 3)
  {
    if (vect.getNum() == 1)
    {
      column(*mat.Sym, r) = vect.getDense();
    }
    else if (vect.getNum() == 2)
    {
      column(*mat.Sym, r) = vect.getSparse();
    }
  }
  else if (num == 4)
  {
    if (vect.getNum() == 1)
    {
      column(*mat.Sparse, r) = vect.getDense();
    }
    else if (vect.getNum() == 2)
    {
      column(*mat.Sparse, r) = vect.getSparse();
    }
  }
  else if (num == 5)
  {
    if (vect.getNum() == 1)
    {
      column(*mat.Banded, r) = vect.getDense();
    }
    else if (vect.getNum() == 2)
    {
      column(*mat.Banded, r) = vect.getSparse();
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

  // \warning FP : Review transpose, with void as return and output as an input arg, to avoid unallocated matrices
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
  else if (m.getNum() == 4)
  {
    SparseMat p = trans(m.getSparse());
    return p;
  }
  else if (m.getNum() == 5) // Boost error
  {
    SiconosMatrixException::selfThrow("transpose of a banded matrix not allowed.");
    BandedMat p = trans(m.getBanded());
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
  unsigned int size1 = (*this).size1();
  unsigned int size2 = (*this).size2();
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
  unsigned int size1 = (*this).size1();
  unsigned int size2 = (*this).size2();
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

double& MySimpleMatrix::operator()(unsigned int row, unsigned int col)
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
    SiconosMatrixException::selfThrow("op() (unsigned int, unsigned int) : invalid type of matrix");
    break;
  }
  return d;
}

double MySimpleMatrix::operator()(unsigned int row, unsigned int col)const
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
    SiconosMatrixException::selfThrow("op() (unsigned int, unsigned int) : invalid type of matrix");
    break;
  }
  return d;
}

const MySimpleMatrix& MySimpleMatrix::operator = (const MySimpleMatrix& m)
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
      SiconosMatrixException::selfThrow("op= (const MySimpleMatrix) : invalid type of matrix");
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
      SiconosMatrixException::selfThrow("assignment of a bad type of matrix into a triangular one.");
      break;
    }
    break;
  case 3:
    if (m.getNum() == 3)
      *mat.Sym = m.getSym();
    else
      SiconosMatrixException::selfThrow("bad assignment of matrix (symetric one = dense or ...)");
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
      SiconosMatrixException::selfThrow("op= (const MySimpleMatrix) : invalid type of matrix");
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
      SiconosMatrixException::selfThrow("op= (const MySiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 2:
    switch (m.getNum())
    {
    case 1:
      *mat.Triang = m.getDense();
      break;
    case 2:
      *mat.Triang = m.getTriang();
      break;
    case 3:
      *mat.Triang = m.getSym();
      break;
    case 4:
      *mat.Triang = m.getSparse();
      break;
    case 5:
      *mat.Triang = m.getBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op= (const MySiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 3:
    switch (m.getNum())
    {
    case 1:
      *mat.Sym = m.getDense();
      break;
    case 2:
      *mat.Sym = m.getTriang();
      break;
    case 3:
      *mat.Sym = m.getSym();
      break;
    case 4:
      *mat.Sym = m.getSparse();
      break;
    case 5:
      *mat.Sym = m.getBanded();
      break;
    default:
      SiconosMatrixException::selfThrow("op= (const MySiconosMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 4:
    switch (m.getNum())
    {
    case 1:
      *mat.Sparse = m.getDense();
      break;
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
      SiconosMatrixException::selfThrow("op= (const MySimpleMatrix) : invalid type of matrix");
      break;
    }
    break;
  case 5:
    switch (m.getNum())
    {
    case 1:
      *mat.Banded = m.getDense();
      break;
    case 2:
      *mat.Banded = m.getTriang();
      break;
    case 3:
      *mat.Banded = m.getSym();
      break;
    case 4:
      *mat.Banded = m.getSparse();
      break;
    case 5:
      *mat.Banded = m.getBanded();
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
      SiconosMatrixException::selfThrow("op+= (const MySiconosMatrix) : invalid type of matrix");
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
      SiconosMatrixException::selfThrow("op+= (const MySiconosMatrix) : invalid type of matrix");
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
      SiconosMatrixException::selfThrow("op+= (const MySiconosMatrix) : invalid type of matrix");
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
      SiconosMatrixException::selfThrow("op= (const MySimpleMatrix) : invalid type of matrix");
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
      SiconosMatrixException::selfThrow("op-= (const MySiconosMatrix) : invalid type of matrix");
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
      SiconosMatrixException::selfThrow("op-= (const MySiconosMatrix) : invalid type of matrix");
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
      SiconosMatrixException::selfThrow("op-= (const MySiconosMatrix) : invalid type of matrix");
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
      SiconosMatrixException::selfThrow("op-= (const MySimpleMatrix) : invalid type of matrix");
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
      SiconosMatrixException::selfThrow("op-= (const MySimpleMatrix) : invalid type of matrix");
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
  double norm = (sub(m, x)).normInf();
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

  if (x.getNum() == m.getNum())
  {
    if (x.getNum() == 1)
    {
      p = x.getDense() + m.getDense();
      return p;
    }
    else if (x.getNum() == 2)
    {
      t = x.getTriang() + m.getTriang();
      return t;
    }
    else if (x.getNum() == 3)
    {
      s = x.getSym() + m.getSym();
      return s;
    }
    else if (x.getNum() == 4)
    {
      sp = x.getSparse() + m.getSparse();
      return sp;
    }
    else if (x.getNum() == 5)
    {
      b.resize(m.size1(), m.size2(), (m.getBanded()).lower(), (m.getBanded()).upper(), false);
      b = x.getBanded() + m.getBanded();
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

  if (x.getNum() == m.getNum())
  {
    if (x.getNum() == 1)
    {
      p = x.getDense() - m.getDense();
      return p;
    }
    else if (x.getNum() == 2)
    {
      t = x.getTriang() - m.getTriang();
      return t;
    }
    else if (x.getNum() == 3)
    {
      s = x.getSym() - m.getSym();
      return s;
    }
    else if (x.getNum() == 4)
    {
      sp = x.getSparse() - m.getSparse();
      return sp;
    }
    else if (x.getNum() == 5)
    {
      b.resize(m.size1(), m.size2(), (m.getBanded()).lower(), (m.getBanded()).upper(), false);
      b = x.getBanded() - m.getBanded();
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

  if (x.getNum() == m.getNum())
  {
    if (x.getNum() == 1)
    {
      p = prod(x.getDense(), m.getDense());
      return p;
    }
    else if (x.getNum() == 2)
    {
      t = prod(x.getTriang(), m.getTriang());
      return t;
    }
    else if (x.getNum() == 3)
    {
      s = prod(x.getSym(), m.getSym());
      return s;
    }
    else if (x.getNum() == 4)
    {
      sp = prod(x.getSparse(), m.getSparse());
      return sp;
    }
    else if (x.getNum() == 5)
    {
      SiconosMatrixException::selfThrow("Matrix product : banded*banded does not return a banded matrix, use prod instead of *");
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

MySimpleMatrix sub(const MySiconosMatrix &x, const MySiconosMatrix& m)
{
  DenseMat p;

  if ((x.size1() != m.size1()) || (x.size2() != m.size2()))
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

MySimpleMatrix prod(const MySiconosMatrix &x, const MySiconosMatrix& m)
{
  DenseMat p;
  TriangMat t;
  SymMat s;
  SparseMat sp;
  BandedMat b;

  if ((x.size2() != m.size1()))
    SiconosMatrixException::selfThrow("Matrix function prod : inconsistent sizes");

  if (m.getNum() == 1)
  {
    DenseMat q;
    q = m.getDense();
    if (x.getNum() == 1)
    {
      p = prod(x.getDense(), q);
    }
    else if (x.getNum() == 2)
    {
      p = prod(x.getTriang(), q);
    }
    else if (x.getNum() == 3)
    {
      p = prod(x.getSym(), q);
    }
    else if (x.getNum() == 4)
    {
      p = prod(x.getSparse(), q);
    }
    else if (x.getNum() == 5)
    {
      p = prod(x.getBanded(), q);
    }
    else
      SiconosMatrixException::selfThrow("Matrix function prod: invalid type of matrix");

    return p;
  }
  else if (m.getNum() == 2)
  {
    TriangMat q;
    q = m.getTriang();
    if (x.getNum() == 1)
    {
      p = prod(x.getDense(), q);
      return p;
    }
    else if (x.getNum() == 2)
    {
      t = prod(x.getTriang(), q);
      return t;
    }
    else if (x.getNum() == 3)
    {
      p = prod(x.getSym(), q);
      return p;
    }
    else if (x.getNum() == 4)
    {
      p = prod(x.getSparse(), q);
      return p;
    }
    else if (x.getNum() == 5)
    {
      p = prod(x.getBanded(), q);
      return p;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function prod: invalid type of matrix");

  }
  else if (m.getNum() == 3)
  {
    SymMat q;
    q = m.getSym();
    if (x.getNum() == 1)
    {
      p = prod(x.getDense(), q);
      return p;
    }
    else if (x.getNum() == 2)
    {
      p = prod(x.getTriang(), q);
      return p;
    }
    else if (x.getNum() == 3)
    {
      s = prod(x.getSym(), q);
      return s;
    }
    else if (x.getNum() == 4)
    {
      p = prod(x.getSparse(), q);
      return p;
    }
    else if (x.getNum() == 5)
    {
      p = prod(x.getBanded(), q);
      return p;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function prod: invalid type of matrix");

  }
  else if (m.getNum() == 4)
  {
    SparseMat q;
    q = m.getSparse();
    if (x.getNum() == 1)
    {
      p = prod(x.getDense(), q);
      return p;
    }
    else if (x.getNum() == 2)
    {
      p = prod(x.getTriang(), q);
      return p;
    }
    else if (x.getNum() == 3)
    {
      p = prod(x.getSym(), q);
      return p;
    }
    else if (x.getNum() == 4)
    {
      sp = prod(x.getSparse(), q);
      return sp;
    }
    else if (x.getNum() == 5)
    {
      p = prod(x.getBanded(), q);
      return p;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function prod: invalid type of matrix");

  }
  else if (m.getNum() == 5)
  {
    BandedMat q;
    q = m.getBanded();
    if (x.getNum() == 1)
    {
      p = prod(x.getDense(), q);
      return p;
    }
    else if (x.getNum() == 2)
    {
      p = prod(x.getTriang(), q);
      return p;
    }
    else if (x.getNum() == 3)
    {
      p = prod(x.getSym(), q);
      return p;
    }
    else if (x.getNum() == 4)
    {
      p = prod(x.getSparse(), q);
      return p;
    }
    else if (x.getNum() == 5)
    {
      b = prod(x.getBanded(), q);
      return b;
    }
    else
      SiconosMatrixException::selfThrow("Matrix function prod: invalid type of matrix");

  }
  else
  {
    SiconosMatrixException::selfThrow("Matrix function prod: invalid type of matrix");
  }

}

MySimpleMatrix multTranspose(const MySiconosMatrix &x, const MySiconosMatrix &m)
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
  else if (m.getNum() == 5)
  {
    b.resize(m.size1(), m.size2(), (m.getBanded()).lower(), (m.getBanded()).upper(), false);
    b = m.getBanded() * d;
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
  else if (m.getNum() == 5)
  {
    b.resize(m.size1(), m.size2(), (m.getBanded()).lower(), (m.getBanded()).upper(), false);
    b = m.getBanded() * d;
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
  if (m.getNum() == 1)
  {
    p = d * m.getDense();
    return p;
  }
  else if (m.getNum() == 2)
  {
    t = d * m.getTriang();
    return t;
  }
  else if (m.getNum() == 3)
  {
    s = d * m.getSym();
    return s;
  }
  else if (m.getNum() == 4)
  {
    sp = d * m.getSparse();
    return sp;
  }
  else if (m.getNum() == 5)
  {
    b.resize(m.size1(), m.size2(), (m.getBanded()).lower(), (m.getBanded()).upper(), false);
    b = d * m.getBanded();
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
  if (m.getNum() == 1)
  {
    p = d * m.getDense();
    return p;
  }
  else if (m.getNum() == 2)
  {
    t = d * m.getTriang();
    return t;
  }
  else if (m.getNum() == 3)
  {
    s = d * m.getSym();
    return s;
  }
  else if (m.getNum() == 4)
  {
    sp = d * m.getSparse();
    return sp;
  }
  else if (m.getNum() == 5)
  {
    b.resize(m.size1(), m.size2(), (m.getBanded()).lower(), (m.getBanded()).upper(), false);
    b = d * m.getBanded();
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
  if (m.getNum() == 1)
  {
    p = m.getDense() / d;
    return p;
  }
  else if (m.getNum() == 2)
  {
    t = m.getTriang() / d;
    return t;
  }
  else if (m.getNum() == 3)
  {
    s = m.getSym() / d;
    return s;
  }
  else if (m.getNum() == 4)
  {
    sp = m.getSparse() / d;
    return sp;
  }
  else if (m.getNum() == 5)
  {
    b.resize(m.size1(), m.size2(), (m.getBanded()).lower(), (m.getBanded()).upper(), false);
    b = m.getBanded() / d;
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
  if (m.getNum() == 1)
  {
    p = m.getDense() / d;
    return p;
  }
  else if (m.getNum() == 2)
  {
    t = m.getTriang() / d;
    return t;
  }
  else if (m.getNum() == 3)
  {
    s = m.getSym() / d;
    return s;
  }
  else if (m.getNum() == 4)
  {
    sp = m.getSparse() / d;
    return sp;
  }
  else if (m.getNum() == 5)
  {
    b.resize(m.size1(), m.size2(), (m.getBanded()).lower(), (m.getBanded()).upper(), false);
    b = m.getBanded() / d;
    return b;
  }
  else
    SiconosMatrixException::selfThrow("Matrix op / (const MySiconosMatrix&, int): invalid type of matrix");
}



