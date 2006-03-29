/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#include "SimpleMatrix.h"
using namespace std;


// --- CONSTRUCTORS ---

// Default (private)
SimpleMatrix::SimpleMatrix(): SiconosMatrix(false), ipiv(NULL), isPLUFactorized(false), isPLUInversed(false)
{}

// copy
SimpleMatrix::SimpleMatrix(const SiconosMatrix& m):
  SiconosMatrix(false), ipiv(NULL), isPLUFactorized(m.isFactorized()), isPLUInversed(m.isInversed())
{
  mat = m.getLaGenMatDouble();
}

SimpleMatrix::SimpleMatrix(const SimpleMatrix& m):
  SiconosMatrix(false), ipiv(NULL), isPLUFactorized(m.isFactorized()), isPLUInversed(m.isInversed())
{
  mat = m.getLaGenMatDouble();
}

// From dimensions
SimpleMatrix::SimpleMatrix(const unsigned int& row, const unsigned int& col):
  SiconosMatrix(false), ipiv(NULL), isPLUFactorized(false), isPLUInversed(false)
{
  mat.resize(row, col);
  zero();
}

// copy a LaGenMatDouble
SimpleMatrix::SimpleMatrix(const LaGenMatDouble& m):
  SiconosMatrix(false), ipiv(NULL), isPLUFactorized(false), isPLUInversed(false)
{
  mat = m;
}

// built from dimensions and a LaVectorDouble
SimpleMatrix::SimpleMatrix(const LaVectorDouble& v, const unsigned int& row, const unsigned int& col):
  SiconosMatrix(false), ipiv(NULL), isPLUFactorized(false), isPLUInversed(false)

{
  if ((unsigned int)v.size() != row * col)
    SiconosMatrixException::selfThrow("constructor(LaVectorDouble,unsigned int,unsigned int) : invalid vector size");

  mat.resize(row, col);
  unsigned int index = 0;
  for (unsigned int i = 0; i < row; i++)
  {
    for (unsigned int j = 0; j < col; j++)
    {
      mat(i, j) = v(index);
      index ++;
    }
  }
}

SimpleMatrix::SimpleMatrix(const string& file, const bool& ascii):
  SiconosMatrix(false), ipiv(NULL), isPLUFactorized(false), isPLUInversed(false)
{
  if (ascii) read(file, "ascii");
  else read(file, "binary");
}

// --- DESTRUCTOR ---

SimpleMatrix::~SimpleMatrix()
{
  if (ipiv != NULL) delete ipiv;
}

// --- FUNCTIONS TO GET INFO ABOUT THE MATRIX ---

unsigned int SimpleMatrix::size(const unsigned int& d) const
{
  if ((d != 0) && (d != 1))
    SiconosMatrixException::selfThrow("function size() : Index out of range");
  return (unsigned int)mat.size(d);
}

bool SimpleMatrix::isSquare() const
{
  return (mat.size(0) == mat.size(1));
}

// --- GETTERS/SETTERS ---

const LaGenMatDouble SimpleMatrix::getLaGenMatDouble(const unsigned int& row, const unsigned int& col) const
{
  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("SimpleMatrix getLaGenMatDouble(row,col), row or col not equal to 0.");
  return mat;
}

void SimpleMatrix::setValue(const LaGenMatDouble& newMat, const unsigned int& row, const unsigned int& col)
{
  if (row != 0 || col != 0)
    SiconosMatrixException::selfThrow("SimpleMatrix getLaGenMatDouble(row,col), row or col not equal to 0.");
  mat = newMat;
  if (ipiv != NULL) delete ipiv;
  isPLUFactorized = false;
  isPLUInversed = false;
}

void SimpleMatrix::setRow(const unsigned int& row, const SiconosVector &v)
{
  if (v.size() != size(1))
    SiconosMatrixException::selfThrow("SimpleMatrix setRow: Index out of range");

  if (row >= size(0))
    SiconosMatrixException::selfThrow("SimpleMatrix setRow: Index out of range");

  for (unsigned int i = 0; i < size(1); i++)
    mat(row, i) = v(i);
}


void SimpleMatrix::getRow(const unsigned int& index, const SimpleVector& vOut) const
{
  unsigned int rowSize = vOut.size();
  if (rowSize != size(1))
    SiconosMatrixException::selfThrow("getRow : inconsistent sizes");
  if (index >= size(0) || index < 0)
    SiconosMatrixException::selfThrow("getRow : Index out of range");

  for (unsigned int i = 0; i < rowSize; i++)
    vOut(i) = mat(index, i);
}

void SimpleMatrix::setCol(const unsigned int& col, const SiconosVector &v)
{
  if (v.size() != size(0))
    SiconosMatrixException::selfThrow("SiconosMatrix setCol: Index out of range");

  if (col >= size(1))
    SiconosMatrixException::selfThrow("SiconosMatrix setCol: Index out of range");

  for (unsigned int i = 0; i < size(0); i++)
    mat(i , col) = v(i);
}

void SimpleMatrix::getCol(const unsigned int& index, const SimpleVector& vOut) const
{
  unsigned int colSize = vOut.size();
  if (colSize != size(0))
    SiconosMatrixException::selfThrow("getCol : inconsistent sizes");

  if (index >= size(1) || index < 0)
    SiconosMatrixException::selfThrow("getCol : Index out of range");

  for (unsigned int i = 0; i < colSize; i++)
    vOut(i) = mat(i, index);
}

void SimpleMatrix::getBlock(const vector<unsigned int>& index_list, SiconosMatrix& block) const
{
  // index_list = {i,j,k,l}
  // get block between lines i-j and columns k-l

  unsigned int i = index_list[0];
  unsigned int j = index_list[1];
  unsigned int k = index_list[2];
  unsigned int l = index_list[3];

  if (i >= size(0) || j >= size(0) || k >= size(1) || l >= size(1))
    SiconosMatrixException::selfThrow("getBlock : Index out of range");
  if (i > j || k > l)
    SiconosMatrixException::selfThrow("getBlock : wrong index_list");

  // we use Lapack++ block matrix getter
  LaIndex * I = new LaIndex(i, j);
  LaIndex * J = new LaIndex(k, l);

  if (mat(*I, *J).size(0) != (int)block.size(0) || mat(*I, *J).size(1) != (int)block.size(1))
    SiconosMatrixException::selfThrow("getBlock : inconsistent sizes");
  block.setValue(mat(*I, *J));
  delete I;
  delete J;
}

void SimpleMatrix::getBlock(const vector<unsigned int>& indexRow, const vector<unsigned int>& indexCol, SiconosMatrix& block) const
{
  unsigned int row = block.size(0);
  unsigned int col = block.size(1);
  if (row != indexRow.size() || col != indexCol.size())
    SiconosMatrixException::selfThrow("getBlock : wrong indexes list");

  unsigned int k = 0, l = 0;
  unsigned int matRow = size(0), matCol = size(1);

  vector<unsigned int>::const_iterator itRow, itCol;

  for (itRow = indexRow.begin(); itRow != indexRow.end(); itRow++)
  {
    if (*itRow >= matRow) SiconosMatrixException::selfThrow("getBlock : row index out of range");
    for (itCol = indexCol.begin(); itCol != indexCol.end(); itCol++)
    {
      if (*itCol >= matCol) SiconosMatrixException::selfThrow("getBlock : column index out of range");
      block(k, l) = mat(*itRow, *itCol);
      l++;
    }
    l = 0;
    k++;
  }
}

// --- READ, WRITE ... ---

bool SimpleMatrix::read(const string& fileName, const string& mode)
{
  bool tmp = false;
  if (mode == "binary")
  {
    FILE * inFile = fopen(fileName.c_str(), "rb");    // open the input file in binary mode
    if (inFile == NULL)
      SiconosMatrixException::selfThrow("function read error : Fail to open \"" + fileName + "\"");

    int m, n;
    fread((char *) &m, sizeof(int), 1, inFile);   // read m
    fread((char *) &n, sizeof(int), 1, inFile);   // read n
    mat.resize(m, n);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
        fread((char*)&mat(i, j), sizeof(double), 1, inFile); // read a double

    fclose(inFile);
    tmp = true;
  }

  else if (mode == "ascii")
  {
    ifstream inFile(fileName.c_str(),  ifstream::in);

    if (inFile == NULL)
    {
      SiconosVectorException::selfThrow("function read error : Fail to open \"" + fileName + "\"");
    }

    int n, m;
    inFile >> m;
    inFile >> n;
    mat.resize(m, n);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
        inFile >> mat(i, j);
    inFile.close();
    tmp = true;
  }
  else
    SiconosMatrixException::selfThrow("Incorrect mode for reading");
  return tmp;
}

bool SimpleMatrix::write(const string& fileName, const string& mode) const
{
  //  if( (size(0) == 0)||(size(1) == 0) ) SiconosMatrixException::selfThrow("write impossible - SiconosMatrix empty");

  if ((mode != "binary") && (mode != "ascii"))
    SiconosMatrixException::selfThrow("Incorrect mode for writing");

  // open the file
  ofstream outFile(fileName.c_str());           // checks that it's opened
  if (!outFile.is_open())
    SiconosMatrixException::selfThrow("function write error : Fail to open \"" + fileName + "\"");

  int m = mat.size(0);
  int n = mat.size(1);
  if (mode == "binary")
  {
    outFile.write((char*)&m, sizeof(int));
    outFile.write((char*)&n, sizeof(int));
  }
  else if (mode == "ascii")
    outFile << m << ' ' << n;

  for (int i = 0; i < m; i++)
  {
    if (mode == "ascii")
    {
      outFile << endl;
    }
    for (int j = 0; j < n; j++)
    {
      if (mode == "binary")
      {
        outFile.write((char*)&mat(i, j), sizeof(double));
      }
      else if (mode == "ascii")
      {
        char buffer[30];
        sprintf(buffer, "%1.17e ", mat(i, j)); // /!\ depends on machine precision
        outFile << buffer;
      }
    }
  }
  outFile.close();
  return true;
}

bool SimpleMatrix::rawWrite(const string& fileName, const string& mode) const
{
  //  if( (size(0) == 0)||(size(1) == 0) ) SiconosMatrixException::selfThrow("write impossible - SiconosMatrix empty");

  if ((mode != "binary") && (mode != "ascii"))
    SiconosMatrixException::selfThrow("Incorrect mode for writing");

  // open the file
  ofstream outFile(fileName.c_str());           // checks that it's opened
  if (!outFile.is_open())
    SiconosMatrixException::selfThrow("function write error : Fail to open \"" + fileName + "\"");

  int m = mat.size(0);
  int n = mat.size(1);

  for (int i = 0; i < m; i++)
  {
    if (mode == "ascii")
    {
      outFile << endl;
    }
    for (int j = 0; j < n; j++)
    {
      if (mode == "binary")
      {
        outFile.write((char*)&mat(i, j), sizeof(double));
      }
      else if (mode == "ascii")
      {
        char buffer[30];
        sprintf(buffer, "%1.17e ", mat(i, j)); // /!\ depends on machine precision
        outFile << buffer;
      }
    }
  }
  outFile.close();
  return true;
}

void SimpleMatrix::zero()
{
  double* array = mat.addr();
  int sizeV = size(0) * size(1);
  for (int i = 0; i < sizeV; i++)
    array[i] = 0.0;
}

void SimpleMatrix::eye()
{
  unsigned int sizeMin;
  sizeMin = size(0);
  if (sizeMin > size(1)) sizeMin = size(1);
  zero();
  for (unsigned int j = 0; j < sizeMin; j++) mat(j, j) = 1.0;

}

void SimpleMatrix::display() const
{
  cout << "=======> SimpleMatrix, with " << size(0) << " lines and " << size(1) << " columns." << endl;
  cout << " --> isPLUInversed: " << isPLUInversed << endl;
  cout << " --> isPLUFactorized: " << isPLUFactorized << endl;
  if ((size(0) <= MAXSIZEFORDISPLAY) || (size(1) <= MAXSIZEFORDISPLAY))
    cout << mat << endl;
  else
    cout << "Sorry, can not display this SimpleMatrix, size > MaxSizeForDisplay." << endl;
  cout << "===== End of SimpleMatrix display ====" << endl;
}

// --- MATRICES HANDLING AND OPERATORS ---

SimpleMatrix SimpleMatrix::multTranspose(const SiconosMatrix & B)
{
  if (size(1) != B.size(1))
    SiconosMatrixException::selfThrow("Incompatible matrix dimension. Operation multTranspose is impossible");

  LaGenMatDouble matResult(mat.size(0), B.size(0));
  Blas_Mat_Mat_Trans_Mult(mat, B.getLaGenMatDouble(), matResult);

  return SimpleMatrix(matResult);
}

void SimpleMatrix::blockMatrixCopy(const SiconosMatrix &blockMat, const unsigned int& xPos, const unsigned int& yPos)
{
  if ((int)xPos > mat.size(0) || (int)yPos > mat.size(1))
    SiconosMatrixException::selfThrow("ERROR. SimpleMatrix::blockMatrixCopy : Cannot copy block matrix into specified matrix [block matrix to copy is too big]");
  else if ((int)(xPos + blockMat.size(0)) > mat.size(0) || (int)(yPos + blockMat.size(1)) > mat.size(1))
    SiconosMatrixException::selfThrow("ERROR. SimpleMatrix::blockMatrixCopy : Cannot copy block matrix into specified matrix [bad position for the copy of the block matrix]");
  else
  {
    for (unsigned int i = 0; i < blockMat.size(0); i++)
      for (unsigned int j = 0; j < blockMat.size(1); j++)
        mat(i + xPos, j + yPos) = blockMat(i, j);
  }
}

/****************** () ******************/
// subscript operator to get/set individual elements
double& SimpleMatrix::operator()(const int& row, const int& col)
{
  if ((row >= mat.size(0)) || (col >= mat.size(1)))
    SiconosMatrixException::selfThrow("operator() : Index out of range");

  return mat(row, col);
}

// subscript operator to get/set individual elements
double& SimpleMatrix::operator()(const unsigned int& row, const unsigned int& col)
{
  if (((int)row >= mat.size(0)) || ((int)col >= mat.size(1)))
    SiconosMatrixException::selfThrow("operator() : Index out of range");

  return mat(row, col);
}

// subscript operator to get/set individual elements
double& SimpleMatrix::operator()(const unsigned int& row, const unsigned int& col) const
{
  if (((int)row >= mat.size(0)) || ((int)col >= mat.size(1)))
    SiconosMatrixException::selfThrow("operator() : Index out of range");

  return mat(row, col);
}

/*************************************************/
SimpleMatrix& SimpleMatrix::operator = (const SiconosMatrix& m)
{
  if (m.isBlock())
    SiconosMatrixException::selfThrow("operator = : can not assign a block matrix to a simple, not yet implemented.");

  if (&m == this) return *this;

  if (m.size(0) != size(0) || m.size(1) != size(1))
    SiconosMatrixException::selfThrow("operator = : left and right value have inconsistent sizes.");

  mat = m.getLaGenMatDouble();
  if (ipiv != NULL) delete ipiv;
  ipiv = NULL;
  isPLUFactorized = m.isFactorized();
  isPLUInversed   = m.isInversed();
  return *this;
}

SimpleMatrix& SimpleMatrix::operator = (const SimpleMatrix& m)
{
  if (&m == this) return *this;
  if (m.size(0) != size(0) || m.size(1) != size(1))
    SiconosMatrixException::selfThrow("operator = : left and right value have inconsistent sizes.");
  mat = m.getLaGenMatDouble();
  if (ipiv != NULL) delete ipiv;
  ipiv = NULL;
  isPLUFactorized = m.isFactorized();
  isPLUInversed   = m.isInversed();
  return *this;
}

SimpleMatrix& SimpleMatrix::operator+=(const SiconosMatrix &m)
{
  unsigned int row = m.size(0);
  unsigned int col = m.size(1);
  if (row != size(0) || col != size(1))
    SiconosMatrixException::selfThrow("SimpleMatrix::operator+=  inconsistent sizes");

  for (unsigned int i = 0; i < row; i++)
    for (unsigned int j = 0; j < col; j++)
      mat(i, j) += (m.getLaGenMatDouble())(i, j);
  return *this;
}

SimpleMatrix& SimpleMatrix::operator-=(const SiconosMatrix &m)
{
  unsigned int row = m.size(0);
  unsigned int col = m.size(1);
  if (row != size(0) || col != size(1))
    SiconosMatrixException::selfThrow("SimpleMatrix::operator-=  inconsistent sizes");

  for (unsigned int i = 0; i < row; i++)
    for (unsigned int j = 0; j < col; j++)
      mat(i, j) -= (m.getLaGenMatDouble())(i, j);
  return *this;
}

SimpleMatrix& SimpleMatrix::operator*=(const double& d)
{
  mat *= d;
  return *this;
}

bool operator==(const SiconosMatrix& m1, const SiconosMatrix& m2)
{
  if (m1.isBlock() || m2.isBlock())
    SiconosMatrixException::selfThrow("SimpleMatrix comparison with a block matrix: not yet implemented.");
  double norm = (m1 - m2).normInf();
  return(norm < tolerance);
}

/*************************************************/
SimpleMatrix operator * (const SimpleMatrix& m1, const SimpleMatrix& m2)
{
  if (m1.size(1) != m2.size(0))
    SiconosMatrixException::selfThrow("SimpleMatrix product: Inconsistent matrices dimensions.");

  LaGenMatDouble matm3(m1.size(0), m2.size(1));
  Blas_Mat_Mat_Mult(m1.getLaGenMatDouble(), m2.getLaGenMatDouble(), matm3);

  return SimpleMatrix(matm3);
}

SimpleMatrix operator * (const SiconosMatrix& m1, const double& d)
{
  SimpleMatrix tmp(m1);
  tmp *= d;
  return tmp;
}

SimpleMatrix operator * (const double& d, const SiconosMatrix& m1)
{
  SimpleMatrix tmp(m1);
  tmp *= d;
  return tmp;
}

SimpleMatrix operator / (const SiconosMatrix& m1, const double& d)
{
  if (d == 0.0)
    SiconosMatrixException::selfThrow("SimpleMatrix/double : division by 0");
  SimpleMatrix tmp(m1);
  tmp *= (1.0 / d);
  return tmp;
}

SimpleMatrix operator + (const SiconosMatrix& m1, const SiconosMatrix& m2)
{
  if ((m1.size(0) != m2.size(0)) || (m1.size(1) != m2.size(1)))
    SiconosMatrixException::selfThrow("Matrix addition: inconsistent sizes");
  return (m1.getLaGenMatDouble() + m2.getLaGenMatDouble());
}

SimpleMatrix operator - (const SiconosMatrix& m1, const SiconosMatrix& m2)
{
  if ((m1.size(0) != m2.size(0)) || (m1.size(1) != m2.size(1)))
    SiconosMatrixException::selfThrow("Matrices substraction : inconsistent sizes");
  return (m1.getLaGenMatDouble() - m2.getLaGenMatDouble());
}

SimpleMatrix pow(const SimpleMatrix& m, const unsigned int& power)
{
  SimpleMatrix temp(m);
  unsigned int size = m.size(0);

  if (!m.isSquare())
    SiconosMatrixException::selfThrow("pow(SimpleMatrix), matrix is not square.");

  if (power < 0)
    SiconosMatrixException::selfThrow("pow(SimpleMatrix,n) with negative value is not supported");

  if (power > 0)
    for (unsigned int i = 1; i < power; i++)
      temp = temp * m.getLaGenMatDouble();

  else if (power == 0)
    for (unsigned int i = 0; i < size; i++)
      for (unsigned int j = 0; j < size; j++)
        temp.mat(i, j) = i == j ? 1 : 0;

  return temp;
}

const double SimpleMatrix::normInf() const
{
  return Blas_Norm_Inf(mat);
}

// --- COMPUTING WITH MATRICES  ---

void SimpleMatrix::linearSolve(const SiconosMatrix &B, SiconosMatrix &X)
{
  LaGenMatDouble tmp = X.getLaGenMatDouble();
  LaLinearSolve(mat, tmp, B.getLaGenMatDouble());
  X.setValue(tmp);
};

void  SimpleMatrix::PLUFactorizationInPlace()
{
  long int info;

  int M = size(0);
  long nbRow = M;
  long int lda = (mat).inc(0) * (mat).gdim(0);

  ipiv = new LaVectorLongInt(M);

  if (mat.inc(0) != 1 || mat.inc(1) != 1)
    SiconosMatrixException::selfThrow(" SimpleMatrix::PLUFactorizationInPlace : The Matrix is non-contiguous. ");

  if (size(0) != size(1))
    SiconosMatrixException::selfThrow(" SimpleMatrix::PLUFactorizationInPlace : Square matrix expected.\n");


  if (M <= 0)  // Test like in dgesv
    SiconosMatrixException::selfThrow("SimpleMatrix::PLUFactorizationInPlace :  Problem in Matrix Size");

  //   if (lda <= max(1,M))
  //    SiconosMatrixException::selfThrow("SimpleMatrix::PLUFactorizationInPlace :  Problem in Matrix Leading Size lda");

  F77NAME(dgetrf)(&nbRow,  &nbRow, &(mat(0, 0)), &lda, &((*(ipiv))(0)), &info);

  if (info != 0)
    SiconosMatrixException::selfThrow(" SimpleMatrix::PLUFactorizationInPlace : Internal error in LAPACK: DGETRF ");

  isPLUFactorized = true;
}

void  SimpleMatrix::PLUInverseInPlace()
{
  if ((isPLUFactorized))
  {
  }
  else
    SiconosMatrixException::selfThrow(" SimpleMatrix::PLUInverseInPlace : This Matrix is not LU Factorized   with Partial pivoting");


  long int info;

  int M = size(0);
  long nbRow = M;
  long int lda = (mat).inc(0) * (mat).gdim(0);

  if (mat.inc(0) != 1 || mat.inc(1) != 1)
    SiconosMatrixException::selfThrow(" SimpleMatrix::PLUFactorizationInPlace : The Matrix is non-contiguous. ");

  if (size(0) != size(1))
    SiconosMatrixException::selfThrow(" SimpleMatrix::PLUFactorizationInPlace : Square matrix expected.\n");

  if (M <= 0)  // Test like in dgesv
    SiconosMatrixException::selfThrow("SimpleMatrix::PLUFactorizationInPlace :  Problem in Matrix Size");

  // First step : Query for the optimal size for the workspace work(0,0)
  long int lwork = -1;
  SimpleMatrix work(1, 1);

  F77NAME(dgetri)(&nbRow, &(mat(0, 0)), &lda, &((*(ipiv))(0)), &(work(0, 0)), &lwork, &info);

  // Second step :  allocation of the Workspace and computtuaion of the inverse.
  lwork = static_cast<long int>(work(0, 0));
  work.mat.resize(lwork, lwork);

  F77NAME(dgetri)(&nbRow, &(mat(0, 0)), &lda, &((*(ipiv))(0)), &(work(0, 0)), &lwork, &info);

  //SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )

  if (info != 0)
    SiconosMatrixException::selfThrow(" SimpleMatrix::PLUInverseInPlace : Internal error in LAPACK: DGETRI ");

  isPLUInversed = true ;

}

void  SimpleMatrix::PLUForwardBackwardInPlace(SiconosMatrix &B)
{
  if (!isPLUFactorized)
    SiconosMatrixException::selfThrow(" SimpleMatrix::PLUForwardBackwardInPlace : This Matrix is not LU Factorized   with Partial pivoting");

  if (mat.inc(0) != 1 || mat.inc(1) != 1)
    SiconosMatrixException::selfThrow(" SimpleMatrix::PLUForwardBackwardInPlace : This Matrix is non-contiguous. ");

  if (size(0) != size(1))
    SiconosMatrixException::selfThrow(" SimpleMatrix::PLUForwardBackwardInPlace :Square matrix expected.\n");

  if (size(1) != B.size(0))
    SiconosMatrixException::selfThrow(" SimpleMatrix::PLUForwardBackwardInPlace : This Matrix and X are non-comformant. ");


  long int info;
  int M = size(0);
  long  nbRow = M;
  //long int N = A.size(1);
  long int K = B.size(1);
  long int lda = mat.inc(0) * mat.gdim(0);
  long int ldx = B.getLaGenMatDouble().inc(0) * B.getLaGenMatDouble().gdim(0);


  char nt = 'N';

  if (M <= 0)  // Test like in dgesv
    SiconosMatrixException::selfThrow("SimpleMatrix::PLUForwardBackwardInPlace :  Problem in Matrix Size");
  //   else if (lda <= max(1,M))
  //    SiconosMatrixException::selfThrow("SimpleMatrix::PLUFactorizationInPlace :  Problem in Matrix Leading Size lda");


  F77NAME(dgetrs)(&nt,  &nbRow, &K,  &(mat(0, 0)), &lda, &((*(ipiv))(0)),  &B(0, 0), &ldx,  &info);


  if (info != 0)
    SiconosMatrixException::selfThrow(" SimpleMatrix::PLUForwardBackwardInPlace() : Internal error in LAPACK: DGETRS ");
}

void SimpleMatrix::PLUForwardBackwardInPlace(SiconosVector &B)
{
  if (isFactorized())
  {
  }
  else
    SiconosMatrixException::selfThrow(" SimpleMatrix::PLUForwardBackwardInPlace : This Matrix is not LU Factorized with Partial pivoting");


  if (mat.inc(0) != 1 || mat.inc(1) != 1)
    SiconosMatrixException::selfThrow(" SimpleMatrix::PLUForwardBackwardInPlace : This Matrix is non-contiguous. ");

  if (size(0) != size(1))
    SiconosMatrixException::selfThrow(" SimpleMatrix::PLUForwardBackwardInPlace :Square matrix expected.\n");

  if (size(1) != B.size())
    SiconosMatrixException::selfThrow(" SimpleMatrix::PLUForwardBackwardInPlace : This Matrix and X are non-comformant. ");
  long int info;
  int M = size(0);
  long  nbRow = M;
  //long int N = A.size(1);
  long int K = 1;
  long int lda = mat.inc(0) * mat.gdim(0);

  //long int ldx = X.inc(0) * X.gdim(0);
  long int ldx = B.size();

  char nt = 'N';

  if (M <= 0)  // Test like in dgesv
    SiconosMatrixException::selfThrow("SimpleMatrix::PLUForwardBackwardInPlace :  Problem in Matrix Size");

  F77NAME(dgetrs)(&nt,  &nbRow, &K,  &(mat(0, 0)), &lda, &((*(ipiv))(0)),  &B(0), &ldx,  &info);


  if (info != 0)
    SiconosMatrixException::selfThrow(" SimpleMatrix::PLUForwardBackwardInPlace() : Internal error in LAPACK: DGETRS ");
}

SimpleMatrix BlockMatrixAssemble(const vector<SiconosMatrix*>& VM)
{
  // compute the size of the result
  unsigned int sizeRes = 0;

  vector<SiconosMatrix*>::const_iterator iter;
  for (iter = VM.begin(); iter != VM.end(); ++iter)
    sizeRes += (*iter)->size(0); // we assume that the matrices contained in VM are squared

  SimpleMatrix res(sizeRes, sizeRes);
  res.zero();

  // assemble the blocks
  unsigned int start = 0;
  unsigned int sizeOfBlock = 0;

  for (unsigned int k = 0; k < VM.size(); k++)
  {
    sizeOfBlock = VM[k]->size(0);
    for (unsigned int i = start; i < start + sizeOfBlock; i++)
    {
      for (unsigned int j = start; j < start + sizeOfBlock; j++)
        res(i, j) = (*VM[k])(i - start, j - start);
    }
    start += sizeOfBlock;
  }

  return res;
}



