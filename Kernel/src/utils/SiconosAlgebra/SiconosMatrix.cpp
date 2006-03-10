/* Siconos-Kernel version 1.1.3, Copyright INRIA 2005-2006.
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
#include "SiconosMatrix.h"
using namespace std;

// LAPACK F77 function not declared in lapack++.h
Extern  "C"
{
  void F77NAME(dgetri)(integer * m, doublereal * A, integer * lda, integer * ipiv,
  doublereal * WORK, integer * LWORK, integer * info);
}

// Private functions (why?)

void SiconosMatrix::verbose(const string& msg)
{
  if (printVerbose)
    cout << msg << endl;
}

// --- CONSTRUCTORS ---

// Default (private)
SiconosMatrix::SiconosMatrix(): ipiv(NULL), isPLUFactorized(false), isPLUInversed(false)
{}

// copy
SiconosMatrix::SiconosMatrix(const SiconosMatrix& m):
  ipiv(NULL), isPLUFactorized(m.isFactorized()), isPLUInversed(m.isInversed())
{
  mat = m.getLaGenMatDouble();
}

// From dimensions
SiconosMatrix::SiconosMatrix(const int& row, const int& col):
  ipiv(NULL), isPLUFactorized(false), isPLUInversed(false)
{
  mat.resize(row, col);
  zero();
}

// copy a LaGenMatDouble
SiconosMatrix::SiconosMatrix(const LaGenMatDouble& m):
  ipiv(NULL), isPLUFactorized(false), isPLUInversed(false)
{
  mat = m;
}

// built from dimensions and a LaVectorDouble
SiconosMatrix::SiconosMatrix(const LaVectorDouble& v, const int& row, const int& col):
  ipiv(NULL), isPLUFactorized(false), isPLUInversed(false)

{
  if (v.size() != row * col)
    SiconosMatrixException::selfThrow("constructor(LaVectorDouble,int,int) : invalid vector size");

  mat.resize(row, col);
  int index = 0;
  for (int i = 0; i < row; i++)
  {
    for (int j = 0; j < col; j++)
    {
      mat(i, j) = v(index);
      index ++;
    }
  }
}

SiconosMatrix::SiconosMatrix(const string& file, const bool& ascii):
  ipiv(NULL), isPLUFactorized(false), isPLUInversed(false)

{
  if (ascii) read(file, "ascii");
  else read(file, "binary");
}

// --- DESTRUCTOR ---

SiconosMatrix::~SiconosMatrix()
{
  IN("SiconosMatrix::~SiconosMatrix \n");
  if (ipiv != NULL) delete ipiv;
  OUT("SiconosMatrix::~SiconosMatrix \n");

}

// --- FUNCTIONS TO GET INFO ABOUT THE MATRIX ---

unsigned int SiconosMatrix::size(const unsigned int& d) const
{
  if ((d != 0) && (d != 1))
    SiconosMatrixException::selfThrow("function size() : Index out of range");
  return (unsigned int)mat.size(d);
}

bool SiconosMatrix::isSquare() const
{
  return (mat.size(0) == mat.size(1));
}

// --- GETTERS/SETTERS ---

// \Warning (FP): double def, why??
LaGenMatDouble SiconosMatrix::getLaGenMatDouble() const
{
  return mat;
}

//LaGenMatDouble& SiconosMatrix::getLaGenMatDouble() { return mat;}

void SiconosMatrix::setValue(const LaGenMatDouble& newMat)
{
  mat = newMat;
  if (ipiv != NULL) delete ipiv;
  isPLUFactorized = false;
  isPLUInversed = false;
}

void SiconosMatrix::setRow(const unsigned int& row, const SiconosVector &v)
{
  if ((int)v.size() != mat.size(1))
    SiconosMatrixException::selfThrow("SiconosMatrix setRow: Index out of range");

  if ((int)row >= mat.size(0))
    SiconosMatrixException::selfThrow("SiconosMatrix setRow: Index out of range");

  for (unsigned int i = 0; (int)i < mat.size(1); i++)
    mat(row, i) = v(i);
}


void SiconosMatrix::getRow(const int& index, const SimpleVector& vOut) const
{
  unsigned int rowSize = vOut.size();
  if ((int)rowSize != mat.size(1))
    SiconosMatrixException::selfThrow("getRow : inconsistent sizes");
  if (index >= mat.size(0) || index < 0)
    SiconosMatrixException::selfThrow("getRow : Index out of range");

  for (unsigned int i = 0; i < rowSize; i++)
    vOut(i) = mat(index, i);
}

void SiconosMatrix::getCol(const int& index, const SimpleVector& vOut) const
{
  unsigned int colSize = vOut.size();
  if ((int)colSize != mat.size(0))
    SiconosMatrixException::selfThrow("getCol : inconsistent sizes");

  if (index >= mat.size(1) || index < 0)
    SiconosMatrixException::selfThrow("getCol : Index out of range");

  for (unsigned int i = 0; i < colSize; i++)
    vOut(i) = mat(i, index);
}

void SiconosMatrix::getBlock(const vector<unsigned int>& index_list, SiconosMatrix& block) const
{
  // index_list = {i,j,k,l}
  // get block between lines i-j and columns k-l

  unsigned int i = index_list[0];
  unsigned int j = index_list[1];
  unsigned int k = index_list[2];
  unsigned int l = index_list[3];

  if ((int)i >= mat.size(0) || (int)j >= mat.size(0) || (int)k >= mat.size(1) || (int)l >= mat.size(1))
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

void SiconosMatrix::getBlock(const vector<unsigned int>& indexRow, const vector<unsigned int>& indexCol, SiconosMatrix& block) const
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

bool SiconosMatrix::read(const string& fileName, const string& mode)
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

bool SiconosMatrix::write(const string& fileName, const string& mode) const
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

bool SiconosMatrix::rawWrite(const string& fileName, const string& mode) const
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

void SiconosMatrix::display() const
{
  cout << "| size : " << size(0) << ", " << size(1) << endl;
  cout << "| isPLUInversed : " << isPLUInversed << endl;

  if ((size(0) <= MAXSIZEFORDISPLAY) || (size(1) <= MAXSIZEFORDISPLAY))
    cout << mat << endl;
  else
    cout << "Display SiconosMatrix : matrix too large" << endl;
}

void SiconosMatrix::zero()
{
  double* array = mat.addr();
  int sizeV = size(0) * size(1);
  for (int i = 0; i < sizeV; i++)
    array[i] = 0.0;
}

void SiconosMatrix::eye()
{
  unsigned int sizeMin;
  sizeMin = size(0);
  if (sizeMin > size(1)) sizeMin = size(1);
  zero();
  for (unsigned int j = 0; j < sizeMin; j++) mat(j, j) = 1.0;

}

// --- MATRICES HANDLING AND OPERATORS ---

SiconosMatrix SiconosMatrix::multTranspose(const SiconosMatrix & B)
{
  if (size(1) != B.size(1))
    SiconosMatrixException::selfThrow("Incompatible matrix dimension. Operation multTranspose is impossible");

  LaGenMatDouble matResult(mat.size(0), B.mat.size(0));
  Blas_Mat_Mat_Trans_Mult(mat, B.mat, matResult);

  //SiconosMatrix result(matResult);

  return SiconosMatrix(matResult);
}

void SiconosMatrix::blockMatrixCopy(const SiconosMatrix &blockMat, const unsigned int& xPos, const unsigned int& yPos)
{
  if ((int)xPos > mat.size(0) || (int)yPos > mat.size(1))
    SiconosMatrixException::selfThrow("ERROR. SiconosMatrix::blockMatrixCopy : Cannot copy block matrix into specified matrix [block matrix to copy is too big]");
  else if ((int)(xPos + blockMat.size(0)) > mat.size(0) || (int)(yPos + blockMat.size(1)) > mat.size(1))
    SiconosMatrixException::selfThrow("ERROR. SiconosMatrix::blockMatrixCopy : Cannot copy block matrix into specified matrix [bad position for the copy of the block matrix]");
  else
  {
    for (unsigned int i = 0; i < blockMat.size(0); i++)
      for (unsigned int j = 0; j < blockMat.size(1); j++)
        mat(i + xPos, j + yPos) = blockMat(i, j);
  }
}

/****************** output stream ******************/
ostream& operator << (ostream &ostrm, SiconosMatrix& m)
{
  if (m.mat.size(0) == 0)
  {
    cout << "Display SiconosMatrix : empty matrix" << endl;
  }
  //SiconosMatrixException::selfThrow("Try to display a 0-size matrix");

  else cout << m.getLaGenMatDouble() ;
  return ostrm;
}

/****************** input stream ******************/
istream& operator >> (istream &istrm, SiconosMatrix& m)
{
  LaGenMatDouble matTmp = m.getLaGenMatDouble();
  for (int i = 0; i < matTmp.size(0); i++)
    for (int j = 0; j < matTmp.size(1); j++)
    {
      cout << '[' << i + 1 << ',' << j + 1 << "] = ";
      cin >> m(i, j);
      if (cin.fail())
        cout << "INPUT ERROR";
    }
  return istrm;
}

/****************** () ******************/
// subscript operator to get/set individual elements
double& SiconosMatrix::operator()(const int& row, const int& col)
{
  if ((row >= mat.size(0)) || (col >= mat.size(1)))
    SiconosMatrixException::selfThrow("operator() : Index out of range");

  return mat(row, col);
}

// subscript operator to get/set individual elements
double& SiconosMatrix::operator()(const unsigned int& row, const unsigned int& col)
{
  if (((int)row >= mat.size(0)) || ((int)col >= mat.size(1)))
    SiconosMatrixException::selfThrow("operator() : Index out of range");

  return mat(row, col);
}

// subscript operator to get/set individual elements
double& SiconosMatrix::operator()(const unsigned int& row, const unsigned int& col) const
{
  if (((int)row >= mat.size(0)) || ((int)col >= mat.size(1)))
    SiconosMatrixException::selfThrow("operator() : Index out of range");

  return mat(row, col);
}

/*************************************************/
SiconosMatrix& SiconosMatrix::operator = (const SiconosMatrix& m)
{
  if (&m != this)
  {
    mat = m.getLaGenMatDouble();
    if (ipiv != NULL) delete ipiv;
    ipiv = NULL;
    isPLUFactorized = m.isFactorized();
    isPLUInversed   = m.isInversed();
  }
  return *this;
}

SiconosMatrix& SiconosMatrix::operator+=(const SiconosMatrix &m)
{
  IN(" SiconosMatrix::operator+= \n");
  int row = m.size(0);
  int col = m.size(1);
  if (row != mat.size(0) || col != mat.size(1))
    SiconosMatrixException::selfThrow("SiconosMatrix::operator+=  inconsistent sizes");

  for (int i = 0; i < row; i++)
    for (int j = 0; j < col; j++)
      mat(i, j) += (m.getLaGenMatDouble())(i, j);

  OUT(" SiconosMatrix::operator+=\n");
  return *this;
}




/*************************************************/
// logical equal-to operator
bool operator == (const SiconosMatrix& m1, const SiconosMatrix& m2)
{
  SiconosMatrix::verbose("WARNING : operator == and != not performed by Blas.");

  unsigned int m = m1.size(0);
  unsigned int n = m1.size(1);

  if ((m != m2.size(0)) || (n != m2.size(1)))
    return false;
  for (unsigned int i = 0; i < m; i++)
    for (unsigned int j = 0; j < n; j++)
      if ((double)m1.mat(i, j) != (double)m2.mat(i, j))
        return false;
  return true;
}


/*************************************************/
// logical no-equal-to operator
bool operator != (const SiconosMatrix& m1, const SiconosMatrix& m2)
{
  return !(m1 == m2);
}



/*************************************************/
SiconosMatrix operator * (const SiconosMatrix& m1, const SiconosMatrix& m2)
{
  if (m1.mat.size(1) != m2.mat.size(0))
    SiconosMatrixException::selfThrow("Incompatible matrix dimension. Multiplication is impossible");

  return (m1.mat * m2.mat);
}


/*************************************************/
SiconosMatrix operator + (const SiconosMatrix& m1, const SiconosMatrix& m2)
{
  if ((m1.mat.size(0) != m2.mat.size(0)) || (m1.mat.size(1) != m2.mat.size(1)))
    SiconosMatrixException::selfThrow("Matrix addition: inconsistent sizes");
  return (m1.mat + m2.mat);
}


/*************************************************/
SiconosMatrix operator - (const SiconosMatrix& m1, const SiconosMatrix& m2)
{
  if ((m1.mat.size(0) != m2.mat.size(0)) || (m1.mat.size(1) != m2.mat.size(1)))
    SiconosMatrixException::selfThrow("Matrices substraction : inconsistent sizes");
  return (m1.mat - m2.mat);
}


/*************************************************/
SiconosMatrix operator * (const SiconosMatrix& m1, const double& d)
{

  //SiconosMatrix matTmp(m1);

  int m = m1.size(0);
  int n = m1.size(1);
  SiconosMatrix matTmp(m, n);
  matTmp.mat = m1.mat;

  long int M = matTmp.size(0) * matTmp.size(1);
  long int incx = 1;
  double alpha = d;

  F77NAME(dscal)(&M , &alpha, &(matTmp.mat(0, 0)), &incx);

  return SiconosMatrix(matTmp);
}

/*************************************************/
SiconosMatrix operator / (const SiconosMatrix& m, const double d)
{
  if (d == 0.0)
    SiconosMatrixException::selfThrow("Operator '/' : division by 0");
  return (m * (1 / d));
}

/*************************************************/
SiconosMatrix operator ^ (const SiconosMatrix& m, const int pow)
{
  SiconosMatrix temp(m);
  int size = m.mat.size(0);

  if (size != m.mat.size(1))
    SiconosMatrixException::selfThrow("Incompatible matrix dimension. Operation ^ is impossible");

  SiconosMatrix::verbose("WARNING : operator ^ not performed by Blas.");
  if (pow < 0)
    SiconosMatrixException::selfThrow("Operator '^' with negative value is not supported");

  if (pow > 0)
    for (int i = 1; i < pow; i++)
      temp = temp * m.mat;

  else if (pow == 0)
    for (int i = 0; i < size; i++)
      for (int j = 0; j < size; j++)
        temp.mat(i, j) = i == j ? 1 : 0;

  return temp;
}

const double SiconosMatrix::normInf() const
{
  // get number or rows and columns
  unsigned int row = mat.size(0);
  unsigned int col = mat.size(1);
  double norm = 0, normMax = 0 ;
  for (unsigned int i = 0; i < row; ++i)
  {
    for (unsigned int j = 0; j < col; ++j)
      norm += fabs(mat(i, j));
    if (norm > normMax) normMax = norm;
    norm = 0;
  }

  return normMax;

}

// --- COMPUTING WITH MATRICES  ---

SiconosMatrix SiconosMatrix::linearSolve(const SiconosMatrix &B)
{
  SiconosMatrix X(B);
  LaLinearSolve(mat, X.mat, B.mat);
  return X;
};

// LU factorization with partial pivoting

SiconosMatrix  SiconosMatrix::PLUFactorization()
{
  IN(" SiconosMatrix::PLUFactorization()\n");
  SiconosMatrix Plu(*(this));
  Plu.PLUFactorizationInPlace();
  return Plu;
  OUT(" SiconosMatrix::PLUFactorization()\n");
}

void  SiconosMatrix::PLUFactorizationInPlace()
{
  IN(" SiconosMatrix::PLUFactorizationInPlace()\n");

  long int info;

  int M = size(0);
  long nbRow = M;
  long int lda = (mat).inc(0) * (mat).gdim(0);

  ipiv = new LaVectorLongInt(M);

  if (mat.inc(0) != 1 || mat.inc(1) != 1)
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUFactorizationInPlace : The Matrix is non-contiguous. ");

  if (size(0) != size(1))
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUFactorizationInPlace : Square matrix expected.\n");


  if (M <= 0)  // Test like in dgesv
    SiconosMatrixException::selfThrow("SiconosMatrix::PLUFactorizationInPlace :  Problem in Matrix Size");

  //   if (lda <= max(1,M))
  //    SiconosMatrixException::selfThrow("SiconosMatrix::PLUFactorizationInPlace :  Problem in Matrix Leading Size lda");

  F77NAME(dgetrf)(&nbRow,  &nbRow, &(mat(0, 0)), &lda, &((*(ipiv))(0)), &info);

  if (info != 0)
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUFactorizationInPlace : Internal error in LAPACK: DGETRF ");

  isPLUFactorized = true;

  OUT(" SiconosMatrix::PLUFactorizationInPlace()\n");
}

SiconosMatrix  SiconosMatrix::PLUInverse()
{
  IN(" SiconosMatrix::PLUInverse()\n");
  SiconosMatrix PluInv(*(this));
  PluInv.PLUInverseInPlace();
  OUT(" SiconosMatrix::PLUInverse()\n");
  return PluInv;
}

void  SiconosMatrix::PLUInverseInPlace()
{
  IN(" SiconosMatrix::PLUInverseInPlace()\n");

  if ((isPLUFactorized))
  {
  }
  else
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUInverseInPlace : This Matrix is not LU Factorized   with Partial pivoting");


  long int info;

  int M = size(0);
  long nbRow = M;
  long int lda = (mat).inc(0) * (mat).gdim(0);

  if (mat.inc(0) != 1 || mat.inc(1) != 1)
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUFactorizationInPlace : The Matrix is non-contiguous. ");

  if (size(0) != size(1))
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUFactorizationInPlace : Square matrix expected.\n");

  if (M <= 0)  // Test like in dgesv
    SiconosMatrixException::selfThrow("SiconosMatrix::PLUFactorizationInPlace :  Problem in Matrix Size");

  // First step : Query for the optimal size for the workspace work(0,0)
  long int lwork = -1;
  SiconosMatrix work(1, 1);

  F77NAME(dgetri)(&nbRow, &(mat(0, 0)), &lda, &((*(ipiv))(0)), &(work(0, 0)), &lwork, &info);

  // Second step :  allocation of the Workspace and computtuaion of the inverse.
  lwork = static_cast<long int>(work(0, 0));
  work.mat.resize(lwork, lwork);

  F77NAME(dgetri)(&nbRow, &(mat(0, 0)), &lda, &((*(ipiv))(0)), &(work(0, 0)), &lwork, &info);

  //SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )

  if (info != 0)
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUInverseInPlace : Internal error in LAPACK: DGETRI ");

  isPLUInversed = true ;

  OUT(" SiconosMatrix::PLUInverseInPlace()\n");
}

SiconosMatrix  SiconosMatrix::PLUForwardBackward(SiconosMatrix &B)
{
  IN(" SiconosMatrix::PLUForwardBackward()\n");
  SiconosMatrix X(B);
  PLUForwardBackwardInPlace(X);
  return X;
  OUT(" SiconosMatrix::PLUForwardBackward()\n");
}

void  SiconosMatrix::PLUForwardBackwardInPlace(SiconosMatrix &B)
{
  IN(" SiconosMatrix::PLUForwardBackwardInPlace()\n");

  if (!isPLUFactorized)
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUForwardBackwardInPlace : This Matrix is not LU Factorized   with Partial pivoting");

  if (mat.inc(0) != 1 || mat.inc(1) != 1)
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUForwardBackwardInPlace : This Matrix is non-contiguous. ");

  if (size(0) != size(1))
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUForwardBackwardInPlace :Square matrix expected.\n");

  if (size(1) != B.size(0))
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUForwardBackwardInPlace : This Matrix and X are non-comformant. ");


  long int info;
  int M = size(0);
  long  nbRow = M;
  //long int N = A.size(1);
  long int K = B.size(1);
  long int lda = mat.inc(0) * mat.gdim(0);
  long int ldx = B.mat.inc(0) * B.mat.gdim(0);


  char nt = 'N';

  if (M <= 0)  // Test like in dgesv
    SiconosMatrixException::selfThrow("SiconosMatrix::PLUForwardBackwardInPlace :  Problem in Matrix Size");
  //   else if (lda <= max(1,M))
  //    SiconosMatrixException::selfThrow("SiconosMatrix::PLUFactorizationInPlace :  Problem in Matrix Leading Size lda");


  F77NAME(dgetrs)(&nt,  &nbRow, &K,  &(mat(0, 0)), &lda, &((*(ipiv))(0)),  &B(0, 0), &ldx,  &info);


  if (info != 0)
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUForwardBackwardInPlace() : Internal error in LAPACK: DGETRS ");

  OUT(" SiconosMatrix::PLUForwardBackwardInPlace()\n");
}

SimpleVector SiconosMatrix::PLUForwardBackward(SiconosVector &B)
{
  IN(" SiconosMatrix::PLUForwardBackward()\n");
  SimpleVector X(B);
  PLUForwardBackwardInPlace(X);
  return X;
  OUT(" SiconosMatrix::PLUForwardBackward()\n");
}

void SiconosMatrix::PLUForwardBackwardInPlace(SiconosVector &B)
{
  IN(" SiconosMatrix::PLUForwardBackwardInPlace(SiconosVector B)\n");

  if ((isPLUFactorized))
  {
  }
  else
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUForwardBackwardInPlace : This Matrix is not LU Factorized with Partial pivoting");


  if (mat.inc(0) != 1 || mat.inc(1) != 1)
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUForwardBackwardInPlace : This Matrix is non-contiguous. ");

  if (size(0) != size(1))
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUForwardBackwardInPlace :Square matrix expected.\n");

  if (size(1) != B.size())
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUForwardBackwardInPlace : This Matrix and X are non-comformant. ");
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
    SiconosMatrixException::selfThrow("SiconosMatrix::PLUForwardBackwardInPlace :  Problem in Matrix Size");

  F77NAME(dgetrs)(&nt,  &nbRow, &K,  &(mat(0, 0)), &lda, &((*(ipiv))(0)),  &B(0), &ldx,  &info);


  if (info != 0)
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUForwardBackwardInPlace() : Internal error in LAPACK: DGETRS ");


  OUT(" SiconosMatrix::PLUForwardBackwardInPlace(SiconosVector B)\n");
}


SiconosMatrix BlockMatrixAssemble(const vector<SiconosMatrix*>& VM)
{
  IN("SiconosMatrix BlockMatrixAssemble\n");
  // compute the size of the result
  unsigned int sizeRes = 0;

  vector<SiconosMatrix*>::const_iterator iter;
  for (iter = VM.begin(); iter != VM.end(); ++iter)
    sizeRes += (*iter)->size(0); // we assume that the matrices contained in VM are squared

  SiconosMatrix res(sizeRes, sizeRes);
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

  return res

         OUT("SiconosMatrix BlockMatrixAssemble\n");
}



