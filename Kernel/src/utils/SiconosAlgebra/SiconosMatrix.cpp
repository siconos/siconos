#include "SiconosMatrix.h"

// LAPACK F77 function not declared in lapack++.h
Extern  "C"
{
  void F77NAME(dgetri)(integer * m, doublereal * A, integer * lda, integer * ipiv,
  doublereal * WORK, integer * LWORK, integer * info);
}


/****************** constructor ******************/
SiconosMatrix::SiconosMatrix()
{
  this->mat.resize(0, 0);
  this->ipiv = NULL;
  this ->isPLUFactorized = false;
  this ->isPLUInversed = false;
}


/****************** constructor ******************/
SiconosMatrix::SiconosMatrix(int row, int col)
{
  this->mat.resize(row, col);
  this->ipiv = NULL;
  this ->isPLUFactorized = false;
  this ->isPLUInversed = false;
}

/****************** constructor ******************/
SiconosMatrix::SiconosMatrix(const LaGenMatDouble m)
{
  this->mat.resize(m.size(0), m.size(1));
  this->mat = m;
  this->ipiv = NULL;
  this->isPLUFactorized = false;
  this->isPLUInversed = false;

}

/****************** constructor ******************/
SiconosMatrix::SiconosMatrix(const LaVectorDouble v, int row, int col)
{
  if (v.size() != row * col)
    SiconosMatrixException::selfThrow("constructor(LaVectorDouble,int,int) : invalid vector size");

  this ->mat.resize(row, col);
  int index = 0;
  for (int i = 0; i < this->mat.size(0); i++)
  {
    for (int j = 0; j < this->mat.size(1); j++)
    {
      this->mat(i, j) = v(index);
      index ++;
    }
  }
  this->ipiv = NULL;
  this -> isPLUFactorized = false;
  this ->isPLUInversed = false;

}

/****************** constructor ******************/
SiconosMatrix::SiconosMatrix(string file, bool ascii)
{
  if (ascii) this->read(file, "ascii");
  else this->read(file, "binary");
  this->ipiv = NULL;
  this ->isPLUFactorized = false;
  this ->isPLUInversed = false;
}


/****************** destructor ******************/
SiconosMatrix::~SiconosMatrix(void)
{
  IN("SiconosMatrix::~SiconosMatrix \n");
  if (ipiv != NULL) delete ipiv;
  OUT("SiconosMatrix::~SiconosMatrix \n");

}

/****************** size(d) ******************/
int SiconosMatrix::size(int d) const
{
  if ((d != 0) && (d != 1))
    SiconosMatrixException::selfThrow("function size() : Index out of range");
  return this->mat.size(d);
}


/****************** isSquare ******************/
bool SiconosMatrix::isSquare()
{
  return (this->mat.size(0) == this->mat.size(1));
}

/****************** addRow ******************/
bool SiconosMatrix::addRow(int row, SiconosVector &v)
{
  bool res = true;
  if (v.size() != this->mat.size(1))
    res = false;

  else
    for (int i = 0; i < this->mat.size(1); i++)
      this->mat(row, i) = v(i);

  return res;
}

/****************** getMatrix ******************/
LaGenMatDouble SiconosMatrix::getLaGenMatDouble() const
{
  return this->mat;
}

LaGenMatDouble& SiconosMatrix::getLaGenMatDouble()
{
  return this->mat;
}


/****************** () ******************/
// subscript operator to get/set individual elements
double& SiconosMatrix::operator()(int row, int col)
{
  if ((row >= mat.size(0)) || (col >= mat.size(1)))
    SiconosMatrixException::selfThrow("operator() : Index out of range");

  return mat(row, col);
}



/****************** getRow ******************/
SimpleVector SiconosMatrix::getRow(const int index) const
{
  //  cout<<"SiconosMatrix::getRow(int index)"<<endl;

  //  SiconosVector* res = new SimpleVector();
  //  if (index < mat.size(0) )
  //  {
  //    vector<double> v(mat.size(1));
  //    for (int i=0; i < mat.size(1); i++)
  //    {
  //      v[i] =
  //    }
  //    res->setValues(v);
  //  }
  //  return *res;

  if (index >= this->mat.size(0))
    SiconosMatrixException::selfThrow("getRow : Index out of range");
  else
  {
    const int rowSize = this->mat.size(1);
    SimpleVector res(rowSize);
    for (int i = 0; i < rowSize; i++)
      res(i) = mat(index, i);

    return res;
  }
}



/****************** multTranspose ******************/
SiconosMatrix SiconosMatrix::multTranspose(SiconosMatrix &B)
{
  LaGenMatDouble matResult(mat.size(0), B.mat.size(0));

  if (this->size(1) != B.size(1))
    SiconosMatrixException::selfThrow("Incompatible matrix dimension. Operation multTranspose is impossible");


  Blas_Mat_Mat_Trans_Mult(mat, B.mat, matResult);

  SiconosMatrix result(matResult);

  return result;
}


/****************** block Matrix Copy ******************/
void SiconosMatrix::blockMatrixCopy(SiconosMatrix &blockMat, int xPos, int yPos)
{
  if (xPos > this->mat.size(0) || yPos > this->mat.size(1))
  {
    SiconosMatrixException::selfThrow("ERROR. SiconosMatrix::blockMatrixCopy : Cannot copy block matrix into specified matrix [block matrix to copy is too big]");
  }
  else if ((xPos + blockMat.size(0)) > this->mat.size(0) || (yPos + blockMat.size(1)) > this->mat.size(1))
  {
    SiconosMatrixException::selfThrow("ERROR. SiconosMatrix::blockMatrixCopy : Cannot copy block matrix into specified matrix [bad position for the copy of the block matrix]");
  }
  else
  {
    for (int i = 0; i < blockMat.size(0); i++)
      for (int j = 0; j < blockMat.size(1); j++)
        this->mat(i + xPos, j + yPos) = blockMat(i, j);
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


/*************************************************/
SiconosMatrix& SiconosMatrix::operator = (const SiconosMatrix& m)
{
  //this->verbose("Matrix = operator ");

  this->mat.resize(m.mat.size(0), m.mat.size(1));
  this->mat = m.mat;
  if (m.ipiv == NULL)
  {
    if (this->ipiv != NULL)
    {
      delete ipiv;
      this->ipiv == NULL;
    }
  }
  else
  {

    if (this->ipiv != NULL)
    {
      *(this->ipiv) = *(m.ipiv);
    }
    else
    {
      ipiv = new LaVectorLongInt(*(m.ipiv));

    }

  }
  this->isPLUFactorized = m.isPLUFactorized;
  this->isPLUInversed = m.isPLUInversed;
  return *this;
}



/*************************************************/
// logical equal-to operator
bool operator == (const SiconosMatrix& m1, const SiconosMatrix& m2)
{
  SiconosMatrix::verbose("WARNING : operator == and != not performed by Blas.");
  int m = m1.size(0);
  int n = m1.size(1);

  if ((m != m2.size(0)) || (n != m2.size(1)))
    return false;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
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
  if ((m1.mat.size(0) == m2.mat.size(0)) && (m1.mat.size(1) == m2.mat.size(1)))
    return (m1.mat + m2.mat);
  else
    SiconosMatrixException::selfThrow("Incompatible matrix dimension. Addition is impossible");
}


/*************************************************/
SiconosMatrix operator - (const SiconosMatrix& m1, const SiconosMatrix& m2)
{
  if ((m1.mat.size(0) == m2.mat.size(0)) && (m1.mat.size(1) == m2.mat.size(1)))
    return (m1.mat - m2.mat);
  else
    SiconosMatrixException::selfThrow("Incompatible matrix dimension. Subtraction is impossible");
}


/*************************************************/
SiconosMatrix operator * (const SiconosMatrix& mat, const double& d)
{
  SiconosMatrix::verbose("WARNING : SiconosMatrix * double is not compute by Blas.");
  int m = mat.size(0);
  int n = mat.size(1);
  SiconosMatrix matTmp(m, n);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      matTmp(i, j) = mat.mat(i, j) * d;

  return SiconosMatrix(matTmp);
}


/*************************************************/
SiconosMatrix operator / (const SiconosMatrix& m, const double d)
{
  if (d == 0.0)
    SiconosMatrixException::selfThrow("Operator '/' : try to divide by 0");
  else
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




/*************************************************/
bool SiconosMatrix::read(string fileName, string mode)
{
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
    return true;
  }

  if (mode == "ascii")
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
    return true;
  }
  else
    SiconosMatrixException::selfThrow("Incorrect mode for reading");
}



/*************************************************/
bool SiconosMatrix::write(string fileName, string mode)
{
  //  if( (this->size(0) == 0)||(this->size(1) == 0) ) SiconosMatrixException::selfThrow("write impossible - SiconosMatrix empty");

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

//void SiconosMatrix::setVerbose(bool set)
//{
//  isVerbose = set;
//}


void SiconosMatrix::verbose(string msg)
{
  if (printVerbose)
  {
    cout << msg << endl;
  }
}

void SiconosMatrix::display() const
{
  cout << "| size : " << this->size(0) << ", " << this->size(1) << endl;
  cout << "| isPLUInversed : " << this->isPLUInversed << endl;

  if ((this->size(0) <= MAXSIZEFORDISPLAY) || (this->size(1) <= MAXSIZEFORDISPLAY))
    cout << this->mat << endl;
  else
    cout << "Display SiconosMatrix : matrix too large" << endl;
}

void SiconosMatrix::zero()
{
  /** VERY SLOW !

  for (int i = 0; i < this->size(0); i++)
  {
  for (int j = 0; j < this->size(1); j++)
  {
  (*this)(i, j) = 0.0;
  }
  }
  *
  **/

  double* array = this->mat.addr();
  int size = this->size(0) * this->size(1);
  for (int i = 0; i < size; i++)
    array[i] = 0.0;
}

SiconosMatrix  BlockMatrixAssemble(vector<SiconosMatrix*> VM)
{
  IN("SiconosMatrix BlockMatrixAssemble\n");
  // compute the size of the result
  int i, size = 0;
  SiconosMatrix res;


  for (i = 0; i < VM.size(); i++)
  {
    size += VM[i]->size(0); // we assume that the matrices contained in VM are squared
  }
  res.mat.resize(size, size);
  res.zero();

  // assemble the blocks
  int start = 0;
  int sizeOfBlock = 0;

  for (int k = 0; k < VM.size(); k++)
  {
    sizeOfBlock = VM[k]->size(0);
    for (int i = start; i < start + sizeOfBlock; i++)
    {
      for (int j = start; j < start + sizeOfBlock; j++)
      {
        res(i, j) = (*VM[k])(i - start, j - start);
        //cout<<i<<" "<<j<<endl;
      }
    }
    start += sizeOfBlock;
  }
  OUT("SiconosMatrix BlockMatrixAssemble\n");
  return res;
}


/****************************** linearSolve *************************************************/

SiconosMatrix SiconosMatrix::linearSolve(SiconosMatrix &B)
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

  int M = this->size(0);
  long nbRow = M;
  long int lda = (this->mat).inc(0) * (this->mat).gdim(0);

  this->ipiv = new LaVectorLongInt(M);

  if (this->mat.inc(0) != 1 || this->mat.inc(1) != 1)
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUFactorizationInPlace : The Matrix is non-contiguous. ");

  if (this->size(0) != this->size(1))
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUFactorizationInPlace : Square matrix expected.\n");


  if (M <= 0)  // Test like in dgesv
    SiconosMatrixException::selfThrow("SiconosMatrix::PLUFactorizationInPlace :  Problem in Matrix Size");

  //   if (lda <= max(1,M))
  //    SiconosMatrixException::selfThrow("SiconosMatrix::PLUFactorizationInPlace :  Problem in Matrix Leading Size lda");

  F77NAME(dgetrf)(&nbRow,  &nbRow, &(this->mat(0, 0)), &lda, &((*(this->ipiv))(0)), &info);

  if (info != 0)
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUFactorizationInPlace : Internal error in LAPACK: DGETRF ");

  this->isPLUFactorized = true;

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

  if ((this->isPLUFactorized))
  {
  }
  else
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUInverseInPlace : This Matrix is not LU Factorized   with Partial pivoting");


  long int info;

  int M = this->size(0);
  long nbRow = M;
  long int lda = (this->mat).inc(0) * (this->mat).gdim(0);

  if (this->mat.inc(0) != 1 || this->mat.inc(1) != 1)
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUFactorizationInPlace : The Matrix is non-contiguous. ");

  if (this->size(0) != this->size(1))
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUFactorizationInPlace : Square matrix expected.\n");

  if (M <= 0)  // Test like in dgesv
    SiconosMatrixException::selfThrow("SiconosMatrix::PLUFactorizationInPlace :  Problem in Matrix Size");

  // First step : Query for the optimal size for the workspace work(0,0)
  long int lwork = -1;
  SiconosMatrix work(1, 1);

  F77NAME(dgetri)(&nbRow, &(this->mat(0, 0)), &lda, &((*(this->ipiv))(0)), &(work(0, 0)), &lwork, &info);

  // Second step :  allocation of the Workspace and computtuaion of the inverse.
  lwork = work(0, 0);
  work.mat.resize(lwork, lwork);

  F77NAME(dgetri)(&nbRow, &(this->mat(0, 0)), &lda, &((*(this->ipiv))(0)), &(work(0, 0)), &lwork, &info);

  //SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )

  if (info != 0)
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUInverseInPlace : Internal error in LAPACK: DGETRI ");

  this->isPLUInversed = true ;

  OUT(" SiconosMatrix::PLUInverseInPlace()\n");
}

SiconosMatrix  SiconosMatrix::PLUForwardBackward(SiconosMatrix &B)
{
  IN(" SiconosMatrix::PLUForwardBackward()\n");
  SiconosMatrix X(B);
  this->PLUForwardBackwardInPlace(X);
  return X;
  OUT(" SiconosMatrix::PLUForwardBackward()\n");
}

void  SiconosMatrix::PLUForwardBackwardInPlace(SiconosMatrix &B)
{
  IN(" SiconosMatrix::PLUForwardBackwardInPlace()\n");

  if ((this->isPLUFactorized))
  {
  }
  else
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUForwardBackwardInPlace : This Matrix is not LU Factorized   with Partial pivoting");


  if (this->mat.inc(0) != 1 || this->mat.inc(1) != 1)
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUForwardBackwardInPlace : This Matrix is non-contiguous. ");



  if (this->size(0) != this->size(1))
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUForwardBackwardInPlace :Square matrix expected.\n");

  if (this->size(1) != B.size(0))
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUForwardBackwardInPlace : This Matrix and X are non-comformant. ");


  long int info;
  int M = this->size(0);
  long  nbRow = M;
  //long int N = A.size(1);
  long int K = B.size(1);
  long int lda = this->mat.inc(0) * this->mat.gdim(0);
  long int ldx = B.mat.inc(0) * B.mat.gdim(0);


  char nt = 'N';

  if (M <= 0)  // Test like in dgesv
    SiconosMatrixException::selfThrow("SiconosMatrix::PLUForwardBackwardInPlace :  Problem in Matrix Size");
  //   else if (lda <= max(1,M))
  //    SiconosMatrixException::selfThrow("SiconosMatrix::PLUFactorizationInPlace :  Problem in Matrix Leading Size lda");


  F77NAME(dgetrs)(&nt,  &nbRow, &K,  &(this->mat(0, 0)), &lda, &((*(this->ipiv))(0)),  &B(0, 0), &ldx,  &info);


  if (info != 0)
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUForwardBackwardInPlace() : Internal error in LAPACK: DGETRS ");

  OUT(" SiconosMatrix::PLUForwardBackwardInPlace()\n");
}

SimpleVector SiconosMatrix::PLUForwardBackward(SiconosVector &B)
{
  IN(" SiconosMatrix::PLUForwardBackward()\n");
  SimpleVector X(B);
  this->PLUForwardBackwardInPlace(X);
  return X;
  OUT(" SiconosMatrix::PLUForwardBackward()\n");
}

void SiconosMatrix::PLUForwardBackwardInPlace(SiconosVector &B)
{
  IN(" SiconosMatrix::PLUForwardBackwardInPlace(SiconosVector B)\n");

  if ((this->isPLUFactorized))
  {
  }
  else
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUForwardBackwardInPlace : This Matrix is not LU Factorized with Partial pivoting");


  if (this->mat.inc(0) != 1 || this->mat.inc(1) != 1)
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUForwardBackwardInPlace : This Matrix is non-contiguous. ");

  if (this->size(0) != this->size(1))
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUForwardBackwardInPlace :Square matrix expected.\n");

  if (this->size(1) != B.size())
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUForwardBackwardInPlace : This Matrix and X are non-comformant. ");
  long int info;
  int M = this->size(0);
  long  nbRow = M;
  //long int N = A.size(1);
  long int K = 1;
  long int lda = this->mat.inc(0) * this->mat.gdim(0);

  //long int ldx = X.inc(0) * X.gdim(0);
  long int ldx = B.size();

  char nt = 'N';

  if (M <= 0)  // Test like in dgesv
    SiconosMatrixException::selfThrow("SiconosMatrix::PLUForwardBackwardInPlace :  Problem in Matrix Size");

  F77NAME(dgetrs)(&nt,  &nbRow, &K,  &(this->mat(0, 0)), &lda, &((*(this->ipiv))(0)),  &B(0), &ldx,  &info);


  if (info != 0)
    SiconosMatrixException::selfThrow(" SiconosMatrix::PLUForwardBackwardInPlace() : Internal error in LAPACK: DGETRS ");


  OUT(" SiconosMatrix::PLUForwardBackwardInPlace(SiconosVector B)\n");
}



