#include "SimpleVector.h"
using namespace std;

// default
SimpleVector::SimpleVector(): SiconosVector()
{
  IN("SimpleVector() \n");
  composite = false;
  lavd = LaVectorDouble();
  OUT("SimpleVector() \n");
}

SimpleVector::SimpleVector(const string& file, const bool& ascii):
  SiconosVector()
{
  IN("SimpleVector(const string file, const bool ascii) \n");
  composite = false;
  string mode;
  if (ascii) mode = "ascii";
  else mode = "binary";
  read(file, mode);
  OUT("SimpleVector(const string file, const bool ascii) \n");
}

// copy from a std vector
SimpleVector::SimpleVector(const vector<double>& v):
  SiconosVector()
{
  composite = false;
  setValues(v);
}

// copy
SimpleVector::SimpleVector(const SimpleVector& v):
  SiconosVector()
{
  IN("SimpleVector(const SiconosVector& v) \n");
  composite = false;
  lavd = v.getValues();
  OUT("SimpleVector(const SiconosVector& v) \n");
}

// copy from SiconosVector
SimpleVector::SimpleVector(const SiconosVector& v):
  SiconosVector()
{
  if (!v.isComposite())
  {
    composite = false;
    lavd = v.getValues();
  }
  else
  {
    composite = false;
    unsigned int size = v.size();
    lavd.resize(size, 1);
    for (unsigned int i = 0; i < size; i++)
      lavd(i) = v.getValue(i);
  }
}

// with size
SimpleVector::SimpleVector(const int unsigned& size):
  SiconosVector()
{
  IN("SimpleVector (const int size) \n");
  composite = false;
  if (size < 0)
    SiconosVectorException::selfThrow(" SimpleVector:: constructor, negative size");
  lavd.resize(size, 1);
  OUT("SimpleVector (const int size) \n");
}

SimpleVector::~SimpleVector()
{
  IN("~SimpleVector() \n");
  OUT("~SimpleVector() \n");
}

/***********************************************************************************************/

void SimpleVector::display() const
{
  IN("SimpleVector::display() \n");
  cout << "=== Simple vector display === " << "of size " << size() << endl;
  if (size() <= M_MAXSIZEFORDISPLAY)
  {
    for (unsigned int i = 0; i < size(); i++)
      cout << lavd(i) << " ";
  }
  else cout << "Display SiconosVector : vector too large" << endl;
  cout << endl << " ========== " << endl;
  OUT("SimpleVector::display() \n");
}

double& SimpleVector::operator()(const int unsigned& index) const
{
  IN("SimpleVector::operator() \n");

  if ((int)index >= lavd.size())
    SiconosVectorException::selfThrow(" SimpleVector::operator()   -- index greater than vector size");

  OUT("SimpleVector::operator()\n");

  return const_cast<double&>(lavd(index));
}

vector<SiconosVector*> SimpleVector::getSvref() const
{
  SiconosVectorException::selfThrow(" SimpleVector::getSvref()  -- do not use this function for simple vector");
  // return to avoid warning ...
  vector<SiconosVector*> tmp;
  tmp.push_back(NULL);
  return tmp;
}

std::vector<int> SimpleVector::getTabIndex() const
{
  vector<int> tmp;
  tmp.push_back(size());
  return tmp;
}


void SimpleVector::zero()
{
  double *array = lavd.addr();
  const int sizeV = size();
  for (int i = 0; i < sizeV; i++)
    array[i] = 0.0;
}

string SimpleVector::toString() const
{
  char element[100];
  unsigned int i = 0;
  string vectorContent = "";
  while (i < size())
  {
    strcpy(element, "");
    sprintf(element, N_DOUBLE_PRECISION, (*this)(i));
    if (i > 0)
    {
      vectorContent += " ";
      vectorContent += element;
    }
    else vectorContent = element;
    i++;
  }
  return vectorContent;
}

void SimpleVector::setValue(const int unsigned& index, const double& newVal)
{
  if (index < 0 || (int)index > (lavd.size() - 1))
    SiconosVectorException::selfThrow("SimpleVector::setValue - wrong index value");
  lavd(index) = newVal;
}

void SimpleVector::setValues(const vector<double>& v, const int unsigned& index)
{
  if (v.size() < 0)
    SiconosVectorException::selfThrow("SimpleVector::setValues - negative vector size");

  lavd.resize(v.size(), 1);
  for (unsigned int i = 0; i < v.size(); i++)
    lavd(i) = v[i];
}

const double SimpleVector::getValue(const int unsigned& index) const
{
  if (index < 0 || (int)index > (lavd.size() - 1))
    SiconosVectorException::selfThrow("SimpleVector::setValue - wrong index value");
  return lavd(index);
}

void SimpleVector::getBlock(const vector<unsigned int>& index, SimpleVector& block) const
{
  unsigned int sizeBlock = block.size();
  if (sizeBlock != index.size())
    SiconosVectorException::selfThrow("getBlock : wrong indexes list");

  unsigned int k = 0;

  vector<unsigned int>::const_iterator it;

  for (it = index.begin(); it != index.end(); it++)
  {
    if (*it >= size()) SiconosVectorException::selfThrow("getBlock : index out of range");
    block(k) = lavd(*it);
    k++;
  }
}


unsigned int SimpleVector::size(const unsigned int& i) const
{
  if (i != 0)  SiconosVectorException::selfThrow("SimpleVector::size(i) - i!=0 -> not available");
  return lavd.size();
}

bool SimpleVector::read(const string& fileName, const string& mode)
{
  IN(" SimpleVector::read(string fileName, string mode) \n");

  bool res = false;
  int size, nbVectors;

  // Binary file
  if (mode == "binary")
  {
    FILE * inFile = fopen(fileName.c_str(), "rb");    // open the input file in binary mode
    if (inFile == NULL)
      SiconosVectorException::selfThrow(" SimpleVector::read : Fail to open file \"" + fileName + "\"");
    fread((char *) &nbVectors, sizeof(int), 1, inFile);   // read nbVectors
    if (nbVectors <= 0)
      SiconosVectorException::selfThrow(" SimpleVector::read : negative number of vectors!");
    if (nbVectors != 1)
      cout << "SimpleVector:: read; Warning: only the first vector will be read in the file" << endl;
    fread((char *) &size, sizeof(int), 1, inFile);   // read size
    if (size < 0)
      SiconosVectorException::selfThrow(" SimpleVector::read : try to read a vector with a negative size");
    lavd.resize(size, 1);
    for (int i = 0; i < size; i++)
      fread((char*) & (lavd)(i), sizeof(double), 1, inFile); // read a double

    fclose(inFile);
    res = true;
  }
  // Ascii file
  else if (mode == "ascii")
  {
    ifstream inFile(fileName.c_str(),  ifstream::in);
    if (inFile == NULL)
      SiconosVectorException::selfThrow(" SimpleVector::read : Fail to open file \"" + fileName + "\"");
    double tmp;
    inFile >> nbVectors;
    if (nbVectors <= 0)
      SiconosVectorException::selfThrow(" SimpleVector::read : negative number of vectors!");
    if (nbVectors != 1)
      cout << "SimpleVector:: read; Warning: only the first vector will be read in the file" << endl;
    // get values
    vector<double> vect;
    inFile >> size; // read size
    if (size < 0)
      SiconosVectorException::selfThrow(" SimpleVector::read : try to read a vector with a negative size");
    LaVectorDouble lvd(size);
    for (int i = 0; i < size; i++)
    {
      inFile >> tmp;
      vect.push_back(tmp);
    }

    setValues(vect);
    inFile.close();
    res = true;
  }
  OUT("SimpleVector::read(string fileName, string mode) \n");
  return res;
}


bool SimpleVector::write(const string& fileName, const string& mode) const
{
  IN("SimpleVector::write(string fileName, string mode) \n");

  bool res = false;
  if ((mode != "binary") && (mode != "ascii"))
    SiconosVectorException::selfThrow("SimpleVector::write : unknown mode");

  // open the file
  ofstream outFile(fileName.c_str());

  if (!outFile.is_open())
    SiconosVectorException::selfThrow("SimpleVector::write : : Fail to open file \"" + fileName + "\"");

  unsigned int nbVectors = 1;

  if (mode == "binary")
  {
    outFile.write((char*)&nbVectors, sizeof(int));
    unsigned int sizeV = size();
    outFile.write((char*)&sizeV, sizeof(int));
    for (unsigned int i = 0; i < sizeV; i++)
      outFile.write((char*) & (lavd)(i), sizeof(double));
    res = true;
  }
  else if (mode == "ascii")
  {
    outFile << nbVectors << endl;;

    int sizeV = size();
    outFile << sizeV << endl;
    //for (int i = 0; i < sizeV; i++)
    //{
    /* WARNING this buffer is dangerous, the size is machine dependant*/
    // char buffer[30];
    // sprintf(buffer,"%1.17e ",(*this)(i));
    // outFile << buffer;
    outFile << lavd;
    cout << endl;
    //}
    res = true;
  }
  outFile.close();

  OUT("SimpleVector::write(string fileName, string mode) \n");
  return res;
}

double* SimpleVector::getArray() const
{
  return lavd.addr();
}

double SimpleVector::norm() const
{
  return  Blas_Norm2(lavd);
}

/*******************************************************************************
 *          GENERIC INTERNAL OPERATORS                                 *
 *******************************************************************************/

SimpleVector &SimpleVector::operator+=(const SiconosVector &v)
{
  IN(" SimpleVector::operator+=(const SiconosVector &) \n");
  unsigned int sizeV = size();
  if (sizeV != v.size())
    SiconosVectorException::selfThrow(" SimpleVector::operator+=   -- the vectors have not the same size");

  for (unsigned int i = 0; i < sizeV; i++)
    (*this)(i) += v(i);
  OUT(" SimpleVector::operator-=(const SiconosVector &) \n");
  return *this;
}


SimpleVector &SimpleVector::operator-=(const SiconosVector &v)
{
  IN(" SimpleVector::operator-=(const SiconosVector &) \n");
  unsigned int sizeV = size();
  if (sizeV != v.size())
    SiconosVectorException::selfThrow(" SimpleVector::operator-=   -- the vectors have not the same size");

  for (unsigned int i = 0; i < sizeV; i++)
    (*this)(i) -= v(i);
  OUT(" SimpleVector::operator-=(const SiconosVector &) \n");
  return *this;
}


SimpleVector& SimpleVector::operator = (const SiconosVector& v)
{
  IN("SimpleVector::operator = GENERIC\n");
  if (this != &v)
  {
    unsigned int sizeV = size();
    if (sizeV != v.size())
      SiconosVectorException::selfThrow(" SimpleVector::operator = GENERIC  -- inconsistent sizes");

    else
    {
      for (unsigned int i = 0; i < sizeV; i++)
        (*this)(i) = v(i);
    }
  }
  OUT("SimpleVector::operator = GENERIC\n");
  return *this;
}

bool SimpleVector::operator == (const SiconosVector& v) const
{
  IN(" SimpleVector::operator == \n");
  bool res = true;

  if (v.isComposite()) res = false;
  else
  {
    if (this != &v)
    {
      unsigned int sizeV = size();
      if (sizeV == v.size())
      {
        if (sizeV == 0)
          res = true;
        else for (unsigned int i = 0; i < sizeV; i++)
          {
            if ((*this)(i) != v(i))
            {
              res = false;
              break;
            }
          }
      }
      else res =  false;
    }
  }
  OUT(" SimpleVector::operator == \n");
  return res;
}

bool SimpleVector::operator == (const SimpleVector& v) const
{
  IN(" SimpleVector::operator == \n");
  bool res = true;

  if (this != &v)
  {
    unsigned int sizeV = size();
    if (sizeV == v.size())
    {
      if (sizeV == 0)
        res = true;
      else for (unsigned int i = 0; i < sizeV; i++)
        {
          if ((*this)(i) != v(i))
          {
            res = false;
            break;
          }
        }
    }
    else res =  false;
  }
  OUT(" SimpleVector::operator == \n");
  return res;
}

bool SimpleVector::operator != (const SiconosVector& v) const
{
  return !(*this == v);
}

bool SimpleVector::operator != (const SimpleVector& v) const
{
  return !(*this == v);
}

/*******************************************************************************
 *          SPECIFIC INTERNAL OPERATORS                                *
 *******************************************************************************/
SimpleVector& SimpleVector::operator = (const SimpleVector& v)
{
  IN("SimpleVector::operator = SPECIFIC\n");
  if (v.size() != size())
    SiconosVectorException::selfThrow("SimpleVector::operator= inconsistent sizes");

  if (this != &v)
    lavd.copy(v.lavd);

  OUT("SimpleVector::operator = SPECIFIC\n");
  return *this;
}


SimpleVector &SimpleVector::operator+=(const SimpleVector &v)
{
  IN(" SimpleVector::operator+=(const SimpleVector &) \n");

  if (size() != v.size())
    SiconosVectorException::selfThrow(" SimpleVector::operator+=   -- the vectors have not the same size");

  Blas_Add_Mult(lavd, 1.0, v.lavd);

  OUT(" SimpleVector::operator+=(const SimpleVector &) \n");
  return *this;
}


SimpleVector &SimpleVector::operator-=(const SimpleVector &v)
{
  IN(" SimpleVector::operator-=(const SimpleVector &) \n");
  if (size() != v.size())
    SiconosVectorException::selfThrow(" SimpleVector::operator-=   -- the vectors have not the same size");

  Blas_Add_Mult(lavd, -1.0, v.lavd);

  OUT(" SimpleVector::operator-=(const SimpleVector &) \n");
  return *this;
}


SimpleVector &SimpleVector::operator*=(const double& d)
{
  IN(" SimpleVector::operator*=(const double d) \n");
  lavd = lavd * d;
  OUT(" SimpleVector::operator*=(const double d) \n");
  return *this;
}


SimpleVector &SimpleVector::operator/=(const double& d)
{
  IN(" SimpleVector::operator/=(const double d) \n");
  if (d == 0)
    SiconosVectorException::selfThrow(" SimpleVector::operator/ =   --division by 0");
  lavd = lavd * (1.0 / d);
  OUT(" SimpleVector::operator/=(const double d) \n");
  return *this;
}

//==================================================================================
//      GENERIC EXTERNAL OPERATORS
//==================================================================================


SimpleVector operator + (const SiconosVector& v1, const SiconosVector& v2)
{
  IN(" SimpleVector operator +  (const SimpleVector& , const SiconosVector&) GENERIC \n");

  if (v1.size() != v2.size())
    SiconosVectorException::selfThrow(" SimpleVector operator + --- the vectors have not the same size");

  SimpleVector sv(v1);
  sv += v2;

  OUT(" SimpleVector operator +  (const SimpleVector& , const SiconosVector&) GENERIC \n");
  return sv;
}


SimpleVector operator - (const SiconosVector& v1, const SiconosVector& v2)
{
  IN(" SimpleVector operator -  (const SimpleVector& , const SiconosVector&) GENERIC \n");

  if (v1.size() != v2.size())
    SiconosVectorException::selfThrow(" SimpleVector operator - --- the vectors have not the same size");

  SimpleVector sv(v1);
  sv -= v2;
  OUT(" SimpleVector operator -  (const SimpleVector& , const SiconosVector&) GENERIC \n");
  return sv;
}

//==================================================================================
//      GENERIC INTERNAL FUNCTIONS FOR MIXED OPERATIONS
//==================================================================================

SimpleVector SimpleVector::addition(const SiconosVector& v) const
{
  IN(" SimpleVector::addition(const SiconosVector& v) GENERIC \n");

  if (size() != v.size())
    SiconosVectorException::selfThrow(" SimpleVector::addition(const SiconosVector& v) --- the vectors have not the same size");

  SimpleVector sv(*this);
  sv += v;
  OUT(" SimpleVector::addition(const SiconosVector& v) GENERIC \n");
  return sv;
}


SimpleVector SimpleVector::subtraction(const SiconosVector& v) const
{
  IN(" SimpleVector::subtraction(const SiconosVector& v) GENERIC \n");

  if (size() != v.size())
    SiconosVectorException::selfThrow(" SimpleVector::subtraction(const SiconosVector& v) --- the vectors have not the same size");

  SimpleVector sv(*this);
  sv -= v;

  OUT(" SimpleVector::subtraction(const SiconosVector& v) GENERIC \n");
  return sv;
}


/*******************************************************************************
 *          SPECIFIC EXTERNAL OPERATORS                                *
 *******************************************************************************/
SimpleVector operator * (const SimpleVector& v, const double& d)
{
  IN(" SimpleVector operator * (const SimpleVector& v, const double d)  \n");

  SimpleVector sv(v);
  Blas_Scale(d, sv.lavd);
  OUT(" SimpleVector operator * (const SimpleVector& v, const double d)  \n");
  return sv;
}


SimpleVector operator * (const double& d, const SimpleVector& v)
{
  IN(" SimpleVector::operator *  (const SimpleVector& v, const double d) \n");

  SimpleVector sv(v);
  Blas_Scale(d, sv.lavd);

  OUT(" SimpleVector::operator * (const SimpleVector& v, const double d)  \n");
  return sv;
}


SimpleVector operator / (const SimpleVector& v, const double& d)
{
  IN(" SimpleVector operator /  (const SimpleVector& v, const double d) \n");

  if (d == 0.0)
    SiconosVectorException::selfThrow(" SimpleVector operator / --- division by 0");

  SimpleVector sv(v);
  Blas_Scale(1 / d, sv.lavd);

  OUT(" SimpleVector operator / (const SimpleVector& v, const double d)  \n");
  return sv;
}

/*
SimpleVector operator + (const SimpleVector& v1, const SimpleVector& v2)
{
  IN(" SimpleVector operator +  (const SimpleVector& , const SimpleVector&) \n");

  if (v1.size() != v2.size())
    SiconosVectorException::selfThrow(" SimpleVector operator + --- vectors have not the same size");

  SimpleVector sv(v1);
  Blas_Add_Mult(sv.lavd, 1.0, v2.lavd);

  return sv;

  OUT("SimpleVector operator +  (const SimpleVector& , const SimpleVector&)  \n");
}


SimpleVector operator - (const SimpleVector& v1, const SimpleVector& v2)
{
  IN(" SimpleVector operator -  (const SimpleVector& , const SimpleVector&) \n");

  if (v1.size() != v2.size())
    SiconosVectorException::selfThrow(" SimpleVector operator - --- the 2 vectors have not the same size");

  SimpleVector sv(v1);
  Blas_Add_Mult(sv.lavd, -1.0, v2.lavd);

  OUT("SimpleVector operator -  (const SimpleVector& , const SimpleVector&)  \n");
  return sv;
}
*/

SimpleVector operator * (const SiconosMatrix &m, const SimpleVector &v)
{
  SimpleVector sv(m.size(0));
  Blas_Mat_Vec_Mult(m.getLaGenMatDouble(), v.lavd, sv.lavd, 1.0, 0.0);
  return sv;
}

SimpleVector operator * (const SiconosMatrix &m, const SiconosVector &v)
{
  SimpleVector sv(m.size(0));
  // "transform" v (which is simple or composite) into a simple
  // not good -> to be optimized
  SimpleVector tmp(v);
  Blas_Mat_Vec_Mult(m.getLaGenMatDouble(), tmp.lavd, sv.lavd, 1.0, 0.0);
  return sv;
}

SimpleVector operator * (const SimpleVector &v, const SiconosMatrix &m)
{
  SimpleVector sv(m.size(1));
  Blas_Mat_Trans_Vec_Mult(m.getLaGenMatDouble(), v.lavd, sv.lavd, 1.0, 0.0);
  return sv;
}

SimpleVector operator * (const SiconosVector &v, const SiconosMatrix &m)
{
  SimpleVector sv(m.size(1));
  // "transform" v (which is simple or composite) into a simple
  // not good -> to be optimized
  SimpleVector tmp(v);
  Blas_Mat_Trans_Vec_Mult(m.getLaGenMatDouble(), tmp.lavd, sv.lavd, 1.0, 0.0);
  return sv;
}

SimpleVector matTransVecMult(SiconosMatrix &m, SimpleVector &v)
{
  SimpleVector sv(m.size(1));
  Blas_Mat_Trans_Vec_Mult(m.getLaGenMatDouble(), v.lavd, sv.lavd, 1.0, 0.0);
  return sv;
}

SimpleVector matTransVecMult(SiconosMatrix &m, SiconosVector &v)
{
  SimpleVector sv(m.size(1));
  if (v.isComposite())
  {
    SimpleVector tmp(v);
    Blas_Mat_Trans_Vec_Mult(m.getLaGenMatDouble(), tmp.lavd, sv.lavd, 1.0, 0.0);
  }
  else
    Blas_Mat_Trans_Vec_Mult(m.getLaGenMatDouble(), v.getValues(), sv.lavd, 1.0, 0.0);

  return sv;
}

