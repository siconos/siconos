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
#include "SimpleVector.h"
using namespace std;

// default (private)
SimpleVector::SimpleVector(): SiconosVector()
{
  isBlockVector = false;
  lavd = LaVectorDouble();
}

SimpleVector::SimpleVector(const string& file, const bool& ascii):
  SiconosVector()
{
  isBlockVector = false;
  string mode;
  if (ascii) mode = "ascii";
  else mode = "binary";
  read(file, mode);
}

// copy from a std vector
SimpleVector::SimpleVector(const vector<double>& v):
  SiconosVector()
{
  isBlockVector = false;
  setValues(v);
}

// copy
SimpleVector::SimpleVector(const SimpleVector& v):
  SiconosVector()
{
  isBlockVector = false;
  lavd = v.getValues();
}

// copy from SiconosVector
SimpleVector::SimpleVector(const SiconosVector& v):
  SiconosVector()
{
  if (!v.isBlock())
  {
    isBlockVector = false;
    lavd = v.getValues();
  }
  else
  {
    isBlockVector = false;
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
  isBlockVector = false;
  // resize and init lavd
  lavd.resize(size, 1);
  for (unsigned int i = 0; i < size; ++i)
    lavd(i) = 0.0;
}

SimpleVector::~SimpleVector()
{
}

/***********************************************************************************************/

void SimpleVector::display() const
{
  cout << "=== Simple vector of size " << size() << ":" << endl;
  if (size() <= M_MAXSIZEFORDISPLAY)
  {
    for (unsigned int i = 0; i < size(); i++)
      cout << lavd(i) << " ";
  }
  else cout << "Display SiconosVector : vector too large" << endl;
  cout << endl << "=================================== " << endl;
}

double& SimpleVector::operator()(const int unsigned& index) const
{
  if ((int)index >= lavd.size())
    SiconosVectorException::selfThrow(" SimpleVector::operator()   -- index greater than vector size");
  return const_cast<double&>(lavd(index));
}

std::vector<unsigned int> SimpleVector::getTabIndex() const
{
  vector<unsigned int> tmp;
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

void SimpleVector::setValues(const vector<double>& v, const int unsigned& index) // index default value = 0, useless for simpleVector BUT required for block vector
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
    SiconosVectorException::selfThrow("SimpleVector::getBlock : wrong indexes list");

  unsigned int k = 0;

  vector<unsigned int>::const_iterator it;

  for (it = index.begin(); it != index.end(); it++)
  {
    if (*it >= size()) SiconosVectorException::selfThrow("SimpleVector::getBlock : index out of range");
    block(k) = lavd(*it);
    k++;
  }
}

void SimpleVector::getBlock(const unsigned int& index , const unsigned int& nbval , SimpleVector& block) const
{
  if ((block.size() != nbval) || (index + nbval > size()))
    SiconosVectorException::selfThrow("SimpleVector::getBlock : wrong index or sizes");

  for (unsigned int i = 0; i < nbval; i++)
    block(i) = lavd(i + index);

}

void SimpleVector::setBlock(const unsigned int& posi, const SiconosVector &v)
{
  if (v.isBlock())
    SiconosVectorException::selfThrow("SimpleVector::setBlock : argument vector is composite");

  if (posi + v.size() > (unsigned int)lavd.size())
    SiconosVectorException::selfThrow("SimpleVector::setBlock : filling out of range");

  for (unsigned int i = 0; i < v.size(); i++)
    lavd(i + posi) = v.getValue(i);
}

unsigned int SimpleVector::size(const unsigned int& i) const
{
  if (i != 0)  SiconosVectorException::selfThrow("SimpleVector::size(i) - i!=0 -> not available");
  return lavd.size();
}

bool SimpleVector::read(const string& fileName, const string& mode)
{
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
  return res;
}


bool SimpleVector::write(const string& fileName, const string& mode) const
{
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
  unsigned int sizeV = size();
  if (sizeV != v.size())
    SiconosVectorException::selfThrow(" SimpleVector::operator+=   -- the vectors have not the same size");

  for (unsigned int i = 0; i < sizeV; i++)
    (*this)(i) += v(i);
  return *this;
}


SimpleVector &SimpleVector::operator-=(const SiconosVector &v)
{
  unsigned int sizeV = size();
  if (sizeV != v.size())
    SiconosVectorException::selfThrow(" SimpleVector::operator-=   -- the vectors have not the same size");

  for (unsigned int i = 0; i < sizeV; i++)
    (*this)(i) -= v(i);
  return *this;
}


SimpleVector& SimpleVector::operator = (const SiconosVector& v)
{
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
  return *this;
}

/*******************************************************************************
 *          SPECIFIC INTERNAL OPERATORS                                *
 *******************************************************************************/
SimpleVector& SimpleVector::operator = (const SimpleVector& v)
{
  if (v.size() != size())
    SiconosVectorException::selfThrow("SimpleVector::operator= inconsistent sizes");

  if (this != &v)
    lavd.copy(v.lavd);

  return *this;
}


SimpleVector &SimpleVector::operator+=(const SimpleVector &v)
{
  if (size() != v.size())
    SiconosVectorException::selfThrow(" SimpleVector::operator+=   -- the vectors have not the same size");

  Blas_Add_Mult(lavd, 1.0, v.lavd);

  return *this;
}


SimpleVector &SimpleVector::operator-=(const SimpleVector &v)
{
  if (size() != v.size())
    SiconosVectorException::selfThrow(" SimpleVector::operator-=   -- the vectors have not the same size");

  Blas_Add_Mult(lavd, -1.0, v.lavd);
  return *this;
}


SimpleVector &SimpleVector::operator*=(const double& d)
{
  lavd = lavd * d;
  return *this;
}


SimpleVector &SimpleVector::operator/=(const double& d)
{
  if (d == 0)
    SiconosVectorException::selfThrow(" SimpleVector::operator/ =   --division by 0");
  lavd = lavd * (1.0 / d);
  return *this;
}

bool operator==(const SiconosVector& v1, const SiconosVector& v2)
{
  double norm = (v1 - v2).norm();
  return(norm < tolerance);
}

//==================================================================================
//      GENERIC INTERNAL FUNCTIONS FOR MIXED OPERATIONS
//==================================================================================

SimpleVector SimpleVector::addition(const SiconosVector& v) const
{
  if (size() != v.size())
    SiconosVectorException::selfThrow(" SimpleVector::addition(const SiconosVector& v) --- the vectors have not the same size");

  SimpleVector sv(*this);
  sv += v;
  return sv;
}


SimpleVector SimpleVector::subtraction(const SiconosVector& v) const
{
  if (size() != v.size())
    SiconosVectorException::selfThrow(" SimpleVector::subtraction(const SiconosVector& v) --- the vectors have not the same size");

  SimpleVector sv(*this);
  sv -= v;
  return sv;
}

SimpleVector operator + (const SiconosVector& v1, const SiconosVector& v2)
{
  if (v1.size() != v2.size())
    SiconosVectorException::selfThrow(" SiconosVector operator + --- the vectors have not the same size");

  SimpleVector sv(v1);
  sv += v2;

  return sv;
}


SimpleVector operator - (const SiconosVector& v1, const SiconosVector& v2)
{
  if (v1.size() != v2.size())
    SiconosVectorException::selfThrow(" SiconosVector operator - --- the vectors have not the same size");

  SimpleVector sv(v1);
  sv -= v2;
  return sv;
}

/*******************************************************************************
 *          SPECIFIC EXTERNAL OPERATORS                                *
 *******************************************************************************/
SimpleVector operator * (const SimpleVector& v, const double& d)
{
  SimpleVector sv(v);
  Blas_Scale(d, sv.lavd);
  return sv;
}


SimpleVector operator * (const double& d, const SimpleVector& v)
{
  SimpleVector sv(v);
  Blas_Scale(d, sv.lavd);
  return sv;
}


SimpleVector operator / (const SimpleVector& v, const double& d)
{
  if (d == 0.0)
    SiconosVectorException::selfThrow(" SimpleVector operator / --- division by 0");

  SimpleVector sv(v);
  Blas_Scale(1 / d, sv.lavd);
  return sv;
}

// SimpleVector operator + (const SimpleVector& v1, const SimpleVector& v2)
// {
//   if (v1.size() != v2.size())
//     SiconosVectorException::selfThrow(" SimpleVector operator + --- vectors have not the same size");

//   SimpleVector sv(v1);
//   Blas_Add_Mult(sv.lavd, 1.0, v2.lavd);

//   return sv;
// }


// SimpleVector operator - (const SimpleVector& v1, const SimpleVector& v2)
// {
//   if (v1.size() != v2.size())
//     SiconosVectorException::selfThrow(" SimpleVector operator - --- the 2 vectors have not the same size");

//   SimpleVector sv(v1);
//   Blas_Add_Mult(sv.lavd, -1.0, v2.lavd);
//   return sv;
// }


SimpleVector operator * (const SiconosMatrix &m, const SimpleVector &v)
{
  SimpleVector sv(m.size(0));
  Blas_Mat_Vec_Mult(*m.getLaGenMatDoubleRef(), v.lavd, sv.lavd, 1.0, 0.0);
  return sv;
}

SimpleVector operator * (const SiconosMatrix &m, const SiconosVector &v)
{
  SimpleVector sv(m.size(0));
  // "transform" v (which is simple or block vector) into a simple
  // not good -> to be optimized
  SimpleVector tmp(v);
  Blas_Mat_Vec_Mult(*m.getLaGenMatDoubleRef(), tmp.lavd, sv.lavd, 1.0, 0.0);
  return sv;
}

SimpleVector operator * (const SimpleVector &v, const SiconosMatrix &m)
{
  SimpleVector sv(m.size(1));
  Blas_Mat_Trans_Vec_Mult(*m.getLaGenMatDoubleRef(), v.lavd, sv.lavd, 1.0, 0.0);
  return sv;
}

SimpleVector operator * (const SiconosVector &v, const SiconosMatrix &m)
{
  SimpleVector sv(m.size(1));
  // "transform" v (which is simple or block vector) into a simple
  // not good -> to be optimized
  SimpleVector tmp(v);
  Blas_Mat_Trans_Vec_Mult(*m.getLaGenMatDoubleRef(), tmp.lavd, sv.lavd, 1.0, 0.0);
  return sv;
}

SimpleVector matTransVecMult(SiconosMatrix &m, SimpleVector &v)
{
  SimpleVector sv(m.size(1));
  Blas_Mat_Trans_Vec_Mult(*m.getLaGenMatDoubleRef(), v.lavd, sv.lavd, 1.0, 0.0);
  return sv;
}

SimpleVector matTransVecMult(SiconosMatrix &m, SiconosVector &v)
{
  SimpleVector sv(m.size(1));
  if (v.isBlock())
  {
    SimpleVector tmp(v);
    Blas_Mat_Trans_Vec_Mult(*m.getLaGenMatDoubleRef(), tmp.lavd, sv.lavd, 1.0, 0.0);
  }
  else
    Blas_Mat_Trans_Vec_Mult(*m.getLaGenMatDoubleRef(), v.getValues(), sv.lavd, 1.0, 0.0);

  return sv;
}

