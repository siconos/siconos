/* Siconos-Kernel version 1.1.1, Copyright INRIA 2005-2006.
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
#include "CompositeVector.h"
using namespace std;

// CONSTRUCTORS
// Default
CompositeVector::CompositeVector(): SiconosVector()
{
  composite = true;
}

// From a file
CompositeVector::CompositeVector(const string& file, const bool& ascii): SiconosVector()
{
  string mode;
  composite = true;
  if (ascii) mode = "ascii";
  else mode = "binary";
  read(file, mode);
}

// copy from a Simple
CompositeVector::CompositeVector(const SimpleVector& v): SiconosVector()
{
  composite = true;
  svref.reserve(1);
  tabindex.reserve(1);
  svref.push_back(new SimpleVector(v));
  isSvrefAllocatedIn.push_back(true);
  tabindex.push_back(v.size());
}


// copy
CompositeVector::CompositeVector(const CompositeVector& v): SiconosVector()
{
  composite = true;
  if (this != &v)
  {
    unsigned int sizeV = v.size(1);
    svref.reserve(sizeV);
    vector<SimpleVector*> tmpVector = v.getSvref();

    for (unsigned int i = 0; i < sizeV; i++)
    {
      svref.push_back(new SimpleVector(*(tmpVector[i])));
      isSvrefAllocatedIn.push_back(true);
    }
    tabindex.reserve(sizeV);
    tabindex = v.getTabIndex();
  }
  else
    SiconosVectorException::selfThrow("CompositeVector: copy constructor, auto copy");
}

CompositeVector::~CompositeVector()
{
  for (unsigned int i = 0; i < svref.size(); i++)
  {
    if (isSvrefAllocatedIn[i]) delete svref[i];
    svref[i] = NULL;
  }
}

void CompositeVector::display() const
{
  cout << "=== Composite Display === " << endl;
  const int sizeV = svref.size();
  cout << "with " << sizeV << " simple vectors" << endl;
  for (int i = 0; i < sizeV; i ++)
  {
    cout << "Vector number " << i << endl;
    svref[i]->display();
  }
  cout << endl << " ======================== " << endl;
}

double& CompositeVector::operator()(const int unsigned& index) const
{
  if ((int)index > tabindex[tabindex.size() - 1])
    SiconosVectorException::selfThrow(" CompositeVector::operator() -- index out of range");

  // locate which vector of svref corresponds to index
  unsigned int numVect = 0;
  while ((int)index >= tabindex[numVect] && numVect < tabindex.size()) numVect++;
  // get position of required index in vector numVect
  unsigned int pos = 0;
  if (numVect == 0)
    pos = index;
  else
    pos = index - tabindex[numVect - 1];

  return (*svref[numVect])(pos);
}

void CompositeVector::add(const SimpleVector& v)
{
  // copy the vector into svref
  svref.push_back(new SimpleVector(v));
  isSvrefAllocatedIn.push_back(true);

  if (tabindex.size() > 0)
    tabindex.push_back(tabindex[tabindex.size() - 1] + v.size());
  else
    tabindex.push_back(v.size());
}

void CompositeVector::addPtr(SimpleVector*v)
{
  svref.push_back(v);
  if (tabindex.size() > 0)
    tabindex.push_back(tabindex[tabindex.size() - 1] + v->size());
  else
    tabindex.push_back(v->size());
  isSvrefAllocatedIn.push_back(false);
}

void CompositeVector::setValues(const vector<double>& v, const int unsigned& i)
{
  if (svref.size() <= i)
    SiconosVectorException::selfThrow("CompositeVector::setValues -- index out of range");
  unsigned int oldSize;
  if (isSvrefAllocatedIn[i])
  {
    oldSize = svref[i]->size();
    delete svref[i];
  }
  svref[i] = new SimpleVector(v.size());
  isSvrefAllocatedIn[i] = true;
  svref[i]->setValues(v);
  // update tabindex
  for (unsigned int j = i; j < tabindex.size(); j++)
    tabindex[j] = tabindex[j] - oldSize + svref[i]->size();
}


const LaVectorDouble CompositeVector::getValues(const unsigned int& i) const
{
  if (svref.size() <= i)
    SiconosVectorException::selfThrow("CompositeVector::getValues -- index out of range");
  if (svref[i] == NULL)
    SiconosVectorException::selfThrow("CompositeVector::getValues -- NULL vector");
  return svref[i]->getValues();
}

unsigned int CompositeVector::size(const unsigned int& i) const
{
  if (i > 1) SiconosVectorException::selfThrow("CompositeVector::size(i) -- index i out of range");
  if (i == 0)
  {
    int sizeV = 0;
    if (tabindex.size() > 0) sizeV = tabindex[tabindex.size() - 1];
    return sizeV;
  }
  else
    return tabindex.size();
}

bool CompositeVector::read(const string& fileName, const string& mode)
{
  // warning: this function reads the whole vector and overwrites existing values !!
  // to read only one part of the std::vector, use read for a SimpleVector and then
  // add it at the right place of the present composite
  bool res = false;
  unsigned int i, j, size, nbVectors;
  double tmp;
  vector<double> vect;
  if (mode == "binary")
  {
    FILE * inFile = fopen(fileName.c_str(), "rb");    // open the input file in binary mode
    if (inFile == NULL)
      SiconosVectorException::selfThrow(" CompositeVector::read : Fail to open file \"" + fileName + "\"");
    fread((char *) &nbVectors, sizeof(int), 1, inFile);   // read nbVectors
    if (nbVectors <= 0)
      SiconosVectorException::selfThrow(" CompositeVector::read : try to read a vector with a negative size");
    svref.clear();
    tabindex.clear();
    isSvrefAllocatedIn.clear();
    // get values
    for (i = 0; i < nbVectors; i++)
    {
      fread((char *) &size, sizeof(int), 1, inFile);  // read size of vector i
      if (size <= 0)
        SiconosVectorException::selfThrow(" CompositeVector::read : try to read a vector with a negative size");
      // Memory allocation for vector i
      svref.push_back(new SimpleVector(size));
      isSvrefAllocatedIn.push_back(true);
      for (j = 0; j < size; j++)
      {
        fread((char*) &tmp, sizeof(double), 1, inFile);  // read a double
        vect.push_back(tmp);
      }
      svref[i]->setValues(vect);
      if (tabindex.size() == 0)
        tabindex.push_back(vect.size());
      else
        tabindex.push_back(vect.size() + tabindex[i - 1]);
      vect.clear();
    }
    fclose(inFile);
    res = true;
  }
  else if (mode == "ascii")
  {
    ifstream inFile(fileName.c_str(),  ifstream::in);
    if (inFile == NULL)
      SiconosVectorException::selfThrow(" CompositeVector::read : Failed to open file \"" + fileName + "\"");

    // get the number of vectors to load
    inFile >> nbVectors;
    if (nbVectors <= 0)
      SiconosVectorException::selfThrow(" CompositeVector::read : try to read a vector with a negative size");
    svref.clear();
    tabindex.clear();
    isSvrefAllocatedIn.clear();
    // get values
    for (i = 0; i < nbVectors; i++)
    {
      inFile >> size; // read size of vector i
      if (size <= 0)
        SiconosVectorException::selfThrow(" CompositeVector::read : try to read a vector with a negative size");
      // Memory allocation for vector i
      svref.push_back(new SimpleVector(size));
      isSvrefAllocatedIn.push_back(true);
      for (j = 0; j < size; j++)
      {
        inFile >> tmp; // read a double
        vect.push_back(tmp);
      }
      svref[i]->setValues(vect);
      if (tabindex.size() == 0)
        tabindex.push_back(vect.size());
      else
        tabindex.push_back(vect.size() + tabindex[i - 1]);
      vect.clear();
    }
    inFile.close();
    res = true;
  }
  return res;
}


bool CompositeVector::write(const string& fileName, const  string& mode) const
{
  bool res = false;
  if ((mode != "binary") && (mode != "ascii"))
    SiconosVectorException::selfThrow("CompositeVector::write : unknown mode");

  // open the file
  ofstream outFile(fileName.c_str());

  if (!outFile.is_open())
    SiconosVectorException::selfThrow("CompositeVector::write : : Fail to open file \"" + fileName + "\"");

  unsigned int nbVectors = svref.size();

  if (mode == "binary")
  {
    outFile.write((char*)&nbVectors, sizeof(int));
    unsigned int sizeV = 0;
    for (unsigned int i = 0; i < nbVectors; i++)
    {
      sizeV = svref[i]->size();
      outFile.write((char*)&sizeV, sizeof(int));
      for (unsigned int j = 0; j < sizeV; j++)
        outFile.write((char*) & (*(svref[i]))(j), sizeof(double));
    }
    res = true;
  }
  else if (mode == "ascii")
  {
    outFile << nbVectors << endl;

    int sizeV = 0;
    for (unsigned int i = 0; i < nbVectors; i++)
    {
      sizeV = svref[i]->size();
      outFile << sizeV << endl;
      outFile << (getValues(i));
    }
    res = true;
  }
  outFile.close();
  return res;
}


double* CompositeVector::getArray() const
{
  SiconosVectorException::selfThrow(" CompositeVector::getArray, should not be used for composite vector");
  return NULL;
}

void CompositeVector::zero()
{
  vector<SimpleVector*>::iterator it;
  for (it = svref.begin(); it != svref.end(); ++it)
    (*it)->zero();
}

//===================================================================================
//          GENERIC INTERNAL OPERATORS
//===================================================================================

CompositeVector &CompositeVector::operator+=(const SiconosVector &v)
{
  unsigned int sizeV = size();
  if (sizeV != v.size())
    SiconosVectorException::selfThrow(" CompositeVector::operator+=   -- the vectors have not the same size");

  for (unsigned int i = 0; i < sizeV; i++)
    (*this)(i) = (*this)(i) + v(i);

  return *this;
}

CompositeVector &CompositeVector::operator-=(const SiconosVector &v)
{
  unsigned int sizeV = size();
  if (sizeV != v.size())
    SiconosVectorException::selfThrow(" CompositeVector::operator-=   -- the vectors have not the same size");

  for (unsigned int i = 0; i < sizeV; i++)
    (*this)(i) -= v(i);
  return *this;
}

CompositeVector& CompositeVector::operator = (const SiconosVector& v)
{
  if (this != &v)
  {
    unsigned int sizeV = size();
    if (sizeV != v.size())
      SiconosVectorException::selfThrow(" CompositeVector::operator = -- the vectors have not the same size");
    else
    {
      for (unsigned int i = 0; i < sizeV; i++)
        (*this)(i) = v(i);
      if (v.isComposite())
        tabindex = v.getTabIndex();
    }
  }
  return *this;
}


bool CompositeVector::operator == (const SiconosVector& v) const
{
  bool res = true;
  if (this != &v)
  {
    unsigned int sizeV = size();
    if (sizeV == v.size())
    {
      if (sizeV != 0)
      {
        for (unsigned int i = 0; i < sizeV; i++)
        {
          if ((*this)(i) != v(i))
          {
            res = false;
            break;
          }
        }
      }
    }
    else res =  false;
  }
  return res;
}

bool CompositeVector::operator == (const CompositeVector& v) const
{
  bool res = true;
  if (this != &v)
  {
    unsigned int sizeV = size();
    if (v.getTabIndex() == getTabIndex())
    {
      if (sizeV != 0)
      {
        for (unsigned int i = 0; i < sizeV; i++)
        {
          if ((*this)(i) != v(i))
          {
            res = false;
            break;
          }
        }
      }
    }
    else res =  false;
  }
  return res;
}

bool CompositeVector::operator != (const SiconosVector& v) const
{
  return !(*this == v);
}

bool CompositeVector::operator != (const CompositeVector& v) const
{
  return !(*this == v);
}

//==============================================================================
//          SPECIFIC INTERNAL OPERATORS
//==============================================================================

CompositeVector &CompositeVector::operator*=(const double& d)
{
  *this = *this * d;
  return *this;
}

CompositeVector &CompositeVector::operator/=(const double& d)
{
  *this = *this / d;
  return *this;
}

CompositeVector &CompositeVector::operator+=(const CompositeVector& v)
{
  if (v.size() != size())
    SiconosVectorException::selfThrow(" CompositeVector::+=, the vectors have not the same size");
  *this = *this + v;
  return *this;
}


CompositeVector &CompositeVector::operator-=(const CompositeVector& v)
{
  if (size() != v.size())
    SiconosVectorException::selfThrow(" CompositeVector::-=, the vectors have not the same size");
  *this = *this - v;
  return *this;
}


CompositeVector& CompositeVector::operator = (const CompositeVector& v)
{
  if (this != &v)
  {
    unsigned int sizeV = size();
    if (sizeV != v.size())
      SiconosVectorException::selfThrow(" CompositeVector::operator = GENERIC  -- the vectors have not the same size");
    else
    {
      for (unsigned int i = 0; i < sizeV; i++)
        (*this)(i) = v(i);
      tabindex = v.getTabIndex();
    }
  }
  return *this;
}

//=============================================================================
//          GENERIC EXTERNAL OPERATORS
//=============================================================================

CompositeVector CompositeVector::addition(const SiconosVector& v1) const
{
  unsigned int sizeV = size();

  if (sizeV != v1.size())
    SiconosVectorException::selfThrow(" CompositeVector::addition(const CompositeVector& v1)  --- the vectors have not the same size");

  CompositeVector cv(*this);
  cv += v1;
  return cv;
}


CompositeVector CompositeVector::subtraction(const SiconosVector& v1) const
{
  unsigned int sizeV = size();

  if (sizeV != v1.size())
    SiconosVectorException::selfThrow(" CompositeVector::subtraction(const CompositeVector& v1)  --- the vectors have not the same size");

  CompositeVector cv(*this);
  cv -= v1;
  return cv;
}


//=============================================================================
//          SPECIFIC EXTERNAL OPERATORS
//=============================================================================

CompositeVector operator * (const CompositeVector& v, const double& d)
{
  CompositeVector tmp(v);
  unsigned int size = v.size();
  for (unsigned int i = 0; i < size; i++)
    tmp(i) = d * v(i);
  return tmp;
}


CompositeVector operator * (const double& d, const CompositeVector& v)
{
  CompositeVector tmp(v);
  unsigned int size = v.size();
  for (unsigned int i = 0; i < size; i++)
    tmp(i) = d * v(i);
  return tmp;
}


CompositeVector operator / (const CompositeVector& v, const double& d)
{
  if (d == 0.0)
    SiconosVectorException::selfThrow(" CompositeVector operator/   --- division by 0");
  CompositeVector tmp(v);
  unsigned int size = v.size();
  for (unsigned int i = 0; i < size; i++)
    tmp(i) = v(i) / d;
  return tmp;
}


CompositeVector operator + (const CompositeVector& v1, const CompositeVector& v2)
{
  if (v1.getTabIndex() != v2.getTabIndex())
    SiconosVectorException::selfThrow(" CompositeVector operator+  --- the vectors have not the same size or structure");

  CompositeVector tmp(v1);
  unsigned int size = v1.size();
  for (unsigned int i = 0; i < size; i++)
    tmp(i) += v2(i);

  return tmp;
}

CompositeVector operator - (const CompositeVector& v1, const CompositeVector& v2)
{
  if (v1.getTabIndex() != v2.getTabIndex())
    SiconosVectorException::selfThrow(" CompositeVector operator+  --- the vectors have not the same size or structure");

  CompositeVector tmp(v1);
  unsigned int size = v1.size();
  for (unsigned int i = 0; i < size; i++)
    tmp(i) -= v2(i);

  return tmp;
}

SimpleVector operator * (const SiconosMatrix &m, const CompositeVector &v)
{
  // Check sizes
  if (v.size() != m.size(1))
    SiconosVectorException::selfThrow(" CompositeVector, SiconosMatrix*CompositeVector: inconsistent sizes");
  SimpleVector tmp(v);
  return m * tmp;
}

SimpleVector matTransVecMult(SiconosMatrix &m, CompositeVector &v)
{
  // Check sizes
  if (v.size() != m.size(0))
    SiconosVectorException::selfThrow(" CompositeVector, SiconosMatrix*CompositeVector: inconsistent sizes");

  SimpleVector tmp(v);
  return matTransVecMult(m, tmp);
}


