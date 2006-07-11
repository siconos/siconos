/* Siconos-Kernel version 1.2.0, Copyright INRIA 2005-2006.
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
#include "BlockVector.h"
using namespace std;

// CONSTRUCTORS
// Default
BlockVector::BlockVector(): SiconosVector()
{
  isBlockVector = true;
}

// From a file
BlockVector::BlockVector(const string file, const bool ascii): SiconosVector()
{
  string mode;
  isBlockVector = true;
  if (ascii) mode = "ascii";
  else mode = "binary";
  read(file, mode);
}

// copy from a Simple
BlockVector::BlockVector(const SimpleVector& v): SiconosVector()
{
  isBlockVector = true;
  svref.reserve(1);
  tabindex.reserve(1);
  svref.push_back(new SimpleVector(v));
  isSvrefAllocatedIn.push_back(true);
  tabindex.push_back(v.size());
}

// fill with a list of SimpleVector* (Warning: links between pointers!!)
BlockVector::BlockVector(vector<SimpleVector*> v): SiconosVector()
{
  isBlockVector = true;
  unsigned int numberOfVectors = v.size();
  svref.reserve(numberOfVectors);
  tabindex.reserve(numberOfVectors);
  for (unsigned int i = 0; i < numberOfVectors; ++i)
  {
    svref.push_back(v[i]);
    tabindex.push_back(v[i]->size());
    if (i > 0) tabindex[i] += tabindex[i - 1];
    isSvrefAllocatedIn.push_back(false);
  }
}

// link with 2 SimpleVector
BlockVector::BlockVector(SimpleVector* v1, SimpleVector* v2): SiconosVector()
{
  isBlockVector = true;
  svref.reserve(2);
  tabindex.reserve(2);
  if ((v1 == NULL) && (v2 == NULL))
    SiconosVectorException::selfThrow("BlockVector:constructor(SimpleVector*,SimpleVector*), both vectors are NULL.");

  if (v1 != NULL)
  {
    svref.push_back(v1);
    tabindex.push_back(v1->size());
    isSvrefAllocatedIn.push_back(false);
  }
  else
    // If first parameter is a NULL pointer, then set this(1) to a SimpleVector of the same size as v2, and equal to 0.
  {
    // This case is usefull to set xDot in LagrangianDS.
    svref.push_back(new SimpleVector(v2->size()));
    tabindex.push_back(v2->size());
    isSvrefAllocatedIn.push_back(true);
  }

  if (v2 != NULL)
  {
    svref.push_back(v2);
    tabindex.push_back(v2->size());
    tabindex[1] += tabindex[0];
    isSvrefAllocatedIn.push_back(false);
  }
  else // If second parameter is a NULL pointer, then set this(2) to a SimpleVector of the same size as v1, and equal to 0.
  {
    // This case is usefull to set xDot in LagrangianDS.
    svref.push_back(new SimpleVector(v1->size()));
    tabindex.push_back(v1->size());
    tabindex[1] += tabindex[0];
    isSvrefAllocatedIn.push_back(true);
  }
}

// copy
BlockVector::BlockVector(const BlockVector& v): SiconosVector()
{
  isBlockVector = true;
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
    SiconosVectorException::selfThrow("BlockVector: copy constructor, auto copy");
}

// with two int: the number of blocks and the dim of each block
BlockVector::BlockVector(unsigned int numberOfBlocks, unsigned int dimOfaBlock): SiconosVector()
{

  isBlockVector = true;
  svref.reserve(numberOfBlocks);
  tabindex.reserve(numberOfBlocks);
  for (unsigned int i = 0; i < numberOfBlocks; ++i)
  {
    svref.push_back(new SimpleVector(dimOfaBlock));
    tabindex.push_back((i + 1)*dimOfaBlock);
    isSvrefAllocatedIn.push_back(true);
  }
}

BlockVector::~BlockVector()
{
  for (unsigned int i = 0; i < svref.size(); i++)
  {
    if (isSvrefAllocatedIn[i]) delete svref[i];
    svref[i] = NULL;
  }
}

SimpleVector* BlockVector::getVectorPtr(const unsigned int i)
{
  if (i > svref.size())
    SiconosVectorException::selfThrow("BlockVector:getVectorPtr(i), i out of range.");

  return svref[i];
}

void BlockVector::display() const
{
  cout << "=== Block vector display === " << endl;
  const int sizeV = svref.size();
  cout << "It handles " << sizeV << " simple vector(s)" << endl;
  for (int i = 0; i < sizeV; i ++)
  {
    cout << "Vector number " << i << ":" << endl;
    svref[i]->display();
  }
  cout << endl << " ======================== " << endl;
}

double& BlockVector::operator()(const int unsigned index) const
{
  if (index > tabindex[tabindex.size() - 1])
    SiconosVectorException::selfThrow(" BlockVector::operator() -- index out of range");

  // locate which vector of svref corresponds to index
  unsigned int numVect = 0;
  while (index >= tabindex[numVect] && numVect < tabindex.size())
    numVect++;

  // get position of required index in vector numVect
  unsigned int pos = 0;
  if (numVect == 0)
    pos = index;
  else
    pos = index - tabindex[numVect - 1];

  return (*svref[numVect])(pos);
}

void BlockVector::add(const SimpleVector& v)
{
  // copy the vector into svref
  svref.push_back(new SimpleVector(v));
  isSvrefAllocatedIn.push_back(true);

  if (tabindex.size() > 0)
    tabindex.push_back(tabindex[tabindex.size() - 1] + v.size());
  else
    tabindex.push_back(v.size());
}

void BlockVector::addPtr(SimpleVector*v)
{
  svref.push_back(v);
  if (tabindex.size() > 0)
    tabindex.push_back(tabindex[tabindex.size() - 1] + v->size());
  else
    tabindex.push_back(v->size());
  isSvrefAllocatedIn.push_back(false);
}

void BlockVector::setValues(const vector<double>& v, const int unsigned i)
{
  if (svref.size() <= i)
    SiconosVectorException::selfThrow("BlockVector::setValues -- index out of range");
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


const LaVectorDouble BlockVector::getValues(const unsigned int i) const
{
  if (svref.size() <= i)
    SiconosVectorException::selfThrow("BlockVector::getValues -- index out of range");
  if (svref[i] == NULL)
    SiconosVectorException::selfThrow("BlockVector::getValues -- NULL vector");
  return svref[i]->getValues();
}

unsigned int BlockVector::size(const unsigned int i) const
{
  if (i > 1) SiconosVectorException::selfThrow("BlockVector::size(i) -- index i out of range");
  if (i == 0)
  {
    int sizeV = 0;
    if (tabindex.size() > 0) sizeV = tabindex[tabindex.size() - 1];
    return sizeV;
  }
  else
    return tabindex.size();
}

bool BlockVector::read(const string fileName, const string mode)
{
  // warning: this function reads the whole vector and overwrites existing values !!
  // to read only one part of the std::vector, use read for a SimpleVector and then
  // add it at the right place of the present block vector
  bool res = false;
  unsigned int i, j, size, nbVectors;
  double tmp;
  vector<double> vect;
  if (mode == "binary")
  {
    FILE * inFile = fopen(fileName.c_str(), "rb");    // open the input file in binary mode
    if (inFile == NULL)
      SiconosVectorException::selfThrow(" BlockVector::read : Fail to open file \"" + fileName + "\"");
    fread((char *) &nbVectors, sizeof(int), 1, inFile);   // read nbVectors
    if (nbVectors <= 0)
      SiconosVectorException::selfThrow(" BlockVector::read : try to read a vector with a negative size");
    svref.clear();
    tabindex.clear();
    isSvrefAllocatedIn.clear();
    // get values
    for (i = 0; i < nbVectors; i++)
    {
      fread((char *) &size, sizeof(int), 1, inFile);  // read size of vector i
      if (size <= 0)
        SiconosVectorException::selfThrow(" BlockVector::read : try to read a vector with a negative size");
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
      SiconosVectorException::selfThrow(" BlockVector::read : Failed to open file \"" + fileName + "\"");

    // get the number of vectors to load
    inFile >> nbVectors;
    if (nbVectors <= 0)
      SiconosVectorException::selfThrow(" BlockVector::read : try to read a vector with a negative size");
    svref.clear();
    tabindex.clear();
    isSvrefAllocatedIn.clear();
    // get values
    for (i = 0; i < nbVectors; i++)
    {
      inFile >> size; // read size of vector i
      if (size <= 0)
        SiconosVectorException::selfThrow(" BlockVector::read : try to read a vector with a negative size");
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


bool BlockVector::write(const string fileName, const  string mode) const
{
  bool res = false;
  if ((mode != "binary") && (mode != "ascii"))
    SiconosVectorException::selfThrow("BlockVector::write : unknown mode");

  // open the file
  ofstream outFile(fileName.c_str());

  if (!outFile.is_open())
    SiconosVectorException::selfThrow("BlockVector::write : : Fail to open file \"" + fileName + "\"");

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


double* BlockVector::getArray() const
{
  SiconosVectorException::selfThrow(" BlockVector::getArray, should not be used for a block vector.");
  return NULL;
}

void BlockVector::zero()
{
  vector<SimpleVector*>::iterator it;
  for (it = svref.begin(); it != svref.end(); ++it)
    (*it)->zero();
}

//===================================================================================
//          GENERIC INTERNAL OPERATORS
//===================================================================================

BlockVector &BlockVector::operator+=(const SiconosVector &v)
{
  unsigned int sizeV = size();
  if (sizeV != v.size())
    SiconosVectorException::selfThrow(" BlockVector::operator+=   -- the vectors have not the same size");

  for (unsigned int i = 0; i < sizeV; i++)
    (*this)(i) = (*this)(i) + v(i);

  return *this;
}

BlockVector &BlockVector::operator-=(const SiconosVector &v)
{
  unsigned int sizeV = size();
  if (sizeV != v.size())
    SiconosVectorException::selfThrow(" BlockVector::operator-=   -- the vectors have not the same size");

  for (unsigned int i = 0; i < sizeV; i++)
    (*this)(i) -= v(i);
  return *this;
}

BlockVector& BlockVector::operator = (const SiconosVector& v)
{
  if (this != &v)
  {
    unsigned int sizeV = size();
    if (sizeV != v.size())
      SiconosVectorException::selfThrow(" BlockVector::operator = -- the vectors have not the same size");
    else
    {
      for (unsigned int i = 0; i < sizeV; i++)
        (*this)(i) = v(i);
      if (v.isBlock())
        tabindex = v.getTabIndex();
    }
  }
  return *this;
}

//==============================================================================
//          SPECIFIC INTERNAL OPERATORS
//==============================================================================

BlockVector &BlockVector::operator*=(const double d)
{
  *this = *this * d;
  return *this;
}

BlockVector &BlockVector::operator/=(const double d)
{
  *this = *this / d;
  return *this;
}

BlockVector &BlockVector::operator+=(const BlockVector& v)
{
  if (v.size() != size())
    SiconosVectorException::selfThrow(" BlockVector::+=, the vectors have not the same size");
  *this = *this + v;
  return *this;
}


BlockVector &BlockVector::operator-=(const BlockVector& v)
{
  if (size() != v.size())
    SiconosVectorException::selfThrow(" BlockVector::-=, the vectors have not the same size");
  *this = *this - v;
  return *this;
}


BlockVector& BlockVector::operator = (const BlockVector& v)
{
  if (this != &v)
  {
    unsigned int sizeV = size();
    if (sizeV != v.size())
      SiconosVectorException::selfThrow(" BlockVector::operator = GENERIC  -- the vectors have not the same size");
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

BlockVector BlockVector::addition(const SiconosVector& v1) const
{
  unsigned int sizeV = size();

  if (sizeV != v1.size())
    SiconosVectorException::selfThrow(" BlockVector::addition(const BlockVector& v1)  --- the vectors have not the same size");

  BlockVector cv(*this);
  cv += v1;
  return cv;
}


BlockVector BlockVector::subtraction(const SiconosVector& v1) const
{
  unsigned int sizeV = size();

  if (sizeV != v1.size())
    SiconosVectorException::selfThrow(" BlockVector::subtraction(const BlockVector& v1)  --- the vectors have not the same size");

  BlockVector cv(*this);
  cv -= v1;
  return cv;
}


//=============================================================================
//          SPECIFIC EXTERNAL OPERATORS
//=============================================================================

BlockVector operator * (const BlockVector& v, const double d)
{
  BlockVector tmp(v);
  unsigned int size = v.size();
  for (unsigned int i = 0; i < size; i++)
    tmp(i) = d * v(i);
  return tmp;
}


BlockVector operator * (const double d, const BlockVector& v)
{
  BlockVector tmp(v);
  unsigned int size = v.size();
  for (unsigned int i = 0; i < size; i++)
    tmp(i) = d * v(i);
  return tmp;
}


BlockVector operator / (const BlockVector& v, const double d)
{
  if (d == 0.0)
    SiconosVectorException::selfThrow(" BlockVector operator/   --- division by 0");
  BlockVector tmp(v);
  unsigned int size = v.size();
  for (unsigned int i = 0; i < size; i++)
    tmp(i) = v(i) / d;
  return tmp;
}

BlockVector operator + (const BlockVector& v1, const BlockVector& v2)
{
  if (v1.getTabIndex() != v2.getTabIndex())
    SiconosVectorException::selfThrow(" BlockVector operator+  --- the vectors have not the same size or structure");

  BlockVector tmp(v1);
  unsigned int size = v1.size();
  for (unsigned int i = 0; i < size; i++)
    tmp(i) += v2(i);

  return tmp;
}

BlockVector operator - (const BlockVector& v1, const BlockVector& v2)
{
  if (v1.getTabIndex() != v2.getTabIndex())
    SiconosVectorException::selfThrow(" BlockVector operator+  --- the vectors have not the same size or structure");

  BlockVector tmp(v1);
  unsigned int size = v1.size();
  for (unsigned int i = 0; i < size; i++)
    tmp(i) -= v2(i);

  return tmp;
}



SimpleVector operator * (const SiconosMatrix &m, const BlockVector &v)
{
  // Check sizes
  if (v.size() != m.size(1))
    SiconosVectorException::selfThrow(" BlockVector, SiconosMatrix*BlockVector: inconsistent sizes");
  SimpleVector tmp(v);
  return m * tmp;
}

SimpleVector matTransVecMult(SiconosMatrix &m, BlockVector &v)
{
  // Check sizes
  if (v.size() != m.size(0))
    SiconosVectorException::selfThrow(" BlockVector, SiconosMatrix*BlockVector: inconsistent sizes");

  SimpleVector tmp(v);
  return matTransVecMult(m, tmp);
}


