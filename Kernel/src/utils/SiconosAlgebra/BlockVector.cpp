/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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
#include "SimpleVector.h"
#include "SiconosVectorException.h"

/***************************** CONSTRUCTORS ****************************/
// Default (private)
BlockVector::BlockVector(): SiconosVector(true)
{}

BlockVector::BlockVector(const std::string & file, bool ascii): SiconosVector(true)
{
  SiconosVectorException::selfThrow(" BlockVector::constructor from a file : read BlockVector is not implemented");
}

BlockVector::BlockVector(const BlockVector &v): SiconosVector(true)
{
  unsigned int numberOfBlocks;

  tabIndex = v.getTabIndex();
  numberOfBlocks = tabIndex.size();
  for (unsigned int i = 0; i < numberOfBlocks; ++i)
  {
    // Call copy-constructor of SimpleVector
    vect.push_back(new SimpleVector(*v[i])) ;
    isBlockAllocatedIn.push_back(true);
  }
}

BlockVector::BlockVector(const SiconosVector &v): SiconosVector(true)
{
  if (v.isBlock())
  {
    //       unsigned int numberOfBlocks;

    //       SiconosVector* tmp1 = const_cast<SiconosVector*>(&v);
    //       BlockVector * tmp = static_cast<BlockVector*>(tmp1);
    //       tabIndex = tmp->getTabIndex();
    //       numberOfBlocks = tabIndex.size();
    //       for (unsigned int i=0;i < numberOfBlocks; ++i)
    //  {
    //    // Call copy-constructor of SimpleVector
    //    vect.push_back(new SimpleVector(*(tmp->getVectorPtr(i)))) ;
    //    isBlockAllocatedIn.push_back(true);
    //  }
    unsigned int numberOfBlocks = v.getNumberOfBlocks();;
    tabIndex.resize(numberOfBlocks);
    for (unsigned int i = 0; i < numberOfBlocks; ++i)
    {
      // Call copy-constructor of SimpleVector
      vect.push_back(new SimpleVector(*v[i])) ;
      isBlockAllocatedIn.push_back(true);
      tabIndex[i] = v[i]->size();
      if (i > 0) tabIndex[i] += tabIndex[i - 1];
    }
  }
  else
  {
    // Call copy-constructor of SimpleVector
    vect.push_back(new SimpleVector(v)) ;
    isBlockAllocatedIn.push_back(true);
    tabIndex.push_back(v.size());
  }
}

BlockVector::BlockVector(SiconosVector* v1, SiconosVector* v2): SiconosVector(true)
{
  // Add the two vectors in the container
  // NO COPY !!
  if ((v1 == NULL) && (v2 == NULL))
    SiconosVectorException::selfThrow("BlockVector:constructor(SimpleVector*,SimpleVector*), both vectors are NULL.");

  if (v1 != NULL)
  {
    if (v1->isBlock())
      SiconosVectorException::selfThrow("BlockVector::constructor(v1,v2) : v1 is a block vector and cannot be added into another block vector.");
    vect.push_back(v1);
    tabIndex.push_back(v1->size());
    isBlockAllocatedIn.push_back(false);
  }
  else
    // If first parameter is a NULL pointer, then set this(1) to a SimpleVector of the same size as v2, and equal to 0.
  {
    // This case is usefull to set xDot in LagrangianDS.
    vect.push_back(new SimpleVector(v2->size()));
    tabIndex.push_back(v2->size());
    isBlockAllocatedIn.push_back(true);
  }
  if (v2 != NULL)
  {
    if (v2->isBlock())
      SiconosVectorException::selfThrow("BlockVector::constructor(v1,v2) : v2 is a block vector and cannot be added into another block vector.");
    vect.push_back(v2);
    tabIndex.push_back(v2->size());
    tabIndex[1] += tabIndex[0];
    isBlockAllocatedIn.push_back(false);
  }
  else // If second parameter is a NULL pointer, then set this(2) to a SimpleVector of the same size as v1, and equal to 0.
  {
    // This case is usefull to set xDot in LagrangianDS.
    vect.push_back(new SimpleVector(v1->size()));
    tabIndex.push_back(v1->size());
    tabIndex[1] += tabIndex[0];
    isBlockAllocatedIn.push_back(true);
  }
}

BlockVector::BlockVector(unsigned int numberOfBlocks, unsigned int dim): SiconosVector(true)
{
  vect.resize(numberOfBlocks, NULL);
  BlockVectIterator it;
  unsigned int i = 1;
  for (it = vect.begin(); it != vect.end(); ++it)
  {
    *it = new SimpleVector(dim);
    tabIndex.push_back(dim * i);
    i++;
    isBlockAllocatedIn.push_back(true);
  }
}

/****************************** DESTRUCTOR  ****************************/
BlockVector::~BlockVector()
{
  unsigned int numberOfBlocks = tabIndex.size();
  for (unsigned int i = 0; i < numberOfBlocks; ++i)
    if (isBlockAllocatedIn[i]) delete vect[i];
  vect.clear();
  tabIndex.clear();
  isBlockAllocatedIn.clear();
}

/******************************** METHODS ******************************/
unsigned int BlockVector::getNum() const
{
  return 0;
}

const DenseVect BlockVector::getDense(unsigned int i) const
{
  if (vect[i]->getNum() != 1)
    SiconosVectorException::selfThrow("BlockVector::getDense(unsigned int num) : the vector[num] is not a Dense vector");
  return (vect[i])->getDense();
}

const SparseVect BlockVector::getSparse(unsigned int i)const
{
  if (vect[i]->getNum() != 4)
    SiconosVectorException::selfThrow("BlockVector::getSparse(unsigned int num) : the vector[num] is not a Sparse vector");
  return (vect[i])->getSparse();
}

DenseVect* BlockVector::getDensePtr(unsigned int i) const
{
  if (vect[i]->getNum() != 1)
    SiconosVectorException::selfThrow("BlockVector::getDensePtr(unsigned int num) : the vector[num] is not a Dense vector");
  return (vect[i])->getDensePtr();
}

SparseVect* BlockVector::getSparsePtr(unsigned int i)const
{

  if (vect[i]->getNum() != 4)
    SiconosVectorException::selfThrow("BlockVector::getSparsePtr(unsigned int num) : the vector[num] is not a Sparse vector");
  return (vect[i])->getSparsePtr();
}

double* BlockVector::getArray(unsigned int i) const
{
  return (vect[i])->getArray();
}

void BlockVector::getBlock(unsigned int, SiconosVector&) const
{
  SiconosVectorException::selfThrow("BlockVector::getBlock(int,SiconosVector): not yet implemented.");
}

void BlockVector::zero()
{
  BlockVectIterator it;
  for (it = vect.begin(); it != vect.end(); ++it)
    (*it)->zero();
}

bool BlockVector::check() const
{
  bool res = true;

  for (unsigned int i = 0; i < tabIndex.size(); ++i)
  {
    if (vect[i] == NULL)
    {
      res = false;
      std::cout << "BlockVector warning: one of the block is a NULL pointer, at position:" << i << std::endl;
    }
  }
  return res;
}

void BlockVector::resize(unsigned int, bool)
{
  SiconosVectorException::selfThrow("BlockVector::resize, not allowed for block vectors.");
}

const double BlockVector::normInf() const
{
  double d = 0, tmp;
  ConstBlockVectIterator it;
  for (it = vect.begin(); it != vect.end(); ++it)
  {
    tmp = (*it)->normInf();
    if (tmp > d) d = tmp;
  }
  return d;
}


const double BlockVector::norm() const
{
  double d = 0;
  ConstBlockVectIterator it;
  for (it = vect.begin(); it != vect.end(); ++it)
  {
    if ((*it) != NULL)
      d += pow((*it)->norm(), 2);
    else
      SiconosVectorException::selfThrow("BlockVector::norm, one of the blocks is equal to NULL pointer.");
  }
  return sqrt(d);
}

void BlockVector::display() const
{
  ConstBlockVectIterator it;
  std::cout << "=======> Block Vector Display (" << tabIndex.size() << " block(s)): " << std::endl;
  for (it = vect.begin(); it != vect.end(); ++it)
    (*it)->display();
}

SimpleVector BlockVector::getVector(unsigned int num) const
{
  if (vect[num] == NULL)
    SiconosVectorException::selfThrow("BlockVector::getVector(num), vector[num] == NULL pointer.");

  return *(vect[num]);
}

SiconosVector * BlockVector::getVectorPtr(unsigned int num)
{
  return vect[num];
}

void BlockVector::fill(double value)
{
  BlockVectIterator it;
  for (it = vect.begin(); it != vect.end(); ++it)
    if ((*it) != NULL)(*it)->fill(value);
}

std::string BlockVector::toString() const
{
  SiconosVectorException::selfThrow("BlockVector::toString, not yet implemented.");
  return "BlockVector";
}

/***************************** OPERATORS ******************************/

double BlockVector::getValue(unsigned int pos)
{
  unsigned int blockNum = 0;

  while (pos >= tabIndex[blockNum] && blockNum < tabIndex.size())
    blockNum ++;

  unsigned int relativePos = pos;

  if (blockNum != 0)
    relativePos -= tabIndex[blockNum - 1];

  return (*vect[blockNum])(relativePos);
}

void BlockVector::setValue(unsigned int pos, double value)
{
  unsigned int blockNum = 0;

  while (pos >= tabIndex[blockNum] && blockNum < tabIndex.size())
    blockNum ++;

  unsigned int relativePos = pos;

  if (blockNum != 0)
    relativePos -= tabIndex[blockNum - 1];

  (*vect[blockNum])(relativePos) = value;
}

double& BlockVector::operator()(unsigned int pos)
{
  unsigned int blockNum = 0;

  while (pos >= tabIndex[blockNum] && blockNum < tabIndex.size())
    blockNum ++;

  unsigned int relativePos = pos;

  if (blockNum != 0)
    relativePos -= tabIndex[blockNum - 1];

  return (*vect[blockNum])(relativePos);
}

double BlockVector::operator()(unsigned int pos) const
{
  unsigned int blockNum = 0;

  while (pos >= tabIndex[blockNum] && blockNum < tabIndex.size())
    blockNum ++;
  unsigned int relativePos = pos;

  if (blockNum != 0)
    relativePos -= tabIndex[blockNum - 1];

  return (*vect[blockNum])(relativePos);
}

SiconosVector* BlockVector::operator [](unsigned int pos)
{
  return  vect[pos];
}

const SiconosVector* BlockVector::operator [](unsigned int pos) const
{
  return  vect[pos];
}

BlockVector& BlockVector::operator = (const SiconosVector& m)
{
  if (&m == this) return *this;
  if (!check())
    SiconosVectorException::selfThrow("BlockVector:operator = m, one of the block is a NULL pointer.");

  if (m.size() != size())
    SiconosVectorException::selfThrow("BlockVector:operator = m, inconsistent size between this and m.");

  if (m.isBlock())
  {
    double numberOfBlocks = tabIndex.size();
    for (unsigned int i = 0; i < numberOfBlocks; ++i)
    {
      // Call = of SimpleVector
      *(vect[i]) = *(m[i]) ;
    }
  }
  else // if m is a Simple, use of () operator
  {
    for (unsigned int i = 0; i < size(); ++i)
      (*this)(i) = m(i);
  }

  return *this;
}

BlockVector& BlockVector::operator = (const BlockVector& m)
{
  if (&m == this) return *this;
  if (!check())
    SiconosVectorException::selfThrow("BlockVector:operator = m, one of the block is a NULL pointer.");

  if (m.size() != size())
    SiconosVectorException::selfThrow("BlockVector:operator = m, inconsistent size between this and m.");

  double numberOfBlocks = tabIndex.size();
  for (unsigned int i = 0; i < numberOfBlocks; ++i)
  {
    // Call = of SimpleVector
    *(vect[i]) = *(m[i]) ;
  }
  return *this;
}

BlockVector& BlockVector::operator =(const DenseVect&)
{
  SiconosVectorException::selfThrow("BlockVector:operator = DenseVect - Not implemented.");
  return *this;
}

BlockVector& BlockVector::operator =(const SparseVect&)
{
  SiconosVectorException::selfThrow("BlockVector:operator = SparseVect - Not implemented.");
  return *this;
}

BlockVector& BlockVector::operator += (const SiconosVector& m)
{
  if (m.size() != size())
    SiconosVectorException::selfThrow("BlockVector:operator += m, inconsistent size between this and m.");
  if (m.isBlock())
  {
    unsigned int n = getNumberOfBlocks();
    for (unsigned int i = 0; i < n; ++i)
      *vect[i] += *(m[i]);
  }
  else // if m is a Simple, use of () operator
  {
    for (unsigned int i = 0; i < size(); i++)
      (*this)(i) += m(i);
  }
  return *this;
}

BlockVector& BlockVector::operator -= (const SiconosVector& m)
{
  if (m.size() != size())
    SiconosVectorException::selfThrow("BlockVector:operator -= m, inconsistent size between this and m.");
  if (m.isBlock())
  {
    unsigned int n = getNumberOfBlocks();
    for (unsigned int i = 0; i < n; ++i)
      *vect[i] -= *(m[i]);
  }
  else // if m is a Simple, use of () operator
  {
    for (unsigned int i = 0; i < size(); i++)
      (*this)(i) -= m(i);
  }
  return *this;
}

BlockVector& BlockVector::operator *= (double m)
{
  BlockVectIterator it;
  for (it = vect.begin(); it != vect.end(); ++it)
    (**it) *= m;

  return *this;
}

BlockVector& BlockVector::operator *= (int m)
{
  BlockVectIterator it;
  for (it = vect.begin(); it != vect.end(); ++it)
    (**it) *= m;

  return *this;
}

BlockVector& BlockVector::operator /= (double m)
{
  BlockVectIterator it;
  for (it = vect.begin(); it != vect.end(); ++it)
    (**it) /= m;

  return *this;
}

BlockVector& BlockVector::operator /= (int m)
{
  BlockVectIterator it;
  for (it = vect.begin(); it != vect.end(); ++it)
    (**it) /= m;

  return *this;
}

void BlockVector::add(const  SiconosVector& v)
{
  if (v.isBlock())
    SiconosVectorException::selfThrow("BlockVector::add(v) : v is a block vector and cannot be added into another block vector.");

  unsigned int tabVal = v.size();
  if (tabIndex.size() != 0)
    tabVal += size();
  vect.push_back(new SimpleVector(v)); // Copy
  isBlockAllocatedIn.push_back(true);
  tabIndex.push_back(tabVal);
}

void BlockVector::addPtr(SiconosVector* v)
{
  if (v->isBlock())
    SiconosVectorException::selfThrow("BlockVector::addPtr(v) : v is a block vector and cannot be added into another block vector.");

  unsigned int tabVal = v->size();
  if (tabIndex.size() != 0)
    tabVal += size();
  vect.push_back(v);
  isBlockAllocatedIn.push_back(false);
  tabIndex.push_back(tabVal);
}
