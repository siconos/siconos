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
#include <vector>

/***************************** CONSTRUCTORS ****************************/
// Default
BlockVector::BlockVector(): SiconosVector(true)
{}

BlockVector::BlockVector(const std::string & file, bool ascii): SiconosVector(true)
{
  SiconosVectorException::selfThrow(" BlockVector::constructor from a file : read BlockVector is not implemented");
}

BlockVector::BlockVector(const BlockVector &v): SiconosVector(true)
{
  tabIndex = v.getTabIndex();
  ConstBlockVectIterator it;
  for (it = v.begin(); it != v.end(); ++it)
  {
    if (!(*it)->isBlock())  // Call copy-constructor of SimpleVector
      vect.push_back(new SimpleVector(**it)) ;
    else
      vect.push_back(new BlockVector(**it)) ;
    isBlockAllocatedIn.push_back(true);
  }
  sizeV = v.size();
}

BlockVector::BlockVector(const SiconosVector &v): SiconosVector(true)
{
  if (v.isBlock())
  {
    unsigned int numberOfBlocks = v.getNumberOfBlocks();;
    tabIndex.resize(numberOfBlocks);
    for (unsigned int i = 0; i < numberOfBlocks; ++i)
    {
      if (!(v[i]->isBlock()))   // Call copy-constructor of SimpleVector
        vect.push_back(new SimpleVector(*v[i])) ;
      else
        vect.push_back(new BlockVector(*v[i])) ;
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
  sizeV = v.size();
}

BlockVector::BlockVector(SiconosVector* v1, SiconosVector* v2): SiconosVector(true)
{
  // Add the two vectors in the container
  // NO COPY !!
  if ((v1 == NULL) && (v2 == NULL))
    SiconosVectorException::selfThrow("BlockVector:constructor(SimpleVector*,SimpleVector*), both vectors are NULL.");

  if (v1 != NULL)
  {
    vect.push_back(v1);
    tabIndex.push_back(v1->size());
    isBlockAllocatedIn.push_back(false);
    sizeV += v1->size();
  }
  else
    // If first parameter is a NULL pointer, then set this(1) to a SimpleVector of the same size as v2, and equal to 0.
  {
    // This case is usefull to set xDot in LagrangianDS.
    vect.push_back(new SimpleVector(v2->size()));
    tabIndex.push_back(v2->size());
    isBlockAllocatedIn.push_back(true);
    sizeV += v2->size();
  }
  if (v2 != NULL)
  {
    vect.push_back(v2);
    tabIndex.push_back(v2->size());
    tabIndex[1] += tabIndex[0];
    isBlockAllocatedIn.push_back(false);
    sizeV += v2->size();
  }
  else // If second parameter is a NULL pointer, then set this(2) to a SimpleVector of the same size as v1, and equal to 0.
  {
    // This case is usefull to set xDot in LagrangianDS.
    vect.push_back(new SimpleVector(v1->size()));
    tabIndex.push_back(v1->size());
    tabIndex[1] += tabIndex[0];
    isBlockAllocatedIn.push_back(true);
    sizeV += v1->size();
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
  sizeV = dim * numberOfBlocks;
}

BlockVector::BlockVector(unsigned int numberOfBlocks, unsigned int numberOfSubBlocks, BlocksVect blocks) : SiconosVector(true)
{
  unsigned int dim = numberOfSubBlocks * numberOfBlocks;
  if (dim != blocks.size())
    SiconosVectorException::selfThrow("BlockVector:constructor(unsigned int numberOfBlocks, unsigned int numberOfSubBlocks, std::vector<SiconosVector*> blocks), unconsistent dimension between blocks and numberOf...");

  vect.resize(numberOfBlocks, NULL);
  BlockVectIterator it, it2 = blocks.begin();
  for (it = vect.begin(); it != vect.end(); it++)
  {
    (*it) = new BlockVector() ;
    for (unsigned int i = 0; i < numberOfSubBlocks; ++i)
    {
      if ((*it2) == NULL)
        SiconosVectorException::selfThrow("BlockVector::constructor(int, int, vector<SiconosVector*>), one of the input vectors is NULL.");
      (*it)->addPtr(*it2);
      it2++; // Step to next element in blocks.
    }
    isBlockAllocatedIn.push_back(true);
    sizeV += (*it)->size();
    tabIndex.push_back(sizeV);
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
const unsigned int BlockVector::getNum() const
{
  SiconosVectorException::selfThrow("BlockVector::getNum of a block is forbidden.");
  return 1;
}

const DenseVect BlockVector::getDense(unsigned int i) const
{
  if (vect[i]->isBlock())
    SiconosVectorException::selfThrow("BlockVector::getDense(unsigned int num) : the vector[num] is a Block vector");

  if (vect[i]->getNum() != 1)
    SiconosVectorException::selfThrow("BlockVector::getDense(unsigned int num) : the vector[num] is not a Dense vector");
  return (vect[i])->getDense();
}

const SparseVect BlockVector::getSparse(unsigned int i)const
{
  if (vect[i]->isBlock())
    SiconosVectorException::selfThrow("BlockVector::getSparse(unsigned int num) : the vector[num] is a Block vector");

  if (vect[i]->getNum() != 4)
    SiconosVectorException::selfThrow("BlockVector::getSparse(unsigned int num) : the vector[num] is not a Sparse vector");
  return (vect[i])->getSparse();
}

DenseVect* BlockVector::getDensePtr(unsigned int i) const
{
  if (vect[i]->isBlock())
    SiconosVectorException::selfThrow("BlockVector::getDensePtr(unsigned int num) : the vector[num] is a Block vector");

  if (vect[i]->getNum() != 1)
    SiconosVectorException::selfThrow("BlockVector::getDensePtr(unsigned int num) : the vector[num] is not a Dense vector");
  return (vect[i])->getDensePtr();
}

SparseVect* BlockVector::getSparsePtr(unsigned int i) const
{
  if (vect[i]->isBlock())
    SiconosVectorException::selfThrow("BlockVector::getSparsePtr(unsigned int num) : the vector[num] is a Block vector");

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


const double BlockVector::norm2() const
{
  double d = 0;
  ConstBlockVectIterator it;
  for (it = vect.begin(); it != vect.end(); ++it)
  {
    if ((*it) != NULL)
      d += pow((*it)->norm2(), 2);
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

  if (vect[num]->isBlock())
    SiconosVectorException::selfThrow("BlockVector::getVector(num), vector[num] is a Block. Use getVectorPtr()");

  return *(vect[num]);
}

SiconosVector * BlockVector::getVectorPtr(unsigned int num)
{
  return vect[num];
}

void BlockVector::setVector(unsigned int num, const SiconosVector& newV)
{
  if (vect[num] == NULL)
    SiconosVectorException::selfThrow("BlockVector::setVector(num,v), this[num] == NULL pointer.");

  if (newV.size() != (vect[num])->size())
    SiconosVectorException::selfThrow("BlockVector::setVector(num,v), this[num] and v have unconsistent sizes.");

  *vect[num] = newV ;
}

void BlockVector::setVectorPtr(unsigned int num, SiconosVector* newV)
{
  if (newV->size() != (vect[num])->size())
    SiconosVectorException::selfThrow("BlockVector::setVectorPtr(num,v), this[num] and v have unconsistent sizes.");

  if (isBlockAllocatedIn[num]) delete vect[num];
  vect[num] = newV;
  isBlockAllocatedIn[num] = false;
}

void BlockVector::fill(double value)
{
  BlockVectIterator it;
  for (it = vect.begin(); it != vect.end(); ++it)
    if ((*it) != NULL)(*it)->fill(value);
}

const std::string BlockVector::toString() const
{
  SiconosVectorException::selfThrow("BlockVector::toString, not yet implemented.");
  return "BlockVector";
}

/***************************** OPERATORS ******************************/

const double BlockVector::getValue(unsigned int pos)
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

const double BlockVector::operator()(unsigned int pos) const
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

  if (m.size() != sizeV)
    SiconosVectorException::selfThrow("BlockVector:operator = m, inconsistent size between this and m.");

  if (m.isBlock())
  {
    unsigned int i = 0;
    BlockVectIterator it1;
    for (it1 = vect.begin(); it1 != vect.end(); ++it1)
      **it1 = *(m[i++]);
  }
  //       double numberOfBlocks = tabIndex.size();
  //       for (unsigned int i=0;i < numberOfBlocks; ++i)
  //  {
  //    // Call = of SimpleVector
  //    *(vect[i]) = *(m[i]) ;
  //  }
  else // if m is a Simple, use of () operator
  {
    unsigned int pos = 0;
    BlockVectIterator it1;
    for (it1 = vect.begin(); it1 != vect.end(); ++it1)
    {
      m.getBlock(pos, **it1);
      pos += (*it1)->size();
    }
  }

  //for(unsigned int i = 0; i<sizeV; ++i)
  //(*this)(i) = m(i);
  return *this;
}

BlockVector& BlockVector::operator = (const BlockVector& m)
{
  if (&m == this) return *this;

  if (m.size() != sizeV)
    SiconosVectorException::selfThrow("BlockVector:operator = m, inconsistent size between this and m.");
  BlockVectIterator it1;
  unsigned int i = 0;
  for (it1 = vect.begin(); it1 != vect.end(); ++it1)
    **it1 = *(m[i++]);
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
  if (m.size() != sizeV)
    SiconosVectorException::selfThrow("BlockVector:operator += m, inconsistent size between this and m.");
  if (m.isBlock())
  {
    unsigned int n = getNumberOfBlocks();
    for (unsigned int i = 0; i < n; ++i)
      *vect[i] += *(m[i]);
  }
  else // if m is a Simple, use of () operator
  {
    for (unsigned int i = 0; i < sizeV; i++)
      (*this)(i) += m(i);
  }
  return *this;
}

BlockVector& BlockVector::operator -= (const SiconosVector& m)
{
  if (m.size() != sizeV)
    SiconosVectorException::selfThrow("BlockVector:operator -= m, inconsistent size between this and m.");
  if (m.isBlock())
  {
    unsigned int n = getNumberOfBlocks();
    for (unsigned int i = 0; i < n; ++i)
      *vect[i] -= *(m[i]);
  }
  else // if m is a Simple, use of () operator
  {
    for (unsigned int i = 0; i < sizeV; i++)
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

  unsigned int tabVal = v.size();
  if (tabIndex.size() != 0)
    tabVal += sizeV;
  if (!v.isBlock())
    vect.push_back(new SimpleVector(v)); // Copy
  else
    vect.push_back(new BlockVector(v)); // Copy

  isBlockAllocatedIn.push_back(true);
  tabIndex.push_back(tabVal);
  sizeV += v.size();
}

void BlockVector::addPtr(SiconosVector* v)
{
  if (v == NULL)
    SiconosVectorException::selfThrow("BlockVector:addPtr(v), v is a NULL vector.");

  unsigned int tabVal = v->size();
  if (tabIndex.size() != 0)
    tabVal += sizeV;
  vect.push_back(v);
  isBlockAllocatedIn.push_back(false);
  tabIndex.push_back(tabVal);
  sizeV += v->size();
}
