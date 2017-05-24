/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#include "BlockVector.hpp"

#include <boost/numeric/ublas/vector_proxy.hpp>  // for project
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <vector>

#include "SiconosVector.hpp"

#include "SiconosAlgebra.hpp"
#include "Tools.hpp"


//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"


// =================================================
//                CONSTRUCTORS
// =================================================
BlockVector::BlockVector()
{
  _sizeV = 0;
  _tabIndex.reset(new Index());
}

BlockVector::BlockVector(const BlockVector &v)
{
  _sizeV = 0;
  unsigned int nbBlocks = v.numberOfBlocks();
  _tabIndex.reset(new Index());
  _tabIndex->reserve(nbBlocks);
  _vect.reserve(nbBlocks);
  VectorOfVectors::const_iterator it;
  for(it = v.begin(); it != v.end(); ++it)
  {
    _vect.push_back(std11::shared_ptr<SiconosVector>(new SiconosVector(**it))) ;
    _sizeV += (*it)->size();
    _tabIndex->push_back(_sizeV);
  }
}


BlockVector::BlockVector(SP::SiconosVector v1, SP::SiconosVector v2)
{
  _sizeV = 0;

  // Insert the two vectors in the container
  // NO COPY !!
  if(! v1  && ! v2)
    SiconosVectorException::selfThrow("BlockVector:constructor(SiconosVector*,SiconosVector*), both vectors are NULL.");

  _tabIndex.reset(new Index());

  _tabIndex->reserve(2);
  _vect.reserve(2);

  if(v1)
  {
    _vect.push_back(v1);
    _sizeV = v1->size();
    _tabIndex->push_back(_sizeV);

  }
  else
    // If first parameter is a NULL pointer, then set this(1) to a SiconosVector of the same size as v2, and equal to 0.
  {
    // This case is usefull to set xDot in LagrangianDS.
    _sizeV = v2->size();

    _vect.push_back(std11::shared_ptr<SiconosVector>(new SiconosVector(_sizeV)));
    _tabIndex->push_back(_sizeV);

  }
  if(v2)
  {
    _vect.push_back(v2);
    _sizeV += v2->size();
    _tabIndex->push_back(_sizeV);

  }
  else // If second parameter is a NULL pointer, then set this(2) to a SiconosVector of the same size as v1, and equal to 0.
  {
    // This case is usefull to set xDot in LagrangianDS.

    _vect.push_back(std11::shared_ptr<SiconosVector>(new SiconosVector(v1->size())));
    _sizeV += v1->size();
    _tabIndex->push_back(_sizeV);
  }
}

BlockVector::BlockVector(unsigned int numberOfBlocks, unsigned int dim)
{
  _sizeV = 0;

  _tabIndex.reset(new Index());
  _tabIndex->reserve(numberOfBlocks);
  _vect.reserve(numberOfBlocks);
  for(unsigned int i = 0; i < numberOfBlocks; ++i)
  {
    _vect.push_back(std11::shared_ptr<SiconosVector>(new SiconosVector(dim)));
    _tabIndex->push_back(dim * (i + 1));
  }
  _sizeV = dim * numberOfBlocks;
}

BlockVector::BlockVector(unsigned int numberOfBlocks)
{
  _sizeV = 0;
  _tabIndex.reset(new Index());
  _tabIndex->resize(numberOfBlocks);
  _vect.resize(numberOfBlocks);
}

BlockVector::~BlockVector()
{}

// ===========================
//      private method
// ===========================

void BlockVector::updateSizeV()
{
  _sizeV=0;
  VectorOfVectors::iterator it;
  for(it = _vect.begin(); it != _vect.end(); ++it)
  {
    if (*it)
      _sizeV += (*it)->size();
  }
}
void BlockVector::updateTabIndex()
{
  unsigned int cumulated_size=0;
  _tabIndex.reset(new Index());

  VectorOfVectors::iterator it;
  for(it = _vect.begin(); it != _vect.end(); ++it)
  {
    if (*it)
    {
      cumulated_size += (*it)->size();
    }
    _tabIndex->push_back(cumulated_size);
  }
}

// ===========================
//       fill vector
// ===========================

void BlockVector::zero()
{
  VectorOfVectors::iterator it;
  for(it = _vect.begin(); it != _vect.end(); ++it)
    (*it)->zero();
}

void BlockVector::fill(double value)
{
  VectorOfVectors::iterator it;
  for(it = _vect.begin(); it != _vect.end(); ++it)
    if((*it))(*it)->fill(value);
}

//=====================
// screen display
//=====================

void BlockVector::display() const
{
  VectorOfVectors::const_iterator it;
  std::cout << "=======> Block Vector Display (" << _tabIndex->size() << " block(s)): " << std::endl;
  for(it = _vect.begin(); it != _vect.end(); ++it)
  {
    DEBUG_EXPR(std::cout <<"(*it)" << (*it) << std::endl;);
    if (*it)
      (*it)->display();
    else
      std::cout << "(*it)-> NULL" <<std::endl;
  }
}

//=====================
// convert to a string
//=====================

std::string BlockVector::toString() const
{
  return ::toString(*this);
}

//=====================
// convert to an ostream
//=====================

std::ostream& operator<<(std::ostream& os, const BlockVector& bv)
{
  VectorOfVectors::const_iterator it;
  os << "[" << bv._vect.size() << "](";
  for(it = bv._vect.begin(); it != bv._vect.end(); ++it) {
    if (it != bv._vect.begin()) os << ",";
    if (*it) os << **it; else os << "(nil)";
  }
  os << ")";
  return os;
}

//=============================
// Elements access (get or set)
//=============================

double BlockVector::getValue(unsigned int pos) const
{
  unsigned int blockNum = 0;

  while(pos >= (*_tabIndex)[blockNum] && blockNum < _tabIndex->size())
    blockNum ++;

  unsigned int relativePos = pos;

  if(blockNum != 0)
    relativePos -= (*_tabIndex)[blockNum - 1];

  return (*_vect[blockNum])(relativePos);
}

void BlockVector::setValue(unsigned int pos, double value)
{
  unsigned int blockNum = 0;

  while(pos >= (*_tabIndex)[blockNum] && blockNum < _tabIndex->size())
    blockNum ++;

  unsigned int relativePos = pos;

  if(blockNum != 0)
    relativePos -= (*_tabIndex)[blockNum - 1];

  (*_vect[blockNum])(relativePos) = value;
}

double& BlockVector::operator()(unsigned int pos)
{
  unsigned int blockNum = 0;

  while(pos >= (*_tabIndex)[blockNum] && blockNum < _tabIndex->size())
    blockNum ++;

  unsigned int relativePos = pos;

  if(blockNum != 0)
    relativePos -= (*_tabIndex)[blockNum - 1];

  return (*_vect[blockNum])(relativePos);
}

double BlockVector::operator()(unsigned int pos) const
{
  unsigned int blockNum = 0;

  while(pos >= (*_tabIndex)[blockNum] && blockNum < _tabIndex->size())
    blockNum ++;
  unsigned int relativePos = pos;

  if(blockNum != 0)
    relativePos -= (*_tabIndex)[blockNum - 1];

  return (*_vect[blockNum])(relativePos);
}

//============================================
// Access (get or set) to blocks of elements
//============================================


void BlockVector::setVector(unsigned int pos, const SiconosVector& v)
{
  assert(pos < _vect.size() && "insertion out of vector size");
  if(! _vect[pos])
    SiconosVectorException::selfThrow("BlockVector::setVector(pos,v), this[pos] == NULL pointer.");

  // if(v.size() != (_vect[pos])->size())
  //   SiconosVectorException::selfThrow("BlockVector::setVector(pos,v), this[pos] and v have unconsistent sizes.");

  *_vect[pos] = v ;
}

void BlockVector::setVectorPtr(unsigned int pos, SP::SiconosVector v)
{
  assert(pos < _vect.size() && "insertion out of vector size");
  // if(v->size() != (_vect[pos])->size())
  //   SiconosVectorException::selfThrow("BlockVector::setVectorPtr(pos,v), this[pos] and v have unconsistent sizes.");
  _vect[pos] = v;
  updateSizeV();
  updateTabIndex();
}

SP::SiconosVector BlockVector::operator [](unsigned int pos)
{
  return  _vect[pos];
}

SPC::SiconosVector BlockVector::operator [](unsigned int pos) const
{
  return  _vect[pos];
}

unsigned int BlockVector::getNumVectorAtPos(unsigned int pos) const
{
  unsigned int blockNum = 0;

  while(pos >= (*_tabIndex)[blockNum] && blockNum < _tabIndex->size() - 1)
    blockNum ++;
  return blockNum;
}


BlockVector& BlockVector::operator = (const BlockVector& vIn)
{
  if(&vIn == this) return *this;
  else
  {
    if(isComparableTo(*this, vIn))  // if vIn and this are "block-consistent"
    {
      VectorOfVectors::iterator it1;
      VectorOfVectors::const_iterator it2 = vIn.begin();

      for(it1 = _vect.begin(); it1 != _vect.end(); ++it1)
      {
        (**it1) = (**it2);
        it2++;
      }
    }
    else
    {
      for(unsigned int i = 0; i < _sizeV; ++i)
        (*this)(i) = vIn(i);
    }
    return *this;
  }
}

BlockVector& BlockVector::operator = (const double* data)
{
  VectorOfVectors::iterator it1;
  unsigned indxPos = 0;

  for(it1 = _vect.begin(); it1 != _vect.end(); ++it1)
  {
    SiconosVector& v = **it1;
    v = &data[indxPos];
    indxPos += v.size();
  }
  return *this;
}


BlockVector& BlockVector::operator -= (const BlockVector& vIn)
{
  if(isComparableTo(*this, vIn))  // if vIn and this are "block-consistent"
  {
    unsigned int i = 0;
    VectorOfVectors::iterator it1;

    for(it1 = _vect.begin(); it1 != _vect.end(); ++it1)
      **it1 -= *(vIn[i++]);
  }
  else // use of a temporary SimpleVector... bad way, to be improved. But this case happens rarely ...
  {
    for(unsigned int i = 0; i < _sizeV; ++i)
      (*this)(i) -= vIn(i);
  }
  return *this;
}

BlockVector& BlockVector::operator -= (const SiconosVector& vIn)
{
  // Add a part of vIn (starting from index) to the current vector.
  // vIn must be a SimpleVector.

  // At the end of the present function, index is equal to index + the dim. of the added sub-vector.

  unsigned int dim = vIn.size(); // size of the block to be added.
  if(dim > _sizeV) SiconosVectorException::selfThrow("BlockVector::addSimple : invalid ranges");

  VectorOfVectors::const_iterator it;
  unsigned int numVIn = vIn.num();
  unsigned int currentSize, currentNum;
  unsigned int index = 0;
  for(it = _vect.begin(); it != _vect.end(); ++it)
  {
    currentSize = (*it)->size();
    currentNum = (*it)->num();
    if(numVIn != currentNum) SiconosVectorException::selfThrow("BlockVector::addSimple : inconsistent types.");
    if(numVIn == 1)
      noalias(*(*it)->dense()) -=  ublas::subrange(*vIn.dense(), index, index + currentSize) ;
    else
      noalias(*(*it)->sparse()) -=  ublas::subrange(*vIn.sparse(), index, index + currentSize) ;
    index += currentSize;
  }
  return *this;
}

BlockVector& BlockVector::operator += (const BlockVector& vIn)
{
  if(isComparableTo(*this, vIn))  // if vIn and this are "block-consistent"
  {
    unsigned int i = 0;
    VectorOfVectors::iterator it1;

    for(it1 = _vect.begin(); it1 != _vect.end(); ++it1)
      **it1 += *(vIn[i++]);
  }
  else // use of a temporary SimpleVector... bad way, to be improved. But this case happens rarely ...
  {
    for(unsigned int i = 0; i < _sizeV; ++i)
      (*this)(i) += vIn(i);
  }
  return *this;
}

BlockVector& BlockVector::operator += (const SiconosVector& vIn)
{
  // Add a part of vIn (starting from index) to the current vector.
  // vIn must be a SimpleVector.

  // At the end of the present function, index is equal to index + the dim. of the added sub-vector.

  unsigned int dim = vIn.size(); // size of the block to be added.
  if(dim > _sizeV) SiconosVectorException::selfThrow("BlockVector::addSimple : invalid ranges");

  VectorOfVectors::const_iterator it;
  unsigned int numVIn = vIn.num();
  unsigned int currentSize, currentNum;
  unsigned int index = 0;

  for(it = _vect.begin(); it != _vect.end(); ++it)
  {
    currentSize = (*it)->size();
    currentNum = (*it)->num();
    if(numVIn != currentNum) SiconosVectorException::selfThrow("BlockVector::addSimple : inconsistent types.");
    if(numVIn == 1)
      noalias(*(*it)->dense()) += ublas::subrange(*vIn.dense(), index, index + currentSize) ;
    else
      noalias(*(*it)->sparse()) += ublas::subrange(*vIn.sparse(), index, index + currentSize) ;
    index += currentSize;
  }
  return *this;
}

void BlockVector::insert(const  SiconosVector& v)
{
  _sizeV += v.size();

  _vect.push_back(std11::shared_ptr<SiconosVector>(new SiconosVector(v))); // Copy

  _tabIndex->push_back(_sizeV);
}

void BlockVector::insertPtr(SP::SiconosVector v)
{
  if(!v)
    SiconosVectorException::selfThrow("BlockVector:insertPtr(v), v is a NULL vector.");

  _sizeV += v->size();
  _vect.push_back(v);
  _tabIndex->push_back(_sizeV);
}

void BlockVector::setBlock(const SiconosVector& vIn, unsigned int sizeB, unsigned int startIn, unsigned int startOut)
{
  // Check dim ...
  unsigned int endOut = startOut + sizeB;

  assert(startIn < vIn.size());
  assert(startOut < size());
  assert((startIn + sizeB) <= vIn.size());
  assert(endOut <= size());

  // We look for the block of vOut that include index startOut
  unsigned int blockOutStart = 0;
  while(startOut >= (*_tabIndex)[blockOutStart] && blockOutStart < _tabIndex->size())
    blockOutStart++;
  // Relative position in the block blockOutStart.
  unsigned int posOut = startOut;
  if(blockOutStart != 0)
    posOut -= (*_tabIndex)[blockOutStart - 1];

  // We look for the block of vOut that include index endOut
  unsigned int blockOutEnd = blockOutStart;
  while(endOut > (*_tabIndex)[blockOutEnd] && blockOutEnd < _tabIndex->size())
    blockOutEnd ++;

  // => the block to be set runs from block number blockOutStart to block number blockOutEnd.

  if(blockOutEnd == blockOutStart)  //
  {
    vIn.toBlock(*_vect[blockOutStart], sizeB, startIn, posOut);
  }
  else // More that one block of vOut are concerned
  {

    // The current considered block ...
    SP::SiconosVector currentBlock = _vect[blockOutStart];

    // Size of the subBlock of vOut to be set.
    unsigned int subSizeB = currentBlock->size() - posOut;
    unsigned int posIn = startIn;

    // Set first sub-block (currentBlock) values, between index posOut and posOut+subSizeB,
    // with vIn values from posIn to posIn+subSizeB.
    vIn.toBlock(*currentBlock, subSizeB, posIn, posOut);

    // Other blocks, except number blockOutEnd.
    unsigned int currentBlockNum = blockOutStart + 1;
    while(currentBlockNum != blockOutEnd)
    {
      posIn += subSizeB;
      currentBlock = _vect[currentBlockNum];
      subSizeB = currentBlock->size();
      vIn.toBlock(*currentBlock, subSizeB, posIn, 0);
      currentBlockNum++;
    }
    // set last subBlock ...
    currentBlock = _vect[blockOutEnd];

    posIn += subSizeB;

    // Size of the considered sub-block
    subSizeB = endOut - (*_tabIndex)[blockOutEnd - 1];

    vIn.toBlock(*currentBlock, subSizeB, posIn, 0);
  }

}

bool BlockVector::isComparableTo(const BlockVector& v1, const BlockVector& v2)
{
  // return:
  //  - true if both are block but with blocks which are facing each other of the same size.
  //  - false in other cases
  //
  const Index& I1 = *v1.tabIndex();
  const Index& I2 = *v2.tabIndex();

  return (I1 == I2);

}

double BlockVector::norm2() const
{
  double d = 0;
  VectorOfVectors::const_iterator it;
  for(it = _vect.begin(); it != _vect.end(); ++it)
  {
    assert(*it);
    d += pow((*it)->norm2(), 2);
  }
  return sqrt(d);
}

double BlockVector::normInf() const
{
  double d = 0;
  VectorOfVectors::const_iterator it;
  for(it = _vect.begin(); it != _vect.end(); ++it)
  {
    assert(*it);
    d = fmax((*it)->normInf(), d);
  }
  return d;
}

BlockVector& BlockVector::operator =(const SiconosVector& vIn)
{
  setBlock(vIn, _sizeV, 0, 0);
  return *this;
}

BlockVector& BlockVector::operator *= (double s)
{
  VectorOfVectors::iterator it;
  for(it = begin(); it != end(); ++it)
    (**it) *= s;
  return *this;
}

BlockVector& BlockVector::operator /= (double s)
{
  VectorOfVectors::iterator it;
  for(it = begin(); it != end(); ++it)
    (**it) /= s;
  return *this;
}
