/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

#include "BlockVector.hpp"

#include <boost/numeric/ublas/vector_proxy.hpp>  // for project
#include <vector>

#include "SiconosVector.hpp"

#include "SiconosAlgebra.hpp"

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
  unsigned int nbBlocks = v.getNumberOfBlocks();
  _tabIndex.reset(new Index());
  _tabIndex->reserve(nbBlocks);
  vect.reserve(nbBlocks);
  VectorOfVectors::const_iterator it;
  for (it = v.begin(); it != v.end(); ++it)
  {
    vect.push_back(cpp11ns::shared_ptr<SiconosVector>(new SiconosVector(**it))) ;
    _sizeV += (*it)->size();
    _tabIndex->push_back(_sizeV);
  }
}


BlockVector::BlockVector(SP::SiconosVector v1, SP::SiconosVector v2)
{
  _sizeV = 0;

  // Insert the two vectors in the container
  // NO COPY !!
  if (! v1  && ! v2)
    SiconosVectorException::selfThrow("BlockVector:constructor(SiconosVector*,SiconosVector*), both vectors are NULL.");

  _tabIndex.reset(new Index());

  _tabIndex->reserve(2);
  vect.reserve(2);

  if (v1)
  {
    vect.push_back(v1);
    _sizeV = v1->size();
    _tabIndex->push_back(_sizeV);

  }
  else
    // If first parameter is a NULL pointer, then set this(1) to a SiconosVector of the same size as v2, and equal to 0.
  {
    // This case is usefull to set xDot in LagrangianDS.
    _sizeV = v2->size();

    vect.push_back(cpp11ns::shared_ptr<SiconosVector>(new SiconosVector(_sizeV)));
    _tabIndex->push_back(_sizeV);

  }
  if (v2)
  {
    vect.push_back(v2);
    _sizeV += v2->size();
    _tabIndex->push_back(_sizeV);

  }
  else // If second parameter is a NULL pointer, then set this(2) to a SiconosVector of the same size as v1, and equal to 0.
  {
    // This case is usefull to set xDot in LagrangianDS.

    vect.push_back(cpp11ns::shared_ptr<SiconosVector>(new SiconosVector(v1->size())));
    _sizeV += v1->size();
    _tabIndex->push_back(_sizeV);
  }
}

BlockVector::BlockVector(unsigned int numberOfBlocks, unsigned int dim)
{
  _sizeV = 0;

  _tabIndex.reset(new Index());
  _tabIndex->reserve(numberOfBlocks);
  vect.reserve(numberOfBlocks);
  for (unsigned int i = 0; i < numberOfBlocks; ++i)
  {
    vect.push_back(cpp11ns::shared_ptr<SiconosVector>(new SiconosVector(dim)));
    _tabIndex->push_back(dim * (i + 1));
  }
  _sizeV = dim * numberOfBlocks;
}


BlockVector::~BlockVector()
{}

// ===========================
//       fill vector
// ===========================

void BlockVector::zero()
{
  VectorOfVectors::iterator it;
  for (it = vect.begin(); it != vect.end(); ++it)
    (*it)->zero();
}

void BlockVector::fill(double value)
{
  VectorOfVectors::iterator it;
  for (it = vect.begin(); it != vect.end(); ++it)
    if ((*it))(*it)->fill(value);
}

//=====================
// screen display
//=====================

void BlockVector::display() const
{
  VectorOfVectors::const_iterator it;
  std::cout << "=======> Block Vector Display (" << _tabIndex->size() << " block(s)): " << std::endl;
  for (it = vect.begin(); it != vect.end(); ++it)
    (*it)->display();
}


//=============================
// Elements access (get or set)
//=============================

double BlockVector::getValue(unsigned int pos) const
{
  unsigned int blockNum = 0;

  while (pos >= (*_tabIndex)[blockNum] && blockNum < _tabIndex->size())
    blockNum ++;

  unsigned int relativePos = pos;

  if (blockNum != 0)
    relativePos -= (*_tabIndex)[blockNum - 1];

  return (*vect[blockNum])(relativePos);
}

void BlockVector::setValue(unsigned int pos, double value)
{
  unsigned int blockNum = 0;

  while (pos >= (*_tabIndex)[blockNum] && blockNum < _tabIndex->size())
    blockNum ++;

  unsigned int relativePos = pos;

  if (blockNum != 0)
    relativePos -= (*_tabIndex)[blockNum - 1];

  (*vect[blockNum])(relativePos) = value;
}

double& BlockVector::operator()(unsigned int pos)
{
  unsigned int blockNum = 0;

  while (pos >= (*_tabIndex)[blockNum] && blockNum < _tabIndex->size())
    blockNum ++;

  unsigned int relativePos = pos;

  if (blockNum != 0)
    relativePos -= (*_tabIndex)[blockNum - 1];

  return (*vect[blockNum])(relativePos);
}

double BlockVector::operator()(unsigned int pos) const
{
  unsigned int blockNum = 0;

  while (pos >= (*_tabIndex)[blockNum] && blockNum < _tabIndex->size())
    blockNum ++;
  unsigned int relativePos = pos;

  if (blockNum != 0)
    relativePos -= (*_tabIndex)[blockNum - 1];

  return (*vect[blockNum])(relativePos);
}

//============================================
// Access (get or set) to blocks of elements
//============================================


void BlockVector::setVector(unsigned int pos, const SiconosVector& newV)
{
  if (! vect[pos])
    SiconosVectorException::selfThrow("BlockVector::setVector(pos,v), this[pos] == NULL pointer.");

  if (newV.size() != (vect[pos])->size())
    SiconosVectorException::selfThrow("BlockVector::setVector(pos,v), this[pos] and v have unconsistent sizes.");

  *vect[pos] = newV ;
}

void BlockVector::setVectorPtr(unsigned int pos, SP::SiconosVector newV)
{
  if (newV->size() != (vect[pos])->size())
    SiconosVectorException::selfThrow("BlockVector::setVectorPtr(pos,v), this[pos] and v have unconsistent sizes.");

  vect[pos] = newV;
}

SP::SiconosVector BlockVector::operator [](unsigned int pos)
{
  return  vect[pos];
}

SPC::SiconosVector BlockVector::operator [](unsigned int pos) const
{
  return  vect[pos];
}

unsigned int BlockVector::getNumVectorAtPos(unsigned int pos) const
{
  unsigned int blockNum = 0;

  while (pos >= (*_tabIndex)[blockNum] && blockNum < _tabIndex->size() - 1)
    blockNum ++;
  return blockNum;
}


BlockVector& BlockVector::operator = (const BlockVector& vIn)
{
  if (&vIn == this) return *this;
  else
  {
    if (isComparableTo(*this, vIn)) // if vIn and this are "block-consistent"
    {
      VectorOfVectors::iterator it1;
      VectorOfVectors::const_iterator it2 = vIn.begin();

      for (it1 = vect.begin(); it1 != vect.end(); ++it1)
      {
        (**it1) = (**it2);
        it2++;
      }
    }
    else
    {
      for (unsigned int i = 0; i < _sizeV; ++i)
        (*this)(i) = vIn(i);
    }
    return *this;
  }
  return *this;
}

BlockVector& BlockVector::operator -= (const BlockVector& vIn)
{
  if (isComparableTo(*this, vIn)) // if vIn and this are "block-consistent"
  {
    unsigned int i = 0;
    VectorOfVectors::iterator it1;

    for (it1 = vect.begin(); it1 != vect.end(); ++it1)
      **it1 -= *(vIn[i++]);
  }
  else // use of a temporary SimpleVector... bad way, to be improved. But this case happens rarely ...
  {
    for (unsigned int i = 0; i < _sizeV; ++i)
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
  if (dim > _sizeV) SiconosVectorException::selfThrow("BlockVector::addSimple : invalid ranges");

  VectorOfVectors::const_iterator it;
  unsigned int numVIn = vIn.getNum();
  unsigned int currentSize, currentNum;
  unsigned int index = 0;
  for (it = vect.begin(); it != vect.end(); ++it)
  {
    currentSize = (*it)->size();
    currentNum = (*it)->getNum();
    if (numVIn != currentNum) SiconosVectorException::selfThrow("BlockVector::addSimple : inconsistent types.");
    if (numVIn == 1)
      noalias(*(*it)->dense()) -=  ublas::subrange(*vIn.dense(), index, index + currentSize) ;
    else
      noalias(*(*it)->sparse()) -=  ublas::subrange(*vIn.sparse(), index, index + currentSize) ;
    index += currentSize;
  }
  return *this;
}

BlockVector& BlockVector::operator += (const BlockVector& vIn)
{
  if (isComparableTo(*this, vIn)) // if vIn and this are "block-consistent"
  {
    unsigned int i = 0;
    VectorOfVectors::iterator it1;

    for (it1 = vect.begin(); it1 != vect.end(); ++it1)
      **it1 += *(vIn[i++]);
  }
  else // use of a temporary SimpleVector... bad way, to be improved. But this case happens rarely ...
  {
    for (unsigned int i = 0; i < _sizeV; ++i)
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
  if (dim > _sizeV) SiconosVectorException::selfThrow("BlockVector::addSimple : invalid ranges");

  VectorOfVectors::const_iterator it;
  unsigned int numVIn = vIn.getNum();
  unsigned int currentSize, currentNum;
  unsigned int index = 0;

  for (it = vect.begin(); it != vect.end(); ++it)
  {
    currentSize = (*it)->size();
    currentNum = (*it)->getNum();
    if (numVIn != currentNum) SiconosVectorException::selfThrow("BlockVector::addSimple : inconsistent types.");
    if (numVIn == 1)
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

  vect.push_back(cpp11ns::shared_ptr<SiconosVector>(new SiconosVector(v))); // Copy

  _tabIndex->push_back(_sizeV);
}

void BlockVector::insertPtr(SP::SiconosVector v)
{
  if (!v)
    SiconosVectorException::selfThrow("BlockVector:insertPtr(v), v is a NULL vector.");

  _sizeV += v->size();
  vect.push_back(v);
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
  while (startOut >= (*_tabIndex)[blockOutStart] && blockOutStart < _tabIndex->size())
    blockOutStart++;
  // Relative position in the block blockOutStart.
  unsigned int posOut = startOut;
  if (blockOutStart != 0)
    posOut -= (*_tabIndex)[blockOutStart - 1];

  // We look for the block of vOut that include index endOut
  unsigned int blockOutEnd = blockOutStart;
  while (endOut > (*_tabIndex)[blockOutEnd] && blockOutEnd < _tabIndex->size())
    blockOutEnd ++;

  // => the block to be set runs from block number blockOutStart to block number blockOutEnd.

  if (blockOutEnd == blockOutStart) //
  {
    vIn.toBlock(*vect[blockOutStart], sizeB, startIn, posOut);
  }
  else // More that one block of vOut are concerned
  {

    // The current considered block ...
    SP::SiconosVector currentBlock = vect[blockOutStart];

    // Size of the subBlock of vOut to be set.
    unsigned int subSizeB = currentBlock->size() - posOut;
    unsigned int posIn = startIn;

    // Set first sub-block (currentBlock) values, between index posOut and posOut+subSizeB,
    // with vIn values from posIn to posIn+subSizeB.
    vIn.toBlock(*currentBlock, subSizeB, posIn, posOut);

    // Other blocks, except number blockOutEnd.
    unsigned int currentBlockNum = blockOutStart + 1;
    while (currentBlockNum != blockOutEnd)
    {
      posIn += subSizeB;
      currentBlock = vect[currentBlockNum];
      subSizeB = currentBlock->size();
      vIn.toBlock(*currentBlock, subSizeB, posIn, 0);
      currentBlockNum++;
    }
    // set last subBlock ...
    currentBlock = vect[blockOutEnd];

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
  for (it = vect.begin(); it != vect.end(); ++it)
  {
    assert(*it);
    d += pow((*it)->norm2(), 2);
  }
  return sqrt(d);
}

BlockVector& BlockVector::operator =(const SiconosVector& vIn)
{
  setBlock(vIn, _sizeV, 0, 0);
  return *this;
}

BlockVector& BlockVector::operator *= (double s)
{
  VectorOfVectors::iterator it;
  for (it = begin(); it != end(); ++it)
    (**it) *= s;
  return *this;
}

BlockVector& BlockVector::operator /= (double s)
{
  VectorOfVectors::iterator it;
  for (it = begin(); it != end(); ++it)
    (**it) /= s;
  return *this;
}
