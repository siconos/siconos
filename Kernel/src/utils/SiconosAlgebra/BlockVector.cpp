/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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

#include <boost/numeric/ublas/vector_proxy.hpp>  // for project
#include <vector>

#include "BlockVector.h"
#include "SimpleVector.h"


// =================================================
//                CONSTRUCTORS
// =================================================
BlockVector::BlockVector(): SiconosVector(0), tabIndex(NULL)
#ifndef WithSmartPtr
  , isBlockAllocatedIn(NULL)
#endif
{
  tabIndex = new Index();

#ifndef WithSmartPtr
  isBlockAllocatedIn = new std::vector<bool> ();
#endif

}

BlockVector::BlockVector(const std::string & file, bool ascii): SiconosVector(0), tabIndex(NULL)
#ifndef WithSmartPtr
  , isBlockAllocatedIn(NULL)
#endif
{
  SiconosVectorException::selfThrow(" BlockVector::constructor from a file : read BlockVector is not implemented");
}

BlockVector::BlockVector(const BlockVector &v): SiconosVector(0), tabIndex(NULL)
#ifndef WithSmartPtr
  , isBlockAllocatedIn(NULL)
#endif
{
  unsigned int nbBlocks = v.getNumberOfBlocks();
  tabIndex = new Index();
  tabIndex->reserve(nbBlocks);

#ifndef WithSmartPtr
  isBlockAllocatedIn = new std::vector<bool> ();
  isBlockAllocatedIn->reserve(nbBlocks);
#endif

  vect.reserve(nbBlocks);
  ConstBlockVectIterator it;
  for (it = v.begin(); it != v.end(); ++it)
  {

#ifndef WithSmartPtr
    if (!(*it)->isBlock())  // Call copy-constructor of SimpleVector
      vect.push_back(new SimpleVector(**it)) ;
    else
      vect.push_back(new BlockVector(**it)) ;
    isBlockAllocatedIn->push_back(true);
#else
    if (!(*it)->isBlock())  // Call copy-constructor of SimpleVector
      vect.push_back(boost::shared_ptr<SimpleVector>(new SimpleVector(**it))) ;
    else
      vect.push_back(boost::shared_ptr<BlockVector>(new BlockVector(**it))) ;
#endif

    sizeV += (*it)->size();
    tabIndex->push_back(sizeV);
  }
}

BlockVector::BlockVector(const SiconosVector &v): SiconosVector(0), tabIndex(NULL)
#ifndef WithSmartPtr
  , isBlockAllocatedIn(NULL)
#endif
{
  unsigned int nbBlocks = v.getNumberOfBlocks();
  tabIndex = new Index();

#ifndef WithSmartPtr
  isBlockAllocatedIn = new std::vector<bool> ();
#endif

  tabIndex->reserve(nbBlocks);
  vect.reserve(nbBlocks);

#ifndef WithSmartPtr
  isBlockAllocatedIn->reserve(nbBlocks);
#endif

  if (v.isBlock())
  {
    ConstBlockVectIterator it;
    for (it = v.begin(); it != v.end(); ++it)
    {

#ifndef WithSmartPtr
      if (!(*it)->isBlock())  // Call copy-constructor of SimpleVector
        vect.push_back(new SimpleVector(**it)) ;
      else
        vect.push_back(new BlockVector(**it)) ;

      isBlockAllocatedIn->push_back(true);
#else
      if (!(*it)->isBlock())  // Call copy-constructor of SimpleVector
        vect.push_back(boost::shared_ptr<SimpleVector>(new SimpleVector(**it))) ;
      else
        vect.push_back(boost::shared_ptr<BlockVector>(new BlockVector(**it))) ;
#endif

      sizeV += (*it)->size();
      tabIndex->push_back(sizeV);
    }
  }
  else
  {
    // Call copy-constructor of SimpleVector
#ifndef WithSmartPtr
    vect.push_back(new SimpleVector(v));
    isBlockAllocatedIn->push_back(true);
#else
    vect.push_back(boost::shared_ptr<SimpleVector>(new SimpleVector(v)));

#endif

    sizeV = v.size();
    tabIndex->push_back(sizeV);
  }
}

BlockVector::BlockVector(SiconosVectorSPtr v1, SiconosVectorSPtr v2): SiconosVector(0), tabIndex(NULL)
#ifndef WithSmartPtr
  , isBlockAllocatedIn(NULL)
#endif
{
  // Insert the two vectors in the container
  // NO COPY !!
  if (! v1  && ! v2)
    SiconosVectorException::selfThrow("BlockVector:constructor(SimpleVector*,SimpleVector*), both vectors are NULL.");

  tabIndex = new Index();

#ifndef WithSmartPtr
  isBlockAllocatedIn = new std::vector<bool> ();
#endif

  tabIndex->reserve(2);
  vect.reserve(2);

#ifndef WithSmartPtr
  isBlockAllocatedIn->reserve(2);
#endif

  if (v1)
  {
    vect.push_back(v1);
    sizeV = v1->size();
    tabIndex->push_back(sizeV);

#ifndef WithSmartPtr
    isBlockAllocatedIn->push_back(false);
#endif

  }
  else
    // If first parameter is a NULL pointer, then set this(1) to a SimpleVector of the same size as v2, and equal to 0.
  {
    // This case is usefull to set xDot in LagrangianDS.
    sizeV = v2->size();

#ifndef WithSmartPtr
    vect.push_back(new SimpleVector(sizeV));
#else
    vect.push_back(boost::shared_ptr<SimpleVector>(new SimpleVector(sizeV)));
#endif

    tabIndex->push_back(sizeV);

#ifndef WithSmartPtr
    isBlockAllocatedIn->push_back(true);
#endif

  }
  if (v2)
  {
    vect.push_back(v2);
    sizeV += v2->size();
    tabIndex->push_back(sizeV);

#ifndef WithSmartPtr
    isBlockAllocatedIn->push_back(false);
#endif

  }
  else // If second parameter is a NULL pointer, then set this(2) to a SimpleVector of the same size as v1, and equal to 0.
  {
    // This case is usefull to set xDot in LagrangianDS.

#ifndef WithSmartPtr
    vect.push_back(new SimpleVector(v1->size()));
#else
    vect.push_back(boost::shared_ptr<SimpleVector>(new SimpleVector(v1->size())));
#endif

    sizeV += v1->size();
    tabIndex->push_back(sizeV);

#ifndef WithSmartPtr
    isBlockAllocatedIn->push_back(true);
#endif

  }
}

BlockVector::BlockVector(unsigned int numberOfBlocks, unsigned int dim): SiconosVector(0), tabIndex(NULL)
#ifndef WithSmartPtr
  , isBlockAllocatedIn(NULL)
#endif
{
  tabIndex = new Index();

#ifndef WithSmartPtr
  isBlockAllocatedIn = new std::vector<bool> ();
#endif

  tabIndex->reserve(numberOfBlocks);
  vect.reserve(numberOfBlocks);

#ifndef WithSmartPtr
  isBlockAllocatedIn->reserve(numberOfBlocks);
#endif

  for (unsigned int i = 0; i < numberOfBlocks; ++i)
  {
#ifndef WithSmartPtr
    vect.push_back(new SimpleVector(dim));
#else
    vect.push_back(boost::shared_ptr<SimpleVector>(new SimpleVector(dim)));
#endif

    tabIndex->push_back(dim * (i + 1));

#ifndef WithSmartPtr
    isBlockAllocatedIn->push_back(true);
#endif

  }
  sizeV = dim * numberOfBlocks;
}

BlockVector::~BlockVector()
{

#ifndef WithSmartPtr
  purge(vect, *isBlockAllocatedIn);
#endif

  tabIndex->clear();
  delete tabIndex;

#ifndef WithSmartPtr
  isBlockAllocatedIn->clear();
  delete isBlockAllocatedIn;
#endif

}

// =================================================
//        get Ublas component (dense or sparse)
// =================================================

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

SparseVect* BlockVector::getSparsePtr(unsigned int i) const
{
  if (vect[i]->getNum() != 4)
    SiconosVectorException::selfThrow("BlockVector::getSparsePtr(unsigned int num) : the vector[num] is not a Sparse vector");
  return (vect[i])->getSparsePtr();
}

double* BlockVector::getArray(unsigned int i) const
{
  if (vect[i]->isBlock())
    SiconosVectorException::selfThrow("BlockVector::getArray(unsigned int num) : the vector[num] is a Block vector");
  return (vect[i])->getArray();
}

// ===========================
//       fill vector
// ===========================

void BlockVector::zero()
{
  BlockVectIterator it;
  for (it = vect.begin(); it != vect.end(); ++it)
    (*it)->zero();
}

void BlockVector::fill(double value)
{
  BlockVectIterator it;
  for (it = vect.begin(); it != vect.end(); ++it)
    if ((*it))(*it)->fill(value);
}

//=======================
// set vector dimension
//=======================

void BlockVector::resize(unsigned int, bool)
{
  SiconosVectorException::selfThrow("BlockVector::resize, not allowed for block vectors.");
}

//=======================
//       get norm
//=======================

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
    if ((*it))
      d += pow((*it)->norm2(), 2);

    else
      SiconosVectorException::selfThrow("BlockVector::norm, one of the blocks is equal to NULL pointer.");
  }
  return sqrt(d);
}

//=====================
// screen display
//=====================

void BlockVector::display() const
{
  ConstBlockVectIterator it;
  std::cout << "=======> Block Vector Display (" << tabIndex->size() << " block(s)): " << std::endl;
  for (it = vect.begin(); it != vect.end(); ++it)
    (*it)->display();
}

//============================
// Convert vector to a string
//============================

const std::string BlockVector::toString() const
{
  SiconosVectorException::selfThrow("BlockVector::toString, not yet implemented.");
  return "BlockVector";
}

//=============================
// Elements access (get or set)
//=============================

const double BlockVector::getValue(unsigned int pos) const
{
  unsigned int blockNum = 0;

  while (pos >= (*tabIndex)[blockNum] && blockNum < tabIndex->size())
    blockNum ++;

  unsigned int relativePos = pos;

  if (blockNum != 0)
    relativePos -= (*tabIndex)[blockNum - 1];

  return (*vect[blockNum])(relativePos);
}

void BlockVector::setValue(unsigned int pos, double value)
{
  unsigned int blockNum = 0;

  while (pos >= (*tabIndex)[blockNum] && blockNum < tabIndex->size())
    blockNum ++;

  unsigned int relativePos = pos;

  if (blockNum != 0)
    relativePos -= (*tabIndex)[blockNum - 1];

  (*vect[blockNum])(relativePos) = value;
}

double& BlockVector::operator()(unsigned int pos)
{
  unsigned int blockNum = 0;

  while (pos >= (*tabIndex)[blockNum] && blockNum < tabIndex->size())
    blockNum ++;

  unsigned int relativePos = pos;

  if (blockNum != 0)
    relativePos -= (*tabIndex)[blockNum - 1];

  return (*vect[blockNum])(relativePos);
}

const double BlockVector::operator()(unsigned int pos) const
{
  unsigned int blockNum = 0;

  while (pos >= (*tabIndex)[blockNum] && blockNum < tabIndex->size())
    blockNum ++;
  unsigned int relativePos = pos;

  if (blockNum != 0)
    relativePos -= (*tabIndex)[blockNum - 1];

  return (*vect[blockNum])(relativePos);
}

//============================================
// Access (get or set) to blocks of elements
//============================================

SimpleVector BlockVector::getVector(unsigned int pos) const
{

  if (! vect[pos])
    SiconosVectorException::selfThrow("BlockVector::getVector(pos), vector[pos] == NULL pointer.");

  if (vect[pos]->isBlock())
    SiconosVectorException::selfThrow("BlockVector::getVector(pos), vector[pos] is a Block. Use getVectorPtr()");

  return *(vect[pos]);
}

void BlockVector::setVector(unsigned int pos, const SiconosVector& newV)
{
  if (! vect[pos])
    SiconosVectorException::selfThrow("BlockVector::setVector(pos,v), this[pos] == NULL pointer.");

  if (newV.size() != (vect[pos])->size())
    SiconosVectorException::selfThrow("BlockVector::setVector(pos,v), this[pos] and v have unconsistent sizes.");

  *vect[pos] = newV ;
}

void BlockVector::setVectorPtr(unsigned int pos, SiconosVectorSPtr newV)
{
  if (newV->size() != (vect[pos])->size())
    SiconosVectorException::selfThrow("BlockVector::setVectorPtr(pos,v), this[pos] and v have unconsistent sizes.");

#ifndef WithSmartPtr
  if ((*isBlockAllocatedIn)[pos]) delete vect[pos];
  (*isBlockAllocatedIn)[pos] = false;
#endif

  vect[pos] = newV;
}

SiconosVectorSPtr BlockVector::operator [](unsigned int pos)
{
  return  vect[pos];
}

#ifndef WithSmartPtr
const SiconosVector * BlockVector::operator [](unsigned int pos) const
#else
SiconosVectorSPtrConst BlockVector::operator [](unsigned int pos) const
#endif
{
  return  vect[pos];
}

unsigned int BlockVector::getNumVectorAtPos(unsigned int pos) const
{
  unsigned int blockNum = 0;

  while (pos >= (*tabIndex)[blockNum] && blockNum < tabIndex->size() - 1)
    blockNum ++;
  return blockNum;
}

void BlockVector::addSimple(unsigned int& index, const SiconosVector& vIn)
{
  // Add a part of vIn (starting from index) to the current vector.
  // vIn must be a SimpleVector.

  // At the end of the present function, index is equal to index + the dim. of the added sub-vector.

  unsigned int dim = vIn.size() - index; // size of the block to be added.
  if (dim > sizeV) SiconosVectorException::selfThrow("BlockVector::addSimple : invalid ranges");

  ConstBlockVectIterator it;
  unsigned int numVIn = vIn.getNum();
  unsigned int currentSize, currentNum;

  for (it = vect.begin(); it != vect.end(); ++it)
  {
    if ((*it)->isBlock())
#ifndef WithSmartPtr
      (static_cast<BlockVector*>(*it))->addSimple(index, vIn);
#else
      (boost::static_pointer_cast<BlockVector>(*it))->addSimple(index, vIn);
#endif
    else // the current block is a SimpleVector
    {
      currentSize = (*it)->size();
      currentNum = (*it)->getNum();
      if (numVIn != currentNum) SiconosVectorException::selfThrow("BlockVector::addSimple : inconsistent types.");
      if (numVIn == 1)
        noalias(*(*it)->getDensePtr()) +=  ublas::subrange(*vIn.getDensePtr(), index, index + currentSize) ;
      else
        noalias(*(*it)->getSparsePtr()) +=  ublas::subrange(*vIn.getSparsePtr(), index, index + currentSize) ;
      index += currentSize;
    }
  }
}

void BlockVector::subSimple(unsigned int& index, const SiconosVector& vIn)
{
  // subtract a part of vIn (starting from index) to the current vector.
  // vIn must be a SimpleVector.

  // At the end of the present function, index is equal to index + the dim. of the added sub-vector.

  unsigned int dim = vIn.size() - index; // size of the block to be added.
  if (dim > sizeV) SiconosVectorException::selfThrow("BlockVector::addSimple : invalid ranges");

  ConstBlockVectIterator it;
  unsigned int numVIn = vIn.getNum();
  unsigned int currentSize, currentNum;

  for (it = vect.begin(); it != vect.end(); ++it)
  {
    if ((*it)->isBlock())
#ifndef WithSmartPtr
      (static_cast<BlockVector*>(*it))->subSimple(index, vIn);
#else
      (boost::static_pointer_cast<BlockVector>(*it))->subSimple(index, vIn);
#endif

    else // the current block is a SimpleVector
    {
      currentSize = (*it)->size();
      currentNum = (*it)->getNum();
      if (numVIn != currentNum) SiconosVectorException::selfThrow("BlockVector::addSimple : inconsistent types.");
      if (numVIn == 1)
        noalias(*(*it)->getDensePtr()) -=  ublas::subrange(*vIn.getDensePtr(), index, index + currentSize) ;
      else
        noalias(*(*it)->getSparsePtr()) -=  ublas::subrange(*vIn.getSparsePtr(), index, index + currentSize) ;
      index += currentSize;
    }
  }
}

//===============
//  Assignment
//===============
BlockVector& BlockVector::operator = (const SiconosVector& vIn)
{
  if (&vIn == this) return *this;

  if (vIn.size() != sizeV)
    SiconosVectorException::selfThrow("BlockVector:operator = vIn, inconsistent size between this and vIn.");

  if (vIn.isBlock())
  {
    if (isComparableTo(this, &vIn)) // if vIn and this are "block-consistent"
    {
      BlockVectIterator it1;
      ConstBlockVectIterator it2 = vIn.begin();

      for (it1 = vect.begin(); it1 != vect.end(); ++it1)
      {
        (**it1) = (**it2);
        it2++;
      }
    }
    else // use of a temporary SimpleVector... bad way, to be improved. But this case happens rarely ...
    {
      //    std::cout << "WARNING: BlockVector::operator = BlockVector, v=w, with v and w having blocks of different sizes. This operation may be time-consuming and inappropriate." << std::endl;
      //    SiconosVector * tmp = new SimpleVector(vIn);
      //    *this = *tmp;
      //    delete tmp;

      // component by component copy ...
      for (unsigned int i = 0; i < sizeV; ++i)
        (*this)(i) = vIn(i);

    }
  }
  else // if vIn is a SimpleVector
  {
#ifndef WithSmartPtr
    setBlock(&vIn, this, sizeV, 0, 0);
#else
    setBlock(boost::shared_ptr<const SiconosVector>(&vIn), shared_from_this(), sizeV, 0, 0);
#endif
    //       unsigned int pos = 0;
    //       BlockVectIterator it1;
    //       for(it1=vect.begin(); it1!=vect.end();++it1)
    //  {

    //    vIn.getBlock(pos,*it1);
    //    pos += (*it1)->size();
    //  }
  }
  return *this;
}

BlockVector& BlockVector::operator = (const BlockVector& vIn)
{
  if (&vIn == this) return *this;

  if (vIn.size() != sizeV)
    SiconosVectorException::selfThrow("BlockVector:operator = vIn, inconsistent size between this and vIn.");

  if (isComparableTo(this, &vIn)) // if vIn and this are "block-consistent"
  {
    BlockVectIterator it1;
    ConstBlockVectIterator it2 = vIn.begin();

    for (it1 = vect.begin(); it1 != vect.end(); ++it1)
    {
      **it1 = **it2;
      it2++;
    }
  }
  else // use of a temporary SimpleVector... bad way, to be improved. But this case happens rarely ...
  {
    //       std::cout << "WARNING: BlockVector::operator = BlockVector, v=w, with v and w having blocks of different sizes. This operation may be time-consuming and inappropriate." << std::endl;
    //       SiconosVector * tmp = new SimpleVector(vIn);
    //       *this = *tmp;
    //       delete tmp;

    // component by component copy ...
    for (unsigned int i = 0; i < sizeV; ++i)
      (*this)(i) = vIn(i);

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

//=================================
// Op. and assignment (+=, -= ... )
//=================================

BlockVector& BlockVector::operator += (const SiconosVector& vIn)
{
  if (&vIn == this)
  {
    BlockVectIterator it1;
    for (it1 = vect.begin(); it1 != vect.end(); ++it1)
    {
      **it1 += **it1;
    }
    return *this;
  }

  if (vIn.size() != sizeV)
    SiconosVectorException::selfThrow("BlockVector:operator += vIn, inconsistent size between this and vIn.");

  if (vIn.isBlock())
  {
    if (isComparableTo(this, &vIn)) // if vIn and this are "block-consistent"
    {
      BlockVectIterator it1;
      ConstBlockVectIterator it2 = vIn.begin();

      for (it1 = vect.begin(); it1 != vect.end(); ++it1)
      {
        **it1 += **it2;
        it2++;
      }
    }
    else // use of a temporary SimpleVector... bad way, to be improved. But this case happens rarely ...
    {
      //      SiconosVector * tmp = new SimpleVector(vIn);
      //      *this += *tmp;
      //      delete tmp;
      for (unsigned int i = 0; i < sizeV; ++i)
        (*this)(i) += vIn(i);
    }
  }
  else // if vIn is a Simple
  {
    unsigned int index = 0;
    addSimple(index, vIn);
  }
  // a call to a subfunction is required since the current vector
  // may be a Block of Block (...of Blocks...)

  return *this;
}

BlockVector& BlockVector::operator -= (const SiconosVector& vIn)
{
  if (&vIn == this)
  {
    this->zero();
    return *this;
  }
  if (vIn.size() != sizeV)
    SiconosVectorException::selfThrow("BlockVector:operator -= vIn, inconsistent size between this and vIn.");
  if (vIn.isBlock())
  {
    if (isComparableTo(this, &vIn)) // if vIn and this are "block-consistent"
    {
      unsigned int i = 0;
      BlockVectIterator it1;

      for (it1 = vect.begin(); it1 != vect.end(); ++it1)
        **it1 -= *(vIn[i++]);
    }
    else // use of a temporary SimpleVector... bad way, to be improved. But this case happens rarely ...
    {
      //    SiconosVector * tmp = new SimpleVector(vIn);
      //    *this -= *tmp;
      //    delete tmp;
      for (unsigned int i = 0; i < sizeV; ++i)
        (*this)(i) -= vIn(i);
    }
  }
  else // if vIn is a Simple
  {
    unsigned int index = 0;
    subSimple(index, vIn);
  }
  return *this;
}

void BlockVector::insert(const  SiconosVector& v)
{
  sizeV += v.size();

#ifndef WithSmartPtr
  if (!v.isBlock())
    vect.push_back(new SimpleVector(v)); // Copy
  else
    vect.push_back(new BlockVector(v)); // Copy

  isBlockAllocatedIn->push_back(true);
#else
  if (!v.isBlock())
    vect.push_back(boost::shared_ptr<SiconosVector>(new SimpleVector(v))); // Copy
  else
    vect.push_back(boost::shared_ptr<BlockVector>(new BlockVector(v))); // Copy
#endif

  tabIndex->push_back(sizeV);
}

void BlockVector::insertPtr(SiconosVectorSPtr v)
{
  if (!v)
    SiconosVectorException::selfThrow("BlockVector:insertPtr(v), v is a NULL vector.");

  sizeV += v->size();
  vect.push_back(v);

#ifndef WithSmartPtr
  isBlockAllocatedIn->push_back(false);
#endif

  tabIndex->push_back(sizeV);
}
