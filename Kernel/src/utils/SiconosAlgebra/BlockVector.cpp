/* Siconos-Kernel, Copyright INRIA 2005-2010.
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

#include <boost/numeric/ublas/vector_proxy.hpp>  // for project
#include <vector>

#include "BlockVector.hpp"
#include "SimpleVector.hpp"


// =================================================
//                CONSTRUCTORS
// =================================================
BlockVector::BlockVector(): SiconosVector()
{
  _sizeV = 0;
  _tabIndex.reset(new Index());
}

BlockVector::BlockVector(const std::string & file, bool ascii): SiconosVector()
{
  _sizeV = 0;

  SiconosVectorException::selfThrow(" BlockVector::constructor from a file : read BlockVector is not implemented");
}

BlockVector::BlockVector(const BlockVector &v): SiconosVector()
{
  _sizeV = 0;
  unsigned int nbBlocks = v.getNumberOfBlocks();
  _tabIndex.reset(new Index());
  _tabIndex->reserve(nbBlocks);
  vect.reserve(nbBlocks);
  VectorOfVectors::const_iterator it;
  for (it = v.begin(); it != v.end(); ++it)
  {
    if (!(*it)->isBlock())  // Call copy-constructor of SimpleVector
      vect.push_back(boost::shared_ptr<SimpleVector>(new SimpleVector(**it))) ;
    else
      vect.push_back(boost::shared_ptr<BlockVector>(new BlockVector(**it))) ;
    _sizeV += (*it)->size();
    _tabIndex->push_back(_sizeV);
  }
}

BlockVector::BlockVector(const SiconosVector &v): SiconosVector()
{
  _sizeV = 0;

  unsigned int nbBlocks = v.getNumberOfBlocks();
  _tabIndex.reset(new Index());

  _tabIndex->reserve(nbBlocks);
  vect.reserve(nbBlocks);

  if (v.isBlock())
  {
    VectorOfVectors::const_iterator it;
    for (it = v.begin(); it != v.end(); ++it)
    {

      if (!(*it)->isBlock())  // Call copy-constructor of SimpleVector
        vect.push_back(boost::shared_ptr<SimpleVector>(new SimpleVector(**it))) ;
      else
        vect.push_back(boost::shared_ptr<BlockVector>(new BlockVector(**it))) ;
      _sizeV += (*it)->size();
      _tabIndex->push_back(_sizeV);
    }
  }
  else
  {
    // Call copy-constructor of SimpleVector
    vect.push_back(boost::shared_ptr<SimpleVector>(new SimpleVector(v)));
    _sizeV = v.size();
    _tabIndex->push_back(_sizeV);
  }
}

BlockVector::BlockVector(SP::SiconosVector v1, SP::SiconosVector v2): SiconosVector()
{
  _sizeV = 0;

  // Insert the two vectors in the container
  // NO COPY !!
  if (! v1  && ! v2)
    SiconosVectorException::selfThrow("BlockVector:constructor(SimpleVector*,SimpleVector*), both vectors are NULL.");

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
    // If first parameter is a NULL pointer, then set this(1) to a SimpleVector of the same size as v2, and equal to 0.
  {
    // This case is usefull to set xDot in LagrangianDS.
    _sizeV = v2->size();

    vect.push_back(boost::shared_ptr<SimpleVector>(new SimpleVector(_sizeV)));
    _tabIndex->push_back(_sizeV);

  }
  if (v2)
  {
    vect.push_back(v2);
    _sizeV += v2->size();
    _tabIndex->push_back(_sizeV);

  }
  else // If second parameter is a NULL pointer, then set this(2) to a SimpleVector of the same size as v1, and equal to 0.
  {
    // This case is usefull to set xDot in LagrangianDS.

    vect.push_back(boost::shared_ptr<SimpleVector>(new SimpleVector(v1->size())));
    _sizeV += v1->size();
    _tabIndex->push_back(_sizeV);
  }
}

BlockVector::BlockVector(unsigned int numberOfBlocks, unsigned int dim): SiconosVector()
{
  _sizeV = 0;

  _tabIndex.reset(new Index());
  _tabIndex->reserve(numberOfBlocks);
  vect.reserve(numberOfBlocks);
  for (unsigned int i = 0; i < numberOfBlocks; ++i)
  {
    vect.push_back(boost::shared_ptr<SimpleVector>(new SimpleVector(dim)));
    _tabIndex->push_back(dim * (i + 1));
  }
  _sizeV = dim * numberOfBlocks;
}

BlockVector::~BlockVector()
{}

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

DenseVect* BlockVector::dense(unsigned int i) const
{
  if (vect[i]->getNum() != 1)
    SiconosVectorException::selfThrow("BlockVector::dense(unsigned int num) : the vector[num] is not a Dense vector");
  return (vect[i])->dense();
}

SparseVect* BlockVector::sparse(unsigned int i) const
{
  if (vect[i]->getNum() != 4)
    SiconosVectorException::selfThrow("BlockVector::sparse(unsigned int num) : the vector[num] is not a Sparse vector");
  return (vect[i])->sparse();
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
  VectorOfVectors::const_iterator it;
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
  VectorOfVectors::const_iterator it;
  for (it = vect.begin(); it != vect.end(); ++it)
  {
    assert(*it);
    d += pow((*it)->norm2(), 2);
  }
  return sqrt(d);
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

const double BlockVector::operator()(unsigned int pos) const
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

SimpleVector BlockVector::getVector(unsigned int pos) const
{

  if (! vect[pos])
    SiconosVectorException::selfThrow("BlockVector::getVector(pos), vector[pos] == NULL pointer.");

  if (vect[pos]->isBlock())
    SiconosVectorException::selfThrow("BlockVector::getVector(pos), vector[pos] is a Block. Use vector()");

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

void BlockVector::addSimple(unsigned int& index, const SiconosVector& vIn)
{
  // Add a part of vIn (starting from index) to the current vector.
  // vIn must be a SimpleVector.

  // At the end of the present function, index is equal to index + the dim. of the added sub-vector.

  unsigned int dim = vIn.size() - index; // size of the block to be added.
  if (dim > _sizeV) SiconosVectorException::selfThrow("BlockVector::addSimple : invalid ranges");

  VectorOfVectors::const_iterator it;
  unsigned int numVIn = vIn.getNum();
  unsigned int currentSize, currentNum;

  for (it = vect.begin(); it != vect.end(); ++it)
  {
    if ((*it)->isBlock())
      (boost::static_pointer_cast<BlockVector>(*it))->addSimple(index, vIn);
    else // the current block is a SimpleVector
    {
      currentSize = (*it)->size();
      currentNum = (*it)->getNum();
      if (numVIn != currentNum) SiconosVectorException::selfThrow("BlockVector::addSimple : inconsistent types.");
      if (numVIn == 1)
        noalias(*(*it)->dense()) +=  ublas::subrange(*vIn.dense(), index, index + currentSize) ;
      else
        noalias(*(*it)->sparse()) +=  ublas::subrange(*vIn.sparse(), index, index + currentSize) ;
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
  if (dim > _sizeV) SiconosVectorException::selfThrow("BlockVector::addSimple : invalid ranges");

  VectorOfVectors::const_iterator it;
  unsigned int numVIn = vIn.getNum();
  unsigned int currentSize, currentNum;

  for (it = vect.begin(); it != vect.end(); ++it)
  {
    if ((*it)->isBlock())
      (boost::static_pointer_cast<BlockVector>(*it))->subSimple(index, vIn);

    else // the current block is a SimpleVector
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
  }
}

//===============
//  Assignment
//===============
BlockVector& BlockVector::operator = (const SiconosVector& vIn)
{
  if (&vIn == this) return *this;

  if (vIn.size() != _sizeV)
    SiconosVectorException::selfThrow("BlockVector:operator = vIn, inconsistent size between this and vIn.");

  if (vIn.isBlock())
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
    else // use of a temporary SimpleVector... bad way, to be improved. But this case happens rarely ...
    {
      //    std::cout << "WARNING: BlockVector::operator = BlockVector, v=w, with v and w having blocks of different sizes. This operation may be time-consuming and inappropriate." << std::endl;
      //    SiconosVector * tmp = new SimpleVector(vIn);
      //    *this = *tmp;

      // component by component copy ...
      for (unsigned int i = 0; i < _sizeV; ++i)
        (*this)(i) = vIn(i);

    }
  }
  else // if vIn is a SimpleVector
  {
    setBlock(vIn, shared_from_this(), _sizeV, 0, 0);

    //       unsigned int pos = 0;
    //       VectorOfVectors::iterator it1;
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

  if (vIn.size() != _sizeV)
    SiconosVectorException::selfThrow("BlockVector:operator = vIn, inconsistent size between this and vIn.");

  if (isComparableTo(*this, vIn)) // if vIn and this are "block-consistent"
  {
    VectorOfVectors::iterator it1;
    VectorOfVectors::const_iterator it2 = vIn.begin();

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

    // component by component copy ...
    for (unsigned int i = 0; i < _sizeV; ++i)
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
    VectorOfVectors::iterator it1;
    for (it1 = vect.begin(); it1 != vect.end(); ++it1)
    {
      **it1 += **it1;
    }
    return *this;
  }

  if (vIn.size() != _sizeV)
    SiconosVectorException::selfThrow("BlockVector:operator += vIn, inconsistent size between this and vIn.");

  if (vIn.isBlock())
  {
    if (isComparableTo(*this, vIn)) // if vIn and this are "block-consistent"
    {
      VectorOfVectors::iterator it1;
      VectorOfVectors::const_iterator it2 = vIn.begin();

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
      for (unsigned int i = 0; i < _sizeV; ++i)
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
  if (vIn.size() != _sizeV)
    SiconosVectorException::selfThrow("BlockVector:operator -= vIn, inconsistent size between this and vIn.");
  if (vIn.isBlock())
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
      //    SiconosVector * tmp = new SimpleVector(vIn);
      //    *this -= *tmp;
      for (unsigned int i = 0; i < _sizeV; ++i)
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
  _sizeV += v.size();

  if (!v.isBlock())
    vect.push_back(boost::shared_ptr<SiconosVector>(new SimpleVector(v))); // Copy
  else
    vect.push_back(boost::shared_ptr<BlockVector>(new BlockVector(v))); // Copy

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
