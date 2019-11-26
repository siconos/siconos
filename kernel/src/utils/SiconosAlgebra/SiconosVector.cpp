/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

#include "SiconosConfig.h"

#include <boost/numeric/ublas/io.hpp>    
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>                 // for matri...
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/vector_proxy.hpp>
#include <boost/numeric/bindings/blas.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/std/vector.hpp>
#include "SiconosAlgebraTypeDef.hpp"

#include "SiconosAlgebra.hpp"

#include "SiconosVectorIterator.hpp"
#include "Tools.hpp"

namespace siconosBindings = boost::numeric::bindings::blas;

#include "SimpleMatrix.hpp"
#include "BlockVector.hpp"
#include "ioVector.hpp"
#include "SiconosVector.hpp"
#include "SiconosAlgebra.hpp"
#include <cmath>        // std::exp(double)
#include <algorithm>    // std::transform

//#define DEBUG_MESSAGES
#include "debug.h"

// Do not document
/// @cond
#include "Question.hpp"

struct IsDense : public Question<bool>
{
  using SiconosVisitor::visit;

  void visit(const SiconosVector& v)
  {
    answer = v._dense;
  }

  void visit(const BlockVector& v)
  {
    answer = false;
  }
};

struct IsSparse : public Question<bool>
{

  using SiconosVisitor::visit;

  void visit(const SiconosVector& v)
  {
    answer = !v._dense;
  }

  void visit(const BlockVector& v)
  {
    answer = false;
  }
};

struct IsBlock : public Question<bool>
{
  using SiconosVisitor::visit;

  void visit(const SiconosVector& v)
  {
    answer = false;
  }

  void visit(const BlockVector& v)
  {
    answer = true;
  }
};

/// @endcond


// =================================================
//                CONSTRUCTORS
// =================================================

// Default
SiconosVector::SiconosVector()
{
  _dense = true;
  vect.Dense = new DenseVect(ublas::zero_vector<double>());
}

// parameters: dimension and type.
SiconosVector::SiconosVector(unsigned row, Siconos::UBLAS_TYPE type)
{
  if (type == Siconos::SPARSE)
  {
    _dense = false;
    vect.Sparse = new SparseVect(ublas::zero_vector<double>(row));
  }
  else if (type == Siconos::DENSE)
  {
    _dense = true;
    vect.Dense = new DenseVect(ublas::zero_vector<double>(row));
  }
  else
  {
    SiconosVectorException::selfThrow("SiconosVector::constructor(Siconos::UBLAS_TYPE, unsigned int) failed, invalid type given");
  }
}

// parameters: dimension, default value for all components and type.
SiconosVector::SiconosVector(unsigned row, double val, Siconos::UBLAS_TYPE type)
{
  if (type == Siconos::SPARSE)
  {
    _dense = false;
    vect.Sparse = new SparseVect(row);
    fill(val);
  }
  else if (type == Siconos::DENSE)
  {
    _dense = true;
    vect.Dense = new DenseVect(ublas::scalar_vector<double>(row, val));
  }
  else
  {
    SiconosVectorException::selfThrow("SiconosVector::constructor(Siconos::UBLAS_TYPE, unsigned int) : invalid type given");
  }
}

// parameters: a vector (stl) of double and the type.
SiconosVector::SiconosVector(const std::vector<double>& v, Siconos::UBLAS_TYPE typ)
{
  if (typ != Siconos::DENSE)
    SiconosVectorException::selfThrow("SiconosVector::constructor(Siconos::UBLAS_TYPE, std::vector<double>, unsigned int) : invalid type given");

  _dense = true;
  vect.Dense = new DenseVect(v.size());
  std::copy(v.begin(), v.end(), (vect.Dense)->begin());
}

// Copy
SiconosVector::SiconosVector(const SiconosVector &svect) : std11::enable_shared_from_this<SiconosVector>()
{
  if (ask<IsDense>(svect)) // dense
  {
    _dense = true;
    vect.Dense = new DenseVect(svect.size());
    noalias(*vect.Dense) = (*svect.dense());
    // std::copy((vect.Dense)->begin(), (vect.Dense)->end(), (svect.dense())->begin());
  }
  else //sparse
  {
    _dense = false;
    vect.Sparse = new SparseVect(svect.size());
    noalias(*vect.Sparse) = (*svect.sparse());
    //std::copy((vect.Sparse)->begin(), (vect.Sparse)->end(), (svect.sparse())->begin());
  }

  // Note FP: using constructor + noalias = (or std::copy) is more
  // efficient than a call to ublas::vector copy constructor, this for
  // large or small vectors.
}

// Copy from BlockVector
SiconosVector::SiconosVector(const BlockVector & vIn) : std11::enable_shared_from_this<SiconosVector>()
{
  if (ask<IsDense>(**(vIn.begin()))) // dense
  {
    _dense = true;
    vect.Dense = new DenseVect(vIn.size());
  }
  else
  {
    _dense = false;
    vect.Sparse = new SparseVect(vIn.size());
  }

  VectorOfVectors::const_iterator it;
  unsigned int pos = 0;
  for (it = vIn.begin(); it != vIn.end(); ++it)
  {
    setBlock(pos, **it);
    pos += (*it)->size();
  }

}

SiconosVector::SiconosVector(const DenseVect& m)
{
  _dense = true;
  vect.Dense = new DenseVect(m.size());
  noalias(*vect.Dense) = m;

}

SiconosVector::SiconosVector(const SparseVect& m)
{
  _dense = false;
  vect.Sparse = new SparseVect(m.size());
  noalias(*vect.Sparse) = m;
}

SiconosVector::SiconosVector(const std::string &file, bool ascii)
{
  _dense = true;
  vect.Dense = new DenseVect();
  if (ascii)
  {
    ioVector::read(file, *this, ioVector::ASCII_IN);
   }
  else
  {
    ioVector::read(file, *this, ioVector::BINARY_IN);
  }
}

SiconosVector::SiconosVector(const SiconosVector& v1, const SiconosVector& v2)
{
  unsigned int size1 = v1.size();
  if (ask<IsDense>(v1) && ask<IsDense>(v2))
  {
    _dense = true;
    vect.Dense = new DenseVect(size1 + v2.size());
  }
  else if(ask<IsSparse>(v1) && ask<IsSparse>(v2))
  {
    _dense = false;
    vect.Sparse = new SparseVect(size1 + v2.size());
  }
  else
  {
    SiconosVectorException::selfThrow("SiconosVector::SiconosVector :: mixed dense and sparse vector detected");
  }
  setBlock(0, v1);
  setBlock(size1, v2);
}

SiconosVector::~SiconosVector()
{
  if(_dense)
    delete(vect.Dense);
  else delete(vect.Sparse);
}


// =================================================
//        get Ublas component (dense or sparse)
// =================================================

const DenseVect SiconosVector::getDense(unsigned int) const
{
  if(!_dense)
    SiconosVectorException::selfThrow("SiconosVector::getDense(unsigned int row, unsigned int col) : the current vector is not a Dense vector");

  return *vect.Dense;
}

const SparseVect SiconosVector::getSparse(unsigned int)const
{

  if(_dense)
    SiconosVectorException::selfThrow("SiconosVector::getSparse(unsigned int row, unsigned int col) : the current vector is not a Sparse vector");

  return *vect.Sparse;
}

SparseVect* SiconosVector::sparse(unsigned int)const
{

  if(_dense)
    SiconosVectorException::selfThrow("SiconosVector::sparse(unsigned int row, unsigned int col) : the current vector is not a Sparse vector");

  return vect.Sparse;
}

double* SiconosVector::getArray() const
{
  assert(vect.Dense && "SiconosVector::getArray() : not yet implemented for sparse vector.");

  return &(((*vect.Dense).data())[0]);
}

// ===========================
//       fill vector
// ===========================

void SiconosVector::zero()
{
  if(_dense)
    siconosBindings::scal(0.0, *vect.Dense);

  else
  {
    assert(vect.Sparse);
    *vect.Sparse *= 0.0;
  }

}

void SiconosVector::setVector(unsigned int, const SiconosVector& newV)
{
  if(newV.size() != size())
    SiconosVectorException::selfThrow("SiconosVector::setVector(num,v), unconsistent sizes.");

  *this = newV ;
}

void SiconosVector::fill(double value)
{
  if(!_dense)
  {
    for(unsigned int i = 0; i < (vect.Sparse)->size(); ++i)
      (vect.Sparse)->push_back(i, value);
  }
  else
    siconosBindings::set(value, *vect.Dense);


}

//=======================
// set vector dimension
//=======================

void SiconosVector::resize(unsigned int n, bool preserve)
{
  if(_dense)
    (vect.Dense)->resize(n, preserve);
  else
    (vect.Sparse)->resize(n, preserve);
}

//=======================
//       get norm
//=======================

double SiconosVector::normInf() const
{
  if(_dense)
    return norm_inf(*vect.Dense);
  else //if(num==4)
    return norm_inf(*vect.Sparse);
}

double SiconosVector::norm2() const
{
  if(_dense)
    return ublas::norm_2(*vect.Dense);
  else //if(num==4)
    return ublas::norm_2(*vect.Sparse);
}
//======================================
// get sum of all elements of the vector
//=====================================
double SiconosVector::vector_sum() const
{
  if(_dense)
    return ublas::sum(*vect.Dense);
  else
    return ublas::sum(*vect.Sparse);
}

//=====================
// screen display
//=====================

void SiconosVector::display()const
{
  std::cout.setf(std::ios::scientific);
  std::cout.precision(6);
  if(_dense)
    std::cout << *vect.Dense << std::endl;
  else if(vect.Sparse)
    std::cout << *vect.Sparse << std::endl;
}

//============================
// Convert vector to a std::string
//============================

std::string SiconosVector::toString() const
{
  return ::toString(*this);
}

//=============================
// Elements access (get or set)
//=============================

double SiconosVector::getValue(unsigned int row) const
{
  assert(row < size() && "SiconosVector::getValue(index) : Index out of range");

  if(_dense)
    return (*vect.Dense)(row);
  else
    return (*vect.Sparse)(row);
}

void SiconosVector::setValue(unsigned int row, double value)
{
  assert(row < size() && "SiconosVector::setValue(index, value) : Index out of range");
  if(_dense)
    (*vect.Dense)(row) = value ;
  else
    (*vect.Sparse)(row) = value;
}

double& SiconosVector::operator()(unsigned int row)
{
  assert(row < size() && "SiconosVector::operator ( index ): Index out of range");

  if(_dense)
    return (*vect.Dense)(row);
  else
    return (*vect.Sparse)(row).ref();
}

double SiconosVector::operator()(unsigned int row) const
{
  assert(row < size() && "SiconosVector::operator ( index ): Index out of range");

  if(_dense)
    return (*vect.Dense)(row);
  else
    return ((*vect.Sparse)(row)).ref();
}

//============================================
// Access (get or set) to blocks of elements
//============================================

void SiconosVector::setBlock(unsigned int index, const SiconosVector& vIn)
{
  // Set current vector elements, starting from position "index", to the values of vector vIn

  // Exceptions ...
  assert(&vIn != this && "SiconosVector::this->setBlock(pos,vIn): vIn = this.");

  assert(index < size() && "SiconosVector::setBlock : invalid ranges");

  unsigned int end = vIn.size() + index;
  assert(end <= size() && "SiconosVector::setBlock : invalid ranges");

  assert(vIn.num() == num() && "SiconosVector::setBlock: inconsistent types.");

  if(_dense)
    noalias(ublas::subrange(*vect.Dense, index, end)) = *vIn.dense();
  else
    noalias(ublas::subrange(*vect.Sparse, index, end)) = *vIn.sparse();
}

void SiconosVector::toBlock(SiconosVector& vOut, unsigned int sizeB, unsigned int startIn, unsigned int startOut) const
{
  // To copy a subBlock of the vector (from position startIn to startIn+sizeB) into vOut (from pos. startOut to startOut+sizeB).
  // Check dim ...
  assert(startIn < size() && "vector toBlock(v1,v2,...): start position in input vector is out of range.");

  assert(startOut < vOut.size() && "vector toBlock(v1,v2,...): start position in output vector is out of range.");

  assert(startIn + sizeB <= size() && "vector toBlock(v1,v2,...): end position in input vector is out of range.");
  assert(startOut + sizeB <= vOut.size() && "vector toBlock(v1,v2,...): end position in output vector is out of range.");

  unsigned int endOut = startOut + sizeB;
  unsigned int numIn = num();
  unsigned int numOut = vOut.num();

  if(numIn == numOut)
  {
    if(numIn == 1)  // vIn / vOut are Dense
      noalias(ublas::subrange(*vOut.dense(), startOut, endOut)) = ublas::subrange(*vect.Dense, startIn, startIn + sizeB);
    else // if(numIn == 4)// vIn / vOut are Sparse
      noalias(ublas::subrange(*vOut.sparse(), startOut, endOut)) = ublas::subrange(*vect.Sparse, startIn, startIn + sizeB);
  }
  else // vIn and vout of different types ...
  {
    if(numIn == 1)  // vIn Dense
      noalias(ublas::subrange(*vOut.sparse(), startOut, endOut)) = ublas::subrange(*vect.Dense, startIn, startIn + sizeB);
    else // if(numIn == 4)// vIn Sparse
      noalias(ublas::subrange(*vOut.dense(), startOut, endOut)) = ublas::subrange(*vect.Sparse, startIn, startIn + sizeB);
  }
}

void SiconosVector::addBlock(unsigned int index, const SiconosVector& vIn)
{
  // Add vIn to the current vector, starting from position "index".
  // vIn may be a BlockVector.

  //if ( num != 1 ) SiconosVectorException::selfThrow("SiconosVector::addBlock : vector should be dense");

  if(&vIn == this)
    SiconosVectorException::selfThrow("SiconosVector::this->addBlock(pos,vIn): vIn = this.");

  unsigned int end = vIn.size();
  if((index + end) > size()) SiconosVectorException::selfThrow("SiconosVector::addBlock : invalid ranges");

  unsigned int numVin = vIn.num();

  if(numVin != num()) SiconosVectorException::selfThrow("SiconosVector::addBlock : inconsistent types.");

  if(_dense)
    noalias(ublas::subrange(*vect.Dense, index, index + end)) += *vIn.dense();
  else
    noalias(ublas::subrange(*vect.Sparse, index, index + end)) += *vIn.sparse();
}

void SiconosVector::subBlock(unsigned int index, const SiconosVector& vIn)
{
  // Add vIn from the current vector, starting from position "index".
  // vIn may be a BlockVector.

  //  if ( num != 1 ) SiconosVectorException::selfThrow("SiconosVector::subBlock : vector should be dense");

  unsigned int end = vIn.size();
  if((index + end) > size()) SiconosVectorException::selfThrow("SiconosVector::subBlock : invalid ranges");

  unsigned int numVin = vIn.num();
  if(numVin != num()) SiconosVectorException::selfThrow("SiconosVector::subBlock : inconsistent types.");

  if(_dense)
    noalias(ublas::subrange(*vect.Dense, index, index + end)) -= *vIn.dense();
  else
    noalias(ublas::subrange(*vect.Sparse, index, index + end)) -= *vIn.sparse();
}

//===============
//  Assignment
//===============

SiconosVector& SiconosVector::operator = (const SiconosVector& vIn)
{
  if(&vIn == this) return *this;  // auto-assignment.

  assert(size() == vIn.size() && "SiconosVector::operator = failed: inconsistent sizes.");

  unsigned int vInNum = vIn.num();
  {
    switch(num())
    {
    case 1:
      switch(vInNum)
      {
      case 1:
        //siconosBindings::copy(*vIn.dense(),*vect.Dense);
        noalias(*vect.Dense) = *vIn.dense();
        break;
      case 4:
        noalias(*vect.Dense) = *vIn.sparse();
        break;
      default:
        SiconosVectorException::selfThrow("SiconosVector::operator = : invalid type given");
        break;
      }
      break;
    case 4:
      if(vInNum == 4)
        noalias(*vect.Sparse) = *vIn.sparse();
      else
        SiconosVectorException::selfThrow("SiconosVector::operator = : can not set sparse = dense.");
      break;
    default:
      SiconosVectorException::selfThrow("SiconosVector::operator = : invalid type given");
      break;
    }
  }
  return *this;
}

SiconosVector& SiconosVector::operator = (const BlockVector& vIn)
{
  VectorOfVectors::const_iterator it;
  unsigned int pos = 0;
  for(it = vIn.begin(); it != vIn.end(); ++it)
  {
    setBlock(pos, **it);
    pos += (*it)->size();
  }
  return *this;
}


SiconosVector& SiconosVector::operator = (const DenseVect& d)
{
  if(!_dense)
    SiconosVectorException::selfThrow("SiconosVector::operator = DenseVect : forbidden: the current vector is not dense.");
  if(d.size() != size())
    SiconosVectorException::selfThrow("SiconosVector::operator = DenseVect : inconsistent size.");

  siconosBindings::copy(d, *vect.Dense);
  return *this;
}

SiconosVector& SiconosVector::operator = (const SparseVect& sp)
{
  if(_dense)
    SiconosVectorException::selfThrow("SiconosVector::operator = SparseVect : current vector is not sparse.");
  if(sp.size() != size())
    SiconosVectorException::selfThrow("SiconosVector::operator = SparseVect : inconsistent size.");

  noalias(*vect.Sparse) = sp;

  return *this;
}

SiconosVector& SiconosVector::operator = (const double* d)
{
  assert(_dense && "SiconosVector::operator = double* : forbidden: the current vector is not dense.");

  siconosBindings::detail::copy(vect.Dense->size(), d, 1, getArray(), 1);
  return *this;
}

unsigned SiconosVector::copyData(double* data) const
{
  assert(_dense && "SiconosVector::copyData : forbidden: the current vector is not dense.");

  unsigned size = vect.Dense->size();
  siconosBindings::detail::copy(vect.Dense->size(), getArray(), 1, data, 1);
  return size;
}


//=================================
// Op. and assignment (+=, -= ... )
//=================================

SiconosVector& SiconosVector::operator += (const SiconosVector& vIn)
{
  if(&vIn == this)  // alias
  {
    // Note: using this *= 2.0 is much more time-consuming.
    switch(num())
    {
    case 1:
      *vect.Dense += *vect.Dense;
      break;
    case 4:
      *vect.Sparse += *vect.Sparse;
      break;
    default:
      SiconosVectorException::selfThrow("SiconosVector::operator += : invalid type given");
      break;
    }
    return *this;
  }

  unsigned int vInNum = vIn.num();
  {
    switch(num())
    {
    case 1:
      switch(vInNum)
      {
      case 1:
        noalias(*vect.Dense) += *vIn.dense();
        break;
      case 4:
        noalias(*vect.Dense) += *vIn.sparse();
        break;
      default:
        SiconosVectorException::selfThrow("SiconosVector::operator += : invalid type given");
        break;
      }
      break;
    case 4:
      if(vInNum == 4)
        noalias(*vect.Sparse) += *vIn.sparse();
      else SiconosVectorException::selfThrow("SiconosVector::operator += : can not add a dense to a sparse.");
      break;
    default:
      SiconosVectorException::selfThrow("SiconosVector::operator += : invalid type given");
      break;
    }
  }
  return *this;
}
SiconosVector& SiconosVector::operator += (const BlockVector& vIn)
{
  VectorOfVectors::const_iterator it;
  unsigned int pos = 0;
  for(it = vIn.begin(); it != vIn.end(); ++it)
  {
    addBlock(pos, **it);
    pos += (*it)->size();
  }
  return *this;
}

SiconosVector& SiconosVector::operator -= (const SiconosVector& vIn)
{
  if(&vIn == this)
  {
    this->zero();
    return *this;
  }

  unsigned int vInNum = vIn.num();
  {
    switch(num())
    {
    case 1:
      switch(vInNum)
      {
      case 1:
        noalias(*vect.Dense) -= *vIn.dense();
        break;
      case 4:
        noalias(*vect.Dense) -= *vIn.sparse();
        break;
      default:
        SiconosVectorException::selfThrow("SiconosVector::operator -= : invalid type given");
        break;
      }
      break;
    case 4:
      if(vInNum == 4)
        noalias(*vect.Sparse) -= *vIn.sparse();
      else SiconosVectorException::selfThrow("SiconosVector::operator -= : can not sub a dense to a sparse.");
      break;
    default:
      SiconosVectorException::selfThrow("SiconosVector::operator -= : invalid type given");
      break;
    }
  }
  return *this;
}

SiconosVector& SiconosVector::operator -= (const BlockVector& vIn)
{
  VectorOfVectors::const_iterator it;
  unsigned int pos = 0;
  for(it = vIn.begin(); it != vIn.end(); ++it)
  {
    subBlock(pos, **it);
    pos += (*it)->size();
  }
  return *this;
}



struct exp_op
{
  double operator()(double d) const
  {
    return std::exp(d);
  }
};

void SiconosVector::exp_in_place()
{
  // struct exp_op { double operator() (double d) const { return std::exp(d); } };
  // assert(num() == 1);
  // std::transform(vect.Dense->begin(), vect.Dense->end(), vect.Dense->begin(), exp_op);
}

void SiconosVector::exp(SiconosVector& input)
{
  // assert(num() == 1 && input.num()==1);
  // std::transform(input.dense()->begin(), input.dense()->end(), vect.Dense->begin(), exp_op);
}

unsigned int SiconosVector::size(void) const
{
  if(!_dense)
  {
    return (vect.Sparse->size());
  }
  else
  {
    return (vect.Dense->size());
  }
}

SiconosVector::iterator SiconosVector::begin()
{
  return SiconosVector::iterator(*this, 0);
}

SiconosVector::const_iterator SiconosVector::begin() const
{
  return SiconosVector::const_iterator(*this, 0);
}

SiconosVector::iterator SiconosVector::end()
{
  return SiconosVector::iterator(*this, size());
}

SiconosVector::const_iterator SiconosVector::end() const
{
  return SiconosVector::const_iterator(*this, size());
}

SiconosVector::operator std::vector<double>()
{
  std::vector<double> v;
  v.resize(size());
  std::copy(begin(), end(), v.begin());
  return v;
}


void SiconosVector::private_prod(SPC::SiconosVector x, SPC::SiconosMatrix A, unsigned int startCol, bool init)
{
  assert(!(A->isPLUFactorized()) && "A is PLUFactorized in prod !!");

  // Computes y = subA *x (or += if init = false), subA being a sub-matrix of trans(A), between el. of A of index (col) startCol and startCol + sizeY
  if (init) // y = subA * x , else y += subA * x
    zero();
  private_addprod(x, A, startCol, 0);

}

void SiconosVector::private_addprod(SPC::SiconosVector x, SPC::SiconosMatrix A, unsigned int startRow, unsigned int startCol)
{
  assert(!(A->isPLUFactorized()) && "A is PLUFactorized in prod !!");

  if (A->isBlock())
    SiconosMatrixException::selfThrow("private_addprod(x,A,start,y) error: not yet implemented for block matrix.");

  // we take a submatrix subA of A, starting from row startRow to row (startRow+sizeY) and between columns startCol and (startCol+sizeX).
  // Then computation of y = subA*x + y.
  unsigned int numA = A->num();
  unsigned int numY = num();
  unsigned int numX = x->num();
  unsigned int sizeX = x->size();
  unsigned int sizeY = size();

  if (numX != numY)
    SiconosMatrixException::selfThrow("private_addprod(x,A,start) error: not yet implemented for x and y of different types.");

  if (num() == 1 && numX == 1)
  {

    assert(dense() != x->dense());

    if (numA == 1)
      noalias(*dense()) += prod(ublas::subrange(trans(*A->dense()), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
    else if (numA == 2)
      noalias(*dense()) += prod(ublas::subrange(trans(*A->triang()), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
    else if (numA == 3)
      noalias(*dense()) += prod(ublas::subrange(trans(*A->sym()), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
    else if (numA == 4)
      noalias(*dense()) += prod(ublas::subrange(trans(*A->sparse()), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
    else //if(numA==5)
      noalias(*dense()) += prod(ublas::subrange(trans(*A->banded()), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->dense());
  }
  else // x and y sparse
  {
    if (numA == 4)
      *sparse() += prod(ublas::subrange(trans(*A->sparse()), startRow, startRow + sizeY, startCol, startCol + sizeX), *x->sparse());
    else
      SiconosMatrixException::selfThrow("private_addprod(x,A,start,y) error: not yet implemented for x, y  sparse and A not sparse.");
  }
}

