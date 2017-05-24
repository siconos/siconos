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

#include "SiconosConfig.h"

#include "SiconosAlgebraTypeDef.hpp"

#include <boost/numeric/ublas/io.hpp>            // for >>
//#include <boost/numeric/ublas/vector_proxy.hpp>  // for project
#include <boost/numeric/ublas/vector_sparse.hpp>


#include <boost/numeric/bindings/ublas/vector_proxy.hpp>
#include <boost/numeric/bindings/blas.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/std/vector.hpp>

#include "SiconosVectorStorage.hpp"
#include "SiconosVectorOperators.hpp"
#include "SimpleMatrix.hpp"

#include "ioVector.hpp"
#include "SiconosVector.hpp"
#include "SiconosAlgebra.hpp"
#include "BlockVector.hpp"
//#define DEBUG_MESSAGES
#include "debug.h"
#include "TypeName.hpp"

namespace bindings = boost::numeric::bindings::blas;

template<typename C>
struct IOuterProd : public ParamVisitor<const C&>
{

  SimpleMatrix answer;

  IOuterProd(const C& p) : ParamVisitor<const C&>(p) {};

  template<typename T>
  void operator() (const T& storage)
  {
    answer = outer_prod(storage.internal_data, this->param().internal_data);
  };
};


struct OuterProd : public ParamVisitor<const SiconosVectorStorage&>
{

  SimpleMatrix answer;

  OuterProd(const SiconosVectorStorage& st) : ParamVisitor(st) {};

  template<typename T>
  void operator() (const T& storage)
  {
    answer = apply_visitor<IOuterProd<T>, SimpleMatrix >(this->param(), storage);
  };
};
// =================================================
//                CONSTRUCTORS
// =================================================

// Default
SiconosVector::SiconosVector() : _storage(new DenseVectStorage()) {};

// parameters: dimension and type.
SiconosVector::SiconosVector(unsigned row, Siconos::UBLAS_TYPE type)
{
  if (type == Siconos::SPARSE)
  {
    _storage = new SparseVectStorage(row);
  }
  else if (type == Siconos::DENSE)
  {
    _storage = new DenseVectStorage(row);
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
    _storage = new SparseVectStorage(row);
    fill(val);
  }
  else if (type == Siconos::DENSE)
  {
    _storage = new DenseVectStorage(row, val);
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

 _storage = new DenseVectStorage(v.size());
 std::copy(v.begin(), v.end(), static_cast<DenseVectStorage*>(_storage)->internal_data.begin());
}

// Copy
SiconosVector::SiconosVector(const SiconosVector &svect) : boost::enable_shared_from_this<SiconosVector>()
{
  this->_storage = apply_visitor<StorageAllocator, SiconosVectorStorage*>(storage(svect), svect.size());
  // *this <- svect
  apply_visitor<Copy>(storage(svect), storage(*this));
}

// Copy from BlockVector
SiconosVector::SiconosVector(const BlockVector & vIn) : boost::enable_shared_from_this<SiconosVector>()
{
  this->_storage = apply_visitor<StorageAllocator, SiconosVectorStorage*>(storage(**vIn.begin()), vIn.size());
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
  this->_storage = new DenseVectStorage(m.size());
  noalias(this->dense()) = m;

}

SiconosVector::SiconosVector(const SparseVect& m)
{
  this->_storage = new SparseVectStorage(m.size());
  noalias(this->sparse()) = m;
}

SiconosVector::SiconosVector(const std::string &file, bool ascii)
{
  this->_storage = new DenseVectStorage();
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

  if (Type::value(storage(v1)) == Type::value(storage(v2)))
  {
    this->_storage = apply_visitor<StorageAllocator, SiconosVectorStorage*>(storage(v1), v1.size() + v2.size()) ;
  }
  else
  {
    SiconosVectorException::selfThrow("SiconosVector::SiconosVector :: mixed dense and sparse vector detected");
  }
  setBlock(0, v1);
  setBlock(size1, v2);
}

SiconosVector::SiconosVector(SiconosVectorStorage& storage) : _storage(&storage) {};

SiconosVector::~SiconosVector()
{
  if (_storage) delete(_storage);
}


// =================================================
//        get Ublas component (dense or sparse)
// =================================================

unsigned int SiconosVector::num() const
{
  if (Type::value(*_storage) == Type::DenseVectStorage)
  {
    return 1;
  }
  else
    if (Type::value(*_storage) == Type::SparseVectStorage)
    {
      return 4;
    }
    else
    {
      return 0;
    }
}

DenseVect& SiconosVector::dense(unsigned int) const
{
  if (Type::value(*_storage) != Type::DenseVectStorage)
  {
    SiconosVectorException::selfThrow("SiconosVector::Densevect(unsigned int) : cannot get dense storage.");
  }
  return static_cast<DenseVectStorage*>(_storage)->internal_data;
}

SparseVect& SiconosVector::sparse(unsigned int)const
{
  if (Type::value(*_storage) != Type::SparseVectStorage)
  {
    SiconosVectorException::selfThrow("SiconosVector::Sparsevect(unsigned int) : cannot get sparse storage.");
  }
  return static_cast<SparseVectStorage*>(_storage)->internal_data;
}

double* SiconosVector::getArray() const
{
  if (Type::value(*_storage) != Type::DenseVectStorage)
  {
    SiconosVectorException::selfThrow("SiconosVector::getArray() : cannot get array for this kind of vector.");
  }
  return &(((static_cast<DenseVectStorage*>(_storage)->internal_data).data())[0]);
}

// ===========================
//       fill vector
// ===========================

void SiconosVector::zero()
{
  apply_visitor<Zero>(storage(*this));
}

void SiconosVector::setVector(unsigned int , const SiconosVector& newV)
{
  if (newV.size() != size())
    SiconosVectorException::selfThrow("SiconosVector::setVector(num,v), unconsistent sizes.");

  *this = newV ;
}

void SiconosVector::fill(const double& value)
{
  apply_visitor<Fill>(storage(*this), value);
}

//=======================
// set vector dimension
//=======================

void SiconosVector::resize(unsigned int n, bool preserve)
{
  apply_visitor<Resize>(storage(*this), n);
}

//=======================
//       get norm
//=======================

double SiconosVector::normInf() const
{
  return apply_visitor<NormInf, double>(storage(*this));
}

double SiconosVector::norm2() const
{
  return apply_visitor<Norm2, double>(storage(*this));
}
//======================================
// get sum of all elements of the vector
//=====================================
double SiconosVector::vector_sum() const
{
  return apply_visitor<Sum, double>(storage(*this));
}

//=====================
// screen display
//=====================

void SiconosVector::display(unsigned int n)const
{
  apply_visitor<Display>(storage(*this), n);
}

//============================
// Convert vector to a std::string
//============================

const std::string SiconosVector::toString() const
{
  return apply_visitor<ToString, std::string>(storage(*this));
}

//=============================
// Elements access (get or set)
//=============================

double SiconosVector::getValue(unsigned int row) const
{
  return apply_visitor<GetValue, double>(storage(*this), row);
}

void SiconosVector::setValue(unsigned int row, const double value)
{
  return apply_visitor<SetValue>(storage(*this), row, value);
}

double& SiconosVector::operator()(unsigned int row)
{
  return *apply_visitor<GetRValue, double* >(storage(*this), row);
};

double SiconosVector::operator()(unsigned int row) const
{
  return getValue(row);
}

// //============================================
// // Access (get or set) to blocks of elements
// //============================================

void SiconosVector::setBlock(unsigned int index, const SiconosVector& vIn)
{
  apply_visitor<SetBlock>(storage(vIn), index, storage(*this));
}

void SiconosVector::toBlock(SiconosVector& vOut, unsigned int sizeB, unsigned int startIn, unsigned int startOut) const
{
  apply_visitor<ToBlock>(storage(*this), sizeB, startIn, startOut, storage(vOut));
}

void SiconosVector::addBlock(unsigned int index, const SiconosVector& vIn)
{
  apply_visitor<AddBlock>(storage(vIn), index, storage(*this));
}

void SiconosVector::subBlock(unsigned int index, const SiconosVector& vIn)
{
  apply_visitor<SubBlock>(storage(vIn), index, storage(*this));
}

// //===============
// //  Assignment
// //===============

SiconosVector& SiconosVector::operator = (const SiconosVector& vIn)
{
  if (&vIn == this)
  {
    return *this; // auto-assignment.
  }
  else
  {
    apply_visitor<Copy>(storage(vIn), storage(*this));
    return *this;
  }
}

SiconosVector& SiconosVector::operator = (const BlockVector& vIn)
{
  VectorOfVectors::const_iterator it;
  unsigned int pos = 0;
  for (it = vIn.begin(); it != vIn.end(); ++it)
  {
    setBlock(pos, **it);
    pos += (*it)->size();
  }
  return *this;
}


SiconosVector& SiconosVector::operator = (const DenseVect& d)
{
  if (Type::value(storage(*this)) != Type::DenseVectStorage)
  {
    SiconosVectorException::selfThrow("SiconosVector::operator = DenseVect : current vector is not dense.");
  }
  if (this->size() != d.size())
  {
    SiconosVectorException::selfThrow("SiconosVector::operator = DenseVect : inconsistent size.");
  }
  bindings::copy(d, this->dense());
  return *this;
}

SiconosVector& SiconosVector::operator = (const SparseVect& sp)
{
  if (Type::value(storage(*this)) != Type::SparseVectStorage)
  {
    SiconosVectorException::selfThrow("SiconosVector::operator = SparseVect : current vector is not sparse.");
  }
  if (this->size() != sp.size())
  {
    SiconosVectorException::selfThrow("SiconosVector::operator = SparseVect : inconsistent size.");
  }

  noalias(this->sparse()) = sp;

  return *this;
}

SiconosVector& SiconosVector::operator = (const double* d)
{
  if (Type::value(storage(*this)) == Type::SparseVectStorage)
  {
    SiconosVectorException::selfThrow("SiconosVector::operator = double* : forbidden: the current vector is not dense.");
  }
  bindings::detail::copy(this->size(), d, 1, getArray(), 1);
  return *this;
}

unsigned SiconosVector::copyData(double* data) const
{
  if (Type::value(storage(*this)) == Type::SparseVectStorage)
  {
    SiconosVectorException::selfThrow("SiconosVector::copyData : forbidden: the current vector is not dense.");
  }
  unsigned size = this->size();
  bindings::detail::copy(size, getArray(), 1, data, 1);
  return size;
}


//=================================
// Op. and assignment (+=, -= ... )
//=================================

SiconosVector& SiconosVector::operator += (const SiconosVector& vIn)
{
  apply_visitor<Plus>(storage(vIn), storage(*this));
  return *this;
}
SiconosVector& SiconosVector::operator += (const BlockVector& vIn)
{
  VectorOfVectors::const_iterator it;
  unsigned int pos = 0;
  for (it = vIn.begin(); it != vIn.end(); ++it)
  {
    addBlock(pos, **it);
    pos += (*it)->size();
  }
  return *this;
}

SiconosVector& SiconosVector::operator -= (const SiconosVector& vIn)
{
  apply_visitor<Minus>(storage(vIn), storage(*this));
  return *this;
}

SiconosVector& SiconosVector::operator -= (const BlockVector& vIn)
{
  VectorOfVectors::const_iterator it;
  unsigned int pos = 0;
  for (it = vIn.begin(); it != vIn.end(); ++it)
  {
    subBlock(pos, **it);
    pos += (*it)->size();
  }
  return *this;
}


//===============
// Comparison
//===============

bool operator == (const SiconosVector &m, const SiconosVector &x)
{
  DEBUG_PRINTF("norm = %12.8e \n", (m - x).normInf() );
  DEBUG_PRINTF("std::numeric_limits<double>::epsilon() = %12.8e \n", std::numeric_limits<double>::epsilon() );
  DEBUG_EXPR(std::cout << std::boolalpha << ( (m - x).normInf() <= std::numeric_limits<double>::epsilon()) <<std::endl;);
  return ((m - x).normInf() <= std::numeric_limits<double>::epsilon());
}

//==================
// y = scalar * x
//==================

SiconosVector operator * (const  SiconosVector&m, double d)
{
  SiconosVector tmp = m;
  apply_visitor<Scal>(storage(tmp), d);
  return tmp;
}

SiconosVector operator * (double d, const  SiconosVector&m)
{
  SiconosVector tmp = m;
  apply_visitor<Scal>(storage(tmp), d);
  return tmp;
}

SiconosVector operator / (const SiconosVector &m, double d)
{
  return m * (1.0/d);
}

//====================
//  Vectors addition
//====================

SiconosVector operator + (const  SiconosVector& x, const  SiconosVector& y)
{
  SiconosVector tmp = x;
  tmp += y;
  return tmp;
}

void add(const SiconosVector& x, const SiconosVector& y, SiconosVector& z)
{
  apply_visitor<Copy>(storage(x), storage(z));
  z += y;
}

//======================
//  Vectors subtraction
//======================

SiconosVector operator - (const  SiconosVector& x, const  SiconosVector& y)
{
  SiconosVector tmp = x;
  tmp -= y;
  return tmp;
}

void sub(const SiconosVector& x, const SiconosVector& y, SiconosVector& z)
{
  apply_visitor<Copy>(storage(x), storage(z));
  z -= y;
}




void axpby(double a, const SiconosVector& x, double b, SiconosVector& y)
{
  apply_visitor<Axpby>(storage(x), a, b, storage(y));
}

void axpy(double a, const SiconosVector& x, SiconosVector& y)
{
  apply_visitor<Axpy>(storage(x), a, storage(y));
}


double inner_prod(const SiconosVector &x, const SiconosVector &m)
{
  return apply_visitor<Dot, double>(storage(x), storage(m));
}

//// outer_prod(v,w) = trans(v)*w
SimpleMatrix outer_prod(const SiconosVector &x, const SiconosVector& m)
{
  return apply_visitor<OuterProd, SimpleMatrix >(storage(x), storage(m));
}

void scal(double a, const SiconosVector & x, SiconosVector & y, bool init)
{
  if(init)
  {
    apply_visitor<Copy>(storage(x), storage(y));
    apply_visitor<Scal>(storage(y), a);
  }
  else
  {
    apply_visitor<Axpy>(storage(x), a, storage(y));
  }
}

void subscal(double a, const SiconosVector & x, SiconosVector & y, const Index& coord, bool init)
{
  apply_visitor<Subscal>(storage(x), a, coord, init, storage(y));
}
void cross_product(const SiconosVector& V1, const SiconosVector& V2, SiconosVector& VOUT)
{
  if (V1.size() != 3 || V2.size() != 3 || VOUT.size() != 3)
    SiconosVectorException::selfThrow("SiconosVector::cross_product allowed only with dim 3.");

  double aux = V1.getValue(1) * V2.getValue(2) - V1.getValue(2) * V2.getValue(1);
  VOUT.setValue(0, aux);

  aux = V1.getValue(2) * V2.getValue(0) - V1.getValue(0) * V2.getValue(2);
  VOUT.setValue(1, aux);

  aux = V1.getValue(0) * V2.getValue(1) - V1.getValue(1) * V2.getValue(0);
  VOUT.setValue(2, aux);

}

//

void abs_wise(const SiconosVector& V, SiconosVector& Vabs)
{
  for (unsigned int it = 0; it < V.size(); ++it)
  {
    Vabs.setValue(it, std::abs(V.getValue(it)));
  };
}

//

void getMax(const SiconosVector& V, double& maxvalue, unsigned int& idmax)
{
  maxvalue = V.getValue(0);
  idmax = 0;
  for (unsigned int it = 1; it < V.size(); ++it)
  {
    if (V.getValue(it) > maxvalue)
    {
      maxvalue = V.getValue(it);
      idmax = it;
    };
  };
}

//

void getMin(const SiconosVector& V, double& minvalue, unsigned int& idmin)
{
  minvalue = V.getValue(0);
  idmin = 0;
  for (unsigned int it = 1; it < V.size(); ++it)
  {
    if (V.getValue(it) < minvalue)
    {
      minvalue = V.getValue(it);
      idmin = it;
    };
  };
}


struct exp_op { double operator() (double d) const { return std::exp(d); } };

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


//
/*
SiconosVector abs_wise(const SiconosVector& V){
  SiconosVector Vabs(V.size());
  for (int it = 0; it < V.size(); ++it){
    Vabs.setValue(it,std::abs(V.getValue(it)));
  };
  return Vabs;
}
//
void getMin(const SiconosVector& V, double& minvalue, unsigned int& idmin){
  minvalue = V.getValue(0);
  idmin = 0;
  for (unsigned int it = 1; it < V.size(); ++it){
    if (V.getValue(it) < minvalue){
      minvalue = V.getValue(it);
      idmin = it;
    };
  };
}
*/
void setBlock(const SiconosVector& vIn, SP::SiconosVector vOut, unsigned int sizeB,
              unsigned int startIn, unsigned int startOut)
{
  apply_visitor<ToBlock>(storage(vIn), sizeB, startIn, startOut, storage(*vOut));
}

unsigned int SiconosVector::size(void) const
{
  return apply_visitor<Size, unsigned int>(storage(*this));
}

SiconosVector& operator *= (SiconosVector& v, const double& s)
{

  apply_visitor<Scal>(storage(v), s);
  return v;
}


SiconosVector& operator /= (SiconosVector& v, const double& s)
{
  apply_visitor<Scal>(storage(v), 1.0/s);
  return v;
}
