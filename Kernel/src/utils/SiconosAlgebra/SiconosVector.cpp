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


#include <boost/numeric/ublas/io.hpp>            // for >> 
#include <boost/numeric/ublas/vector_proxy.hpp>  // for project

#include "KernelConfig.h"

#define BIND_FORTRAN_LOWERCASE_UNDERSCORE

#if defined(HAVE_ATLAS)
#include <boost/numeric/bindings/atlas/cblas1.hpp>
#define SB_IS_ATLAS
namespace siconosBindings = boost::numeric::bindings::atlas;
#else
#include <boost/numeric/bindings/blas/blas1.hpp>
#define SB_IS_BLAS
namespace siconosBindings = boost::numeric::bindings::blas;
#endif

#include "SimpleMatrix.hpp"
#include "BlockVector.hpp"
#include "ioVector.hpp"
#include "SiconosVector.hpp"
#include "SiconosAlgebra.hpp"


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
SiconosVector::SiconosVector(unsigned int row, Siconos::UBLAS_TYPE typ)
{
  if (typ == Siconos::SPARSE)
  {
    _dense = false;
    vect.Sparse = new SparseVect(ublas::zero_vector<double>(row));
  }
  else if (typ == Siconos::DENSE)
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
SiconosVector::SiconosVector(unsigned int row, double val, Siconos::UBLAS_TYPE typ)
{
  if (typ == Siconos::SPARSE)
  {
    _dense = false;
    vect.Sparse = new SparseVect(row);
    fill(val);
  }
  else if (typ == Siconos::DENSE)
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
SiconosVector::SiconosVector(const SiconosVector &svect) : boost::enable_shared_from_this<SiconosVector>()
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
SiconosVector::SiconosVector(const BlockVector & vIn) : boost::enable_shared_from_this<SiconosVector>()
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
    ioVector::read(file, "ascii", *this);
  }
  else
  {
    ioVector::read(file, "binary", *this);
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
  else if (ask<IsSparse>(v1) && ask<IsSparse>(v2))
  {
    _dense = false;
    vect.Sparse = new SparseVect(size1 + v2.size());
  }
  setBlock(0, v1);
  setBlock(size1, v2);
}

SiconosVector::~SiconosVector()
{
  if (_dense)
    delete(vect.Dense);
  else delete(vect.Sparse);
}


// =================================================
//        get Ublas component (dense or sparse)
// =================================================

const DenseVect SiconosVector::getDense(unsigned int) const
{
  if (!_dense)
    SiconosVectorException::selfThrow("SiconosVector::getDense(unsigned int row, unsigned int col) : the current vector is not a Dense vector");

  return *vect.Dense;
}

const SparseVect SiconosVector::getSparse(unsigned int)const
{

  if (_dense)
    SiconosVectorException::selfThrow("SiconosVector::getSparse(unsigned int row, unsigned int col) : the current vector is not a Sparse vector");

  return *vect.Sparse;
}

SparseVect* SiconosVector::sparse(unsigned int)const
{

  if (_dense)
    SiconosVectorException::selfThrow("SiconosVector::sparse(unsigned int row, unsigned int col) : the current vector is not a Sparse vector");

  return vect.Sparse;
}

double* SiconosVector::getArray(unsigned int) const
{
  if (!_dense)
    SiconosVectorException::selfThrow("SiconosVector::getArray() : not yet implemented for sparse vector.");

  assert(vect.Dense);

  return &(((*vect.Dense).data())[0]);
}

// ===========================
//       fill vector
// ===========================

void SiconosVector::zero()
{
  if (_dense)
    siconosBindings::scal(0.0, *vect.Dense);

  else
  {
    assert(vect.Sparse);
    *vect.Sparse *= 0.0;
  }

}

void SiconosVector::setVector(unsigned int, const SiconosVector& newV)
{
  if (newV.size() != size())
    SiconosVectorException::selfThrow("SiconosVector::setVector(num,v), unconsistent sizes.");

  *this = newV ;
}

void SiconosVector::fill(double value)
{
  if (!_dense)
  {
    for (unsigned int i = 0; i < (vect.Sparse)->size(); ++i)
      (vect.Sparse)->push_back(i, value);
  }
  else
  {
    // only atlas provides set functions
    // they are called catlas_*set for a reason
#if defined(SB_IS_ATLAS)
    siconosBindings::set(value, *vect.Dense);
#else
    for (unsigned int i = 0; i < (vect.Dense)->size(); ++i)
      (*vect.Dense)(i) = value;
#endif
  }

}

//=======================
// set vector dimension
//=======================

void SiconosVector::resize(unsigned int n, bool preserve)
{
  if (_dense)
    (vect.Dense)->resize(n, preserve);
  else
    (vect.Sparse)->resize(n, preserve);
}

//=======================
//       get norm
//=======================

double SiconosVector::normInf() const
{
  if (_dense)
    return norm_inf(*vect.Dense);
  else //if(num==4)
    return norm_inf(*vect.Sparse);
}

double SiconosVector::norm2() const
{
  if (_dense)
    return ublas::norm_2(*vect.Dense);
  else //if(num==4)
    return ublas::norm_2(*vect.Sparse);
}
//======================================
// get sum of all elements of the vector
//=====================================
double SiconosVector::sum() const
{
  if (_dense)
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
  if (_dense)
    std::cout << *vect.Dense << std::endl;
  else if (vect.Sparse)
    std::cout << *vect.Sparse << std::endl;
}

//============================
// Convert vector to a string
//============================

const std::string SiconosVector::toString() const
{
  std::stringstream sstr;
  std::string s;
  if (_dense)
    sstr << *vect.Dense;
  else
    sstr << *vect.Sparse;
  sstr >> s;
  s = s.substr(4, s.size() - 5); // Remove "[size](" at the beginning of the string
  std::string::size_type pos;
  while ((pos = s.find(",")) != std::string::npos) // Replace "," by " " in the string
    s[pos] = ' ';
  return s;
}

//=============================
// Elements access (get or set)
//=============================

double SiconosVector::getValue(unsigned int row) const
{
  if (row >= size())
    SiconosVectorException::selfThrow("SiconosVector::getValue(index) : Index out of range");

  if (_dense)
    return (*vect.Dense)(row);
  else
    return (*vect.Sparse)(row);
}

void SiconosVector::setValue(unsigned int row, double value)
{
  if (row >= size())
    SiconosVectorException::selfThrow("SiconosVector::setValue(index, value) : Index out of range");
  if (_dense)
    (*vect.Dense)(row) = value ;
  else
    (*vect.Sparse)(row) = value;
}

double& SiconosVector::operator()(unsigned int row)
{
  if (row >= size())
    SiconosVectorException::selfThrow("SiconosVector::operator ( index ): Index out of range");

  if (_dense)
    return (*vect.Dense)(row);
  else
    return (*vect.Sparse)(row).ref();
}

double SiconosVector::operator()(unsigned int row) const
{
  if (row >= size())
    SiconosVectorException::selfThrow("SiconosVector::operator ( index ): Index out of range");

  if (_dense)
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
  if (&vIn == this)
    SiconosVectorException::selfThrow("SiconosVector::this->setBlock(pos,vIn): vIn = this.");

  if (index > size())
    SiconosVectorException::selfThrow("SiconosVector::setBlock : invalid ranges");

  unsigned int end = vIn.size() + index;
  if (end > size())
    SiconosVectorException::selfThrow("SiconosVector::setBlock : invalid ranges");

  unsigned int numVin = vIn.getNum();
  if (numVin != getNum()) SiconosVectorException::selfThrow("SiconosVector::setBlock: inconsistent types.");

  if (_dense)
    noalias(ublas::subrange(*vect.Dense, index, end)) = *vIn.dense();
  else
    noalias(ublas::subrange(*vect.Sparse, index, end)) = *vIn.sparse();
}

void SiconosVector::toBlock(SiconosVector& vOut, unsigned int sizeB, unsigned int startIn, unsigned int startOut) const
{
  // To copy a subBlock of the vector (from position startIn to startIn+sizeB) into vOut (from pos. startOut to startOut+sizeB).
  // Check dim ...
  unsigned int sizeIn = size();
  unsigned int sizeOut = vOut.size();

  if (startIn >= sizeIn)
    SiconosVectorException::selfThrow("vector toBlock(v1,v2,...): start position in input vector is out of range.");

  if (startOut >= sizeOut)
    SiconosVectorException::selfThrow("vector toBlock(v1,v2,...): start position in output vector is out of range.");

  unsigned int endIn = startIn + sizeB;
  unsigned int endOut = startOut + sizeB;

  if (endIn > sizeIn)
    SiconosVectorException::selfThrow("vector toBlock(v1,v2,...): end position in input vector is out of range.");
  if (endOut > sizeOut)
    SiconosVectorException::selfThrow("vector toBlock(v1,v2,...): end position in output vector is out of range.");

  unsigned int numIn = getNum();
  unsigned int numOut = vOut.getNum();

  if (numIn == numOut)
  {
    if (numIn == 1) // vIn / vOut are Dense
      noalias(ublas::subrange(*vOut.dense(), startOut, endOut)) = ublas::subrange(*vect.Dense, startIn, startIn + sizeB);
    else // if(numIn == 4)// vIn / vOut are Sparse
      noalias(ublas::subrange(*vOut.sparse(), startOut, endOut)) = ublas::subrange(*vect.Sparse, startIn, startIn + sizeB);
  }
  else // vIn and vout of different types ...
  {
    if (numIn == 1) // vIn Dense
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

  if (&vIn == this)
    SiconosVectorException::selfThrow("SiconosVector::this->addBlock(pos,vIn): vIn = this.");

  unsigned int end = vIn.size();
  if ((index + end) > size()) SiconosVectorException::selfThrow("SiconosVector::addBlock : invalid ranges");

  unsigned int numVin = vIn.getNum();

  if (numVin != getNum()) SiconosVectorException::selfThrow("SiconosVector::addBlock : inconsistent types.");

  if (_dense)
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
  if ((index + end) > size()) SiconosVectorException::selfThrow("SiconosVector::subBlock : invalid ranges");

  unsigned int numVin = vIn.getNum();
  if (numVin != getNum()) SiconosVectorException::selfThrow("SiconosVector::subBlock : inconsistent types.");

  if (_dense)
    noalias(ublas::subrange(*vect.Dense, index, index + end)) -= *vIn.dense();
  else
    noalias(ublas::subrange(*vect.Sparse, index, index + end)) -= *vIn.sparse();
}

//===============
//  Assignment
//===============

SiconosVector& SiconosVector::operator = (const SiconosVector& vIn)
{
  if (&vIn == this) return *this; // auto-assignment.

  if (size() != vIn.size())
    SiconosVectorException::selfThrow("SiconosVector::operator = failed: inconsistent sizes.");

  unsigned int vInNum = vIn.getNum();

  {
    switch (getNum())
    {
    case 1:
      switch (vInNum)
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
      if (vInNum == 4)
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
  for (it = vIn.begin(); it != vIn.end(); ++it)
  {
    setBlock(pos, **it);
    pos += (*it)->size();
  }
  return *this;
}


SiconosVector& SiconosVector::operator = (const DenseVect& d)
{
  if (!_dense)
    SiconosVectorException::selfThrow("SiconosVector::operator = DenseVect : forbidden: the current vector is not dense.");
  if (d.size() != size())
    SiconosVectorException::selfThrow("SiconosVector::operator = DenseVect : inconsistent size.");

  siconosBindings::copy(d, *vect.Dense);
  return *this;
}

SiconosVector& SiconosVector::operator = (const SparseVect& sp)
{
  if (_dense)
    SiconosVectorException::selfThrow("SiconosVector::operator = SparseVect : current vector is not sparse.");
  if (sp.size() != size())
    SiconosVectorException::selfThrow("SiconosVector::operator = SparseVect : inconsistent size.");

  noalias(*vect.Sparse) = sp;

  return *this;
}

//=================================
// Op. and assignment (+=, -= ... )
//=================================

SiconosVector& SiconosVector::operator += (const SiconosVector& vIn)
{
  if (&vIn == this) // alias
  {
    // Note: using this *= 2.0 is much more time-consuming.
    switch (getNum())
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

  unsigned int vInNum = vIn.getNum();
  {
    switch (getNum())
    {
    case 1:
      switch (vInNum)
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
      if (vInNum == 4)
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
  for (it = vIn.begin(); it != vIn.end(); ++it)
  {
    addBlock(pos, **it);
    pos += (*it)->size();
  }
  return *this;
}

SiconosVector& SiconosVector::operator -= (const SiconosVector& vIn)
{
  if (&vIn == this)
  {
    this->zero();
    return *this;
  }

  unsigned int vInNum = vIn.getNum();
  {
    switch (getNum())
    {
    case 1:
      switch (vInNum)
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
      if (vInNum == 4)
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
  return ((m - x).norm2() < tolerance);
}

//==================
// y = scalar * x
//==================

SiconosVector operator * (const  SiconosVector&m, double d)
{
  unsigned int numM = m.getNum();

  if (numM == 1)
  {
    // Copy m into p and call siconosBindings::scal(d,p), p = d*p.
    DenseVect p = *m.dense();
    siconosBindings::scal(d, p);
    return p;
  }
  else// if(numM==4)
  {
    return (SparseVect)(*m.sparse() * d);
  }
}

SiconosVector operator * (double d, const  SiconosVector&m)
{
  unsigned int numM = m.getNum();

  if (numM == 1)
  {
    // Copy m into p and call siconosBindings::scal(d,p), p = d*p.
    DenseVect p = *m.dense();
    siconosBindings::scal(d, p);
    return p;
  }
  else// if(numM==4)
  {
    return (SparseVect)(*m.sparse() * d);
  }
}

SiconosVector operator / (const SiconosVector &m, double d)
{
  unsigned int numM = m.getNum();

  if (numM == 1)
  {
    DenseVect p = *m.dense();
    siconosBindings::scal((1.0 / d), p);
    return p;
  }

  else// if(numM==4){
    return (SparseVect)(*m.sparse() / d);
}

//====================
//  Vectors addition
//====================

SiconosVector operator + (const  SiconosVector& x, const  SiconosVector& y)
{
  if (x.size() != y.size())
    SiconosVectorException::selfThrow("SiconosVector, x + y: inconsistent sizes");

  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();

  if (numX == numY)  // x, y SiconosVector of the same type
  {
    if (numX == 1)
    {
      //    siconosBindings::xpy(*x.dense(),p);
      //    return p;
      return (DenseVect)(*x.dense() + *y.dense());
    }
    else
      return (SparseVect)(*x.sparse() + *y.sparse());
  }

  else  // x, y SiconosVector with y and x of different types
  {
    if (numX == 1)
      return (DenseVect)(*x.dense() + *y.sparse());
    else
      return (DenseVect)(*x.sparse() + *y.dense());
  }

}

void add(const SiconosVector& x, const SiconosVector& y, SiconosVector& z)
{
  // Computes z = x + y in an "optimized" way (in comparison with operator +)

  if (x.size() != y.size() || x.size() != z.size())
    SiconosVectorException::selfThrow("add(x,y,z): inconsistent sizes");

  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();
  unsigned int numZ = z.getNum();

  if (&z == &x) // x, and z are the same object.
  {
    z += y;
  }
  else if (&z == &y) // y and z are the same object, different from x
  {
    z += x;
  }
  else // No common memory between x,y and z
  {

    if (numZ != 0) // z is a SiconosVector
    {
      if (numX == numY && numX != 0) // x, y SiconosVector of the same type
      {
        if (numX == 1)
        {
          if (numZ != 1)
            SiconosVectorException::selfThrow("SiconosVector addition, add(x,y,z) failed - Addition of two dense vectors into a sparse.");
          noalias(*z.dense()) = *x.dense() + *y.dense() ;
        }
        else
        {
          if (numZ == 1)
            noalias(*z.dense()) = *x.sparse() + *y.sparse() ;
          else
            noalias(*z.sparse()) = *x.sparse() + *y.sparse() ;
        }
      }
      else if (numX != 0 && numY != 0) // x and y of different types => z must be dense.
      {
        if (numZ != 1)
          SiconosVectorException::selfThrow("SiconosVector addition, add(x,y,z) failed - z can not be sparse.");
        if (numX == 1)
          noalias(*z.dense()) = *x.dense() + *y.sparse();
        else
          noalias(*z.dense()) = *x.sparse() + *y.dense() ;
      }
    }
  }
}

//======================
//  Vectors subtraction
//======================

SiconosVector operator - (const  SiconosVector& x, const  SiconosVector& y)
{
  if (x.size() != y.size())
    SiconosVectorException::selfThrow("SiconosVector, x - y: inconsistent sizes");

  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();

  if (numX == numY) // x, y SiconosVector of the same type
  {
    if (numX == 1)
    {
      //    siconosBindings::xpy(*x.dense(),p);
      //    return p;
      return (DenseVect)(*x.dense() - *y.dense());
    }
    else
      return (SparseVect)(*x.sparse() - *y.sparse());
  }
  else  // x, y SiconosVector with y and x of different types
  {
    if (numX == 1)
      return (DenseVect)(*x.dense() - *y.sparse());
    else
      return (DenseVect)(*x.sparse() - *y.dense());
  }
}

void sub(const SiconosVector& x, const SiconosVector& y, SiconosVector& z)
{
  // Computes z = x - y in an "optimized" way (in comparison with operator +)

  if (x.size() != y.size() || x.size() != z.size())
    SiconosVectorException::selfThrow("sub(x,y,z): inconsistent sizes");

  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();
  unsigned int numZ = z.getNum();

  if (&z == &x) // x and z are the same object.
  {
    z -= y;
  }
  else if (&z == &y) // y and z are the same object
  {
    {
      if (numX == 1)
      {
        if (numZ != 1)
          SiconosVectorException::selfThrow("SiconosVector subtraction, sub(x,y,z) failed - Subtraction of two dense vectors into a sparse.");
        *z.dense() = *x.dense() - *y.dense() ;
      }
      else
      {
        if (numZ == 1)
          *z.dense() = *x.sparse() - *y.dense() ;
        else
          *z.sparse() = *x.sparse() - *y.sparse() ;
      }
    }
  }
  else // No common memory between x or y and z
  {

    if (numZ != 0) // z is a SiconosVector
    {
      if (numX == numY && numX != 0) // x, y SiconosVector of the same type
      {
        if (numX == 1)
        {
          if (numZ != 1)
            SiconosVectorException::selfThrow("SiconosVector addition, sub(x,y,z) failed - Addition of two dense vectors into a sparse.");
          noalias(*z.dense()) = *x.dense() - *y.dense() ;
        }
        else
        {
          if (numZ == 1)
            noalias(*z.dense()) = *x.sparse() - *y.sparse() ;
          else
            noalias(*z.sparse()) = *x.sparse() - *y.sparse() ;
        }
      }
      else if (numX != 0 && numY != 0) // x and y of different types => z must be dense.
      {
        if (numZ != 1)
          SiconosVectorException::selfThrow("SiconosVector addition, sub(x,y,z) failed - z can not be sparse.");
        if (numX == 1)
          noalias(*z.dense()) = *x.dense() - *y.sparse();
        else
          noalias(*z.dense()) = *x.sparse() - *y.dense() ;
      }
    }
  }
}

void axpby(double a, const SiconosVector& x, double b, SiconosVector& y)
{
  // Computes y = ax + by

  if (x.size() != y.size())
    SiconosVectorException::selfThrow("axpby(x,y,z): inconsistent sizes");

  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();

  if (numX == numY) // x and y of the same type
  {
    if (numX == 1) // all dense
    {
      siconosBindings::scal(b, *y.dense());
      siconosBindings::axpy(a, *x.dense(), *y.dense());
    }
    else // all sparse
    {
      *y.sparse() *= b;
      if (&y != &x)
        noalias(*y.sparse()) += a**x.sparse();
      else
        *y.sparse() += a**x.sparse();
    }
  }

  else // x and y of different types
  {
    y *= b;
    {
      if (numX == 1)
        *y.sparse() += a**x.dense();
      else
        *y.dense() +=  a**x.sparse();
    }
  }
}

void axpy(double a, const SiconosVector& x, SiconosVector& y)
{
  // Computes y = ax + y

  if (x.size() != y.size())
    SiconosVectorException::selfThrow("axpy(x,y,z): inconsistent sizes");

  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();

  if (numX == numY) // x and y of the same type
  {
    if (numX == 1) // all dense
      siconosBindings::axpy(a, *x.dense(), *y.dense());

    else // all sparse
    {
      if (&y != &x)
        noalias(*y.sparse()) += a**x.sparse();
      else
        *y.sparse() += a**x.sparse();
    }
  }

  else // x and y of different types
  {
    {
      if (numX == 1)
        *y.sparse() += a**x.dense();
      else
        *y.dense() +=  a**x.sparse();
    }
  }
}

double inner_prod(const SiconosVector &x, const SiconosVector &m)
{
  if (x.size() != m.size())
    SiconosVectorException::selfThrow("inner_prod: inconsistent sizes");

  unsigned int numM = m.getNum();
  unsigned int numX = x.getNum();

  if (numX == numM)
  {
    if (numM == 1)
      return siconosBindings::dot(*x.dense(), *m.dense());
    else
      return inner_prod(*x.sparse(), *m.sparse());
  }
  else if (numM == 1)
    return inner_prod(*x.sparse(), *m.dense());
  else
    return inner_prod(*x.dense(), *m.sparse());
}

// outer_prod(v,w) = trans(v)*w
SimpleMatrix outer_prod(const SiconosVector &x, const SiconosVector& m)
{
  unsigned int numM = m.getNum();
  unsigned int numX = x.getNum();

  if (numM == 1)
  {
    if (numX == 1)
      return (DenseMat)(outer_prod(*x.dense(), *m.dense()));

    else// if(numX == 4)
      return (DenseMat)(outer_prod(*x.sparse(), *m.dense()));
  }
  else // if(numM == 4)
  {
    if (numX == 1)
      return (DenseMat)(outer_prod(*x.dense(), *m.sparse()));

    else //if(numX == 4)
      return (DenseMat)(outer_prod(*x.sparse(), *m.sparse()));
  }
}

void scal(double a, const SiconosVector & x, SiconosVector & y, bool init)
{
  // To compute y = a *x (init = true) or y += a*x (init = false)

  if (&x == &y)
  {
    if (init)
      y *= a;
    else
    {
      y *= (1.0 + a);
    }
  }
  else
  {
    unsigned int sizeX = x.size();
    unsigned int sizeY = y.size();

    if (sizeX != sizeY)
      SiconosVectorException::selfThrow("scal(a,SiconosVector,SiconosVector) failed, sizes are not consistent.");

    unsigned int numY = y.getNum();
    unsigned int numX = x.getNum();
    if (numX == numY)
    {

      if (numX == 1) // ie if both are Dense
      {
        if (init)
          //siconosBindings::axpby(a,*x.dense(),0.0,*y.dense());
          noalias(*y.dense()) = a * *x.dense();
        else
          noalias(*y.dense()) += a * *x.dense();
      }
      else  // if both are sparse
      {
        if (init)
          noalias(*y.sparse()) = a**x.sparse();
        else
          noalias(*y.sparse()) += a**x.sparse();
      }
    }
    else
    {
      if (numY == 0 || numX == 0) // if y or x is block
      {
        if (init)
        {
          y = x;
          y *= a;
        }
        else
        {
          SiconosVector tmp(x);
          tmp *= a;
          y += tmp;
        }
      }
      else
      {
        if (numY == 1) // if y is dense
        {
          if (init)
            noalias(*y.dense()) = a**x.sparse();
          else
            noalias(*y.dense()) += a**x.sparse();

        }
        else
          SiconosVectorException::selfThrow("SiconosVector::scal(a,dense,sparse) not allowed.");
      }
    }
  }
}

void subscal(double a, const SiconosVector & x, SiconosVector & y, const Index& coord, bool init)
{
  // To compute sub_y = a *sub_x (init = true) or sub_y += a*sub_x (init = false)
  // Coord  = [r0x r1x r0y r1y];
  // subX is the sub-vector of x, for row numbers between r0x and r1x-1.
  // The same for y with riy.


  // Check dimensions
  unsigned int dimX = coord[1] - coord[0];
  unsigned int dimY = coord[3] - coord[2];
  if (dimY != dimX)
    SiconosVectorException::selfThrow("subscal(a,x,y,...) error: inconsistent sizes between (sub)x and (sub)y.");
  if (dimY > y.size() || dimX > x.size())
    SiconosVectorException::selfThrow("subscal(a,x,y,...) error: input index too large.");

  unsigned int numY = y.getNum();
  unsigned int numX = x.getNum();

  if (&x == &y) // if x and y are the same object
  {
    if (numX == 1) // Dense
    {
      ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[2], coord[3]));
      if (coord[0] == coord[2])
      {
        if (init)
          subY *= a;
        else
          subY *= (1.0 + a);
      }
      else
      {
        ublas::vector_range<DenseVect> subX(*x.dense(), ublas::range(coord[0], coord[1]));
        if (init)
          subY = a * subX;
        else
          subY += a * subX;
      }
    }
    else //if (numX == 4) // Sparse
    {
      ublas::vector_range<SparseVect> subY(*y.sparse(), ublas::range(coord[2], coord[3]));
      if (coord[0] == coord[2])
      {
        if (init)
          subY *= a;
        else
          subY *= (1.0 + a);
      }
      else
      {
        ublas::vector_range<SparseVect> subX(*x.sparse(), ublas::range(coord[0], coord[1]));
        if (init)
          subY = a * subX;
        else
          subY += a * subX;
      }
    }
  }
  else
  {
    if (numX == numY)
    {
      if (numX == 1) // ie if both are Dense
      {
        ublas::vector_range<DenseVect> subX(*x.dense(), ublas::range(coord[0], coord[1]));
        ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[2], coord[3]));

        if (init)
          noalias(subY) = a * subX;
        else
          noalias(subY) += a * subX;
      }
      else  // if both are sparse
      {
        ublas::vector_range<SparseVect> subX(*x.sparse(), ublas::range(coord[0], coord[1]));
        ublas::vector_range<SparseVect> subY(*y.sparse(), ublas::range(coord[2], coord[3]));

        if (init)
          noalias(subY) = a * subX;
        else
          noalias(subY) += a * subX;
      }
    }
    else // x and y of different types ...
    {
      if (numY == 1) // y dense, x sparse
      {
        ublas::vector_range<DenseVect> subY(*y.dense(), ublas::range(coord[2], coord[3]));
        ublas::vector_range<SparseVect> subX(*x.sparse(), ublas::range(coord[0], coord[1]));

        if (init)
          noalias(subY) = a * subX;
        else
          noalias(subY) += a * subX;
      }
      else // y sparse, x dense => fails
        SiconosVectorException::selfThrow("SiconosVector::subscal(a,dense,sparse) not allowed.");
    }
  }
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
void getMax(const SiconosVector& V, double& maxvalue, unsigned int& idmax){
  maxvalue = V.getValue(0);
  idmax = 0;
  for (unsigned int it = 1; it < V.size(); ++it){
    if (V.getValue(it) > maxvalue){
    maxvalue = V.getValue(it);
    idmax = it;
    };
  };
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
void setBlock(const SiconosVector& vIn, SP::SiconosVector vOut, unsigned int sizeB, unsigned int startIn, unsigned int startOut)
{
  unsigned int endOut = startOut + sizeB;
  unsigned int numIn = vIn.getNum();
  unsigned int numOut = vOut->getNum();
  if (numIn == numOut)
  {
    if (numIn == 1) // vIn / vOut are Dense
      noalias(ublas::subrange(*vOut->dense(), startOut, endOut)) = ublas::subrange(*vIn.dense(), startIn, startIn + sizeB);
    else // if(numIn == 4)// vIn / vOut are Sparse
      noalias(ublas::subrange(*vOut->sparse(), startOut, endOut)) = ublas::subrange(*vIn.sparse(), startIn, startIn + sizeB);
  }
  else // vIn and vout of different types ...
  {
    if (numIn == 1) // vIn Dense
      noalias(ublas::subrange(*vOut->sparse(), startOut, endOut)) = ublas::subrange(*vIn.dense(), startIn, startIn + sizeB);
    else // if(numIn == 4)// vIn Sparse
      noalias(ublas::subrange(*vOut->dense(), startOut, endOut)) = ublas::subrange(*vIn.sparse(), startIn, startIn + sizeB);
  }
}

unsigned int SiconosVector::size(void) const
{
  if (!_dense)
  {
    return (vect.Sparse->size());
  }
  else
  {
    return (vect.Dense->size());
  }
}

SiconosVector& operator *= (SiconosVector& v, const double& s)
{
  if (v._dense)
    *v.dense() *= s;
  else
    *v.sparse() *= s;
  return v;
}


SiconosVector& operator /= (SiconosVector& v, const double& s)
{
  if (v._dense)
    *v.dense() /= s;
  else
    *v.sparse() /= s;
  return v;
}

