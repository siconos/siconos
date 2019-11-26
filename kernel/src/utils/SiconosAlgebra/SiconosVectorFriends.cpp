/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2019 INRIA.
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

#include "SiconosVectorFriends.hpp"
#include <boost/numeric/ublas/io.hpp>    
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/bindings/ublas/vector_proxy.hpp>
#include <boost/numeric/bindings/blas.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/std/vector.hpp>
#include "SiconosAlgebra.hpp"
#include "BlockVector.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "SiconosMatrixException.hpp"
//#define DEBUG_MESSAGES
#include "debug.h"

namespace siconosBindings = boost::numeric::bindings::blas;

void setBlock(const SiconosVector& vIn, SP::SiconosVector vOut, unsigned int sizeB,
              unsigned int startIn, unsigned int startOut)
{
  unsigned int endOut = startOut + sizeB;
  unsigned int numIn = vIn.num();
  unsigned int numOut = vOut->num();
  assert(vOut->size() >= endOut && "The output vector is too small");
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

//===============
// Comparison
//===============

bool operator == (const SiconosVector &m, const SiconosVector &x)
{
  DEBUG_PRINTF("norm = %12.8e \n", (m - x).normInf() );
  DEBUG_PRINTF("std::numeric_limits<double>::epsilon() = %12.8e \n", std::numeric_limits<double>::epsilon() );
  DEBUG_EXPR(std::cout << std::boolalpha << ( (m - x).normInf() <= std::numeric_limits<double>::epsilon()) <<std::endl;);
  double atol = 1e-14;
  double rtol = std::numeric_limits<double>::epsilon();
  return ((m - x).normInf() <= atol + rtol * x.normInf()) ;
}

//==================
// y = scalar * x
//==================

SiconosVector operator * (const  SiconosVector&m, double d)
{
  unsigned int numM = m.num();

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
  unsigned int numM = m.num();

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
  unsigned int numM = m.num();

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

  unsigned int numX = x.num();
  unsigned int numY = y.num();

  if (numX == numY) // x, y SiconosVector of the same type
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

  else // x, y SiconosVector with y and x of different types
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

  unsigned int numX = x.num();
  unsigned int numY = y.num();
  unsigned int numZ = z.num();

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

  unsigned int numX = x.num();
  unsigned int numY = y.num();

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
  else // x, y SiconosVector with y and x of different types
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

  unsigned int numX = x.num();
  unsigned int numY = y.num();
  unsigned int numZ = z.num();

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

  unsigned int numX = x.num();
  unsigned int numY = y.num();

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

  unsigned int numX = x.num();
  unsigned int numY = y.num();

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

  unsigned int numM = m.num();
  unsigned int numX = x.num();

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
  unsigned int numM = m.num();
  unsigned int numX = x.num();

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

    unsigned int numY = y.num();
    unsigned int numX = x.num();
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

  unsigned int numY = y.num();
  unsigned int numX = x.num();

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

void abs_wise(const SiconosVector& V, SiconosVector& Vabs)
{
  for (unsigned int it = 0; it < V.size(); ++it)
  {
    Vabs.setValue(it, std::abs(V.getValue(it)));
  };
}

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


void prod(const SiconosVector& x, const SiconosMatrix& A, BlockVector& y, bool init)
{

  assert(!(A.isPLUFactorized()) && "A is PLUFactorized in prod !!");

  if (A.size(0) != x.size())
    SiconosMatrixException::selfThrow("prod(x,A,y) error: inconsistent sizes between A and x.");

  if (A.size(1) != y.size())
    SiconosMatrixException::selfThrow("prod(x,A,y) error: inconsistent sizes between A and y.");
  unsigned int pos = 0;
  VectorOfVectors::const_iterator it;
  // For Each subvector of y, y[i], private_prod computes y[i] = subA x, subA being a submatrix of A corresponding to y[i] position.
  for (it = y.begin(); it != y.end(); ++it)
  {
    (*it)->private_prod(createSPtrConstSiconosVector(x), createSPtrConstSiconosMatrix(A), pos, init);
    pos += (*it)->size();
  }
}

//=====================
// convert to an ostream
//=====================

std::ostream& operator<<(std::ostream& os, const SiconosVector& sv)
{
  if (sv._dense)
    os << *sv.vect.Dense;
  else
    os << *sv.vect.Sparse;
  return os;
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


