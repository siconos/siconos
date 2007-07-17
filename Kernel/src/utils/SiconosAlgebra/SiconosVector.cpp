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

#include "SiconosVector.h"

// Default (protected) constructor
SiconosVector::SiconosVector(unsigned int newNum): sizeV(0), num(newNum)
{}

// Basic (protected) constructor
SiconosVector::SiconosVector(unsigned int newNum, unsigned int size): sizeV(size), num(newNum)
{}

Index SiconosVector::getTabIndex() const
{
  SiconosVectorException::selfThrow("SiconosVector::getTabIndex() : not implemented for this type of vector (Simple?) reserved to BlockVectors.");
  // fake to avoid error on warning.
  Index tmp;
  return tmp;
}

const Index * SiconosVector::getTabIndexPtr() const
{
  SiconosVectorException::selfThrow("SiconosVector::getTabIndexPtr() : not implemented for this type of vector (Simple?) reserved to BlockVectors.");
  // fake to avoid error on warning.
  Index * tmp = NULL;
  return tmp;
}

BlockVectIterator SiconosVector::begin()
{
  SiconosVectorException::selfThrow("SiconosVector::begin(): reserved to BlockVectors");
  BlockVectIterator it;
  return it;
};

BlockVectIterator SiconosVector::end()
{
  SiconosVectorException::selfThrow("SiconosVector::end(): reserved to BlockVectors");
  BlockVectIterator it;
  return it;
};

ConstBlockVectIterator SiconosVector::begin() const
{
  SiconosVectorException::selfThrow("SiconosVector::begin(): reserved to BlockVectors");
  ConstBlockVectIterator it;
  return it;
};

ConstBlockVectIterator SiconosVector::end() const
{
  SiconosVectorException::selfThrow("SiconosVector::end(): reserved to BlockVectors");
  ConstBlockVectIterator it;
  return it;
};

//=====================
// vectors comparison
//=====================
const bool isComparableTo(const SiconosVector * v1, const SiconosVector * v2)
{
  // return:
  // - true if one of the vectors is a Simple and if they have the same size
  // - true if both are block but with blocks which are facing each other of the same size.
  // - false in other cases

  if ((!v1->isBlock() || !v2->isBlock()) && (v1->size() == v2->size())) return true;

  const Index * I1 = v1->getTabIndexPtr();
  const Index * I2 = v2->getTabIndexPtr();

  return (*I1 == *I2);
}

void swap(SiconosVector& x, SiconosVector& y)
{
  unsigned int numX = x.getNum();
  unsigned int numY = y.getNum();


  if (numX == numY)
  {
    if (numX == 1) // ublas swap seems to be much more efficient than atlas::swap (bindings) ...
      swap(*x.getDensePtr(), *y.getDensePtr());
    else if (numX == 4)
      swap(*x.getSparsePtr(), *y.getSparsePtr());
    else
      SiconosVectorException::selfThrow("Swap(x,y): not yet implemented for block vectors.");
  }
  else
    SiconosVectorException::selfThrow("Swap(x,y): can not swap two vectors of different types.");
}


