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
#include "SiconosMatrix.hpp"

// Constructor with the type-number
SiconosMatrix::SiconosMatrix(unsigned int newNum): dimRow(0), dimCol(0), num(newNum)
{}

// Constructor with the dimensions and the type-number
SiconosMatrix::SiconosMatrix(unsigned int newNum, unsigned int row, unsigned int col): dimRow(row), dimCol(col), num(newNum)
{}

BlockIterator1 SiconosMatrix::begin()
{
  SiconosMatrixException::selfThrow("SiconosMatrix::begin(): reserved to BlockMatrix");
  BlockIterator1 it;
  return it;
};

BlockIterator1 SiconosMatrix::end()
{
  SiconosMatrixException::selfThrow("SiconosMatrix::end(): reserved to BlockMatrix");
  BlockIterator1 it;
  return it;
};

ConstBlockIterator1 SiconosMatrix::begin() const
{
  SiconosMatrixException::selfThrow("SiconosMatrix::begin(): reserved to BlockMatrix");
  ConstBlockIterator1 it;
  return it;
};

ConstBlockIterator1 SiconosMatrix::end() const
{
  SiconosMatrixException::selfThrow("SiconosMatrix::end(): reserved to BlockMatrix");
  ConstBlockIterator1 it;
  return it;
};

const SP::Index SiconosMatrix::tabRow() const
{
  SiconosMatrixException::selfThrow("SiconosMatrix::tabRow() : not implemented for this type of matrix (Simple?) reserved to BlockMatrix.");
  // fake to avoid error on warning.
  return SP::Index();
}

const SP::Index SiconosMatrix::tabCol() const
{
  SiconosMatrixException::selfThrow("SiconosMatrix::tabCol() : not implemented for this type of matrix (Simple?) reserved to BlockMatrix.");
  // fake to avoid error on warning.
  return SP::Index();
}

//=====================
// matrices comparison
//=====================
bool isComparableTo(const  SiconosMatrix& m1, const  SiconosMatrix& m2)
{
  // return:
  // - true if one of the matrices is a Simple and if they have the same dimensions.
  // - true if both are block but with blocks which are facing each other of the same size.
  // - false in other cases

  if ((!m1.isBlock() || !m2.isBlock()) && (m1.size(0) == m2.size(0)) && (m1.size(1) == m2.size(1)))
    return true;

  const SP::Index I1R = m1.tabRow();
  const SP::Index I2R = m2.tabRow();
  const SP::Index I1C = m1.tabCol();
  const SP::Index I2C = m2.tabCol();

  return ((*I1R == *I2R) && (*I1C == *I2C));
}

