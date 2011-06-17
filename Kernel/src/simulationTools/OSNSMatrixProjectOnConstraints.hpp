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
/*! \file OSNSMatrix.h
  Specific storage for matrices used in OneStepNSProblem
*/

#ifndef OSNSMPROJECTONCONSTRAINT_H
#define OSNSMPROJECTONCONSTRAINT_H

#include "OSNSMatrix.hpp"


/** Interface to some specific storage types for matrices used in
 * OneStepNSProblem
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) 05/02/2010
 *
 * This class is used to define an interface for various storage used
 * for matrices in OneStepNSProblem. Its aim is to fill the
 * Numerics structure NumericsMatrix, required in many XXX_problem
 * structures of Numerics as input argument for drivers. \n
 *
 * The idea is to remove all matrix storage management problems from
 * OSNS classes (LCP ...) and to leave it into this class. \n
 *
 * Two main functions:
 * - fill(indexSet, unitaryBlocks): fill the matrix using a list of
 *   "active" UnitaryRelation, in indexSet, and a
 *   MapOfMapOfUnitaryMatrices, unitaryBlocks, which determines
 *   which UR are connected or not (ie have common DynamicalSystem).
 *   - convert(): fill the NumericsMatrix structure (indeed only
 *   pointers links to the components of the present class)
 *
 * Note that OSNSMatrix are square.
 *
 *  For example, if in a LCP, constraints of interest are
 *  indexSet={UR2,UR3,UR8,UR12}, whith common DynamicalSystem between
 *  2 and 3, 2 and 8 and 8 and 12.

 *  unitaryBlocks contains matrices for all (URi,URj) which have
 *  common DS, for (URi,URj) in I0, the set of all UnitaryRelation.

 *  (for details on how unitaryBlocks is computed see OneStepNSProblem.h).
 *
 * We denote unitaryBlocks[URi][URj] = mij \n Then, a call to
 * fill(indexSet, unitaryBlock) results in a matrix which looks like:
 *
 * \f{eqnarray*}
 M=\left\lbrace\begin{array}{cccc}
 m22 & m23 & m28 &  0 \\
 m32 & m33 & 0   &  0 \\
 0  &  0  & m88 & m812 \\
 0  &  0  & m128& m1212
 \end{array}\right.
 \f}
 *
 *
 * Note: at the time the available storage types are:
 *
 *  - full matrix in a SiconosMatrix (storageType = 0). In this case,
 *  for each call to fill(), the SiconosMatrix M is resized
 *  according  to the sizes of the UR present in indexSet and then
 *  all the required unitaryBlocks mij are COPIED into M.
 *
 *  - Sparse Block Storage (storageType = 1): corresponds to
 *  SparseBlockStructuredMatrix structure of Numerics. Only non-null
 *  unitaryBlocks are saved in the matrix M and there is no copy of
 *  sub-unitaryBlocks, only links thanks to pointers.
 *
 */


class OSNSMatrixProjectOnConstraints : public OSNSMatrix
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(OSNSMatrixProjectOnConstraints);

  virtual void updateSizeAndPositions(unsigned int&, SP::UnitaryRelationsGraph);
public:


  /** Constructor with dimRow and DimColumn of the matrix
      \param n and m sizes of the rectangle matrix
      \param stor storage type (0:dense, 1:sparse unitaryBlock)
  */
  OSNSMatrixProjectOnConstraints(unsigned int, unsigned int, int);



  /** destructor
   */
  virtual ~OSNSMatrixProjectOnConstraints();


  /** fill the current class using an index set and a map of unitaryBlocks
      \param UnitaryRelationsGraph*, the index set of the active constraints
  */
  virtual void fill(SP::UnitaryRelationsGraph, bool updateSize = true);

  virtual unsigned int getPositionOfUnitaryBlock(SP::UnitaryRelation) const;
};

DEFINE_SPTR(OSNSMatrixProjectOnConstraints);

#endif
