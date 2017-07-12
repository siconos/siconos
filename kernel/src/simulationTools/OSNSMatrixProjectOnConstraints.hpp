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
/*! \file OSNSMatrixProjectOnConstraints.hpp
  Specific storage for matrices used in OneStepNSProblem with a projection scheme
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
 * - fill(indexSet, interactionBlocks): fill the matrix using a list of
 *   "active" Interaction, in indexSet, and a
 *   MapOfMapOfInteractionMatrices, interactionBlocks, which determines
 *   which Interaction are connected or not (ie have common DynamicalSystem).
 *   - convert(): fill the NumericsMatrix structure (indeed only
 *   pointers links to the components of the present class)
 *
 * Note that OSNSMatrix are square.
 *
 *  For example, if in a LCP, constraints of interest are
 *  indexSet={inter2,inter3,inter8,inter12}, whith common DynamicalSystem between
 *  2 and 3, 2 and 8 and 8 and 12.

 *  interactionBlocks contains matrices for all (interi,interj) which have
 *  common DS, for (interi,interj) in I0, the set of all Interaction.

 *  (for details on how interactionBlocks is computed see OneStepNSProblem.h).
 *
 * We denote interactionBlocks[interi][interj] = mij \n Then, a call to
 * fill(indexSet, interactionBlock) results in a matrix which looks like:
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
 *  - full matrix in a SiconosMatrix (_storageType = 0). In this case,
 *  for each call to fill(), the SiconosMatrix M is resized
 *  according  to the sizes of the Interaction present in indexSet and then
 *  all the required interactionBlocks mij are COPIED into M.
 *
 *  - Sparse Block Storage (_storageType = 1): corresponds to
 *  SparseBlockStructuredMatrix structure of Numerics. Only non-null
 *  interactionBlocks are saved in the matrix M and there is no copy of
 *  sub-interactionBlocks, only links thanks to pointers.
 *
 */


class OSNSMatrixProjectOnConstraints : public OSNSMatrix
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(OSNSMatrixProjectOnConstraints);

  /* default constructor
   */
  OSNSMatrixProjectOnConstraints() {};

  using OSNSMatrix::updateSizeAndPositions;
  virtual unsigned updateSizeAndPositions(InteractionsGraph& indexSet);

public:


  /** Constructor with dimRow and DimColumn of the matrix
   * \param n row size of the rectangle matrix
   * \param m column size of the rectangle matrix
   * \param stor storage type (0:dense, 1:sparse interactionBlock)
   */
  OSNSMatrixProjectOnConstraints(unsigned int n, unsigned int m, int stor);

  /** compute the size of the vector to project for a given Interaction.
   * \param inter the corresponding interaction
   * \return  unsigned int
   */
  unsigned int computeSizeForProjection(SP::Interaction inter);


  /** destructor
   */
  virtual ~OSNSMatrixProjectOnConstraints();


  /** fill the current class using an index set and a map of interactionBlocks
      \param indexSet the index set of the active constraints
      \param update if true update the size and position
  */
  void fillW(InteractionsGraph& indexSet, bool update = true);

};

DEFINE_SPTR(OSNSMatrixProjectOnConstraints)

#endif
