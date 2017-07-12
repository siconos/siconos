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
/*! \file OSNSMatrix.hpp
  Specific storage for matrices used in OneStepNSProblem
*/

#ifndef OSNSM_H
#define OSNSM_H

#include "SiconosFwd.hpp"
#include "SimulationTypeDef.hpp"

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
 * - Sparse matrix (_storageType = 2): at the time of writting, only csc (compressed-sparse column).
 *   Could also be triplet (coo or coordinate) or csr (compressed-sparse row).
 */
class OSNSMatrix
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(OSNSMatrix);


  /** number of rows */
  unsigned int _dimRow;

  /** number of columns */
  unsigned int _dimColumn;

  /** Storage type used for the present matrix */
  int _storageType;

  /** Numerics structure to be filled  */
  SP::NumericsMatrix _numericsMatrix;

  /** Matrix used for default storage type (_storageType = 0) */
  SP::SiconosMatrix _M1;

  /** Matrix which corresponds to Numerics SparseBlockStructuredMatrix
      (_storageType = 1) */
  SP::BlockCSRMatrix _M2;

  /** For each Interaction in the graph, compute its absolute position
   *  \param indexSet the index set ot the concerned interactios.
   * \return the dimension of the problem (or size of the matrix),
   * computed as the sum of the nslaw of all the Interaction in indexSet
   */
  virtual unsigned updateSizeAndPositions(InteractionsGraph & indexSet);

  /** For each DynamicalSystem in the graph, compute its absolute position
   * \param indexSet the index set of the dynamical systems
   * \return the dimension of the problem (or size of the matrix),
   * computed as the sum of the nslaw of all the Interaction in indexSet
   */

  virtual unsigned updateSizeAndPositions(DynamicalSystemsGraph & DSG);

private:
  /** Private copy constructor => no copy nor pass by value */
  OSNSMatrix(const OSNSMatrix&);

  /** Private assignment -> forbidden
   * \return  OSNSMatrix&
   */
  OSNSMatrix& operator=(const OSNSMatrix&);

public:

  /** Default constructor -> empty matrix
   */
  OSNSMatrix();

  /** Constructor with _dimRow. of the matrix
   *   \param n size of the square matrix
   *   \param stor storage type (0:dense, 1:sparse interactionBlock)
   */
  OSNSMatrix(unsigned int n, int stor);

  /** Constructor with _dimRow and DimColumn of the matrix
   * \param n row sizes of the rectangle matrix
   * \param m column size of the rectangle matrix
   * \param stor storage type (0:dense, 1:sparse interactionBlock)
   */
  OSNSMatrix(unsigned int n, unsigned int m, int stor);

  /** Constructor from index set and map
   * \param indexSet InteractionsGraph* the index set of the active constraints
   * \param stor storage type
   */
  OSNSMatrix(InteractionsGraph& indexSet, int stor);

  /** Constructor with copy of a SiconosMatrix => _storageType = 0
   * \param MSource matrix to be copied
   */
  OSNSMatrix(const SiconosMatrix& MSource);

  /** destructor
   */
  virtual ~OSNSMatrix();

  /** get dimension of the square matrix
   * \return unsigned int
   */
  inline unsigned int size() const
  {
    return _dimRow;
  };

  /** get dimension of the square matrix
   * \return unsigned int
   */
  inline unsigned int sizeColumn() const
  {
    return _dimColumn;
  };

  /** get the type of storage for current matrix
   * \return unsigned int
   */
  inline int storagetype() const
  {
    return _storageType;
  };

  /** set which type of storage will be used for current matrix
   * \param i the type of storage
   */
  inline void setStorageType(int i)
  {
    _storageType = i;
  };

  /** get the numerics-readable structure
   * \return SP::NumericsMatrix
   */
  inline SP::NumericsMatrix numericsMatrix()
  {
    return _numericsMatrix;
  };

  /** get the matrix used for default storage
   * \return SP::NumericsMatrix
   */
  inline SP::SiconosMatrix defaultMatrix()
  {
    return _M1;
  };

  /** fill the current class using an index set
   * \param indexSet the index set of the active constraints
   * \param update if true update the size of the Matrix (default true)
   */
  virtual void fillW(InteractionsGraph&indexSet, bool update = true);

  /** fill the current class using an index set
   * \param indexSet the index set of the dynamicalSystems
   * \param update if true update the size of the Matrix (default true)
   */
  virtual void fillM(DynamicalSystemsGraph& DSG, bool update = true);

  /** fill the current class using an index set
   * \param indexSet the index set of the dynamicalSystems
   * \param update if true update the size of the Matrix (default true)
   */
  virtual void fillH(DynamicalSystemsGraph& DSG, InteractionsGraph& indexSet,  bool update = true);

  /** fill the numerics structure _numericsMatSparse using MBlockCSR */
  void convert();

  /** display the current matrix
   */
  void display() const;
};

#endif
