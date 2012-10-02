/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
/*! \file OSNSMatrix.hpp
  Specific storage for matrices used in OneStepNSProblem
*/

#ifndef OSNSM_H
#define OSNSM_H

#include <boost/shared_ptr.hpp>
#include "SimpleMatrix.hpp"
#include "SiconosNumerics.h"
#include "SimulationTypeDef.hpp"
#include "Topology.hpp"
#include "BlockCSRMatrix.hpp"

TYPEDEF_SPTR(NumericsMatrix)

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
 *  - full matrix in a SiconosMatrix (storageType = 0). In this case,
 *  for each call to fill(), the SiconosMatrix M is resized
 *  according  to the sizes of the Interaction present in indexSet and then
 *  all the required interactionBlocks mij are COPIED into M.
 *
 *  - Sparse Block Storage (storageType = 1): corresponds to
 *  SparseBlockStructuredMatrix structure of Numerics. Only non-null
 *  interactionBlocks are saved in the matrix M and there is no copy of
 *  sub-interactionBlocks, only links thanks to pointers.
 *
 */
enum SICONOS_STORAGE_TYPE
{
  SICONOS_DENSE = 0,
  SICONOS_SPARSE = 1
};

class OSNSMatrix
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(OSNSMatrix);


  /** number of rows */
  unsigned int dimRow;

  /** number of columns */
  unsigned int dimColumn;

  /** Storage type used for the present matrix */
  int storageType;

  /** map that links each Interaction with an int that gives the
   * position (in number of scalar elements, not interactionBlocks) \n of
   * the corresponding interactionBlock matrix in the full matrix (M in
   * LCP case) - Warning: it depends on the considered index set \n
   * (ie on which constraints are "active")
   */
  //  SP::Interaction_int interactionBlocksPositions;

  /** map that links each DynamicalSystem with an int that gives the
   * position (in number of scalar elements, not DSBlocks) of the
   * corresponding DSBlock matrix in the full matrix (M in
   * PrimalFrictionalCase case) - Warning: it depends on the
   * considered index set (ie on which constraints are "active")
   */
  SP::DS_int DSBlocksPositions;

  /** Numerics structure to be filled  */
  SP::NumericsMatrix numericsMat;

  /** Matrix used for default storage type (storageType = 0) */
  SP::SiconosMatrix M1;
  SP::SiconosMatrix Mt;

  /** Matrix which corresponds to Numerics SparseBlockStructuredMatrix
      (storageType = 1) */
  SP::BlockCSRMatrix M2;

  /** Private copy constructor => no copy nor pass by value */
  OSNSMatrix(const OSNSMatrix&);

  /** Private assignment -> forbidden */
  OSNSMatrix& operator=(const OSNSMatrix&);

  /** To update dim and interactionBlocksPositions for a new set of
      Interaction
      \param InteractionsGraph* the index set of
      the active constraints
  */
  virtual void updateSizeAndPositions(unsigned int&, SP::InteractionsGraph);

  /** To update dim and DSBlocksPositions for a new set of DynamicalSystem
      \param DynamicalSystemsSet* the DyncamicalSystemsSet
  */
  void updateSizeAndPositions(unsigned int&, SP::DynamicalSystemsSet);

  /** To update dim, DSBlocksPositions and interactionBlocksPositions for
      a new set of DynamicalSystem and a new set of Interaction
      \param DynamicalSystemsSet* the DynamicalSystemsSet
      \param InteractionsGraph* the index set of the active constraints
  */
  void updateSizeAndPositions(unsigned int&, SP::DynamicalSystemsSet, SP::InteractionsGraph);

public:

  /** Default constructor -> empty matrix
   */
  OSNSMatrix();

  /** Constructor with dimRow. of the matrix
      \param n size of the square matrix
      \param stor storage type (0:dense, 1:sparse interactionBlock)
  */
  OSNSMatrix(unsigned int, int);

  /** Constructor with dimRow and DimColumn of the matrix
      \param n and m sizes of the rectangle matrix
      \param stor storage type (0:dense, 1:sparse interactionBlock)
  */
  OSNSMatrix(unsigned int, unsigned int, int);

  /** Constructor from index set and map
      \param InteractionsGraph* the index set of the active constraints
      \param storage type
  */
  OSNSMatrix(SP::InteractionsGraph, int);

  /** Constructor from DynamicalSystemsSet and map
      \param InteractionsGraph* the index set of the active constraints
      \param MapOfDSMatrices the list of matrices linked to a couple of SP::Interaction
      \param storage type
  */
  OSNSMatrix(SP::DynamicalSystemsSet, MapOfDSMatrices&, int);

  /** Constructor from DynamicalSystemsSet and indexSet and map
      \param InteractionsGraph* the index set of the active constraints
      \param MapOfMapOfInteractionMatrices the list of matrices linked to a couple of SP::Interaction
      \param storage type
  */
  OSNSMatrix(SP::DynamicalSystemsSet, SP::InteractionsGraph, MapOfDSMapOfInteractionMatrices&, int);

  /** Constructor from DynamicalSystemsSet and indexSet and map
      \param InteractionsGraph* the index set of the active constraints
      \param MapOfMapOfInteractionMatrices the list of matrices linked to a couple of SP::Interaction
      \param storage type
  */
  OSNSMatrix(SP::InteractionsGraph, SP::DynamicalSystemsSet , MapOfInteractionMapOfDSMatrices&, int);

  /** Constructor from DynamicalSystemsSet and indexSet and maps of Blocks
      \param InteractionsGraph* the index set of the active constraints
      \param storage type
  */
  OSNSMatrix(SP::InteractionsGraph, SP::DynamicalSystemsSet,  MapOfDSMatrices&, MapOfDSMapOfInteractionMatrices&,  MapOfInteractionMapOfDSMatrices&, int);

  /** Constructor with copy of a SiconosMatrix => storageType = 0
      \param MSource matrix to be copied
  */
  OSNSMatrix(const SiconosMatrix&);

  /** destructor
   */
  virtual ~OSNSMatrix();

  /** get dimension of the square matrix */
  inline unsigned int size() const
  {
    return dimRow;
  };

  /** get dimension of the square matrix */
  inline unsigned int sizeColumn() const
  {
    return dimColumn;
  };

  /** get the type of storage for current matrix  */
  inline int getStorageType() const
  {
    return storageType;
  };

  /** set which type of storage will be used for current matrix */
  inline void setStorageType(int i)
  {
    storageType = i;
  };

  /** get the position (real, not interactionBlock) of the first element of the interactionBlock which corresponds to the Interaction
      \param Interaction Interaction from which position is required
  */
  unsigned int getPositionOfDSBlock(SP::DynamicalSystem) const;

  /** get the position (real, not DSBlock) of the first element of the DSBlock which corresponds to DS
      \param DS DynamicalSystem  from which position is required
  */
  virtual unsigned int getPositionOfInteractionBlock(SP::Interaction) const;

  /** get the numerics-readable structure */
  inline SP::NumericsMatrix getNumericsMatrix()
  {
    return numericsMat;
  };

  /** get the matrix used for default storage */
  inline SP::SiconosMatrix defaultMatrix()
  {
    return M1;
  };

  /** fill the current class using an index set and a map of interactionBlocks
      \param InteractionsGraph*, the index set of the active constraints
  */
  virtual void fill(SP::InteractionsGraph, bool updateSize = true);
  /** fill diagonal of thecurrent class using an index set and a map of interactionBlocks
      \param InteractionsGraph*, the index set of the active constraints
  */
  void fillDiagonal(SP::InteractionsGraph, bool updateSize = true);

  /** fill the current class using an DynamicalSystemsSet and a map of DSBlocks
      \param DynamicalSystemsSet*, the Dynamical set
      \param MapOfDSMatrices, the list of matrices linked to a DynamicalSystems
  */
  void fill(SP::DynamicalSystemsSet, MapOfDSMatrices&, bool updateSize = true);
  /** fill the current class using an index set , a DynamicalSystemsSet and a map of interactionBlocks
      \param InteractionsGraph*, the index set of the active constraints
      \param DynamicalSystemsSet*, the Dynamical set
      \param MapOfMapOfInteractionMatrices, the list of matrices linked to a couple of SP::Interaction
  */
  void fill(SP::DynamicalSystemsSet, SP::InteractionsGraph, MapOfDSMapOfInteractionMatrices&, bool updateSize = true);

  /** fill the current class using an index set and a map of interactionBlocks
      \param InteractionsGraph*, the index set of the active constraints
      \param DynamicalSystemsSet*, the Dynamical set
      \param MapOfMapOfInteractionMatrices, the list of matrices linked to a couple of SP::Interaction
  */
  void fill(SP::InteractionsGraph, SP::DynamicalSystemsSet, MapOfInteractionMapOfDSMatrices&, bool updateSize = true);

  /** fill the current class using an index set and  maps of Blocks
      \param InteractionsGraph*, the index set of the active constraints
      \param DynamicalSystemsSet*, the Dynamical set
  */
  void fill(SP::InteractionsGraph, SP::DynamicalSystemsSet,  MapOfDSMatrices&, MapOfDSMapOfInteractionMatrices&,  MapOfInteractionMapOfDSMatrices&, bool updateSize = true);

  /** fill the numerics structure numericsMatSparse using MBlockCSR */
  void convert();

  /** display the current matrix
   */
  void display() const;
};

DEFINE_SPTR(OSNSMatrix)

#endif
