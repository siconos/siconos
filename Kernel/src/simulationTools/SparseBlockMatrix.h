/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
/*! \file SparseBlockMatrix.h
Definition of a sparse block matrix of SiconosMatrix*
*/

#ifndef SBM_H
#define SBM_H

#include <boost/shared_ptr.hpp>
#include "SiconosNumerics.h"
#include "SimulationTypeDef.h"

typedef  boost::numeric::ublas::compressed_matrix<double*> CompressedRowMat;
TYPEDEF_SPTR(CompressedRowMat);


/** Definition of a sparse block matrix of SiconosMatrix*, used in OneStepNSProblem to save matrix M.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) 29/11/2007
 *
 * This class defines a specific sparse storage for blocks matrices, each block being a SiconosMatrix*.
 *
 * It handles:
 * - a SparseMat (boost-ublas) of SiconosMatrix*
 * - a vector<SiconosMatrix*> which handles the non-null blocks
 * - three vector<int> (IndexInt) to save non-null blocks position in row, columns and the list of the sizes of diagonal blocks.
 * - two int, the number of blocks in a row and the number of non null blocks.
 *
 * Each block of the current object represents the connection between two coupled Unitary Relations, \n
 * (for example for Lagrangian systems, a single \f$ H W^{-1} H^t \f$ block or for first order systems \f$ hCW^{-1}B \f$ ...) \n
 *
 * This objects is built using an index set of SP::UnitaryRelation, that represents the "active" constraints in the OSNS problem and
 * a map<UR* u1, <UR* u2, SiconosMatrix* block> >, block being the link between u1 and u2. Only UR present in the index set are picked out in the map.
 *
 *  A convert method is also implemented to create a SparseBlockStructuredMatrix which is Numerics-readable.
 *
 * As an example, consider the index set I={u1, u3, u5, u8} and the map where non null blocks are (ui,ui), (u1,u3), (u1,u8), (u3,u1), (u8,u1).\n
 * Each block being a pointer to a 3x3 matrix.\n
 * Then the resulting matrix has 4 X 4 blocks, with 8 non-null blocks and looks like:
 *
 * \f{eqnarray*}
   M=\left\lbrace\begin{array}{cccc}
    b11 & b13 & 0 & b18 \\
    b31 & b22 & 0 & 0 \\
    0   & 0   & b33&0 \\
    b81 & 0   & 0 & b44
    \end{array}\right.
    \f}
 *
 * with nc = 4, nbNonNullBlocks = 8, RowPos = [0 0 0 1 1 2 3 3], RowCol = [0 1 3 0 1 2 0 3]\n
 * and diagSizes = [3 6 9 12].
 *
 * We use stl::vector (which may seems redundent with the double* of the numerics SparseBlockStructuredMatrix) because memory can be
 * reserved during construction or initialized and then vectors are resized when the object is filled in. This avoid some call to
 * malloc/free at each iteration.
 *
 * Note FP: MSparseBlock boost-ublas matrix is not used at the time, because not readable by Numerics. But I leave it for future improvements
 * (Numerics with boost or implement a way to get the address of the embedded double** in boost sparse matrix).
 *
 */
class SparseBlockMatrix
{
private:
  /** Number of blocks in  row*/
  unsigned int nr;

  /** Number of blocks in  col*/
  unsigned int nc;

  /** Specific structure required when a (Numerics) solver block is used */
  boost::shared_ptr<SparseBlockStructuredMatrix> numericsMatSparse;

  /** Sparse-Block Boost Matrix. Each block is a SiconosMatrix**/
  boost::shared_ptr<CompressedRowMat> MSparseBlock;

  /** Vector used to save the sum of dim of diagonal blocks of M: diagSizes[i] = diagSizes[i-1] + ni, ni being the size of the diagonal block at row(block) i */
  boost::shared_ptr<IndexInt> diagSizes;

  /** List of non null blocks positions (in row) */
  boost::shared_ptr<IndexInt> rowPos;

  /** List of non null blocks positions (in col) */
  boost::shared_ptr<IndexInt> colPos;

  /** Private copy constructor => no copy nor pass by value */
  SparseBlockMatrix(const SparseBlockMatrix&);

  /** Private assignment -> forbidden */
  SparseBlockMatrix& operator=(const SparseBlockMatrix&);

public:

  /** Default constructor -> empty matrix
   */
  SparseBlockMatrix();

  /** Constructor with dimension (number of blocks)
      \param n number of blocks in a row/column (only square matrices allowed)
  */
  SparseBlockMatrix(unsigned int);

  /** Constructor from index set and map
      \param SP::UnitaryRelation, the index set of the active constraints
      \param MapOfMapOfUnitaryMatrices, the list of matrices linked to a couple of UR*
  */
  SparseBlockMatrix(SP::UnitaryRelationsSet, MapOfMapOfUnitaryMatrices&);

  /** Constructor from DynamicalSystemsSet and map
      \param DynamicalSystemsSet*, the index set of the active constraints
      \param MapOfDSMatrices, the list of matrices linked to a couple of UR*
  */
  SparseBlockMatrix(SP::DynamicalSystemsSet, MapOfDSMatrices&);

  /** Constructor from DynamicalSystemsSet and map
    \param SP::UnitaryRelation, the index set of the active constraints
     \param DynamicalSystemsSet*, the index set of the active constraints
     \param MapOfDSMatrices, the list of matrices linked to a couple of UR*
  */
  SparseBlockMatrix(SP::UnitaryRelationsSet, SP::DynamicalSystemsSet, MapOfUnitaryMapOfDSMatrices&);

  /** destructor
   */
  ~SparseBlockMatrix();

  /** get size (in block-components) */
  inline const unsigned int getNumberOfBlocksInARow() const
  {
    return nr;
  };

  /** get total number of non-null blocks */
  inline const unsigned int getNbNonNullBlocks() const
  {
    return MSparseBlock->nnz();
  };

  /** get the numerics-readable structure */
  inline SP::SparseBlockStructuredMatrix getNumericsMatSparse()
  {
    return numericsMatSparse;
  };

  /** get the ublas sparse mat*/
  inline SP::CompressedRowMat getMSparse()
  {
    return MSparseBlock;
  };

  /** get the dimension of the square-diagonal block number num
      \param num block position
  */
  unsigned int getSizeOfDiagonalBlock(int i) const
  {
    if (i == 0) return diagSizes->at(0);
    else return (diagSizes->at(i) - diagSizes->at(i - 1));
  };

  /** get the index of blocks position (i=0 -> rows, i=1 -> columns)
      \param unsigned int, 0 for rows, 1 for columns
  */
  inline SP::IndexInt getPositionsIndex(bool i)
  {
    if (i) return rowPos;
    else return colPos;
  };

  /** fill the current class using an index set and a map of blocks
      \param UnitaryRelationsSet*, the index set of the active constraints
      \param MapOfMapOfUnitaryMatrices, the list of matrices linked to a couple of UR*
  */
  void fill(SP::UnitaryRelationsSet, MapOfMapOfUnitaryMatrices&);

  /** fill the current class using an index set and a map of DSblocks
       \param DynamicalSystemsSet*, the  set of DynamicalSystem
       \param MapOfDSMatrices, the list of matrices linked to a DynamicalSystem
   */
  void fill(SP::DynamicalSystemsSet, MapOfDSMatrices&);

  /** fill the current class using an index set and a map of DSblocks
       \param DynamicalSystemsSet*, the  set of DynamicalSystem
       \param UnitaryRelationsSet*, the index set of the active constraints
       \param MapOfUnitaryMapOfDSMatrices, the list of matrices linked to a DynamicalSystem
   */
  void fill(SP::UnitaryRelationsSet, SP::DynamicalSystemsSet, MapOfUnitaryMapOfDSMatrices&);

  /** fill the numerics structure numericsMatSparse using MSparseBlock */
  void convert();

  /** display the current matrix
   */
  void display() const;
};

#endif
