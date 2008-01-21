/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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

#include "SimpleMatrix.h"
#include "SiconosNumerics.h"
#include "SimulationTypeDef.h"

typedef std::vector<int> IndexInt;
typedef  boost::numeric::ublas::compressed_matrix<SiconosMatrix*> SparseMat2;

/** Definition of a sparse block matrix of SiconosMatrix*, used in OneStepNSProblem to save matrix M.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
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
 * This objects is built using an index set of UnitaryRelation*, that represents the "active" constraints in the OSNS problem and
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

  /** Number of blocks in a row/col*/
  unsigned int nc;

  /** Total number of non-null blocks in the matrix */
  unsigned int nbNonNullBlocks;

  /** Specific structure required when a (Numerics) solver block is used */
  SparseBlockStructuredMatrix* numericsMatSparse;

  /** Sparse-Block Boost Matrix. Each block is a SiconosMatrix**/
  SparseMat2 * MSparseBlock;

  /** vector of the addresses of the non-null blocks. */
  std::vector<double*> * blocksList;

  /** Vector used to save the sum of dim of diagonal blocks of M: diagSizes[i] = diagSizes[i-1] + ni, ni being the size of the diagonal block at row(block) i */
  IndexInt* diagSizes;

  /** List of non null blocks positions (in row) */
  IndexInt* rowPos;

  /** List of non null blocks positions (in col) */
  IndexInt* colPos;

  /** Private copy constructor => no copy nor pass by value */
  SparseBlockMatrix(const SparseBlockMatrix&);

  /** Private assignment -> forbidden */
  SparseBlockMatrix& operator=(const SparseBlockMatrix&);

public:

  /** Default constructor -> empty matrix
   */
  SparseBlockMatrix();

  /** Constructor with dim.
      \param nrow, int, number of blocks in a row
      \param ncol, int, number of blocks in a column (note that at the time we only authorize nrow = ncol).
  */
  SparseBlockMatrix(unsigned int, unsigned int);

  /** Constructor from index set and map
      \param UnitaryRelation*, the index set of the active constraints
      \param MapOfMapOfUnitaryMatrices, the list of matrices linked to a couple of UR*
  */
  SparseBlockMatrix(UnitaryRelationsSet*, MapOfMapOfUnitaryMatrices&);

  /** destructor
   */
  ~SparseBlockMatrix();

  /** get number of blocks in a row */
  inline const unsigned int size() const
  {
    return nc;
  };

  /** get total number of non-null blocks */
  inline const unsigned int getNbNonNullBlocks() const
  {
    return nbNonNullBlocks;
  };

  /** get the numerics-readable structure */
  inline SparseBlockStructuredMatrix* getNumericsMatSparse()
  {
    return numericsMatSparse;
  };

  /** get the ublas sparse mat*/
  inline SparseMat2 * getMSparse()
  {
    return MSparseBlock;
  };

  /** get the list of pointer to non null blocks */
  inline std::vector<double*>* getBlocksList()
  {
    return blocksList;
  };

  /** get the index of dimension of diagonale blocks */
  inline IndexInt * getDiagSizes()
  {
    return diagSizes;
  };

  /** get the index of blocks position (i=0 -> rows, i=1 -> columns)
      \param unsigned int, 0 for rows, 1 for columns
  */
  inline IndexInt * getPositionsIndex(bool i)
  {
    if (i) return rowPos;
    else return colPos;
  };

  /** fill the current class using an index set and a map of blocks
      \param UnitaryRelationsSet*, the index set of the active constraints
      \param MapOfMapOfUnitaryMatrices, the list of matrices linked to a couple of UR*
  */
  void fill(UnitaryRelationsSet*, MapOfMapOfUnitaryMatrices&);

  /** fill the numerics structure numericsMatSparse using MSparseBlock */
  void convert();

  /** display the current matrix
   */
  void display() const;
};

#endif
