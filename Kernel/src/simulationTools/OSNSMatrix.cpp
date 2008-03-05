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
#include "OSNSMatrix.h"
#include "SparseBlockMatrix.h"
#include "Tools.h"

using namespace std;

void OSNSMatrix::updateSizeAndPositions(UnitaryRelationsSet* indexSet)
{
  // === Description ===
  // For a block (diagonal or extra diagonal) corresponding to a Unitary Relation, we need to know the position of its first element
  // in the full-matrix M. This position depends on the previous blocks sizes.
  //
  // positions are saved in a map<UnitaryRelation*, unsigned int>, named blocksPositions.
  //

  // Computes real size of the current matrix = sum of the dim. of all UR in indexSet
  dim = 0;
  for (UnitaryRelationsIterator it = indexSet->begin(); it != indexSet->end(); ++it)
  {
    (*blocksPositions)[*it] = dim;
    dim += (*it)->getNonSmoothLawSize();
  }
}

// Default constructor: empty matrix, default storage
// No allocation for M1 or M2
OSNSMatrix::OSNSMatrix():
  dim(0), storageType(0), blocksPositions(NULL), numericsMat(NULL), M1(NULL), M2(NULL)
{
  blocksPositions = new UR_int();
  numericsMat = new NumericsMatrix;
}

// Constructor with dimensions (one input: square matrix only)
OSNSMatrix::OSNSMatrix(unsigned int n, int stor):
  dim(n), storageType(stor), blocksPositions(NULL), numericsMat(NULL), M1(NULL), M2(NULL)
{
  // Note:
  // for storageType = 0 (dense) n represents the real dimension of the matrix
  // and for sparse storage (storageType == 1) the number of blocks in a row or column.

  if (storageType == 0)
  {
    // A zero matrix M of size nXn is built.
    // blocksPositions remains empty (=NULL) since we have no information concerning the UR.
    M1 = new SimpleMatrix(n, n);
  }
  else // if(storageType == 1)
    M2 = new SparseBlockMatrix(n);
  blocksPositions = new UR_int();
  numericsMat = new NumericsMatrix;
}

// Basic constructor
OSNSMatrix::OSNSMatrix(UnitaryRelationsSet* indexSet, MapOfMapOfUnitaryMatrices& blocks, int stor):
  dim(0), storageType(stor), blocksPositions(NULL), numericsMat(NULL), M1(NULL), M2(NULL)
{
  blocksPositions = new UR_int();
  numericsMat = new NumericsMatrix;

  fill(indexSet, blocks);
}

// Copy of a SiconosMatrix (used when OSNS xml constructor is called with M input in XML file)
OSNSMatrix::OSNSMatrix(const SiconosMatrix& MSource):
  dim(MSource.size(0)), storageType(0), blocksPositions(NULL), numericsMat(NULL), M1(NULL), M2(NULL)
{
  blocksPositions = new UR_int();
  numericsMat = new NumericsMatrix;
  M1 = new SimpleMatrix(MSource);
  // Warning: blocksPositions remains empty since we have no information concerning indexSet and blocks in MSource
}


// Destructor
OSNSMatrix::~OSNSMatrix()
{
  blocksPositions->clear();
  delete blocksPositions;
  blocksPositions = NULL;
  delete numericsMat;
  numericsMat = NULL;
  if (M1 != NULL) delete M1;
  M1 = NULL;
  if (M2 != NULL) delete M2;
  M2 = NULL;
}

unsigned int OSNSMatrix::getPositionOfBlock(UnitaryRelation* UR) const
{
  UR_int::const_iterator it = blocksPositions->find(UR);
  if (it == blocksPositions->end())
    RuntimeException::selfThrow("OSNSMatrix::getPositionOfBlock(UnitaryRelation ur), ur does not belong to the index set used to built the OSNS matrix.");
  return it->second;
}

// Fill the SparseMat
void OSNSMatrix::fill(UnitaryRelationsSet* indexSet, MapOfMapOfUnitaryMatrices& blocks)
{
  if (indexSet == NULL)
    RuntimeException::selfThrow("OSNSMatrix::fill(IndexInt* i, ...), i is a null pointer");

  // Computes dim and blocksPositions according to indexSet
  updateSizeAndPositions(indexSet);

  if (storageType == 0)
  {
    // === Memory allocation, if required ===
    // Mem. is allocate only if M==NULL or if its size has changed.
    if (M1 == NULL)
      M1 = new SimpleMatrix(dim, dim);
    else
    {
      if (M1->size(0) != dim || M1->size(1) != dim)
        M1->resize(dim, dim);
      M1->zero();
    }

    // ======>  Aim: find UR1 and UR2 both in indexSet and which have common DynamicalSystems.
    // Then get the corresponding matrix from map blocks, and copy it into M

    unsigned int pos = 0, col = 0; // index position used for block copy into M, see below.
    // === Loop through "active" Unitary Relations (ie present in indexSets[level]) ===
    for (UnitaryRelationsIterator itRow = indexSet->begin(); itRow != indexSet->end(); ++itRow)
    {
      // Look for UR connected (ie with common DS) to the current one
      // The connection is checked thanks to blocks map
      for (UnitaryMatrixColumnIterator itCol = blocks[*itRow].begin(); itCol != blocks[*itRow].end(); ++itCol)
      {
        // We keep connected UR only if they are also in indexSets[level]
        if ((indexSet->find((*itCol).first)) != indexSet->end())
        {
          // UR1 = *itRow
          // UR2 = *itCol
          // corresponding matrix = blocks[*itRow][(*itCol).first]

          // Case 1: basic storage
          pos = (*blocksPositions)[*itRow];
          col = (*blocksPositions)[(*itCol).first];
          static_cast<SimpleMatrix*>(M1)->setBlock(pos, col, *(blocks[*itRow][(*itCol).first]));
        }
      }
    }
  }
  else // if storageType == 1
  {
    if (M2 == NULL)
      M2 = new SparseBlockMatrix(indexSet, blocks);
    else
      M2->fill(indexSet, blocks);
  }
  convert();
}

// convert current matrix to NumericsMatrix structure
void OSNSMatrix::convert()
{
  numericsMat->storageType = storageType;
  numericsMat->size0 = dim;
  numericsMat->size1 = dim;
  if (storageType == 0)
  {
    numericsMat->matrix0 = M1->getArray(); // Pointer link
    // numericsMat->matrix1 = NULL;
    // matrix1 is not set to NULL: we keep previous allocation. May be usefull if we switch between different storages during simu
  }
  else
  {
    M2->convert();
    numericsMat->matrix1 = M2->getNumericsMatSparse();
  }
}

// Display data
void OSNSMatrix::display() const
{
  if (storageType == 0)
  {
    cout << "----- OSNS Matrix using default storage type for Numerics structure (SiconosMatrix -> double*)" << endl;
    if (M1 == NULL)
      cout << " matrix = NULL pointer" << endl;
    else M1->display();
  }
  else
  {
    cout << "----- OSNS Matrix using Sparse Block storage type for Numerics (SparseBlockStructuredMatrix)" << endl;
    if (M2 == NULL)
      cout << " matrix = NULL pointer" << endl;
    else M2->display();
  }
}
