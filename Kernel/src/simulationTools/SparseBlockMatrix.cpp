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
#include "SparseBlockMatrix.h"

using namespace std;

// Default constructor: empty matrix
SparseBlockMatrix::SparseBlockMatrix():
  nr(0)
{
  MSparseBlock.reset(new CompressedRowMat());
  numericsMatSparse.reset(new SparseBlockStructuredMatrix());
  diagSizes.reset(new IndexInt());
  rowPos.reset(new IndexInt());
  colPos.reset(new IndexInt());
}

// Constructor with dimensions
SparseBlockMatrix::SparseBlockMatrix(unsigned int nRow):
  nr(nRow)
{
  // Only square-blocks matrices for the moment (ie nRow = nr = nrol)

  // Allocate memory and fill in the matrix
  // rowPos, rowCol ... are initialized with nr to reserve at first step the maximum possible (according to given nr) space in memory.
  // Thus a future resize will not require memory allocation or copy.

  MSparseBlock.reset(new CompressedRowMat(nr, nr));
  numericsMatSparse.reset(new SparseBlockStructuredMatrix);
  diagSizes.reset(new IndexInt());
  rowPos.reset(new IndexInt());
  colPos.reset(new IndexInt());
  diagSizes->reserve(nr);
  rowPos->reserve(nr);
  colPos->reserve(nr);
}

// Basic constructor
SparseBlockMatrix::SparseBlockMatrix(SP::UnitaryRelationsSet I, MapOfMapOfUnitaryMatrices& blocks):
  nr(0)
{
  // Allocate memory and fill in the matrix
  nr = I->size();
  MSparseBlock.reset(new CompressedRowMat(nr, nr));
  numericsMatSparse.reset(new SparseBlockStructuredMatrix());
  diagSizes.reset(new IndexInt());
  diagSizes->reserve(nr);
  rowPos.reset(new IndexInt());
  rowPos->reserve(nr);
  colPos.reset(new IndexInt());
  colPos->reserve(nr);
  fill(I, blocks);
}
SparseBlockMatrix::SparseBlockMatrix(SP::DynamicalSystemsSet DSSet, MapOfDSMatrices& DSblocks):
  nr(0)
{
  // Allocate memory and fill in the matrix
  nr = DSSet->size();
  MSparseBlock.reset(new CompressedRowMat(nr, nr));
  numericsMatSparse.reset(new SparseBlockStructuredMatrix());
  diagSizes.reset(new IndexInt());
  diagSizes->reserve(nr);
  rowPos.reset(new IndexInt());
  rowPos->reserve(nr);
  colPos.reset(new IndexInt());
  colPos->reserve(nr);
  fill(DSSet, DSblocks);
}
SparseBlockMatrix::SparseBlockMatrix(SP::UnitaryRelationsSet I, SP::DynamicalSystemsSet DSSet, MapOfUnitaryMapOfDSMatrices& unitaryDSblocks):
  nr(0)
{
  // Allocate memory and fill in the matrix
  nr = I->size();
  nc = DSSet->size();
  MSparseBlock.reset(new CompressedRowMat(nr, nc));
  numericsMatSparse.reset(new SparseBlockStructuredMatrix());
  diagSizes.reset(new IndexInt());
  diagSizes->reserve(nr);
  rowPos.reset(new IndexInt());
  rowPos->reserve(nr);
  colPos.reset(new IndexInt());
  colPos->reserve(nr);
  fill(I, DSSet, unitaryDSblocks);
}

// Destructor -> see smart pointers
SparseBlockMatrix::~SparseBlockMatrix()
{

}

// Fill the SparseMat
void SparseBlockMatrix::fill(SP::UnitaryRelationsSet indexSet, MapOfMapOfUnitaryMatrices& blocks)
{
  // ======>  Aim: find UR1 and UR2 both in indexSets[level] and which have common DynamicalSystems.
  // Then get the corresponding matrix from map blocks.

  assert(indexSet.get() && "NULL pointer");

  // Number of blocks in a row = number of active constraints.
  nr = indexSet->size();

  // (re)allocate memory for ublas matrix
  MSparseBlock->resize(nr, nr, false);
  diagSizes->resize(nr);

  size_t row = 0; // (block) row position
  size_t col = 0; // (block) col. position

  // reset rowPos, colPos
  rowPos->resize(0);
  colPos->resize(0);
  int sizeV = 0;

  // === Loop through "active" Unitary Relations (ie present in indexSets[level]) ===
  for (UnitaryRelationsIterator itRow = indexSet->begin(); itRow != indexSet->end(); ++itRow)
  {
    // Size of the current diagonal block is added into diagSizes
    sizeV  += (*itRow)->getNonSmoothLawSize();
    (*diagSizes)[row] = sizeV;
    // Look for UR connected (ie with common DS) to the current one
    // The connection is checked thanks to blocks map
    for (UnitaryRelationsIterator itCol = indexSet->begin(); itCol != indexSet->end(); ++itCol)
    {
      //for(UnitaryMatrixColumnIterator itCol = blocks[*itRow].begin(); itCol!=blocks[*itRow].end(); ++itCol)
      //{

      if ((blocks[*itRow]).find(*itCol) != (blocks[*itRow]).end())
      {
        if (blocks[*itRow][*itCol] != NULL)
        {
          // We keep connected UR only if they are also in indexSets[level] and if the block is not NULL
          //   if( (indexSet->find((*itCol).first))!= indexSet->end() && blocks[*itRow][(*itCol).first] != NULL)
          {
            // UR1 = *itRow
            // UR2 = *itCol
            // corresponding matrix = blocks[*itRow][(*itCol).first]

            // Insert block into the list
            (*MSparseBlock)(row, col) = blocks[*itRow][*itCol]->getArray();

            // Insert block positions into rowPos and colPos
            (rowPos)->push_back(row);
            (colPos)->push_back(col);

          }
        }
      }
      col++;// increment col position
    }
    col = 0;
    row++;// increment row position
  }
}

// Fill the SparseMat
void SparseBlockMatrix::fill(SP::DynamicalSystemsSet DSSet, MapOfDSMatrices& DSblocks)
{
  RuntimeException::selfThrow(" SparseBlockMatrix::fill(DynamicalSystemsSet* DSSet, MapOfDSMatrices& DSblocks), Not Yet Implemented");
}
// Fill the SparseMat
void SparseBlockMatrix::fill(SP::UnitaryRelationsSet indexSet, SP::DynamicalSystemsSet DSSet, MapOfUnitaryMapOfDSMatrices& unitaryDSblocks)
{
  RuntimeException::selfThrow(" SparseBlockMatrix::fill(DynamicalSystemsSet* DSSet, MapOfDSMatrices& DSblocks), Not Yet Implemented");
}


// convert MSparseBlock to numerics structure
void SparseBlockMatrix::convert()
{
  numericsMatSparse->size = nr;
  numericsMatSparse->nbblocks = (*MSparseBlock).nnz();
  // Next copies: pointer links!!
  numericsMatSparse->blocksize =  &((*diagSizes)[0]);
  numericsMatSparse->RowIndex = &((*rowPos)[0]);
  numericsMatSparse->ColumnIndex = &((*colPos)[0]);

  // boost
  numericsMatSparse->filled1 = (*MSparseBlock).filled1();
  numericsMatSparse->filled2 = (*MSparseBlock).filled2();
  numericsMatSparse->index1_data = &((*MSparseBlock).index1_data()[0]);
  numericsMatSparse->index2_data = &((*MSparseBlock).index2_data()[0]);
  numericsMatSparse->block =  &((*MSparseBlock).value_data()[0]);

  //   // Loop through the non-null blocks
  //   for (SpMatIt1 i1 = MSparseBlock->begin1(); i1 != MSparseBlock->end1(); ++i1)
  //     {
  //       for (SpMatIt2 i2 = i1.begin(); i2 != i1.end(); ++i2)
  //  {
  //    block[i] = *i2;
  //  }
  //     }
}

// Display data
void SparseBlockMatrix::display() const
{
  cout << "----- Sparse Block Matrix with " << nr << " blocks in a row/col and " << MSparseBlock->nnz() << " non-null blocks" << endl;
  cout << "Non null blocks positions are (row index and columns index):" << endl;
  print(rowPos->begin(), rowPos->end());
  print(colPos->begin(), colPos->end());
  cout << "Sum of sizes of the diagonal blocks:" << endl;
  print(diagSizes->begin(), diagSizes->end());
}


