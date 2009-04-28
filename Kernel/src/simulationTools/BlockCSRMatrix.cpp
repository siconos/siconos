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
#include "BlockCSRMatrix.h"

using namespace std;

// Default constructor: empty matrix
BlockCSRMatrix::BlockCSRMatrix():
  nr(0)
{
  MBlockCSR.reset(new CompressedRowMat());
  numericsMatSparse.reset(new SparseBlockStructuredMatrix());
  diagSizes.reset(new IndexInt());
  rowPos.reset(new IndexInt());
  colPos.reset(new IndexInt());
}

// Constructor with dimensions
BlockCSRMatrix::BlockCSRMatrix(unsigned int nRow):
  nr(nRow)
{
  // Only square-blocks matrices for the moment (ie nRow = nr = nrol)

  // Allocate memory and fill in the matrix rowPos, rowCol ... are
  // initialized with nr to reserve at first step the maximum possible
  // (according to given nr) space in memory.  Thus a future resize
  // will not require memory allocation or copy.

  MBlockCSR.reset(new CompressedRowMat(nr, nr));
  numericsMatSparse.reset(new SparseBlockStructuredMatrix);
  diagSizes.reset(new IndexInt());
  rowPos.reset(new IndexInt());
  colPos.reset(new IndexInt());
  diagSizes->reserve(nr);
  rowPos->reserve(nr);
  colPos->reserve(nr);
}

// Basic constructor
BlockCSRMatrix::BlockCSRMatrix(SP::UnitaryRelationsGraph indexSet,
                               MapOfMapOfUnitaryMatrices& blocks):
  nr(0)
{
  // Allocate memory and fill in the matrix
  nr = indexSet->size();
  MBlockCSR.reset(new CompressedRowMat(nr, nr));
  numericsMatSparse.reset(new SparseBlockStructuredMatrix());
  diagSizes.reset(new IndexInt());
  diagSizes->reserve(nr);
  rowPos.reset(new IndexInt());
  rowPos->reserve(nr);
  colPos.reset(new IndexInt());
  colPos->reserve(nr);
  fill(indexSet, blocks);
}
BlockCSRMatrix::BlockCSRMatrix(SP::DynamicalSystemsSet DSSet,
                               MapOfDSMatrices& DSblocks):
  nr(0)
{
  // Allocate memory and fill in the matrix
  nr = DSSet->size();
  MBlockCSR.reset(new CompressedRowMat(nr, nr));
  numericsMatSparse.reset(new SparseBlockStructuredMatrix());
  diagSizes.reset(new IndexInt());
  diagSizes->reserve(nr);
  rowPos.reset(new IndexInt());
  rowPos->reserve(nr);
  colPos.reset(new IndexInt());
  colPos->reserve(nr);
  fill(DSSet, DSblocks);
}
BlockCSRMatrix::BlockCSRMatrix(SP::UnitaryRelationsGraph indexSet,
                               SP::DynamicalSystemsSet DSSet,
                               MapOfUnitaryMapOfDSMatrices& unitaryDSblocks):
  nr(0)
{
  // Allocate memory and fill in the matrix
  nr = indexSet->size();
  nc = DSSet->size();
  MBlockCSR.reset(new CompressedRowMat(nr, nc));
  numericsMatSparse.reset(new SparseBlockStructuredMatrix());
  diagSizes.reset(new IndexInt());
  diagSizes->reserve(nr);
  rowPos.reset(new IndexInt());
  rowPos->reserve(nr);
  colPos.reset(new IndexInt());
  colPos->reserve(nr);
  fill(indexSet, DSSet, unitaryDSblocks);
}

// Destructor -> see smart pointers
BlockCSRMatrix::~BlockCSRMatrix()
{

}

// Fill the SparseMat
void BlockCSRMatrix::fill(SP::UnitaryRelationsGraph indexSet,
                          MapOfMapOfUnitaryMatrices& blocks)
{
  // ======> Aim: find UR1 and UR2 both in indexSets[level] and which
  // have common DynamicalSystems.  Then get the corresponding matrix
  // from map blocks.

  assert(indexSet);

  // Number of blocks in a row = number of active constraints.
  nr = indexSet->size();

  // (re)allocate memory for ublas matrix
  MBlockCSR->resize(nr, nr, false);
  diagSizes->resize(nr);

  int sizeV = 0;

  // === Loop through "active" Unitary Relations (ie present in
  // indexSets[level]) ===
  UnitaryRelationsGraph::VIterator ui1, ui1end;
  for (boost::tie(ui1, ui1end) = indexSet->vertices();
       ui1 != ui1end; ++ui1)
  {
    SP::UnitaryRelation ur1 = indexSet->bundle(*ui1);
    sizeV  += ur1->getNonSmoothLawSize();
    (*diagSizes)[indexSet->index(*ui1)] = sizeV;

    assert((*diagSizes)[indexSet->index(*ui1)] > 0);
    assert(blocks[ur1][ur1]);
    assert(blocks[ur1][ur1]->size(0) > 0 && blocks[ur1][ur1]->size(1) > 0);
    assert(blocks[ur1][ur1]->size(0) == ur1->getNonSmoothLawSize());
    assert(blocks[ur1][ur1]->size(1) == ur1->getNonSmoothLawSize());
    assert(indexSet->index(*ui1) < nr);

    (*MBlockCSR)(indexSet->index(*ui1), indexSet->index(*ui1)) =
      blocks[ur1][ur1]->getArray();

    UnitaryRelationsGraph::AVIterator ui2, ui2end;
    for (boost::tie(ui2, ui2end) =
           indexSet->adjacent_vertices(*ui1);
         ui2 != ui2end; ++ui2)
    {
      SP::UnitaryRelation ur2 = indexSet->bundle(*ui2);

      assert(indexSet->is_vertex(ur2));
      assert(ur2 != ur1);
      assert(blocks[ur1][ur2]);
      assert(blocks[ur1][ur2]->size(0) > 0 && blocks[ur1][ur2]->size(1) > 0);
      assert(blocks[ur1][ur2]->size(0) == ur1->getNonSmoothLawSize());
      assert(blocks[ur1][ur2]->size(1) == ur2->getNonSmoothLawSize());

      assert(indexSet->index(*ui1) < nr);
      assert(indexSet->index(*ui2) < nr);

      assert(*ui2 == indexSet->descriptor(ur2));
      assert(indexSet->index(*ui2) == indexSet->index(indexSet->descriptor(ur2)));

      (*MBlockCSR)(indexSet->index(*ui1), indexSet->index(*ui2)) =
        blocks[ur1][ur2]->getArray();

    }
  }
}

// Fill the SparseMat
void BlockCSRMatrix::fill(SP::DynamicalSystemsSet DSSet,
                          MapOfDSMatrices& DSblocks)
{
  RuntimeException::selfThrow
  (" BlockCSRMatrix::fill(DynamicalSystemsSet* DSSet, MapOfDSMatrices& DSblocks), Not Yet Implemented");
}
// Fill the SparseMat
void BlockCSRMatrix::fill(SP::UnitaryRelationsGraph indexSet,
                          SP::DynamicalSystemsSet DSSet,
                          MapOfUnitaryMapOfDSMatrices& unitaryDSblocks)
{
  RuntimeException::selfThrow
  (" BlockCSRMatrix::fill(DynamicalSystemsSet* DSSet, MapOfDSMatrices& DSblocks), Not Yet Implemented");
}


// convert MBlockCSR to numerics structure
void BlockCSRMatrix::convert()
{
  numericsMatSparse->blocknumber0 = nr;
  numericsMatSparse->blocknumber1 = nr;  // nc not always set
  numericsMatSparse->nbblocks = (*MBlockCSR).nnz();
  // Next copies: pointer links!!
  numericsMatSparse->blocksize0 =  &((*diagSizes)[0]);
  numericsMatSparse->blocksize1 =  &((*diagSizes)[0]);  // nr = nc

  // boost
  numericsMatSparse->filled1 = (*MBlockCSR).filled1();
  numericsMatSparse->filled2 = (*MBlockCSR).filled2();
  numericsMatSparse->index1_data = &((*MBlockCSR).index1_data()[0]);
  if (nr > 0)
  {
    numericsMatSparse->index2_data = &((*MBlockCSR).index2_data()[0]);
    numericsMatSparse->block =  &((*MBlockCSR).value_data()[0]);
  };

  //   // Loop through the non-null blocks
  //   for (SpMatIt1 i1 = MBlockCSR->begin1(); i1 != MBlockCSR->end1(); ++i1)
  //     {
  //       for (SpMatIt2 i2 = i1.begin(); i2 != i1.end(); ++i2)
  //  {
  //    block[i] = *i2;
  //  }
  //     }
}

// Display data
void BlockCSRMatrix::display() const
{
  cout << "----- Sparse Block Matrix with "
       << nr << " blocks in a row/col and "
       << MBlockCSR->nnz()
       << " non-null blocks" << endl;
  cout << "filled1:" << MBlockCSR->filled1() << endl;
  cout << "filled2:" << MBlockCSR->filled2() << endl;
  cout << "index1_data:\n";
  print(MBlockCSR->index1_data().begin(), MBlockCSR->index1_data().end());
  cout << endl;
  cout << "index2_data:\n";
  print(MBlockCSR->index2_data().begin(), MBlockCSR->index2_data().end());
  cout << endl;
  cout << "Sum of sizes of the diagonal blocks:"
       << endl;
  print(diagSizes->begin(), diagSizes->end());
}







