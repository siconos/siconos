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
#include <assert.h>
#include "OSNSMatrix.hpp"
#include "Tools.hpp"
//#define OSNSM_DEBUG
using namespace std;

void OSNSMatrix::updateSizeAndPositions(unsigned int& dim,
                                        SP::InteractionsGraph indexSet)
{
  // === Description ===

  // For a interactionBlock (diagonal or extra diagonal) corresponding to
  // a Interaction, we need to know the position of its first
  // element in the full-matrix M. This position dep`ends on the
  // previous interactionBlocks sizes.
  //
  // positions are saved in a map<SP::Interaction, unsigned int>,
  // named interactionBlocksPositions.
  //

  // Computes real size of the current matrix = sum of the dim. of all
  // Interactionin indexSet
  dim = 0;
  InteractionsGraph::VIterator vd, vdend;
  for (boost::tie(vd, vdend) = indexSet->vertices(); vd != vdend; ++vd)
  {
    assert(indexSet->descriptor(indexSet->bundle(*vd)) == *vd);

    //    (*interactionBlocksPositions)[indexSet->bundle(*vd)] = dim;
    indexSet->bundle(*vd)->setAbsolutePosition(dim);
    dim += (indexSet->bundle(*vd)->getNonSmoothLawSize());

    assert(indexSet->bundle(*vd)->absolutePosition() < dim);
  }
}

void OSNSMatrix::updateSizeAndPositions(unsigned int& dim,
                                        SP::DynamicalSystemsSet DSSet)
{
  // === Description ===

  // For a DSBlock (diagonal or extra diagonal) corresponding to a
  // DynamicalSet, we need to know the position of its first element
  // in the full-matrix M. This position depends on the previous
  // DSBlocks sizes.
  //
  // positions are saved in a map<SP::DynamicalSystem, unsigned int>,
  // named DSBlocksPositions.
  //

  // Computes real size of the current matrix = sum of the dim. of all
  // Interactionin indexSet
  dim = 0;
  for (DSIterator it = DSSet->begin(); it != DSSet->end(); ++it)
  {
    (*DSBlocksPositions)[*it] = dim;
    dim += (*it)->getDim();
  }
}

void OSNSMatrix::updateSizeAndPositions(unsigned int& dim,
                                        SP::DynamicalSystemsSet DSSet,
                                        SP::InteractionsGraph indexSet)
{
  // === Description ===

  // positions are saved in a map<SP::Interaction, unsigned int>,
  // named interactionBlocksPositions.  positions are saved in a
  // map<SP::DynamicalSystem, unsigned int>, named DSBlocksPositions.
  //

  // Computes real size of the current matrix = sum of the dim. of all
  // Interactionin indexSet
  dim = 0;
  for (DSIterator it = DSSet->begin(); it != DSSet->end(); ++it)
  {
    (*DSBlocksPositions)[*it] = dim;
    dim += (*it)->getDim();
  }
  InteractionsGraph::VIterator vd, vdend;
  for (boost::tie(vd, vdend) = indexSet->vertices(); vd != vdend; ++vd)
  {
    indexSet->bundle(*vd)->setAbsolutePosition(dim);
    dim += indexSet->bundle(*vd)->getNonSmoothLawSize();
  }
}

// Default constructor: empty matrix, default storage
// No allocation for M1 or M2
OSNSMatrix::OSNSMatrix():
  dimRow(0),  dimColumn(0), storageType(0)
{
  //  interactionBlocksPositions.reset(new Interaction_int());
  DSBlocksPositions.reset(new DS_int());
  numericsMat.reset(new NumericsMatrix);
}

// Constructor with dimensions (one input: square matrix only)
OSNSMatrix::OSNSMatrix(unsigned int n, int stor):
  dimRow(n),  dimColumn(n), storageType(stor)
{
  // Note:

  // for storageType = 0 (dense) n represents the real dimRowension of
  // the matrix and for sparse storage (storageType == 1) the number
  // of interactionBlocks in a row or column.

  if (storageType == 0)
  {
    // A zero matrix M of size nXn is built.  interactionBlocksPositions
    // remains empty (=NULL) since we have no information concerning
    // the Interaction.
    M1.reset(new SimpleMatrix(n, n));
  }
  else // if(storageType == 1)
  {
    M2.reset(new BlockCSRMatrix(n));
  }
  //  interactionBlocksPositions.reset(new Interaction_int());
  DSBlocksPositions.reset(new DS_int());
  numericsMat.reset(new NumericsMatrix);
}
OSNSMatrix::OSNSMatrix(unsigned int n, unsigned int m, int stor):
  dimRow(n),  dimColumn(m), storageType(stor)
{
  // Note:

  // for storageType = 0 (dense) n represents the real dimension of
  // the matrix and for sparse storage (storageType == 1) the number
  // of interactionBlocks in a row or column.

  if (storageType == 0)
  {
    // A zero matrix M of size nXn is built.  interactionBlocksPositions
    // remains empty (=NULL) since we have no information concerning
    // the Interaction.
    M1.reset(new SimpleMatrix(n, m));
  }
  else // if(storageType == 1)
    M2.reset(new BlockCSRMatrix(n));

  //  interactionBlocksPositions.reset(new Interaction_int());
  DSBlocksPositions.reset(new DS_int());
  numericsMat.reset(new NumericsMatrix);
}

// Basic constructor
OSNSMatrix::OSNSMatrix(SP::InteractionsGraph indexSet, int stor):
  dimRow(0), dimColumn(0), storageType(stor)
{
  //  interactionBlocksPositions.reset(new Interaction_int());
  DSBlocksPositions.reset(new DS_int());
  numericsMat.reset(new NumericsMatrix);

  fill(indexSet);
}

OSNSMatrix::OSNSMatrix(SP::DynamicalSystemsSet DSSet, MapOfDSMatrices& DSBlocks, int stor):
  dimRow(0), dimColumn(0), storageType(stor)
{
  //  interactionBlocksPositions.reset(new Interaction_int());
  //  DSBlocksPositions.reset(new DS_int());
  numericsMat.reset(new NumericsMatrix);

  fill(DSSet, DSBlocks);
}
OSNSMatrix::OSNSMatrix(SP::DynamicalSystemsSet DSSet, SP::InteractionsGraph indexSet, MapOfDSMapOfInteractionMatrices& DSInteractionBlocks, int stor):
  dimRow(0), dimColumn(0), storageType(stor)
{
  //  interactionBlocksPositions.reset(new Interaction_int());
  DSBlocksPositions.reset(new DS_int());
  numericsMat.reset(new NumericsMatrix);

  fill(DSSet, indexSet, DSInteractionBlocks);
}

OSNSMatrix::OSNSMatrix(SP::InteractionsGraph indexSet, SP::DynamicalSystemsSet DSSet,   MapOfDSMatrices& DSBlocks, MapOfDSMapOfInteractionMatrices& DSInteractionBlocks, MapOfInteractionMapOfDSMatrices& interactionDSBlocks, int stor):
  dimRow(0), dimColumn(0), storageType(stor)
{
  //  interactionBlocksPositions.reset(new Interaction_int());
  DSBlocksPositions.reset(new DS_int());
  numericsMat.reset(new NumericsMatrix);

  fill(indexSet, DSSet, DSBlocks, DSInteractionBlocks, interactionDSBlocks);
}

// Copy of a SiconosMatrix (used when OSNS xml constructor is called with M input in XML file)
OSNSMatrix::OSNSMatrix(const SiconosMatrix& MSource):
  dimRow(MSource.size(0)), dimColumn(MSource.size(1)), storageType(0)
{
  //  interactionBlocksPositions.reset(new Interaction_int());
  DSBlocksPositions.reset(new DS_int());
  numericsMat.reset(new NumericsMatrix);
  M1.reset(new SimpleMatrix(MSource));

  // Warning: interactionBlocksPositions remains empty since we have no
  // information concerning indexSet and interactionBlocks in MSource
}


// Destructor : pointers are smart
OSNSMatrix::~OSNSMatrix()
{
}

unsigned int OSNSMatrix::getPositionOfInteractionBlock(SP::Interaction inter) const
{
  return inter->absolutePosition();
  /*  Interaction_int::const_iterator it = interactionBlocksPositions->find(inter);
  if (it== interactionBlocksPositions->end())
    RuntimeException::selfThrow
    ("OSNSMatrix::getPositionOfInteractionBlock(Interaction inter), inter does not belong to the index set used to built the OSNS matrix.");
    return it->second;*/
}
unsigned int OSNSMatrix::getPositionOfDSBlock(SP::DynamicalSystem DS) const
{
  DS_int::const_iterator it = DSBlocksPositions->find(DS);
  if (it == DSBlocksPositions->end())
    RuntimeException::selfThrow
    ("OSNSMatrix::getPositionOfDSBlock(DynamicalSystems ds), ds does not belong to the DynamicalSet used to built the OSNS matrix.");
  return it->second;
}

// Fill the matrix
void OSNSMatrix::fill(SP::InteractionsGraph indexSet, bool update)
{

  assert(indexSet);

  if (update)
  {
    // Computes dimRow and interactionBlocksPositions according to indexSet
    updateSizeAndPositions(dimColumn, indexSet);
    dimRow = dimColumn;
  }

  if (storageType == 0)
  {

    // === Memory allocation, if required ===
    // Mem. is allocate only if !M or if its size has changed.
    if (update)
    {
      if (! M1)
        M1.reset(new SimpleMatrix(dimRow, dimColumn));
      else
      {
        if (M1->size(0) != dimRow || M1->size(1) != dimColumn)
          M1->resize(dimRow, dimColumn);
        M1->zero();
      }
    }

    // ======> Aim: find inter1 and inter2 both in indexSet and which have
    // common DynamicalSystems.  Then get the corresponding matrix
    // from map interactionBlocks, and copy it into M

    unsigned int pos = 0, col = 0; // index position used for
    // interactionBlock copy into M, see
    // below.
    // === Loop through "active" Interactions (ie present in
    // indexSets[level]) ===
    InteractionsGraph::VIterator vi, viend;
    for (boost::tie(vi, viend) = indexSet->vertices();
         vi != viend; ++vi)
    {
      SP::Interaction inter = indexSet->bundle(*vi);
      pos = inter->absolutePosition();

      boost::static_pointer_cast<SimpleMatrix>(M1)
      ->setBlock(pos, pos, *indexSet->properties(*vi).block);
#ifdef OSNSMPROJ_DEBUG
      printf("OSNSMatrix M1: %i %i\n", M1->size(0), M1->size(1));
      printf("OSNSMatrix block: %i %i\n", indexSet->properties(*vi).block->size(0), indexSet->properties(*vi).block->size(1));
#endif
    }


    InteractionsGraph::EIterator ei, eiend;
    for (boost::tie(ei, eiend) = indexSet->edges();
         ei != eiend; ++ei)
    {
      InteractionsGraph::VDescriptor vd1 = indexSet->source(*ei);
      InteractionsGraph::VDescriptor vd2 = indexSet->target(*ei);

      SP::Interaction inter1 = indexSet->bundle(vd1);
      SP::Interaction inter2 = indexSet->bundle(vd2);

      pos = inter1->absolutePosition();

      assert(indexSet->is_vertex(inter2));

      col = inter2->absolutePosition();


      assert(pos < dimRow);
      assert(col < dimColumn);

#ifdef OSNSM_DEBUG
      printf("OSNSMatrix M1: %i %i\n", M1->size(0), M1->size(1));
      printf("OSNSMatrix upper: %i %i\n", indexSet->properties(*ei).upper_block->size(0), indexSet->properties(*ei).upper_block->size(1));
      printf("OSNSMatrix lower: %i %i\n", indexSet->properties(*ei).lower_block->size(0), indexSet->properties(*ei).lower_block->size(1));
#endif

      assert(indexSet->properties(*ei).lower_block);
      assert(indexSet->properties(*ei).upper_block);
      boost::static_pointer_cast<SimpleMatrix>(M1)
      ->setBlock(std::min(pos, col), std::max(pos, col),
                 *indexSet->properties(*ei).upper_block);

      boost::static_pointer_cast<SimpleMatrix>(M1)
      ->setBlock(std::max(pos, col), std::min(pos, col),
                 *indexSet->properties(*ei).lower_block);
    }

  }
  else // if storageType == 1
  {
    if (! M2)
      M2.reset(new BlockCSRMatrix(indexSet));
    else
      M2->fill(indexSet);
  }
  if (update)
    convert();

}
void OSNSMatrix::fillDiagonal(SP::InteractionsGraph IG, bool update)
{

  assert(0);

}
/*
  assert(URSet);

  if (update)
  {
    // Computes dimRow and interactionBlocksPositions according to indexSet
    updateSizeAndPositions(dimColumn, URSet);
    updateSizeAndPositions(dimRow, URSet);
  }

  if (storageType==0)
  {
    // === Memory allocation, if required ===
    // Mem. is allocate only if !M or if its size has changed.
    if (update)
    {
      if (! M1)
        M1.reset(new SimpleMatrix(dimRow,dimColumn));
      else
      {
        if (M1->size(0) != dimRow || M1->size(1)!= dimColumn)
          M1->resize(dimRow,dimColumn);
        M1->zero();
      }
    }

    // ======> Aim: find inter1 and inter2 both in indexSet and which have
    // common DynamicalSystems.  Then get the corresponding matrix
    // from map interactionBlocks, and copy it into M

    unsigned int pos = 0, col =0; // index position used for
    // interactionBlock copy into M, see
    // below.
    // === Loop through "active" Interactions (ie present in
    // indexSets[level]) ===

    InteractionsGraph::VIterator ui,uiend;
    for (boost::tie(ui,uiend)=URSet->vertices(); ui!=uiend; ++ui)
    {
      SP::Interaction inter = URSet->bundle(*ui);
      pos = (*interactionBlocksPositions)[inter];
      boost::static_pointer_cast<SimpleMatrix>(M1)->setBlock(pos,pos,*(interactionBlocks[inter][inter]));
    }
  }
  else // if storageType == 1
  {
    if (! M2)
      M2.reset(new BlockCSRMatrix(URSet,interactionBlocks));
    else
      M2->fill(URSet,interactionBlocks);
  }

  if (update)
    convert();
    }*/

void OSNSMatrix::fill(SP::DynamicalSystemsSet DSSet, MapOfDSMatrices& DSBlocks, bool update)
{

  assert(DSSet && "NULL pointer");

  if (update)
  {
    // Computes dimRow and interactionBlocksPositions according to indexSet
    updateSizeAndPositions(dimRow, DSSet);
    updateSizeAndPositions(dimColumn, DSSet);
  }
  if (storageType == 0)
  {
    // === Memory allocation, if required ===
    // Mem. is allocate only if !M or if its size has changed.
    if (update)
    {
      if (! M1)
        M1.reset(new SimpleMatrix(dimRow, dimColumn));
      else
      {
        if (M1->size(0) != dimRow || M1->size(1) != dimColumn)
          M1->resize(dimRow, dimColumn);
        M1->zero();
      }
    }

    // get the corresponding matrix from map interactionBlocks, and copy it into M

    unsigned int pos = 0; // index position used for interactionBlock copy into M, see below.
    // === Loop through "active" Interactions (ie present in indexSets[level]) ===
    for (DSIterator itDS = DSSet->begin(); itDS != DSSet->end(); ++itDS)
    {
      // DS = *itDS
      // corresponding matrix = DSBlocks[*itDS]

      // Case 1: basic storage
      pos = (*DSBlocksPositions)[*itDS];
      boost::static_pointer_cast<SimpleMatrix>(M1)->setBlock(pos, pos, *(DSBlocks[(*itDS)->number()]));
    }
  }
  else // if storageType == 1
  {
    if (! M2)
      M2.reset(new BlockCSRMatrix(DSSet, DSBlocks));
    else
      M2->fill(DSSet, DSBlocks);
  }
  if (update)
    convert();

}
void OSNSMatrix::fill(SP::DynamicalSystemsSet DSSet, SP::InteractionsGraph IG, MapOfDSMapOfInteractionMatrices& DSInteractionBlocks, bool update)
{

  assert(IG && "NULL pointer");
  assert(DSSet && "NULL pointer");

  if (update)
  {
    // Computes dimRow and interactionBlocksPositions and  DSBlocksPositions according to indexSet
    updateSizeAndPositions(dimColumn, IG);


    updateSizeAndPositions(dimRow, DSSet);
  }

  if (storageType == 0)
  {
    // === Memory allocation, if required ===
    // Mem. is allocate only if !M or if its size has changed.
    if (update)
    {
      if (! M1)
        M1.reset(new SimpleMatrix(dimRow, dimColumn));
      else
      {
        if (M1->size(0) != dimRow || M1->size(1) != dimColumn)
          M1->resize(dimRow, dimColumn);
        M1->zero();
      }
    }


    // Get the  matrix from map interactionDSBlocks which corresponds to Interactionand DS, and copy it into M

    unsigned int pos = 0, col = 0; // index position used for interactionBlock copy into M, see below.
    // === Loop through "active" Interactions (ie present in indexSets[level]) ===


    for (DSIterator itCol = DSSet->begin(); itCol != DSSet->end(); ++itCol)
    {
      // Look for DS
      // The connection is checked thanks to interactionBlocks map
      for (InteractionMatrixColumnIterator itRow = DSInteractionBlocks[*itCol].begin(); itRow != DSInteractionBlocks[*itCol].end(); ++itRow)
      {
        // DS = *itCol
        // Interaction= *itRow
        // corresponding matrix = DSInteractionBlocks[*itCol][(*itRow).first]

        // Case 1: basic storage
        pos = (*DSBlocksPositions)[*itCol];
        col = ((*itRow).first)->absolutePosition();
        boost::static_pointer_cast<SimpleMatrix>(M1)->setBlock(pos, col, *(DSInteractionBlocks[*itCol][(*itRow).first]));
      }
    }

  }
  else // if storageType == 1
  {
    RuntimeException::selfThrow("Not yet Implemented case storageType == 1:OSNSMatrix::fill(DynamicalSystemsSet* DSSet, SP::InteractionsGraph IG, MapOfDSMapOfInteractionMatrices& DSInteractionBlocks,bool update) ");

    //       if(M2==NULL)
    //  M2 = new BlockCSRMatrix(DSSet, IG, DSInteractionBlocks);
    //       else
    //  M2->fill(DSSet, IG ,DSInteractionBlocks);
  }
  if (update)
    convert();

}




void OSNSMatrix::fill(SP::InteractionsGraph indexSet, SP::DynamicalSystemsSet DSSet,  MapOfInteractionMapOfDSMatrices& interactionDSBlocks, bool update)
{
  assert(false);
}
void OSNSMatrix::fill(SP::InteractionsGraph IG, SP::DynamicalSystemsSet DSSet,  MapOfDSMatrices& DSBlocks, MapOfDSMapOfInteractionMatrices& DSInteractionBlocks ,  MapOfInteractionMapOfDSMatrices& interactionDSBlocks , bool update)
{
  assert(false);
}

// convert current matrix to NumericsMatrix structure
void OSNSMatrix::convert()
{
  numericsMat->storageType = storageType;
  numericsMat->size0 = dimRow;
  numericsMat->size1 = dimColumn;
  if (storageType == 0)
  {
    numericsMat->matrix0 = M1->getArray(); // Pointer link
    // numericsMat->matrix1 = NULL; matrix1 is not set to NULL: we
    // keep previous allocation. May be usefull if we switch between
    // different storages during simu
  }
  else
  {
    M2->convert();
    numericsMat->matrix1 = &*M2->getNumericsMatSparse();
  }
}

// Display data
void OSNSMatrix::display() const
{
  if (storageType == 0)
  {
    cout << "----- OSNS Matrix using default storage type for Numerics structure (SiconosMatrix -> double*)" << endl;
    if (! M1)
      cout << " matrix = NULL pointer" << endl;
    else M1->display();
  }
  else
  {
    cout << "----- OSNS Matrix using Sparse InteractionBlock storage type for Numerics (SparseBlockStructuredMatrix)" << endl;
    if (! M2)
      cout << " matrix = NULL pointer" << endl;
    else M2->display();
  }
}
