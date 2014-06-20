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
#include <assert.h>
#include "OSNSMatrix.hpp"
#include "NonSmoothLaw.hpp"
#include "Tools.hpp"
#include "BlockCSRMatrix.hpp"
//#define OSNSM_DEBUG


void OSNSMatrix::updateSizeAndPositions(unsigned int& dim,
                                        SP::InteractionsGraph indexSet)
{
  // === Description ===

  // For an interactionBlock (diagonal or extra diagonal) corresponding to
  // an Interaction, we need to know the position of its first
  // element in the full-matrix M. This position depends on the
  // previous interactionBlocks sizes.
  //
  // Note FP: at the time positions are saved in the Interaction
  // but this is wrong (I think) since it prevents the inter
  // to be present in several different osns.
  //

  // Computes real size of the current matrix = sum of the dim. of all
  // Interactionin indexSet
  dim = 0;
  InteractionsGraph::VIterator vd, vdend;
  for (std11::tie(vd, vdend) = indexSet->vertices(); vd != vdend; ++vd)
  {
    assert(indexSet->descriptor(indexSet->bundle(*vd)) == *vd);

    indexSet->bundle(*vd)->setAbsolutePosition(dim); 
    dim += (indexSet->bundle(*vd)->nonSmoothLaw()->size());

    assert(indexSet->bundle(*vd)->absolutePosition() < dim);
  }
}

// Default constructor: empty matrix, default storage
// No allocation for M1 or M2
OSNSMatrix::OSNSMatrix():
  dimRow(0),  dimColumn(0), storageType(0)
{
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

  numericsMat.reset(new NumericsMatrix);
}

// Basic constructor
OSNSMatrix::OSNSMatrix(SP::InteractionsGraph indexSet, int stor):
  dimRow(0), dimColumn(0), storageType(stor)
{
  numericsMat.reset(new NumericsMatrix);
  fill(indexSet);
}


// construct by copy of SiconosMatrix
OSNSMatrix::OSNSMatrix(const SiconosMatrix& MSource):
  dimRow(MSource.size(0)), dimColumn(MSource.size(1)), storageType(0)
{
  numericsMat.reset(new NumericsMatrix);
  M1.reset(new SimpleMatrix(MSource));
}


// Destructor : pointers are smart
OSNSMatrix::~OSNSMatrix()
{
}

unsigned int OSNSMatrix::getPositionOfInteractionBlock(SP::Interaction inter) const
{
  // Note FP: I think the return value below is not the right one :
  // this position does not depend on the interaction but on
  // the OSNS and the corresponding indexSet. 
  // One Interaction may have different absolute positions if it is present
  // in several OSNS. ==> add this pos as a property on vertex in Interactions Graph
  // 
  return inter->absolutePosition();
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
    for (std11::tie(vi, viend) = indexSet->vertices();
         vi != viend; ++vi)
    {
      SP::Interaction inter = indexSet->bundle(*vi);
      pos = inter->absolutePosition();

      std11::static_pointer_cast<SimpleMatrix>(M1)
      ->setBlock(pos, pos, *indexSet->properties(*vi).block);
#ifdef OSNSMPROJ_DEBUG
      printf("OSNSMatrix M1: %i %i\n", M1->size(0), M1->size(1));
      printf("OSNSMatrix block: %i %i\n", indexSet->properties(*vi).block->size(0), indexSet->properties(*vi).block->size(1));
#endif
    }

    InteractionsGraph::EIterator ei, eiend;
    for (std11::tie(ei, eiend) = indexSet->edges();
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
      std11::static_pointer_cast<SimpleMatrix>(M1)
      ->setBlock(std::min(pos, col), std::max(pos, col),
                 *indexSet->properties(*ei).upper_block);

      std11::static_pointer_cast<SimpleMatrix>(M1)
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
    std::cout << "----- OSNS Matrix using default storage type for Numerics structure (SiconosMatrix -> double*)" <<std::endl;
    if (! M1)
      std::cout << " matrix = NULL pointer" <<std::endl;
    else M1->display();
  }
  else
  {
    std::cout << "----- OSNS Matrix using Sparse InteractionBlock storage type for Numerics (SparseBlockStructuredMatrix)" <<std::endl;
    if (! M2)
      std::cout << " matrix = NULL pointer" <<std::endl;
    else M2->display();
  }
}
