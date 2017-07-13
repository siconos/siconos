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
#include <assert.h>
#include "OSNSMatrix.hpp"
#include "NonSmoothLaw.hpp"
#include "Tools.hpp"
#include "BlockCSRMatrix.hpp"
#include "SimulationGraphs.hpp"
#include "SimpleMatrix.hpp"

#include "DynamicalSystem.hpp"


#include "NumericsSparseMatrix.h"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"

// Default constructor: empty matrix, default storage
// No allocation for _M1 or _M2
OSNSMatrix::OSNSMatrix():
  _dimRow(0),  _dimColumn(0), _storageType(NM_DENSE)
{
  _numericsMatrix.reset(new NumericsMatrix);
}

// Constructor with dimensions (one input: square matrix only)
OSNSMatrix::OSNSMatrix(unsigned int n, int stor):
  _dimRow(n),  _dimColumn(n), _storageType(stor)
{
  // Note:

  // for _storageType = NM_DENSE (dense) n represents the real _dimRowension of
  // the matrix and for sparse storage (_storageType == 1) the number
  // of interactionBlocks in a row or column.
  DEBUG_BEGIN("OSNSMatrix::OSNSMatrix(unsigned int n, int stor) \n");
  switch(_storageType)
  {
  case NM_DENSE:
  {
    // A zero matrix M of size nXn is built.  interactionBlocksPositions
    // remains empty (=NULL) since we have no information concerning
    // the Interaction.
    _M1.reset(new SimpleMatrix(n, n));
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    DEBUG_PRINTF(" _M2 is reset with a matrix of size = %i\n", n);
    _M2.reset(new BlockCSRMatrix(n));
    break;
  }
  default: {} // do nothing here
  }
  _numericsMatrix.reset(new NumericsMatrix);
  NM_null(_numericsMatrix.get());
  DEBUG_END("OSNSMatrix::OSNSMatrix(unsigned int n, int stor) \n");
}

OSNSMatrix::OSNSMatrix(unsigned int n, unsigned int m, int stor):
  _dimRow(n),  _dimColumn(m), _storageType(stor)
{
  // Note:

  // for _storageType = NM_DENSE (dense) n represents the real dimension of
  // the matrix and for sparse storage (_storageType == 1) the number
  // of interactionBlocks in a row or column.
  DEBUG_BEGIN("OSNSMatrix::OSNSMatrix(unsigned int n, unsigned int m, int stor)\n");
  switch(_storageType)
  {
  case NM_DENSE:
  {
    // A zero matrix M of size nXn is built.  interactionBlocksPositions
    // remains empty (=NULL) since we have no information concerning
    // the Interaction.
    _M1.reset(new SimpleMatrix(n, n));
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    _M2.reset(new BlockCSRMatrix(n));
    break;
  }
  default: {} // do nothing here
  }

  _numericsMatrix.reset(new NumericsMatrix);
  NM_null(_numericsMatrix.get());
  DEBUG_END("OSNSMatrix::OSNSMatrix(unsigned int n, unsigned int m, int stor)\n");

}

// Basic constructor
OSNSMatrix::OSNSMatrix(InteractionsGraph& indexSet, int stor):
  _dimRow(0), _dimColumn(0), _storageType(stor)
{
  _numericsMatrix.reset(new NumericsMatrix);
  NM_null(_numericsMatrix.get());
  fillW(indexSet);
}


// construct by copy of SiconosMatrix
OSNSMatrix::OSNSMatrix(const SiconosMatrix& MSource):
  _dimRow(MSource.size(0)), _dimColumn(MSource.size(1)), _storageType(NM_DENSE)
{
  _numericsMatrix.reset(new NumericsMatrix);
  NM_null(_numericsMatrix.get());
  _M1.reset(new SimpleMatrix(MSource));
}


// Destructor : pointers are smart
OSNSMatrix::~OSNSMatrix()
{
}

unsigned OSNSMatrix::updateSizeAndPositions(InteractionsGraph& indexSet)
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
  unsigned dim = 0;
  InteractionsGraph::VIterator vd, vdend;
  for (std11::tie(vd, vdend) = indexSet.vertices(); vd != vdend; ++vd)
  {
    assert(indexSet.descriptor(indexSet.bundle(*vd)) == *vd);
    indexSet.properties(*vd).absolute_position = dim;
    dim += (indexSet.bundle(*vd)->nonSmoothLaw()->size());
    DEBUG_PRINTF("Position = %i for interaction %i\n",dim, indexSet.bundle(*vd)->number()  );
    assert(indexSet.properties(*vd).absolute_position < dim);
  }

  return dim;
}
unsigned OSNSMatrix::updateSizeAndPositions(DynamicalSystemsGraph & DSG)
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
  unsigned dim = 0;
  DynamicalSystemsGraph::VIterator dsi, dsend;
  // first loop to compute sizeM and nnz
  for(std11::tie(dsi, dsend) = DSG.vertices(); dsi != dsend; ++dsi)
  {
    SP::DynamicalSystem ds = DSG.bundle(*dsi);
    DSG.properties(*dsi).absolute_position = dim;
    dim += ds->dimension();
  }

  return dim;
}

// Fill the matrix W
void OSNSMatrix::fillW(InteractionsGraph& indexSet, bool update)
{
  DEBUG_BEGIN("void OSNSMatrix::fillW(SP::InteractionsGraph indexSet, bool update)\n");

  if (update)
  {
    // Computes _dimRow and interactionBlocksPositions according to indexSet
    _dimColumn = updateSizeAndPositions(indexSet);
    _dimRow = _dimColumn;
  }

  if (_storageType == NM_DENSE)
  {

    // === Memory allocation, if required ===
    // Mem. is allocate only if !M or if its size has changed.
    if (update)
    {
      if (! _M1)
        _M1.reset(new SimpleMatrix(_dimRow, _dimColumn));
      else
      {
        if (_M1->size(0) != _dimRow || _M1->size(1) != _dimColumn)
          _M1->resize(_dimRow, _dimColumn);
        _M1->zero();
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
    for (std11::tie(vi, viend) = indexSet.vertices();
         vi != viend; ++vi)
    {
      SP::Interaction inter = indexSet.bundle(*vi);
      pos = indexSet.properties(*vi).absolute_position;

      std11::static_pointer_cast<SimpleMatrix>(_M1)
      ->setBlock(pos, pos, *indexSet.properties(*vi).block);
      DEBUG_PRINTF("OSNSMatrix _M1: %i %i\n", _M1->size(0), _M1->size(1));
      DEBUG_PRINTF("OSNSMatrix block: %i %i\n", indexSet.properties(*vi).block->size(0), indexSet.properties(*vi).block->size(1));
    }

    InteractionsGraph::EIterator ei, eiend;
    for (std11::tie(ei, eiend) = indexSet.edges();
         ei != eiend; ++ei)
    {
      InteractionsGraph::VDescriptor vd1 = indexSet.source(*ei);
      InteractionsGraph::VDescriptor vd2 = indexSet.target(*ei);

      SP::Interaction inter1 = indexSet.bundle(vd1);
      SP::Interaction inter2 = indexSet.bundle(vd2);
      pos = indexSet.properties(vd1).absolute_position;

      assert(indexSet.is_vertex(inter2));

      col = indexSet.properties(vd2).absolute_position;

      assert(pos < _dimRow);
      assert(col < _dimColumn);

      DEBUG_PRINTF("OSNSMatrix _M1: %i %i\n", _M1->size(0), _M1->size(1));
      DEBUG_PRINTF("OSNSMatrix upper: %i %i\n", indexSet.properties(*ei).upper_block->size(0), indexSet.properties(*ei).upper_block->size(1));
      DEBUG_PRINTF("OSNSMatrix lower: %i %i\n", indexSet.properties(*ei).lower_block->size(0), indexSet.properties(*ei).lower_block->size(1));

      assert(indexSet.properties(*ei).lower_block);
      assert(indexSet.properties(*ei).upper_block);
      std11::static_pointer_cast<SimpleMatrix>(_M1)
      ->setBlock(std::min(pos, col), std::max(pos, col),
                 *indexSet.properties(*ei).upper_block);

      std11::static_pointer_cast<SimpleMatrix>(_M1)
      ->setBlock(std::max(pos, col), std::min(pos, col),
                 *indexSet.properties(*ei).lower_block);
    }

  }
  else if (_storageType == NM_SPARSE_BLOCK)
  {
    if (! _M2)
    {
      DEBUG_PRINT("Reset _M2 shared pointer using new BlockCSRMatrix(indexSet) \n ");
      _M2.reset(new BlockCSRMatrix(indexSet));

    }
    else
    {
      DEBUG_PRINT("fill existing _M2\n");
      _M2->fill(indexSet);
    }
  }

  if (update)
    convert();
  DEBUG_END("void OSNSMatrix::fill(SP::InteractionsGraph indexSet, bool update)\n");
}

// convert current matrix to NumericsMatrix structure
void OSNSMatrix::convert()
{
  DEBUG_BEGIN("OSNSMatrix::convert()\n");
  DEBUG_PRINTF("_storageType = %i\n", _storageType);
  _numericsMatrix->storageType = _storageType;
  _numericsMatrix->size0 = _dimRow;
  _numericsMatrix->size1 = _dimColumn;
  switch (_storageType)
  {
  case NM_DENSE:
  {
    _numericsMatrix->matrix0 = _M1->getArray(); // Pointer link
    // _numericsMatrix->matrix1 = NULL; matrix1 is not set to NULL: we
    // keep previous allocation. May be usefull if we switch between
    // different storages during simu
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    _M2->convert();
    _numericsMatrix->matrix1 = &*_M2->getNumericsMatSparse();
    break;
  }
  case NM_SPARSE:
  {
    // we already filled the matrix
    break;
  }
  default:
  {
     RuntimeException::selfThrow("OSNSMatrix::convert unknown _storageType");
  }
  }
  DEBUG_END("OSNSMatrix::convert()\n");
}

// Fill the matrix M
void OSNSMatrix::fillM(DynamicalSystemsGraph & DSG, bool update)
{
  DEBUG_BEGIN("void OSNSMatrix::fillM(SP::DynamicalSystemsGraph DSG, bool update)\n");

  if (update)
    {
      _dimColumn = updateSizeAndPositions(DSG);
      _dimRow = _dimColumn;
    }

  switch (_storageType)
  {
  case NM_SPARSE:
  {
    if (update)
    {
      size_t sizeM = _dimRow;
      DEBUG_PRINTF("sizeM = %lu \n", sizeM);

      // We choose a triplet matrix format for inserting values.
      // This simplifies the memory manipulation.
      NumericsMatrix& M_NM = *numericsMatrix();
      NM_clearSparse(&M_NM);
      M_NM.storageType = NM_SPARSE;
      M_NM.size0 = sizeM;
      M_NM.size1 = sizeM;
      NM_triplet_alloc(&M_NM, sizeM); // At least one element per row
      M_NM.matrix2->origin = NS_TRIPLET;
      CSparseMatrix* Mtriplet = NM_triplet(&M_NM);

      unsigned int pos =0;
      // Loop over the DS for filling M
      DynamicalSystemsGraph::VIterator dsi, dsend;
      for(std11::tie(dsi, dsend) = DSG.vertices(); dsi != dsend; ++dsi)
      {
        SP::DynamicalSystem ds = DSG.bundle(*dsi);
        SiconosMatrix* W = DSG.properties(*dsi).W.get();
        pos = DSG.properties(*dsi).absolute_position;
        W->fillTriplet(Mtriplet, pos, pos);
        DEBUG_PRINTF("pos = %u \n", pos);
      }
    }
    DEBUG_EXPR(NM_display(numericsMatrix().get() ););
    break;
  }
  default:
  {
     RuntimeException::selfThrow("OSNSMatrix::convert unknown _storageType");
  }
  }


  DEBUG_END("void OSNSMatrix::fillM(SP::DynamicalSystemsGraph DSG, bool update)\n");
}

// Fill the matrix H
void OSNSMatrix::fillH(DynamicalSystemsGraph & DSG, InteractionsGraph& indexSet, bool update)
{
  DEBUG_BEGIN("void OSNSMatrix::fillH(SP::DynamicalSystemsGraph DSG, InteractionsGraph& indexSet, bool update)\n");
  if (update)
    {

      _dimColumn = updateSizeAndPositions(indexSet);
      _dimRow = updateSizeAndPositions(DSG);
    }

  switch (_storageType)
  {
  case NM_SPARSE:
  {
    if (update)
    {
      // We choose a triplet matrix format for inserting values.
      // This simplifies the memory manipulation.
      NumericsMatrix& H_NM = *numericsMatrix();
      NM_clearSparse(&H_NM);
      H_NM.storageType = NM_SPARSE;
      H_NM.size0 = _dimRow;
      H_NM.size1 = _dimColumn;
      NM_triplet_alloc(&H_NM, _dimColumn);// At least one element per column
      H_NM.matrix2->origin = NS_TRIPLET;
      CSparseMatrix* Htriplet= NM_triplet(&H_NM);
      unsigned int pos = 0, pos_ds=0;


      SP::SiconosMatrix leftInteractionBlock;
      InteractionsGraph::VIterator ui, uiend;
      for (std11::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui)
      {
        Interaction& inter = *indexSet.bundle(*ui);
        VectorOfSMatrices& workMInter = *indexSet.properties(*ui).workMatrices;

        SP::DynamicalSystem ds1 = indexSet.properties(*ui).source;
        SP::DynamicalSystem ds2 = indexSet.properties(*ui).target;

        bool endl = false;
        size_t posBlock = indexSet.properties(*ui).source_pos;
        size_t pos2 = indexSet.properties(*ui).target_pos;

        pos =  indexSet.properties(*ui).absolute_position;
        for (SP::DynamicalSystem ds = ds1; !endl; ds = ds2, posBlock = pos2)
        {
          endl = (ds == ds2);
          size_t sizeDS = ds->dimension();
          // this whole part is a hack. Just should just get the rightblock
          leftInteractionBlock.reset(new SimpleMatrix(3, sizeDS));
          inter.getLeftInteractionBlockForDS(posBlock, leftInteractionBlock, workMInter);
          leftInteractionBlock->trans();
          pos_ds =  DSG.properties(DSG.descriptor(ds)).absolute_position;
          DEBUG_PRINTF("pos = %u", pos);
          DEBUG_PRINTF("pos_ds = %u", pos_ds);
          leftInteractionBlock->fillTriplet(Htriplet, pos_ds, pos);
        }
      }
      break;

    }
  }
  default:
  {
     RuntimeException::selfThrow("OSNSMatrix::convert unknown _storageType");
  }
  }
  DEBUG_END("void OSNSMatrix::fillH(SP::DynamicalSystemsGraph DSG, InteractionsGraph& indexSet, bool update)\n");
}



// Display data
void OSNSMatrix::display() const
{
  if (_storageType == NM_DENSE)
  {
    std::cout << "----- OSNS Matrix using default storage type for Numerics structure (SiconosMatrix -> double*)" <<std::endl;
    if (! _M1)
      std::cout << " matrix = NULL pointer" <<std::endl;
    else _M1->display();
  }
  else if (_storageType == NM_SPARSE_BLOCK)
  {
    std::cout << "----- OSNS Matrix using Sparse InteractionBlock storage type for Numerics (SparseBlockStructuredMatrix)" <<std::endl;
    if (! _M2)
      std::cout << " matrix = NULL pointer" <<std::endl;
    else _M2->display();
  }
  else if (_storageType == NM_SPARSE)
  {
    std::cout << "----- OSNS Matrix using sparse storage, nothing to show" << std::endl;
  }
}
