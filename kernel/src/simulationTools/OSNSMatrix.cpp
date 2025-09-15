/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#include "NumericsMatrix.h"
#include "OSNSMatrix.hpp"
#include "NonSmoothLaw.hpp"
#include "Tools.hpp"
#include "BlockCSRMatrix.hpp"
#include "SimulationGraphs.hpp"
#include "SimpleMatrix.hpp"
#include "Interaction.hpp"
#include "DynamicalSystem.hpp"
#include "NumericsSparseMatrix.h"
#include "CSparseMatrix.h"
#include "OneStepIntegrator.hpp"
#include "MoreauJeanOSI.hpp"
#include "MoreauJeanGOSI.hpp"
#include "SecondOrderDS.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

// Default constructor: empty matrix, default storage
// No allocation for _M1 or _M2
OSNSMatrix::OSNSMatrix():
  _dimRow(0),  _dimColumn(0), _storageType(NM_DENSE)
{
  //_numericsMatrix.reset(new NumericsMatrix);
}

// Constructor with dimensions (one input: square matrix only)
OSNSMatrix::OSNSMatrix(unsigned int n, NM_types stor):
  _dimRow(n),  _dimColumn(n), _storageType(stor)
{
  // Note:
  // * Dense matrix (_storageType = NM_DENSE), n represents the real dimension of
  // the matrix
  // * Sparse matrix (_storageType == 1) n represents the number of blocks in a row or column.

  DEBUG_BEGIN("OSNSMatrix::OSNSMatrix(unsigned int n, NM_types stor) \n");
  switch(_storageType)
  {
  case NM_DENSE:
  {
    // A zero matrix M of size nXn is built.
    _M1.reset(new SimpleMatrix(n, n));
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    DEBUG_PRINTF(" _M2 is reset with a matrix of size = %i\n", n);
    _M2.reset(new BlockCSRMatrix(n));
    break;
  }
  default:
  {
    _triplet_nzmax = _dimRow ; /* at least a non zero element per row */
  } // do nothing here
  }

  DEBUG_END("OSNSMatrix::OSNSMatrix(unsigned int n, int stor) \n");

}

OSNSMatrix::OSNSMatrix(unsigned int n, unsigned int m, NM_types stor):
  _dimRow(n),  _dimColumn(m), _storageType(stor)
{
  // Note:

  // for _storageType = NM_DENSE (dense) n represents the real dimension of
  // the matrix and for sparse storage (_storageType == 1) the number
  // of interactionBlocks in a row or column.
  DEBUG_BEGIN("OSNSMatrix::OSNSMatrix(unsigned int n, unsigned int m, NM_types stor)\n");
  switch(_storageType)
  {
  case NM_DENSE:
  {
    // A zero matrix M of size nXn is built.  interactionBlocksPositions
    // remains empty (=nullptr) since we have no information concerning
    // the Interaction.
    _M1.reset(new SimpleMatrix(n, n));
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    _M2.reset(new BlockCSRMatrix(n));
    break;
  }
  default:
  {} // do nothing here
  }


  DEBUG_END("OSNSMatrix::OSNSMatrix(unsigned int n, unsigned int m, NM_types stor)\n");

}

// Build from index set (i.e. get size from number of interactions in the set)
OSNSMatrix::OSNSMatrix(InteractionsGraph& indexSet, NM_types stor):
  _dimRow(0), _dimColumn(0), _storageType(stor)
{
  DEBUG_BEGIN("OSNSMatrix::OSNSMatrix(InteractionsGraph& indexSet, NM_types stor)\n");
//  _numericsMatrix.reset(new NumericsMatrix);
//  NM_null(_numericsMatrix.get());
  fillM(indexSet);
  DEBUG_END("OSNSMatrix::OSNSMatrix(InteractionsGraph& indexSet, NM_types stor)\n");
}


// construct by copy of SiconosMatrix
OSNSMatrix::OSNSMatrix(const SiconosMatrix& MSource):
  _dimRow(MSource.size(0)), _dimColumn(MSource.size(1)), _storageType(NM_DENSE)
{
//  _numericsMatrix.reset(new NumericsMatrix);
//  NM_null(_numericsMatrix.get());
  _M1.reset(new SimpleMatrix(MSource));
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
  for(std::tie(vd, vdend) = indexSet.vertices(); vd != vdend; ++vd)
  {
    assert(indexSet.descriptor(indexSet.bundle(*vd)) == *vd);
    indexSet.properties(*vd).absolute_position = dim;
    dim += (indexSet.bundle(*vd)->nonSmoothLaw()->size());
    DEBUG_PRINTF("Position = %i for interaction %zu\n",dim, indexSet.bundle(*vd)->number());
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

  // Computes real size of the current matrix = sum of the dim. of all
  // Interactionin indexSet
  unsigned dim = 0;
  DynamicalSystemsGraph::VIterator dsi, dsend;
  // first loop to compute sizeM and nnz
  for(std::tie(dsi, dsend) = DSG.vertices(); dsi != dsend; ++dsi)
  {
    SP::DynamicalSystem ds = DSG.bundle(*dsi);
    DSG.properties(*dsi).absolute_position = dim;
    dim += ds->dimension();
  }

  return dim;
}

// Fill the matrix W
void OSNSMatrix::fillM(InteractionsGraph& indexSet, bool update)
{
  DEBUG_BEGIN("void OSNSMatrix::fillM(SP::InteractionsGraph indexSet, bool update)\n");
  DEBUG_PRINTF(" update = %i\n", update);
  if(update)  // If index set vertices list has changed
  {
    // Computes _dimRow and interactionBlocksPositions according to indexSet
    _dimColumn = updateSizeAndPositions(indexSet);
    _dimRow = _dimColumn;
  }

  if(_storageType == NM_DENSE)
  {

    // === Memory allocation, if required ===
    // Mem. is allocate only if !M or if its size has changed.
    if(update)
    {
      if(! _M1)
        _M1.reset(new SimpleMatrix(_dimRow, _dimColumn));
      else
      {
        if(_M1->size(0) != _dimRow || _M1->size(1) != _dimColumn)
          _M1->resize(_dimRow, _dimColumn);
        _M1->zero();
      }
    }

    unsigned int pos = 0, col = 0; // index position used for

    // === Loop through "active" Interactions (ie vertices present in indexSets[level]) ===
    InteractionsGraph::VIterator vi, viend;
    for(std::tie(vi, viend) = indexSet.vertices();
        vi != viend; ++vi)
    {
      SP::Interaction inter = indexSet.bundle(*vi);
      pos = indexSet.properties(*vi).absolute_position;

      std::static_pointer_cast<SimpleMatrix>(_M1)
      ->setBlock(pos, pos, *indexSet.properties(*vi).block);
      DEBUG_PRINTF("OSNSMatrix _M1: %i %i\n", _M1->size(0), _M1->size(1));
      DEBUG_PRINTF("OSNSMatrix block: %i %i\n", indexSet.properties(*vi).block->size(0), indexSet.properties(*vi).block->size(1));
    }

    // == Loop through all edges (ds) in active index set ==
    // Computation of extra-diagonal blocks.
    InteractionsGraph::EIterator ei, eiend;
    for(std::tie(ei, eiend) = indexSet.edges();
        ei != eiend; ++ei)
    {
      // For current edge (ds) get source and target vertices (interactions)
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
      std::static_pointer_cast<SimpleMatrix>(_M1)
      ->setBlock(std::min(pos, col), std::max(pos, col),
                 *indexSet.properties(*ei).upper_block);

      std::static_pointer_cast<SimpleMatrix>(_M1)
      ->setBlock(std::max(pos, col), std::min(pos, col),
                 *indexSet.properties(*ei).lower_block);
    }
  }
  else if(_storageType == NM_SPARSE_BLOCK)
  {
    if(! _M2)
    {
      DEBUG_PRINT("Reset _M2 shared pointer using new BlockCSRMatrix(indexSet) \n ");
      _M2.reset(new BlockCSRMatrix(indexSet));

    }
    else
    {
      DEBUG_PRINT("fill existing _M2\n");
      _M2->fill(indexSet);
      DEBUG_EXPR(_M2->display(););
    }
  }
  if(update)
    convert();
  DEBUG_END("void OSNSMatrix::fillM(SP::InteractionsGraph indexSet, bool update)\n");
}

// convert current matrix to NumericsMatrix structure
void OSNSMatrix::convert()
{
  DEBUG_BEGIN("OSNSMatrix::convert()\n");
  DEBUG_PRINTF("_storageType = %i\n", _storageType);


  switch(_storageType)
  {
  case NM_DENSE:
  {
    _numericsMatrix.reset(NM_new(),NM_free_not_dense);
    _numericsMatrix.get()->storageType = _storageType ;
    _numericsMatrix.get()->size0 = _dimRow ;
    _numericsMatrix.get()->size1 = _dimColumn ;
    _numericsMatrix->matrix0 = _M1->getArray(); // Pointer link, be careful when freed.
    DEBUG_EXPR(NM_display(_numericsMatrix.get()););
    DEBUG_EXPR(_M1->display(););
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    _M2->convert();
    _numericsMatrix.reset(NM_new(), NM_free_not_SBM);
    _numericsMatrix.get()->storageType = _storageType ;
    _numericsMatrix.get()->size0 = _dimRow ;
    _numericsMatrix.get()->size1 = _dimColumn ;
    _numericsMatrix->matrix1 = &*_M2->getNumericsMatSparse(); // Pointer link, be careful when freed.
    break;
  }
  case NM_SPARSE:
  {
    // we already filled the matrix
    break;
  }
  default:
  {
    THROW_EXCEPTION("OSNSMatrix::convert unknown _storageType");
  }
  }
  DEBUG_END("OSNSMatrix::convert()\n");
}

// Fill the matrix W
// Used only in GlobalFrictionContact
void OSNSMatrix::fillW(DynamicalSystemsGraph & DSG, bool update)
{
  DEBUG_BEGIN("void OSNSMatrix::fillW(SP::DynamicalSystemsGraph DSG, bool update)\n");

  if(update)
  {
    _dimColumn = updateSizeAndPositions(DSG);
    _dimRow = _dimColumn;
  }

  switch(_storageType)
  {
  case NM_SPARSE:
  {
    if(update)
    {
      size_t sizeM = _dimRow;
      DEBUG_PRINTF("sizeM = %lu \n", sizeM);

      // We choose a triplet matrix format for inserting values.
      // This simplifies the memory manipulation.
      _numericsMatrix.reset(NM_create(NM_SPARSE, sizeM, sizeM),NM_free);

      NumericsMatrix& M_NM = *numericsMatrix();
      NM_triplet_alloc(&M_NM, _triplet_nzmax);
      CSparseMatrix* Mtriplet = NM_triplet(&M_NM);

      unsigned int pos =0;
      // Loop over the DS for filling M
      DynamicalSystemsGraph::VIterator dsi, dsend;
      for(std::tie(dsi, dsend) = DSG.vertices(); dsi != dsend; ++dsi)
      {
        SP::DynamicalSystem ds = DSG.bundle(*dsi);
        SiconosMatrix* W = DSG.properties(*dsi).W.get();
        pos = DSG.properties(*dsi).absolute_position;
        W->fillTriplet(Mtriplet, pos, pos);
        DEBUG_PRINTF("pos = %u \n", pos);
      }
      _triplet_nzmax =  NM_nnz(&M_NM);
    }
    break;
  }
  default:
  {
    THROW_EXCEPTION("OSNSMatrix::fillW unknown _storageType");
  }
  }


  DEBUG_END("void OSNSMatrix::fillW(SP::DynamicalSystemsGraph DSG, bool update)\n");
}


// Fill the matrix Winverse
// Used only in GlobalFrictionContact
void OSNSMatrix::fillWinverse(DynamicalSystemsGraph & DSG, bool update)
{
  DEBUG_BEGIN("void OSNSMatrix::fillWinverse(SP::DynamicalSystemsGraph DSG, bool update)\n");

  if(update)
  {
    _dimColumn = updateSizeAndPositions(DSG);
    _dimRow = _dimColumn;
  }

  switch(_storageType)
  {
  case NM_SPARSE:
  {
    if(update)
    {
      size_t sizeM = _dimRow;
      DEBUG_PRINTF("sizeM = %lu \n", sizeM);

      // We choose a triplet matrix format for inserting values.
      // This simplifies the memory manipulation.
      _numericsMatrix.reset(NM_create(NM_SPARSE, sizeM, sizeM),NM_free);

      NumericsMatrix& M_NM = *numericsMatrix();
      NM_triplet_alloc(&M_NM, _triplet_nzmax);
      CSparseMatrix* Mtriplet = NM_triplet(&M_NM);

      unsigned int pos =0;
      // Loop over the DS for filling M
      DynamicalSystemsGraph::VIterator dsi, dsend;
      for(std::tie(dsi, dsend) = DSG.vertices(); dsi != dsend; ++dsi)
      {
        SP::SimpleMatrix Winverse;

        OneStepIntegrator& osi = *DSG.properties(*dsi).osi;
        SP::DynamicalSystem ds =  DSG.bundle(*dsi);
        SP::SecondOrderDS sods =  std::static_pointer_cast<SecondOrderDS> (ds);

        if (typeid(osi) == typeid(MoreauJeanGOSI))
        {
          Winverse =  static_cast<MoreauJeanGOSI&>(osi).Winverse(sods, true);
        }
        else if (typeid(osi) == typeid(MoreauJeanOSI))
        {
          Winverse =  static_cast<MoreauJeanOSI&>(osi).Winverse(sods);
        }
        else
          THROW_EXCEPTION("OSNSMatrix::fillWinverse not yet implemented for this type of OSI  ");

        pos = DSG.properties(*dsi).absolute_position;
        Winverse->fillTriplet(Mtriplet, pos, pos);
        DEBUG_PRINTF("pos = %u \n", pos);
        //W->display();
      }

      // // Ugly inversion
      // for (unsigned int k =0 ; k < M_NM.size0; k++)
      // {
      //   Mtriplet->x[k] = 1.0/Mtriplet->x[k];
      // }
      //NM_display(numericsMatrix().get());
      _triplet_nzmax =  NM_nnz(&M_NM);
      //getchar();
    }
    break;
  }
  default:
  {
    THROW_EXCEPTION("OSNSMatrix::fillWinverse unknown _storageType");
  }
  }


  DEBUG_END("void OSNSMatrix::fillW(SP::DynamicalSystemsGraph DSG, bool update)\n");
}


#include <float.h>
// Fill the matrix H
void OSNSMatrix::fillH(DynamicalSystemsGraph & DSG, InteractionsGraph& indexSet, bool update)
{
  DEBUG_BEGIN("void OSNSMatrix::fillH(SP::DynamicalSystemsGraph DSG, InteractionsGraph& indexSet, bool update)\n");

  fillHtrans(DSG, indexSet, update);

  SP::NumericsMatrix Htrans = _numericsMatrix;
  _numericsMatrix.reset(NM_transpose(Htrans.get()), NM_free);
  _dimColumn = updateSizeAndPositions(indexSet);
  _dimRow = updateSizeAndPositions(DSG);

  DEBUG_END("void OSNSMatrix::fillH(SP::DynamicalSystemsGraph DSG, InteractionsGraph& indexSet, bool update)\n");
}

// Fill the matrix Htrans
void OSNSMatrix::fillHtrans(DynamicalSystemsGraph & DSG, InteractionsGraph& indexSet, bool update)
{
  DEBUG_BEGIN("void OSNSMatrix::fillHtrans(SP::DynamicalSystemsGraph DSG, InteractionsGraph& indexSet, bool update)\n");
  if(update)
  {

    _dimRow = updateSizeAndPositions(indexSet);
    _dimColumn = updateSizeAndPositions(DSG);
  }

  switch(_storageType)
  {
  case NM_SPARSE:
  {
    if(update)
    {
      // We choose a triplet matrix format for inserting values.
      // This simplifies the memory manipulation.
      _numericsMatrix.reset(NM_create(NM_SPARSE, _dimRow, _dimColumn),NM_free);
      NumericsMatrix& H_NM = *numericsMatrix();
      NM_triplet_alloc(&H_NM, _triplet_nzmax);
      CSparseMatrix* Htriplet= NM_triplet(&H_NM);

      unsigned int pos = 0, abs_pos_ds=0;
      SP::SiconosMatrix leftInteractionBlock;


      InteractionsGraph::VIterator ui, uiend;
      for(std::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui)
      {
        Interaction& inter = *indexSet.bundle(*ui);
        size_t sizeY = inter.dimension();
        leftInteractionBlock = inter.getLeftInteractionBlock();

        double * array = &*leftInteractionBlock->getArray();
        //double * array_with_bc= nullptr;

        SP::DynamicalSystem ds1 = indexSet.properties(*ui).source;
        SP::DynamicalSystem ds2 = indexSet.properties(*ui).target;

        bool endl = false;
        size_t posBlock = indexSet.properties(*ui).source_pos;
        size_t pos_ds2 = indexSet.properties(*ui).target_pos;

        pos =  indexSet.properties(*ui).absolute_position;

        for(SP::DynamicalSystem ds = ds1; !endl; ds = ds2, posBlock = pos_ds2)
        {
          endl = (ds == ds2);
          size_t sizeDS = ds->dimension();

          SecondOrderDS* sods = dynamic_cast<SecondOrderDS*> (ds.get());

          if (sods)
          {
            SP::BoundaryCondition bc;
            if(sods->boundaryConditions())
            {
              // bc = sods->boundaryConditions();
              // NM_dense_display(array,sizeY,sizeDS,sizeY);
              // array_with_bc = (double *) calloc(sizeY*sizeDS,sizeof(double));
              // memcpy(array_with_bc, array ,sizeY*sizeDS,sizeof(double));
              // NM_dense_display(array_with_bc,sizeY,sizeDS,sizeY);
              // for(std::vector<unsigned int>::iterator itindex = bc->velocityIndices()->begin() ;
              //     itindex != bc->velocityIndices()->end();
              //     ++itindex)
              // {


              //   for (unsigned int row; row < sizeY; row++  )
              //   {
              //     array_with_bc[row + (sizeY) * (posBlock + *itindex)] = 0.0
              //   }
              //     // (nslawSize,sizeDS));
              //   //SP::SiconosVector coltmp(new SiconosVector(nslawSize));
              //   //coltmp->zero();
              //   std::cout <<  "bc indx "<< *itindex << std::endl;
              // }


              // //getchar();
              THROW_EXCEPTION("OSNSMatrix::fillHtrans boundary conditions not yet implemented.");
            }
          }


          abs_pos_ds =  DSG.properties(DSG.descriptor(ds)).absolute_position;
          CSparseMatrix_block_dense_zentry(Htriplet,  pos, abs_pos_ds, array+posBlock*sizeY, sizeY, sizeDS, DBL_EPSILON);
        }
      }
      _triplet_nzmax =  NM_nnz(&H_NM);
    }
    break;
  }
  default:
  {
    THROW_EXCEPTION("OSNSMatrix::fillHtrans unknown _storageType");
  }
  }
  DEBUG_END("void OSNSMatrix::fillHtrans(SP::DynamicalSystemsGraph DSG, InteractionsGraph& indexSet, bool update)\n");
}

void OSNSMatrix::computeM(SP::NumericsMatrix Winverse, SP::NumericsMatrix Htrans)
{
   // Compute M = H^T * Winverse * H
  NumericsMatrix *   H_NM = NM_transpose(Htrans.get());


  NumericsMatrix *   NM1 = NM_multiply(Winverse.get(), H_NM);


  _numericsMatrix.reset(NM_multiply(Htrans.get(), NM1), NM_free);


    // NumericsMatrix *   NM1 = NM_multiply(Winverse.get(), H.get());
    // NumericsMatrix *   Htrans_NM = NM_transpose(H.get());

    // _numericsMatrix.reset(NM_multiply(Htrans_NM, NM1), NM_free);



    _dimRow = _numericsMatrix->size0;
    _dimColumn = _numericsMatrix->size1;

    NM_free(NM1);
    NM_free(H_NM);
}



// Display data
void OSNSMatrix::display() const
{
  if(_storageType == NM_DENSE)
  {
    std::cout << "----- OSNS Matrix ( "<< this <<") using default storage type for Numerics structure (SiconosMatrix -> double*)" <<std::endl;
    if(! _M1)
      std::cout << " matrix = nullptr pointer" <<std::endl;
    else _M1->display();
  }
  else if(_storageType == NM_SPARSE_BLOCK)
  {
    std::cout << "----- OSNS Matrix using Sparse InteractionBlock storage type for Numerics (SparseBlockStructuredMatrix)" <<std::endl;
    if(! _M2)
      std::cout << " matrix = nullptr pointer" <<std::endl;
    else _M2->display();
  }
  else if(_storageType == NM_SPARSE)
  {
    std::cout << "----- OSNS Matrix using sparse storage, nothing to show" << std::endl;
  }
}
