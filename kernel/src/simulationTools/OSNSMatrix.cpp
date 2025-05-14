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
#include "OSNSMatrix.hpp"

#include <iostream>
// #include <assert.h>

#include "BlockCSRMatrix.hpp"
#include "Interaction.hpp"
#include "MoreauJeanGOSI.hpp"
#include "MoreauJeanOSI.hpp"
#include "NonSmoothLaw.hpp"
#include "NumericsToolsNamespace.h"  // For NumericsMatrix
#include "SecondOrderDS.hpp"
#include "SiconosAlgebraAddons.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

siconos::nonsmooth_formulations::OSNSMatrix::OSNSMatrix(unsigned int n, unsigned int m,
                                                        NM_types stor)
    : _dimRow(n), _dimColumn(m), _storageType(stor) {
  // Note:

  // for _storageType = NM_DENSE (dense) n represents the real dimension of
  // the matrix and for sparse storage (_storageType == 1) the number
  // of interactionBlocks in a row or column.
  DEBUG_BEGIN(
      "siconos:simulation::OSNSMatrix::OSNSMatrix(unsigned int n, unsigned int m, "
      "NM_types stor)\n");
  switch (_storageType) {
    case NM_DENSE: {
      // A zero matrix M of size nXn is built.  interactionBlocksPositions
      // remains empty (=nullptr) since we have no information concerning
      // the Interaction.
      _M1 = std::make_shared<siconos::algebra::SiconosMatrix>(n, n);
      break;
    }
    case NM_SPARSE_BLOCK: {
      _M2 = std::make_shared<siconos::simulation::BlockCSRMatrix>(n);
      break;
    }
    default: {
      _triplet_nzmax = _dimRow; /* at least a non zero element per row */
    }
  }

  DEBUG_END(
      "siconos:simulation::OSNSMatrix::OSNSMatrix(unsigned int n, unsigned int m, "
      "NM_types stor)\n");
}

// Build from index set (i.e. get size from number of interactions in the set)
siconos::nonsmooth_formulations::OSNSMatrix::OSNSMatrix(
    siconos::graphs::InteractionsGraph& indexSet, NM_types stor)
    : _dimRow(0), _dimColumn(0), _storageType(stor) {
  DEBUG_BEGIN(
      "siconos:simulation::OSNSMatrix::OSNSMatrix(siconos::graphs::InteractionsGraph& "
      "indexSet, NM_types stor)\n");
  //  _numericsMatrix = std::make_shared<NumericsMatrix>()
  //  NM_null(_numericsMatrix.get());
  fillM(indexSet);
  DEBUG_END(
      "siconos:simulation::OSNSMatrix::OSNSMatrix(siconos::graphs::InteractionsGraph& "
      "indexSet, NM_types stor)\n");
}

// construct by copy of SiconosMatrix
siconos::nonsmooth_formulations::OSNSMatrix::OSNSMatrix(
    const siconos::algebra::SiconosMatrix& MSource)
    : _dimRow(MSource.rows()), _dimColumn(MSource.cols()), _storageType(NM_DENSE) {
  //  _numericsMatrix = std::make_shared<NumericsMatrix>()
  //  NM_null(_numericsMatrix.get());
  _M1 = std::make_shared<siconos::algebra::SiconosMatrix>(MSource);
}

unsigned siconos::nonsmooth_formulations::OSNSMatrix::updateSizeAndPositions(
    siconos::graphs::InteractionsGraph& indexSet) {
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
  siconos::graphs::InteractionsGraph::VIterator vd, vdend;
  for (std::tie(vd, vdend) = indexSet.vertices(); vd != vdend; ++vd) {
    assert(indexSet.descriptor(indexSet.bundle(*vd)) == *vd);
    indexSet.properties(*vd).absolute_position = dim;
    dim += (indexSet.bundle(*vd)->nonSmoothLaw()->size());
    DEBUG_PRINTF("Position = %i for interaction %zu\n", dim, indexSet.bundle(*vd)->number());
    assert(indexSet.properties(*vd).absolute_position < dim);
  }

  return dim;
}

unsigned siconos::nonsmooth_formulations::OSNSMatrix::updateSizeAndPositions(
    siconos::graphs::DynamicalSystemsGraph& DSG) {
  // === Description ===

  // For an interactionBlock (diagonal or extra diagonal) corresponding to
  // an Interaction, we need to know the position of its first
  // element in the full-matrix M. This position depends on the
  // previous interactionBlocks sizes.

  // Computes real size of the current matrix = sum of the dim. of all
  // Interactionin indexSet
  unsigned dim = 0;
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = DSG.vertices(); dsi != dsend; ++dsi) {
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds = DSG.bundle(*dsi);
    DSG.properties(*dsi).absolute_position = dim;
    dim += ds->dimension();
  }

  return dim;
}

// Fill the matrix W
void siconos::nonsmooth_formulations::OSNSMatrix::fillM(
    siconos::graphs::InteractionsGraph& indexSet, bool update) {
  DEBUG_BEGIN("siconos:simulation::OSNSMatrix::fillM()\n");
  DEBUG_PRINTF(" update = %i\n", update);
  if (update)  // If index set vertices list has changed
  {
    // Computes _dimRow and interactionBlocksPositions according to indexSet
    _dimColumn = updateSizeAndPositions(indexSet);
    _dimRow = _dimColumn;
  }

  if (_storageType == NM_DENSE) {
    // === Memory allocation, if required ===
    // Mem. is allocate only if !M or if its size has changed.
    if (update) {
      if (!_M1)
        _M1 = std::make_shared<siconos::algebra::SiconosMatrix>(_dimRow, _dimColumn);
      else {
        if (_M1->rows() != _dimRow || _M1->cols() != _dimColumn)
          _M1->resize(_dimRow, _dimColumn);
        _M1->setZero();
      }
    }

    unsigned int pos = 0, col = 0;  // index position used for

    // === Loop through "active" Interactions (ie vertices present in indexSets[level]) ===
    siconos::graphs::InteractionsGraph::VIterator vi, viend;
    for (std::tie(vi, viend) = indexSet.vertices(); vi != viend; ++vi) {
      std::shared_ptr<siconos::modeling::Interaction> inter = indexSet.bundle(*vi);
      pos = indexSet.properties(*vi).absolute_position;

      siconos::algebra::setBlock(
          pos, pos, *indexSet.properties(*vi).block,
          *std::static_pointer_cast<siconos::algebra::SiconosMatrix>(_M1));
      DEBUG_PRINTF("OSNSMatrix _M1: %i %i\n", _M1->rows(), _M1->cols());
      DEBUG_PRINTF("OSNSMatrix block: %i %i\n", indexSet.properties(*vi).block->rows(),
                   indexSet.properties(*vi).block->cols());
    }

    // == Loop through all edges (ds) in active index set ==
    // Computation of extra-diagonal blocks.
    siconos::graphs::InteractionsGraph::EIterator ei, eiend;
    for (std::tie(ei, eiend) = indexSet.edges(); ei != eiend; ++ei) {
      // For current edge (ds) get source and target vertices (interactions)
      auto vd1 = indexSet.source(*ei);
      auto vd2 = indexSet.target(*ei);

      auto inter1 = indexSet.bundle(vd1);
      auto inter2 = indexSet.bundle(vd2);
      pos = indexSet.properties(vd1).absolute_position;

      assert(indexSet.is_vertex(inter2));

      col = indexSet.properties(vd2).absolute_position;

      assert(pos < _dimRow);
      assert(col < _dimColumn);

      DEBUG_PRINTF("OSNSMatrix _M1: %i %i\n", _M1->rows(), _M1->cols());
      DEBUG_PRINTF("OSNSMatrix upper: %i %i\n", indexSet.properties(*ei).upper_block->rows(),
                   indexSet.properties(*ei).upper_block->cols());
      DEBUG_PRINTF("OSNSMatrix lower: %i %i\n", indexSet.properties(*ei).lower_block->rows(),
                   indexSet.properties(*ei).lower_block->cols());

      assert(indexSet.properties(*ei).lower_block);
      assert(indexSet.properties(*ei).upper_block);
      siconos::algebra::setBlock(
          std::min(pos, col), std::max(pos, col), *indexSet.properties(*ei).upper_block,
          *std::static_pointer_cast<siconos::algebra::SiconosMatrix>(_M1));
      siconos::algebra::setBlock(
          std::max(pos, col), std::min(pos, col), *indexSet.properties(*ei).lower_block,
          *std::static_pointer_cast<siconos::algebra::SiconosMatrix>(_M1));
    }
  } else if (_storageType == NM_SPARSE_BLOCK) {
    if (!_M2) {
      DEBUG_PRINT("Reset _M2 shared pointer with make_shared<BlockCSRMatrix>(indexSet) \n ");
      _M2 = std::make_shared<siconos::simulation::BlockCSRMatrix>(indexSet);
    } else {
      DEBUG_PRINT("fill existing _M2\n");
      _M2->fill(indexSet);
      DEBUG_EXPR(_M2->display(););
    }
  }
  if (update) convert();
  DEBUG_END(
      "void "
      "siconos:simulation::OSNSMatrix::fillM(std::shared_ptr<siconos::graphs::"
      "InteractionsGraph> indexSet, "
      "bool update)\n");
}

// convert current matrix to NumericsMatrix structure
void siconos::nonsmooth_formulations::OSNSMatrix::convert() {
  DEBUG_BEGIN("siconos:simulation::OSNSMatrix::convert()\n");
  DEBUG_PRINTF("_storageType = %i\n", _storageType);

  switch (_storageType) {
    case NM_DENSE: {
      _numericsMatrix.reset(NM_new(), NM_free_not_dense);
      _numericsMatrix.get()->storageType = _storageType;
      _numericsMatrix.get()->size0 = _dimRow;
      _numericsMatrix.get()->size1 = _dimColumn;
      _numericsMatrix->matrix0 = _M1->data();  // Pointer link, be careful when freed.
      DEBUG_EXPR(NM_display(_numericsMatrix.get()););
      DEBUG_EXPR(_M1->display(););
      break;
    }
    case NM_SPARSE_BLOCK: {
      _M2->convert();
      _numericsMatrix.reset(NM_new(), NM_free_not_SBM);
      _numericsMatrix.get()->storageType = _storageType;
      _numericsMatrix.get()->size0 = _dimRow;
      _numericsMatrix.get()->size1 = _dimColumn;
      _numericsMatrix->matrix1 =
          &*_M2->getNumericsMatSparse();  // Pointer link, be careful when freed.
      break;
    }
    case NM_SPARSE: {
      // we already filled the matrix
      break;
    }
    default: {
      THROW_EXCEPTION("siconos:simulation::OSNSMatrix::convert unknown _storageType");
    }
  }
  DEBUG_END("siconos:simulation::OSNSMatrix::convert()\n");
}

// Fill the matrix W
// Used only in GlobalFrictionContact
void siconos::nonsmooth_formulations::OSNSMatrix::fillW(
    siconos::graphs::DynamicalSystemsGraph& DSG, bool update) {
  DEBUG_BEGIN(
      "void siconos:simulation::OSNSMatrix::fillW(std::shared_ptr<DynamicalSystemsGraph> DSG, "
      "bool update)\n");

  if (update) {
    _dimColumn = updateSizeAndPositions(DSG);
    _dimRow = _dimColumn;
  }

  switch (_storageType) {
    case NM_SPARSE: {
      if (update) {
        size_t sizeM = _dimRow;
        DEBUG_PRINTF("sizeM = %lu \n", sizeM);

        // We choose a triplet matrix format for inserting values.
        // This simplifies the memory manipulation.
        _numericsMatrix.reset(NM_create(NM_SPARSE, sizeM, sizeM), NM_free);

        auto& M_NM = *numericsMatrix();
        NM_triplet_alloc(&M_NM, _triplet_nzmax);
        auto Mtriplet = NM_triplet(&M_NM);

        unsigned int pos = 0;
        // Loop over the DS for filling M
        siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
        for (std::tie(dsi, dsend) = DSG.vertices(); dsi != dsend; ++dsi) {
          auto ds = DSG.bundle(*dsi);
          auto W = DSG.properties(*dsi).iterationMatrix.get();
          pos = DSG.properties(*dsi).absolute_position;
          siconos::algebra::fillTriplet(*W, Mtriplet, pos, pos);
          DEBUG_PRINTF("pos = %u \n", pos);
        }
        _triplet_nzmax = NM_nnz(&M_NM);
      }
      break;
    }
    default: {
      THROW_EXCEPTION("siconos:simulation::OSNSMatrix::fillW unknown _storageType");
    }
  }

  DEBUG_END(
      "void "
      "siconos:simulation::OSNSMatrix::fillW(std::shared_ptr<siconos::graphs::"
      "DynamicalSystemsGraph> DSG, "
      "bool update)\n");
}

// Fill the matrix Winverse
// Used only in GlobalFrictionContact
void siconos::nonsmooth_formulations::OSNSMatrix::fillWinverse(
    siconos::graphs::DynamicalSystemsGraph& DSG, bool update) {
  DEBUG_BEGIN(
      "void "
      "siconos:simulation::OSNSMatrix::fillWinverse(std::shared_ptr<siconos::graphs::"
      "DynamicalSystemsGraph> "
      "DSG, bool update)\n");

  if (update) {
    _dimColumn = updateSizeAndPositions(DSG);
    _dimRow = _dimColumn;
  }

  switch (_storageType) {
    case NM_SPARSE: {
      if (update) {
        size_t sizeM = _dimRow;
        DEBUG_PRINTF("sizeM = %lu \n", sizeM);

        // We choose a triplet matrix format for inserting values.
        // This simplifies the memory manipulation.
        _numericsMatrix.reset(NM_create(NM_SPARSE, sizeM, sizeM), NM_free);

        auto& M_NM = *numericsMatrix();
        NM_triplet_alloc(&M_NM, _triplet_nzmax);
        auto Mtriplet = NM_triplet(&M_NM);

        unsigned int pos = 0;
        // Loop over the DS for filling M
        siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
        for (std::tie(dsi, dsend) = DSG.vertices(); dsi != dsend; ++dsi) {
          std::shared_ptr<siconos::algebra::SiconosMatrix> iterationMatrixInverse;

          auto osi = DSG.properties(*dsi).osi;
          std::shared_ptr<siconos::modeling::DynamicalSystem> ds = DSG.bundle(*dsi);
          auto sods = std::static_pointer_cast<siconos::modeling::SecondOrderDS>(ds);

          // if (auto mosi = dynamic_pointer_cast<siconos::integrators::MoreauJeanGOSI>(osi)) {
          //   iterationMatrixInverse = mosi->iterationMatrixInverse(sods, true);
          //} else
          if (auto mosi = dynamic_pointer_cast<siconos::integrators::MoreauJeanOSI>(osi)) {
            auto iterationMatrixInverse =
                mosi->iterationMatrixInverse(sods, *mosi->LUiterationMatrix(ds));

            pos = DSG.properties(*dsi).absolute_position;
            siconos::algebra::fillTriplet(*iterationMatrixInverse, Mtriplet, pos, pos);
            DEBUG_PRINTF("pos = %u \n", pos);

          } else
            THROW_EXCEPTION(
                "siconos:simulation::OSNSMatrix::fillWinverse not yet implemented for this "
                "type of OSI  ");

          // siconos::algebra::print(*W);
        }

        // // Ugly inversion
        // for (unsigned int k =0 ; k < M_NM.size0; k++)
        // {
        //   Mtriplet->x[k] = 1.0/Mtriplet->x[k];
        // }
        // NM_display(numericsMatrix().get());
        _triplet_nzmax = NM_nnz(&M_NM);
        // getchar();
      }
      break;
    }
    default: {
      THROW_EXCEPTION("siconos:simulation::OSNSMatrix::fillWinverse unknown _storageType");
    }
  }

  DEBUG_END(
      "void "
      "siconos:simulation::OSNSMatrix::fillW(std::shared_ptr<siconos::graphs::"
      "DynamicalSystemsGraph> DSG, "
      "bool update)\n");
}

#include <float.h>
// Fill the matrix H
void siconos::nonsmooth_formulations::OSNSMatrix::fillH(
    siconos::graphs::DynamicalSystemsGraph& DSG, siconos::graphs::InteractionsGraph& indexSet,
    bool update) {
  DEBUG_BEGIN(
      "void "
      "siconos:simulation::OSNSMatrix::fillH(std::shared_ptr<siconos::graphs::"
      "DynamicalSystemsGraph> DSG, "
      "siconos::graphs::InteractionsGraph& indexSet, bool update)\n");

  fillHtrans(DSG, indexSet, update);

  auto Htrans = _numericsMatrix;
  _numericsMatrix.reset(NM_transpose(Htrans.get()), NM_free);
  _dimColumn = updateSizeAndPositions(indexSet);
  _dimRow = updateSizeAndPositions(DSG);

  DEBUG_END(
      "void "
      "siconos:simulation::OSNSMatrix::fillH(std::shared_ptr<siconos::graphs::"
      "DynamicalSystemsGraph> DSG, "
      "siconos::graphs::InteractionsGraph& indexSet, bool update)\n");
}

// Fill the matrix Htrans
void siconos::nonsmooth_formulations::OSNSMatrix::fillHtrans(
    siconos::graphs::DynamicalSystemsGraph& DSG, siconos::graphs::InteractionsGraph& indexSet,
    bool update) {
  DEBUG_BEGIN(
      "void "
      "siconos:simulation::OSNSMatrix::fillHtrans(std::shared_ptr<siconos::graphs::"
      "DynamicalSystemsGraph> "
      "DSG, siconos::graphs::InteractionsGraph& indexSet, bool update)\n");
  if (update) {
    _dimRow = updateSizeAndPositions(indexSet);
    _dimColumn = updateSizeAndPositions(DSG);
  }

  switch (_storageType) {
    case NM_SPARSE: {
      if (update) {
        // We choose a triplet matrix format for inserting values.
        // This simplifies the memory manipulation.
        _numericsMatrix.reset(NM_create(NM_SPARSE, _dimRow, _dimColumn), NM_free);
        auto& H_NM = *numericsMatrix();
        NM_triplet_alloc(&H_NM, _triplet_nzmax);
        auto Htriplet = NM_triplet(&H_NM);

        unsigned int pos = 0, abs_pos_ds = 0;

        siconos::graphs::InteractionsGraph::VIterator ui, uiend;
        for (std::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui) {
          auto& inter = *indexSet.bundle(*ui);
          size_t sizeY = inter.dimension();

          auto ds1 = indexSet.properties(*ui).source;
          auto ds2 = indexSet.properties(*ui).target;

          bool endl = false;
          auto posBlock = indexSet.properties(*ui).source_pos;
          auto pos_ds2 = indexSet.properties(*ui).target_pos;

          pos = indexSet.properties(*ui).absolute_position;

          for (auto ds = ds1; !endl; ds = ds2, posBlock = pos_ds2) {
            endl = (ds == ds2);
            size_t sizeDS = ds->dimension();

            auto sods = dynamic_cast<siconos::modeling::SecondOrderDS*>(ds.get());

            if (sods) {
              // std::shared_ptr<siconos::modeling::BoundaryCondition> bc;
              if (sods->boundaryConditions()) {
                // bc = sods->boundaryConditions();
                // NM_dense_display(array,sizeY,sizeDS,sizeY);
                // array_with_bc = (double *) calloc(sizeY*sizeDS,sizeof(double));
                // memcpy(array_with_bc, array ,sizeY*sizeDS,sizeof(double));
                // NM_dense_display(array_with_bc,sizeY,sizeDS,sizeY);
                // for(const auto itindex: bc->velocityIndices())
                // {

                //   for (unsigned int row; row < sizeY; row++  )
                //   {
                //     array_with_bc[row + (sizeY) * (posBlock + itindex)] = 0.0
                //   }
                //     // (nslawSize,sizeDS));
                //   //std::shared_ptr<siconos::algebra::SiconosVector> coltmp =
                //   std::make_shared<siconos::algebra::SiconosVector>(nslawSize));
                //   //coltmp->setZero();
                //   std::cout <<  "bc indx "<< itindex << std::endl;
                // }

                // //getchar();
                THROW_EXCEPTION(
                    "siconos:simulation::OSNSMatrix::fillHtrans boundary conditions not yet "
                    "implemented.");
              }
            }

            abs_pos_ds = DSG.properties(DSG.descriptor(ds)).absolute_position;
            auto leftInteractionBlock = inter.getLeftInteractionBlock();
            const double* raw_array = leftInteractionBlock.data();
            CSparseMatrix_block_dense_zentry(Htriplet, pos, abs_pos_ds,
                                             raw_array + posBlock * sizeY, sizeY, sizeDS,
                                             DBL_EPSILON);
          }
        }
        _triplet_nzmax = NM_nnz(&H_NM);
      }
      break;
    }
    default: {
      THROW_EXCEPTION("siconos:simulation::OSNSMatrix::fillHtrans unknown _storageType");
    }
  }
  DEBUG_END(
      "void "
      "siconos:simulation::OSNSMatrix::fillHtrans(std::shared_ptr<siconos::graphs::"
      "DynamicalSystemsGraph> "
      "DSG, siconos::graphs::InteractionsGraph& indexSet, bool update)\n");
}

void siconos::nonsmooth_formulations::OSNSMatrix::computeM(
    std::shared_ptr<NumericsMatrix> Winverse, std::shared_ptr<NumericsMatrix> Htrans) {
  // Compute M = H^T * Winverse * H
  auto H_NM = NM_transpose(Htrans.get());

  auto NM1 = NM_multiply(Winverse.get(), H_NM);

  _numericsMatrix.reset(NM_multiply(Htrans.get(), NM1), NM_free);

  // auto NM1 = NM_multiply(Winverse.get(), H.get());
  // auto Htrans_NM = NM_transpose(H.get());

  // _numericsMatrix.reset(NM_multiply(Htrans_NM, NM1),
  // NM_free);

  _dimRow = _numericsMatrix->size0;
  _dimColumn = _numericsMatrix->size1;

  NM_free(NM1);
  NM_free(H_NM);
}

// Display data
void siconos::nonsmooth_formulations::OSNSMatrix::display() const {
  if (_storageType == NM_DENSE) {
    std::cout
        << "----- OSNS Matrix ( " << this
        << ") using default storage type for Numerics structure (SiconosMatrix -> double*)"
        << std::endl;
    if (!_M1)
      std::cout << " matrix = nullptr pointer" << std::endl;
    else
      siconos::algebra::print(*_M1);
  } else if (_storageType == NM_SPARSE_BLOCK) {
    std::cout << "----- OSNS Matrix using Sparse InteractionBlock storage type for Numerics "
                 "(SparseBlockStructuredMatrix)"
              << std::endl;
    if (!_M2)
      std::cout << " matrix = nullptr pointer" << std::endl;
    else
      _M2->display();
  } else if (_storageType == NM_SPARSE) {
    std::cout << "----- OSNS Matrix using sparse storage, nothing to show" << std::endl;
  }
}
