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
#include "OSNSMatrixProjectOnConstraints.hpp"

#include <assert.h>

#include "BlockCSRMatrix.hpp"
#include "Interaction.hpp"
#include "NonSmoothLaw.hpp"
#include "Relation.hpp"
#include "SiconosMatrix.hpp"
#include "SimulationGraphs.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

siconos::nonsmooth_formulations::OSNSMatrixProjectOnConstraints::
    OSNSMatrixProjectOnConstraints(siconos::algebra::Index n, siconos::algebra::Index m,
                                   NM_types stor)
    : OSNSMatrix(n, m, stor) {}

siconos::algebra::Index
siconos::nonsmooth_formulations::OSNSMatrixProjectOnConstraints::updateSizeAndPositions(
    siconos::graphs::InteractionsGraph& indexSet) {
  // === Description ===

  // For a interactionBlock (diagonal or extra diagonal) corresponding to
  // a Interaction, we need to know the position of its first
  // element in the full-matrix M. This position depends on the
  // previous interactionBlocks sizes.
  //
  // positions are saved in a map<std::shared_ptr<siconos::modeling::Interaction>, unsigned
  // int>, named interactionBlocksPositions.
  //

  // Computes real size of the current matrix = sum of the dim. of all
  // Interactionin indexSet
  siconos::algebra::Index dim = 0;
  siconos::graphs::InteractionsGraph::VIterator vd, vdend;
  DEBUG_EXPR_WE(std::cout << "indexSet :" << &indexSet << std::endl;
                siconos::algebra::print(indexSet););
  for (std::tie(vd, vdend) = indexSet.vertices(); vd != vdend; ++vd) {
    assert(indexSet.descriptor(indexSet.bundle(*vd)) == *vd);

    //    (*interactionBlocksPositions)[indexSet.bundle(*vd)] = dim;
    DEBUG_EXPR_WE(std::cout << " dim :" << dim << std::endl;
                  std::cout << "vd :" << *vd << std::endl;);

    assert(indexSet.blockProj[*vd]);
    indexSet.properties(*vd).absolute_position_proj = dim;
    auto inter = indexSet.bundle(*vd);
    auto nslawSize = computeSizeForProjection(inter);
    dim += nslawSize;
    assert(indexSet.properties(*vd).absolute_position_proj < dim);
  }

  return dim;
}

void siconos::nonsmooth_formulations::OSNSMatrixProjectOnConstraints::fillM(
    siconos::graphs::InteractionsGraph& indexSet, bool update) {
  if (update) {
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

    // ======> Aim: find inter1 and inter2 both in indexSet and which have
    // common DynamicalSystems.  Then get the corresponding matrix
    // from map interactionBlocks, and copy it into M

    siconos::algebra::Index pos = 0, col = 0;  // index position used for
    // interactionBlock copy into M, see
    // below.
    // === Loop through "active" Interactions (ie present in
    // indexSets[level]) ===
    siconos::graphs::InteractionsGraph::VIterator vi, viend;
    for (std::tie(vi, viend) = indexSet.vertices(); vi != viend; ++vi) {
      auto inter = indexSet.bundle(*vi);
      pos = indexSet.properties(*vi).absolute_position_proj;
      assert(indexSet.blockProj[*vi]);
      siconos::algebra::setBlock(
          pos, pos, *(indexSet.blockProj[*vi]),
          *std::static_pointer_cast<siconos::algebra::SiconosMatrix>(_M1));
    }

    siconos::graphs::InteractionsGraph::EIterator ei, eiend;
    for (std::tie(ei, eiend) = indexSet.edges(); ei != eiend; ++ei) {
      auto vd1 = indexSet.source(*ei);
      auto vd2 = indexSet.target(*ei);

      auto inter1 = indexSet.bundle(vd1);
      auto inter2 = indexSet.bundle(vd2);

      pos = indexSet.properties(vd1).absolute_position_proj;
      assert(indexSet.is_vertex(inter2));
      col = indexSet.properties(vd2).absolute_position_proj;

      assert(pos < _dimRow);
      assert(col < _dimColumn);

      DEBUG_PRINTF("OSNSMatrix _M1: %i %i\n", _M1->rows(), _M1->cols());
      DEBUG_PRINTF("OSNSMatrix upper: %i %i\n", (indexSet.upper_blockProj[*ei])->rows(),
                   (indexSet.upper_blockProj[*ei])->cols());
      DEBUG_PRINTF("OSNSMatrix lower: %i %i\n", (indexSet.lower_blockProj[*ei])->rows(),
                   (indexSet.upper_blockProj[*ei])->cols());

      siconos::algebra::setBlock(
          std::min(pos, col), std::max(pos, col), *(indexSet.upper_blockProj[*ei]),
          *std::static_pointer_cast<siconos::algebra::SiconosMatrix>(_M1));

      siconos::algebra::setBlock(
          std::max(pos, col), std::min(pos, col), *(indexSet.lower_blockProj[*ei]),
          *std::static_pointer_cast<siconos::algebra::SiconosMatrix>(_M1));
    }
  } else  // if _storageType == NM_SPARSE_BLOCK
  {
    if (!_M2)
      _M2 = std::make_shared<siconos::simulation::BlockCSRMatrix>(indexSet);
    else
      _M2->fill(indexSet);
  }
  if (update) convert();
}

siconos::algebra::Index
siconos::nonsmooth_formulations::OSNSMatrixProjectOnConstraints::computeSizeForProjection(
    std::shared_ptr<siconos::modeling::Interaction> inter) {
  DEBUG_BEGIN(
      "siconos::nonsmooth_formulations::OSNSMatrixProjectOnConstraints::"
      "computeSizeForProjection(std::"
      "shared_ptr<siconos::"
      "modeling::Interaction> inter)\n");
  auto relationType = inter->relation()->getType();
  auto nslawSize = inter->nonSmoothLaw()->size();

  auto size = nslawSize;

  if (siconos::types::type_value(*(inter->nonSmoothLaw())) ==
          siconos::modeling::Type::NewtonImpactFrictionNSL ||
      siconos::types::type_value(*(inter->nonSmoothLaw())) ==
          siconos::modeling::Type::NewtonImpactNSL) {
    if (relationType == siconos::modeling::RelationType::NewtonEuler) {
      // auto ri = std::static_pointer_cast<NewtonEuler1DR> (inter->relation());
      // if(ri->_isOnContact)
      //   equalitySize = 1;
      size = 1;
      DEBUG_EXPR_WE(
          std::cout << "siconos::nonsmooth_formulations::OSNSMatrixProjectOnConstraints::"
                       "computeSizeForProjection : "
                       "NewtonImpact * nslaw and  relationType NewtonEuler. size=1"
                    << std::endl;);
    } else if (relationType == siconos::modeling::RelationType::Lagrangian) {
      size = 1;
      DEBUG_EXPR_WE(
          std::cout << "siconos::nonsmooth_formulations::OSNSMatrixProjectOnConstraints::"
                       "computeSizeForProjection : "
                       "NewtonImpact * nslaw and relationType Lagrangian. size=1"
                    << std::endl;);
    } else {
      THROW_EXCEPTION(
          "MLCPProjectOnConstraints::computeSizeForProjection. relation is not of the right "
          "type. neither Lagrangian nor NewtonEuler ");
    }
  }
  DEBUG_END(
      "siconos::nonsmooth_formulations::OSNSMatrixProjectOnConstraints::"
      "computeSizeForProjection(std::"
      "shared_ptr<siconos::"
      "modeling::Interaction> inter)\n");
  return size;
}
