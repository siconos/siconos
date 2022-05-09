/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#include "OSNSMatrixProjectOnConstraints.hpp"
#include "Tools.hpp"
#include "RelationTypes.hpp"
#include "NewtonEulerR.hpp"
#include "LagrangianR.hpp"
#include "NonSmoothLaw.hpp"
#include "OSNSMatrix.hpp"
#include "BlockCSRMatrix.hpp"
#include "SimulationGraphs.hpp"
#include "SimpleMatrix.hpp"
using namespace RELATION;
using namespace Siconos;


// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"


OSNSMatrixProjectOnConstraints::OSNSMatrixProjectOnConstraints(unsigned int n, unsigned int m, NM_types stor):
  OSNSMatrix(n, m, stor)
{
}

OSNSMatrixProjectOnConstraints::~OSNSMatrixProjectOnConstraints()
{
}

unsigned OSNSMatrixProjectOnConstraints::updateSizeAndPositions(InteractionsGraph& indexSet)
{
  // === Description ===

  // For a interactionBlock (diagonal or extra diagonal) corresponding to
  // a Interaction, we need to know the position of its first
  // element in the full-matrix M. This position depends on the
  // previous interactionBlocks sizes.
  //
  // positions are saved in a map<SP::Interaction, unsigned int>,
  // named interactionBlocksPositions.
  //

  // Computes real size of the current matrix = sum of the dim. of all
  // Interactionin indexSet
  unsigned dim = 0;
  InteractionsGraph::VIterator vd, vdend;
  DEBUG_EXPR_WE(std::cout << "indexSet :" << &indexSet << std::endl;
             indexSet.display(););
  for(std::tie(vd, vdend) = indexSet.vertices(); vd != vdend; ++vd)
  {
    assert(indexSet.descriptor(indexSet.bundle(*vd)) == *vd);

    //    (*interactionBlocksPositions)[indexSet.bundle(*vd)] = dim;
    DEBUG_EXPR_WE(std::cout << " dim :" << dim << std::endl;
                  std::cout << "vd :" << *vd << std::endl;);


    assert(indexSet.blockProj[*vd]);
    indexSet.properties(*vd).absolute_position_proj = dim;
    SP::Interaction inter = indexSet.bundle(*vd);
    unsigned int nslawSize = computeSizeForProjection(inter);
    dim += nslawSize;
    assert(indexSet.properties(*vd).absolute_position_proj  < dim);
  }

  return dim;
}

void OSNSMatrixProjectOnConstraints::fillM(InteractionsGraph& indexSet, bool update)
{

  if(update)
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

    // ======> Aim: find inter1 and inter2 both in indexSet and which have
    // common DynamicalSystems.  Then get the corresponding matrix
    // from map interactionBlocks, and copy it into M

    unsigned int pos = 0, col = 0; // index position used for
    // interactionBlock copy into M, see
    // below.
    // === Loop through "active" Interactions (ie present in
    // indexSets[level]) ===
    InteractionsGraph::VIterator vi, viend;
    for(std::tie(vi, viend) = indexSet.vertices();
        vi != viend; ++vi)
    {
      SP::Interaction inter = indexSet.bundle(*vi);
      pos = indexSet.properties(*vi).absolute_position_proj;
      assert(indexSet.blockProj[*vi]);
      std::static_pointer_cast<SimpleMatrix>(_M1)
      ->setBlock(pos, pos, *(indexSet.blockProj[*vi]));
    }


    InteractionsGraph::EIterator ei, eiend;
    for(std::tie(ei, eiend) = indexSet.edges();
        ei != eiend; ++ei)
    {
      InteractionsGraph::VDescriptor vd1 = indexSet.source(*ei);
      InteractionsGraph::VDescriptor vd2 = indexSet.target(*ei);

      SP::Interaction inter1 = indexSet.bundle(vd1);
      SP::Interaction inter2 = indexSet.bundle(vd2);

      pos = indexSet.properties(vd1).absolute_position_proj;
      assert(indexSet.is_vertex(inter2));
      col = indexSet.properties(vd2).absolute_position_proj;


      assert(pos < _dimRow);
      assert(col < _dimColumn);


      DEBUG_PRINTF("OSNSMatrix _M1: %i %i\n", _M1->size(0), _M1->size(1));
      DEBUG_PRINTF("OSNSMatrix upper: %i %i\n", (indexSet.upper_blockProj[*ei])->size(0), (indexSet.upper_blockProj[*ei])->size(1));
      DEBUG_PRINTF("OSNSMatrix lower: %i %i\n", (indexSet.lower_blockProj[*ei])->size(0), (indexSet.upper_blockProj[*ei])->size(1));


      std::static_pointer_cast<SimpleMatrix>(_M1)
      ->setBlock(std::min(pos, col), std::max(pos, col),
                 *(indexSet.upper_blockProj[*ei]));

      std::static_pointer_cast<SimpleMatrix>(_M1)
      ->setBlock(std::max(pos, col), std::min(pos, col),
                 *(indexSet.lower_blockProj[*ei]));
    }

  }
  else // if _storageType == NM_SPARSE_BLOCK
  {
    if(! _M2)
      _M2.reset(new BlockCSRMatrix(indexSet));
    else
      _M2->fill(indexSet);
  }
  if(update)
    convert();
}


unsigned int OSNSMatrixProjectOnConstraints::computeSizeForProjection(SP::Interaction inter)
{
  DEBUG_BEGIN( "OSNSMatrixProjectOnConstraints::computeSizeForProjection(SP::Interaction inter)\n");
  RELATION::TYPES relationType;
  relationType = inter->relation()->getType();
  unsigned int nslawSize = inter->nonSmoothLaw()->size();

  unsigned int size =  nslawSize;

  if(Type::value(*(inter->nonSmoothLaw())) == Type::NewtonImpactFrictionNSL ||
      Type::value(*(inter->nonSmoothLaw())) == Type::NewtonImpactNSL)
  {
    if(relationType == NewtonEuler)
    {
      // SP::NewtonEuler1DR ri = std::static_pointer_cast<NewtonEuler1DR> (inter->relation());
      // if(ri->_isOnContact)
      //   equalitySize = 1;
      size = 1;
      DEBUG_EXPR_WE(std::cout << "OSNSMatrixProjectOnConstraints::computeSizeForProjection : NewtonImpact * nslaw and  relationType NewtonEuler. size=1" << std::endl;);
    }
    else if(relationType == Lagrangian)
    {
      size = 1;
      DEBUG_EXPR_WE(std::cout <<
                 "OSNSMatrixProjectOnConstraints::computeSizeForProjection : NewtonImpact * nslaw and relationType Lagrangian. size=1"
                 << std::endl;);
    }
    else
    {
      THROW_EXCEPTION("MLCPProjectOnConstraints::computeSizeForProjection. relation is not of the right type. neither Lagrangian nor NewtonEuler ");
    }
  }
  DEBUG_END( "OSNSMatrixProjectOnConstraints::computeSizeForProjection(SP::Interaction inter)\n");
  return size;

}
