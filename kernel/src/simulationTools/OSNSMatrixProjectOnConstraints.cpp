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
//#define OSNSMPROJ_DEBUG


OSNSMatrixProjectOnConstraints::OSNSMatrixProjectOnConstraints(unsigned int n, unsigned int m, int stor):
  OSNSMatrix(n, m, stor)
{
}



// Destructor : pointers are smart
OSNSMatrixProjectOnConstraints::~OSNSMatrixProjectOnConstraints()
{
}

unsigned OSNSMatrixProjectOnConstraints::updateSizeAndPositions(SP::InteractionsGraph indexSet)
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
#ifdef OSNSMPROJ_DEBUG
  std::cout << "indexSet :" << indexSet << std::endl;
  indexSet->display();
#endif
  for (std11::tie(vd, vdend) = indexSet->vertices(); vd != vdend; ++vd)
  {
    assert(indexSet->descriptor(indexSet->bundle(*vd)) == *vd);

    //    (*interactionBlocksPositions)[indexSet->bundle(*vd)] = dim;
#ifdef OSNSMPROJ_DEBUG
    std::cout << " dim :" << dim << std::endl;
    std::cout << "vd :" << *vd << std::endl;

#endif
    assert(indexSet->blockProj[*vd]);

    indexSet->bundle(*vd)->setAbsolutePositionProj(dim);
    SP::Interaction inter = indexSet->bundle(*vd);

    unsigned int nslawSize = computeSizeForProjection(inter);

    dim += nslawSize;
    assert(indexSet->bundle(*vd)->absolutePositionProj() < dim);
  }

  return dim;
}

void OSNSMatrixProjectOnConstraints::fill(SP::InteractionsGraph indexSet, bool update)
{
  assert(indexSet);

  if (update)
  {
    // Computes _dimRow and interactionBlocksPositions according to indexSet
    _dimColumn = updateSizeAndPositions(indexSet);
    _dimRow = _dimColumn;
  }

  if (_storageType == 0)
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
#ifdef OSNSMPROJ_DEBUG
    std::cout << "indexSet :" << indexSet << std::endl;
    indexSet->display();
#endif
    for (std11::tie(vi, viend) = indexSet->vertices();
         vi != viend; ++vi)
    {
      SP::Interaction inter = indexSet->bundle(*vi);
      pos = inter->absolutePositionProj();
      assert(indexSet->blockProj[*vi]);
      std11::static_pointer_cast<SimpleMatrix>(_M1)
      ->setBlock(pos, pos, *(indexSet->blockProj[*vi]));
#ifdef OSNSMPROJ_DEBUG
      printf("OSNSMatrix _M1: %i %i\n", _M1->size(0), _M1->size(1));
      printf("OSNSMatrix block: %i %i\n", indexSet->blockProj[*vi]->size(0), indexSet->blockProj[*vi]->size(1));
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

      pos =  inter1->absolutePositionProj();//(*interactionBlocksPositions)[inter1];

      assert(indexSet->is_vertex(inter2));

      col =  inter2->absolutePositionProj();//(*interactionBlocksPositions)[inter2];


      assert(pos < _dimRow);
      assert(col < _dimColumn);

#ifdef OSNSMPROJ_DEBUG
      printf("OSNSMatrix _M1: %i %i\n", _M1->size(0), _M1->size(1));
      printf("OSNSMatrix upper: %i %i\n", (indexSet->upper_blockProj[*ei])->size(0), (indexSet->upper_blockProj[*ei])->size(1));
      printf("OSNSMatrix lower: %i %i\n", (indexSet->lower_blockProj[*ei])->size(0), (indexSet->upper_blockProj[*ei])->size(1));
#endif

      std11::static_pointer_cast<SimpleMatrix>(_M1)
      ->setBlock(std::min(pos, col), std::max(pos, col),
                 *(indexSet->upper_blockProj[*ei]));

      std11::static_pointer_cast<SimpleMatrix>(_M1)
      ->setBlock(std::max(pos, col), std::min(pos, col),
                 *(indexSet->lower_blockProj[*ei]));
    }

  }
  else // if _storageType == 1
  {
    if (! _M2)
      _M2.reset(new BlockCSRMatrix(indexSet));
    else
      _M2->fill(indexSet);
  }
  if (update)
    convert();
}
unsigned int OSNSMatrixProjectOnConstraints::getPositionOfInteractionBlock(Interaction& inter) const
{
  return inter.absolutePositionProj();
}

unsigned int OSNSMatrixProjectOnConstraints::computeSizeForProjection(SP::Interaction inter)
{
#ifdef OSNSMPROJ_DEBUG
  std::cout << "OSNSMatrixProjectOnConstraints::computeSizeForProjection(SP::Interaction inter)" << std::endl;
#endif




  RELATION::TYPES relationType;
  relationType = inter->relation()->getType();
  unsigned int nslawSize = inter->nonSmoothLaw()->size();

  unsigned int size =  nslawSize;

  if (Type::value(*(inter->nonSmoothLaw())) == Type::NewtonImpactFrictionNSL ||
      Type::value(*(inter->nonSmoothLaw())) == Type::NewtonImpactNSL)
  {
    if (relationType == NewtonEuler)
    {
      // SP::NewtonEulerFrom1DLocalFrameR ri = std11::static_pointer_cast<NewtonEulerFrom1DLocalFrameR> (inter->relation());
      // if(ri->_isOnContact)
      //   equalitySize = 1;
      size = 1;
#ifdef OSNSMPROJ_DEBUG
      std::cout << "OSNSMatrixProjectOnConstraints::computeSizeForProjection : NewtonImpact * nslaw and  relationType NewtonEuler. size=1" << std::endl;
#endif
    }
    else if (relationType == Lagrangian)
    {
      size = 1;
#ifdef OSNSMPROJ_DEBUG
      std::cout << "OSNSMatrixProjectOnConstraints::computeSizeForProjection : NewtonImpact * nslaw and relationType Lagrangian. size=1" << std::endl;
#endif
    }
    else
    {
      RuntimeException::selfThrow("MLCPProjectOnConstraints::computeSizeForProjection. relation is not of the right type. neither Lagrangian nor NewtonEuler ");
    }
  }
#ifdef OSNSMPROJ_DEBUG
  std::cout << "OSNSMatrixProjectOnConstraints::computeSizeForProjection : size= " << size << std::endl;
#endif
  return size;

}
