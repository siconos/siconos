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
#include "MLCPProjectOnConstraints.hpp"
#include "MixedComplementarityConditionNSL.hpp"
#include "EqualityConditionNSL.hpp"
#include "Simulation.hpp"
#include "LagrangianR.hpp"
#include "NewtonEulerDS.hpp"
#include "NewtonEulerR.hpp"
#include "NewtonEulerFrom1DLocalFrameR.hpp"
#include "OSNSMatrixProjectOnConstraints.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "Model.hpp"
#include "NonSmoothDynamicalSystem.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"
#include <mlcp_cst.h>

using namespace RELATION;
using namespace Siconos;
//#define MLCPPROJ_DEBUG
//#define MLCPPROJ_WITH_CT
void MLCPProjectOnConstraints::initOSNSMatrix()
{
  _M.reset(new OSNSMatrixProjectOnConstraints(0, 0, _numericsMatrixStorageType));
  _n = 0;
  _m = 0;
  _curBlock = 0;
  _doProjOnEquality = false;
  _useMassNormalization = false;
}


// Constructor from a set of data
MLCPProjectOnConstraints::MLCPProjectOnConstraints(const int numericsSolverId, double alphaval):
  MLCP(numericsSolverId), _alpha(alphaval)
{
  _indexSetLevel = 2;
  _inputOutputLevel = 0;
}


void MLCPProjectOnConstraints::display() const
{
  std::cout << "======= MLCPProjectOnConstraints of size " << _sizeOutput << " with: " <<std::endl;
  std::cout << "======= m " << _m << " _n " << _n <<std::endl;
  LinearOSNS::display();
}
void MLCPProjectOnConstraints::updateInteractionBlocks()
{
  // The present functions checks various conditions and possibly
  // compute interactionBlocks matrices.
  //
  // Let interi and interj be two Interactions.
  //
  // Things to be checked are:
  //  1 - is the topology time invariant?
  //  2 - does interactionBlocks[interi][interj] already exists (ie has been
  //  computed in a previous time step)?
  //  3 - do we need to compute this interactionBlock? A interactionBlock is
  //  to be computed if interi and interj are in IndexSet1 AND if interi and
  //  interj have common DynamicalSystems.
  //
  // The possible cases are:
  //
  //  - If 1 and 2 are true then it does nothing. 3 is not checked.
  //  - If 1 == true, 2 == false, 3 == false, it does nothing.
  //  - If 1 == true, 2 == false, 3 == true, it computes the
  //    interactionBlock.
  //  - If 1==false, 2 is not checked, and the interactionBlock is
  //    computed if 3==true.
  //

  DEBUG_BEGIN(" MLCPProjectOnConstraints::updateInteractionBlocks()\n");

  // Get index set from Simulation
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());

  // It seems that index() in not update in Index(0)
  // see comment in void Simulation::updateIndexSets()
  if (indexSetLevel() == 0)
  {
    indexSet->update_vertices_indices();
    indexSet->update_edges_indices();
  }

  bool isLinear = simulation()->nonSmoothDynamicalSystem()->isLinear();





  // we put diagonal informations on vertices
  // self loops with bgl are a *nightmare* at the moment
  // (patch 65198 on standard boost install)

  if (indexSet->properties().symmetric)
  {
    RuntimeException::selfThrow
      (" MLCPProjectOnConstraints::updateInteractionBlocks() - not yet implemented for symmetric case");
  }
  else // not symmetric => follow out_edges for each vertices
  {
    if (!_hasBeenUpdated)
    {
      //      printf("MLCPProjectOnConstraints::updateInteractionBlocks must be updated.\n");
      _n = 0;
      _m = 0;
      _curBlock = 0;
    }
    InteractionsGraph::VIterator vi, viend;
    for (std11::tie(vi, viend) = indexSet->vertices();
         vi != viend; ++vi)
    {




      SP::Interaction inter = indexSet->bundle(*vi);
      unsigned int nslawSize = std11::static_pointer_cast<OSNSMatrixProjectOnConstraints>
        (_M)->computeSizeForProjection(inter);

      DEBUG_PRINTF("Start to work on Interaction %i of vertex %p\n", inter->number(), *vi);
      if (! indexSet->blockProj[*vi])
      {
        DEBUG_PRINTF("Allocation of blockProj of size %i x %i for interaction %i \n",nslawSize, nslawSize,inter->number());
        indexSet->blockProj[*vi].reset(new SimpleMatrix(nslawSize, nslawSize));
      }

      if (!isLinear || !_hasBeenUpdated)
      {
        computeDiagonalInteractionBlock(*vi);
      }

      /* on a undirected graph, out_edges gives all incident edges */
      InteractionsGraph::OEIterator oei, oeiend;
      /* interactionBlock must be zeroed at init */
      std::map<SP::SiconosMatrix, bool> initialized;
      for (std11::tie(oei, oeiend) = indexSet->out_edges(*vi);
           oei != oeiend; ++oei)
      {
        /* on adjoint graph there is at most 2 edges between source and target */
        InteractionsGraph::EDescriptor ed1, ed2;
        std11::tie(ed1, ed2) = indexSet->edges(indexSet->source(*oei), indexSet->target(*oei));
        if (indexSet->upper_blockProj[ed1])
        {
          initialized[indexSet->upper_blockProj[ed1]] = false;
        }
        // if(indexSet->upper_blockProj[ed2])
        // {
        //   initialized[indexSet->upper_blockProj[ed1]] = false;
        // }

        if (indexSet->lower_blockProj[ed1])
        {
          initialized[indexSet->lower_blockProj[ed2]] = false;
        }
        // if(indexSet->lower_blockProj[ed2])
        // {
        //   initialized[indexSet->lower_blockProj[ed2]] = false;
        // }
      }


      for (std11::tie(oei, oeiend) = indexSet->out_edges(*vi);
           oei != oeiend; ++oei)
      {

        /* on adjoint graph there is at most 2 edges between source and target */
        InteractionsGraph::EDescriptor ed1, ed2;
        std11::tie(ed1, ed2) = indexSet->edges(indexSet->source(*oei), indexSet->target(*oei));

        assert(*oei == ed1 || *oei == ed2);

        /* the first edge as the lower index */
        assert(indexSet->index(ed1) <= indexSet->index(ed2));

        SP::Interaction inter1 = indexSet->bundle(indexSet->source(*oei));
        SP::Interaction inter2 = indexSet->bundle(indexSet->target(*oei));

        // Memory allocation if needed
        unsigned int nslawSize1 = std11::static_pointer_cast<OSNSMatrixProjectOnConstraints>
          (_M)->computeSizeForProjection(inter1);
        unsigned int nslawSize2 = std11::static_pointer_cast<OSNSMatrixProjectOnConstraints>
          (_M)->computeSizeForProjection(inter2);
        unsigned int isrc = indexSet->index(indexSet->source(*oei));
        unsigned int itar = indexSet->index(indexSet->target(*oei));

        SP::SiconosMatrix currentInteractionBlock;

        if (itar > isrc) // upper block
        {
          if (! indexSet->upper_blockProj[ed1])
          {
            indexSet->upper_blockProj[ed1].reset(new SimpleMatrix(nslawSize1, nslawSize2));
            initialized[indexSet->upper_blockProj[ed1]] = false;
            
#ifdef MLCPPROJ_DEBUG
            std::cout <<  "Allocation of upper_blockProj " <<  indexSet->upper_blockProj[ed1].get() << " of edge " << ed1 << " of size " << nslawSize1 << " x " << nslawSize2 << " for interaction " << inter1->number() << " and interaction " <<  inter2->number() <<  std::endl;
#endif

            if (ed2 != ed1)
              indexSet->upper_blockProj[ed2] = indexSet->upper_blockProj[ed1];
          }
#ifdef MLCPPROJ_DEBUG
          else
            std::cout <<  "No Allocation of upper_blockProj of size " << nslawSize1 << " x " << nslawSize2 <<  std::endl;
#endif
          currentInteractionBlock = indexSet->upper_blockProj[ed1];
#ifdef MLCPPROJ_DEBUG
          std::cout << "currentInteractionBlock->size(0)" << currentInteractionBlock->size(0) << std::endl;
          std::cout << "currentInteractionBlock->size(1)" << currentInteractionBlock->size(1) << std::endl;

          std::cout << "inter1->display() " << inter1->number() << std::endl;
          //inter1->display();

          std::cout << "inter2->display() " << inter2->number() << std::endl;
          //inter2->display();
#endif
        }
        else  // lower block
        {
          if (! indexSet->lower_blockProj[ed1])
          {

#ifdef MLCPPROJ_DEBUG
            std::cout <<  "Allocation of lower_blockProj of size " << nslawSize1 << " x " << nslawSize2 << " for interaction " << inter1->number() << " and interaction " <<  inter2->number() <<  std::endl;
#endif
            indexSet->lower_blockProj[ed1].reset(new SimpleMatrix(nslawSize1, nslawSize2));
            initialized[indexSet->lower_blockProj[ed1]] = false;
            if (ed2 != ed1)
              indexSet->lower_blockProj[ed2] = indexSet->lower_blockProj[ed1];
          }
#ifdef MLCPPROJ_DEBUG
          else
            std::cout <<  "No Allocation of lower_blockProj of size " << nslawSize1 << " x " << nslawSize2 <<  std::endl;
#endif
          currentInteractionBlock = indexSet->lower_blockProj[ed1];

#ifdef MLCPPROJ_DEBUG
          std::cout << "currentInteractionBlock->size(0)" << currentInteractionBlock->size(0) << std::endl;
          std::cout << "currentInteractionBlock->size(1)" << currentInteractionBlock->size(1) << std::endl;


          std::cout << "inter1->display() " << inter1->number() << std::endl;
          //inter1->display();

          std::cout << "inter2->display() " << inter2->number() << std::endl;
          //inter2->display();
#endif

        }



        //assert(indexSet->index(ed1));

        if (!initialized[currentInteractionBlock])
        {
          initialized[currentInteractionBlock] = true;
          currentInteractionBlock->zero();
        }


        if (!isLinear || !_hasBeenUpdated)
        {
          if (isrc != itar)
            computeInteractionBlock(*oei);
        }

      }
    }
  }
  DEBUG_EXPR(displayBlocks(indexSet););
  DEBUG_END(" MLCPProjectOnConstraints::updateInteractionBlocks()\n");
}
void MLCPProjectOnConstraints::displayBlocks(SP::InteractionsGraph indexSet)
{

  std::cout <<  "MLCPProjectOnConstraints::displayBlocks(SP::InteractionsGraph indexSet) " << std::endl;
  std::cout << "                          indexSet :" << indexSet << std::endl;


  InteractionsGraph::VIterator vi, viend;
  for (std11::tie(vi, viend) = indexSet->vertices();
       vi != viend; ++vi)
  {
    SP::Interaction inter = indexSet->bundle(*vi);
    std::cout << "                          vertex :" << *vi << std::endl;
    std::cout << "                          bundle :" << indexSet->bundle(*vi) << std::endl;

    if (indexSet->blockProj[*vi])
    {
      std::cout << "                          blockProj ";
      indexSet->blockProj[*vi]->display();
    }

    InteractionsGraph::OEIterator oei, oeiend;



    for (std11::tie(oei, oeiend) = indexSet->out_edges(*vi);
         oei != oeiend; ++oei)
    {
      unsigned int isrc = indexSet->index(indexSet->source(*oei));
      unsigned int itar = indexSet->index(indexSet->target(*oei));
      std::cout << "                          isrc :" << isrc << std::endl;
      std::cout << "                          itar :" << itar << std::endl;


      InteractionsGraph::EDescriptor ed1, ed2;
      std::cout << "                          outedges :" << *oei << std::endl;
      std11::tie(ed1, ed2) = indexSet->edges(indexSet->source(*oei), indexSet->target(*oei));
      std::cout << "                          edges(ed1,ed2) :" << ed1 << " " << ed2  << std::endl;
      std::cout << "                          (ed1)->upper_blockProj : ";
      if (indexSet->upper_blockProj[ed1])
      {
        std::cout << indexSet->upper_blockProj[ed1] << "   :" ;
        indexSet->upper_blockProj[ed1]->display();
      }
      else
        std::cout << "NULL " << std::endl;

      std::cout << "                          (ed1)->lower_blockProj : ";
      if (indexSet->lower_blockProj[ed1])
      {
        std::cout << indexSet->lower_blockProj[ed1] << "   :" ;
        indexSet->lower_blockProj[ed1]->display();
      }
      else
        std::cout << "NULL " << std::endl;

      std::cout << "                          (ed2)->upper_blockProj : ";
      if (indexSet->upper_blockProj[ed2])
      {
        std::cout << indexSet->upper_blockProj[ed2] << "   :" ;
        indexSet->upper_blockProj[ed2]->display();
      }
      else
        std::cout << "NULL" << std::endl;

      std::cout << "                          (ed2)->lower_blockProj : ";
      if (indexSet->lower_blockProj[ed2])
      {
        std::cout << indexSet->lower_blockProj[ed2] << "   :" ;
        indexSet->lower_blockProj[ed2]->display();
      }
      else
        std::cout << "NULL" << std::endl;
    }

  }
}
void MLCPProjectOnConstraints::updateInteractionBlocksOLD()
{

  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());

  bool isLinear = simulation()->nonSmoothDynamicalSystem()->isLinear();
  //  std::cout<<"isLinear: "<<isLinear<<" hasTopologyChanged: "<<hasTopologyChanged<<"hasBeenUpdated: "<<_hasBeenUpdated<<endl;


  if (indexSet->properties().symmetric)
  {
    RuntimeException::selfThrow
      ("MLCPProjectOnConstraints::updateInteractionBlocks() - symmetric case for the indexSet is not yet implemented");
  }
  else // not symmetric => follow out_edges for each vertices
  {
    if (!_hasBeenUpdated || !isLinear)
    {
      if (!_hasBeenUpdated)
      {
        //      printf("MLCPProjectOnConstraints::updateInteractionBlocks must be updated.\n");
        _n = 0;
        _m = 0;
        _curBlock = 0;
      }
      InteractionsGraph::VIterator vi, viend;
      for (std11::tie(vi, viend) = indexSet->vertices();
           vi != viend; ++vi)
      {
        SP::Interaction inter = indexSet->bundle(*vi);
        unsigned int sizeY = 0;
        sizeY = std11::static_pointer_cast<OSNSMatrixProjectOnConstraints>
          (_M)->computeSizeForProjection(inter);

        // #ifdef MLCPPROJ_DEBUG
        //       std::cout<<"\nMLCPProjectOnConstraints::updateInteractionBlocks()"<<endl;
        //        std::cout << "indexSet :"<< indexSet << std::endl;
        //       indexSet->display();
        //        std::cout << "vi :"<< *vi << std::endl;
        //        std::cout << "indexSet->blockProj[*vi]: before"<< indexSet->blockProj[*vi] << std::endl;
        // #endif

        if (! indexSet->blockProj[*vi])
        {
          indexSet->blockProj[*vi].reset(new SimpleMatrix(sizeY, sizeY));
        }
        // #ifdef MLCPPROJ_DEBUG
        //        std::cout << "indexSet->blockProj[*vi]: after"<< indexSet->blockProj[*vi] << std::endl;
        // #endif

        computeDiagonalInteractionBlock(*vi);
      }





      InteractionsGraph::EIterator ei, eiend;
      for (std11::tie(ei, eiend) = indexSet->edges();
           ei != eiend; ++ei)
      {
        SP::Interaction inter1 = indexSet->bundle(indexSet->source(*ei));
        SP::Interaction inter2 = indexSet->bundle(indexSet->target(*ei));
        unsigned int sizeY1 = 0;
        unsigned int sizeY2 = 0;
        sizeY1 = std11::static_pointer_cast<OSNSMatrixProjectOnConstraints>
          (_M)->computeSizeForProjection(inter1);
        sizeY2 = std11::static_pointer_cast<OSNSMatrixProjectOnConstraints>
          (_M)->computeSizeForProjection(inter2);

        // Memory allocation if needed
        unsigned int isrc = indexSet->index(indexSet->source(*ei));
        unsigned int itar = indexSet->index(indexSet->target(*ei));
        if (itar > isrc) // upper block
        {
          if (! indexSet->upper_blockProj[*ei])
          {
            indexSet->upper_blockProj[*ei].reset(new SimpleMatrix(sizeY1, sizeY2));
          }
        }
        else  // lower block
        {
          if (! indexSet->lower_blockProj[*ei])
          {
            indexSet->lower_blockProj[*ei].reset(new SimpleMatrix(sizeY1, sizeY2));
          }
        }

        // Computation of the diagonal block
        computeInteractionBlock(*ei);

        // allocation for transposed block
        // should be avoided
        if (itar > isrc) // upper block has been computed
        {
          // if (!indexSet->lower_blockProj[*ei])
          //   {
          //     indexSet->lower_blockProj[*ei].
          //  reset(new SimpleMatrix(indexSet->upper_blockProj[*ei]->size(1),
          //             indexSet->upper_blockProj[*ei]->size(0)));
          //   }
          indexSet->lower_blockProj[*ei].reset(new SimpleMatrix(*(indexSet->upper_blockProj[*ei])));
          indexSet->lower_blockProj[*ei]->trans();
          //          indexSet->lower_blockProj[*ei]->trans(*indexSet->upper_blockProj[*ei]);
        }
        else
        {
          assert(itar < isrc);    // lower block has been computed
          // if (!indexSet->upper_blockProj[*ei])
          //   {
          //     indexSet->upper_blockProj[*ei].
          //  reset(new SimpleMatrix(indexSet->lower_blockProj[*ei]->size(1),
          //             indexSet->lower_blockProj[*ei]->size(0)));
          //   }
          indexSet->upper_blockProj[*ei].
            reset(new SimpleMatrix(*(indexSet->lower_blockProj[*ei])));
          indexSet->upper_blockProj[*ei]->trans();
        }
        // #ifdef MLCPPROJ_DEBUG
        //             printf("MLCPP upper: %i %i\n",indexSet->upper_blockProj[*ei]->size(0),indexSet->upper_blockProj[*ei]->size(1));
        //             printf("MLCPP lower: %i %i\n",indexSet->lower_blockProj[*ei]->size(0),indexSet->lower_blockProj[*ei]->size(1));
        // #endif

      }

    }

  }

}

void MLCPProjectOnConstraints::computeDiagonalInteractionBlock(const InteractionsGraph::VDescriptor& vd)
{
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());
  SP::Interaction inter = indexSet->bundle(vd);

  // At most 2 DS are linked by an Interaction
  SP::DynamicalSystem ds1;
  SP::DynamicalSystem ds2;
  unsigned int pos1, pos2;
  // --- Get the dynamical system(s) (edge(s)) connected to the current interaction (vertex) ---
  if (indexSet->properties(vd).source != indexSet->properties(vd).target)
  {
    DEBUG_PRINT("a two DS Interaction\n");
    ds1 = indexSet->properties(vd).source;
    ds2 = indexSet->properties(vd).target;
  }
  else
  {
    DEBUG_PRINT("a single DS Interaction\n");
    ds1 = indexSet->properties(vd).source;
    ds2 = ds1;
    // \warning this looks like some debug code, but it gets executed even with NDEBUG.
    // may be compiler does something smarter, but still it should be rewritten. --xhub
    InteractionsGraph::OEIterator oei, oeiend;
    for (std11::tie(oei, oeiend) = indexSet->out_edges(vd);
         oei != oeiend; ++oei)
    {
      // note : at most 4 edges
      ds2 = indexSet->bundle(*oei);
      if (ds2 != ds1)
      {
        assert(false);
        break;
      }
    }
  }
  assert(ds1);
  assert(ds2);
  pos1 = indexSet->properties(vd).source_pos;
  pos2 = indexSet->properties(vd).target_pos;

  // We assume that all ds in vertex_inter have the same osi.
  // SP::OneStepIntegrator Osi = indexSet->properties(vd).osi;
  // //SP::OneStepIntegrator Osi = simulation()->integratorOfDS(ds);
  // OSI::TYPES  osiType = Osi->getType();
  DynamicalSystemsGraph& DSG0 = *simulation()->nonSmoothDynamicalSystem()->dynamicalSystems();

  OneStepIntegrator& osi1 = *DSG0.properties(DSG0.descriptor(ds1)).osi;
  //OneStepIntegrator& osi2 = *DSG0.properties(DSG0.descriptor(ds2)).osi;


  unsigned int sizeY = 0;
  sizeY = std11::static_pointer_cast<OSNSMatrixProjectOnConstraints>
    (_M)->computeSizeForProjection(inter);


#ifdef MLCPPROJ_DEBUG
  std::cout << "\nMLCPProjectOnConstraints::computeDiagonalInteractionBlock" <<std::endl;
  std::cout << "indexSetLevel()" << indexSetLevel() << std::endl;
  //   std::cout << "indexSet :"<< indexSet << std::endl;
  //   std::cout << "vd :"<< vd << std::endl;
  //  indexSet->display();
  //  std::cout << "ds1 :" << std::endl;
  // ds1->display();
  //  std::cout << "ds2 :" << std::endl;
  // ds2->display();
#endif
  assert(indexSet->blockProj[vd]);
  SP::SiconosMatrix currentInteractionBlock = indexSet->blockProj[vd];

#ifdef MLCPPROJ_DEBUG
  //     std::cout<<"MLCPProjectOnConstraints::computeDiagonalInteractionBlock  "<<std::endl;
  //    currentInteractionBlock->display();
  std::cout << "sizeY " << sizeY  << std::endl;
  std::cout <<  "blockProj " <<  indexSet->blockProj[vd].get() << " of edge " << vd << " of size " << currentInteractionBlock->size(0) << " x " << currentInteractionBlock->size(0) << " for interaction " << inter->number() <<  std::endl;
  // std::cout<<"inter1->display() "<< inter1->number()<< std::endl;
  //inter1->display();
  // std::cout<<"inter2->display() "<< inter2->number()<< std::endl;
  //inter2->display();

#endif

  assert(currentInteractionBlock->size(0) == sizeY);
  assert(currentInteractionBlock->size(1) == sizeY);

  if (!_hasBeenUpdated)
    computeOptions(inter, inter);
  // Computes matrix _interactionBlocks[inter1][inter2] (and allocates memory if
  // necessary) if inter1 and inter2 have commond DynamicalSystem.  How
  // _interactionBlocks are computed depends explicitely on the type of
  // Relation of each Interaction.

  // Warning: we suppose that at this point, all non linear
  // operators (G for lagrangian relation for example) have been
  // computed through plug-in mechanism.

  // Get the W and Theta maps of one of the Interaction -
  // Warning: in the current version, if OSI!=MoreauJeanOSI, this fails.
  // If OSI = MOREAU, centralInteractionBlocks = W if OSI = LSODAR,
  // centralInteractionBlocks = M (mass matrices)
  SP::SiconosMatrix leftInteractionBlock, rightInteractionBlock, leftInteractionBlock1;


  // General form of the interactionBlock is : interactionBlock =
  // a*extraInteractionBlock + b * leftInteractionBlock * centralInteractionBlocks
  // * rightInteractionBlock a and b are scalars, centralInteractionBlocks a
  // matrix depending on the integrator (and on the DS), the
  // simulation type ...  left, right and extra depend on the relation
  // type and the non smooth law.


  VectorOfSMatrices& workMInter = *indexSet->properties(vd).workMatrices;

  currentInteractionBlock->zero();

  // loop over the common DS
  bool endl = false;
  unsigned int pos = pos1;
  for (SP::DynamicalSystem ds = ds1; !endl; ds = ds2)
  {
    assert(ds == ds1 || ds == ds2);
    endl = (ds == ds2);

    if (Type::value(*ds) == Type::LagrangianLinearTIDS ||
        Type::value(*ds) == Type::LagrangianDS)
    {
      if (inter->relation()->getType() != Lagrangian)
      {
        RuntimeException::selfThrow(
          "MLCPProjectOnConstraints::computeDiagonalInteractionBlock - relation is not of type Lagrangian with a LagrangianDS.");
      }


      SP::LagrangianDS lds = (std11::static_pointer_cast<LagrangianDS>(ds));
      unsigned int sizeDS = lds->dimension();
      leftInteractionBlock.reset(new SimpleMatrix(sizeY, sizeDS));
      inter->getLeftInteractionBlockForDS(pos, leftInteractionBlock, workMInter);

      if (lds->boundaryConditions()) // V.A. Should we do that ?
      {
        for (std::vector<unsigned int>::iterator itindex =
               lds->boundaryConditions()->velocityIndices()->begin() ;
             itindex != lds->boundaryConditions()->velocityIndices()->end();
             ++itindex)
        {
          // (sizeY,sizeDS));
          SP::SiconosVector coltmp(new SiconosVector(sizeY));
          coltmp->zero();
          leftInteractionBlock->setCol(*itindex, *coltmp);
        }
      }
      // (inter1 == inter2)
      SP::SiconosMatrix work(new SimpleMatrix(*leftInteractionBlock));
      //
      //        std::cout<<"LinearOSNS : leftUBlock\n";
      //        work->display();
      work->trans();
      //        std::cout<<"LinearOSNS::computeInteractionBlock leftInteractionBlock"<<endl;
      //        leftInteractionBlock->display();



      if (_useMassNormalization)
      {
        SP::SiconosMatrix centralInteractionBlock = getOSIMatrix(osi1, ds);
        centralInteractionBlock->PLUForwardBackwardInPlace(*work);
        prod(*leftInteractionBlock, *work, *currentInteractionBlock, false);
        //      gemm(CblasNoTrans,CblasNoTrans,1.0,*leftInteractionBlock,*work,1.0,*currentInteractionBlock);
      }
      else
      {
        prod(*leftInteractionBlock, *work, *currentInteractionBlock, false);
      }


      //*currentInteractionBlock *=h;
    }
    else if (Type::value(*ds) == Type::NewtonEulerDS)
    {

      if (inter->relation()->getType() != NewtonEuler)
      {
        RuntimeException::selfThrow("MLCPProjectOnConstraints::computeDiagonalInteractionBlock - relation is not from NewtonEulerR.");
      }
      SP::NewtonEulerDS neds = (std11::static_pointer_cast<NewtonEulerDS>(ds));
#ifdef MLCPPROJ_WITH_CT
      unsigned int sizeDS = neds->dimension();
      SP::SimpleMatrix T = neds->T();
      SP::SimpleMatrix workT(new SimpleMatrix(*T));
      workT->trans();
      SP::SimpleMatrix workT2(new SimpleMatrix(6, 6));
      prod(*workT, *T, *workT2, true);
      leftInteractionBlock.reset(new SimpleMatrix(sizeY, sizeDS));
      inter->getLeftInteractionBlockForDS(pos, leftInteractionBlock);
      SP::SiconosMatrix work(new SimpleMatrix(*leftInteractionBlock));
      std::cout << "LinearOSNS : leftUBlock\n";
      work->display();
      work->trans();
      std::cout << "LinearOSNS::computeInteractionBlock workT2" <<std::endl;
      workT2->display();
      workT2->PLUForwardBackwardInPlace(*work);
      prod(*leftInteractionBlock, *work, *currentInteractionBlock, false);
#else
      if (0) //(std11::static_pointer_cast<NewtonEulerR> inter->relation())->_isConstact){
      {
//        unsigned int sizeDS = neds->dimension();
//        SP::SimpleMatrix T = neds->T();
//        SP::SimpleMatrix workT(new SimpleMatrix(*T));
//        workT->trans();
//        SP::SimpleMatrix workT2(new SimpleMatrix(6, 6));
//        prod(*workT, *T, *workT2, true);
//        leftInteractionBlock1.reset(new SimpleMatrix(sizeY, sizeDS));
//        inter->getLeftInteractionBlockForDS(pos, leftInteractionBlock);
//        leftInteractionBlock.reset(new SimpleMatrix(1, sizeDS));
//        for (unsigned int ii = 0; ii < sizeDS; ii++)
//          leftInteractionBlock->setValue(1, ii, leftInteractionBlock1->getValue(1, ii));
//
//        SP::SiconosMatrix work(new SimpleMatrix(*leftInteractionBlock));
//        //cout<<"LinearOSNS : leftUBlock\n";
//        //work->display();
//        work->trans();
//        //cout<<"LinearOSNS::computeInteractionBlock workT2"<<endl;
//        //workT2->display();
//        workT2->PLUForwardBackwardInPlace(*work);
//        prod(*leftInteractionBlock, *work, *currentInteractionBlock, false);
      }
      else
      {
        unsigned int sizeDS = (std11::static_pointer_cast<NewtonEulerDS>(ds))->getqDim();
        leftInteractionBlock.reset(new SimpleMatrix(sizeY, sizeDS));
        inter->getLeftInteractionBlockForDSProjectOnConstraints(pos, leftInteractionBlock);
        // #ifdef MLCPPROJ_DEBUG
        //          std::cout << "MLCPProjectOnConstraints::computeDiagonalInteractionBlock - NewtonEuler case leftInteractionBlock : " << std::endl;
        //         leftInteractionBlock->display();
        // #endif

        SP::SiconosMatrix work(new SimpleMatrix(*leftInteractionBlock));
        //cout<<"LinearOSNS sizeY="<<sizeY<<": leftUBlock\n";
        //work->display();
        work->trans();
        prod(*leftInteractionBlock, *work, *currentInteractionBlock, false);
        // #ifdef MLCPPROJ_DEBUG
        //          std::cout << "MLCPProjectOnConstraints::computeDiagonalInteractionBlock - NewtonEuler case currentInteractionBlock : "<< std::endl;
        //         currentInteractionBlock->display();
        // #endif


      }

    }
    else
      RuntimeException::selfThrow("MLCPProjectOnConstraints::computeDiagonalInteractionBlock - ds is not from NewtonEulerDS neither a LagrangianDS.");



#endif
#ifdef MLCPPROJ_DEBUG
    std::cout << "MLCPProjectOnConstraints::computeDiagonalInteractionBlock DiaginteractionBlock " << std::endl;
    currentInteractionBlock->display();
#endif
    // Set pos for next loop.
    pos = pos2;

  }

}

void MLCPProjectOnConstraints::computeInteractionBlock(const InteractionsGraph::EDescriptor& ed)
{

  // Computes matrix _interactionBlocks[inter1][inter2] (and allocates memory if
  // necessary) if inter1 and inter2 have commond DynamicalSystem.  How
  // _interactionBlocks are computed depends explicitely on the type of
  // Relation of each Interaction.

  // Warning: we suppose that at this point, all non linear
  // operators (G for lagrangian relation for example) have been
  // computed through plug-in mechanism.

#ifdef MLCPPROJ_DEBUG
  std::cout << "MLCPProjectOnConstraints::computeInteractionBlock currentInteractionBlock start " << std::endl;
#endif
  // Get dimension of the NonSmoothLaw (ie dim of the interactionBlock)
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());

  SP::DynamicalSystem ds = indexSet->bundle(ed);
  SP::Interaction inter1 = indexSet->bundle(indexSet->source(ed));
  SP::Interaction inter2 = indexSet->bundle(indexSet->target(ed));
  DynamicalSystemsGraph& DSG0 = *simulation()->nonSmoothDynamicalSystem()->dynamicalSystems();
  OneStepIntegrator& Osi = *DSG0.properties(DSG0.descriptor(ds)).osi;
  // For the edge 'ds', we need to find relative position of this ds
  // in inter1 and inter2 relation matrices (--> pos1 and pos2 below)
  // - find if ds is source or target in inter_i
  InteractionsGraph::VDescriptor vertex_inter;
  // - get the corresponding position
  unsigned int pos1, pos2;
  // source of inter1 :
  vertex_inter = indexSet->source(ed);
  VectorOfSMatrices& workMInter1 = *indexSet->properties(vertex_inter).workMatrices;
  SP::DynamicalSystem tmpds = indexSet->properties(vertex_inter).source;
  if (tmpds == ds)
    pos1 =  indexSet->properties(vertex_inter).source_pos;
  else
  {
    tmpds  = indexSet->properties(vertex_inter).target;
    pos1 =  indexSet->properties(vertex_inter).target_pos;
  }
  // now, inter2
  vertex_inter = indexSet->target(ed);
  VectorOfSMatrices& workMInter2 = *indexSet->properties(vertex_inter).workMatrices;
  tmpds = indexSet->properties(vertex_inter).source;
  if (tmpds == ds)
    pos2 =  indexSet->properties(vertex_inter).source_pos;
  else
  {
    tmpds  = indexSet->properties(vertex_inter).target;
    pos2 =  indexSet->properties(vertex_inter).target_pos;
  }

  unsigned int index1 = indexSet->index(indexSet->source(ed));
  unsigned int index2 = indexSet->index(indexSet->target(ed));

  unsigned int sizeY1 = 0;
  sizeY1 = std11::static_pointer_cast<OSNSMatrixProjectOnConstraints>
    (_M)->computeSizeForProjection(inter1);
  unsigned int sizeY2 = 0;
  sizeY2 = std11::static_pointer_cast<OSNSMatrixProjectOnConstraints>
    (_M)->computeSizeForProjection(inter2);

  SP::SiconosMatrix currentInteractionBlock;

  assert(index1 != index2);

  if (index2 > index1) // upper block
  {
    //     if (! indexSet->properties(ed).upper_block)
    //     {
    //       indexSet->properties(ed).upper_block.reset(new SimpleMatrix(sizeY1, sizeY2));
    //     }

    currentInteractionBlock = indexSet->upper_blockProj[ed];
#ifdef MLCPPROJ_DEBUG
    std::cout << "MLCPProjectOnConstraints::computeInteractionBlock currentInteractionBlock " << std::endl;
    //    currentInteractionBlock->display();
    std::cout << "sizeY1 " << sizeY1  << std::endl;
    std::cout << "sizeY2 " << sizeY2  << std::endl;
    std::cout <<  "upper_blockProj " <<  indexSet->upper_blockProj[ed].get() << " of edge " << ed << " of size " << currentInteractionBlock->size(0) << " x " << currentInteractionBlock->size(0) << " for interaction " << inter1->number() << " and interaction " <<  inter2->number() <<  std::endl;
    // std::cout<<"inter1->display() "<< inter1->number()<< std::endl;
    //inter1->display();
    // std::cout<<"inter2->display() "<< inter2->number()<< std::endl;
    //inter2->display();

#endif
    assert(currentInteractionBlock->size(0) == sizeY1);
    assert(currentInteractionBlock->size(1) == sizeY2);

  }
  else  // lower block
  {
    //     if (! indexSet->properties(ed).lower_block)
    //     {
    //       indexSet->properties(ed).lower_block.reset(new SimpleMatrix(sizeY1, sizeY2));
    //     }

    assert(indexSet->lower_blockProj[ed]->size(0) == sizeY1);
    assert(indexSet->lower_blockProj[ed]->size(1) == sizeY2);

    currentInteractionBlock = indexSet->lower_blockProj[ed];
  }


  SP::SiconosMatrix leftInteractionBlock, rightInteractionBlock;

  RELATION::TYPES relationType1, relationType2;

  // General form of the interactionBlock is : interactionBlock =
  // a*extraInteractionBlock + b * leftInteractionBlock * centralInteractionBlocks
  // * rightInteractionBlock a and b are scalars, centralInteractionBlocks a
  // matrix depending on the integrator (and on the DS), the
  // simulation type ...  left, right and extra depend on the relation
  // type and the non smooth law.
  relationType1 = inter1->relation()->getType();
  relationType2 = inter2->relation()->getType();
  if (relationType1 == NewtonEuler &&
      relationType2 == NewtonEuler)
  {
    assert(inter1 != inter2);
    currentInteractionBlock->zero();
#ifdef MLCPPROJ_WITH_CT
    unsigned int sizeDS = (std11::static_pointer_cast<NewtonEulerDS>(ds))->dimension();
    leftInteractionBlock.reset(new SimpleMatrix(sizeY1, sizeDS));
    inter1->getLeftInteractionBlockForDS(pos1, leftInteractionBlock);
    SP::NewtonEulerDS neds = (std11::static_pointer_cast<NewtonEulerDS>(ds));
    SP::SimpleMatrix T = neds->T();
    SP::SimpleMatrix workT(new SimpleMatrix(*T));
    workT->trans();
    SP::SimpleMatrix workT2(new SimpleMatrix(6, 6));
    prod(*workT, *T, *workT2, true);
    rightInteractionBlock.reset(new SimpleMatrix(sizeY2, sizeDS));
    inter2->getLeftInteractionBlockForDS(pos2, rightInteractionBlock);
    rightInteractionBlock->trans();
    workT2->PLUForwardBackwardInPlace(*rightInteractionBlock);
    prod(*leftInteractionBlock, *rightInteractionBlock, *currentInteractionBlock, false);

#else

    unsigned int sizeDS = (std11::static_pointer_cast<NewtonEulerDS>(ds))->getqDim();
    leftInteractionBlock.reset(new SimpleMatrix(sizeY1, sizeDS));
    inter1->getLeftInteractionBlockForDSProjectOnConstraints(pos1, leftInteractionBlock);
    SP::NewtonEulerDS neds = (std11::static_pointer_cast<NewtonEulerDS>(ds));
    rightInteractionBlock.reset(new SimpleMatrix(sizeY2, sizeDS));
    inter2->getLeftInteractionBlockForDSProjectOnConstraints(pos2, rightInteractionBlock);
    rightInteractionBlock->trans();
    prod(*leftInteractionBlock, *rightInteractionBlock, *currentInteractionBlock, false);
  }
#endif
  else if (relationType1 == Lagrangian &&
           relationType2 == Lagrangian)
  {
    unsigned int sizeDS =  ds->dimension();
    leftInteractionBlock.reset(new SimpleMatrix(sizeY1, sizeDS));
    inter1->getLeftInteractionBlockForDS(pos1, leftInteractionBlock, workMInter1);

    Type::Siconos dsType = Type::value(*ds);
    if (dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianDS)
    {
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);

      if (d->boundaryConditions()) // V.A. Should we do that ?
      {
        for (std::vector<unsigned int>::iterator itindex =
               d->boundaryConditions()->velocityIndices()->begin() ;
             itindex != d->boundaryConditions()->velocityIndices()->end();
             ++itindex)
        {
          // (sizeY1,sizeDS));
          SP::SiconosVector coltmp(new SiconosVector(sizeY1));
          coltmp->zero();
          leftInteractionBlock->setCol(*itindex, *coltmp);
        }
      }
    }
#ifdef MLCPPROJ_DEBUG
    std::cout << "MLCPProjectOnConstraints::computeInteractionBlock : leftInteractionBlock" << std::endl;
    leftInteractionBlock->display();
#endif
    // inter1 != inter2
    rightInteractionBlock.reset(new SimpleMatrix(sizeY2, sizeDS));
    inter2->getLeftInteractionBlockForDS(pos2, rightInteractionBlock, workMInter2);
#ifdef MLCPPROJ_DEBUG
    std::cout << "MLCPProjectOnConstraints::computeInteractionBlock : rightInteractionBlock" << std::endl;
    rightInteractionBlock->display();
#endif
    // Warning: we use getLeft for Right interactionBlock
    // because right = transpose(left) and because of
    // size checking inside the getBlock function, a
    // getRight call will fail.
    SP::SiconosMatrix centralInteractionBlock = getOSIMatrix(Osi, ds);
#ifdef MLCPPROJ_DEBUG
    std::cout << "MLCPProjectOnConstraints::computeInteractionBlock : centralInteractionBlocks " << std::endl;
    centralInteractionBlock->display();
#endif
    rightInteractionBlock->trans();

    if (_useMassNormalization)
    {
      centralInteractionBlock->PLUForwardBackwardInPlace(*rightInteractionBlock);
      //*currentInteractionBlock +=  *leftInteractionBlock ** work;
      prod(*leftInteractionBlock, *rightInteractionBlock, *currentInteractionBlock, false);
    }
    else
    {
      prod(*leftInteractionBlock, *rightInteractionBlock, *currentInteractionBlock, false);
    }
#ifdef MLCPPROJ_DEBUG
    std::cout << "MLCPProjectOnConstraints::computeInteractionBlock : currentInteractionBlock" << std::endl;
    currentInteractionBlock->display();
#endif
  }

  else
    RuntimeException::selfThrow("MLCPProjectOnConstraints::computeInteractionBlock not yet implemented for relation of type " + relationType1);

}

void MLCPProjectOnConstraints::computeqBlock(InteractionsGraph::VDescriptor& vertex_inter, unsigned int pos)
{
  DEBUG_BEGIN("MLCPProjectOnConstraints::computeqBlock(InteractionsGraph::VDescriptor& vertex_inter, unsigned int pos)\n");
  
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());
  SP::Interaction inter = indexSet->bundle(vertex_inter);
  unsigned int sizeY = std11::static_pointer_cast<OSNSMatrixProjectOnConstraints>
    (_M)->computeSizeForProjection(inter);
  DEBUG_PRINTF("pos = %i",pos);
  for (unsigned int i = 0; i < sizeY; i++)
    _q->setValue(pos + i, inter->y(0)->getValue(0 + i));

  DEBUG_EXPR(_q->display(););
  DEBUG_END("MLCPProjectOnConstraints::computeqBlock(InteractionsGraph::VDescriptor& vertex_inter, unsigned int pos)\n");
 
}

void MLCPProjectOnConstraints::computeq(double time)
{
  if (_q->size() != _sizeOutput)
    _q->resize(_sizeOutput);
  _q->zero();

  // === Get index set from Simulation ===
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());
  // === Loop through "active" Interactions (ie present in
  // indexSets[level]) ===

  unsigned int pos = 0;
  InteractionsGraph::VIterator ui, uiend;
  for (std11::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    // Compute q, this depends on the type of non smooth problem, on
    // the relation type and on the non smooth law
    pos = indexSet->properties(*ui).absolute_position_proj;
    computeqBlock(*ui, pos); // free output is saved in y
  }
}



void MLCPProjectOnConstraints::postCompute()
{
  _hasBeenUpdated = true;
  // This function is used to set y/lambda values using output from
  // lcp_driver (w,z).  Only Interactions (ie Interactions) of
  // indexSet(leveMin) are concerned.

  // === Get index set from Topology ===
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());

  // y and lambda vectors
  SP::SiconosVector lambda;
  SP::SiconosVector y;

  // === Loop through "active" Interactions (ie present in
  // indexSets[1]) ===
  /** We chose to do a small step _alpha in view of stabilized the algorithm.*/
#ifdef MLCPPROJ_DEBUG
  printf("MLCPProjectOnConstraints::postCompute damping value = %f\n", _alpha);
#endif
  (*_z) *= _alpha;
  unsigned int pos = 0;
#ifdef MLCPPROJ_DEBUG
  printf("MLCPProjectOnConstraints::postCompute _z\n");
  _z->display();
  display();
#endif



  InteractionsGraph::VIterator ui, uiend;

  for (std11::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {

    SP::Interaction inter = indexSet->bundle(*ui);
    // Get the relative position of inter-interactionBlock in the vector w
    // or z
    pos = indexSet->properties(*ui).absolute_position_proj;
    RELATION::TYPES relationType = inter->relation()->getType();
    if (relationType == NewtonEuler)
    {
      postComputeNewtonEulerR(inter, pos);
    }
    else if (relationType == Lagrangian)
    {
      postComputeLagrangianR(inter, pos);
    }
    else
    {
      RuntimeException::selfThrow("MLCPProjectOnConstraints::computeInteractionBlock - relation type is not from Lagrangian type neither NewtonEuler.");
    }

  }



}

void MLCPProjectOnConstraints::postComputeLagrangianR(SP::Interaction inter, unsigned int pos)
{
  SP::LagrangianR  lr = std11::static_pointer_cast<LagrangianR>(inter->relation());
#ifdef MLCPPROJ_DEBUG
  printf("MLCPProjectOnConstraints::postComputeLagrangian inter->y(0)\n");
  inter->y(0)->display();
  printf("MLCPProjectOnConstraints::postComputeLagrangian lr->jachq \n");
  lr->jachq()->display();
  printf("MLCPProjectOnConstraints::postComputeLagrangianR q before update\n");


  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());
  InteractionsGraph::VDescriptor ui = indexSet->descriptor(inter);
  InteractionsGraph::OEIterator oei, oeiend;
    for(std11::tie(oei, oeiend) = indexSet->out_edges(ui);
        oei != oeiend; ++oei)
    {

      SP::LagrangianDS lds =  std11::static_pointer_cast<LagrangianDS>(indexSet->bundle(*oei));
      lds->q()->display();
  }
#endif



  //unsigned int sizeY = inter->nonSmoothLaw()->size();

  // y and lambda vectors
  SP::SiconosVector lambda = inter->lambda(0);
  SP::SiconosVector y = inter->y(0);
  unsigned int sizeY = std11::static_pointer_cast<OSNSMatrixProjectOnConstraints>
    (_M)->computeSizeForProjection(inter);
  // Copy _w/_z values, starting from index pos into y/lambda.

  //setBlock(*_w, y, sizeY, pos, 0);
  setBlock(*_z, lambda, sizeY, pos, 0);

#ifdef MLCPPROJ_DEBUG
  printf("MLCPP lambda of Interaction is pos =%i :\n", pos);
  //  aBuff->display();
  lambda->display();
  unsigned int nslawsize = inter->nonSmoothLaw()->size();
  SP::SiconosVector aBuff(new SiconosVector(nslawsize));
  setBlock(*_z, aBuff, sizeY, pos, 0);
  SP::SiconosMatrix J = lr->jachq();
  SP::SimpleMatrix aux(new SimpleMatrix(*J));
  aux->trans();
  // SP::SiconosVector tmp(new SiconosVector(*(lr->q())));
  // prod(*aux, *aBuff, *(tmp), false);
  // //prod(*aux,*lambda,*(lr->q()),false);
  // std:: std::cout << " tmp =  tmp + J^T * lambda" << std::endl;
  // tmp->display();
#endif



  // // WARNING : Must not be done here. and should be called with the correct time.
  // // compute p(0)
  // inter->computeInput(0.0 ,0);

  // // \warning aBuff should normally be in lambda[0]
  // // The update of the position in DS should be made
  // //  in MoreauJeanOSI::upateState or ProjectedMoreauJeanOSI::updateState
  // SP::SiconosMatrix J=lr->jachq();
  // SP::SimpleMatrix aux(new SimpleMatrix(*J));
  // aux->trans();

  // SP::SiconosVector tmp (new SiconosVector(*(lr->q())));
  // std:: std::cout << " tmp ="<<std::endl;
  // tmp->display();
  // std:: std::cout << " lr->q() ="<<std::endl;
  // lr->q()->display();

  // //prod(*aux,*lambda,*(lr->q()),false);
  // prod(*aux,*aBuff,*(tmp),false);
  // std:: std::cout << " tmp =  tmp + J * lambda"<<std::endl;
  // tmp->display();


  // // The following step should be done on MoreauJeanOSI::upateState or ProjectedMoreauJeanOSI::updateState
  // DSIterator itDS = inter->dynamicalSystemsBegin();
  // while(itDS!=inter->dynamicalSystemsEnd())
  // {
  //   Type::Siconos dsType = Type::value(**itDS);
  //   if((dsType !=Type::LagrangianDS) and
  //      (dsType !=Type::LagrangianLinearTIDS) )
  //   {
  //     RuntimeException::selfThrow("MLCPProjectOnConstraint::postCompute- ds is not of Lagrangian DS type.");
  //   }

  //   SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (*itDS);
  //   SP::SiconosVector q = d->q();

  //   *q +=  *d->p(0);
  //    std::cout << " q=" << std::endl;
  //   q->display();
  //   itDS++;
  // }

  // if ((*lr->q() - *tmp).normInf() > 1e-12)
  // {
  //   RuntimeException::selfThrow("youyou");
  // }

#ifdef MLCPPROJ_DEBUG
  printf("MLCPProjectOnConstraints::postComputeLagrangianR _z\n");
  _z->display();
  printf("MLCPProjectOnConstraints::postComputeLagrangianR updated\n");

  VectorOfBlockVectors& DSlink = *(indexSet->properties(ui)).DSlink;
//  (*DSlink[LagrangianR::q0]).display();
//  (lr->q())->display();
#endif



  //RuntimeException::selfThrow("MLCPProjectOnConstraints::postComputeLagrangianR() - not yet implemented");
}

void MLCPProjectOnConstraints::postComputeNewtonEulerR(SP::Interaction inter, unsigned int pos)
{
  SP::NewtonEulerR ner = (std11::static_pointer_cast<NewtonEulerR>(inter->relation()));
  SP::SiconosVector lambda = inter->lambda(0);
  SP::SiconosVector y = inter->y(0);
  unsigned int sizeY = std11::static_pointer_cast<OSNSMatrixProjectOnConstraints>
    (_M)->computeSizeForProjection(inter);
  // Copy _w/_z values, starting from index pos into y/lambda.

  //setBlock(*_w, y, sizeY, pos, 0);
  setBlock(*_z, lambda, sizeY, pos, 0);

}

void MLCPProjectOnConstraints::computeOptions(SP::Interaction inter1, SP::Interaction inter2)
{
  //  printf("MLCPProjectOnConstraints::computeOptions\n");
  // Get dimension of the NonSmoothLaw (ie dim of the interactionBlock)
  // Retrieve size of Y (projected variable)
  unsigned int sizeY1;
  sizeY1 = std11::static_pointer_cast<OSNSMatrixProjectOnConstraints>
    (_M)->computeSizeForProjection(inter1);

  // Compute the number of equalities
  unsigned int equalitySize1 =  sizeY1; //default behavior

  if (Type::value(*(inter1->nonSmoothLaw())) == Type::NewtonImpactFrictionNSL ||
      Type::value(*(inter1->nonSmoothLaw())) == Type::NewtonImpactNSL)
  {
    if (_doProjOnEquality)
    {
      equalitySize1 = sizeY1;
    }
    else
    {
      equalitySize1 = 0;
    }
  }
  else if (Type::value(*(inter1->nonSmoothLaw()))
           == Type::MixedComplementarityConditionNSL)
  {
    equalitySize1 = std11::static_pointer_cast<MixedComplementarityConditionNSL>(inter1->nonSmoothLaw())->equalitySize();
  }

  // Compute the number of inequalities
  unsigned int inequalitySize1 =  sizeY1 - equalitySize1;



  if (inter1 == inter2)
  {
    //inter1->getExtraInteractionBlock(currentInteractionBlock);
    _m += inequalitySize1;
    _n += equalitySize1;
    //    _m=0;
    //_n=6;
    if (_curBlock > MLCP_NB_BLOCKS - 2)
      printf("MLCP.cpp : number of block to small, memory crach below!!!\n");
    /*add an equality block.*/

    // #ifdef MLCPPROJ_DEBUG
    //   printf("MLCPProjectOnConstraints::computeOptions()\n");
    // #endif

    if (equalitySize1 > 0)
    {
      _numerics_problem.blocksRows[_curBlock + 1] = _numerics_problem.blocksRows[_curBlock] + equalitySize1;
      _numerics_problem.blocksIsComp[_curBlock] = 0;
      // #ifdef MLCPPROJ_DEBUG
      //        std::cout << "_curBlock : " << _curBlock <<std::endl;
      //        std::cout << "_numerics_problem.blocksRows["<<_curBlock+1 <<" ] : " << _numerics_problem.blocksRows[_curBlock+1] <<std::endl;
      //        std::cout << "_numerics_problem.blocksIsComp["<<_curBlock <<" ] : " << _numerics_problem.blocksIsComp[_curBlock] <<std::endl;
      // #endif

      _curBlock++;
    }
    /*add a complementarity block.*/
    if (inequalitySize1 > 0)
    {
      _numerics_problem.blocksRows[_curBlock + 1] = _numerics_problem.blocksRows[_curBlock] + inequalitySize1;
      _numerics_problem.blocksIsComp[_curBlock] = 1;
      // #ifdef MLCPPROJ_DEBUG
      //        std::cout << "_curBlock : " << _curBlock <<std::endl;
      //        std::cout << "_numerics_problem.blocksRows["<<_curBlock+1<< "] : " << _numerics_problem.blocksRows[_curBlock+1] <<std::endl;
      //        std::cout << "_numerics_problem.blocksIsComp["<<_curBlock<< "] : " << _numerics_problem.blocksIsComp[_curBlock] <<std::endl;
      // #endif

      _curBlock++;

    }
  }
  // #ifdef MLCPPROJ_DEBUG
  //    std::cout << "_m : " << _m <<std::endl;
  //    std::cout << "_n : " << _n <<std::endl;
  // #endif
}
