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
#include "MLCPProjectOnConstraints.hpp"

#include "BoundaryCondition.hpp"
#include "Interaction.hpp"
#include "LagrangianDS.hpp"
#include "LagrangianR.hpp"
#include "MixedComplementarityConditionNSL.hpp"
#include "NewtonEulerDS.hpp"
#include "NewtonEulerR.hpp"
#include "NumericsSolversNamespace.h"  // solver_options stuff
#include "OSNSMatrixProjectOnConstraints.hpp"
#include "Relation.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixOp.hpp"  // mat prod
#include "SiconosVector.hpp"
#include "SiconosVectorOp.hpp"  // for setBlock
#include "Simulation.hpp"
#include "Tools.hpp"  // enum_to_string
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES

#include "siconos_debug.h"

// #define MLCPPROJ_DEBUG
// #define MLCPPROJ_WITH_CT
void siconos::nonsmooth_formulations::MLCPProjectOnConstraints::initOSNSMatrix() {
  _M = std::make_shared<OSNSMatrixProjectOnConstraints>(0, 0, _numericsMatrixStorageType);
  _n = 0;
  _m = 0;
  _curBlock = 0;
  _doProjOnEquality = false;
  _useMassNormalization = false;
}

// Constructor from a set of data
siconos::nonsmooth_formulations::MLCPProjectOnConstraints::MLCPProjectOnConstraints(
    std::shared_ptr<SolverOptions> options, double alphaval)
    : MLCP(options), _alpha(alphaval) {
  _indexSetLevel = 2;
  _inputOutputLevel = 0;
}

// Constructor from a set of data
siconos::nonsmooth_formulations::MLCPProjectOnConstraints::MLCPProjectOnConstraints(
    const int numericsSolverId, double alphaval)
    : MLCP(numericsSolverId), _alpha(alphaval) {
  _indexSetLevel = 2;
  _inputOutputLevel = 0;
}

void siconos::nonsmooth_formulations::MLCPProjectOnConstraints::display() const {
  std::cout << "======= MLCPProjectOnConstraints of size " << _sizeOutput << " with: \n";
  std::cout << "======= m " << _m << " _n " << _n << "\n";
  LinearOSNS::display();
}
void siconos::nonsmooth_formulations::MLCPProjectOnConstraints::updateInteractionBlocks() {
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

  DEBUG_BEGIN(
      " siconos::nonsmooth_formulations::MLCPProjectOnConstraints::updateInteractionBlocks()"
      "\n");

  // Get index set from Simulation
  auto indexSet = simulation()->indexSet(indexSetLevel());

  // It seems that index() in not update in Index(0)
  // see comment in void Simulation::updateIndexSets()
  if (indexSetLevel() == 0) {
    indexSet->update_vertices_indices();
    indexSet->update_edges_indices();
  }

  bool isLinear = simulation()->nonSmoothDynamicalSystem()->isLinear();

  // we put diagonal information on vertices
  // self loops with bgl are a *nightmare* at the moment
  // (patch 65198 on standard boost install)

  if (indexSet->properties().symmetric) {
    THROW_EXCEPTION(
        " siconos::nonsmooth_formulations::MLCPProjectOnConstraints::updateInteractionBlocks()"
        " - not yet "
        "implemented for "
        "symmetric case");
  } else  // not symmetric => follow out_edges for each vertices
  {
    if (!_hasBeenUpdated) {
      //      printf("siconos::nonsmooth_formulations::MLCPProjectOnConstraints::updateInteractionBlocks
      //      must be updated.\n");
      _n = 0;
      _m = 0;
      _curBlock = 0;
    }
    siconos::graphs::InteractionsGraph::VIterator vi, viend;
    for (std::tie(vi, viend) = indexSet->vertices(); vi != viend; ++vi) {
      auto inter = indexSet->bundle(*vi);
      auto nslawSize = std::static_pointer_cast<OSNSMatrixProjectOnConstraints>(_M)
                           ->computeSizeForProjection(inter);

      DEBUG_PRINTF("Start to work on Interaction %i of vertex %p\n", inter->number(), *vi);
      if (!indexSet->blockProj[*vi]) {
        DEBUG_PRINTF("Allocation of blockProj of size %i x %i for interaction %i \n",
                     nslawSize, nslawSize, inter->number());
        indexSet->blockProj[*vi] =
            std::make_shared<siconos::algebra::SiconosMatrix>(nslawSize, nslawSize);
      }

      if (!isLinear || !_hasBeenUpdated) {
        computeDiagonalInteractionBlock(*vi);
      }

      /* on a undirected graph, out_edges gives all incident edges */
      siconos::graphs::InteractionsGraph::OEIterator oei, oeiend;
      /* interactionBlock must be zeroed at init */
      std::map<std::shared_ptr<siconos::algebra::SiconosMatrix>, bool> initialized;
      for (std::tie(oei, oeiend) = indexSet->out_edges(*vi); oei != oeiend; ++oei) {
        /* on adjoint graph there is at most 2 edges between source and target */
        siconos::graphs::InteractionsGraph::EDescriptor ed1, ed2;
        std::tie(ed1, ed2) = indexSet->edges(indexSet->source(*oei), indexSet->target(*oei));
        if (indexSet->upper_blockProj[ed1]) {
          initialized[indexSet->upper_blockProj[ed1]] = false;
        }
        // if(indexSet->upper_blockProj[ed2])
        // {
        //   initialized[indexSet->upper_blockProj[ed1]] = false;
        // }

        if (indexSet->lower_blockProj[ed1]) {
          initialized[indexSet->lower_blockProj[ed2]] = false;
        }
        // if(indexSet->lower_blockProj[ed2])
        // {
        //   initialized[indexSet->lower_blockProj[ed2]] = false;
        // }
      }

      for (std::tie(oei, oeiend) = indexSet->out_edges(*vi); oei != oeiend; ++oei) {
        /* on adjoint graph there is at most 2 edges between source and target */
        siconos::graphs::InteractionsGraph::EDescriptor ed1, ed2;
        std::tie(ed1, ed2) = indexSet->edges(indexSet->source(*oei), indexSet->target(*oei));

        assert(*oei == ed1 || *oei == ed2);

        /* the first edge as the lower index */
        assert(indexSet->index(ed1) <= indexSet->index(ed2));

        auto inter1 = indexSet->bundle(indexSet->source(*oei));
        auto inter2 = indexSet->bundle(indexSet->target(*oei));

        // Memory allocation if needed
        auto nslawSize1 = std::static_pointer_cast<OSNSMatrixProjectOnConstraints>(_M)
                              ->computeSizeForProjection(inter1);
        auto nslawSize2 = std::static_pointer_cast<OSNSMatrixProjectOnConstraints>(_M)
                              ->computeSizeForProjection(inter2);
        auto isrc = indexSet->index(indexSet->source(*oei));
        auto itar = indexSet->index(indexSet->target(*oei));

        std::shared_ptr<siconos::algebra::SiconosMatrix> currentInteractionBlock;

        if (itar > isrc)  // upper block
        {
          if (!indexSet->upper_blockProj[ed1]) {
            indexSet->upper_blockProj[ed1] =
                std::make_shared<siconos::algebra::SiconosMatrix>(nslawSize1, nslawSize2);
            initialized[indexSet->upper_blockProj[ed1]] = false;

#ifdef MLCPPROJ_DEBUG
            std::cout << "Allocation of upper_blockProj "
                      << indexSet->upper_blockProj[ed1].get() << " of edge " << ed1
                      << " of size " << nslawSize1 << " x " << nslawSize2
                      << " for interaction " << inter1->number() << " and interaction "
                      << inter2->number() << "\n";
#endif

            if (ed2 != ed1) indexSet->upper_blockProj[ed2] = indexSet->upper_blockProj[ed1];
          }
#ifdef MLCPPROJ_DEBUG
          else
            std::cout << "No Allocation of upper_blockProj of size " << nslawSize1 << " x "
                      << nslawSize2 << "\n";
#endif
          currentInteractionBlock = indexSet->upper_blockProj[ed1];
#ifdef MLCPPROJ_DEBUG
          std::cout << "currentInteractionBlock->size(0)" << currentInteractionBlock->size(0)
                    << "\n";
          std::cout << "currentInteractionBlock->size(1)" << currentInteractionBlock->size(1)
                    << "\n";

          std::cout << "inter1->display() " << inter1->number() << "\n";
          // inter1->display();

          std::cout << "inter2->display() " << inter2->number() << "\n";
          // inter2->display();
#endif
        } else  // lower block
        {
          if (!indexSet->lower_blockProj[ed1]) {
#ifdef MLCPPROJ_DEBUG
            std::cout << "Allocation of lower_blockProj of size " << nslawSize1 << " x "
                      << nslawSize2 << " for interaction " << inter1->number()
                      << " and interaction " << inter2->number() << "\n";
#endif
            indexSet->lower_blockProj[ed1] =
                std::make_shared<siconos::algebra::SiconosMatrix>(nslawSize1, nslawSize2);
            initialized[indexSet->lower_blockProj[ed1]] = false;
            if (ed2 != ed1) indexSet->lower_blockProj[ed2] = indexSet->lower_blockProj[ed1];
          }
#ifdef MLCPPROJ_DEBUG
          else
            std::cout << "No Allocation of lower_blockProj of size " << nslawSize1 << " x "
                      << nslawSize2 << "\n";
#endif
          currentInteractionBlock = indexSet->lower_blockProj[ed1];

#ifdef MLCPPROJ_DEBUG
          std::cout << "currentInteractionBlock->size(0)" << currentInteractionBlock->size(0)
                    << "\n";
          std::cout << "currentInteractionBlock->size(1)" << currentInteractionBlock->size(1)
                    << "\n";

          std::cout << "inter1->display() " << inter1->number() << "\n";
          // inter1->display();

          std::cout << "inter2->display() " << inter2->number() << "\n";
          // inter2->display();
#endif
        }

        // assert(indexSet->index(ed1));

        if (!initialized[currentInteractionBlock]) {
          initialized[currentInteractionBlock] = true;
          currentInteractionBlock->setZero();
        }

        if (!isLinear || !_hasBeenUpdated) {
          if (isrc != itar) computeInteractionBlock(*oei);
        }
      }
    }
  }
  DEBUG_EXPR(displayBlocks(indexSet););
  DEBUG_END(
      " siconos::nonsmooth_formulations::MLCPProjectOnConstraints::updateInteractionBlocks()"
      "\n");
}
void siconos::nonsmooth_formulations::MLCPProjectOnConstraints::displayBlocks(
    std::shared_ptr<siconos::graphs::InteractionsGraph> indexSet) {
  std::cout << "siconos::nonsmooth_formulations::MLCPProjectOnConstraints::displayBlocks(std::"
               "shared_ptr<"
               "siconos::graphs::"
               "InteractionsGraph> indexSet) "
            << "\n";
  std::cout << "                          indexSet :" << indexSet << "\n";

  siconos::graphs::InteractionsGraph::VIterator vi, viend;
  for (std::tie(vi, viend) = indexSet->vertices(); vi != viend; ++vi) {
    auto inter = indexSet->bundle(*vi);
    std::cout << "                          vertex :" << *vi << "\n";
    std::cout << "                          bundle :" << indexSet->bundle(*vi) << "\n";

    if (indexSet->blockProj[*vi]) {
      std::cout << "                          blockProj ";
      indexSet->blockProj[*vi]->display();
    }

    siconos::graphs::InteractionsGraph::OEIterator oei, oeiend;

    for (std::tie(oei, oeiend) = indexSet->out_edges(*vi); oei != oeiend; ++oei) {
      auto isrc = indexSet->index(indexSet->source(*oei));
      auto itar = indexSet->index(indexSet->target(*oei));
      std::cout << "                          isrc :" << isrc << "\n";
      std::cout << "                          itar :" << itar << "\n";

      siconos::graphs::InteractionsGraph::EDescriptor ed1, ed2;
      std::cout << "                          outedges :" << *oei << "\n";
      std::tie(ed1, ed2) = indexSet->edges(indexSet->source(*oei), indexSet->target(*oei));
      std::cout << "                          edges(ed1,ed2) :" << ed1 << " " << ed2 << "\n";
      std::cout << "                          (ed1)->upper_blockProj : ";
      if (indexSet->upper_blockProj[ed1]) {
        std::cout << indexSet->upper_blockProj[ed1] << "   :";
        indexSet->upper_blockProj[ed1]->display();
      } else
        std::cout << "nullptr \n";

      std::cout << "                          (ed1)->lower_blockProj : ";
      if (indexSet->lower_blockProj[ed1]) {
        std::cout << indexSet->lower_blockProj[ed1] << "   :";
        indexSet->lower_blockProj[ed1]->display();
      } else
        std::cout << "nullptr \n";

      std::cout << "                          (ed2)->upper_blockProj : ";
      if (indexSet->upper_blockProj[ed2]) {
        std::cout << indexSet->upper_blockProj[ed2] << "   :";
        indexSet->upper_blockProj[ed2]->display();
      } else
        std::cout << "nullptr\n";

      std::cout << "                          (ed2)->lower_blockProj : ";
      if (indexSet->lower_blockProj[ed2]) {
        std::cout << indexSet->lower_blockProj[ed2] << "   :";
        indexSet->lower_blockProj[ed2]->display();
      } else
        std::cout << "nullptr\n";
    }
  }
}

void siconos::nonsmooth_formulations::MLCPProjectOnConstraints::
    computeDiagonalInteractionBlock(
        const siconos::graphs::InteractionsGraph::VDescriptor& vd) {
  auto indexSet = simulation()->indexSet(indexSetLevel());
  auto inter = indexSet->bundle(vd);

  // At most 2 DS are linked by an Interaction
  std::shared_ptr<siconos::modeling::DynamicalSystem> ds1;
  std::shared_ptr<siconos::modeling::DynamicalSystem> ds2;
  unsigned int pos1, pos2;
  // --- Get the dynamical system(s) (edge(s)) connected to the current interaction (vertex)
  // ---
  if (indexSet->properties(vd).source != indexSet->properties(vd).target) {
    DEBUG_PRINT("a two DS Interaction\n");
    ds1 = indexSet->properties(vd).source;
    ds2 = indexSet->properties(vd).target;
  } else {
    DEBUG_PRINT("a single DS Interaction\n");
    ds1 = indexSet->properties(vd).source;
    ds2 = ds1;
    // \warning this looks like some debug code, but it gets executed even with NDEBUG.
    // may be compiler does something smarter, but still it should be rewritten. --xhub
    siconos::graphs::InteractionsGraph::OEIterator oei, oeiend;
    for (std::tie(oei, oeiend) = indexSet->out_edges(vd); oei != oeiend; ++oei) {
      // note : at most 4 edges
      ds2 = indexSet->bundle(*oei);
      if (ds2 != ds1) {
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
  auto& DSG0 = *simulation()->nonSmoothDynamicalSystem()->dynamicalSystems();

  auto& osi1 = *DSG0.properties(DSG0.descriptor(ds1)).osi;
  // OneStepIntegrator& osi2 = *DSG0.properties(DSG0.descriptor(ds2)).
  auto sizeY =
      std::static_pointer_cast<OSNSMatrixProjectOnConstraints>(_M)->computeSizeForProjection(
          inter);

#ifdef MLCPPROJ_DEBUG
  std::cout << "\nsiconos::nonsmooth_formulations::MLCPProjectOnConstraints::"
               "computeDiagonalInteractionBlock"
            << "\n";
  std::cout << "indexSetLevel()" << indexSetLevel() << "\n";
  //   std::cout << "indexSet :"<< indexSet << "\n";
  //   std::cout << "vd :"<< vd << "\n";
  //  indexSet->display();
  //  std::cout << "ds1 :\n";
  // ds1->display();
  //  std::cout << "ds2 :\n";
  // ds2->display();
#endif
  assert(indexSet->blockProj[vd]);
  auto currentInteractionBlock = indexSet->blockProj[vd];

#ifdef MLCPPROJ_DEBUG
  //     std::cout<<"siconos::nonsmooth_formulations::MLCPProjectOnConstraints::computeDiagonalInteractionBlock
  //     "<<std::endl;
  //    currentInteractionBlock->display();
  std::cout << "sizeY " << sizeY << "\n";
  std::cout << "blockProj " << indexSet->blockProj[vd].get() << " of edge " << vd
            << " of size " << currentInteractionBlock->size(0) << " x "
            << currentInteractionBlock->size(0) << " for interaction " << inter->number()
            << "\n";
  // std::cout<<"inter1->display() "<< inter1->number()<< "\n";
  // inter1->display();
  // std::cout<<"inter2->display() "<< inter2->number()<< "\n";
  // inter2->display();

#endif

  assert(currentInteractionBlock->size(0) == sizeY);
  assert(currentInteractionBlock->size(1) == sizeY);

  if (!_hasBeenUpdated) computeOptions(inter, inter);
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

  // General form of the interactionBlock is : interactionBlock =
  // a*extraInteractionBlock + b * leftInteractionBlock * centralInteractionBlocks
  // * rightInteractionBlock a and b are scalars, centralInteractionBlocks a
  // matrix depending on the integrator (and on the DS), the
  // simulation type ...  left, right and extra depend on the relation
  // type and the non smooth law.

  currentInteractionBlock->setZero();

  // loop over the common DS
  bool endl = false;
  unsigned int pos = pos1;
  for (auto ds = ds1; !endl; ds = ds2) {
    assert(ds == ds1 || ds == ds2);
    endl = (ds == ds2);

    if (auto lds = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      if (inter->relation()->getType() != siconos::modeling::RelationType::Lagrangian) {
        THROW_EXCEPTION(
            "siconos::nonsmooth_formulations::MLCPProjectOnConstraints::"
            "computeDiagonalInteractionBlock - "
            "relation is not of type Lagrangian with a LagrangianDS.");
      }

      auto sizeDS = lds->dimension();
      auto leftInteractionBlock = inter->getLeftInteractionBlockForDS(pos, sizeY, sizeDS);

      if (lds->boundaryConditions())  // V.A. Should we do that ?
      {
        for (const auto itindex : lds->boundaryConditions()->velocityIndices()) {
          auto coltmp = std::make_shared<siconos::algebra::SiconosVector>(sizeY);
          coltmp->setZero();
          leftInteractionBlock->col(itindex) = *coltmp;
        }
      }
      // (inter1 == inter2)
      auto work = std::make_shared<siconos::algebra::SiconosMatrix>(*leftInteractionBlock);
      //
      //        std::cout<<"LinearOSNS : leftUBlock\n";
      //        work->display();
      work->transposeInPlace();
      //        std::cout<<"LinearOSNS::computeInteractionBlock leftInteractionBlock"<<endl;
      //        leftInteractionBlock->display();

      if (_useMassNormalization) {
        auto centralInteractionBlock = getOSIMatrix(osi1, ds);
        *work = centralInteractionBlock->solve(*work);
        *currentInteractionBlock += *leftInteractionBlock * *work;
        //      gemm(CblasNoTrans,CblasNoTrans,1.0,*leftInteractionBlock,*work,1.0,*currentInteractionBlock);
      } else {
        siconos::algebra::prod(*leftInteractionBlock, *work, *currentInteractionBlock, false);
      }

      //*currentInteractionBlock *=h;
    }

    else if (auto neds = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
      if (inter->relation()->getType() != siconos::modeling::RelationType::NewtonEuler) {
        THROW_EXCEPTION(
            "siconos::nonsmooth_formulations::MLCPProjectOnConstraints::"
            "computeDiagonalInteractionBlock - "
            "relation is not from NewtonEulerR.");
      }
#ifdef MLCPPROJ_WITH_CT
      auto sizeDS = neds->dimension();
      auto T = neds->T();
      auto workT = std::make_shared<siconos::algebra::SiconosMatrix>(T);
      workT->trans();
      auto workT2 = std::make_shared<siconos::algebra::SiconosMatrix>(6, 6);
      siconos::algebra::prod(*workT, *T, *workT2, true);
      auto leftInteractionBlock =
          std::make_shared<siconos::algebra::SiconosMatrix>(sizeY, sizeDS);
      inter->getLeftInteractionBlockForDS(pos, leftInteractionBlock);
      auto work = std::make_shared<siconos::algebra::SiconosMatrix>(*leftInteractionBlock);
      std::cout << "LinearOSNS : leftUBlock\n";
      work->display();
      work->trans();
      std::cout << "LinearOSNS::computeInteractionBlock workT2\n";
      workT2->display();
      workT2->Solve(*work);
      siconos::algebra::prod(*leftInteractionBlock, *work, *currentInteractionBlock, false);
#else
      if (0)  //(std::static_pointer_cast<NewtonEulerR> inter->relation())->_isConstact){
      {
        //        auto sizeDS = neds->dimension();
        //        auto T = neds->T();
        //        auto workT = std::make_shared<siconos::algebra::SiconosMatrix>(*T));
        //        workT->trans();
        //        auto workT2 = std::make_shared<siconos::algebra::SiconosMatrix>(6, 6));
        //        siconos::algebra::prod(*workT, *T, *workT2, true);
        //        leftInteractionBlock1 =
        //        std::make_shared<siconos::algebra::SiconosMatrix>(sizeY, sizeDS));
        //        inter->getLeftInteractionBlockForDS(pos, leftInteractionBlock);
        //        leftInteractionBlock = std::make_shared<siconos::algebra::SiconosMatrix>(1,
        //        sizeDS)); for (auto ii = 0; ii < sizeDS; ii++)
        //          leftInteractionBlock->setValue(1, ii, leftInteractionBlock1->getValue(1,
        //          ii));
        //
        //        auto work =
        //        std::make_shared<siconos::algebra::SiconosMatrix>(*leftInteractionBlock));
        //        //cout<<"LinearOSNS : leftUBlock\n";
        //        //work->display();
        //        work->trans();
        //        //cout<<"LinearOSNS::computeInteractionBlock workT2"<<endl;
        //        //workT2->display();
        //        workT2->Solve(*work);
        //        siconos::algebra::prod(*leftInteractionBlock, *work,
        //        *currentInteractionBlock, false);
      } else {
        auto sizeDS =
            (std::static_pointer_cast<siconos::modeling::NewtonEulerDS>(ds))->getqDim();
        auto leftInteractionBlock =
            std::make_shared<siconos::algebra::SiconosMatrix>(sizeY, sizeDS);
        inter->getLeftInteractionBlockForDSProjectOnConstraints(pos, leftInteractionBlock);
        // #ifdef MLCPPROJ_DEBUG
        //          std::cout <<
        //          "siconos::nonsmooth_formulations::MLCPProjectOnConstraints::computeDiagonalInteractionBlock
        //          - NewtonEuler case leftInteractionBlock : \n";
        //         leftInteractionBlock->display();
        // #endif

        auto work = std::make_shared<siconos::algebra::SiconosMatrix>(*leftInteractionBlock);
        // cout<<"LinearOSNS sizeY="<<sizeY<<": leftUBlock\n";
        // work->display();
        work->transposeInPlace();
        siconos::algebra::prod(*leftInteractionBlock, *work, *currentInteractionBlock, false);
        // #ifdef MLCPPROJ_DEBUG
        //          std::cout <<
        //          "siconos::nonsmooth_formulations::MLCPProjectOnConstraints::computeDiagonalInteractionBlock
        //          - NewtonEuler case currentInteractionBlock : "<< "\n";
        //         currentInteractionBlock->display();
        // #endif
      }
#endif
    } else {
      THROW_EXCEPTION(
          "siconos::nonsmooth_formulations::MLCPProjectOnConstraints::"
          "computeDiagonalInteractionBlock - "
          "ds must be either NewtonEulerDS or LagrangianDS.");
    }
#ifdef MLCPPROJ_DEBUG
    std::cout << "siconos::nonsmooth_formulations::MLCPProjectOnConstraints:: "
                 "computeDiagonalInteractionBlock DiaginteractionBlock\n";
    currentInteractionBlock->display();
#endif
    // Set pos for next loop.
    pos = pos2;
  }
}

void siconos::nonsmooth_formulations::MLCPProjectOnConstraints::computeInteractionBlock(
    const siconos::graphs::InteractionsGraph::EDescriptor& ed) {
  // Computes matrix _interactionBlocks[inter1][inter2] (and allocates memory if
  // necessary) if inter1 and inter2 have commond DynamicalSystem.  How
  // _interactionBlocks are computed depends explicitely on the type of
  // Relation of each Interaction.

  // Warning: we suppose that at this point, all non linear
  // operators (G for lagrangian relation for example) have been
  // computed through plug-in mechanism.

#ifdef MLCPPROJ_DEBUG
  std::cout
      << "siconos::nonsmooth_formulations::MLCPProjectOnConstraints::computeInteractionBlock "
         "currentInteractionBlock start "
      << "\n";
#endif
  // Get dimension of the siconos::modeling::NonSmoothLaw (ie dim of the interactionBlock)
  auto indexSet = simulation()->indexSet(indexSetLevel());

  auto ds = indexSet->bundle(ed);
  auto inter1 = indexSet->bundle(indexSet->source(ed));
  auto inter2 = indexSet->bundle(indexSet->target(ed));
  auto& DSG0 = *simulation()->nonSmoothDynamicalSystem()->dynamicalSystems();
  auto& Osi = *DSG0.properties(DSG0.descriptor(ds)).osi;
  // For the edge 'ds', we need to find relative position of this ds
  // in inter1 and inter2 relation matrices (--> pos1 and pos2 below)
  // - find if ds is source or target in inter_i
  siconos::graphs::InteractionsGraph::VDescriptor vertex_inter;
  // - get the corresponding position
  unsigned int pos1, pos2;
  // source of inter1 :
  vertex_inter = indexSet->source(ed);
  auto tmpds = indexSet->properties(vertex_inter).source;
  if (tmpds == ds)
    pos1 = indexSet->properties(vertex_inter).source_pos;
  else {
    tmpds = indexSet->properties(vertex_inter).target;
    pos1 = indexSet->properties(vertex_inter).target_pos;
  }
  // now, inter2
  vertex_inter = indexSet->target(ed);
  tmpds = indexSet->properties(vertex_inter).source;
  if (tmpds == ds)
    pos2 = indexSet->properties(vertex_inter).source_pos;
  else {
    tmpds = indexSet->properties(vertex_inter).target;
    pos2 = indexSet->properties(vertex_inter).target_pos;
  }

  auto index1 = indexSet->index(indexSet->source(ed));
  auto index2 = indexSet->index(indexSet->target(ed));

  auto sizeY1 =
      std::static_pointer_cast<OSNSMatrixProjectOnConstraints>(_M)->computeSizeForProjection(
          inter1);
  auto sizeY2 =
      std::static_pointer_cast<OSNSMatrixProjectOnConstraints>(_M)->computeSizeForProjection(
          inter2);

  std::shared_ptr<siconos::algebra::SiconosMatrix> currentInteractionBlock;

  assert(index1 != index2);

  if (index2 > index1)  // upper block
  {
    //     if (! indexSet->properties(ed).upper_block)
    //     {
    //       indexSet->properties(ed).upper_block =
    //       std::make_shared<siconos::algebra::SiconosMatrix>(sizeY1, sizeY2));
    //     }

    currentInteractionBlock = indexSet->upper_blockProj[ed];
#ifdef MLCPPROJ_DEBUG
    std::cout << "siconos::nonsmooth_formulations::MLCPProjectOnConstraints::"
                 "computeInteractionBlock "
                 "currentInteractionBlock "
              << "\n";
    //    currentInteractionBlock->display();
    std::cout << "sizeY1 " << sizeY1 << "\n";
    std::cout << "sizeY2 " << sizeY2 << "\n";
    std::cout << "upper_blockProj " << indexSet->upper_blockProj[ed].get() << " of edge " << ed
              << " of size " << currentInteractionBlock->size(0) << " x "
              << currentInteractionBlock->size(0) << " for interaction " << inter1->number()
              << " and interaction " << inter2->number() << "\n";
    // std::cout<<"inter1->display() "<< inter1->number()<< "\n";
    // inter1->display();
    // std::cout<<"inter2->display() "<< inter2->number()<< "\n";
    // inter2->display();

#endif
    assert(currentInteractionBlock->size(0) == sizeY1);
    assert(currentInteractionBlock->size(1) == sizeY2);
  } else  // lower block
  {
    //     if (! indexSet->properties(ed).lower_block)
    //     {
    //       indexSet->properties(ed).lower_block =
    //       std::make_shared<siconos::algebra::SiconosMatrix>(sizeY1, sizeY2));
    //     }

    assert(indexSet->lower_blockProj[ed]->size(0) == sizeY1);
    assert(indexSet->lower_blockProj[ed]->size(1) == sizeY2);

    currentInteractionBlock = indexSet->lower_blockProj[ed];
  }
  // General form of the interactionBlock is : interactionBlock =
  // a*extraInteractionBlock + b * leftInteractionBlock * centralInteractionBlocks
  // * rightInteractionBlock a and b are scalars, centralInteractionBlocks a
  // matrix depending on the integrator (and on the DS), the
  // simulation type ...  left, right and extra depend on the relation
  // type and the non smooth law.
  auto relationType1 = inter1->relation()->getType();
  auto relationType2 = inter2->relation()->getType();
  if (relationType1 == siconos::modeling::RelationType::NewtonEuler &&
      relationType2 == siconos::modeling::RelationType::NewtonEuler) {
    assert(inter1 != inter2);
    currentInteractionBlock->setZero();
#ifdef MLCPPROJ_WITH_CT
    auto sizeDS = (std::static_pointer_cast<NewtonEulerDS>(ds))->dimension();
    auto leftInteractionBlock =
        std::make_shared<siconos::algebra::SiconosMatrix>(sizeY1, sizeDS);
    inter1->getLeftInteractionBlockForDS(pos1, leftInteractionBlock);
    auto neds = (std::static_pointer_cast<NewtonEulerDS>(ds));
    auto T = neds->T();
    auto workT = std::make_shared<siconos::algebra::SiconosMatrix>(T);
    workT->trans();
    auto workT2 = std::make_shared<siconos::algebra::SiconosMatrix>(6, 6);
    siconos::algebra::prod(*workT, *T, *workT2, true);
    auto rightInteractionBlock =
        std::make_shared<siconos::algebra::SiconosMatrix>(sizeY2, sizeDS);
    inter2->getLeftInteractionBlockForDS(pos2, rightInteractionBlock);
    rightInteractionBlock->trans();
    workT2->Solve(*rightInteractionBlock);
    siconos::algebra::prod(*leftInteractionBlock, *rightInteractionBlock,
                           *currentInteractionBlock, false);

#else

    auto sizeDS = (std::static_pointer_cast<siconos::modeling::NewtonEulerDS>(ds))->getqDim();
    auto leftInteractionBlock =
        std::make_shared<siconos::algebra::SiconosMatrix>(sizeY1, sizeDS);
    inter1->getLeftInteractionBlockForDSProjectOnConstraints(pos1, leftInteractionBlock);
    auto neds = (std::static_pointer_cast<siconos::modeling::NewtonEulerDS>(ds));
    auto rightInteractionBlock =
        std::make_shared<siconos::algebra::SiconosMatrix>(sizeY2, sizeDS);
    inter2->getLeftInteractionBlockForDSProjectOnConstraints(pos2, rightInteractionBlock);
    rightInteractionBlock->transposeInPlace();
    siconos::algebra::prod(*leftInteractionBlock, *rightInteractionBlock,
                           *currentInteractionBlock, false);
#endif
  } else if (relationType1 == siconos::modeling::RelationType::Lagrangian &&
             relationType2 == siconos::modeling::RelationType::Lagrangian) {
    auto sizeDS = ds->dimension();
    auto leftInteractionBlock = inter1->getLeftInteractionBlockForDS(pos1, sizeY1, sizeDS);
    if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      if (d->boundaryConditions())  // V.A. Should we do that ?
      {
        for (const auto itindex : d->boundaryConditions()->velocityIndices()) {
          auto coltmp = std::make_shared<siconos::algebra::SiconosVector>(sizeY1);
          coltmp->setZero();
          leftInteractionBlock->col(itindex) = *coltmp;
        }
      }
    }
#ifdef MLCPPROJ_DEBUG
    std::cout << "siconos::nonsmooth_formulations::MLCPProjectOnConstraints::"
                 "computeInteractionBlock "
                 ": leftInteractionBlock"
              << "\n";
    leftInteractionBlock->display();
#endif
    // inter1 != inter2
    auto rightInteractionBlock = inter2->getLeftInteractionBlockForDS(pos2, sizeY2, sizeDS);
#ifdef MLCPPROJ_DEBUG
    std::cout << "siconos::nonsmooth_formulations::MLCPProjectOnConstraints::"
                 "computeInteractionBlock "
                 ": rightInteractionBlock"
              << "\n";
    rightInteractionBlock->display();
#endif
    // Warning: we use getLeft for Right interactionBlock
    // because right = transpose(left) and because of
    // size checking inside the getBlock function, a
    // getRight call will fail.
    auto centralInteractionBlock = getOSIMatrix(Osi, ds);
#ifdef MLCPPROJ_DEBUG
    std::cout << "siconos::nonsmooth_formulations::MLCPProjectOnConstraints::"
                 "computeInteractionBlock "
                 ": centralInteractionBlocks "
              << "\n";
    centralInteractionBlock->display();
#endif
    rightInteractionBlock->transposeInPlace();

    if (_useMassNormalization) {
      *rightInteractionBlock = centralInteractionBlock->solve(*rightInteractionBlock);
      //*currentInteractionBlock +=  *leftInteractionBlock ** work;
      *currentInteractionBlock += *leftInteractionBlock * *rightInteractionBlock;
    } else {
      *currentInteractionBlock += *leftInteractionBlock * *rightInteractionBlock;
    }
#ifdef MLCPPROJ_DEBUG
    std::cout << "siconos::nonsmooth_formulations::MLCPProjectOnConstraints::"
                 "computeInteractionBlock "
                 ": currentInteractionBlock"
              << "\n";
    currentInteractionBlock->display();
#endif
  }

  else
    THROW_EXCEPTION(
        "siconos::nonsmooth_formulations::MLCPProjectOnConstraints::computeInteractionBlock "
        "not yet "
        "implemented for relation of type " +
        siconos::tools::enum_to_string(relationType1));
}

void siconos::nonsmooth_formulations::MLCPProjectOnConstraints::computeqBlock(
    siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter, unsigned int pos) {
  DEBUG_BEGIN(
      "siconos::nonsmooth_formulations::MLCPProjectOnConstraints::computeqBlock(siconos::"
      "graphs::"
      "InteractionsGraph::"
      "VDescriptor& vertex_inter, unsigned int pos)\n");

  auto indexSet = simulation()->indexSet(indexSetLevel());
  auto inter = indexSet->bundle(vertex_inter);
  auto sizeY =
      std::static_pointer_cast<OSNSMatrixProjectOnConstraints>(_M)->computeSizeForProjection(
          inter);
  DEBUG_PRINTF("pos = %i", pos);
  for (decltype(sizeY) i = 0; i < sizeY; i++)
    _q->setValue(pos + i, inter->y(0)->getValue(0 + i));

  DEBUG_EXPR(_q->display(););
  DEBUG_END(
      "siconos::nonsmooth_formulations::MLCPProjectOnConstraints::computeqBlock(siconos::"
      "graphs::"
      "InteractionsGraph::"
      "VDescriptor& vertex_inter, unsigned int pos)\n");
}

void siconos::nonsmooth_formulations::MLCPProjectOnConstraints::computeq(double time) {
  if (_q->size() != _sizeOutput) _q->resize(_sizeOutput);
  _q->setZero();

  // === Get index set from Simulation ===
  auto indexSet = simulation()->indexSet(indexSetLevel());
  // === Loop through "active" Interactions (ie present in
  // indexSets[level]) ===

  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  for (std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui) {
    // Compute q, this depends on the type of non smooth problem, on
    // the relation type and on the non smooth law
    auto pos = indexSet->properties(*ui).absolute_position_proj;
    computeqBlock(*ui, pos);  // free output is saved in y
  }
}

void siconos::nonsmooth_formulations::MLCPProjectOnConstraints::postCompute() {
  _hasBeenUpdated = true;
  // This function is used to set y/lambda values using output from
  // lcp_driver (w,z).  Only Interactions (ie Interactions) of
  // indexSet(leveMin) are concerned.

  // === Get index set from Topology ===
  auto indexSet = simulation()->indexSet(indexSetLevel());

  // === Loop through "active" Interactions (ie present in
  // indexSets[1]) ===
  /** We chose to do a small step _alpha in view of stabilized the algorithm.*/
#ifdef MLCPPROJ_DEBUG
  printf(
      "siconos::nonsmooth_formulations::MLCPProjectOnConstraints::postCompute damping value = "
      "%f\n",
      _alpha);
#endif
  (*_z) *= _alpha;
#ifdef MLCPPROJ_DEBUG
  printf("siconos::nonsmooth_formulations::MLCPProjectOnConstraints::postCompute _z\n");
  _z->display();
  display();
#endif

  siconos::graphs::InteractionsGraph::VIterator ui, uiend;

  for (std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui) {
    auto inter = indexSet->bundle(*ui);
    // Get the relative position of inter-interactionBlock in the vector w
    // or z
    auto pos = indexSet->properties(*ui).absolute_position_proj;
    auto relationType = inter->relation()->getType();
    if (relationType == siconos::modeling::RelationType::NewtonEuler) {
      postComputeNewtonEulerR(inter, pos);
    } else if (relationType == siconos::modeling::RelationType::Lagrangian) {
      postComputeLagrangianR(inter, pos);
    } else {
      THROW_EXCEPTION(
          "siconos::nonsmooth_formulations::MLCPProjectOnConstraints::computeInteractionBlock "
          "- "
          "relation type is not from "
          "Lagrangian type neither NewtonEuler.");
    }
  }
}

void siconos::nonsmooth_formulations::MLCPProjectOnConstraints::postComputeLagrangianR(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int pos) {
  auto lr = std::static_pointer_cast<siconos::modeling::LagrangianR>(inter->relation());
#ifdef MLCPPROJ_DEBUG
  printf(
      "siconos::nonsmooth_formulations::MLCPProjectOnConstraints::postComputeLagrangian "
      "inter->y(0)\n");
  inter->y(0)->display();
  printf(
      "siconos::nonsmooth_formulations::MLCPProjectOnConstraints::postComputeLagrangian "
      "lr->jachq \n");
  lr->jachq()->display();
  printf(
      "siconos::nonsmooth_formulations::MLCPProjectOnConstraints::postComputeLagrangianR q "
      "before "
      "update\n");

  auto indexSet = simulation()->indexSet(indexSetLevel());
  auto ui = indexSet->descriptor(inter);
  siconos::graphs::InteractionsGraph::OEIterator oei, oeiend;
  for (std::tie(oei, oeiend) = indexSet->out_edges(ui); oei != oeiend; ++oei) {
    auto lds = std::static_pointer_cast<LagrangianDS>(indexSet->bundle(*oei));
    lds->q()->display();
  }
#endif

  // auto sizeY =inter->nonSmoothLaw()->size();

  // y and lambda vectors
  auto lambda = inter->lambda(0);
  auto y = inter->y(0);
  auto sizeY =
      std::static_pointer_cast<OSNSMatrixProjectOnConstraints>(_M)->computeSizeForProjection(
          inter);
  // Copy _w/_z values, starting from index pos into y/lambda.

  *lambda = _z->segment(pos, sizeY);

#ifdef MLCPPROJ_DEBUG
  printf("MLCPP lambda of Interaction is pos =%i :\n", pos);
  //  aBuff->display();
  lambda->display();
  auto nslawsize = inter->nonSmoothLaw()->size();
  auto aBuff = std::make_shared<siconos::algebra::SiconosVector>(nslawsize);
  *aBuff = _z->segment(pos, sizeY);
  auto J = lr->jachq();
  auto aux = std::make_shared<siconos::algebra::SiconosMatrix>(*J);
  aux->trans();
  // std::shared_ptr<siconos::algebra::SiconosVector> tmp =
  // std::make_shared<siconos::algebra::SiconosVector>(*(lr->q())));
  // siconos::algebra::prod(*aux, *aBuff, *(tmp), false);
  // //siconos::algebra::prod(*aux,*lambda,*(lr->q()),false);
  // std:: std::cout << " tmp =  tmp + J^T * lambda\n";
  // tmp->display();
#endif

  // // WARNING : Must not be done here. and should be called with the correct time.
  // // compute p(0)
  // inter->computeInput(0.0 ,0);

  // // \warning aBuff should normally be in lambda[0]
  // // The update of the position in DS should be made
  // //  in MoreauJeanOSI::upateState or ProjectedMoreauJeanOSI::updateState
  // auto J=lr->jachq();
  // auto aux = std::make_shared<siconos::algebra::SiconosMatrix>(*J));
  // aux->trans();

  // std::shared_ptr<siconos::algebra::SiconosVector> tmp  =
  // std::make_shared<siconos::algebra::SiconosVector>(*(lr->q()))); std:: std::cout << "
  // tmp ="<<std::endl; tmp->display(); std:: std::cout << " lr->q() ="<<std::endl;
  // lr->q()->display();

  // //siconos::algebra::prod(*aux,*lambda,*(lr->q()),false);
  // siconos::algebra::prod(*aux,*aBuff,*(tmp),false);
  // std:: std::cout << " tmp =  tmp + J * lambda"<<std::endl;
  // tmp->display();

  // // The following step should be done on MoreauJeanOSI::upateState or
  // ProjectedMoreauJeanOSI::updateState DSIterator itDS = inter->dynamicalSystemsBegin();
  // while(itDS!=inter->dynamicalSystemsEnd())
  // {
  //   Type::Siconos dsType = Type::value(**itDS);
  //   if((dsType !=Type::LagrangianDS) and
  //      (dsType !=Type::LagrangianLinearTIDS) )
  //   {
  //     THROW_EXCEPTION("MLCPProjectOnConstraint::postCompute- ds is not of Lagrangian DS
  //     type.");
  //   }

  //   auto d = std::static_pointer_cast<LagrangianDS> (*itDS);
  //   std::shared_ptr<siconos::algebra::SiconosVector> q = d->q();

  //   *q +=  *d->p(0);
  //    std::cout << " q=\n";
  //   q->display();
  //   itDS++;
  // }

  // if ((*lr->q() - *tmp).normInf() > 1e-12)
  // {
  //   THROW_EXCEPTION("youyou");
  // }

#ifdef MLCPPROJ_DEBUG
  printf(
      "siconos::nonsmooth_formulations::MLCPProjectOnConstraints::postComputeLagrangianR "
      "_z\n");
  _z->display();
  printf(
      "siconos::nonsmooth_formulations::MLCPProjectOnConstraints::postComputeLagrangianR "
      "updated\n");

  auto& DSlink = *(indexSet->properties(ui)).DSlink;
//  (*DSlink[siconos::modeling::LagrangianR::q0]).display();
//  (lr->q())->display();
#endif

  // THROW_EXCEPTION("siconos::nonsmooth_formulations::MLCPProjectOnConstraints::postComputeLagrangianR()
  // - not yet implemented");
}

void siconos::nonsmooth_formulations::MLCPProjectOnConstraints::postComputeNewtonEulerR(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int pos) {
  auto ner = (std::static_pointer_cast<siconos::modeling::NewtonEulerR>(inter->relation()));
  auto lambda = inter->lambda(0);
  auto y = inter->y(0);
  auto sizeY =
      std::static_pointer_cast<OSNSMatrixProjectOnConstraints>(_M)->computeSizeForProjection(
          inter);
  // Copy _w/_z values, starting from index pos into y/lambda.

  *lambda = _z->segment(pos, sizeY);
}

void siconos::nonsmooth_formulations::MLCPProjectOnConstraints::computeOptions(
    std::shared_ptr<siconos::modeling::Interaction> inter1,
    std::shared_ptr<siconos::modeling::Interaction> inter2) {
  //  printf("siconos::nonsmooth_formulations::MLCPProjectOnConstraints::computeOptions\n");
  // Get dimension of the siconos::modeling::NonSmoothLaw (ie dim of the interactionBlock)
  // Retrieve size of Y (projected variable)
  auto sizeY1 =
      std::static_pointer_cast<OSNSMatrixProjectOnConstraints>(_M)->computeSizeForProjection(
          inter1);

  // Compute the number of equalities
  auto equalitySize1 = sizeY1;  // default behavior

  if (siconos::types::type_value(*(inter1->nonSmoothLaw())) ==
          siconos::modeling::Type::NewtonImpactFrictionNSL ||
      siconos::types::type_value(*(inter1->nonSmoothLaw())) ==
          siconos::modeling::Type::NewtonImpactNSL) {
    if (_doProjOnEquality) {
      equalitySize1 = sizeY1;
    } else {
      equalitySize1 = 0;
    }
  } else if (siconos::types::type_value(*(inter1->nonSmoothLaw())) ==
             siconos::modeling::Type::MixedComplementarityConditionNSL) {
    equalitySize1 =
        std::static_pointer_cast<siconos::modeling::MixedComplementarityConditionNSL>(
            inter1->nonSmoothLaw())
            ->equalitySize();
  }

  // Compute the number of inequalities
  auto inequalitySize1 = sizeY1 - equalitySize1;

  if (inter1 == inter2) {
    // inter1->getExtraInteractionBlock(currentInteractionBlock);
    _m += inequalitySize1;
    _n += equalitySize1;
    //    _m=0;
    //_n=6;
    if (_curBlock > MLCP_NB_BLOCKS_MAX - 2)
      printf("MLCP.cpp : number of block to small, memory crach below!!!\n");
    /*add an equality block.*/

    // #ifdef MLCPPROJ_DEBUG
    //   printf("siconos::nonsmooth_formulations::MLCPProjectOnConstraints::computeOptions()\n");
    // #endif

    if (equalitySize1 > 0) {
      _numerics_problem->blocksRows[_curBlock + 1] =
          _numerics_problem->blocksRows[_curBlock] + equalitySize1;
      _numerics_problem->blocksIsComp[_curBlock] = 0;
      // #ifdef MLCPPROJ_DEBUG
      //        std::cout << "_curBlock : " << _curBlock <<std::endl;
      //        std::cout << "_numerics_problem->blocksRows["<<_curBlock+1 <<" ] : " <<
      //        _numerics_problem->blocksRows[_curBlock+1] <<std::endl; std::cout <<
      //        "_numerics_problem->blocksIsComp["<<_curBlock <<" ] : " <<
      //        _numerics_problem->blocksIsComp[_curBlock] <<std::endl;
      // #endif

      _curBlock++;
    }
    /*add a complementarity block.*/
    if (inequalitySize1 > 0) {
      _numerics_problem->blocksRows[_curBlock + 1] =
          _numerics_problem->blocksRows[_curBlock] + inequalitySize1;
      _numerics_problem->blocksIsComp[_curBlock] = 1;
      // #ifdef MLCPPROJ_DEBUG
      //        std::cout << "_curBlock : " << _curBlock <<std::endl;
      //        std::cout << "_numerics_problem->blocksRows["<<_curBlock+1<< "] : " <<
      //        _numerics_problem->blocksRows[_curBlock+1] <<std::endl; std::cout <<
      //        "_numerics_problem->blocksIsComp["<<_curBlock<< "] : " <<
      //        _numerics_problem->blocksIsComp[_curBlock] <<std::endl;
      // #endif

      _curBlock++;
    }
  }
  // #ifdef MLCPPROJ_DEBUG
  //    std::cout << "_m : " << _m <<std::endl;
  //    std::cout << "_n : " << _n <<std::endl;
  // #endif
}
