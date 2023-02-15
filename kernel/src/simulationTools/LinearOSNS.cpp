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
#include "LinearOSNS.hpp"

#include "BoundaryCondition.hpp"
#include "D1MinusLinearOSI.hpp"
#include "DynamicalSystem.hpp"
#include "EulerMoreauOSI.hpp"
#include "Interaction.hpp"
#include "LsodarOSI.hpp"
#include "MoreauJeanBilbaoOSI.hpp"
#include "MoreauJeanOSI.hpp"
#include "NewMarkAlphaOSI.hpp"
#include "NonSmoothLaw.hpp"
#include "OneStepIntegrator.hpp"
#include "Relation.hpp"
#include "SchatzmanPaoliOSI.hpp"
#include "SecondOrderDS.hpp"
#include "SiconosConst.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrixOp.hpp"  // for prod
#include "SiconosVector.hpp"
#include "SiconosVectorOp.hpp"  // for setBlock
#include "SimpleMatrix.hpp"
#include "Simulation.hpp"
#include "Topology.hpp"
#include "ZeroOrderHoldOSI.hpp"
// #include "NewtonEulerR.hpp"
// #include "FirstOrderLinearR.hpp"
// #include "FirstOrderLinearTIR.hpp"
#include "LagrangianCompliantLinearTIR.hpp"
#include "LagrangianLinearDiagonalDS.hpp"
// #include "NewtonImpactNSL.hpp"
// #include "MultipleImpactNSL.hpp"
// #include "NewtonImpactFrictionNSL.hpp"
// #include "NonSmoothDynamicalSystem.hpp"
// #include "LagrangianRheonomousR.hpp"
// #include "LagrangianScleronomousR.hpp"
// #include "LagrangianLinearTIDS.hpp"
// #include "NewtonEulerDS.hpp"
#include "OSNSMatrix.hpp"

// #include "Tools.hpp"
// #include <chrono>

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"
// #define WITH_TIMER

void siconos::nonsmooth_formulations::LinearOSNS::initVectorsMemory() {
  // Memory allocation for _w, M, z and q.
  // If one of them has already been allocated, nothing is done.
  // We suppose that user has chosen a correct size.
  if (!_w)
    _w = std::make_shared<siconos::algebra::SiconosVector>(maxSize());
  else {
    if (_w->size() != maxSize()) _w->resize(maxSize());
  }

  if (!_z)
    _z = std::make_shared<siconos::algebra::SiconosVector>(maxSize());
  else {
    if (_z->size() != maxSize()) _z->resize(maxSize());
  }

  if (!_q)
    _q = std::make_shared<siconos::algebra::SiconosVector>(maxSize());
  else {
    if (_q->size() != maxSize()) _q->resize(maxSize());
  }
}
void siconos::nonsmooth_formulations::LinearOSNS::initOSNSMatrix() {
  // Default size for M = maxSize()
  if (_assemblyType == LinearOSNSAssemblyType::REDUCED_BLOCK or
      _assemblyType == LinearOSNSAssemblyType::REDUCED_DIRECT) {
    if (!_M) {
      switch (_numericsMatrixStorageType) {
        case NM_DENSE:
        case NM_SPARSE: {
          _M = std::make_shared<OSNSMatrix>(0, _numericsMatrixStorageType);
          break;
        }
        case NM_SPARSE_BLOCK: {
          // = number of Interactionin the largest considered indexSet
          if (indexSetLevel() != siconos::internal::LEVELMAX &&
              simulation()->nonSmoothDynamicalSystem()->topology()->indexSetsSize() >
                  indexSetLevel()) {
            _M = std::make_shared<OSNSMatrix>(simulation()->indexSet(indexSetLevel())->size(),
                                              _numericsMatrixStorageType);
          } else {
            _M = std::make_shared<OSNSMatrix>(1, _numericsMatrixStorageType);
          }
          break;
        }
          {
            default:
              THROW_EXCEPTION(
                  "siconos::nonsmooth_formulations::LinearOSNS::initOSNSMatrix unknown "
                  "_storageType");
          }
      }
    }
  }

  if (_assemblyType == LinearOSNSAssemblyType::GLOBAL) {
    // Default size for M = _maxSize
    if (!_W) {
      // if (_numericsMatrixStorageType == NM_DENSE)
      //   _W = std::make_shared<OSNSMatrix>(_maxSize, NM_DENSE));
      // else // if(MStorageType == 1) size = number of DSBlocks = number of DS in the largest
      // considered graph of ds
      //   _W =
      //   std::make_shared<OSNSMatrix>(simulation()->nonSmoothDynamicalSystem()->dynamicalSystems()->size(),
      //   1));

      switch (_numericsMatrixStorageType) {
        case NM_DENSE: {
          _W = std::make_shared<OSNSMatrix>(_maxSize, NM_DENSE);
          break;
        }
        case NM_SPARSE: {
          _W = std::make_shared<OSNSMatrix>(0, NM_SPARSE);
          break;
        }
        case NM_SPARSE_BLOCK: {
          _W = std::make_shared<OSNSMatrix>(
              simulation()->nonSmoothDynamicalSystem()->dynamicalSystems()->size(),
              NM_SPARSE_BLOCK);
          break;
        }
          {
            default:
              THROW_EXCEPTION(
                  "siconos::nonsmooth_formulations::LinearOSNS::initOSNSMatrix unknown "
                  "_storageType");
          }
      }
    }

    if (!_H) {
      switch (_numericsMatrixStorageType) {
        case NM_DENSE: {
          _H = std::make_shared<OSNSMatrix>(_maxSize, NM_DENSE);
          break;
        }
        case NM_SPARSE: {
          _H = std::make_shared<OSNSMatrix>(0, simulation()->indexSet(_indexSetLevel)->size(),
                                            NM_SPARSE);
          break;
        }
        case NM_SPARSE_BLOCK: {
          _H = std::make_shared<OSNSMatrix>(
              simulation()->nonSmoothDynamicalSystem()->dynamicalSystems()->size(),
              simulation()->indexSet(_indexSetLevel)->size(), NM_SPARSE_BLOCK);
          break;
        }
          {
            default:
              THROW_EXCEPTION(
                  "siconos::nonsmooth_formulations::LinearOSNS::initOSNSMatrix unknown "
                  "_storageType");
          }
      }
    }
  }
  if (_assemblyType == LinearOSNSAssemblyType::REDUCED_DIRECT) {
    // Default size for M = _maxSize
    if (!_W_inverse) {
      // if (_numericsMatrixStorageType == NM_DENSE)
      //   _W_inverse = std::make_shared<OSNSMatrix>(_maxSize, NM_DENSE));
      // else // if(MStorageType == 1) size = number of DSBlocks = number of DS in the largest
      // considered graph of ds
      //   _W_inverse =
      //   std::make_shared<OSNSMatrix>(simulation()->nonSmoothDynamicalSystem()->dynamicalSystems()->size(),
      //   1));

      switch (_numericsMatrixStorageType) {
        case NM_DENSE: {
          _W_inverse = std::make_shared<OSNSMatrix>(_maxSize, NM_DENSE);
          break;
        }
        case NM_SPARSE: {
          _W_inverse = std::make_shared<OSNSMatrix>(0, NM_SPARSE);
          break;
        }
        case NM_SPARSE_BLOCK: {
          _W_inverse = std::make_shared<OSNSMatrix>(
              simulation()->nonSmoothDynamicalSystem()->dynamicalSystems()->size(),
              NM_SPARSE_BLOCK);
          break;
        }
          {
            default:
              THROW_EXCEPTION(
                  "siconos::nonsmooth_formulations::LinearOSNS::initOSNSMatrix unknown "
                  "_storageType");
          }
      }
    }

    if (!_H) {
      switch (_numericsMatrixStorageType) {
        case NM_DENSE: {
          _H = std::make_shared<OSNSMatrix>(_maxSize, NM_DENSE);
          break;
        }
        case NM_SPARSE: {
          _H = std::make_shared<OSNSMatrix>(0, simulation()->indexSet(_indexSetLevel)->size(),
                                            NM_SPARSE);
          break;
        }
        case NM_SPARSE_BLOCK: {
          _H = std::make_shared<OSNSMatrix>(
              simulation()->nonSmoothDynamicalSystem()->dynamicalSystems()->size(),
              simulation()->indexSet(_indexSetLevel)->size(), NM_SPARSE_BLOCK);
          break;
        }
          {
            default:
              THROW_EXCEPTION(
                  "siconos::nonsmooth_formulations::LinearOSNS::initOSNSMatrix unknown "
                  "_storageType");
          }
      }
    }
  }
  if (_assemblyType == LinearOSNSAssemblyType::GLOBAL_REDUCED) {
    // Default size for M = _maxSize
    if (!_W) {
      // if (_numericsMatrixStorageType == NM_DENSE)
      //   _W = std::make_shared<OSNSMatrix>(_maxSize, NM_DENSE));
      // else // if(MStorageType == 1) size = number of DSBlocks = number of DS in the largest
      // considered graph of ds
      //   _W =
      //   std::make_shared<OSNSMatrix>(simulation()->nonSmoothDynamicalSystem()->dynamicalSystems()->size(),
      //   1));

      switch (_numericsMatrixStorageType) {
        case NM_DENSE: {
          _W = std::make_shared<OSNSMatrix>(_maxSize, NM_DENSE);
          break;
        }
        case NM_SPARSE: {
          _W = std::make_shared<OSNSMatrix>(0, NM_SPARSE);
          break;
        }
        case NM_SPARSE_BLOCK: {
          _W = std::make_shared<OSNSMatrix>(
              simulation()->nonSmoothDynamicalSystem()->dynamicalSystems()->size(),
              NM_SPARSE_BLOCK);
          break;
        }
          {
            default:
              THROW_EXCEPTION(
                  "siconos::nonsmooth_formulations::LinearOSNS::initOSNSMatrix unknown "
                  "_storageType");
          }
      }
    }
    // Default size for M = _maxSize
    if (!_W_inverse) {
      switch (_numericsMatrixStorageType) {
        case NM_DENSE: {
          _W_inverse = std::make_shared<OSNSMatrix>(_maxSize, NM_DENSE);
          break;
        }
        case NM_SPARSE: {
          _W_inverse = std::make_shared<OSNSMatrix>(0, NM_SPARSE);
          break;
        }
        case NM_SPARSE_BLOCK: {
          _W_inverse = std::make_shared<OSNSMatrix>(
              simulation()->nonSmoothDynamicalSystem()->dynamicalSystems()->size(),
              NM_SPARSE_BLOCK);
          break;
        }
          {
            default:
              THROW_EXCEPTION(
                  "siconos::nonsmooth_formulations::LinearOSNS::initOSNSMatrix unknown "
                  "_storageType");
          }
      }
    }

    if (!_H) {
      switch (_numericsMatrixStorageType) {
        case NM_DENSE: {
          _H = std::make_shared<OSNSMatrix>(_maxSize, NM_DENSE);
          break;
        }
        case NM_SPARSE: {
          _H = std::make_shared<OSNSMatrix>(0, simulation()->indexSet(_indexSetLevel)->size(),
                                            NM_SPARSE);
          break;
        }
        case NM_SPARSE_BLOCK: {
          _H = std::make_shared<OSNSMatrix>(
              simulation()->nonSmoothDynamicalSystem()->dynamicalSystems()->size(),
              simulation()->indexSet(_indexSetLevel)->size(), NM_SPARSE_BLOCK);
          break;
        }
          {
            default:
              THROW_EXCEPTION(
                  "siconos::nonsmooth_formulations::LinearOSNS::initOSNSMatrix unknown "
                  "_storageType");
          }
      }
    }
  }
}
void siconos::nonsmooth_formulations::LinearOSNS::initialize(
    std::shared_ptr<siconos::simulation::Simulation> sim) {
  // - Checks memory allocation for main variables (_M,q,_w,z)
  // - Formalizes the problem if the topology is time-invariant

  // This function performs all steps that are time-invariant

  // General initialize for OneStepNSProblem
  OneStepNSProblem::initialize(sim);

  initVectorsMemory();

  // Note that _interactionBlocks is up to date since updateInteractionBlocks
  // has been called during OneStepNSProblem::initialize()

  initOSNSMatrix();
}

void siconos::nonsmooth_formulations::LinearOSNS::computeDiagonalInteractionBlock(
    const siconos::graphs::InteractionsGraph::VDescriptor& vd) {
  DEBUG_BEGIN(
      "siconos::nonsmooth_formulations::LinearOSNS::computeDiagonalInteractionBlock(const "
      "siconos::graphs::InteractionsGraph::VDescriptor& vd)\n");

  // Compute diagonal block (graph property) for a given interaction.
  // - How  blocks are computed depends explicitely on the nonsmooth law, the relation and the
  // dynamical system types.
  // - No memory allocation there : blocks should be previously (updateInteractionBlocks) and
  // properly allocated.

  // Get dimension of the NonSmoothLaw (ie dim of the interactionBlock)
  auto indexSet = simulation()->indexSet(indexSetLevel());
  std::shared_ptr<siconos::modeling::Interaction> inter = indexSet->bundle(vd);
  // Get osi property from interaction

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

  auto& DSG0 = *simulation()->nonSmoothDynamicalSystem()->dynamicalSystems();
  auto nslaw = inter->nonSmoothLaw();

  // --- Check compatible nslaws ----
  checkCompatibleNSLaw(*nslaw);

  // --- Check block size ---
  assert(indexSet->properties(vd).block->size(0) == nslaw->size());
  assert(indexSet->properties(vd).block->size(1) == nslaw->size());

  // --- Compute diagonal block ---
  // Block to be set in OSNS Matrix, corresponding to
  // the current interaction
  auto currentInteractionBlock = indexSet->properties(vd).block;
  std::shared_ptr<siconos::algebra::SiconosMatrix> leftInteractionBlock, rightInteractionBlock;

  auto relationType = inter->relation()->getType();
  auto relationSubType = inter->relation()->getSubType();
  double h = simulation()->currentTimeStep();

  // General form of the interactionBlock is : interactionBlock =
  // a*extraInteractionBlock + b * leftInteractionBlock * centralInteractionBlocks
  // * rightInteractionBlock a and b are scalars, centralInteractionBlocks a
  // matrix depending on the integrator (and on the DS), the
  // simulation type ...  left, right and extra depend on the relation
  // type and the non smooth law.
  relationType = inter->relation()->getType();

  inter->getExtraInteractionBlock(currentInteractionBlock);

  unsigned int nslawSize = nslaw->size();
  // loop over the DS connected to the interaction.
  bool endl = false;
  unsigned int pos = pos1;
  for (std::shared_ptr<siconos::modeling::DynamicalSystem> ds = ds1; !endl; ds = ds2) {
    assert(ds == ds1 || ds == ds2);
    endl = (ds == ds2);

    auto& osi = *DSG0.properties(DSG0.descriptor(ds)).osi;
    auto osiType = osi.getType();
    auto sizeDS = ds->dimension();

    // get _interactionBlocks corresponding to the current DS
    // These _interactionBlocks depends on the relation type.
    leftInteractionBlock = inter->getLeftInteractionBlockForDS(pos, nslawSize, sizeDS);
    DEBUG_EXPR(leftInteractionBlock->display(););
    // Computing depends on relation type -> move this in Interaction method?
    if (relationType == siconos::modeling::RelationType::FirstOrder) {
      ////// TODO : BACK WHEN NAMESPACE OK

      rightInteractionBlock = inter->getRightInteractionBlockForDS(pos, sizeDS, nslawSize);

      if (osiType == siconos::integrators::IntegratorType::EULERMOREAUOSI) {
        if ((static_cast<siconos::integrators::EulerMoreauOSI&>(osi)).useGamma() ||
            (static_cast<siconos::integrators::EulerMoreauOSI&>(osi)).useGammaForRelation()) {
          *rightInteractionBlock *=
              (static_cast<siconos::integrators::EulerMoreauOSI&>(osi)).gamma();
        }
      }

      // for ZOH, we have a different formula ...
      if ((osiType == siconos::integrators::IntegratorType::ZOHOSI) &&
          indexSet->properties(vd).forControl) {
        *rightInteractionBlock =
            static_cast<siconos::integrators::ZeroOrderHoldOSI&>(osi).Bd(ds);
        siconos::algebra::prod(*leftInteractionBlock, *rightInteractionBlock,
                               *currentInteractionBlock, false);
      } else {
        // centralInteractionBlock contains a lu-factorized matrix and we solve
        // centralInteractionBlock * X = rightInteractionBlock with PLU
        auto centralInteractionBlock = getOSIMatrix(osi, ds);
        centralInteractionBlock->Solve(*rightInteractionBlock);
        auto& workMInter = *indexSet->properties(vd).workMatrices;
        static_cast<siconos::integrators::EulerMoreauOSI&>(osi).computeKhat(
            *inter, *rightInteractionBlock, workMInter, h);

        //      integration of r with theta method removed
        //      *currentInteractionBlock += h *Theta[*itDS]* *leftInteractionBlock *
        // (*rightInteractionBlock);  // left = C, right = W.B
        // gemm(h,*leftInteractionBlock,*rightInteractionBlock,1.0,*currentInteractionBlock);
        *leftInteractionBlock *= h;
        siconos::algebra::prod(*leftInteractionBlock, *rightInteractionBlock,
                               *currentInteractionBlock, false);
        // left = C, right = inv(W).B
      }
    } else if (relationType == siconos::modeling::RelationType::Lagrangian ||
               relationType == siconos::modeling::RelationType::NewtonEuler) {
      // Applying boundary conditions
      std::shared_ptr<siconos::modeling::BoundaryCondition> bc;
      auto d = std::static_pointer_cast<siconos::modeling::SecondOrderDS>(ds);
      if (d->boundaryConditions()) bc = d->boundaryConditions();
      if (bc) {
        for (auto itindex = bc->velocityIndices()->begin();
             itindex != bc->velocityIndices()->end(); ++itindex) {
          // (nslawSize,sizeDS));
          auto coltmp = std::make_shared<siconos::algebra::SiconosVector>(nslawSize);
          coltmp->zero();
          leftInteractionBlock->setCol(*itindex, *coltmp);
        }
      }

      if (osiType == siconos::integrators::IntegratorType::MOREAUJEANBILBAOOSI ||
          std::dynamic_pointer_cast<siconos::modeling::LagrangianLinearDiagonalDS>(ds)) {
        auto work = std::make_shared<siconos::algebra::SimpleMatrix>(*leftInteractionBlock);
        // Get inverse of the iteration matrix
        auto& inv_iteration_matrix = *getOSIMatrix(osi, ds);
        // work = HW (remind that W contains the inverse of the iteration matrix)
        siconos::algebra::axpy_prod(*leftInteractionBlock, inv_iteration_matrix, *work, true);
        leftInteractionBlock->trans();
        siconos::algebra::prod(*work, *leftInteractionBlock, *currentInteractionBlock, false);
        if (relationSubType == siconos::modeling::RelationSubType::CompliantLinearTIR) {
          if (osiType == siconos::integrators::IntegratorType::MOREAUJEANOSI) {
            *currentInteractionBlock *=
                (static_cast<siconos::integrators::MoreauJeanOSI&>(osi)).theta();
            *currentInteractionBlock +=
                *std::static_pointer_cast<siconos::modeling::LagrangianCompliantLinearTIR>(
                     inter->relation())
                     ->D() /
                simulation()->timeStep();
          }
        }
      } else {
        DEBUG_PRINT("leftInteractionBlock after application of boundary conditions\n");
        DEBUG_EXPR(leftInteractionBlock->display(););
        // (inter1 == inter2)
        auto work = std::make_shared<siconos::algebra::SimpleMatrix>(*leftInteractionBlock);
        work->trans();
        auto centralInteractionBlock = getOSIMatrix(osi, ds);
        DEBUG_EXPR_WE(std::cout << std::boolalpha
                                << " centralInteractionBlock->isFactorized() = "
                                << centralInteractionBlock->isFactorized() << std::endl;);
        centralInteractionBlock->Solve(*work);
        //*currentInteractionBlock +=  *leftInteractionBlock ** work;
        DEBUG_EXPR(work->display(););
        siconos::algebra::prod(*leftInteractionBlock, *work, *currentInteractionBlock, false);
        //      gemm(CblasNoTrans,CblasNoTrans,1.0,*leftInteractionBlock,*work,1.0,*currentInteractionBlock);
        //*currentInteractionBlock *=h;
        DEBUG_EXPR(currentInteractionBlock->display(););
        // assert(currentInteractionBlock->checkSymmetry(1e-10));
        if (relationSubType == siconos::modeling::RelationSubType::CompliantLinearTIR) {
          if (osiType == siconos::integrators::IntegratorType::MOREAUJEANOSI) {
            *currentInteractionBlock *=
                (static_cast<siconos::integrators::MoreauJeanOSI&>(osi)).theta();
            *currentInteractionBlock +=
                *std::static_pointer_cast<siconos::modeling::LagrangianCompliantLinearTIR>(
                     inter->relation())
                     ->D() /
                simulation()->timeStep();
          }
        }
      }
    } else
      THROW_EXCEPTION(
          "siconos::nonsmooth_formulations::LinearOSNS::computeDiagonalInteractionBlock not "
          "yet "
          "implemented for relation of type " +
          std::to_string(
              static_cast<std::underlying_type<siconos::modeling::RelationType>::type>(
                  relationType)));
    // Set pos for next loop.
    pos = pos2;
  }
  DEBUG_END(
      "siconos::nonsmooth_formulations::LinearOSNS::computeDiagonalInteractionBlock(const "
      "siconos::graphs::InteractionsGraph::VDescriptor& vd) ends \n");
}

void siconos::nonsmooth_formulations::LinearOSNS::computeInteractionBlock(
    const siconos::graphs::InteractionsGraph::EDescriptor& ed) {
  DEBUG_BEGIN(
      "siconos::nonsmooth_formulations::LinearOSNS::computeInteractionBlock(const "
      "siconos::graphs::InteractionsGraph::EDescriptor& ed)\n");

  // Computes matrix _interactionBlocks[inter1][inter2] (and allocates memory if
  // necessary) if inter1 and inter2 have common DynamicalSystem.  How
  // _interactionBlocks are computed depends explicitely on the type of
  // Relation of each Interaction.

  // Warning: we suppose that at this point, all non linear
  // operators (G for lagrangian relation for example) have been
  // computed through plug-in mechanism.

  // Get dimension of the NonSmoothLaw (ie dim of the interactionBlock)
  auto indexSet = simulation()->indexSet(indexSetLevel());
  auto ds = indexSet->bundle(ed);
  auto inter1 = indexSet->bundle(indexSet->source(ed));
  auto inter2 = indexSet->bundle(indexSet->target(ed));
  // Once again we assume that inter1 and inter2 have the same osi ...
  // auto Osi = indexSet->properties(indexSet->source(ed)).osi;
  auto& DSG0 = *simulation()->nonSmoothDynamicalSystem()->dynamicalSystems();
  auto& osi = *DSG0.properties(DSG0.descriptor(ds)).osi;
  auto osiType = osi.getType();

  // For the edge 'ds', we need to find relative position of this ds
  // in inter1 and inter2 relation matrices (--> pos1 and pos2 below)
  // - find if ds is source or target in inter_i
  siconos::graphs::InteractionsGraph::VDescriptor vertex_inter;
  // - get the corresponding position
  unsigned int pos1, pos2;
  // source of inter1 :
  vertex_inter = indexSet->source(ed);
  std::shared_ptr<siconos::modeling::DynamicalSystem> tmpds =
      indexSet->properties(vertex_inter).source;
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

  unsigned int index1 = indexSet->index(indexSet->source(ed));
  unsigned int index2 = indexSet->index(indexSet->target(ed));
  unsigned int nslawSize1 = inter1->nonSmoothLaw()->size();
  unsigned int nslawSize2 = inter2->nonSmoothLaw()->size();

  std::shared_ptr<siconos::algebra::SiconosMatrix> currentInteractionBlock;

  assert(index1 != index2);

  if (index2 > index1)  // upper block
  {
    assert(indexSet->properties(ed).upper_block->size(0) == nslawSize1);
    assert(indexSet->properties(ed).upper_block->size(1) == nslawSize2);

    currentInteractionBlock = indexSet->properties(ed).upper_block;
  } else  // lower block
  {
    assert(indexSet->properties(ed).lower_block->size(0) == nslawSize1);
    assert(indexSet->properties(ed).lower_block->size(1) == nslawSize2);

    currentInteractionBlock = indexSet->properties(ed).lower_block;
  }

  std::shared_ptr<siconos::algebra::SiconosMatrix> leftInteractionBlock, rightInteractionBlock;

  auto h = simulation()->currentTimeStep();

  // General form of the interactionBlock is : interactionBlock =
  // a*extraInteractionBlock + b * leftInteractionBlock * centralInteractionBlocks
  // * rightInteractionBlock a and b are scalars, centralInteractionBlocks a
  // matrix depending on the integrator (and on the DS), the
  // simulation type ...  left, right and extra depend on the relation
  // type and the non smooth law.
  auto relationType1 = inter1->relation()->getType();
  auto relationType2 = inter2->relation()->getType();

  // ==== First Order Relations - Specific treatment for diagonal
  // _interactionBlocks ===
  assert(inter1 != inter2);
  //  currentInteractionBlock->zero();

  // loop over the common DS
  auto sizeDS = ds->dimension();

  // get _interactionBlocks corresponding to the current DS
  // These _interactionBlocks depends on the relation type.
  leftInteractionBlock = inter1->getLeftInteractionBlockForDS(pos1, nslawSize1, sizeDS);

  // Computing depends on relation type -> move this in Interaction method?
  if (relationType1 == siconos::modeling::RelationType::FirstOrder &&
      relationType2 == siconos::modeling::RelationType::FirstOrder) {
    rightInteractionBlock = inter2->getRightInteractionBlockForDS(pos2, sizeDS, nslawSize2);
    // centralInteractionBlock contains a lu-factorized matrix and we solve
    // centralInteractionBlock * X = rightInteractionBlock with PLU
    auto centralInteractionBlock = getOSIMatrix(osi, ds);
    centralInteractionBlock->Solve(*rightInteractionBlock);

    //      integration of r with theta method removed
    //      *currentInteractionBlock += h *Theta[*itDS]* *leftInteractionBlock *
    //      (*rightInteractionBlock); //left = C, right = W.B
    // gemm(h,*leftInteractionBlock,*rightInteractionBlock,1.0,*currentInteractionBlock);
    *leftInteractionBlock *= h;

    siconos::algebra::prod(*leftInteractionBlock, *rightInteractionBlock,
                           *currentInteractionBlock, false);
    // left = C, right = inv(W).B
  } else if (relationType1 == siconos::modeling::RelationType::Lagrangian ||
             relationType2 == siconos::modeling::RelationType::Lagrangian ||
             relationType1 == siconos::modeling::RelationType::NewtonEuler ||
             relationType2 == siconos::modeling::RelationType::NewtonEuler) {
    // Applying boundary conditions
    std::shared_ptr<siconos::modeling::BoundaryCondition> bc;
    auto d = std::static_pointer_cast<siconos::modeling::SecondOrderDS>(ds);
    if (d->boundaryConditions()) bc = d->boundaryConditions();
    if (bc) {
      for (std::vector<unsigned int>::iterator itindex = bc->velocityIndices()->begin();
           itindex != bc->velocityIndices()->end(); ++itindex) {
        // (nslawSize,sizeDS));
        std::shared_ptr<siconos::algebra::SiconosVector> coltmp =
            std::make_shared<siconos::algebra::SiconosVector>(nslawSize1);
        coltmp->zero();
        leftInteractionBlock->setCol(*itindex, *coltmp);
      }
    }

    auto dsType = siconos::types::type_value(*ds);

    if (osiType == siconos::integrators::IntegratorType::MOREAUJEANBILBAOOSI ||
        dsType == siconos::modeling::Type::LagrangianLinearDiagonalDS) {
      // Rightinteractionblock used first as buffer to save left * W-1
      rightInteractionBlock =
          std::make_shared<siconos::algebra::SimpleMatrix>(nslawSize2, sizeDS);
      // auto work(new SimpleMatrix(*leftInteractionBlock));
      //  Get inverse of the iteration matrix
      auto& inv_iteration_matrix = *getOSIMatrix(osi, ds);
      // remind that W contains the inverse of the iteration matrix
      siconos::algebra::axpy_prod(*leftInteractionBlock, inv_iteration_matrix,
                                  *rightInteractionBlock, true);
      // Then save block corresponding to the 'right' interaction into leftInteractionBlock
      leftInteractionBlock = inter2->getLeftInteractionBlockForDS(pos2, nslawSize1, sizeDS);
      leftInteractionBlock->trans();
      // and compute LW-1R == rightInteractionBlock * leftInteractionBlock into
      // currentInteractionBlock
      siconos::algebra::prod(*rightInteractionBlock, *leftInteractionBlock,
                             *currentInteractionBlock, false);
    } else {
      // inter1 != inter2
      rightInteractionBlock = inter2->getLeftInteractionBlockForDS(pos2, nslawSize2, sizeDS);
      rightInteractionBlock->trans();
      // Warning: we use getLeft for Right interactionBlock
      // because right = transpose(left) and because of
      // size checking inside the getBlock function, a
      // getRight call will fail.
      auto centralInteractionBlock = getOSIMatrix(osi, ds);
      centralInteractionBlock->Solve(*rightInteractionBlock);
      //*currentInteractionBlock +=  *leftInteractionBlock ** work;
      siconos::algebra::prod(*leftInteractionBlock, *rightInteractionBlock,
                             *currentInteractionBlock, false);
    }
  } else
    THROW_EXCEPTION(
        "siconos::nonsmooth_formulations::LinearOSNS::computeInteractionBlock not yet "
        "implemented for "
        "relation of type " +
        std::to_string(
            static_cast<std::underlying_type<siconos::modeling::RelationType>::type>(
                relationType1)));
  DEBUG_END(
      "siconos::nonsmooth_formulations::LinearOSNS::computeInteractionBlock(const "
      "siconos::graphs::InteractionsGraph::EDescriptor& ed)\n");
}

void siconos::nonsmooth_formulations::LinearOSNS::computeqBlock(
    siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter, unsigned int pos) {
  DEBUG_BEGIN(
      "siconos::nonsmooth_formulations::LinearOSNS::computeqBlock(std::shared_ptr<siconos::"
      "modeling::"
      "Interaction> inter, unsigned int pos)\n");
  auto indexSet = simulation()->indexSet(indexSetLevel());

  auto& osi1 = *indexSet->properties(vertex_inter).osi1;
  // auto& osi2 = *indexSet->properties(vertex_inter).osi2;

  auto inter = indexSet->bundle(vertex_inter);
  auto sizeY = inter->nonSmoothLaw()->size();

  osi1.computeFreeOutput(vertex_inter, this);
  auto& osnsp_rhs = osi1.osnsp_rhs(vertex_inter, *indexSet);
  setBlock(osnsp_rhs, _q, sizeY, 0, pos);

  DEBUG_EXPR(_q->display());
  DEBUG_END(
      "siconos::nonsmooth_formulations::LinearOSNS::computeqBlock(std::shared_ptr<siconos::"
      "modeling::"
      "Interaction> inter, unsigned int pos)\n");
}

void siconos::nonsmooth_formulations::LinearOSNS::computeq(double time) {
  DEBUG_BEGIN("void siconos::nonsmooth_formulations::LinearOSNS::computeq(double time)\n");
  if (_q->size() != _sizeOutput) _q->resize(_sizeOutput);
  _q->zero();

  // === Get index set from Simulation ===
  auto indexSet = simulation()->indexSet(indexSetLevel());
  // === Loop through "active" Interactions (ie present in
  // indexSets[level]) ===

  unsigned int pos = 0;
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  for (std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui) {
    // Compute q, this depends on the type of non smooth problem, on
    // the relation type and on the non smooth law
    pos = indexSet->properties(*ui).absolute_position;
    computeqBlock(*ui, pos);  // free output is saved in y
  }
  DEBUG_END("void siconos::nonsmooth_formulations::LinearOSNS::computeq(double time)\n");
}

void siconos::nonsmooth_formulations::LinearOSNS::computeM() {
  if (_assemblyType == LinearOSNSAssemblyType::REDUCED_BLOCK) {
    auto& indexSet = *simulation()->indexSet(indexSetLevel());

    // Computes new _interactionBlocks if required
    updateInteractionBlocks();

    //    _M->fill(indexSet);
    _M->fillM(indexSet, !_hasBeenUpdated);
  } else if (_assemblyType == LinearOSNSAssemblyType::REDUCED_DIRECT) {
    auto& indexSet = *simulation()->indexSet(indexSetLevel());
    auto& DSG0 = *simulation()->nonSmoothDynamicalSystem()->dynamicalSystems();
#ifdef WITH_TIMER
    std::chrono::time_point<std::chrono::system_clock> start, end, end_old;
    start = std::chrono::system_clock::now();
#endif
    // fill _Winverse
    _W_inverse->fillWinverse(DSG0);
#ifdef WITH_TIMER
    end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "\nLinearOSNS: fill W inverse " << elapsed << " ms" << std::endl;
#endif
    // fill H
    _H->fillHtrans(DSG0, indexSet);
#ifdef WITH_TIMER
    end_old = end;
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - end_old).count();
    std::cout << "LinearOSNS: FillH " << elapsed << " ms" << std::endl;
#endif
    // ComputeM
    _M->computeM(_W_inverse->numericsMatrix(), _H->numericsMatrix());
#ifdef WITH_TIMER
    end_old = end;
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - end_old).count();
    std::cout << "LinearOSNS: _computeM " << elapsed << " ms" << std::endl;
#endif
  } else
    THROW_EXCEPTION(
        "siconos::nonsmooth_formulations::LinearOSNS::computeM unknown _assemblyTYPE");

  DEBUG_EXPR(_M->display(););
  // NumericsMatrix *   M_NM = _M->numericsMatrix().get();
  // if (M_NM )
  //   NM_display(M_NM);

  // getchar();
}

bool siconos::nonsmooth_formulations::LinearOSNS::preCompute(double time) {
  DEBUG_BEGIN("bool siconos::nonsmooth_formulations::LinearOSNS::preCompute(double time)\n");
  // This function is used to prepare data for the
  // LinearComplementarityProblem

  // - computation of M and q
  // - set _sizeOutput
  // - check dim. for _z,_w

  // If the topology is time-invariant, only q needs to be computed at
  // each time step.  M, _sizeOutput have been computed in initialize
  // and are uptodate.

  // Get topology
  auto topology = simulation()->nonSmoothDynamicalSystem()->topology();
  bool isLinear = simulation()->nonSmoothDynamicalSystem()->isLinear();
  // int elapsed =0;
  //    std::cout << "!b || !isLinear :"  << boolalpha <<  (!b || !isLinear) <<  std::endl;

  // nothing to do
  if (indexSetLevel() == siconos::internal::LEVELMAX) {
    DEBUG_END("bool siconos::nonsmooth_formulations::LinearOSNS::preCompute(double time)\n");
    return false;
  }

  auto& indexSet = *simulation()->indexSet(indexSetLevel());

  if (indexSet.size() == 0) {
    DEBUG_END("bool siconos::nonsmooth_formulations::LinearOSNS::preCompute(double time)\n");
    return false;
  }
#ifdef WITH_TIMER
  std::chrono::time_point<std::chrono::system_clock> start, end, end_old;
  start = std::chrono::system_clock::now();
#endif
  if (!_hasBeenUpdated || !isLinear) {
    computeM();
#ifdef WITH_TIMER
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "\nLinearOSNS: ComputeM " << elapsed << " ms" << std::endl;
#endif
    //      updateOSNSMatrix();
    _sizeOutput = _M->size();

    // Checks z and _w sizes and reset if necessary

    if (_z->size() != _sizeOutput) {
      _z->resize(_sizeOutput, false);
      _z->zero();
    }

    if (_w->size() != _sizeOutput) {
      _w->resize(_sizeOutput);
      _w->zero();
    }

    // Reset _w and _z with previous values of y and lambda
    // VA 31/05/2021. Values was before the one of the Newton Loop
    // (i.e. values saved in yOutputOld and lambdaOld of the interaction)
    // should be better but needs memory comsumption

    // Note : sizeOuput can be unchanged, but positions may have changed. (??)
    if (_keepLambdaAndYState) {
      siconos::graphs::InteractionsGraph::VIterator ui, uiend;
      for (std::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui) {
        auto& inter = *indexSet.bundle(*ui);
        // Get the position of inter-interactionBlock in the vector w
        // or z
        unsigned int pos = indexSet.properties(*ui).absolute_position;
        // VA 30/08/2021  : Warning. the values of y_k and lambda_k that are stored in Memory
        // may be undefined at the first time step.
        const auto& yOutput_k = inter.y_k(inputOutputLevel());
        const auto& lambda_k = inter.lambda_k(inputOutputLevel());

        if (_sizeOutput >= yOutput_k.size() + pos) {
          siconos::algebra::setBlock(yOutput_k, _w, yOutput_k.size(), 0, pos);
          siconos::algebra::setBlock(lambda_k, _z, lambda_k.size(), 0, pos);
        }
      }
    } else {
      _w->zero();
      _z->zero();
    }
  }
  // else
  // nothing to do (IsLinear and not changed)
#ifdef WITH_TIMER
  end_old = end;
  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - end_old).count();
  std::cout << "LinearOSNS: init w and z " << elapsed << " ms" << std::endl;
#endif
  // Computes q of LinearOSNS
  computeq(time);
#ifdef WITH_TIMER
  end_old = end;
  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - end_old).count();
  std::cout << "LinearOSNS: compute q " << elapsed << " ms" << std::endl;
#endif
  DEBUG_END("bool siconos::nonsmooth_formulations::LinearOSNS::preCompute(double time)\n");
  return true;
}

void siconos::nonsmooth_formulations::LinearOSNS::postCompute() {
  DEBUG_BEGIN("void siconos::nonsmooth_formulations::LinearOSNS::postCompute()\n");
  // This function is used to set y/lambda values using output from
  // lcp_driver (w,z).  Only Interactions (ie Interactions) of
  // indexSet(leveMin) are concerned.

  // === Get index set from Topology ===
  auto& indexSet = *simulation()->indexSet(indexSetLevel());

  // y and lambda vectors
  std::shared_ptr<siconos::algebra::SiconosVector> lambda;
  std::shared_ptr<siconos::algebra::SiconosVector> y;

  // === Loop through "active" Interactions (ie present in
  // indexSets[1]) ===

  unsigned int pos = 0;

  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  for (std::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui) {
    auto& inter = *indexSet.bundle(*ui);
    // Get the  position of inter-interactionBlock in the vector w
    // or z
    pos = indexSet.properties(*ui).absolute_position;

    // Get Y and Lambda for the current Interaction
    y = inter.y(inputOutputLevel());
    lambda = inter.lambda(inputOutputLevel());
    // Copy _w/_z values, starting from index pos into y/lambda.

    // siconos::algebra::setBlock(*_w, y, y->size(), pos, 0);// Warning: yEquivalent is
    //  saved in y !!
    siconos::algebra::setBlock(*_z, lambda, lambda->size(), pos, 0);
    DEBUG_EXPR(lambda->display(););

  }
  DEBUG_END("void siconos::nonsmooth_formulations::LinearOSNS::postCompute()\n");
}

void siconos::nonsmooth_formulations::LinearOSNS::display() const {
  std::cout << "==========================" << std::endl;
  std::cout << "this : " << this << std::endl;
  std::cout << "_M  ";
  if (_M)
    _M->display();
  else
    std::cout << "-> nullptr" << std::endl;
  std::cout << std::endl << "q : ";
  if (_q)
    _q->display();
  else
    std::cout << "-> nullptr" << std::endl;
  std::cout << std::endl;
  std::cout << "w : ";
  if (_w)
    _w->display();
  else
    std::cout << "-> nullptr" << std::endl;
  std::cout << std::endl << "z : ";
  if (_z)
    _z->display();
  else
    std::cout << "-> nullptr" << std::endl;
  std::cout << std::endl;
  std::cout << "The linearOSNSP works on the index set of level  " << _indexSetLevel
            << std::endl;
  std::cout << "==========================" << std::endl;
}
