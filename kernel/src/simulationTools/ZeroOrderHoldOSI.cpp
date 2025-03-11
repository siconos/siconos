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

#include "ZeroOrderHoldOSI.hpp"

#include "BlockVector.hpp"
#include "EventsManager.hpp"
#include "FirstOrderLinearDS.hpp"
#include "FirstOrderLinearR.hpp"
#include "FirstOrderLinearTIR.hpp"
#include "FirstOrderR.hpp"
#include "Interaction.hpp"
#include "MatrixIntegrator.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "NewtonImpactNSL.hpp"
#include "OneStepNSProblem.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixVectorOp.hpp"
#include "SiconosVector.hpp"
#include "Simulation.hpp"
#include "Topology.hpp"

// #define DEBUG_WHERE_MESSAGES

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES

#include "siconos_debug.h"

// --- constructor from a minimum set of data ---
siconos::integrators::ZeroOrderHoldOSI::ZeroOrderHoldOSI()
    : OneStepIntegrator{IntegratorType::ZOHOSI, 1, 0, 0, 0, 0} {}

void siconos::integrators::ZeroOrderHoldOSI::initializeWorkVectorsForDS(
    double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  DEBUG_BEGIN(
      "void siconos::integrators::ZeroOrderHoldOSI::initializeWorkVectorsForDS( double t, "
      "std::shared_ptr<siconos::modeling::DynamicalSystem> ds)\n");
  // Get work buffers from the graph
  auto& ds_work_vectors = *_initializeDSWorkVectors(ds);
  ds_work_vectors.resize(siconos::integrators::ZeroOrderHoldOSI::WORK_LENGTH);

  auto& DSG0 = *_dynamicalSystemsGraph;
  auto& IG0 = *_simulation->nonSmoothDynamicalSystem()->topology()->indexSet0();

  if (not std::dynamic_pointer_cast<siconos::modeling::FirstOrderLinearDS>(ds))
    THROW_EXCEPTION(
        "siconos::integrators::ZeroOrderHoldOSI::initialize - the DynamicalSystem does not "
        "have the right type");

  unsigned int indxIter = 0;
  siconos::graphs::DynamicalSystemsGraph::AVIterator avi, aviend;
  auto dsgVD = DSG0.descriptor(ds);
  if (!DSG0.Ad.hasKey(dsgVD)) {
    DSG0.Ad[dsgVD] = std::make_shared<siconos::simulation::MatrixIntegrator>(
        *ds, *_simulation->nonSmoothDynamicalSystem(),
        _simulation->eventsManager()->timeDiscretisation());
    if (DSG0.Ad.at(dsgVD)->isConst()) DSG0.Ad.at(dsgVD)->integrate();
  } else
    THROW_EXCEPTION(
        "siconos::integrators::ZeroOrderHoldOSI::initialize - Ad MatrixIntegrator is already "
        "initialized for ds the DS");

  if ((static_cast<const siconos::modeling::FirstOrderLinearDS&>(*ds)).hasbVector()) {
    siconos::algebra::SiconosMatrix E{ds->dimension(), ds->dimension()};
    E.setIdentity();
    DSG0.AdInt.insert(dsgVD, std::make_shared<siconos::simulation::MatrixIntegrator>(
                                 *ds, *_simulation->nonSmoothDynamicalSystem(),
                                 _simulation->eventsManager()->timeDiscretisation(), E));
    if (DSG0.AdInt.at(dsgVD)->isConst()) DSG0.AdInt.at(dsgVD)->integrate();
  }

  // init extra term, usually to add control terms
  if (_extraAdditionalTerms)
    _extraAdditionalTerms->init(DSG0, *_simulation->nonSmoothDynamicalSystem(),
                                _simulation->eventsManager()->timeDiscretisation());

  // Now we search for an Interaction dedicated to control
  for (std::tie(avi, aviend) = DSG0.adjacent_vertices(dsgVD); avi != aviend; ++avi) {
    siconos::graphs::DynamicalSystemsGraph::EDescriptor ed1, ed2;
    std::tie(ed1, ed2) = DSG0.edges(dsgVD, *avi);

    if (IG0.properties(IG0.descriptor(DSG0.bundle(ed1))).forControl) {
      auto& inter = *DSG0.bundle(ed1);
      auto& rel = *inter.relation();
      if (rel.getType() != siconos::modeling::RelationType::FirstOrder)
        THROW_EXCEPTION(
            "siconos::integrators::ZeroOrderHoldOSI::initialize - the Integrator can only "
            "deal with FirstOrder Relation");
      auto& relR = static_cast<siconos::modeling::FirstOrderR&>(rel);

      if (indxIter == 0) {
        indxIter++;
        if (relR.hasConstantJacobiangOver_lambda()) {
          DSG0.Bd[dsgVD] = std::make_shared<siconos::simulation::MatrixIntegrator>(
              *ds, *_simulation->nonSmoothDynamicalSystem(),
              _simulation->eventsManager()->timeDiscretisation(), relR.jacobiangOver_lambda());
          if (DSG0.Bd.at(dsgVD)->isConst()) DSG0.Bd.at(dsgVD)->integrate();
        } else {  // user defined function for jacobiangOver_lambda
          THROW_EXCEPTION(
              "siconos::integrators::ZeroOrderHoldOSI::initialize - Case not implemented");
          // DSG0.Bd[dsgVD] = std::make_shared<siconos::simulation::MatrixIntegrator>(
          //     *ds, *_simulation->nonSmoothDynamicalSystem(),
          //     _simulation->eventsManager()->timeDiscretisation(), relR.getPluging(),
          //     inter.dimension());
        }
      } else {
        //        THROW_EXCEPTION("siconos::integrators::ZeroOrderHoldOSI::initialize - DS
        //        linked with more that one iteraction");
        DEBUG_PRINTF("number of iteraction attached to the process : %d\n", indxIter);
      }
    }
  }

  // Get work buffers from the graph
  ds_work_vectors[siconos::integrators::ZeroOrderHoldOSI::FREE] =
      std::make_shared<siconos::algebra::SiconosVector>(ds->dimension());
  ds_work_vectors[siconos::integrators::ZeroOrderHoldOSI::DELTA_X_FOR_RELATION] =
      std::make_shared<siconos::algebra::SiconosVector>(ds->dimension());
  DEBUG_END(
      "void siconos::integrators::ZeroOrderHoldOSI::initializeWorkVectorsForDS( double t, "
      "std::shared_ptr<siconos::modeling::DynamicalSystem> ds)\n");
}

void siconos::integrators::ZeroOrderHoldOSI::initializeWorkVectorsForInteraction(
    siconos::modeling::Interaction& inter, siconos::graphs::InteractionProperties& interProp,
    siconos::graphs::DynamicalSystemsGraph& DSG) {
  auto ds1 = interProp.source;
  auto ds2 = interProp.target;
  assert(ds1);
  assert(ds2);

  if (!interProp.workVectors) {
    interProp.workVectors =
        std::make_shared<std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>>(
            siconos::integrators::ZeroOrderHoldOSI::WORK_INTERACTION_LENGTH);
  }

  if (!interProp.workBlockVectors) {
    interProp.workBlockVectors =
        std::make_shared<std::vector<std::shared_ptr<siconos::algebra::BlockVector>>>(
            siconos::integrators::ZeroOrderHoldOSI::BLOCK_WORK_LENGTH);
  }

  auto& inter_work = *interProp.workVectors;
  auto& inter_work_block = *interProp.workBlockVectors;

  auto& relation = *inter.relation();
  auto relationType = relation.getType();
  auto relationSubType = inter.relation()->getSubType();

  if (!inter_work[siconos::integrators::ZeroOrderHoldOSI::OSNSP_RHS])
    inter_work[siconos::integrators::ZeroOrderHoldOSI::OSNSP_RHS] =
        std::make_shared<siconos::algebra::SiconosVector>(inter.dimension());

  // Check if interations levels (i.e. y and lambda sizes) are compliant with the current osi.
  _check_and_update_interaction_levels(inter);
  // Initialize/allocate memory buffers in interaction.
  inter.initializeMemory(_steps);

  /* allocate and set work vectors for the osi */
  if (!(checkOSI(DSG.descriptor(ds1)) && checkOSI(DSG.descriptor(ds2)))) {
    THROW_EXCEPTION(
        "siconos::integrators::ZeroOrderHoldOSI::initializeWorkVectorsForInteraction. The "
        "implementation is not correct for two different OSI for one interaction");
  }

  auto xfree = siconos::integrators::ZeroOrderHoldOSI::xfree;

  auto& workVds1 = *DSG.properties(DSG.descriptor(ds1)).workVectors;
  if (relationType == siconos::modeling::RelationType::FirstOrder) {
    if (relationSubType == siconos::modeling::RelationSubType::NonLinearR ||
        relationSubType == siconos::modeling::RelationSubType::Type2R) {
      inter_work[siconos::integrators::ZeroOrderHoldOSI::H_ALPHA] =
          std::make_shared<siconos::algebra::SiconosVector>(inter.dimension());
    }
    inter_work_block[xfree] = std::make_shared<siconos::algebra::BlockVector>();
    inter_work_block[xfree]->insertPtr(workVds1[siconos::integrators::ZeroOrderHoldOSI::FREE]);
  }

  if (ds1 != ds2) {
    auto& workVds2 = *DSG.properties(DSG.descriptor(ds2)).workVectors;
    if (relationType == siconos::modeling::RelationType::Lagrangian) {
      inter_work_block[siconos::integrators::ZeroOrderHoldOSI::xfree]->insertPtr(
          workVds2[siconos::integrators::ZeroOrderHoldOSI::FREE]);
    }
  }

  if (!inter_work_block[siconos::integrators::ZeroOrderHoldOSI::DELTA_X]) {
    inter_work_block[siconos::integrators::ZeroOrderHoldOSI::DELTA_X] =
        std::make_shared<siconos::algebra::BlockVector>();
    inter_work_block[siconos::integrators::ZeroOrderHoldOSI::DELTA_X]->insertPtr(
        workVds1[siconos::integrators::ZeroOrderHoldOSI::DELTA_X_FOR_RELATION]);
  } else
    inter_work_block[siconos::integrators::ZeroOrderHoldOSI::DELTA_X]->setVectorPtr(
        0, workVds1[siconos::integrators::ZeroOrderHoldOSI::DELTA_X_FOR_RELATION]);
}

double siconos::integrators::ZeroOrderHoldOSI::computeResidu() {
  DEBUG_BEGIN("double siconos::integrators::ZeroOrderHoldOSI::computeResidu()\n");
  // This function is used to compute the residu for each "MoreauJeanOSI-discretized" dynamical
  // system. It then computes the norm of each of them and finally return the maximum value for
  // those norms.
  //
  // The state values used are those saved in the DS, ie the last computed ones.
  //  $\mathcal R(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) - h r$
  //  $\mathcal R_{free}(x) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) $

  // Operators computed at told have index i, and (i+1) at t.

  // Iteration through the set of Dynamical Systems.
  //
  double maxResidu = 0;
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;

  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);

    // 1 - First Order Linear Systems
    if (std::dynamic_pointer_cast<siconos::modeling::FirstOrderLinearDS>(ds)) {
      // No residu with ZOH ...
    } else
      THROW_EXCEPTION(
          "siconos::integrators::ZeroOrderHoldOSI::computeResidu - Only implemented for first "
          "order linear systems");
  }
  DEBUG_END("double siconos::integrators::ZeroOrderHoldOSI::computeResidu()\n");
  return maxResidu;
}

void siconos::integrators::ZeroOrderHoldOSI::computeFreeState() {
  DEBUG_BEGIN("void siconos::integrators::ZeroOrderHoldOSI::computeFreeState()\n");
  // This function computes "free" states of the DS belonging to this Integrator.
  // "Free" means without taking non-smooth effects into account.

  // Operators computed at told have index i, and (i+1) at t.
  double t = _simulation->nextTime();         // End of the time step
  double told = _simulation->startingTime();  // Beginning of the time step
  double h = t - told;                        // time step length

  auto& DSG0 = *_dynamicalSystemsGraph;

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  DEBUG_EXPR(display(););
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    DEBUG_EXPR(ds->display(););
    auto dsgVD = DSG0.descriptor(ds);
    auto& ds_work_vectors = *DSG0.properties(dsgVD).workVectors;
    if (auto d = std::dynamic_pointer_cast<siconos::modeling::FirstOrderLinearDS>(ds)) {
      // Check whether we have to recompute things
      if (!DSG0.Ad.at(dsgVD)->isConst()) DSG0.Ad.at(dsgVD)->integrate();
      if (d->hasbVector() && !DSG0.AdInt.at(dsgVD)->isConst())
        DSG0.AdInt.at(dsgVD)->integrate();

      auto& xfree = *ds_work_vectors[siconos::integrators::ZeroOrderHoldOSI::FREE];
      DEBUG_EXPR(siconos::algebra::print(xfree););

      xfree = DSG0.Ad.at(dsgVD)->mat() * *d->x();  // xfree = Ad*xold
      DEBUG_EXPR(siconos::algebra::print(xfree););
      if (d->hasbVector()) {
        assert(DSG0.AdInt.hasKey(dsgVD));
        xfree += DSG0.AdInt.at(dsgVD)->mat() * d->bVector();  // xfree += AdInt*b
        DEBUG_EXPR(siconos::algebra::print(xfree););
      }

      // add extra term, possible control terms
      if (_extraAdditionalTerms) {
        DEBUG_PRINT("add extra additional terms\n");
        _extraAdditionalTerms->addSmoothTerms(DSG0, dsgVD, h, xfree);
      }
      DEBUG_EXPR(siconos::algebra::print(xfree););
    } else
      THROW_EXCEPTION(
          "siconos::integrators::ZeroOrderHoldOSI::computeFreeState - Only implemented for "
          "first "
          "order linear systems");
  }
  DEBUG_END("void siconos::integrators::ZeroOrderHoldOSI::computeFreeState()\n");
}

void siconos::integrators::ZeroOrderHoldOSI::prepareNewtonIteration(double time) {
  // siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;

  // for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  // {
  //   //if (!checkOSI(dsi)) continue;
  //   //auto ds = _dynamicalSystemsGraph->bundle(*dsi);
  //   //    computeMatrices(time, *ds);
  // }

  if (!_explicitJacobiansOfRelation) {
    _simulation->nonSmoothDynamicalSystem()->computeInteractionJacobians(time);
  }
}

struct siconos::integrators::ZeroOrderHoldOSI::_NSLEffectOnFreeOutput
    : public siconos::modeling::nonsmooth_laws::Visitor {
  using siconos::modeling::nonsmooth_laws::Visitor::visit;

  siconos::nonsmooth_formulations::OneStepNSProblem* _osnsp{nullptr};
  std::shared_ptr<siconos::modeling::Interaction> _inter{nullptr};
  siconos::graphs::InteractionProperties& _interProp;
  _NSLEffectOnFreeOutput(siconos::nonsmooth_formulations::OneStepNSProblem* p,
                         std::shared_ptr<siconos::modeling::Interaction> inter,
                         siconos::graphs::InteractionProperties& interProp)
      : _osnsp(p), _inter(inter), _interProp(interProp) {};

  void visit(const siconos::modeling::NewtonImpactNSL& nslaw) override {
    auto e = nslaw.e();
    siconos::algebra::SiconosVector& osnsp_rhs =
        *(*_interProp.workVectors)[siconos::integrators::ZeroOrderHoldOSI::OSNSP_RHS];
    osnsp_rhs += e * _inter->y_k(_osnsp->inputOutputLevel());
  }

  void visit(const siconos::modeling::NewtonImpactFrictionNSL& nslaw) override {
    double e;
    e = nslaw.en();
    // Only the normal part is multiplied by e
    auto& osnsp_rhs =
        *(*_interProp.workVectors)[siconos::integrators::ZeroOrderHoldOSI::OSNSP_RHS];
    osnsp_rhs(0) += e * _inter->y_k(_osnsp->inputOutputLevel())(0);
  }

  void visit(const siconos::modeling::EqualityConditionNSL& nslaw) override { ; }
  void visit(const siconos::modeling::MixedComplementarityConditionNSL& nslaw) override { ; }
};

void siconos::integrators::ZeroOrderHoldOSI::computeFreeOutput(
    siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter,
    siconos::nonsmooth_formulations::OneStepNSProblem* osnsp) {
  DEBUG_BEGIN("void siconos::integrators::ZeroOrderHoldOSI::computeFreeOutput(...)\n");
  /** \warning: ensures that it can also work with two different osi for two different ds ?
   */
  auto indexSet = osnsp->simulation()->indexSet(osnsp->indexSetLevel());
  auto inter = indexSet->bundle(vertex_inter);
  auto& inter_work = *indexSet->properties(vertex_inter).workVectors;
  auto& inter_work_block = *indexSet->properties(vertex_inter).workBlockVectors;

  // Get relation and non smooth law types
  auto relationType = inter->relation()->getType();
  auto relationSubType = inter->relation()->getSubType();
  auto deltax = inter_work_block[siconos::integrators::ZeroOrderHoldOSI::DELTA_X];

  auto& osnsp_rhs = *(*indexSet->properties(vertex_inter)
                           .workVectors)[siconos::integrators::ZeroOrderHoldOSI::OSNSP_RHS];

  std::shared_ptr<siconos::algebra::BlockVector> Xfree;
  if (relationType == siconos::modeling::RelationType::FirstOrder) {
    Xfree = inter_work_block[siconos::integrators::ZeroOrderHoldOSI::xfree];
  }
  assert(Xfree);

  auto rel = inter->relation();
  assert(inter);
  assert(rel);

  assert(relationType == siconos::modeling::RelationType::FirstOrder);
  auto forel = std::static_pointer_cast<siconos::modeling::FirstOrderR>(rel);

  //  if (!IG0.properties(IG0.descriptor(inter)).forControl) // the integration is not for
  //  control
  {
    if (relationSubType == siconos::modeling::RelationSubType::Type2R) {
      auto lambda = inter->lambda(0);
      assert(lambda);

      if (forel->hasJacobianhOver_lambda()) {
        auto D = rel->jacobianhOver_lambda();
        osnsp_rhs = D * *lambda;
        osnsp_rhs *= -1.0;
      }
      if (forel->hasJacobianhOver_state()) {
        auto C = forel->jacobianhOver_state();
        siconos::algebra::matrixBlockVector_prod(C, *deltax, osnsp_rhs, false);
      }

      if (_useGammaForRelation) {
        THROW_EXCEPTION(
            "siconos::integrators::ZeroOrderHoldOSI::ComputeFreeOutput not yet implemented "
            "with useGammaForRelation() for FirstorderR and Typ2R and H_alpha() "
            "should return the mid-point value");
      }

      auto& hAlpha = *inter_work[siconos::integrators::ZeroOrderHoldOSI::H_ALPHA];
      osnsp_rhs += hAlpha;
    }

    else {
      if (forel->hasJacobianhOver_state()) {
        auto C = forel->jacobianhOver_state();
        assert(Xfree);
        assert(deltax);
        // creates a POINTER link between workX[ds] (xfree) and the
        // corresponding interactionBlock in each Interactionfor each ds of the
        // current Interaction.

        if (_useGammaForRelation) {
          siconos::algebra::matrixBlockVector_prod(C, *deltax, osnsp_rhs, true);
        } else {
          siconos::algebra::matrixBlockVector_prod(C, *Xfree, osnsp_rhs, true);
        }
      }

      if (relationSubType == siconos::modeling::RelationSubType::LinearTIR ||
          relationSubType == siconos::modeling::RelationSubType::LinearR) {
        // In the first order linear case it may be required to add e to q.
        // q = HXfree + e

        if (relationSubType == siconos::modeling::RelationSubType::LinearTIR) {
          auto linrel = std::static_pointer_cast<siconos::modeling::FirstOrderLinearTIR>(rel);
          if (linrel->haseVector()) {
            auto e = linrel->eVector();
            osnsp_rhs += e;
          }
        } else {
          auto linrel = std::static_pointer_cast<siconos::modeling::FirstOrderLinearR>(rel);
          if (linrel->haseVector()) {
            auto e = linrel->eVector();
            osnsp_rhs += e;
          }
        }
      }
    }
  }

  DEBUG_END("void siconos::integrators::ZeroOrderHoldOSI::computeFreeOutput(...)\n");
}
void siconos::integrators::ZeroOrderHoldOSI::integrate(double& tinit, double& tend,
                                                       double& tout, int&) {
  // This function should not be used
  THROW_EXCEPTION("siconos::integrators::ZeroOrderHoldOSI::integrate - should not be used");
}

void siconos::integrators::ZeroOrderHoldOSI::updateState(const unsigned int level) {
  DEBUG_BEGIN(
      "siconos::integrators::ZeroOrderHoldOSI::updateState(const unsigned int level)\n");
  bool useRCC = _simulation->useRelativeConvergenceCriteron();
  if (useRCC) _simulation->setRelativeConvergenceCriterionHeld(true);

  auto& DSG0 = *_dynamicalSystemsGraph;
  siconos::graphs::DynamicalSystemsGraph::OEIterator oei, oeiend;
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    DEBUG_EXPR(ds->display(););
    auto dsgVD = DSG0.descriptor(ds);
    auto& ds_work_vectors = *DSG0.properties(dsgVD).workVectors;
    // 1 - First Order Linear Systems
    if (auto d = std::dynamic_pointer_cast<siconos::modeling::FirstOrderLinearDS>(ds)) {
      auto& x = *d->x();
      // 1 - First Order Linear Time Invariant Systems
      // \Phi is already computed
      x = *ds_work_vectors[siconos::integrators::ZeroOrderHoldOSI::FREE];  // x = xfree =
                                                                           // Phi*xold (+ Bd*u
                                                                           // ) (+  Ld*e)
      if (level != siconos::internal::LEVELMAX) {
        std::shared_ptr<siconos::modeling::Interaction> interC;
        // we have to find the control interaction
        for (std::tie(oei, oeiend) = DSG0.out_edges(dsgVD); oei != oeiend; ++oei) {
          if (DSG0.properties(*oei).forControl) {
            interC = DSG0.bundle(*oei);
            break;
          }
        }
        if (interC) {
          DEBUG_PRINT("A control interaction is found\n");
          auto& Bd = *DSG0.Bd[dsgVD];
          if (!Bd.isConst()) {
            Bd.integrate();
          }
          x += Bd.mat() * *interC->lambda(0);  // x += Bd*\lambda
        }
      }
      DEBUG_EXPR(ds->display(););
    } else
      THROW_EXCEPTION(
          "siconos::integrators::ZeroOrderHoldOSI::updateState - Only implemented for first "
          "order linear DS");
  }
  DEBUG_END("siconos::integrators::ZeroOrderHoldOSI::updateState(const unsigned int level)\n");
}

bool siconos::integrators::ZeroOrderHoldOSI::addInteractionInIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i) {
  assert(i == 1);
  double h = _simulation->timeStep();
  double y = (*inter->y(i - 1))(0);   // for i=1 y(i-1) is the position
  double yDot = (*(inter->y(i)))(0);  // for i=1 y(i) is the velocity
  double gamma = .5;
  DEBUG_PRINTF(
      "siconos::integrators::ZeroOrderHoldOSI::addInteractionInIndexSet yref=%e, yDot=%e, "
      "y_estimated=%e.\n",
      y, yDot, y + gamma * h * yDot);
  y += gamma * h * yDot;
  assert(!std::isnan(y));
  if (y <= 0) {
    DEBUG_PRINT(
        "siconos::integrators::ZeroOrderHoldOSI::addInteractionInIndexSet ACTIVATE.\n");
  }
  return (y <= 0);
}

bool siconos::integrators::ZeroOrderHoldOSI::removeInteractionFromIndexSet(
    std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i) {
  return !(addInteractionInIndexSet(inter, i));
}

const siconos::algebra::SiconosMatrix& siconos::integrators::ZeroOrderHoldOSI::Ad(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) const {
  auto& DSG0 = *_simulation->nonSmoothDynamicalSystem()->topology()->dSG(0);
  auto dsgVD = DSG0.descriptor(ds);
  return DSG0.Ad.at(dsgVD)->mat();
}

const siconos::algebra::SiconosMatrix& siconos::integrators::ZeroOrderHoldOSI::Bd(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) const {
  auto& DSG0 = *_simulation->nonSmoothDynamicalSystem()->topology()->dSG(0);
  auto dsgVD = DSG0.descriptor(ds);
  return DSG0.Bd.at(dsgVD)->mat();
}

void siconos::integrators::ZeroOrderHoldOSI::display() const {
  OneStepIntegrator::display();

  std::cout << "====== ZOH OSI display ======" << std::endl;
  std::cout << "--------------------------------" << std::endl;
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;

    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    if (_dynamicalSystemsGraph->Ad[*dsi]) {
      std::cout << "--> Phi of dynamical system number: \n";
      siconos::algebra::print(_dynamicalSystemsGraph->Ad[*dsi]->mat());
    }

    if (_dynamicalSystemsGraph->Bd[*dsi]) {
      std::cout << "--> Psi of dynamical system number: \n";
      siconos::algebra::print(_dynamicalSystemsGraph->Bd[*dsi]->mat());
    }
  }
  std::cout << "================================" << std::endl;
}
