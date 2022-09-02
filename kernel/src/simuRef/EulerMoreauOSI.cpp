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
#include "EulerMoreauOSI.hpp"
// #include "SiconosAlgebraProd.hpp" // for prod and subprod
// #include "SiconosAlgebraScal.hpp" // for scal
// //#include "siconos::algebra::SiconosVectorFriends.hpp" // for scal
// #include "Simulation.hpp"
// #include "NonSmoothDynamicalSystem.hpp"
// #include "FirstOrderLinearTIDS.hpp"
// #include "FirstOrderLinearTIR.hpp"
// #include "FirstOrderLinearR.hpp"
// #include "FirstOrderNonLinearR.hpp"
// #include "FirstOrderType2R.hpp"
// #include "FirstOrderType1R.hpp"
// #include "NonSmoothLaw.hpp"
// #include "OneStepNSProblem.hpp"
// #include "BlockVector.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
// #define DEBUG_WHERE_MESSAGES
#include "siconos_debug.h"

// --- constructor with theta parameter value  ---
siconos::integrators::EulerMoreauOSI::EulerMoreauOSI(double theta)
    : OneStepIntegrator(IntegratorType::EULERMOREAUOSI)
{
  _levelMinForOutput = 0;
  _levelMaxForOutput = 0;
  _levelMinForInput = 0;
  _levelMaxForInput = 0;
  _steps = 1;
  _theta = theta;
}

// --- constructor from a set of data ---
siconos::integrators::EulerMoreauOSI::EulerMoreauOSI(double theta, double gamma)
    : EulerMoreauOSI{theta}, _theta{theta}, _gamma{gamma}, _useGamma{true}
{
}

const siconos::algebra::SimpleMatrix siconos::integrators::EulerMoreauOSI::getW(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds)
{
  assert(ds && "siconos::integrators::EulerMoreauOSI::getW(ds): ds == nullptr.");
  assert(_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W &&
         "siconos::integrators::EulerMoreauOSI::getW(ds): W[ds] == nullptr.");
  return *(_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds))
               .W);  // Copy !!
}

std::shared_ptr<siconos::algebra::SimpleMatrix> siconos::integrators::EulerMoreauOSI::W(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds)
{
  assert(ds && "siconos::integrators::EulerMoreauOSI::W(ds): ds == nullptr.");
  return _dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W;
}

const siconos::algebra::SimpleMatrix
siconos::integrators::EulerMoreauOSI::getWBoundaryConditions(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds)
{
  assert(ds &&
         "siconos::integrators::EulerMoreauOSI::getWBoundaryConditions(ds): ds == nullptr.");
  //    return *(WBoundaryConditionsMap[0]);
  assert(_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds))
             .WBoundaryConditions &&
         "siconos::integrators::EulerMoreauOSI::getWBoundaryConditions(ds): "
         "WBoundaryConditions[ds] == nullptr.");
  return *(_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds))
               .WBoundaryConditions);  // Copy !!
}

std::shared_ptr<siconos::algebra::SiconosMatrix>
siconos::integrators::EulerMoreauOSI::WBoundaryConditions(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds)
{
  assert(ds &&
         "siconos::integrators::EulerMoreauOSI::WBoundaryConditions(ds): ds == nullptr.");
  return _dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds))
      .WBoundaryConditions;
}

void siconos::integrators::EulerMoreauOSI::initializeWorkVectorsForDS(
    double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds)
{
  auto& ds_work_vectors = *_initializeDSWorkVectors(ds);
  ds_work_vectors.resize(siconos::integrators::EulerMoreauOSI::WORK_LENGTH);

  // Check dynamical system type
  auto fods = std::dynamic_pointer_cast<siconos::modeling::FirstOrderNonLinearDS>(ds);
  assert(fods);

  // Compute W (iteration matrix)
  initializeIterationMatrixW(t, ds);

  // buffers allocation (into the graph)
  ds_work_vectors[siconos::integrators::EulerMoreauOSI::RESIDU] = std::make_shared<siconos::algebra::SiconosVector>(ds->dimension()));
  ds_work_vectors[siconos::integrators::EulerMoreauOSI::RESIDU_FREE] = std::make_shared<siconos::algebra::SiconosVector>(ds->dimension()));
  ds_work_vectors[siconos::integrators::EulerMoreauOSI::FREE] = std::make_shared<siconos::algebra::SiconosVector>(ds->dimension()));
  ds_work_vectors[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS_FOR_RELATION] = std::make_shared<siconos::algebra::SiconosVector>(ds->dimension()));
  ds_work_vectors[siconos::integrators::EulerMoreauOSI::DELTA_X_FOR_RELATION] = std::make_shared<siconos::algebra::SiconosVector>(ds->dimension()));
  ds_work_vectors[siconos::integrators::EulerMoreauOSI::LOCAL_BUFFER] = std::make_shared<siconos::algebra::SiconosVector>(ds->dimension()));

  // Update dynamical system components (for memory swap).
  fods->computef(t, fods->x());  // Only fold is concerned, for FirstOrderNonLinearDS.
  // Update memory buffers
  ds->swapInMemory();
}

void siconos::integrators::EulerMoreauOSI::initializeWorkVectorsForInteraction(
    siconos::modeling::Interaction& inter, siconos::graphs::InteractionProperties& interProp,
    siconos::graphs::DynamicalSystemsGraph& DSG)
{
  auto ds1 = interProp.source;
  auto ds2 = interProp.target;
  assert(ds1);
  assert(ds2);

  if (!interProp.workVectors) {
    interProp.workVectors = std::make_shared<std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>>;
    interProp.workVectors->resize(
        siconos::integrators::EulerMoreauOSI::WORK_INTERACTION_LENGTH);
  }
  if (!interProp.workMatrices) {
    interProp.workMatrices =
        std::make_shared<std::vector<std::shared_ptr<siconos::algebra::SimpleMatrix>>>();
    interProp.workMatrices->resize(siconos::integrators::EulerMoreauOSI::MAT_WORK_LENGTH);
  }
  if (!interProp.workBlockVectors) {
    interProp.workBlockVectors =
        std::make_shared<std::vector<std::shared_ptr<siconos::algebra::BlockVector>>>();
    interProp.workBlockVectors->resize(
        siconos::integrators::EulerMoreauOSI::BLOCK_WORK_LENGTH);
  }

  auto& inter_work = *interProp.workVectors;
  auto& inter_work_mat = *interProp.workMatrices;
  auto& inter_work_block = *interProp.workBlockVectors;

  Relation& relation = *inter.relation();

  unsigned int sizeY = inter.dimension();
  unsigned int sizeOfDS = inter.getSizeOfDS();

  auto relationType = relation.getType();
  auto relationSubType = relation.getSubType();
  inter_work[siconos::integrators::EulerMoreauOSI::OSNSP_RHS] = std::make_shared<siconos::algebra::SiconosVector>(inter.dimension()));

  // Check if interations levels (i.e. y and lambda sizes) are compliant with the current osi.
  _check_and_update_interaction_levels(inter);
  // Initialize/allocate memory buffers in interaction.
  inter.initializeMemory(_steps);

  if (checkOSI(DSG.descriptor(ds1))) {
    DEBUG_PRINTF("ds1->number() %i is taken in to account\n", ds1->number());
    assert(DSG.properties(DSG.descriptor(ds1)).workVectors);
    auto& workVds1 = *DSG.properties(DSG.descriptor(ds1)).workVectors;

    if (relationType == FirstOrder) {
      inter_work[siconos::integrators::EulerMoreauOSI::VEC_X] = std::make_shared<siconos::algebra::SiconosVector>(sizeOfDS));
      inter_work[siconos::integrators::EulerMoreauOSI::VEC_RESIDU_Y] = std::make_shared<siconos::algebra::SiconosVector>(sizeY));
      inter_work[siconos::integrators::EulerMoreauOSI::H_ALPHA] = std::make_shared<siconos::algebra::SiconosVector>(sizeY));
      inter_work[siconos::integrators::EulerMoreauOSI::YOLD] = std::make_shared<siconos::algebra::SiconosVector>(sizeY));
      inter_work[siconos::integrators::EulerMoreauOSI::LAMBDAOLD] = std::make_shared<siconos::algebra::SiconosVector>(sizeY));

      if (relationSubType == NonLinearR || relationSubType == Type2R) {
        // inter_work[siconos::integrators::EulerMoreauOSI::G_ALPHA] =
        // std::make_shared<siconos::algebra::SiconosVector>(sizeOfDS));
        inter_work_block[siconos::integrators::EulerMoreauOSI::G_ALPHA].reset(
            new BlockVector(1, sizeOfDS));
        inter_work[siconos::integrators::EulerMoreauOSI::VEC_RESIDU_R] = std::make_shared<siconos::algebra::SiconosVector>(sizeOfDS));
        inter_work_mat[siconos::integrators::EulerMoreauOSI::MAT_KHAT] = std::make_shared<siconos::algebra::SimpleMatrix>(sizeOfDS, sizeY));
        inter_work_mat[siconos::integrators::EulerMoreauOSI::MAT_KTILDE] = std::make_shared<siconos::algebra::SimpleMatrix>(sizeOfDS, sizeY));
      }

      if (!inter_work_block[siconos::integrators::EulerMoreauOSI::XFREE]) {
        inter_work_block[siconos::integrators::EulerMoreauOSI::XFREE] = std::make_shared<siconos::algebra::BlockVector>();
        inter_work_block[siconos::integrators::EulerMoreauOSI::XFREE]->insertPtr(
            workVds1[siconos::integrators::EulerMoreauOSI::FREE]);
      }
      else
        inter_work_block[siconos::integrators::EulerMoreauOSI::XFREE]->setVectorPtr(
            0, workVds1[siconos::integrators::EulerMoreauOSI::FREE]);

      if (!inter_work_block[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS]) {
        inter_work_block[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS].reset(
            new BlockVector());
        inter_work_block[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS]->insertPtr(
            workVds1[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS_FOR_RELATION]);
      }
      else
        inter_work_block[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS]->setVectorPtr(
            0, workVds1[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS_FOR_RELATION]);
      if (!inter_work_block[siconos::integrators::EulerMoreauOSI::DELTA_X]) {
        inter_work_block[siconos::integrators::EulerMoreauOSI::DELTA_X].reset(
            new BlockVector());
        inter_work_block[siconos::integrators::EulerMoreauOSI::DELTA_X]->insertPtr(
            workVds1[siconos::integrators::EulerMoreauOSI::DELTA_X_FOR_RELATION]);
      }
      else
        inter_work_block[siconos::integrators::EulerMoreauOSI::DELTA_X]->setVectorPtr(
            0, workVds1[siconos::integrators::EulerMoreauOSI::DELTA_X_FOR_RELATION]);
    }
  }
  DEBUG_PRINTF("ds1->number() %i\n", ds1->number());
  DEBUG_PRINTF("ds2->number() %i\n", ds2->number());

  if (ds1 != ds2) {
    DEBUG_PRINT("ds1 != ds2\n");

    if (checkOSI(DSG.descriptor(ds2))) {
      DEBUG_PRINTF("ds2->number() %i is taken in to account\n", ds2->number());
      assert(DSG.properties(DSG.descriptor(ds2)).workVectors);
      auto& workVds2 = *DSG.properties(DSG.descriptor(ds2)).workVectors;
      if (relationType == FirstOrder) {
        if (!inter_work_block[siconos::integrators::EulerMoreauOSI::XFREE]) {
          inter_work_block[siconos::integrators::EulerMoreauOSI::XFREE].reset(
              new BlockVector());
          // dummy insertion to reserve first vector for ds1
          inter_work_block[siconos::integrators::EulerMoreauOSI::XFREE]->insertPtr(
              workVds2[siconos::integrators::EulerMoreauOSI::FREE]);
          inter_work_block[siconos::integrators::EulerMoreauOSI::XFREE]->insertPtr(
              workVds2[siconos::integrators::EulerMoreauOSI::FREE]);
        }
        else
          inter_work_block[siconos::integrators::EulerMoreauOSI::XFREE]->insertPtr(
              workVds2[siconos::integrators::EulerMoreauOSI::FREE]);

        if (!inter_work_block[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS]) {
          inter_work_block[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS].reset(
              new BlockVector());
          // dummy insertion to reserve first vector for ds1
          inter_work_block[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS]->insertPtr(
              workVds2[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS_FOR_RELATION]);
          inter_work_block[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS]->insertPtr(
              workVds2[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS_FOR_RELATION]);
        }
        else
          inter_work_block[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS]->insertPtr(
              workVds2[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS_FOR_RELATION]);

        if (!inter_work_block[siconos::integrators::EulerMoreauOSI::DELTA_X]) {
          inter_work_block[siconos::integrators::EulerMoreauOSI::DELTA_X].reset(
              new BlockVector());
          // dummy insertion to reserve first vector for ds1
          inter_work_block[siconos::integrators::EulerMoreauOSI::DELTA_X]->insertPtr(
              workVds2[siconos::integrators::EulerMoreauOSI::DELTA_X_FOR_RELATION]);
          inter_work_block[siconos::integrators::EulerMoreauOSI::DELTA_X]->insertPtr(
              workVds2[siconos::integrators::EulerMoreauOSI::DELTA_X_FOR_RELATION]);
        }
        else
          inter_work_block[siconos::integrators::EulerMoreauOSI::DELTA_X]->insertPtr(
              workVds2[siconos::integrators::EulerMoreauOSI::DELTA_X_FOR_RELATION]);
      }
    }
  }
}

void siconos::integrators::EulerMoreauOSI::initializeIterationMatrixW(
    double time, std::shared_ptr<siconos::modeling::DynamicalSystem> ds)
{
  // This function:
  // - allocate memory for the matrix W
  // - update its content for the current (initial) state of the dynamical system, depending on
  // its type.

  if (!ds)
    THROW_EXCEPTION(
        "siconos::integrators::EulerMoreauOSI::initializeIterationMatrixW(t,ds) - ds == "
        "nullptr");

  if (!(checkOSI(_dynamicalSystemsGraph->descriptor(ds))))
    THROW_EXCEPTION(
        "siconos::integrators::EulerMoreauOSI::initializeIterationMatrixW(t,ds) - ds does not "
        "belong to the OSI.");

  const DynamicalSystemsGraph::VDescriptor& dsv = _dynamicalSystemsGraph->descriptor(ds);

  if (_dynamicalSystemsGraph->properties(dsv).W)
    THROW_EXCEPTION(
        "siconos::integrators::EulerMoreauOSI::initializeIterationMatrixW(t,ds) - W(ds) is "
        "already in the map and has been initialized.");

  unsigned int sizeW = ds->dimension();  // n for first order systems, ndof for lagrangian.
  // Memory allocation for W
  double h = _simulation->timeStep();
  Type::Siconos dsType = Type::value(*ds);

  // 1 - All 'First order' systems
  if (auto d = std::dynamic_pointer_cast<siconos::modeling::FirstOrderNonLinearDS>(ds)) {
    // W =  M - h*_theta* [jacobian_x f(t,x,z)]

    // Memory allocation for W property of the graph
    if (d->M())  // W = M
    {
      d->computeM(time);
      _dynamicalSystemsGraph->properties(dsv).W = std::make_shared<siconos::algebra::SimpleMatrix>(*d->M()));
    }
    else  // W = I
    {
      _dynamicalSystemsGraph->properties(dsv).W = std::make_shared<siconos::algebra::SimpleMatrix>(sizeW, sizeW));
      _dynamicalSystemsGraph->properties(dsv).W->eye();
    }

    auto W = _dynamicalSystemsGraph->properties(dsv).W;
    // Add -h*_theta*jacobian_XF to W
    if (d->jacobianfx()) {
      d->computeJacobianfx(time, ds->x());
      siconos::algebra::scal(-h * _theta, *d->jacobianfx(), *W, false);
    }
  }
  else
    THROW_EXCEPTION(
        "siconos::integrators::EulerMoreauOSI::initializeIterationMatrixW implemented only "
        "for FirstOrderNonLinearDS (and heirs)\n");

  // Remark: W is not LU-factorized nor inversed here.
  // Function PLUForwardBackward will do that if required.
}

void siconos::integrators::EulerMoreauOSI::initializeIterationMatrixWBoundaryConditions(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds)
{
  // This function:
  // - allocate memory for a matrix WBoundaryConditions
  // - insert this matrix into WBoundaryConditionsMap with ds as a key

  if (!ds)
    THROW_EXCEPTION(
        "siconos::integrators::EulerMoreauOSI::initializeIterationMatrixWBoundaryConditions(t,"
        "ds) - ds == nullptr");

  if (!(checkOSI(_dynamicalSystemsGraph->descriptor(ds))))
    THROW_EXCEPTION(
        "siconos::integrators::EulerMoreauOSI::initializeIterationMatrixW(t,ds) - ds does not "
        "belong to the OSI.");

  THROW_EXCEPTION(
      "siconos::integrators::EulerMoreauOSI::initializeIterationMatrixWBoundaryConditions - "
      "not yet implemented.");
}

void siconos::integrators::EulerMoreauOSI::computeWBoundaryConditions(
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds)
{
  // Compute WBoundaryConditions matrix of the Dynamical System ds, at
  // time t and for the current ds state.

  // When this function is called, WBoundaryConditionsMap[ds] is
  // supposed to exist and not to be null Memory allocation has been
  // done during initializeIterationMatrixWBoundaryConditions.

  assert(ds &&
         "siconos::integrators::EulerMoreauOSI::computeWBoundaryConditions(t,ds) - ds == "
         "nullptr");

  // unsigned int dsN = ds->number();
  THROW_EXCEPTION(
      "siconos::integrators::EulerMoreauOSI::computeWBoundaryConditions - not yet "
      "implemented.");
}

void siconos::integrators::EulerMoreauOSI::computeW(
    double time, DynamicalSystem& ds, siconos::graphs::DynamicalSystemsGraph::VDescriptor& dsv,
    siconos::algebra::SiconosMatrix& W)
{
  DEBUG_BEGIN("siconos::integrators::EulerMoreauOSI::computeW(...)\n");
  // Compute W matrix of the Dynamical System ds, at time t and for the current ds state.

  // When this function is called, W is supposed to exist and not to be null
  // Memory allocation has been done during initializeIterationMatrixW.

  double h = _simulation->timeStep();

  // 1 - First order linear systems
  if (auto fods = std::dynamic_cast<siconos::modeling::FirstOrderLinearDS>(&ds)) {
    if (not std::dynamic_cast<siconos::modeling::FirstOrderLinearTIDS>(&ds)) {
      fods.computeA(time);
      fods.computeM(time);
    }

    if (fods.M())
      W = *fods.M();
    else
      W.eye();

    if (fods.A()) siconos::algebra::scal(-h * _theta, *fods.A(), W, false);
  }
  // 2 - First order non linear systems
  else if (auto d = std::dynamic_cast<siconos::modeling::FirstOrderNonLinearDS>(&ds)) {
    // W =  M - h*_theta* [jacobian_x f(t,x,z)]
    // Copy M or I if M is Null into W
    if (d.M()) {
      d.computeM(time);
      W = *d.M();
    }
    else
      W.eye();

    if (d.jacobianfx()) {
      d.computeJacobianfx(time, d.x());
      // Add -h*_theta*jacobianfx to W
      siconos::algebra::scal(-h * _theta, *d.jacobianfx(), W, false);
    }

    DEBUG_EXPR(W.display(););
  }
  else
    THROW_EXCEPTION(
        "siconos::integrators::EulerMoreauOSI::computeW - not yet implemented for this type "
        "of "
        "Dynamical system.");

  //  if (_useGamma)
  {
    //    InteractionsGraph& indexSet =
    //    *_simulation->nonSmoothDynamicalSystem()->topology()->indexSet(0);

    siconos::graphs::DynamicalSystemsGraph::OEIterator oei, oeiend;
    //    siconos::graphs::InteractionsGraph::VDescriptor ivd;
    std::shared_ptr<siconos::algebra::SiconosMatrix> K;
    std::shared_ptr<siconos::modeling::Interaction> inter;
    for (std::tie(oei, oeiend) = _dynamicalSystemsGraph->out_edges(dsv); oei != oeiend;
         ++oei) {
      inter = _dynamicalSystemsGraph->bundle(*oei);
      auto& relationMat = inter->relationMatrices();
      //      ivd = indexSet.descriptor(inter);
      auto& rel = static_cast<siconos::modeling::FirstOrderR&>(*inter->relation());
      K = rel.K();
      if (!K) K = relationMat[siconos::modeling::FirstOrderR::mat_K];
      if (K) {
        siconos::algebra::scal(-h * _gamma, *K, W, false);
      }
    }
  }
  // Remark: W is not LU-factorized here.
  // Function PLUForwardBackward will do that if required.
  DEBUG_EXPR(W.display());
  DEBUG_END("siconos::integrators::EulerMoreauOSI::computeW(...)\n");
}

void siconos::integrators::EulerMoreauOSI::computeKhat(
    siconos::modeling::Interaction& inter, siconos::algebra::SiconosMatrix& m,
    std::vector<std::shared_ptr<siconos::algebra::SimpleMatrix>>& workM, double h) const
{
  auto relationType = inter.relation()->getType();

  if ((relationType == FirstOrder) &&
      (workM[siconos::integrators::EulerMoreauOSI::MAT_KHAT])) {
    auto K = std::static_pointer_cast<FirstOrderR>(inter.relation())->K();
    if (!K) K = inter.relationMatrices()[FirstOrderR::mat_K];
    siconos::algebra::prod(*K, m, *workM[siconos::integrators::EulerMoreauOSI::MAT_KHAT],
                           true);
    *workM[siconos::integrators::EulerMoreauOSI::MAT_KHAT] *= h;
  }
}

double siconos::integrators::EulerMoreauOSI::computeResidu()
{
  DEBUG_BEGIN("siconos::integrators::EulerMoreauOSI::computeResidu()\n");
  // This function is used to compute the residu for each "EulerMoreauOSI-discretized"
  // dynamical system. It then computes the norm of each of them and finally return the maximum
  // value for those norms.
  //
  // The state values used are those saved in the DS, ie the last computed ones.
  //  $\mathcal R(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) - h r$
  //  $\mathcal R_{free}(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) $

  double time = _simulation->nextTime();      // End of the time step
  double told = _simulation->startingTime();  // Beginning of the time step
  double h = time - told;                     // time step length

  DEBUG_PRINTF("nextTime %f\n", time);
  DEBUG_PRINTF("startingTime %f\n", told);
  DEBUG_PRINTF("time step size %f\n", h);

  // Operators computed at told have index i, and (i+1) at t.

  // Iteration through the set of Dynamical Systems.
  //
  std::shared_ptr<siconos::modeling::DynamicalSystem> ds;  // Current Dynamical System.
  Type::Siconos dsType;                                    // Type of the current DS.

  double maxResidu = 0;
  double normResidu = maxResidu;

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    ds = _dynamicalSystemsGraph->bundle(*dsi);
    auto& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    dsType = Type::value(*ds);  // Its type

    // XXX TMP hack -- xhub
    // we have to iterate over the edges of the DSG0 -> the following won't be necessary
    // anymore Maurice will do that with subgraph :)
    siconos::algebra::SiconosVector& residuFree =
        *ds_work_vectors[siconos::integrators::EulerMoreauOSI::RESIDU_FREE];
    siconos::algebra::SiconosVector& residu =
        *ds_work_vectors[siconos::integrators::EulerMoreauOSI::RESIDU];

    // 1 - First Order Non Linear Systems AND First Order Linear DS
    if (dsType == Type::FirstOrderNonLinearDS || dsType == Type::FirstOrderLinearDS) {
      // ResiduFree = M(x_k,i+1 - x_i) - h*theta*f(t,x_k,i+1) - h*(1-theta)*f(ti,xi)
      //  $\mathcal R(x,r) = M(x - x_{k}) -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) - h
      //  r$
      //  $\mathcal R_{free}(x,r) = M(x - x_{k}) -h\theta f( x , t_{k+1}) -
      //  h(1-\theta)f(x_k,t_k) $

      // Note: indices k/k+1 corresponds to value at the beginning/end of the time step.
      // Newton iterate are x and r
      FirstOrderNonLinearDS& fonlds = *std::static_pointer_cast<FirstOrderNonLinearDS>(ds);

      // 1 - Compute the free residu (purely on the "smooth" dynamics)

      residuFree = *(fonlds.x());  // last saved value for x: could be x_k or x_{k+1}^alpha
      const siconos::algebra::SiconosVector& xold =
          fonlds.xMemory().getsiconos::algebra::SiconosVector(0);
      residuFree -= xold;  // state x_k (at previous time step)

      auto M = fonlds.M();
      if (M) {
        fonlds.computeM(time);
        siconos::algebra::prod(*M, residuFree, residuFree, true);
      }
      // at this step, we have residuFree = M(x - x_k)
      DEBUG_PRINT(
          "siconos::integrators::EulerMoreauOSI::computeResidu residuFree = M(x - x_k)\n");
      DEBUG_EXPR(residuFree.display());
      double coef = -h * (1 - _theta);
      if (dsType == Type::FirstOrderLinearDS) {
        // computes f(t_k,x_k)
        // No fold in FirstOrderLinearDS.
        // residu is used as a tmp buffer to compute Ax + b
        residu.zero();
        FirstOrderLinearDS& folds = *std::static_pointer_cast<FirstOrderLinearDS>(ds);
        if (folds.A()) {
          folds.computeA(told);
          siconos::algebra::prod(*folds.A(), xold, residu);
        }

        if (folds.b()) {
          folds.computeb(told);
          residu += *folds.b();
        }
        DEBUG_EXPR(residuFree.display());
        // residuFree += -h * (1 - _theta) * f(t_k,x_k)
        siconos::algebra::scal(coef, residu, residuFree, false);
        residu.zero();
        if (folds.A()) {
          folds.computeA(time);
          siconos::algebra::prod(*folds.A(), *folds.x(), residu);
        }
        if (folds.b()) {
          folds.computeb(time);
          residu += *folds.b();
        }
        // residuFree += -h * _theta * f(t_{x+1}, x_{k+1}^alpha)
        coef = -h * _theta;
        siconos::algebra::scal(coef, residu, residuFree, false);
        DEBUG_PRINT("- 3 -\n");
        DEBUG_EXPR(residuFree.display());
        DEBUG_EXPR(xold.display());
        DEBUG_EXPR(folds.x()->display());
      }
      else if (dsType == Type::FirstOrderNonLinearDS)  // FirstOrderNonLinearDS
      {
        DEBUG_PRINT("dsType == Type::FirstOrderNonLinearDS\n");
        DEBUG_EXPR(fonlds.f()->display(););
        if (fonlds.f()) {
          coef = -h * (1 - _theta);
          // for these systems, fold is available
          // residuFree += -h * (1 - _theta) * f(t_k,x_k)
          siconos::algebra::scal(coef, *fonlds.fold(), residuFree, false);

          // computes f(t_{x+1}, x_{k+1}^alpha)
          fonlds.computef(time, fonlds.x());
          coef = -h * _theta;
          // residuFree += -h * _theta * f(t_{x+1}, x_{k+1}^alpha)
          siconos::algebra::scal(coef, *(fonlds.f()), residuFree, false);
        }
      }

      // now we compute the residu = residuFree - h*gamma*r - h*(1-gamma)r_k
      residu = residuFree;

      if (!_useGamma)  // no gamma
      {
        DEBUG_EXPR(fonlds.r()->display(););
        DEBUG_EXPR(residu.display());
        siconos::algebra::scal(-h, *fonlds.r(), residu, false);  // residu = residu - h*r
      }
      else {
        siconos::algebra::scal(-h * _gamma, *fonlds.r(), residu, false);
        siconos::algebra::scal(-h * (1 - _gamma),
                               fonlds.rMemory().getsiconos::algebra::SiconosVector(0), residu,
                               false);
      }

      normResidu = residu.norm2();
      DEBUG_EXPR(residu.display());
    }
    // 2 - First Order Linear Systems with Time Invariant coefficients
    else if (dsType == Type::FirstOrderLinearTIDS) {
      FirstOrderLinearTIDS& foltids = *std::static_pointer_cast<FirstOrderLinearTIDS>(ds);
      // Don't use W because it is LU factorized
      // Residu : R_{free} = M(x^{\alpha}_{k+1} - x_{k}) -h( A (\theta x^{\alpha}_{k+1} +
      // (1-\theta)  x_k) +b)

      // 1. R_{free} = -h * b
      if (foltids.b())
        siconos::algebra::scal(-h, *(foltids.b()), residuFree, true);
      else
        residuFree.zero();

      // 2. residuFree += -h * A (\theta x_{k+1}^{\alpha} + (1-\theta) x_k)
      // residu is used as a temp buffer
      if (foltids.A()) {
        auto A = foltids.A();
        siconos::algebra::prod(*A, foltids.xMemory().getsiconos::algebra::SiconosVector(0),
                               residu, true);
        double coef = -h * (1 - _theta);
        siconos::algebra::scal(coef, residu, residuFree, false);

        siconos::algebra::prod(*A, *(foltids.x()), residu, true);
        coef = -h * _theta;
        siconos::algebra::scal(coef, residu, residuFree, false);
      }

      // 3. residuFree += M(x_{k+1}^{\alpha} - x_k)
      residu = *(foltids.x()) - foltids.xMemory().getsiconos::algebra::SiconosVector(0);
      auto M = foltids.M();
      if (M) {
        siconos::algebra::prod(*M, residu, residuFree, false);
      }
      else {
        residuFree += residu;
      }
    }
    else
      THROW_EXCEPTION(
          "siconos::integrators::EulerMoreauOSI::computeResidu - not yet implemented for "
          "Dynamical system type: " +
          std::to_string(dsType));

    DEBUG_PRINT("siconos::integrators::EulerMoreauOSI::computeResidu final residuFree\n");
    DEBUG_EXPR(residuFree.display());

    if (normResidu > maxResidu) maxResidu = normResidu;
  }
  DEBUG_END("siconos::integrators::EulerMoreauOSI::computeResidu()\n");
  return maxResidu;
}

void siconos::integrators::EulerMoreauOSI::computeFreeState()
{
  // This function computes "free" states of the DS belonging to this Integrator.
  // "Free" means without taking non-smooth effects into account.
  DEBUG_BEGIN("siconos::integrators::EulerMoreauOSI::computeFreeState()\n");

  double t = _simulation->nextTime();         // End of the time step
  double told = _simulation->startingTime();  // Beginning of the time step
  double h = t - told;                        // time step length

  // Operators computed at told have index i, and (i+1) at t.

  //  Note: integration of r with a theta method has been removed
  //  siconos::algebra::SiconosVector *rold =
  //  static_cast<siconos::algebra::SiconosVector*>(d.rMemory()->getsiconos::algebra::SiconosVector(0));

  // Iteration through the set of Dynamical Systems.
  //

  std::shared_ptr<siconos::modeling::DynamicalSystem> ds;  // Current Dynamical System.
  Type::Siconos dsType;                                    // Type of the current DS.

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    ds = _dynamicalSystemsGraph->bundle(*dsi);

    // XXX TMP hack -- xhub
    // we have to iterate over the edges of the DSG0 -> the following won't be necessary
    // anymore Maurice will do that with subgraph :)

    auto& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    dsType = Type::value(*ds);  // Its type
    auto& W = *_dynamicalSystemsGraph->properties(*dsi)
                   .W;  // Its W EulerMoreauOSI matrix of iteration.

    // 1 - First Order Non Linear Systems
    if (dsType == Type::FirstOrderNonLinearDS || dsType == Type::FirstOrderLinearDS ||
        dsType == Type::FirstOrderLinearTIDS) {
      // xfree =  x - W^{-1} (ResiduFree - h(1-gamma)*rold)
      // with ResiduFree = = M(x - x_k) - h*theta*f(t_{k+1}, x) - h*(1-theta)*f(t_k, x_k)

      // to be updated at current time: W, f
      // fold is f at t_k
      // not time dependant: M
      FirstOrderNonLinearDS& d = *std::static_pointer_cast<FirstOrderNonLinearDS>(ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      //    std::shared_ptr<siconos::algebra::SiconosVector> xold =
      //    d->xMemory()->getsiconos::algebra::SiconosVector(0); // xi

      siconos::algebra::SiconosVector& x = *d.x();  // x = x_k or x = x_{k+1}^{\alpha}
      // xfree gets ResiduFree at first
      siconos::algebra::SiconosVector& xfree =
          *ds_work_vectors[siconos::integrators::EulerMoreauOSI::FREE];
      xfree = *ds_work_vectors[siconos::integrators::EulerMoreauOSI::RESIDU_FREE];

      DEBUG_PRINT(
          "siconos::integrators::EulerMoreauOSI::computeFreeState xfree <- residuFree\n");
      DEBUG_EXPR(xfree.display());

      if (_useGamma) {
        const siconos::algebra::SiconosVector& rold =
            d.rMemory().getsiconos::algebra::SiconosVector(0);
        double coeff = -h * (1 - _gamma);
        siconos::algebra::scal(coeff, rold, xfree, false);  //  xfree += -h(1-gamma)*rold
      }

      // At this point xfree = (ResiduFree - h(1-gamma)*rold)
      // -> Solve WX = xfree and set xfree = X
      W.Solve(xfree);

      // at this point, xfree = W^{-1} (ResiduFree - h(1-gamma)*rold)
      // -> compute real xfree = x - W^{-1} (ResiduFree - h(1-gamma)*rold)
      xfree *= -1.0;
      xfree += x;

      DEBUG_EXPR(xfree.display());

      // now the crazy intermediate variables
      // xPartialNS was updated before this fonction call
      // It constains either 0 (first Newton iterate)
      // or g(x, \lambda, t_{k+1}) - B_{k+1}^{\alpha} \lambda - K_{k+1}^{\alpha} x
      siconos::algebra::SiconosVector& xPartialNS =
          *ds_work_vectors[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS_FOR_RELATION];
      DEBUG_PRINT(
          "siconos::integrators::EulerMoreauOSI::computeFreeState xPartialNS from "
          "Interaction\n");
      DEBUG_EXPR(xPartialNS.display());

      // -> Solve WX = g(x, \lambda, t_{k+1}) - B_{k+1}^{\alpha} \lambda - K_{k+1}^{\alpha} x
      // and set xPartialNS = X
      W.Solve(xPartialNS);
      siconos::algebra::scal(h, xPartialNS, xPartialNS);

      // compute real xPartialNS = xfree + ...
      xPartialNS += xfree;
      DEBUG_PRINT(
          "siconos::integrators::EulerMoreauOSI::computeFreeState xPartialNS real value\n");
      DEBUG_EXPR(xPartialNS.display());

      // deltaxForRelation = (\widetilde{K}_{k+1}^{\alpha})^{-1} xPartialNS - x
      siconos::algebra::SiconosVector& deltaxForRelation =
          *ds_work_vectors[siconos::integrators::EulerMoreauOSI::DELTA_X_FOR_RELATION];
      deltaxForRelation = xPartialNS;

      deltaxForRelation -= x;

      DEBUG_EXPR(deltaxForRelation.display());

      // have a look at the end of the DevNotes for this part
      if (_useGammaForRelation) {
        if (!(dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS))
          THROW_EXCEPTION(
              "siconos::integrators::EulerMoreauOSI::computeFreeState - _useGammaForRelation "
              "== true is only implemented for FirstOrderLinearDS or FirstOrderLinearTIDS");

        deltaxForRelation = xfree;

        siconos::algebra::scal(_gamma, deltaxForRelation, deltaxForRelation);
        const siconos::algebra::SiconosVector& xold =
            d.xMemory().getsiconos::algebra::SiconosVector(0);

        siconos::algebra::scal(1.0 - _gamma, xold, deltaxForRelation, false);
      }

      // some output
      DEBUG_EXPR(xfree.display(););
      DEBUG_EXPR(xPartialNS.display(););
      DEBUG_EXPR(deltaxForRelation.display(););
    }
    else
      THROW_EXCEPTION(
          "siconos::integrators::EulerMoreauOSI::computeFreeState - not yet implemented for "
          "Dynamical system type: " +
          std::to_string(dsType));
  }
  DEBUG_END("siconos::integrators::EulerMoreauOSI::computeFreeState()\n");
}

void siconos::integrators::EulerMoreauOSI::prepareNewtonIteration(double time)
{
  // XXX TMP hack -- xhub
  // we have to iterate over the edges of the DSG0 -> the following won't be necessary anymore
  // Maurice will do that with subgraph :)
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds =
        _dynamicalSystemsGraph->bundle(*dsi);
    DynamicalSystemsGraph::VDescriptor dsv = _dynamicalSystemsGraph->descriptor(ds);
    auto W = _dynamicalSystemsGraph->properties(*dsi).W;
    computeW(time, *ds, dsv, *W);
  }

  if (!_explicitJacobiansOfRelation) {
    _simulation->nonSmoothDynamicalSystem()->computeInteractionJacobians(time);

    siconos::graphs::InteractionsGraph::VIterator ui, uiend;
    auto indexSet0 = _simulation->nonSmoothDynamicalSystem()->topology()->indexSet0();

    for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
      auto& inter = *indexSet0->bundle(*ui);
      siconos::graphs::InteractionProperties& interProp = indexSet0->properties(*ui);

      auto& inter_work = *interProp.workVectors;
      auto& relationMat = inter.relationMatrices();
      auto& inter_work_block = *interProp.workBlockVectors;

      auto relationType = inter.relation()->getType();
      auto relationSubType = inter.relation()->getSubType();
      if (relationType == FirstOrder) {
        FirstOrderR& relation = static_cast<FirstOrderR&>(*inter.relation());
        BlockVector& xPartialNS =
            *inter_work_block[siconos::integrators::EulerMoreauOSI::X_PARTIAL_NS];

        if (relationSubType == NonLinearR || relationSubType == Type2R) {
          if (relation.B())
            siconos::algebra::prod(*relation.B(), *inter.lambda(0),
                                   *inter_work[siconos::integrators::EulerMoreauOSI::VEC_X],
                                   true);
          else
            siconos::algebra::prod(*relationMat[FirstOrderR::mat_B], *inter.lambda(0),
                                   *inter_work[siconos::integrators::EulerMoreauOSI::VEC_X],
                                   true);

          xPartialNS = *inter_work_block[siconos::integrators::EulerMoreauOSI::G_ALPHA];
          xPartialNS -= *inter_work[siconos::integrators::EulerMoreauOSI::VEC_X];
        }
      }
    }
  }
  // update alpha variables
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  auto indexSet0 = _simulation->nonSmoothDynamicalSystem()->topology()->indexSet0();

  for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
    auto& inter = *indexSet0->bundle(*ui);
    siconos::graphs::InteractionProperties& interProp = indexSet0->properties(*ui);

    auto& inter_work = *interProp.workVectors;
    *inter_work[siconos::integrators::EulerMoreauOSI::YOLD] = *(inter.y(0));
    *inter_work[siconos::integrators::EulerMoreauOSI::LAMBDAOLD] = *(inter.lambda(0));
  }
}

/// @cond

struct siconos::integrators::EulerMoreauOSI::_NSLEffectOnFreeOutput : public SiconosVisitor {
  using SiconosVisitor::visit;

  OneStepNSProblem* _osnsp;
  std::shared_ptr<siconos::modeling::Interaction> _inter;

  _NSLEffectOnFreeOutput(OneStepNSProblem* p,
                         std::shared_ptr<siconos::modeling::Interaction> inter)
      : _osnsp(p), _inter(inter){};

  void visit(const EqualityConditionNSL& nslaw) { ; }
  void visit(const MixedComplementarityConditionNSL& nslaw) { ; }
};

/// @endcond

void siconos::integrators::EulerMoreauOSI::computeFreeOutput(
    siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp)
{
  /** \warning: ensures that it can also work with two different osi for two different ds ?
   */
  DEBUG_BEGIN("siconos::integrators::EulerMoreauOSI::computeFreeOutput(...)\n");
  auto allOSNS = _simulation->oneStepNSProblems();
  auto indexSet = osnsp->simulation()->indexSet(osnsp->indexSetLevel());
  std::shared_ptr<siconos::modeling::Interaction> inter = indexSet->bundle(vertex_inter);

  auto& DSlink = inter->linkToDSVariables();
  auto& relationMat = inter->relationMatrices();
  auto& relationVec = inter->relationVectors();

  auto& inter_work = *indexSet->properties(vertex_inter).workVectors;
  auto& inter_work_block = *(indexSet->properties(vertex_inter)).workBlockVectors;
  // Get relation and non smooth law types
  auto relationType = inter->relation()->getType();
  auto relationSubType = inter->relation()->getSubType();

  unsigned int sizeY = inter->nonSmoothLaw()->size();

  unsigned int relativePosition = 0;

  Index coord(8);
  coord[0] = relativePosition;
  coord[1] = relativePosition + sizeY;
  coord[2] = 0;
  coord[4] = 0;
  coord[6] = 0;
  coord[7] = sizeY;
  std::shared_ptr<siconos::algebra::SiconosMatrix> C;
  std::shared_ptr<siconos::algebra::SiconosMatrix> D;
  std::shared_ptr<siconos::algebra::SiconosMatrix> F;
  std::shared_ptr<siconos::algebra::BlockVector> deltax;
  std::shared_ptr<siconos::algebra::BlockVector> Xfree;

  std::shared_ptr<siconos::algebra::SiconosVector> H_alpha;

  deltax = inter_work_block[siconos::integrators::EulerMoreauOSI::DELTA_X];
  DEBUG_EXPR(deltax->display(););
  siconos::algebra::SiconosVector& osnsp_rhs =
      *(*indexSet->properties(vertex_inter)
             .workVectors)[siconos::integrators::EulerMoreauOSI::OSNSP_RHS];

  Xfree = inter_work_block[siconos::integrators::EulerMoreauOSI::XFREE];
  DEBUG_EXPR(Xfree->display(););
  assert(Xfree);

  std::shared_ptr<siconos::modeling::Interaction> mainInteraction = inter;
  assert(mainInteraction);
  assert(mainInteraction->relation());

  if (relationType == FirstOrder &&
      (relationSubType == Type2R || relationSubType == NonLinearR)) {
    DEBUG_PRINT(
        "relationType == FirstOrder && (relationSubType == Type2R || relationSubType == "
        "NonLinearR)\n")
    siconos::algebra::SiconosVector& lambda = *inter->lambda(0);
    FirstOrderR& rel = *std::static_pointer_cast<FirstOrderR>(mainInteraction->relation());
    C = rel.C();
    if (!C) C = relationMat[FirstOrderR::mat_C];
    D = rel.D();
    if (!D) D = relationMat[FirstOrderR::mat_D];

    if (D) {
      coord[3] = D->size(1);
      coord[5] = D->size(1);
      siconos::algebra::subprod(*D, lambda, osnsp_rhs, coord, true);

      osnsp_rhs *= -1.0;
    }
    else {
      siconos::algebra::subscal(0, osnsp_rhs, osnsp_rhs, coord, true);
    }

    if (C) {
      coord[3] = C->size(1);
      coord[5] = C->size(1);
      siconos::algebra::subprod(*C, *deltax, osnsp_rhs, coord, false);
    }

    if (_useGammaForRelation) {
      THROW_EXCEPTION(
          "siconos::integrators::EulerMoreauOSI::ComputeFreeOutput not yet implemented with "
          "useGammaForRelation() for FirstorderR and Type2R and H_alpha->getValue() should "
          "return the mid-point value");
    }
    siconos::algebra::SiconosVector& hAlpha =
        *inter_work[siconos::integrators::EulerMoreauOSI::H_ALPHA];
    DEBUG_EXPR(hAlpha.display());
    osnsp_rhs += hAlpha;
    DEBUG_EXPR(osnsp_rhs.display(););
  }
  else if (relationType == FirstOrder && relationSubType == Type1R) {
    DEBUG_PRINT("relationType == FirstOrder && relationSubType == Type1R\n");
    FirstOrderType1R& rel =
        *std::static_pointer_cast<FirstOrderType1R>(mainInteraction->relation());
    C = rel.C();
    if (!C) C = relationMat[FirstOrderR::mat_C];
    F = rel.F();
    if (!F) F = relationMat[FirstOrderR::mat_F];
    assert(Xfree);
    assert(deltax);

    if (F) {
      coord[3] = F->size(1);
      coord[5] = F->size(1);
      siconos::algebra::subprod(*F, *DSlink[FirstOrderR::z], osnsp_rhs, coord, true);
    }
    if (C) {
      coord[3] = C->size(1);
      coord[5] = C->size(1);
      siconos::algebra::subprod(*C, *Xfree, osnsp_rhs, coord, false);
    }

    if (_useGammaForRelation) {
      THROW_EXCEPTION(
          "siconos::integrators::EulerMoreauOSI::ComputeFreeOutput not yet implemented with "
          "useGammaForRelation() for FirstorderR and Typ2R and H_alpha->getValue() should "
          "return the mid-point value");
    }
    if (inter_work[siconos::integrators::EulerMoreauOSI::H_ALPHA]) {
      osnsp_rhs += *inter_work[siconos::integrators::EulerMoreauOSI::H_ALPHA];
    }
  }
  else  // First Order Linear Relation
  {
    DEBUG_PRINT("relationType == FirstOrder\n");
    C = mainInteraction->relation()->C();
    if (!C) C = relationMat[FirstOrderR::mat_C];

    if (C) {
      assert(Xfree);
      assert(deltax);

      coord[3] = C->size(1);
      coord[5] = C->size(1);

      if (_useGammaForRelation) {
        siconos::algebra::subprod(*C, *deltax, osnsp_rhs, coord, true);
      }
      else {
        siconos::algebra::subprod(*C, *Xfree, osnsp_rhs, coord, true);
      }
    }
    DEBUG_EXPR(osnsp_rhs.display(););
    if (relationType == FirstOrder &&
        (relationSubType == LinearTIR || relationSubType == LinearR)) {
      // In the first order linear case it may be required to add e + FZ to y.
      // y = CXfree + e + FZ
      std::shared_ptr<siconos::algebra::SiconosVector> e;
      if (relationSubType == LinearTIR) {
        e = std::static_pointer_cast<FirstOrderLinearTIR>(mainInteraction->relation())->e();
        F = std::static_pointer_cast<FirstOrderLinearTIR>(mainInteraction->relation())->F();
      }
      else {
        e = std::static_pointer_cast<FirstOrderLinearR>(mainInteraction->relation())->e();
        if (!e) e = relationVec[FirstOrderR::e];
        F = std::static_pointer_cast<FirstOrderLinearR>(mainInteraction->relation())->F();
        if (!F) F = relationMat[FirstOrderR::mat_F];
      }

      if (e) osnsp_rhs += *e;

      if (F) {
        coord[3] = F->size(1);
        coord[5] = F->size(1);
        siconos::algebra::subprod(*F, *DSlink[FirstOrderR::z], osnsp_rhs, coord, false);
      }
    }
    DEBUG_EXPR(osnsp_rhs.display(););
  }
  DEBUG_END("siconos::integrators::EulerMoreauOSI::computeFreeOutput(...)\n");
}

void siconos::integrators::EulerMoreauOSI::integrate(double& tinit, double& tend, double& tout,
                                                     int&)
{
  // Last parameter is not used (required for LsodarOSI but not for EulerMoreauOSI).

  // double h = tend - tinit;
  tout = tend;

  std::shared_ptr<siconos::algebra::SiconosMatrix> W;
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds =
        _dynamicalSystemsGraph->bundle(*dsi);
    Type::Siconos dsType = Type::value(*ds);
    THROW_EXCEPTION(
        "siconos::integrators::EulerMoreauOSI::integrate - not yet implemented for Dynamical "
        "system type :" +
        std::to_string(dsType));
  }
}

void siconos::integrators::EulerMoreauOSI::updateState(const unsigned int)
{
  DEBUG_PRINT("siconos::integrators::EulerMoreauOSI::updateState\n");

  double h = _simulation->timeStep();

  double RelativeTol = _simulation->relativeConvergenceTol();
  bool useRCC = _simulation->useRelativeConvergenceCriteron();
  if (useRCC) _simulation->setRelativeConvergenceCriterionHeld(true);

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;

  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds =
        _dynamicalSystemsGraph->bundle(*dsi);

    // Get the DS type
    Type::Siconos dsType = Type::value(*ds);
    auto& ds_work_vectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

    SimpleMatrix& W = *_dynamicalSystemsGraph->properties(*dsi).W;

    if (dsType == Type::FirstOrderNonLinearDS || dsType == Type::FirstOrderLinearDS ||
        dsType == Type::FirstOrderLinearTIDS) {
      FirstOrderNonLinearDS& d = *std::static_pointer_cast<FirstOrderNonLinearDS>(ds);
      siconos::algebra::SiconosVector& x = *ds->x();
      DEBUG_PRINT("siconos::integrators::EulerMoreauOSI::updateState Old value of x\n");
      DEBUG_EXPR(x.display());
      DEBUG_PRINT("siconos::integrators::EulerMoreauOSI::updateState residu value\n");
      DEBUG_EXPR(d.r()->display());

      // TODO ???
      bool baux = (useRCC && dsType == Type::FirstOrderNonLinearDS &&
                   _simulation->relativeConvergenceCriterionHeld());

      //    std::shared_ptr<siconos::algebra::SiconosVector> xFree = d->xFree();

      // Save value of q in local_buffer for relative convergence computation
      if (baux) *ds_work_vectors[siconos::integrators::EulerMoreauOSI::LOCAL_BUFFER] = x;

      if (_useGamma) {
        // XXX UseGamma broken ? -- xhub
        siconos::algebra::scal(_gamma * h, *d.r(), x);  // x = gamma*h*r
      }
      else {
        siconos::algebra::scal(h, *d.r(), x);  // x = h*r
      }

      W.Solve(x);  // x = h* W^{-1} *r

      x += *ds_work_vectors[siconos::integrators::EulerMoreauOSI::FREE];  // x+=xfree

      if (baux) {
        double ds_norm_ref = 1. + ds->x0()->norm2();  // Should we save this in the graph?
        *ds_work_vectors[siconos::integrators::EulerMoreauOSI::LOCAL_BUFFER] -= x;
        double aux =
            (ds_work_vectors[siconos::integrators::EulerMoreauOSI::LOCAL_BUFFER]->norm2()) /
            (ds_norm_ref);
        if (aux > RelativeTol) _simulation->setRelativeConvergenceCriterionHeld(false);
      }
      DEBUG_PRINT("siconos::integrators::EulerMoreauOSI::updateState New value of x\n");
      DEBUG_EXPR(x.display());
    }
    else
      THROW_EXCEPTION(
          "siconos::integrators::EulerMoreauOSI::updateState - not yet implemented for "
          "Dynamical system type: " +
          std::to_string(dsType));
  }
}

void siconos::integrators::EulerMoreauOSI::display()
{
  OneStepIntegrator::display();

  std::cout << "====== EulerMoreauOSI OSI display ======" << std::endl;

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds =
        _dynamicalSystemsGraph->bundle(*dsi);
    std::cout << "--------------------------------" << std::endl;
    std::cout << "--> W of dynamical system number " << ds->number() << ": " << std::endl;
    if (_dynamicalSystemsGraph->properties(*dsi).W)
      _dynamicalSystemsGraph->properties(*dsi).W->display();
    else
      std::cout << "-> nullptr" << std::endl;
    std::cout << "--> and corresponding theta is: " << _theta << std::endl;
  }
  std::cout << "================================" << std::endl;
}
void siconos::integrators::EulerMoreauOSI::updateOutput(double time)
{
  /** VA. 16/02/2017 This should normally be done only for interaction managed by the osi */
  for (unsigned int level = _levelMinForOutput; level < _levelMaxForOutput + 1; level++)
    updateOutput(time, level);
}

void siconos::integrators::EulerMoreauOSI::updateInput(double time)
{
  /** VA. 16/02/2017 This should normally be done only for interaction managed by the osi */
  for (unsigned int level = _levelMinForInput; level < _levelMaxForInput + 1; level++)
    updateInput(time, level);
}

void siconos::integrators::EulerMoreauOSI::updateOutput(double time, unsigned int level)
{
  DEBUG_BEGIN(
      "siconos::integrators::EulerMoreauOSI::updateOutput(double time, unsigned int level)\n");
  /** VA. 16/02/2017 This should normally be done only for interaction managed by the osi */
  //_simulation->nonSmoothDynamicalSystem()->updateOutput(time,level);
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  auto indexSet0 = _simulation->nonSmoothDynamicalSystem()->topology()->indexSet0();
  for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
    auto& inter = *indexSet0->bundle(*ui);
    assert(inter.lowerLevelForOutput() <= level);
    assert(inter.upperLevelForOutput() >= level);

    auto& DSlink = inter.linkToDSVariables();
    auto& relationMat = inter.relationMatrices();

    siconos::graphs::InteractionProperties& interProp = indexSet0->properties(*ui);
    auto& inter_work = *interProp.workVectors;
    auto& inter_work_block = *interProp.workBlockVectors;
    auto relationSubType = inter.relation()->getSubType();
    if (relationSubType == Type2R) {
      FirstOrderType2R& r = static_cast<FirstOrderType2R&>(*inter.relation());
      // compute the new y obtained by linearisation (see DevNotes)
      // y_{alpha+1}_{k+1} = h(x_{k+1}^{alpha},lambda_{k+1}^{alpha},t_k+1)
      //                     + C_{k+1}^alpha ( x_{k+1}^{alpha+1}- x_{k+1}^{alpha} )
      //                     + D_{k+1}^alpha ( lambda_{k+1}^{alpha+1} - lambda_{k+1}^{alpha} )
      // or equivalently
      // y_{alpha+1}_{k+1} = y_{alpha}_{k+1} - ResiduY_{k+1}^{alpha}
      //                     + C_{k+1}^alpha ( x_{k+1}^{alpha+1}- x_{k+1}^{alpha} )
      //                     + D_{k+1}^alpha ( lambda_{k+1}^{alpha+1} - lambda_{k+1}^{alpha} )
      siconos::algebra::SiconosVector& y = *inter.y(level);
      DEBUG_EXPR(y.display());

      if (r.D())
        siconos::algebra::prod(
            *r.D(), *inter_work[siconos::integrators::EulerMoreauOSI::LAMBDAOLD], y, true);
      else
        siconos::algebra::prod(*relationMat[FirstOrderR::mat_D],
                               *inter_work[siconos::integrators::EulerMoreauOSI::LAMBDAOLD], y,
                               true);

      y *= -1.0;
      DEBUG_PRINT("FirstOrderType2R::computeOutput : y old(level) \n");
      DEBUG_EXPR(inter_work[siconos::integrators::EulerMoreauOSI::YOLD]->display());

      y += *inter_work[siconos::integrators::EulerMoreauOSI::YOLD]

           DEBUG_PRINT("siconos::integrators::EulerMoreauOSI::updateOutput : ResiduY() \n");
      siconos::algebra::SiconosVector& residuY =
          *inter_work[siconos::integrators::EulerMoreauOSI::VEC_RESIDU_Y];
      DEBUG_EXPR(residuY.display());

      y -= residuY;
      DEBUG_PRINT("siconos::integrators::EulerMoreauOSI::updateOutput : y(level) \n");
      DEBUG_EXPR(y.display());

      BlockVector& deltax = *inter_work_block[siconos::integrators::EulerMoreauOSI::DELTA_X];
      DEBUG_PRINT("siconos::integrators::EulerMoreauOSI::updateOutput : deltax \n");
      DEBUG_EXPR(deltax.display());

      if (r.C())
        siconos::algebra::prod(*r.C(), deltax, y, false);
      else
        siconos::algebra::prod(*relationMat[FirstOrderR::mat_C], deltax, y, false);

      DEBUG_PRINT("siconos::integrators::EulerMoreauOSI::updateOutput : y before osnsM\n");
      DEBUG_EXPR(y.display());
      if (interProp.block) {
        auto& osnsM = *interProp.block;
        siconos::algebra::prod(osnsM, *inter.lambda(level), y, false);
        DEBUG_EXPR(inter.lambda(level)->display());
        DEBUG_EXPR(osnsM.display());
        DEBUG_PRINT(
            "siconos::integrators::EulerMoreauOSI::updateOutput : new linearized y \n");
        DEBUG_EXPR(y.display());
      }

      siconos::algebra::SiconosVector& hAlpha =
          *inter_work[siconos::integrators::EulerMoreauOSI::H_ALPHA];

      r.computeh(time, *DSlink[FirstOrderR::x], *inter.lambda(level), hAlpha);
      DEBUG_PRINT("siconos::integrators::EulerMoreauOSI::updateOutput : new Halpha \n");
      DEBUG_EXPR(hAlpha.display());
    }
    else if (relationSubType == NonLinearR) {
      FirstOrderNonLinearR& r = static_cast<FirstOrderNonLinearR&>(*inter.relation());
      // compute the new y  obtained by linearisation (see DevNotes)
      // y_{alpha+1}_{k+1} = h(x_{k+1}^{alpha},lambda_{k+1}^{alpha},t_k+1)
      //                     + C_{k+1}^alpha ( x_{k+1}^{alpha+1}- x_{k+1}^{alpha} )
      //                     + D_{k+1}^alpha ( lambda_{k+1}^{alpha+1} - lambda_{k+1}^{alpha} )
      // or equivalently
      // y_{alpha+1}_{k+1} = y_{alpha}_{k+1} - ResiduY_{k+1}^{alpha}
      //                     + C_{k+1}^alpha ( x_{k+1}^{alpha+1}- x_{k+1}^{alpha} )
      //                     + D_{k+1}^alpha ( lambda_{k+1}^{alpha+1} - lambda_{k+1}^{alpha} )
      siconos::algebra::SiconosVector& y = *inter.y(level);
      DEBUG_EXPR(y.display());

      if (r.D())
        siconos::algebra::prod(
            *r.D(), *inter_work[siconos::integrators::EulerMoreauOSI::LAMBDAOLD], y, true);
      else
        siconos::algebra::prod(*relationMat[FirstOrderR::mat_D],
                               *inter_work[siconos::integrators::EulerMoreauOSI::LAMBDAOLD], y,
                               true);

      y *= -1.0;
      DEBUG_PRINT("FirstOrderNonLinearR::computeOutput : y Old(level) \n");
      DEBUG_EXPR(inter_work[siconos::integrators::EulerMoreauOSI::YOLD]->display());

      y += *inter_work[siconos::integrators::EulerMoreauOSI::YOLD];

      DEBUG_PRINT("siconos::integrators::EulerMoreauOSI::updateOutput : ResiduY() \n");
      auto& residuY = *inter_work[siconos::integrators::EulerMoreauOSI::VEC_RESIDU_Y];
      DEBUG_EXPR(residuY.display());

      y -= residuY;

      DEBUG_PRINT("siconos::integrators::EulerMoreauOSI::updateOutput : y(level) \n");
      DEBUG_EXPR(y.display());

      auto& deltax = *inter_work_block[siconos::integrators::EulerMoreauOSI::DELTA_X];
      DEBUG_PRINT("siconos::integrators::EulerMoreauOSI::updateOutput : deltax \n");
      DEBUG_EXPR(deltax.display());

      if (r.C())
        siconos::algebra::prod(*r.C(), deltax, y, false);
      else
        siconos::algebra::prod(*relationMat[FirstOrderR::mat_C], deltax, y, false);

      if (interProp.block) {
        auto& osnsM = *interProp.block;
        // osnsM = h * C * W^-1 * B + D
        DEBUG_EXPR(osnsM.display(););
        siconos::algebra::prod(osnsM, *inter.lambda(level), y, false);
      }
      DEBUG_PRINT("siconos::integrators::EulerMoreauOSI::updateOutput : new linearized y \n");
      DEBUG_EXPR(y.display());

      auto& hAlpha = *inter_work[siconos::integrators::EulerMoreauOSI::H_ALPHA];
      r.computeh(time, *DSlink[FirstOrderR::x], *inter.lambda(level), *DSlink[FirstOrderR::z],
                 hAlpha);
      DEBUG_EXPR(x.display(););
      DEBUG_PRINT("siconos::integrators::EulerMoreauOSI::updateOutput : new Halpha \n");
      DEBUG_EXPR(hAlpha.display());
    }
    else
      inter.computeOutput(time, level);
  }
  DEBUG_END(
      "siconos::integrators::EulerMoreauOSI::updateOutput(double time, unsigned int level)\n");
}

void siconos::integrators::EulerMoreauOSI::updateInput(double time, unsigned int level)
{
  /** VA. 16/02/2017 This should normally be done only for interaction managed by the osi */
  //_simulation->nonSmoothDynamicalSystem()->updateInput(time,level);

  // Set dynamical systems non-smooth part to zero.
  _simulation->nonSmoothDynamicalSystem()->reset(level);

  siconos::graphs::InteractionsGraph::VIterator ui, uiend;

  auto indexSet0 = _simulation->nonSmoothDynamicalSystem()->topology()->indexSet0();
  for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
    auto& inter = *indexSet0->bundle(*ui);
    assert(inter.lowerLevelForInput() <= level);
    assert(inter.upperLevelForInput() >= level);

    auto& DSlink = inter.linkToDSVariables();
    auto& relationMat = inter.relationMatrices();

    siconos::graphs::InteractionProperties& interProp = indexSet0->properties(*ui);
    auto& inter_work = *interProp.workVectors;
    auto& inter_work_mat = *interProp.workMatrices;
    auto& inter_work_block = *interProp.workBlockVectors;

    auto relationSubType = inter.relation()->getSubType();
    if (relationSubType == siconos::modeling::RelationSubType::Type2R) {
      auto& r = static_cast<FirstOrderType2R&>(*inter.relation());
      auto lambda = *inter.lambda(level);
      lambda -= *inter_work[siconos::integrators::EulerMoreauOSI::LAMBDAOLD];

      if (r.B())
        siconos::algebra::prod(
            *r.B(), lambda, *inter_work_block[siconos::integrators::EulerMoreauOSI::G_ALPHA],
            false);
      else
        siconos::algebra::prod(
            *relationMat[FirstOrderR::mat_B], lambda,
            *inter_work_block[siconos::integrators::EulerMoreauOSI::G_ALPHA], false);

      *DSlink[siconos::modeling::FirstOrderR::r] +=
          *inter_work_block[siconos::integrators::EulerMoreauOSI::G_ALPHA];
      DEBUG_EXPR(DSlink[FirstOrderR::r]->display(););
      // compute the new g_alpha

      r.computeg(time, *inter.lambda(level),
                 *inter_work_block[siconos::integrators::EulerMoreauOSI::G_ALPHA]);
      DEBUG_EXPR(inter_work_block[siconos::integrators::EulerMoreauOSI::G_ALPHA]->display(););
    }
    else if (relationSubType == NonLinearR) {
      auto& r = static_cast<siconos::modeling::FirstOrderNonLinearR&>(*inter.relation());
      // compute the new r  obtained by linearisation
      // r_{alpha+1}_{k+1} = g(lambda_{k+1}^{alpha},t_k+1)
      //                     + B_{k+1}^alpha ( lambda_{k+1}^{alpha+1}- lambda_{k+1}^{alpha} )

      auto lambda = *inter.lambda(level);
      lambda -= *inter_work[siconos::integrators::EulerMoreauOSI::LAMBDAOLD];

      // Remind that g_alpha has only one block
      auto& g_alpha =
          *(*inter_work_block[siconos::integrators::EulerMoreauOSI::G_ALPHA]).vector(0);

      if (r.B())
        siconos::algebra::prod(*r.B(), lambda, g_alpha, false);
      else
        siconos::algebra::prod(*relationMat[FirstOrderR::mat_B], lambda, g_alpha, false);

      auto& deltax = *inter_work_block[siconos::integrators::EulerMoreauOSI::DELTA_X];
      DEBUG_PRINT("FirstOrderNonLinearR::computeInput : deltax \n");
      DEBUG_EXPR(deltax.display());

      if (r.K())
        siconos::algebra::prod(*r.K(), deltax, g_alpha, false);
      else
        siconos::algebra::prod(*relationMat[FirstOrderR::mat_K], deltax, g_alpha, false);

      // Khat = h * K * W^-1 * B
      siconos::algebra::prod(*inter_work_mat[siconos::integrators::EulerMoreauOSI::MAT_KHAT],
                             *inter.lambda(level), g_alpha, false);

      *DSlink[siconos::modeling::FirstOrderR::r] += g_alpha;

      // compute the new g_alpha
      r.computeg(time, *DSlink[siconos::modeling::FirstOrderR::x], *inter.lambda(level),
                 *DSlink[siconos::modeling::FirstOrderR::z],
                 *inter_work_block[siconos::integrators::EulerMoreauOSI::G_ALPHA]);
    }
    else {
      inter.computeInput(time, level);
    }
  }
}

double siconos::integrators::EulerMoreauOSI::computeResiduOutput(
    double time, std::shared_ptr<siconos::graphs::InteractionsGraph> indexSet)
{
  double residu = 0.0;
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  for (std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui) {
    auto& inter_work = *indexSet->properties(*ui).workVectors;
    auto& residuY = *inter_work[siconos::integrators::EulerMoreauOSI::VEC_RESIDU_Y];
    auto& inter = *indexSet->bundle(*ui);
    residuY = *inter_work[siconos::integrators::EulerMoreauOSI::H_ALPHA];
    siconos::algebra::scal(-1, residuY, residuY);
    residuY += *(inter.y(0));
    DEBUG_EXPR(residuY.display(););
    residu = std::max(residu, residuY.norm2());
  }
  return residu;
}
double siconos::integrators::EulerMoreauOSI::computeResiduInput(
    double time, std::shared_ptr<siconos::graphs::InteractionsGraph> indexSet)
{
  double residu = 0.0;
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  for (std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui) {
    auto& interProp = indexSet->properties(*ui);
    auto& inter_work = *interProp.workVectors;
    auto& inter_work_block = *interProp.workBlockVectors;
    auto inter = indexSet->bundle(*ui);
    auto& DSlink = inter->linkToDSVariables();
    auto& residuR = *inter_work[siconos::integrators::EulerMoreauOSI::VEC_RESIDU_R];
    // Residu_r = r_alpha_k+1 - g_alpha;
    residuR = *DSlink[FirstOrderR::r];
    residuR -= *inter_work_block[siconos::integrators::EulerMoreauOSI::G_ALPHA];
    DEBUG_EXPR(residuR.display(););
    residu = std::max(residu, residuR.norm2());
  }
  return residu;
}
