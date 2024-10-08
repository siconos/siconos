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

#include "LsodarOSI.hpp"

#include "BlockVector.hpp"
#include "EventDriven.hpp"
#include "Interaction.hpp"
#include "LagrangianDS.hpp"
#include "LagrangianR.hpp"  // for LagrangianR::p2
#include "LagrangianRheonomousR.hpp"
#include "LagrangianScleronomousR.hpp"
#include "NewtonImpactNSL.hpp"  // For the visitor
#include "NonSmoothLaw.hpp"
#include "OneStepNSProblem.hpp"
#include "Relation.hpp"
#include "SiconosException.hpp"
#include "SiconosFortran.h"  // for lsodar
#include "SiconosMatrix.hpp"
#include "SiconosMatrixVectorOp.hpp"  // mat-vec prod
#include "SiconosVector.hpp"
#include "SiconosVectorOp.hpp"  // For subscal
#include "SiconosVisitor.hpp"   // for NSLEffectOnFreeOutput visitor
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

int siconos::integrators::LsodarOSI::count_NST = 0;
int siconos::integrators::LsodarOSI::count_NFE = 0;

// ===== Out of class objects and functions =====

namespace {  // Anonymous, local scope only

// global object and wrapping functions -> required for function plug-in and
// call in fortran routine.
std::shared_ptr<siconos::integrators::LsodarOSI> global_object{nullptr};

// This first function must have the same signature as argument F (arg 1) in
// DLSODAR (see opkdmain.f in Numerics) function to compute the righ-hand side
// of xdot = f(x,t) + Tu
extern "C" void LsodarOSI_f_wrapper(int* sizeOfX, double* time, double* x, double* xdot) {
  return global_object->f(sizeOfX, time, x, xdot);
}

// Function to wrap g: same signature as argument G (arg 18) in DLSODAR (see
// opkdmain.f in Numerics)
extern "C" void LsodarOSI_g_wrapper(int* nEq, double* time, double* x, int* ng, double* gOut) {
  return global_object->g(nEq, time, x, ng, gOut);
}

// Function to wrap jacobianf: same signature as argument JAC (arg 16) in
// DLSODAR (see opkdmain.f in Numerics) function to compute the Jacobian/x of
// the rhs.
extern "C" void LsodarOSI_jacobianf_wrapper(int* sizeOfX, double* time, double* x, int* ml,
                                            int* mu, double* jacob, int* nrowpd) {
  return global_object->jacobianfx(sizeOfX, time, x, ml, mu, jacob, nrowpd);
}
}  // namespace

siconos::integrators::LsodarOSI::LsodarOSI()
    : OneStepIntegrator(IntegratorType::LSODAROSI, 1, 0, 2, 1, 2) {
#if !defined(HAS_FORTRAN)
  THROW_EXCEPTION(
      "siconos::integrators::LsodarOSI: the fortran interface is not active. "
      "You can not "
      "create and use Lsodar Solver. Try to configure and reinstall siconos "
      "with "
      "WITH_FORTRAN=ON");

#endif
  _sizeMem = 2;
}

void siconos::integrators::LsodarOSI::setTol(int newItol, std::vector<double>&& newRtol,
                                             std::vector<double>&& newAtol) {
  //            The input parameters ITOL, RTOL, and ATOL determine
  //         the error control performed by the solver.  The solver will
  //         control the vector E = (E(i)) of estimated local errors
  //         in y, according to an inequality of the form
  //                     max-norm of ( E(i)/EWT(i) )   .le.   1,
  //         where EWT = (EWT(i)) is a vector of positive error weights.
  //         The values of RTOL and ATOL should all be non-negative.
  //         The following table gives the types (scalar/array) of
  //         RTOL and ATOL, and the corresponding form of EWT(i).
  //
  //            ITOL    RTOL       ATOL          EWT(i)
  //             1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
  //             2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
  //             3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
  //             4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)
  rtol = std::move(newRtol);
  atol = std::move(newAtol);
  _intData[2] = newItol;  // itol
}

void siconos::integrators::LsodarOSI::setMinMaxStepSizes(double minStep, double maxStep) {
  rwork[5] = minStep;
  rwork[6] = maxStep;
}

void siconos::integrators::LsodarOSI::setMaxNstep(int maxNumberSteps) {
  iwork[5] = maxNumberSteps;
}

void siconos::integrators::LsodarOSI::setTol(int newItol, double newRtol, double newAtol) {
  _intData[2] = newItol;  // itol
  rtol[0] = newRtol;      // rtol
  atol[0] = newRtol;      // atol
}

void siconos::integrators::LsodarOSI::setMaxOrder(int maxorderNonStiff, int maxorderStiff) {
  iwork[7] = maxorderNonStiff;
  iwork[8] = maxorderStiff;
}

void siconos::integrators::LsodarOSI::updateData() {
  // Used to update some data (iwork ...) when _intData is modified.
  // Warning: it only checks sizes and possibly reallocate memory, but no values
  // are set.

  iwork.resize(_intData[5], 0);

  // This is for documentation purposes only
  // Set the flag to generate extra printing at method switches.
  // iwork[4] = 0;
  // Set the maximal number of steps for one call
  // iwork[5] = 0;
  // set  the maximum number of messages printed (per problem)
  // iwork[6] = 0;
  // Set the maximum order to be allowed for the nonstiff (Adams) method
  // iwork[7] = 0;
  // Set   the maximum order to be allowed for the stiff  (BDF) method.
  // iwork[8] = 0;

  rwork.resize(_intData[4], 0.);

  jroot.resize(_intData[1], 0);
}

void siconos::integrators::LsodarOSI::fillXWork(int* sizeOfX, double* x) {
  assert((unsigned int)(*sizeOfX) == _xWork->size() &&
         "siconos::integrators::LsodarOSI::fillXWork xWork and sizeOfX have "
         "different sizes");
  (*_xWork) = x;
}

void siconos::integrators::LsodarOSI::computeRhs(double t) {
  DEBUG_BEGIN(
      "siconos::integrators::LsodarOSI::computeRhs(double t, "
      "DynamicalSystemsGraph& DSG0)\n")
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    // compute standard rhs stored in the dynamical system
    ds->computeRhs(t);
    DEBUG_EXPR(ds->getRhs().display(););
    /* This next line is a good protection  */
    assert(_dynamicalSystemsGraph->properties(*dsi).workVectors);
    auto& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
    if (auto lds = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      auto& free = *workVectors[siconos::integrators::LsodarOSI::FREE];
      // we assume that inverseMass and forces are updated after call of
      // ds->computeRhs(t);
      free = lds->totalForces();
      if (lds->LUMass()) {
        // lds->update_lu_mass();
        lds->LUMass()->solve(free);
      }
      DEBUG_EXPR(free.display(););
    }
    if (_extraAdditionalTerms) {
      auto dsgVD = _dynamicalSystemsGraph->descriptor(ds);
      _extraAdditionalTerms->addSmoothTerms(*_dynamicalSystemsGraph, dsgVD, t, ds->getRhs());
    }
  }
  DEBUG_END(
      "siconos::integrators::LsodarOSI::computeRhs(double t, "
      "DynamicalSystemsGraph& DSG0)\n")
}

void siconos::integrators::LsodarOSI::computeJacobianRhs(
    double t, siconos::graphs::DynamicalSystemsGraph& DSG0) {
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    auto ds = _dynamicalSystemsGraph->bundle(*dsi);
    ds->computeJacobianRhsx(t);
    if (_extraAdditionalTerms) {
      auto dsgVD = DSG0.descriptor(ds);
      _extraAdditionalTerms->addJacobianRhsContribution(
          DSG0, dsgVD, t, *(ds->jacobianRhsx()->toSiconosMatrix()));
    }
  }
}

void siconos::integrators::LsodarOSI::f(int* sizeOfX, double* time, double* x, double* xdot) {
  std::static_pointer_cast<siconos::simulation::EventDriven>(_simulation)
      ->computef(*this, sizeOfX, time, x, xdot);
}

void siconos::integrators::LsodarOSI::g(int* nEq, double* time, double* x, int* ng,
                                        double* gOut) {
  std::static_pointer_cast<siconos::simulation::EventDriven>(_simulation)
      ->computeg(shared_from_this(), nEq, time, x, ng, gOut);
}

void siconos::integrators::LsodarOSI::jacobianfx(int* sizeOfX, double* time, double* x,
                                                 int* ml, int* mu, double* jacob,
                                                 int* nrowpd) {
  std::static_pointer_cast<siconos::simulation::EventDriven>(_simulation)
      ->computeJacobianfx(*this, sizeOfX, time, x, jacob);
}

void siconos::integrators::LsodarOSI::initializeWorkVectorsForDS(
    double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  DEBUG_BEGIN(
      "siconos::integrators::LsodarOSI::initializeWorkVectorsForDS( double t, "
      "std::shared_ptr<siconos::modeling::DynamicalSystem> ds)\n");
  // Get work buffers from the graph
  auto& ds_work_vectors = *_initializeDSWorkVectors(ds);

  ds->initRhs(t);  // This will create p[2] and other required vectors/buffers

  if (auto lds = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
    // TODO FP: use buffer in graph for xWork?
    if (!_xWork) _xWork = std::make_shared<siconos::algebra::BlockVector>();
    _xWork->insertPtr(lds->q());
    _xWork->insertPtr(lds->velocity());
    ds_work_vectors.resize(siconos::integrators::LsodarOSI::WORK_LENGTH);
    ds_work_vectors[siconos::integrators::LsodarOSI::FREE] =
        std::make_shared<siconos::algebra::SiconosVector>(lds->dimension());
  } else {
    if (!_xWork) _xWork = std::make_shared<siconos::algebra::BlockVector>();
    _xWork->insertPtr(ds->x());
  }
  ds->swapInMemory();

  // Update necessary data

  // 1 - Neq; x vector size.
  _intData[0] = _xWork->size();
  // 5 - lrw, size of rwork
  _intData[4] = 22 + _intData[0] * std::max(16, (int)_intData[0] + 9) + 3 * _intData[1];
  // 6 - liw, size of iwork
  _intData[5] = 20 + _intData[0];

  // memory allocation for double*, according to _intData values
  updateData();

  _xtmp = std::make_shared<siconos::algebra::SiconosVector>(_xWork->size());

  computeRhs(t);

  DEBUG_END(
      "siconos::integrators::LsodarOSI::initializeWorkVectorsForDS( double t, "
      "std::shared_ptr<siconos::modeling::DynamicalSystem> ds)\n");
}

void siconos::integrators::LsodarOSI::initializeWorkVectorsForInteraction(
    siconos::modeling::Interaction& inter, siconos::graphs::InteractionProperties& interProp,
    siconos::graphs::DynamicalSystemsGraph& DSG) {
  auto ds1 = interProp.source;
  auto ds2 = interProp.target;
  assert(ds1);
  assert(ds2);

  auto& DSlink = inter.linkToDSVariables();
  if (!interProp.workVectors) {
    interProp.workVectors =
        std::make_shared<std::vector<std::shared_ptr<siconos::algebra::SiconosVector>>>(
            siconos::integrators::LsodarOSI::WORK_INTERACTION_LENGTH);
  }

  if (!interProp.workBlockVectors) {
    interProp.workBlockVectors =
        std::make_shared<std::vector<std::shared_ptr<siconos::algebra::BlockVector>>>();
    interProp.workBlockVectors->resize(siconos::integrators::LsodarOSI::BLOCK_WORK_LENGTH);
  }

  auto& inter_work = *interProp.workVectors;
  auto& inter_work_block = *interProp.workBlockVectors;

  auto& relation = *inter.relation();
  auto relationType = relation.getType();

  inter_work[siconos::integrators::LsodarOSI::OSNSP_RHS] =
      std::make_shared<siconos::algebra::SiconosVector>(inter.dimension());

  auto& nslaw = *inter.nonSmoothLaw();

  auto nslType = siconos::types::type_value(nslaw);

  if (nslType == siconos::modeling::Type::NewtonImpactNSL ||
      nslType == siconos::modeling::Type::MultipleImpactNSL) {
    _levelMinForOutput = 0;
    _levelMaxForOutput = 2;
    _levelMinForInput = 1;
    _levelMaxForInput = 2;
  } else if (nslType == siconos::modeling::Type::NewtonImpactFrictionNSL) {
    _levelMinForOutput = 0;
    _levelMaxForOutput = 4;
    _levelMinForInput = 1;
    _levelMaxForInput = 2;
    THROW_EXCEPTION(
        "siconos::integrators::LsodarOSI::initializeWorkVectorsForInteraction  "
        "not yet "
        "implemented for nonsmooth law of type NewtonImpactFrictionNSL");
  } else
    THROW_EXCEPTION(
        "siconos::integrators::LsodarOSI::initializeWorkVectorsForInteraction "
        "not yet "
        "implemented  for this type of nonsmoothlaw");

  // Check if interations levels (i.e. y and lambda sizes) are compliant with
  // the current osi.
  _check_and_update_interaction_levels(inter);
  // Initialize/allocate memory buffers in interaction.
  inter.initializeMemory(_steps);

  /* allocate and set work vectors for the osi */

  if (!(checkOSI(DSG.descriptor(ds1)) && checkOSI(DSG.descriptor(ds2)))) {
    THROW_EXCEPTION(
        "siconos::integrators::LsodarOSI::initializeWorkVectorsForInteraction. "
        "The "
        "implementation is not correct for two different OSI for one "
        "interaction");
  }

  auto& workVds1 = *DSG.properties(DSG.descriptor(ds1)).workVectors;
  if (relationType == siconos::modeling::RelationType::Lagrangian) {
    auto& lds = *std::static_pointer_cast<siconos::modeling::LagrangianDS>(ds1);
    inter_work_block[siconos::integrators::LsodarOSI::xfree] =
        std::make_shared<siconos::algebra::BlockVector>();
    inter_work_block[siconos::integrators::LsodarOSI::xfree]->insertPtr(
        workVds1[siconos::integrators::LsodarOSI::FREE]);
    DSlink[siconos::modeling::LagrangianR::p2] =
        std::make_shared<siconos::algebra::BlockVector>();
    DSlink[siconos::modeling::LagrangianR::p2]->insertPtr(lds.p(2));
    DSlink[siconos::modeling::LagrangianR::q2] =
        std::make_shared<siconos::algebra::BlockVector>();
    DSlink[siconos::modeling::LagrangianR::q2]->insertPtr(lds.acceleration());
  }
  // else if (relationType == NewtonEuler)
  // {
  //   inter_work_block[::xfree] =
  //   std::make_shared<siconos::algebra::BlockVector>();
  //   inter_work_block[::xfree]->insertPtr(workVds1[siconos::integrators::LsodarOSI::FREE]);
  // }

  if (ds1 != ds2) {
    auto& workVds2 = *DSG.properties(DSG.descriptor(ds2)).workVectors;
    if (relationType == siconos::modeling::RelationType::Lagrangian) {
      auto& lds = *std::static_pointer_cast<siconos::modeling::LagrangianDS>(ds2);
      inter_work_block[siconos::integrators::LsodarOSI::xfree]->insertPtr(
          workVds2[siconos::integrators::LsodarOSI::FREE]);
      DSlink[siconos::modeling::LagrangianR::p2]->insertPtr(lds.p(2));
      DSlink[siconos::modeling::LagrangianR::q2]->insertPtr(lds.acceleration());
    }
    // else if (relationType == NewtonEuler)
    // {
    //   inter_work_block[NewtonEulerR::xfree]->insertPtr(workVds2[siconos::integrators::LsodarOSI::FREE]);
    // }
  }
}

void siconos::integrators::LsodarOSI::initialize() {
  DEBUG_BEGIN("siconos::integrators::LsodarOSI::initialize()\n");
  OneStepIntegrator::initialize();
  // std::string type;
  //  initialize xWork with x values of the dynamical systems present in the set.

  // siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  // for(std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi !=
  // dsend; ++dsi)
  //   {
  //     if(!checkOSI(dsi)) continue;
  //     std::shared_ptr<siconos::modeling::DynamicalSystem> ds =
  //     _dynamicalSystemsGraph->bundle(*dsi); initializeWorkVectorsForDS(m,
  //     m.t0(),ds);
  //     //ds->resetToInitialState();
  //     //ds->swapInMemory();
  //   }

  // auto indexSet0 = m.nonSmoothDynamicalSystem()->topology()->indexSet0();
  // siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  // for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  //   {
  //     auto& inter = *indexSet0->bundle(*ui);
  //     initializeWorkVectorsForInteraction(m.t0(), inter,
  //     indexSet0->properties(*ui),
  //     *_dynamicalSystemsGraph);
  //   }

  //   int parameters for LSODAROSI are saved in vector intParam.
  //   The link with variable names in opkdmain.f is indicated in comments

  // 2 - Ng, number of constraints:
  _intData[1] = std::static_pointer_cast<siconos::simulation::EventDriven>(_simulation)
                    ->computeSizeOfg();
  // 3 - Itol, itask, iopt
  // intData[2,3,4,5] : default values set in class attribute
  // _intData[2] = 1 if ATOL is a scalar, else 2 (ATOL array)

  // 4 - Istate
  _intData[3] = 1;  // istate, an index used for input and output to specify the
                    // state of the calculation.
  // On input:
  //                 1: first call for the problem (initializations will be done).
  //                 2: means this is not the first call, and the calculation is to continue
  //                 normally, with no change in any input
  //                    parameters except possibly TOUT and ITASK.
  //                 3:  means this is not the first call, and the calculation is to continue
  //                 normally, but with
  //                     a change in input parameters other than TOUT and ITASK.
  // On output:
  //                 1: means nothing was done; TOUT = t and ISTATE = 1 on input.
  //                 2: means the integration was performed successfully, and no roots were
  //                 found. 3: means the integration was successful, and one or more roots
  //                 were found before satisfying the stop condition specified by ITASK. See
  //                 JROOT. <0: error. See table below, in integrate function output message.

  // 7 - JT, Jacobian type indicator
  _intData[6] = 1;  // jt, Jacobian type indicator.
  //           1 means a user-supplied full (NEQ by NEQ) Jacobian.
  //           2 means an internally generated (difference quotient) full
  //           Jacobian (using NEQ extra calls to f per df/dx value). 4 means a
  //           user-supplied banded Jacobian. 5 means an internally generated
  //           banded Jacobian (using ML+MU+1 extra calls to f per df/dx
  //           evaluation).

  // set the optional input flags of LSODAROSI to 0
  // LSODAROSI will take the default values

  // === Error handling in LSODAROSI===

  //   parameters: itol, rtol, atol.
  //   Control vector E = (E(i)) of estimated local errors in y:
  //   max-norm of ( E(i)/EWT(i) )< 1
  //   EWT = (EWT(i)) vector of positive error weights.
  //   The values of RTOL and ATOL should all be non-negative.
  //
  //  ITOL    RTOL       ATOL          EWT(i)
  //   1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
  //   2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
  //   3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
  //   4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)
  DEBUG_END("siconos::integrators::LsodarOSI::initialize()\n");
}

void siconos::integrators::LsodarOSI::integrate(double& tinit, double& tend, double& tout,
                                                int& istate) {
  DEBUG_BEGIN(
      "siconos::integrators::LsodarOSI::integrate(double& tinit, double& tend, "
      "double& tout, "
      "int& istate) with \n");
  DEBUG_PRINTF("tinit = %f, tend= %f, tout = %f, istate = %i\n", tinit, tend, tout, istate);

  // For details on DLSODAR parameters, see opkdmain.f in externals/odepack
  double tend_DR = tend;    // next point where output is desired (different from t!)
  double tinit_DR = tinit;  // current (starting) time

  // === Pointers to function ===
  //  --> definition and initialisation thanks to wrapper:
  global_object = std::static_pointer_cast<LsodarOSI>(
      shared_from_this());  // Warning: global object must be initialized to
                            // current one before pointers to function
                            // initialisation.

  // === LSODAR CALL ===
  DEBUG_EXPR(_xWork->display(););
  *_xtmp = *(_xWork->toSiconosVector());
  if (istate == 3) {
    istate = 1;  // restart TEMPORARY
  }

  _intData[3] = istate;

  // call LSODAR to integrate dynamical equation
  siconos::netlib::lsodar(&LsodarOSI_f_wrapper, &(_intData[0]), _xtmp->data(), &tinit_DR,
                          &tend_DR, &(_intData[2]), &rtol.front(), &atol.front(),
                          &(_intData[3]), &rwork.front(), &(_intData[4]), &iwork.front(),
                          &(_intData[5]), &LsodarOSI_jacobianf_wrapper, &(_intData[6]),
                          &LsodarOSI_g_wrapper, &(_intData[1]), &jroot.front());

  // jroot: jroot[i] = 1 if g(i) has a root at t, else jroot[i] = 0.

  // === Post ===
  if (_intData[3] < 0)  // if istate < 0 => LSODAROSI failed
  {
    std::cout << "Lsodar::integrate(...) failed - Istate = " << _intData[3] << std::endl;
    std::cout << " -1 means excess work done on this call (perhaps wrong JT, or so "
                 "small "
                 "tolerance (ATOL and RTOL), or small maximum number of steps for "
                 "one call "
                 "(MXSTEP)). You should increase ATOL or RTOL or increase the MXSTEP"
              << std::endl;
    std::cout << " -2 means excess accuracy requested (tolerances too small)." << std::endl;
    std::cout << " -3 means illegal input detected (see printed message)." << std::endl;
    std::cout << " -4 means repeated error test failures (check all inputs)." << std::endl;
    std::cout << " -5 means repeated convergence failures (perhaps bad "
                 "Jacobian supplied or "
                 "wrong choice of JT or tolerances)."
              << std::endl;
    std::cout << " -6 means error weight became zero during problem. (Solution "
                 "component i "
                 "vanished, and ATOL or ATOL(i) = 0.)"
              << std::endl;
    std::cout << " -7 means work space insufficient to finish (see messages)." << std::endl;
    THROW_EXCEPTION("LsodarOSI, integration failed");
  }

  *_xWork = *_xtmp;
  istate = _intData[3];
  tout = tinit_DR;  // real ouput time
  tend = tend_DR;   // necessary for next start of DLSODAR
  DEBUG_PRINTF("tout = %g, tinit = %g, tend = %g ", tout, tinit, tend);
  DEBUG_EXPR(_xtmp->display(););
  DEBUG_EXPR(_xWork->display(););
  if (istate == 3) {
    //      std:: std::cout << "ok\n";
    assert(true);
  }
  // Update counters
  count_NST = iwork[10];
  count_NFE = iwork[11];
  //  tinit = tinit_DR;
  DEBUG_END("LsodarOSI::integrate(double& tinit, double& tend, double& tout, int& istate)\n");
}

void siconos::integrators::LsodarOSI::updateState(const unsigned int level) {
  // Compute all required (ie time-dependent) data for the DS of the OSI.
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  if (level == 1)  // ie impact case: compute velocity
  {
    for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
      if (!checkOSI(dsi)) continue;
      auto lds = std::static_pointer_cast<siconos::modeling::LagrangianDS>(
          _dynamicalSystemsGraph->bundle(*dsi));
      lds->computePostImpactVelocity();
    }
  } else if (level == 2) {
    auto time = _simulation->nextTime();
    for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
      if (!checkOSI(dsi)) continue;
      {
        auto ds = _dynamicalSystemsGraph->bundle(*dsi);
        ds->update(time);
      }
    }
  } else
    THROW_EXCEPTION(
        "siconos::integrators::LsodarOSI::updateState(index), index is out of "
        "range. Index "
        "= " +
        std::to_string(level));
}

struct siconos::integrators::LsodarOSI::_NSLEffectOnFreeOutput
    : public siconos::internal::SiconosVisitor {
  using SiconosVisitor::visit;

  siconos::nonsmooth_formulations::OneStepNSProblem& _osnsp;
  std::shared_ptr<siconos::modeling::Interaction> _inter{nullptr};
  siconos::graphs::InteractionProperties& _interProp;

  _NSLEffectOnFreeOutput(siconos::nonsmooth_formulations::OneStepNSProblem& p,
                         std::shared_ptr<siconos::modeling::Interaction> inter,
                         siconos::graphs::InteractionProperties& interProp)
      : _osnsp(p), _inter(inter), _interProp(interProp){};

  _NSLEffectOnFreeOutput(const _NSLEffectOnFreeOutput&) = delete;
  _NSLEffectOnFreeOutput(const _NSLEffectOnFreeOutput&&) = delete;
  _NSLEffectOnFreeOutput& operator=(const _NSLEffectOnFreeOutput&) = delete;
  _NSLEffectOnFreeOutput& operator=(const _NSLEffectOnFreeOutput&&) = delete;

  void visit(const siconos::modeling::NewtonImpactNSL& nslaw) const override {
    double e;
    e = nslaw.e();
    std::vector<size_t> subCoord(4);
    subCoord[0] = 0;
    subCoord[1] = _inter->nonSmoothLaw()->size();
    subCoord[2] = 0;
    subCoord[3] = subCoord[1];
    auto& osnsp_rhs = *(*_interProp.workVectors)[siconos::integrators::LsodarOSI::OSNSP_RHS];
    siconos::algebra::subscal(e, _inter->y_k(_osnsp.inputOutputLevel()), osnsp_rhs, subCoord,
                              false);  // q = q + e * q
  }

  // visit function added by Son (9/11/2010)
  void visit(const siconos::modeling::MultipleImpactNSL& nslaw) const override { ; }

  // note : no NewtonImpactFrictionNSL
};

void siconos::integrators::LsodarOSI::computeFreeOutput(
    siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter,
    siconos::nonsmooth_formulations::OneStepNSProblem* osnsp) {
  auto allOSNS = _simulation->oneStepNSProblems();
  auto indexSet = osnsp->simulation()->indexSet(osnsp->indexSetLevel());
  auto inter = indexSet->bundle(vertex_inter);
  auto& DSlink = inter->linkToDSVariables();
  auto& inter_work_block = *indexSet->properties(vertex_inter).workBlockVectors;

  // Get relation and non smooth law types
  auto relationType = inter->relation()->getType();
  auto relationSubType = inter->relation()->getSubType();
  unsigned int sizeY = inter->nonSmoothLaw()->size();

  unsigned int relativePosition = 0;
  auto mainInteraction = inter;
  std::vector<size_t> coord(8);
  coord[0] = relativePosition;
  coord[1] = relativePosition + sizeY;
  coord[2] = 0;
  coord[4] = 0;
  coord[6] = 0;
  coord[7] = sizeY;
  std::shared_ptr<siconos::algebra::SiconosMatrix> C;
  //   std::shared_ptr<siconos::algebra::SiconosMatrix>  D;
  //   std::shared_ptr<siconos::algebra::SiconosMatrix>  F;
  auto& osnsp_rhs = *(*indexSet->properties(vertex_inter)
                           .workVectors)[siconos::integrators::LsodarOSI::OSNSP_RHS];

  std::shared_ptr<siconos::algebra::BlockVector> Xfree;

  // All of these values should be stored in the node corrseponding to the
  // Interactionwhen a LsodarOSI scheme is used.

  /* V.A. 10/10/2010
   * Following the type of OSNS  we need to retrieve the velocity or the
   * acceleration This tricks is not very nice but for the moment the OSNS do
   * not known if it is in accelaration of not
   */

  // auto  allOSNS  = _simulation->oneStepNSProblems();
  if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_ACC]).get() == osnsp) {
    if (relationType == siconos::modeling::RelationType::Lagrangian) {
      Xfree = inter_work_block[siconos::integrators::LsodarOSI::xfree];
      DEBUG_EXPR(Xfree->display(););
    }
    // else if  (relationType == siconos::modeling::RelationType::NewtonEuler)
    // {
    //   Xfree = inter->data(::FREE);
    // }
    assert(Xfree);
    //        std::cout << "Computeqblock Xfree (Gamma)========" << std::endl;
    //       Xfree->display();
  } else if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_ED_IMPACT]).get() == osnsp) {
    Xfree = DSlink[siconos::modeling::LagrangianR::q1];
    //        std::cout << "Computeqblock Xfree (Velocity)========" <<
    //        std::endl;
    //       Xfree->display();
  } else
    THROW_EXCEPTION(" computeqBlock for Event Event-driven is wrong ");

  if (relationType == siconos::modeling::RelationType::Lagrangian) {
    C = mainInteraction->relation()->C();
    if (C) {
      assert(Xfree);

      coord[3] = C->size(1);
      coord[5] = C->size(1);

      siconos::algebra::subprod(*C, *Xfree, osnsp_rhs, coord, true);
    }

    auto ID = std::make_shared<siconos::algebra::SiconosMatrix>(sizeY, sizeY);
    ID->eye();

    std::vector<size_t> xcoord(8);
    xcoord[0] = 0;
    xcoord[1] = sizeY;
    xcoord[2] = 0;
    xcoord[3] = sizeY;
    xcoord[4] = 0;
    xcoord[5] = sizeY;
    xcoord[6] = 0;
    xcoord[7] = sizeY;
    // For the relation of type LagrangianRheonomousR
    if (relationSubType == siconos::modeling::RelationSubType::RheonomousR) {
      if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_ACC]).get() == osnsp) {
        THROW_EXCEPTION(
            "siconos::integrators::LsodarOSI::computeFreeOutput not yet "
            "implemented for LCP "
            "at acceleration level with LagrangianRheonomousR");
      } else if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp) {
        std::static_pointer_cast<siconos::modeling::LagrangianRheonomousR>(inter->relation())
            ->computehDot(simulation()->getTkp1(), *DSlink[siconos::modeling::LagrangianR::q0],
                          *DSlink[siconos::modeling::LagrangianR::z]);
        siconos::algebra::subprod(
            *ID,
            *(std::static_pointer_cast<siconos::modeling::LagrangianRheonomousR>(
                  inter->relation())
                  ->hDot()),
            osnsp_rhs, xcoord, false);  // y += hDot
      } else
        THROW_EXCEPTION(
            "siconos::integrators::LsodarOSI::computeFreeOutput not "
            "implemented for "
            "SICONOS_OSNSP ");
    }
    // For the relation of type LagrangianScleronomousR
    if (relationSubType == siconos::modeling::RelationSubType::ScleronomousR) {
      if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_ACC]).get() == osnsp) {
        std::static_pointer_cast<siconos::modeling::LagrangianScleronomousR>(inter->relation())
            ->computedotjacqhXqdot(simulation()->getTkp1(), *inter);
        siconos::algebra::subprod(
            *ID,
            *(std::static_pointer_cast<siconos::modeling::LagrangianScleronomousR>(
                  inter->relation())
                  ->dotjacqhXqdot()),
            osnsp_rhs, xcoord, false);  // y += NonLinearPart
      }
    }
  } else
    THROW_EXCEPTION(
        "siconos::integrators::LsodarOSI::computeFreeOutput not yet "
        "implemented for Relation "
        "of type " +
        std::to_string(
            static_cast<std::underlying_type<siconos::modeling::RelationType>::type>(
                relationType)));
  if (((*allOSNS)[siconos::simulation::SICONOS_OSNSP_ED_IMPACT]).get() == osnsp) {
    if (inter->relation()->getType() == siconos::modeling::RelationType::Lagrangian ||
        inter->relation()->getType() == siconos::modeling::RelationType::NewtonEuler) {
      _NSLEffectOnFreeOutput nslEffectOnFreeOutput(*osnsp, inter,
                                                   indexSet->properties(vertex_inter));

      inter->nonSmoothLaw()->accept(nslEffectOnFreeOutput);
    }
  }
}
void siconos::integrators::LsodarOSI::display() const {
  OneStepIntegrator::display();
  std::cout << " --- > LsodarOSI specific values: \n";
  std::cout << "Number of equations: " << _intData[0] << "\n";
  std::cout << "Number of constraints: " << _intData[1] << "\n";
  std::cout << "itol, istate, lrw, liw, jt: (for details on what are these "
               "variables see opkdmain.f)\n";
  std::cout << _intData[2] << ", " << _intData[3] << ", " << _intData[4] << ", " << _intData[5]
            << ", " << _intData[6] << "\n";
  std::cout << "====================================\n";
}
