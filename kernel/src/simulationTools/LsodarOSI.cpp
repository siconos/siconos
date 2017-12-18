
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

#include "LsodarOSI.hpp"
#include "EventDriven.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "BlockVector.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Model.hpp"
#include "Topology.hpp"
#include "LagrangianRheonomousR.hpp"
#include "LagrangianScleronomousR.hpp"
#include "NewtonImpactNSL.hpp"
#include "MultipleImpactNSL.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "FirstOrderNonLinearDS.hpp"
#include "ExtraAdditionalTerms.hpp"
#include "OneStepNSProblem.hpp"
#include "TypeName.hpp"
#include <odepack.h>

using namespace RELATION;

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"

int LsodarOSI::count_NST = 0;
int LsodarOSI::count_NFE = 0;

// ===== Out of class objects and functions =====

// global object and wrapping functions -> required for function plug-in and call in fortran routine.
SP::LsodarOSI global_object;

// This first function must have the same signature as argument F (arg 1) in DLSODAR (see opkdmain.f in Numerics)
extern "C" void LsodarOSI_f_wrapper(integer* sizeOfX, doublereal* time, doublereal* x, doublereal* xdot);
extern "C" void LsodarOSI_f_wrapper(integer* sizeOfX, doublereal* time, doublereal* x, doublereal* xdot)
{
  return global_object->f(sizeOfX, time, x, xdot);
}

// Function to wrap g: same signature as argument G (arg 18) in DLSODAR (see opkdmain.f in Numerics)
extern "C" void LsodarOSI_g_wrapper(integer* nEq, doublereal* time, doublereal* x, integer* ng, doublereal* gOut);
extern "C" void LsodarOSI_g_wrapper(integer* nEq, doublereal* time, doublereal* x, integer* ng, doublereal* gOut)
{
  return global_object->g(nEq, time, x, ng, gOut);
}

// Function to wrap jacobianf: same signature as argument JAC (arg 16) in DLSODAR (see opkdmain.f in Numerics)
extern "C" void LsodarOSI_jacobianf_wrapper(integer* sizeOfX, doublereal* time, doublereal* x, integer* ml, integer* mu,  doublereal* jacob, integer* nrowpd);
extern "C" void LsodarOSI_jacobianf_wrapper(integer* sizeOfX, doublereal* time, doublereal* x, integer* ml, integer* mu,  doublereal* jacob, integer* nrowpd)
{
  return global_object->jacobianfx(sizeOfX, time, x, ml, mu, jacob, nrowpd);
}

LsodarOSI::LsodarOSI():
  OneStepIntegrator(OSI::LSODAROSI)
{
  _intData.resize(9);
  for(int i = 0; i < 9; i++) _intData[i] = 0;
  _sizeMem = 2;
  _steps=1;

  // Set levels. This may depend on the nonsmooth law and will be updated during fillDSLinks(...) call.
  _levelMinForOutput=0;
  _levelMaxForOutput=2;
  _levelMinForInput=1;
  _levelMaxForInput=2;

}

void LsodarOSI::setTol(integer newItol, SA::doublereal newRtol, SA::doublereal newAtol)
{
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

  _intData[2] = newItol; // itol

  rtol = newRtol;
  atol = newAtol;
}

void LsodarOSI::setMinMaxStepSizes(doublereal minStep, doublereal maxStep)
{
  _intData[5] = 1; // set IOPT = 1
  rwork[5] = minStep;
  rwork[6] = maxStep;
}

void LsodarOSI::setMaxNstep(integer maxNumberSteps)
{
  _intData[5] = 1; // set IOPT = 1
  iwork[5] = maxNumberSteps;
}

void LsodarOSI::setTol(integer newItol, doublereal newRtol, doublereal newAtol)
{
  _intData[2] = newItol; // itol
  rtol[0] = newRtol; // rtol
  atol[0] = newRtol;  // atol
}

void LsodarOSI::setMaxOrder(integer maxorderNonStiff, integer maxorderStiff)
{
  _intData[5] = 1; // set IOPT = 1
  iwork[7] = maxorderNonStiff;
  iwork[8] = maxorderStiff;
}

void LsodarOSI::updateData()
{
  // Used to update some data (iwork ...) when _intData is modified.
  // Warning: it only checks sizes and possibly reallocate memory, but no values are set.

  unsigned int sizeTol = _intData[0]; // size of rtol, atol ... If itol (_intData[0]) = 1 => scalar else, vector of size neq (_intData[0]).
  //  if(_intData[0]==1) sizeTol = 1;
  //  else sizeTol = _intData[0];

  rtol.reset(new doublereal[sizeTol]) ;    // rtol, relative tolerance

  atol.reset(new doublereal[sizeTol]) ;  // atol, absolute tolerance
  for(unsigned int i = 0; i < sizeTol; i++)
  {
    atol[i] = 0.0;
  }



  iwork.reset(new integer[_intData[7]]);
  for(int i = 0; i < _intData[7]; i++) iwork[i] = 0;

  rwork.reset(new doublereal[_intData[6]]);
  for(int i = 0; i < _intData[6]; i++) rwork[i] = 0.0;

  jroot.reset(new integer[_intData[1]]);
  for(int i = 0; i < _intData[1]; i++) jroot[i] = 0;

}

void LsodarOSI::fillXWork(integer* sizeOfX, doublereal* x)
{
  assert((unsigned int)(*sizeOfX) == _xWork->size() && "LsodarOSI::fillXWork xWork and sizeOfX have different sizes");
  (*_xWork) = x;
}

void LsodarOSI::computeRhs(double t, DynamicalSystemsGraph& DSG0)
{
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    // compute standard rhs stored in the dynamical system
    ds->computeRhs(t);
    VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
    Type::Siconos dsType = Type::value(*ds);
    if(dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianDS)
    {
      SP::LagrangianDS lds = std11::static_pointer_cast<LagrangianDS> (ds);
      SiconosVector &free=*workVectors[OneStepIntegrator::free];
      // we assume that inverseMass and forces are updated after call of ds->computeRhs(t);
      free = *lds->forces();
      if(lds->inverseMass())
        lds->inverseMass()->PLUForwardBackwardInPlace(free);
    }
    if(_extraAdditionalTerms)
    {
      DynamicalSystemsGraph::VDescriptor dsgVD = DSG0.descriptor(ds);
      _extraAdditionalTerms->addSmoothTerms(DSG0, dsgVD, t, ds->getRhs());
    }
  }
}

void LsodarOSI::computeJacobianRhs(double t, DynamicalSystemsGraph& DSG0)
{

  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    ds->computeJacobianRhsx(t);
    if(_extraAdditionalTerms)
    {
      DynamicalSystemsGraph::VDescriptor dsgVD = DSG0.descriptor(ds);
      _extraAdditionalTerms->addJacobianRhsContribution(DSG0, dsgVD, t, *(ds->jacobianRhsx()));
    }
  }
}

void LsodarOSI::f(integer* sizeOfX, doublereal* time, doublereal* x, doublereal* xdot)
{
  std11::static_pointer_cast<EventDriven>(_simulation)->computef(*this, sizeOfX, time, x, xdot);
}

void LsodarOSI::g(integer* nEq, doublereal*  time, doublereal* x, integer* ng, doublereal* gOut)
{
  std11::static_pointer_cast<EventDriven>(_simulation)->computeg(shared_from_this(), nEq, time, x, ng, gOut);
}

void LsodarOSI::jacobianfx(integer* sizeOfX, doublereal* time, doublereal* x, integer* ml, integer* mu,  doublereal* jacob, integer* nrowpd)
{
  std11::static_pointer_cast<EventDriven>(_simulation)->computeJacobianfx(*this, sizeOfX, time, x, jacob);
}


void LsodarOSI::initializeDynamicalSystem(Model& m, double t, SP::DynamicalSystem ds)
{
  DEBUG_BEGIN("LsodarOSI::initializeDynamicalSystem(Model& m, double t, SP::DynamicalSystem ds)\n");
  // Get work buffers from the graph
  VectorOfVectors& workVectors = *_initializeDSWorkVectors(ds);

  Type::Siconos dsType = Type::value(*ds);

  ds->initRhs(t); // This will create p[2] and other required vectors/buffers

  if(dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
  {
    LagrangianDS& lds = *std11::static_pointer_cast<LagrangianDS>(ds);
    // TODO FP: use buffer in graph for xWork?
    _xWork->insertPtr(lds.q());
    _xWork->insertPtr(lds.velocity());
    workVectors.resize(OneStepIntegrator::work_vector_of_vector_size);
    workVectors[OneStepIntegrator::free].reset(new SiconosVector(lds.dimension()));
  }
  else
    _xWork->insertPtr(ds->x());

  ds->swapInMemory();

  DEBUG_END("LsodarOSI::initializeDynamicalSystem(Model& m, double t, SP::DynamicalSystem ds)\n");
}

void LsodarOSI::fillDSLinks(Interaction &inter,
                              InteractionProperties& interProp,
                              DynamicalSystemsGraph & DSG)
{
  SP::DynamicalSystem ds1= interProp.source;
  SP::DynamicalSystem ds2= interProp.target;


  VectorOfVectors& workV = *interProp.workVectors;
  workV.resize(LsodarOSI::WORK_INTERACTION_LENGTH);
  workV[LsodarOSI::OSNSP_RHS].reset(new SiconosVector(inter.getSizeOfY()));

  assert(interProp.DSlink);
  VectorOfBlockVectors& DSlink = *interProp.DSlink;
  // VectorOfVectors& workVInter = *interProp.workVectors;
  // VectorOfSMatrices& workMInter = *interProp.workMatrices;

  Relation &relation =  *inter.relation();
  NonSmoothLaw & nslaw = *inter.nonSmoothLaw();
  RELATION::TYPES relationType = relation.getType();
  Type::Siconos nslType = Type::value(nslaw);

  if (nslType == Type::NewtonImpactNSL || nslType == Type::MultipleImpactNSL)
  {
    _levelMinForOutput = 0;
    _levelMaxForOutput = 2 ;
    _levelMinForInput = 1;
    _levelMaxForInput = 2;
  }
  else if (nslType ==  Type::NewtonImpactFrictionNSL)
  {
    _levelMinForOutput = 0;
    _levelMaxForOutput = 4;
    _levelMinForInput = 1;
    _levelMaxForInput = 2;
    RuntimeException::selfThrow("LsodarOSI::fillDSLinks  not yet implemented for nonsmooth law of type NewtonImpactFrictionNSL");
  }
  else
    RuntimeException::selfThrow("LsodarOSI::fillDSLinks not yet implemented  for nonsmooth of type");

  // Check if interations levels (i.e. y and lambda sizes) are compliant with the current osi.
  _check_and_update_interaction_levels(inter);
  // Initialize/allocate memory buffers in interaction.
  bool computeResidu = relation.requireResidu();
  inter.initializeMemory(computeResidu,_steps);

  /* allocate and set work vectors for the osi */

  if (!(checkOSI(DSG.descriptor(ds1)) && checkOSI(DSG.descriptor(ds2))))
  {
    RuntimeException::selfThrow("LsodarOSI::fillDSLinks. The implementation is not correct for two different OSI for one interaction");
  }



  VectorOfVectors &workVds1 = *DSG.properties(DSG.descriptor(ds1)).workVectors;
  if (relationType == Lagrangian)
  {
    LagrangianDS& lds = *std11::static_pointer_cast<LagrangianDS> (ds1);
    DSlink[LagrangianR::xfree].reset(new BlockVector());
    DSlink[LagrangianR::xfree]->insertPtr(workVds1[OneStepIntegrator::free]);
    DSlink[LagrangianR::p2].reset(new BlockVector());
    DSlink[LagrangianR::p2]->insertPtr(lds.p(2));
    DSlink[LagrangianR::q2].reset(new BlockVector());
    DSlink[LagrangianR::q2]->insertPtr(lds.acceleration());
  }
  // else if (relationType == NewtonEuler)
  // {
  //   DSlink[NewtonEulerR::xfree].reset(new BlockVector());
  //   DSlink[NewtonEulerR::xfree]->insertPtr(workVds1[OneStepIntegrator::free]);
  // }



  if (ds1 != ds2)
  {
    VectorOfVectors &workVds2 = *DSG.properties(DSG.descriptor(ds2)).workVectors;
    if (relationType == Lagrangian)
    {
      LagrangianDS& lds = *std11::static_pointer_cast<LagrangianDS> (ds2);
      DSlink[LagrangianR::xfree]->insertPtr(workVds2[OneStepIntegrator::free]);
      DSlink[LagrangianR::p2]->insertPtr(lds.p(2));
      DSlink[LagrangianR::q2]->insertPtr(lds.acceleration());
    }
    // else if (relationType == NewtonEuler)
    // {
    //   DSlink[NewtonEulerR::xfree]->insertPtr(workVds2[OneStepIntegrator::free]);
    // }
  }
}

void LsodarOSI::initialize(Model& m)
{
  DEBUG_BEGIN("LsodarOSI::initialize(Model& m)\n");
  _xWork.reset(new BlockVector());
  OneStepIntegrator::initialize(m);
  //std::string type;
  // initialize xWork with x values of the dynamical systems present in the set.

  // DynamicalSystemsGraph::VIterator dsi, dsend;
  // for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  //   {
  //     if(!checkOSI(dsi)) continue;
  //     SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
  //     initializeDynamicalSystem(m, m.t0(),ds);
  //     //ds->resetToInitialState();
  //     //ds->swapInMemory();
  //   }

  // SP::InteractionsGraph indexSet0 = m.nonSmoothDynamicalSystem()->topology()->indexSet0();
  // InteractionsGraph::VIterator ui, uiend;
  // for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  //   {
  //     Interaction& inter = *indexSet0->bundle(*ui);
  //     fillDSLinks(m.t0(), inter, indexSet0->properties(*ui), *_dynamicalSystemsGraph);
  //   }

  computeRhs(m.t0(),*_dynamicalSystemsGraph);


  //   Integer parameters for LSODAROSI are saved in vector intParam.
  //   The link with variable names in opkdmain.f is indicated in comments

  // 1 - Neq; x vector size.
  _intData[0] = _xWork->size();
  _xtmp.reset(new SiconosVector(_xWork->size()));

  // 2 - Ng, number of constraints:
  _intData[1] = std11::static_pointer_cast<EventDriven>(_simulation)->computeSizeOfg();
  // 3 - Itol, itask, iopt
  _intData[2] = 1; // itol, 1 if ATOL is a scalar, else 2 (ATOL array)
  _intData[3] = 1; // itask, an index specifying the task to be performed. 1: normal computation.
  _intData[5] = 0; // iopt: 0 if no optional input else 1.

  // 4 - Istate
  _intData[4] = 1; // istate, an index used for input and output to specify the state of the calculation.
  // On input:
  //                 1: first call for the problem (initializations will be done).
  //                 2: means this is not the first call, and the calculation is to continue normally, with no change in any input
  //                    parameters except possibly TOUT and ITASK.
  //                 3:  means this is not the first call, and the calculation is to continue normally, but with
  //                     a change in input parameters other than TOUT and ITASK.
  // On output:
  //                 1: means nothing was done; TOUT = t and ISTATE = 1 on input.
  //                 2: means the integration was performed successfully, and no roots were found.
  //                 3: means the integration was successful, and one or more roots were found before satisfying the stop condition specified by ITASK. See JROOT.
  //                 <0: error. See table below, in integrate function output message.


  // 5 - lrw, size of rwork
  _intData[6] = 22 + _intData[0] * std::max(16, (int)_intData[0] + 9) + 3 * _intData[1];

  // 6 - liw, size of iwork
  _intData[7] = 20 + _intData[0];

  // 7 - JT, Jacobian type indicator
  _intData[8] = 2;   // jt, Jacobian type indicator.
  //           1 means a user-supplied full (NEQ by NEQ) Jacobian.
  //           2 means an internally generated (difference quotient) full Jacobian (using NEQ extra calls to f per df/dx value).
  //           4 means a user-supplied banded Jacobian.
  //           5 means an internally generated banded Jacobian (using ML+MU+1 extra calls to f per df/dx evaluation).

  // memory allocation for doublereal*, according to _intData values ...
  updateData();

  // set the optional input flags of LSODAROSI to 0
  // LSODAROSI will take the default values

  // Set the flag to generate extra printing at method switches.
  iwork[4] = 0;
  // Set the maximal number of steps for one call
  iwork[5] = 0;
  // set  the maximum number of messages printed (per problem)
  iwork[6] = 0;
  // Set the maximum order to be allowed for the nonstiff (Adams) method
  iwork[7] = 0;
  // Set   the maximum order to be allowed for the stiff  (BDF) method.
  iwork[8] = 0;
  // Set atol and rtol values ...
  rtol[0] = RTOL_DEFAULT ; // rtol
  atol[0] = ATOL_DEFAULT ;  // atol

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
  DEBUG_END("LsodarOSI::initialize(Model& m)\n");
}

void LsodarOSI::integrate(double& tinit, double& tend, double& tout, int& istate)
{

  DEBUG_PRINT("LsodarOSI::integrate(double& tinit, double& tend, double& tout, int& istate) with \n");
  DEBUG_PRINTF("tinit = %f, tend= %f, tout = %f, istate = %i\n", tinit, tend,  tout, istate);

  // For details on DLSODAR parameters, see opkdmain.f in externals/odepack
  doublereal tend_DR = tend  ;       // next point where output is desired (different from t!)
  doublereal tinit_DR = tinit;       // current (starting) time

  // === Pointers to function ===
  //  --> definition and initialisation thanks to wrapper:
  global_object = std11::static_pointer_cast<LsodarOSI>(shared_from_this()); // Warning: global object must be initialized to current one before pointers to function initialisation.

  // function to compute the righ-hand side of xdot = f(x,t) + Tu
  fpointer pointerToF = LsodarOSI_f_wrapper;

  // function to compute the Jacobian/x of the rhs.
  jacopointer pointerToJacobianF = LsodarOSI_jacobianf_wrapper; // function to compute the Jacobian/x of the rhs.

  // function to compute the constraints
  gpointer pointerToG;
  pointerToG = LsodarOSI_g_wrapper; // function to compute the constraints

  // === LSODAR CALL ===

  *_xtmp = *_xWork;
  if(istate == 3)
  {
    istate = 1; // restart TEMPORARY
  }

  _intData[4] = istate;

  // call LSODAR to integrate dynamical equation
  CNAME(dlsodar)(pointerToF,
                 &(_intData[0]),
                 _xtmp->getArray(),
                 &tinit_DR,
                 &tend_DR,
                 &(_intData[2]),
                 rtol.get(),
                 atol.get(),
                 &(_intData[3]),
                 &(_intData[4]),
                 &(_intData[5]),
                 rwork.get(),
                 &(_intData[6]),
                 iwork.get(),
                 &(_intData[7]),
                 pointerToJacobianF,
                 &(_intData[8]),
                 pointerToG, &
                 (_intData[1]),
                 jroot.get());

  // jroot: jroot[i] = 1 if g(i) has a root at t, else jroot[i] = 0.

  // === Post ===
  if(_intData[4] < 0)  // if istate < 0 => LSODAROSI failed
  {
    std::cout << "LSodar::integrate(...) failed - Istate = " << _intData[4] <<std::endl;
    std::cout << " -1 means excess work done on this call (perhaps wrong JT, or so small tolerance (ATOL and RTOL), or small maximum number of steps for one call (MXSTEP)). You should increase ATOL or RTOL or increase the MXSTEP" <<std::endl;
    std::cout << " -2 means excess accuracy requested (tolerances too small)." <<std::endl;
    std::cout << " -3 means illegal input detected (see printed message)." <<std::endl;
    std::cout << " -4 means repeated error test failures (check all inputs)." <<std::endl;
    std::cout << " -5 means repeated convergence failures (perhaps bad Jacobian supplied or wrong choice of JT or tolerances)." <<std::endl;
    std::cout << " -6 means error weight became zero during problem. (Solution component i vanished, and ATOL or ATOL(i) = 0.)" <<std::endl;
    std::cout << " -7 means work space insufficient to finish (see messages)." <<std::endl;
    RuntimeException::selfThrow("LsodarOSI, integration failed");
  }

  *_xWork = *_xtmp;
  istate = _intData[4];
  tout  = tinit_DR; // real ouput time
  tend  = tend_DR; // necessary for next start of DLSODAR


  if(istate == 3)
  {
    //      std:: std::cout << "ok\n";
    assert(true);
  }
  // Update counters
  count_NST = iwork[10];
  count_NFE = iwork[11];
  //  tinit = tinit_DR;
}


void LsodarOSI::updateState(const unsigned int level)
{
  // Compute all required (ie time-dependent) data for the DS of the OSI.
  DynamicalSystemsGraph::VIterator dsi, dsend;
  if(level == 1)  // ie impact case: compute velocity
  {
    for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
    {
      if(!checkOSI(dsi)) continue;
      SP::LagrangianDS lds = std11::static_pointer_cast<LagrangianDS>(_dynamicalSystemsGraph->bundle(*dsi));
      lds->computePostImpactVelocity();
    }
  }
  else if(level == 2)
  {
    double time = _simulation->nextTime();
    for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
    {
      if(!checkOSI(dsi)) continue;
      {
        SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
        ds->update(time);
      }
    }
  }
  else RuntimeException::selfThrow("LsodarOSI::updateState(index), index is out of range. Index = " + level);
}

struct LsodarOSI::_NSLEffectOnFreeOutput : public SiconosVisitor
{
  using SiconosVisitor::visit;

  OneStepNSProblem * _osnsp;
  SP::Interaction _inter;
  InteractionProperties& _interProp;
  _NSLEffectOnFreeOutput(OneStepNSProblem *p, SP::Interaction inter, InteractionProperties& interProp) :
    _osnsp(p), _inter(inter), _interProp(interProp) {};

  void visit(const NewtonImpactNSL& nslaw)
  {
    double e;
    e = nslaw.e();
    Index subCoord(4);
    subCoord[0] = 0;
    subCoord[1] = _inter->nonSmoothLaw()->size();
    subCoord[2] = 0;
    subCoord[3] = subCoord[1];
    SiconosVector & osnsp_rhs = *(*_interProp.workVectors)[LsodarOSI::OSNSP_RHS];
    subscal(e, *_inter->yOld(_osnsp->inputOutputLevel()), osnsp_rhs, subCoord, false); // q = q + e * q
  }

  // visit function added by Son (9/11/2010)
  void visit(const MultipleImpactNSL& nslaw)
  {
    ;
  }
  // note : no NewtonImpactFrictionNSL
};


void LsodarOSI::computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp)
{
  SP::OneStepNSProblems  allOSNS  = _simulation->oneStepNSProblems();
  SP::InteractionsGraph indexSet = osnsp->simulation()->indexSet(osnsp->indexSetLevel());
  SP::Interaction inter = indexSet->bundle(vertex_inter);

  VectorOfBlockVectors& DSlink = *indexSet->properties(vertex_inter).DSlink;
  // Get relation and non smooth law types
  RELATION::TYPES relationType = inter->relation()->getType();
  RELATION::SUBTYPES relationSubType = inter->relation()->getSubType();
  unsigned int sizeY = inter->nonSmoothLaw()->size();

  unsigned int relativePosition = 0;
  SP::Interaction mainInteraction = inter;
  Index coord(8);
  coord[0] = relativePosition;
  coord[1] = relativePosition + sizeY;
  coord[2] = 0;
  coord[4] = 0;
  coord[6] = 0;
  coord[7] = sizeY;
  SP::SiconosMatrix  C;
  //   SP::SiconosMatrix  D;
  //   SP::SiconosMatrix  F;
  SiconosVector& osnsp_rhs = *(*indexSet->properties(vertex_inter).workVectors)[LsodarOSI::OSNSP_RHS];

  SP::BlockVector Xfree;


  // All of these values should be stored in the node corrseponding to the Interactionwhen a MoreauJeanOSI scheme is used.

  /* V.A. 10/10/2010
   * Following the type of OSNS  we need to retrieve the velocity or the acceleration
   * This tricks is not very nice but for the moment the OSNS do not known if
   * it is in accelaration of not
   */

  //SP::OneStepNSProblems  allOSNS  = _simulation->oneStepNSProblems();
  if(((*allOSNS)[SICONOS_OSNSP_ED_SMOOTH_ACC]).get() == osnsp)
  {
    if(relationType == Lagrangian)
    {
      Xfree = DSlink[LagrangianR::xfree];
      DEBUG_EXPR(Xfree->display(););
    }
    // else if  (relationType == NewtonEuler)
    // {
    //   Xfree = inter->data(NewtonEulerR::free);
    // }
    assert(Xfree);
    //        std::cout << "Computeqblock Xfree (Gamma)========" << std::endl;
    //       Xfree->display();
  }
  else  if(((*allOSNS)[SICONOS_OSNSP_ED_IMPACT]).get() == osnsp)
  {
    Xfree = DSlink[LagrangianR::q1];
    //        std::cout << "Computeqblock Xfree (Velocity)========" << std::endl;
    //       Xfree->display();

  }
  else
    RuntimeException::selfThrow(" computeqBlock for Event Event-driven is wrong ");

  if(relationType == Lagrangian)
  {
    C = mainInteraction->relation()->C();
    if(C)
    {
      assert(Xfree);

      coord[3] = C->size(1);
      coord[5] = C->size(1);

      subprod(*C, *Xfree, osnsp_rhs, coord, true);
    }

    SP::SiconosMatrix ID(new SimpleMatrix(sizeY, sizeY));
    ID->eye();

    Index xcoord(8);
    xcoord[0] = 0;
    xcoord[1] = sizeY;
    xcoord[2] = 0;
    xcoord[3] = sizeY;
    xcoord[4] = 0;
    xcoord[5] = sizeY;
    xcoord[6] = 0;
    xcoord[7] = sizeY;
    // For the relation of type LagrangianRheonomousR
    if(relationSubType == RheonomousR)
    {
      if(((*allOSNS)[SICONOS_OSNSP_ED_SMOOTH_ACC]).get() == osnsp)
	    {
	      RuntimeException::selfThrow("LsodarOSI::computeFreeOutput not yet implemented for LCP at acceleration level with LagrangianRheonomousR");
	    }
      else if(((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp)
	    {
	      SiconosVector q = *DSlink[LagrangianR::q0];
	      SiconosVector z = *DSlink[LagrangianR::z];

	      std11::static_pointer_cast<LagrangianRheonomousR>(inter->relation())->computehDot(simulation()->getTkp1(), q, z);
	      *DSlink[LagrangianR::z] = z;
	      subprod(*ID, *(std11::static_pointer_cast<LagrangianRheonomousR>(inter->relation())->hDot()), osnsp_rhs, xcoord, false); // y += hDot
	    }
      else
        RuntimeException::selfThrow("LsodarOSI::computeFreeOutput not implemented for SICONOS_OSNSP ");
    }
    // For the relation of type LagrangianScleronomousR
    if(relationSubType == ScleronomousR)
    {
      if(((*allOSNS)[SICONOS_OSNSP_ED_SMOOTH_ACC]).get() == osnsp)
	    {
	      std11::static_pointer_cast<LagrangianScleronomousR>(inter->relation())->computedotjacqhXqdot(simulation()->getTkp1(), *inter, DSlink);
	      subprod(*ID, *(std11::static_pointer_cast<LagrangianScleronomousR>(inter->relation())->dotjacqhXqdot()), osnsp_rhs, xcoord, false); // y += NonLinearPart
	    }
    }
  }
  else
    RuntimeException::selfThrow("LsodarOSI::computeFreeOutput not yet implemented for Relation of type " + relationType);
  if(((*allOSNS)[SICONOS_OSNSP_ED_IMPACT]).get() == osnsp)
  {
    if(inter->relation()->getType() == Lagrangian || inter->relation()->getType() == NewtonEuler)
    {
      SP::SiconosVisitor nslEffectOnFreeOutput(new _NSLEffectOnFreeOutput(osnsp, inter, indexSet->properties(vertex_inter)));
      inter->nonSmoothLaw()->accept(*nslEffectOnFreeOutput);
    }
  }

}
void LsodarOSI::display()
{
  OneStepIntegrator::display();
  std::cout << " --- > LsodarOSI specific values: " <<std::endl;
  std::cout << "Number of equations: " << _intData[0] <<std::endl;
  std::cout << "Number of constraints: " << _intData[1] <<std::endl;
  std::cout << "itol, itask, istate, iopt, lrw, liw, jt: (for details on what are these variables see opkdmain.f)" <<std::endl;
  std::cout << _intData[2] << ", " << _intData[3] << ", " << _intData[4] << ", " << _intData[5] << ", " << _intData[6]  << ", " << _intData[7]  << ", " << _intData[8] <<std::endl;
  std::cout << "====================================" <<std::endl;
}
