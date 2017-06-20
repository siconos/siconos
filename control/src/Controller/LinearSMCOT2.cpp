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

#include "FirstOrderLinearTIDS.hpp"
#include "EventDriven.hpp"

#include "LinearSMCOT2.hpp"

#include "ActuatorFactory.hpp"
#include "SiconosVector.hpp"
#include "LsodarOSI.hpp"
#include "ControlSensor.hpp"
#include "Model.hpp"
#include "TimeDiscretisation.hpp"
#include "NonSmoothDynamicalSystem.hpp"

LinearSMCOT2::LinearSMCOT2(SP::ControlSensor sensor): CommonSMC(LINEAR_SMC_OT2, sensor), _coeff(0.0)
{
}

LinearSMCOT2::~LinearSMCOT2()
{
}

void LinearSMCOT2::initialize(const Model& m)
{
  if (!_Csurface)
  {
    RuntimeException::selfThrow("CommonSMC::initialize - you have to set either _Csurface or h(.) before initializing the Actuator");
  }
  else
  {
    if (_Csurface && !_u)
      _u.reset(new SiconosVector(_Csurface->size(0), 0));
  }

  Actuator::initialize(m);

  // We can only work with FirstOrderNonLinearDS, FirstOrderLinearDS and FirstOrderLinearTIDS
  // We can use the Visitor mighty power to check if we have the right type
  SP::DynamicalSystem DS = _sensor->getDS();
  Type::Siconos dsType;
  dsType = Type::value(*DS);
  if (dsType == Type::FirstOrderLinearDS)
  {
    FirstOrderLinearDS& fods = static_cast<FirstOrderLinearDS&>(*DS);
    SP::SiconosVector x0(new SiconosVector(*fods.x0()));
    _DSPhi.reset(new FirstOrderLinearDS(x0));//(std11::static_pointer_cast<FirstOrderLinearDS>(DS))));
    _DSPred.reset(new FirstOrderLinearDS(x0));//(*(std11::static_pointer_cast<FirstOrderLinearDS>(DS))));
    if(fods.A())
    {
      SP::SiconosMatrix A(new SimpleMatrix(*fods.A()));
      _DSPhi->setA(*fods.A());
      _DSPred->setA(*fods.A());
    }
    if(fods.b())
    {
      SP::SiconosVector b(new SiconosVector(*fods.b()));
      _DSPhi->setb(*fods.b());
      _DSPred->setb(*fods.b());
    }
  }
  else if (dsType == Type::FirstOrderLinearTIDS)
  {
    _DSPhi.reset(new FirstOrderLinearTIDS(*(std11::static_pointer_cast<FirstOrderLinearTIDS>(DS))));
    _DSPred.reset(new FirstOrderLinearTIDS(*(std11::static_pointer_cast<FirstOrderLinearTIDS>(DS))));
  }
  else
  {
    RuntimeException::selfThrow("LinearSMCOT2 is not yet implemented for system of type" + dsType);
  }

  // We have to reset _pluginb
  // _DSPhi->setComputebFunction(NULL);
  // _DSPred->setComputebFunction(NULL);
  // XXX What if there is more than one sensor ...

  _indx = 0;
  //  _Phi.reset(new SimpleMatrix(_nDim, _nDim));
  //  _Phi->eye();
  //  _Xold.reset(new SiconosVector(_nDim));
  //  *_Xold = *(_sensor->y());
  double _t0 = m.t0();
  double _T = m.finalT() + _tdPhi->currentTimeStep(0);

  //  _XPhi.reset(new SiconosVector(_nDim));
  //  (*_XPhi) = _DS->getX0();
  //  _DSPhi->setXPtr(_XPhi);
  _XPhi = _DSPhi->x();

  //  _Xhat.reset(new SiconosVector(_nDim));
  //  *_Xhat = _DS->getX0();
  // _DSPred->setXPtr(_Xhat);
  _Xhat = _DSPred->x();
  SP::SiconosVector dummyb(new SiconosVector(_B->size(0), 0));
  _DSPred->setbPtr(dummyb);
  prod(*_B, *_u, *_DSPred->b());

  //  _Xhat.reset(new SiconosVector(_nDim, 0));
  //  _DSPred->setXPtr(_Xhat);

  _modelPhi.reset(new Model(_t0, _T));
  _PhiOSI.reset(new LsodarOSI());
  _modelPhi->nonSmoothDynamicalSystem()->insertDynamicalSystem(_DSPhi);
  _simulPhi.reset(new EventDriven(_tdPhi, 0));
  _simulPhi->setNonSmoothDynamicalSystemPtr(_modelPhi->nonSmoothDynamicalSystem());
  _simulPhi->prepareIntegratorForDS(_PhiOSI, _DSPhi, _modelPhi, _t0);
  _modelPhi->setSimulation(_simulPhi);
  _modelPhi->initialize();
  // Integration for Gamma
  _modelPred.reset(new Model(_t0, _T));
  _PredOSI.reset(new LsodarOSI());
  _modelPred->nonSmoothDynamicalSystem()->insertDynamicalSystem(_DSPred);
  _simulPred.reset(new EventDriven(_tdPred, 0));
  _simulPred->setNonSmoothDynamicalSystemPtr(_modelPred->nonSmoothDynamicalSystem());
  _simulPred->prepareIntegratorForDS(_PredOSI, _DSPred, _modelPred, _t0);
  _modelPred->setSimulation(_simulPred);
  _modelPred->initialize();

  _X = _sensor->yTk();

}

void LinearSMCOT2::actuate()
{
  double hCurrent = _tdPhi->currentTimeStep(_indx);
  // Get current value of the state
  // Update it
  *_XPhi = *_X;

  // We change the values of the state each time, so we need to change istate to 1
  // See LsodarOSI.cpp for the meaning of istate
  if (_indx > 0)
  {
    _simulPhi->setIstate(1);
    _simulPred->setIstate(1);
  }
  // Compute _XPhi = \Phi*X
  _simulPhi->advanceToEvent();
  _simulPhi->processEvents();
  // XXX small hack here
  SP::SiconosVector CS(new SiconosVector(_B->size(0)));
  _Csurface->getRow(0, *CS);
  _coeff = -1 / (CS->vector_sum() * hCurrent);
  double uEq = inner_prod(*CS, _coeff * (*_XPhi + *_X - *_Xhat));
  double uEqP;
  // We need to project
  // TODO this should work in more than 1D
  uEqP = std::min(uEq, 2.0);
  uEqP = std::max(uEqP, -2.0);
  (*_u)(_u->size() - 1) = uEqP;
  prod(*_B, *_u, *_DSPred->b());
  _indx++;
  *_Xhat = *_X;
  // Compute \hat{x}_k
  _simulPred->advanceToEvent();
  _simulPred->processEvents();
}

void LinearSMCOT2::setTimeDiscretisation(const TimeDiscretisation& td)
{ 
  _tdPhi.reset(new TimeDiscretisation(td));
  _tdPred.reset(new TimeDiscretisation(td));
}


AUTO_REGISTER_ACTUATOR(LINEAR_SMC_OT2, LinearSMCOT2)
