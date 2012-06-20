/* Siconos-Kernel, Copyright INRIA 2005-2011.
* Siconos is a program dedicated to modeling, simulation and control
* of non smooth dynamical systems.
* Siconos is a free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
* Siconos is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Siconos; if not, write to the Free Software
* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*
* Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/

#include "LinearSMCOT2.hpp"
using namespace std;
using namespace ActuatorFactory;

LinearSMCOT2::LinearSMCOT2(SP::TimeDiscretisation t, SP::DynamicalSystem ds): CommonSMC(LINEAR_SMC_OT2, t, ds)
{
}

LinearSMCOT2::LinearSMCOT2(SP::TimeDiscretisation t, SP::DynamicalSystem ds, const Sensors& sensorList): CommonSMC(LINEAR_SMC_OT2, t, ds, sensorList)
{
}

LinearSMCOT2::~LinearSMCOT2()
{
}

void LinearSMCOT2::initialize(SP::Model m)
{
  CommonSMC::initialize(m);

  // We can only work with FirstOrderNonLinearDS, FirstOrderLinearDS and FirstOrderLinearTIDS
  // We can use the Visitor mighty power to check if we have the right type
  Type::Siconos dsType;
  dsType = Type::value(*_DS);
  if (dsType == Type::FirstOrderLinearDS)
  {
    _DSPhi.reset(new FirstOrderLinearDS(*(static_pointer_cast<FirstOrderLinearDS>(_DS))));
    _DSPred.reset(new FirstOrderLinearDS(*(static_pointer_cast<FirstOrderLinearDS>(_DS))));
  }
  else if (dsType == Type::FirstOrderLinearTIDS)
  {
    _DSPhi.reset(new FirstOrderLinearTIDS(*(static_pointer_cast<FirstOrderLinearTIDS>(_DS))));
    _DSPred.reset(new FirstOrderLinearTIDS(*(static_pointer_cast<FirstOrderLinearTIDS>(_DS))));
  }
  else
  {
    RuntimeException::selfThrow("LinearSMCOT2 is not yet implemented for system of type" + dsType);
  }

  // We have to reset _pluginb
  _DSPhi->setComputebFunction(NULL);
  _DSPred->setComputebFunction(NULL);
  // XXX What if there is more than one sensor ...

  _sensor = dynamic_pointer_cast<ControlSensor>(*(_allSensors->begin()));
  if (_sensor == NULL)
  {
    RuntimeException::selfThrow("LinearSMCOT2::initialize - the given sensor is not a ControlSensor");
  }
  else
  {
    _u.reset(new SiconosVector(_nDim, 0));

    // XXX really stupid stuff
    _DS->setzPtr(_u);
  }
  _indx = 0;
  //  _Phi.reset(new SimpleMatrix(_nDim, _nDim));
  //  _Phi->eye();
  //  _Xold.reset(new SiconosVector(_nDim));
  //  *_Xold = *(_sensor->y());
  double _t0 = _model->t0();
  double _T = _model->finalT();

  _timeDPhi.reset(new TimeDiscretisation(*_timeDiscretisation));
  _timeDPred.reset(new TimeDiscretisation(*_timeDiscretisation));

  //  _XPhi.reset(new SiconosVector(_nDim));
  //  (*_XPhi) = _DS->getX0();
  //  _DSPhi->setXPtr(_XPhi);
  _XPhi = _DSPhi->x();

  //  _Xhat.reset(new SiconosVector(_nDim));
  //  *_Xhat = _DS->getX0();
  // _DSPred->setXPtr(_Xhat);
  _Xhat = _DSPred->x();
  _DSPred->setb(_u);

  //  _Xhat.reset(new SiconosVector(_nDim, 0));
  //  _DSPred->setXPtr(_Xhat);

  _modelPhi.reset(new Model(_t0, _T));
  _modelPhi->nonSmoothDynamicalSystem()->insertDynamicalSystem(_DSPhi);
  _PhiOSI.reset(new Lsodar(_DSPhi));
  _simulPhi.reset(new EventDriven(_timeDPhi, 0));
  _simulPhi->insertIntegrator(_PhiOSI);
  _modelPhi->initialize(_simulPhi);
  // Integration for Gamma
  _modelPred.reset(new Model(_t0, _T));
  _modelPred->nonSmoothDynamicalSystem()->insertDynamicalSystem(_DSPred);
  _PredOSI.reset(new Lsodar(_DSPred));
  _simulPred.reset(new EventDriven(_timeDPred, 0));
  _simulPred->insertIntegrator(_PredOSI);
  _modelPred->initialize(_simulPred);

  _X = _sensor->y();

}

void LinearSMCOT2::actuate()
{
  double hCurrent = _timeDiscretisation->currentTimeStep();
  // Get current value of the state
  // Update it
  *_XPhi = *_X;

  // We change the values of the state each time, so we need to change istate to 3
  // The first time, istate has to be 1 for initialization purposes
  // See Lsodar.cpp for the meaning of istate
  if (_indx > 0)
  {
    _simulPhi->setIstate(3);
    _simulPred->setIstate(3);
  }
  // Compute _XPhi = \Phi*X
  _simulPhi->processEvents();
  _simulPhi->advanceToEvent();
  // XXX small hack here
  SP::SiconosVector CS(new SiconosVector(_nDim));
  _Csurface->getRow(0, *CS);
  _coeff = -1 / (CS->sum() * hCurrent);
  double uEq = inner_prod(*CS, _coeff * (*_XPhi + *_X - *_Xhat));
  double uEqP;
  // We need to project
  // TODO this should work in more than 1D
  uEqP = min(uEq, 2.0);
  uEqP = max(uEqP, -2.0);
  (*_u)(_nDim - 1) = uEqP;
  _indx++;
  *_Xhat = *_X;
  // Compute \hat{x}_k
  _simulPred->processEvents();
  _simulPred->advanceToEvent();
}

AUTO_REGISTER_ACTUATOR(LINEAR_SMC_OT2, LinearSMCOT2);
