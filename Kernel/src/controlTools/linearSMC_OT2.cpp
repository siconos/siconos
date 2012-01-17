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

#include "linearSMC_OT2.hpp"
using namespace std;
using namespace ActuatorFactory;

linearSMC_OT2::linearSMC_OT2(int name, SP::TimeDiscretisation t, SP::Model m): commonSMC(name, t, m)
{
}

linearSMC_OT2::linearSMC_OT2(int name, SP::TimeDiscretisation t, SP::Model m, const Sensors& sensorList): commonSMC(name, t, m)
{
}

linearSMC_OT2::~linearSMC_OT2()
{
}

void linearSMC_OT2::initialize()
{
  commonSMC::initialize();

  // We can only work with FirstOrderNonLinearDS, FirstOrderLinearDS and FirstOrderLinearTIDS
  // XXX see how to check this
  _DS = dynamic_pointer_cast<FirstOrderLinearDS>(_model->nonSmoothDynamicalSystem()->dynamicalSystemNumber(0));
  if (_DS == NULL)
  {
    RuntimeException::selfThrow("The control of nonlinear System is not yet implemented");
  }

  // Get the dimension of the output
  // XXX What if there is more than one sensor ...

  _sensor = dynamic_pointer_cast<controlSensor>(*(_allSensors->begin()));
  if (_sensor == NULL)
  {
    RuntimeException::selfThrow("linearSMC_OT2::initialize - the given sensor is not a controlSensor");
  }
  else
  {
    _nDim = _sensor->getYDim();

    _u.reset(new SimpleVector(_nDim, 0));

    // XXX really stupid stuff
    _DS->setzPtr(_u);
  }
  _indx = 0;
  _initDone = true;
  //  _Phi.reset(new SimpleMatrix(_nDim, _nDim));
  //  _Phi->eye();
  //  _Xold.reset(new SimpleVector(_nDim));
  //  *_Xold = *(_sensor->y());
  double _t0 = _model->t0();
  double _T = _model->finalT();
  // XXX and everything breaks if this not a constant ...
  double _hSMC = _timeDiscretisation->currentTimeStep();
  _XPhi.reset(new SimpleVector(_nDim));
  (*_XPhi) = _DS->getX0();
  _DSPhi.reset(new FirstOrderLinearDS(_XPhi, _DS->A()));
  _DSPhi->setXPtr(_XPhi);
  _Xhat.reset(new SimpleVector(_nDim));
  *_Xhat = _DS->getX0();
  _DSPred.reset(new FirstOrderLinearDS(_Xhat, _DS->A()));
  _DSPred->setXPtr(_Xhat);
  _DSPred->setb(_u);

  _Xhat.reset(new SimpleVector(_nDim, 0));
  _DSPred->setXPtr(_Xhat);
  _modelPhi.reset(new Model(_t0, _T));
  _timeDPhi.reset(new TimeDiscretisation(_t0, _hSMC));
  _modelPhi->nonSmoothDynamicalSystem()->insertDynamicalSystem(_DSPhi);
  _PhiOSI.reset(new Lsodar(_DSPhi));
  _simulPhi.reset(new EventDriven(_timeDPhi, 0));
  _simulPhi->insertIntegrator(_PhiOSI);
  _modelPhi->initialize(_simulPhi);
  // Integration for Gamma
  _modelPred.reset(new Model(_t0, _T));
  _PredOSI.reset(new Lsodar(_DSPred));
  _timeDPred.reset(new TimeDiscretisation(_t0, _hSMC));
  _simulPred.reset(new EventDriven(_timeDPred, 0));
  _simulPred->insertIntegrator(_PredOSI);
  _modelPred->initialize(_simulPred);

  _X = _sensor->y();


}

void linearSMC_OT2::actuate()
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
  _coeff = -1 / (_Csurface->sum() * hCurrent);
  double uEq = inner_prod(*_Csurface, _coeff * (*_XPhi + *_X - *_Xhat));
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

void linearSMC_OT2::setCsurface(const SimpleVector& newValue)
{
  // check dimensions ...
  if (newValue.size() != _nDim)
  {
    RuntimeException::selfThrow("linearSMC_OT2::setCsurface - inconstency between the size of the dimension of the state space and Csurface");
  }
  else
  {
    if (_Csurface)
    {
      *_Csurface = newValue;
    }
    else
    {
      _Csurface.reset(new SimpleVector(newValue));
    }
  }
}

void linearSMC_OT2::setCsurfacePtr(SP::SimpleVector newPtr)
{
  // check dimensions ...
  if (newPtr->size() != _nDim)
  {
    RuntimeException::selfThrow("linearSMC_OT2::setCsurfacePtr - inconstency between the size of the dimension of the state space and Csurface");
  }
  else
  {
    _Csurface = newPtr;
  }
}

AUTO_REGISTER_ACTUATOR(104, linearSMC_OT2);
