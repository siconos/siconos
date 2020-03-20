/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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


#include "FirstOrderLinearDS.hpp"
#include "TimeStepping.hpp"
#include "SiconosAlgebraProd.hpp"

#include "LinearSMCimproved.hpp"
#include "SiconosVector.hpp"
#include "ControlSensor.hpp"
#include "ZeroOrderHoldOSI.hpp"
#include "TimeDiscretisation.hpp"
#include "ActuatorFactory.hpp"

#include <boost/circular_buffer.hpp>

LinearSMCimproved::LinearSMCimproved(SP::ControlSensor sensor):
  LinearSMC(sensor, LINEAR_SMC_IMPROVED), _predictionPerturbation(false), _inDisceteTimeSlidingPhase(false), _ubPerturbation(0.0)
{
}

LinearSMCimproved::LinearSMCimproved(SP::ControlSensor sensor, SP::SimpleMatrix B, SP::SimpleMatrix D):
  LinearSMC(sensor, B, D, LINEAR_SMC_IMPROVED), _predictionPerturbation(false), _inDisceteTimeSlidingPhase(false),_ubPerturbation(0.0)
{
}

LinearSMCimproved::~LinearSMCimproved()
{
}

void LinearSMCimproved::initialize(const NonSmoothDynamicalSystem& nsds, const Simulation & s)
{
  LinearSMC::initialize(nsds,s);
  _up.reset(new SiconosVector(_us->size()));
  _measuredPert.reset(new boost::circular_buffer<SP::SiconosVector>(0));
  _predictedPert.reset(new boost::circular_buffer<SP::SiconosVector>(0));
}

void LinearSMCimproved::predictionPerturbation(const SiconosVector& xTk, SimpleMatrix& CBstar)
{
  if(_us->normInf() < _alpha)
  {
    if(_inDisceteTimeSlidingPhase)
    {
      SiconosVector& up = *_up;
      if(_measuredPert->full())
      {
        if(_measuredPert->size() > 1)
        {
          _measuredPert->rotate(_measuredPert->end()-1);
          _predictedPert->rotate(_predictedPert->end()-1);
        }
      }
      else
      {
        // inject new vector in the case where the measurement vector is not full.
        SP::SiconosVector sp1(new SiconosVector(_us->size(), 0));
        SP::SiconosVector sp2(new SiconosVector(_us->size(), 0));
        _measuredPert->push_front(sp1);
        _predictedPert->push_front(sp2);
      }

      // inject new measured value and also perturbation prediction
      SiconosVector& predictedPertC = *(*_predictedPert)[0];
      SiconosVector& measuredPertC = *(*_measuredPert)[0];

      // Cp_k = s_k + Cp_k-tilde
      prod(*_Csurface, xTk, measuredPertC);
      measuredPertC += *(*_predictedPert)[std::min((unsigned int)1, (unsigned int)_predictedPert->size()-1)];

      // compute prediction
      switch(_measuredPert->size()-1)
      {
      case 0:
        predictedPertC = measuredPertC;
        break;
      case 1:
        predictedPertC = 2*measuredPertC - *(*_measuredPert)[1];
        break;
      case 2:
        predictedPertC = 3*measuredPertC - 3*(*(*_measuredPert)[1]) + *(*_measuredPert)[2];
        break;
      default:
        RuntimeException::selfThrow("LinearSMCimproved::predictionPerturbation: unknown order " + std::to_string(_measuredPert->size()));
      }

      // Compute the control to counteract the perturbation
      up = predictedPertC;
      up *= -1;
      CBstar.PLUForwardBackwardInPlace(up);

      // project onto feasible set
      double norm = up.norm2();
      if(norm > _ubPerturbation)
      {
        up *= _ubPerturbation/norm;
        predictedPertC *= _ubPerturbation/norm;
      }
    }
    else
      _inDisceteTimeSlidingPhase = true;
  }
  else if(_inDisceteTimeSlidingPhase)
  {
    _inDisceteTimeSlidingPhase = false;
    _up->zero();
  }
}

void LinearSMCimproved::actuate()
{
  unsigned int sDim = _u->size();
  SP::SimpleMatrix tmpM1(new SimpleMatrix(*_Csurface));
  SP::SimpleMatrix CBstar(new SimpleMatrix(sDim, sDim, 0));
  SP::SiconosVector xTk(new SiconosVector(_sensor->y()));

  ZeroOrderHoldOSI& zoh = *std11::static_pointer_cast<ZeroOrderHoldOSI>(_integratorSMC);

  // equivalent part
  zoh.updateMatrices(_DS_SMC);
  prod(*_Csurface, zoh.Ad(_DS_SMC), *tmpM1);
  *tmpM1 *= -1.0;
  *tmpM1 += *_Csurface;
  prod(*_Csurface, zoh.Bd(_DS_SMC), *CBstar);
  // compute C(I-e^{Ah})x_k
  prod(*tmpM1, *xTk, *_ueq);
  // compute the solution u^eq of the system CB^{*}u^eq = C(I-e^{Ah})x_k
  CBstar->PLUForwardBackwardInPlace(*_ueq);

  *(_DS_SMC->x()) = *xTk;
  prod(*_B, *_ueq, *(std11::static_pointer_cast<FirstOrderLinearDS>(_DS_SMC)->b()));
  _simulationSMC->computeOneStep();
  _simulationSMC->nextStep();

  *_us = *_lambda;

  // inject those in the system
  *_u = *_us;
  *_u += *_ueq;

  // prediction of the perturbation
  if(_predictionPerturbation)
  {
    predictionPerturbation(*xTk, *CBstar);
    *_u += *_up;
  }

  _indx++;

}

void LinearSMCimproved::setPredictionOrder(unsigned int order)
{
  _measuredPert->set_capacity(order+1);
  _predictedPert->set_capacity(order+1);
}

AUTO_REGISTER_ACTUATOR(LINEAR_SMC_IMPROVED, LinearSMCimproved)
