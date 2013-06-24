/* Siconos-Kernel, Copyright INRIA 2005-2012.
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


#include "FirstOrderLinearDS.hpp"
#include "TimeStepping.hpp"

#include "LinearSMCimproved.hpp"
#include "SiconosVector.hpp"
#include "ControlSensor.hpp"
#include "ZeroOrderHold.hpp"
#include "TimeDiscretisation.hpp"
#include "ActuatorFactory.hpp"

LinearSMCimproved::LinearSMCimproved(SP::TimeDiscretisation t):
  LinearSMC(t, LINEAR_SMC_IMPROVED), _predictionPerturbation(false), _inDisceteTimeSlidingPhase(false)
{
}

LinearSMCimproved::LinearSMCimproved(SP::TimeDiscretisation t, SP::SiconosMatrix B, SP::SiconosMatrix D):
  LinearSMC(t, B, D, LINEAR_SMC_IMPROVED), _predictionPerturbation(false), _inDisceteTimeSlidingPhase(false)
{
}

LinearSMCimproved::~LinearSMCimproved()
{
}

void LinearSMCimproved::predictionPerturbation()
{
  if (_us->normInf() < 1)
  {
    if (_inDisceteTimeSlidingPhase)
      *_ueq += *_us;
    else
      _inDisceteTimeSlidingPhase = true;
  }
}

void LinearSMCimproved::actuate()
{
  unsigned int sDim = _u->size();
  SP::SimpleMatrix tmpM1(new SimpleMatrix(*_Csurface));
  SP::SimpleMatrix tmpD(new SimpleMatrix(sDim, sDim, 0));
  SP::SiconosVector xTk(new SiconosVector(_sensor->y()));

  ZeroOrderHold& zoh = *std11::static_pointer_cast<ZeroOrderHold>(_integratorSMC);

  // equivalent part
  zoh.updateMatrices(_DS_SMC);
  prod(*_Csurface, zoh.Ad(_DS_SMC), *tmpM1);
  *tmpM1 *= -1.0;
  *tmpM1 += *_Csurface;
  prod(*_Csurface, zoh.Bd(_DS_SMC), *tmpD);
  // compute C(I-e^{Ah})x_k
  prod(*tmpM1, *xTk, *_ueq);
  // compute the solution u^eq of the system CB^{*}u^eq = C(I-e^{Ah})x_k
  tmpD->PLUForwardBackwardInPlace(*_ueq);

  *(_DS_SMC->x()) = *xTk;
  prod(*_B, *_ueq, *(_DS_SMC->b()));
  _simulationSMC->computeOneStep();
  _simulationSMC->nextStep();


  // discontinous part
  *_us = *_lambda;

  // prediction of the perturbation
  if (_predictionPerturbation)
    predictionPerturbation();

  // inject those in the system
  *_u = *_us;
  *_u += *_ueq;
  _indx++;

}

AUTO_REGISTER_ACTUATOR(LINEAR_SMC_IMPROVED, LinearSMCimproved)
