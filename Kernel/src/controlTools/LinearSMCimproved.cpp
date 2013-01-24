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


#include "LinearSMCimproved.hpp"

using namespace ActuatorFactory;

LinearSMCimproved::LinearSMCimproved(SP::TimeDiscretisation t, SP::DynamicalSystem ds):
  LinearSMC(t, ds, LINEAR_SMC_IMPROVED)
{
}

LinearSMCimproved::LinearSMCimproved(SP::TimeDiscretisation t, SP::DynamicalSystem ds, SP::SiconosMatrix B, SP::SiconosMatrix D):
  LinearSMC(t, ds, B, D, LINEAR_SMC_IMPROVED)
{
}

LinearSMCimproved::LinearSMCimproved(SP::TimeDiscretisation t, SP::DynamicalSystem ds, const Sensors& sensorList):
  LinearSMC(t, ds, sensorList, LINEAR_SMC_IMPROVED)
{
}

LinearSMCimproved::~LinearSMCimproved()
{
}

void LinearSMCimproved::actuate()
{
  unsigned int n = _DS_SMC->A()->size(1);
  // equivalent part, explicit contribution
  SP::SimpleMatrix tmpM1(new SimpleMatrix(*_Csurface));
  SP::SimpleMatrix tmpD(new SimpleMatrix(_sDim, _sDim, 0));
  SP::SiconosVector xTk(new SiconosVector(*(_sensor->y())));

  ZeroOrderHold& zoh = *std11::static_pointer_cast<ZeroOrderHold>(_integratorSMC);

  // equivalent part
  zoh.updateMatrices(*_DS_SMC);
  prod(*_Csurface, zoh.getPhi(*_DS_SMC), *tmpM1);
  *tmpM1 *= -1.0;
  *tmpM1 += *_Csurface;
  prod(*_Csurface, zoh.getPsi(*_DS_SMC), *tmpD);
  // compute C(I-e^{Ah})x_k
  prod(*tmpM1, *xTk, *_ueq);
  // compute the solution x_{k+1} of the system W*X = e^{Ah}x_k
  tmpD->PLUForwardBackwardInPlace(*_ueq);

  if (_indx > 0)
  {
    *(_DS_SMC->x()) = *xTk; // XXX this is sooo wrong
    prod(*_B, *_ueq, *(_DS_SMC->b()));
    _simulationSMC->nextStep();
  }
  _simulationSMC->computeOneStep();


  // discontinous part
  double h = _timeDiscretisation->currentTimeStep();
  prod(h, prod(*_Csurface, *_B), *_lambda, *_us);
  tmpD->PLUForwardBackwardInPlace(*_us);

  // inject those in the system
  prod(1.0, *_B, *_us, *_sampledControl);
  prod(1.0, *_B, *_ueq, *_sampledControl, false);
  _indx++;

}

AUTO_REGISTER_ACTUATOR(LINEAR_SMC_IMPROVED, LinearSMCimproved)
