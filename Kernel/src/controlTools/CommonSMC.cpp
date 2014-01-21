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

#include "ModelingTools.hpp"
#include "SimulationTools.hpp"
#include "CommonSMC.hpp"
#include "ControlSensor.hpp"
#include "Model.hpp"

void CommonSMC::initialize(const Model& m)
{
  if (_Csurface == NULL)
  {
    RuntimeException::selfThrow("CommonSMC::initialize - you have to set _Csurface before initializing the Actuator");
  }
  else
  {
    Actuator::initialize(m);
  }
  // We can only work with FirstOrderNonLinearDS, FirstOrderLinearDS and FirstOrderLinearTIDS
  // We can use the Visitor mighty power to check if we have the right type
  SP::DynamicalSystem DS = _sensor->getDS();
  Type::Siconos dsType;
  dsType = Type::value(*DS);
  // create the DS for the controller
  // if the DS we use is different from the DS we are controlling
  // when we want for instant to see how well the controller behaves
  // if the plant model is not exact, we can use the setSimulatedDS
  // method
  if (dsType == Type::FirstOrderLinearDS)
  {
    _DS_SMC.reset(new FirstOrderLinearDS(*(std11::static_pointer_cast<FirstOrderLinearDS>(DS))));
  }
  else if (dsType == Type::FirstOrderLinearTIDS)
  {
    _DS_SMC.reset(new FirstOrderLinearTIDS(*(std11::static_pointer_cast<FirstOrderLinearTIDS>(DS))));
  }
  else
  {
    RuntimeException::selfThrow("LinearSMC is not yet implemented for system of type" + dsType);
  }
  // We have to reset the _pluginb
  _DS_SMC->setComputebFunction(NULL);
  SP::SiconosVector dummyb(new SiconosVector(_DS_SMC->getN(), 0));
  _DS_SMC->setb(dummyb);
  // Get the dimension of the output
  // XXX What if there is more than one sensor ...
  double t0 = m.t0();
  double T = m.finalT() + _td->currentTimeStep(0);
  // create the SMC Model
  _SMC.reset(new Model(t0, T));
  // create the interaction
  _relationSMC.reset(new FirstOrderLinearTIR(_Csurface, _B));
  std11::static_pointer_cast<FirstOrderLinearTIR>(_relationSMC)->setDPtr(_D);
  unsigned int sDim = _Csurface->size(0);
  _nsLawSMC.reset(new RelayNSL(sDim, -_alpha, _alpha));

  _interactionSMC.reset(new Interaction(sDim, _nsLawSMC, _relationSMC));
  _SMC->nonSmoothDynamicalSystem()->insertDynamicalSystem(_DS_SMC);
  _SMC->nonSmoothDynamicalSystem()->link(_interactionSMC, _DS_SMC);
  _SMC->nonSmoothDynamicalSystem()->setControlProperty(_interactionSMC, true);
  // Set up the simulation
  _simulationSMC.reset(new TimeStepping(_td));
  _simulationSMC->setName("linear sliding mode controller simulation");
  _integratorSMC.reset(new ZeroOrderHoldOSI(_DS_SMC));
  _simulationSMC->insertIntegrator(_integratorSMC);
  // OneStepNsProblem
  _OSNSPB_SMC.reset(new Relay(_numericsSolverId));
  _OSNSPB_SMC->numericsSolverOptions()->dparam[0] = _precision;
  //    std::cout << _OSNSPB_SMC->numericsSolverOptions()->dparam[0] <<std::endl;
  _simulationSMC->insertNonSmoothProblem(_OSNSPB_SMC);
  // Finally we can initialize everything ...
  _SMC->initialize(_simulationSMC);

  // Handy
  _eventsManager = _simulationSMC->eventsManager();
  _lambda.reset(new SiconosVector(sDim));
  _lambda = _interactionSMC->lambda(0);

  SP::SimpleMatrix tmpM(new SimpleMatrix(_Csurface->size(0), _B->size(1)));
  _invCB.reset(new SimpleMatrix(*tmpM));
  prod(*_Csurface, *_B, *tmpM);
  invertMatrix(*tmpM, *_invCB);
  _us.reset(new SiconosVector(sDim));
  _ueq.reset(new SiconosVector(sDim));
}

void CommonSMC::setCsurface(const SiconosMatrix& newValue)
{
  // check dimensions ...
  if (_Csurface)
  {
    *_Csurface = newValue;
  }
  else
  {
    _Csurface.reset(new SimpleMatrix(newValue));
  }
}


void CommonSMC::computeUeq()
{
  unsigned int n = _DS_SMC->A()->size(1);
  // equivalent part, explicit contribution
  SP::SimpleMatrix tmpM1(new SimpleMatrix(_Csurface->size(0), n));
  SP::SimpleMatrix tmpN(new SimpleMatrix(n, n));
  SP::SimpleMatrix quasiProjB_A(new SimpleMatrix(_invCB->size(0), n));
  SP::SimpleMatrix tmpW(new SimpleMatrix(n, n, 0));
  SP::SiconosVector xTk(new SiconosVector(_sensor->y()));
  tmpW->eye();
  prod(*_Csurface, *_DS_SMC->A(), *tmpM1);
  // compute (CB)^{-1}CA
  prod(*_invCB, *tmpM1, *quasiProjB_A);
  prod(_thetaSMC-1, *quasiProjB_A, *xTk, *_ueq);

  // equivalent part, implicit contribution
  // XXX when to call this ?
  ZeroOrderHoldOSI& zoh = *std11::static_pointer_cast<ZeroOrderHoldOSI>(_integratorSMC);
  zoh.updateMatrices(_DS_SMC);

  // tmpN = B^{*}(CB)^{-1}CA
  prod(zoh.Bd(_DS_SMC), *quasiProjB_A, *tmpN, true);
  // W = I + \theta B^{*})CB)^{-1}CA
  scal(_thetaSMC, *tmpN, *tmpW, false);
  // compute e^{Ah}x_k
  prod(zoh.Ad(_DS_SMC), *xTk, *xTk);
  // xTk = (e^{Ah}-(1-\theta)\Psi_k\Pi_B A)x_k
  prod(_thetaSMC-1, *tmpN, _sensor->y(), *xTk, false);
  // compute the solution x_{k+1} of the system W*x_{k+1} = x_k
  tmpW->PLUForwardBackwardInPlace(*xTk);
  // add the contribution from the implicit part to ueq
  prod(-_thetaSMC, *quasiProjB_A, *xTk, *_ueq, false);

}

void CommonSMC::setCsurfacePtr(SP::SiconosMatrix newPtr)
{
  // check dimensions ...
  _Csurface = newPtr;
}

void CommonSMC::setSaturationMatrix(const SiconosMatrix& satM)
{
  // check dimensions ...
  if (satM.size(1) != _B->size(1))
  {
    RuntimeException::selfThrow("CommonSMC::setSaturationMatrix - inconstency between the dimension of the state space and D");
  }
  else
  {
    if (_D)
    {
      *_D = satM;
    }
    else
    {
      _D.reset(new SimpleMatrix(satM));
    }
  }
}

void CommonSMC::setSaturationMatrixPtr(SP::SiconosMatrix newPtr)
{
  // check dimensions ...
  if (newPtr->size(1) != _B->size(1))
  {
    RuntimeException::selfThrow("CommonSMC::setSaturationMatrixPtr - inconstency between the dimension of the state space and D");
  }
  else
  {
    _D = newPtr;
  }
}


void CommonSMC::setTimeDiscretisation(const TimeDiscretisation& td)
{
  _td.reset(new TimeDiscretisation(td));
};
