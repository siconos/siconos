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
void CommonSMC::initialize(SP::Model m)
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
  Type::Siconos dsType;
  dsType = Type::value(*_DS);
  // create the DS for the controller
  // if the DS we use is different from the DS we are controlling
  // when we want for instant to see how well the controller behaves
  // if the plant model is not exact, we can use the setSimulatedDS
  // method
  if (dsType == Type::FirstOrderLinearDS)
  {
    _DS_SMC.reset(new FirstOrderLinearDS(*(std11::static_pointer_cast<FirstOrderLinearDS>(_DS))));
  }
  else if (dsType == Type::FirstOrderLinearTIDS)
  {
    _DS_SMC.reset(new FirstOrderLinearTIDS(*(std11::static_pointer_cast<FirstOrderLinearTIDS>(_DS))));
  }
  else
  {
    RuntimeException::selfThrow("LinearSMC is not yet implemented for system of type" + dsType);
  }
  // We have to reset the _pluginb
  _DS_SMC->setComputebFunction(NULL);
  SP::SiconosVector dummyb(new SiconosVector(_nDim, 0));
  _DS_SMC->setb(dummyb);
  // Get the dimension of the output
  // XXX What if there is more than one sensor ...
  _sensor = std11::dynamic_pointer_cast<ControlSensor>(*(_allSensors->begin()));
  if (_sensor == NULL)
  {
    RuntimeException::selfThrow("LinearSMC::initialize - the given sensor is not a ControlSensor");
  }
  else
  {
    double t0 = _model->t0();
    double T = _model->finalT() + _timeDiscretisation->currentTimeStep();
    // create the SMC Model
    _SMC.reset(new Model(t0, T));
    // create the interaction
    _relationSMC.reset(new FirstOrderLinearTIR(_Csurface, _B));
    std11::static_pointer_cast<FirstOrderLinearTIR>(_relationSMC)->setDPtr(_D);
    _nsLawSMC.reset(new RelayNSL(_sDim));

    std::string id = "interaction for control";
    _interactionSMC.reset(new Interaction(_sDim, _nsLawSMC, _relationSMC));
    _SMC->nonSmoothDynamicalSystem()->insertDynamicalSystem(_DS_SMC);
    std::pair<DynamicalSystemsGraph::EDescriptor, InteractionsGraph::VDescriptor> descr =
      _SMC->nonSmoothDynamicalSystem()->link(_interactionSMC, _DS_SMC);
    _SMC->nonSmoothDynamicalSystem()->setControlProperty(descr.second, descr.first, true);
    // Copy the TD
    _tD_SMC.reset(new TimeDiscretisation(*_timeDiscretisation));
    // Set up the simulation
    _simulationSMC.reset(new TimeStepping(_tD_SMC));
    _simulationSMC->setName("linear sliding mode controller simulation");
//    _integratorSMC.reset(new Moreau(_DS_SMC, _thetaSMC));
    _integratorSMC.reset(new ZeroOrderHold(_DS_SMC));
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
    _lambda.reset(new SiconosVector(_sDim));
    _lambda = _interactionSMC->lambda(0);
    _xController = _DS_SMC->x();
    _u.reset(new SiconosVector(_nDim, 0));

    // XXX really stupid stuff
    _sampledControl.reset(new SiconosVector(_nDim, 0));
    _DS->setzPtr(_sampledControl);
  }
  _indx = 0;
  SP::SimpleMatrix tmpM(new SimpleMatrix(_Csurface->size(0), _B->size(1)));
  _invCB.reset(new SimpleMatrix(*tmpM));
  prod(*_Csurface, *_B, *tmpM);
  invertMatrix(*tmpM, *_invCB);
  _us.reset(new SiconosVector(_sDim));
  _ueq.reset(new SiconosVector(_sDim));
}

void CommonSMC::setCsurface(const SiconosMatrix& newValue)
{
  // check dimensions ...
  if (newValue.size(1) != _nDim)
  {
    RuntimeException::selfThrow("CommonSMC::setCsurface - inconstency between the dimension of the state space and Csurface");
  }
  else
  {
    if (_Csurface)
    {
      *_Csurface = newValue;
    }
    else
    {
      _Csurface.reset(new SimpleMatrix(newValue));
    }
    _sDim = newValue.size(0);
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
  SP::SiconosVector xTk(new SiconosVector(*(_sensor->y())));
  tmpW->eye();
  prod(*_Csurface, *_DS_SMC->A(), *tmpM1);
  prod(*_invCB, *tmpM1, *quasiProjB_A);
  prod(_thetaSMC-1, *quasiProjB_A, *xTk, *_ueq);

  // equivalent part, implicit contribution
  // XXX when to call this ?
  ZeroOrderHold& zoh = *std11::static_pointer_cast<ZeroOrderHold>(_integratorSMC);
  zoh.updateMatrices(*_DS_SMC);

  // tmpN = B^{*})CB)^{-1}CA
  axpy_prod(zoh.getPsi(*_DS_SMC), *quasiProjB_A, *tmpN, true);
  // W = I + \theta B^{*})CB)^{-1}CA
  scal(_thetaSMC, *tmpN, *tmpW, false);
  // compute e^{Ah}x_k
  prod(zoh.getPhi(*_DS_SMC), *xTk, *xTk);
  // xTk = (e^{Ah}-(1-\theta)\Psi_k\Pi_B A)x_k
  prod(_thetaSMC-1, *tmpN, *(_sensor->y()), *xTk, false);
  // compute the solution x_{k+1} of the system W*X = xTk
  tmpW->PLUForwardBackwardInPlace(*xTk);
  // add the contribution from the implicit part to ueq
  prod(-_thetaSMC, *quasiProjB_A, *xTk, *_ueq, false);

}

void CommonSMC::setCsurfacePtr(SP::SiconosMatrix newPtr)
{
  // check dimensions ...
  if (newPtr->size(1) != _nDim)
  {
    RuntimeException::selfThrow("CommonSMC::setCsurfacePtr - inconstency between the dimension of the state space and Csurface");
  }
  else
  {
    _Csurface = newPtr;
    _sDim = newPtr->size(0);
  }
}

void CommonSMC::setSaturationMatrix(const SiconosMatrix& newValue)
{
  // check dimensions ...
  if (newValue.size(1) != _nDim)
  {
    RuntimeException::selfThrow("CommonSMC::setSaturationMatrix - inconstency between the dimension of the state space and D");
  }
  else
  {
    if (_D)
    {
      *_D = newValue;
    }
    else
    {
      _D.reset(new SimpleMatrix(newValue));
    }
  }
}

void CommonSMC::setSaturationMatrixPtr(SP::SiconosMatrix newPtr)
{
  // check dimensions ...
  if (newPtr->size(1) != _sDim)
  {
    RuntimeException::selfThrow("CommonSMC::setSaturationMatrixPtr - inconstency between the dimension of the state space and D");
  }
  else
  {
    _D = newPtr;
  }
}

void CommonSMC::setB(const SiconosMatrix& B)
{
  _B.reset(new SimpleMatrix(B));
}
