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

#include "LinearSMC.hpp"
using namespace std;
using namespace ActuatorFactory;

LinearSMC::LinearSMC(SP::TimeDiscretisation t, SP::DynamicalSystem ds): CommonSMC(LINEAR_SMC, t, ds), _thetaSMC(0.5)
{
}

LinearSMC::LinearSMC(SP::TimeDiscretisation t, SP::DynamicalSystem ds, SP::SiconosMatrix B, SP::SiconosMatrix D):
  CommonSMC(LINEAR_SMC, t, ds, B, D), _thetaSMC(0.5)
{
}

LinearSMC::LinearSMC(SP::TimeDiscretisation t, SP::DynamicalSystem ds, const Sensors& sensorList): CommonSMC(LINEAR_SMC, t, ds, sensorList),
  _thetaSMC(0.5)
{
}

LinearSMC::~LinearSMC()
{
}

void LinearSMC::initialize(SP::Model m)
{
  CommonSMC::initialize(m);

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
    _DS_SMC.reset(new FirstOrderLinearDS(*(static_pointer_cast<FirstOrderLinearDS>(_DS))));
  }
  else if (dsType == Type::FirstOrderLinearTIDS)
  {
    _DS_SMC.reset(new FirstOrderLinearTIDS(*(static_pointer_cast<FirstOrderLinearTIDS>(_DS))));
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
  _sensor = dynamic_pointer_cast<ControlSensor>(*(_allSensors->begin()));
  if (_sensor == NULL)
  {
    RuntimeException::selfThrow("LinearSMC::initialize - the given sensor is not a ControlSensor");
  }
  else
  {
    double t0 = _model->t0();
    double T = _model->finalT();
    // create the SMC Model
    _SMC.reset(new Model(t0, T));
    // create the interaction
    _relationSMC.reset(new FirstOrderLinearR(_Csurface, _B));
    _relationSMC->setDPtr(_D);
    _nsLawSMC.reset(new RelayNSL(_sDim));

    _interactionSMC.reset(new Interaction(_sDim, _nsLawSMC, _relationSMC));
    //_interactionSMC.reset(new Interaction("SMC Interation", _DS_SMC, 0, _sDim, _nsLawSMC, _relationSMC));
    _SMC->nonSmoothDynamicalSystem()->insertDynamicalSystem(_DS_SMC);
    //    _SMC->nonSmoothDynamicalSystem()->insertInteraction(_interactionSMC);
    _SMC->nonSmoothDynamicalSystem()->link(_interactionSMC, _DS_SMC);
    //    SP::NonSmoothDynamicalSystem myNSDS(new NonSmoothDynamicalSystem(_DS_SMC, _interactionSMC));
    //    _SMC->setNonSmoothDynamicalSystemPtr(myNSDS);
    // Copy the TD
    _tD_SMC.reset(new TimeDiscretisation(*_timeDiscretisation));
    // Set up the simulation
    _simulationSMC.reset(new TimeStepping(_tD_SMC));
    _simulationSMC->setName("linear sliding mode controller simulation");
    _integratorSMC.reset(new Moreau(_DS_SMC, _thetaSMC));
    _simulationSMC->insertIntegrator(_integratorSMC);
    // OneStepNsProblem
    _OSNSPB_SMC.reset(new Relay(_numericsSolverId));
    _OSNSPB_SMC->numericsSolverOptions()->dparam[0] = _precision;
    //    cout << _OSNSPB_SMC->numericsSolverOptions()->dparam[0] << endl;
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
}

void LinearSMC::actuate()
{
  if (_indx > 0)
  {
    *(_DS_SMC->x()) = *(_sensor->y()); // XXX this is sooo wrong
    _simulationSMC->nextStep();
  }
  _simulationSMC->computeOneStep();
  prod(1.0, *_B, *_lambda, *_sampledControl, true);
  _indx++;

}

AUTO_REGISTER_ACTUATOR(LINEAR_SMC, LinearSMC);
