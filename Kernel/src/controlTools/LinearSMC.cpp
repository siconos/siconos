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

LinearSMC::LinearSMC(SP::TimeDiscretisation t, SP::DynamicalSystem ds): CommonSMC(101, t, ds)
{
}

LinearSMC::LinearSMC(SP::TimeDiscretisation t, SP::DynamicalSystem ds, const Sensors& sensorList): CommonSMC(101, t, ds, sensorList)
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
  if (dsType != Type::FirstOrderLinearDS && dsType != Type::FirstOrderLinearTIDS)
  {
    RuntimeException::selfThrow("LinearSMC is not yet implemented for system of type" + dsType);
  }

  // Get the dimension of the output
  // XXX What if there is more than one sensor ...

  _sensor = dynamic_pointer_cast<ControlSensor>(*(_allSensors->begin()));
  if (_sensor == NULL)
  {
    RuntimeException::selfThrow("LinearSMC::initialize - the given sensor is not a ControlSensor");
  }
  else
  {
    _t0 = _model->t0();
    _T = _model->finalT();
    // create the SMC Model
    _SMC.reset(new Model(_t0, _T));
    // create the DS for the controller
    SP::SimpleMatrix _A(new SimpleMatrix(_nDim, _nDim));
    _A->zero();
    _DS_SMC.reset(new FirstOrderLinearDS(static_pointer_cast<SiconosVector>(_DS->x0()), _A));
    //    _DS_SMC->setComputebFunction("RelayPlugin.so","computeB");
    // create the interaction
    _B.reset(new SimpleMatrix(_nDim, _sDim));
    _B->eye();
    (*_B) *= 2;
    SP::SimpleMatrix _D(new SimpleMatrix(_sDim, _sDim));
    _D->zero();

    _relationSMC.reset(new FirstOrderLinearR(_Csurface, _B));
    _relationSMC->setDPtr(_D);
    // XXX Where to put this ? What does this variable mean ?
    unsigned int nslawSize = 2;
    _nsLawSMC.reset(new RelayNSL(nslawSize));
    _interactionSMC.reset(new Interaction(_sDim, _nsLawSMC, _relationSMC));
    _SMC->nonSmoothDynamicalSystem()->insertDynamicalSystem(_DS_SMC);
    _SMC->nonSmoothDynamicalSystem()->link(_interactionSMC, _DS_SMC);
    double _hSMC = 0.01;
    _tD_SMC.reset(new TimeDiscretisation(_t0, _hSMC));
    // Set up the simulation
    _simulationSMC.reset(new TimeStepping(_tD_SMC));
    _simulationSMC->setName("linear sliding mode controller simulation");
    _thetaSMC = 0.5;
    _integratorSMC.reset(new Moreau(_DS_SMC, _thetaSMC));
    _simulationSMC->insertIntegrator(_integratorSMC);
    // OneStepNsProblem
    _LCP_SMC.reset(new LCP());
    _OSNSPB_SMC.reset(new Relay(SICONOS_RELAY_PGS));
    _simulationSMC->insertNonSmoothProblem(_OSNSPB_SMC);

    // Finally we can initialize everything ...
    _SMC->initialize(_simulationSMC);

    // Handy
    _eventsManager = _simulationSMC->eventsManager();
    _lambda.reset(new SimpleVector(nslawSize));
    _lambda = _interactionSMC->lambda(0);
    _xController = _DS_SMC->x();
    _u.reset(new SimpleVector(_nDim, 0));

    // XXX really stupid stuff
    _sampledControl.reset(new SimpleVector(2));
    _sampledControl->zero();
    _DS->setzPtr(_sampledControl);
    //    _DS_SMC->setzPtr(sampledControl);
  }
  _indx = 0;
  _initDone = true;
}

void LinearSMC::actuate()
{
  if (_indx > 0)
  {
    *(_DS_SMC->x()) = *(_sensor->y()); // XXX this is sooo wrong
    _simulationSMC->nextStep();
    //    cout << "current time: " << _tD_SMC->currentTime() << endl;
  }
  //    cout << "output from sensor: ";
  //    _sensor->y()->display();
  //    cout << endl;

  //  if ((*_sensor->y())(0) < 0 )
  //    cout << " Change of control";
  _simulationSMC->computeOneStep();
  //  _interactionSMC->display();
  // compute the new control and update it
  prod(1.0, *_B, *_lambda, *_sampledControl, true);  // XXX bad
  //  prod(1.0, *_B, *_lambda, *_u, true);
  //  _DS->setb(_u);
  //  _sampledControl->display();
  _indx++;

}



AUTO_REGISTER_ACTUATOR(101, LinearSMC);
