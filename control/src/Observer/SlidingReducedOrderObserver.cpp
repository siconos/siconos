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

#include "ModelingTools.hpp"
#include "SimulationTools.hpp"
#include "Model.hpp"

#include "SlidingReducedOrderObserver.hpp"
#include "ControlSensor.hpp"
#include "ObserverFactory.hpp"
#include "ControlZOHAdditionalTerms.hpp"


void SlidingReducedOrderObserver::initialize(const Model& m)
{
  if (!_C)
  {
    RuntimeException::selfThrow("SlidingReducedOrderObserver::initialize - you have to set C before initializing the Observer");
  }
  else
  {
    Observer::initialize(m);
  }
  bool isDSinDSG0 = true;
  DynamicalSystemsGraph& originalDSG0 = *m.nonSmoothDynamicalSystem()->topology()->dSG(0);
  DynamicalSystemsGraph::VDescriptor originaldsgVD;
  if (!_DS) // No DynamicalSystem was given
  {
    // We can only work with FirstOrderNonLinearDS, FirstOrderLinearDS and FirstOrderLinearTIDS
    // We can use the Visitor mighty power to check if we have the right type
    DynamicalSystem& observedDS = *_sensor->getDS();
    Type::Siconos dsType;
    dsType = Type::value(observedDS);
    // create the DS for the controller
    // if the DS we use is different from the DS we are controlling
    // when we want for instant to see how well the controller behaves
    // if the plant model is not exact, we can use the setSimulatedDS
    // method
    if (dsType == Type::FirstOrderLinearDS)
    {
      _DS.reset(new FirstOrderLinearDS(static_cast<FirstOrderLinearDS&>(observedDS)));
    }
    else if (dsType == Type::FirstOrderLinearTIDS)
    {
      _DS.reset(new FirstOrderLinearTIDS(static_cast<FirstOrderLinearTIDS&>(observedDS)));
    }
    else
      RuntimeException::selfThrow("SlidingReducedOrderObserver is not yet implemented for system of type" + dsType);

    // is it controlled ?
    originaldsgVD = originalDSG0.descriptor(_sensor->getDS());
  }
  else
  {
    // is it controlled ?
    if (originalDSG0.is_vertex(_DS))
      originaldsgVD = originalDSG0.descriptor(_DS);
    else
      isDSinDSG0 = false;
  }

  // Initialize with the guessed state
  _DS->setX0Ptr(_xHat);

  _e.reset(new SiconosVector(_C->size(0)));
  _y.reset(new SiconosVector(_C->size(0)));

  double t0 = m.t0();
  double h = m.simulation()->currentTimeStep();
  double T = m.finalT() + h;
  _model.reset(new Model(t0, T));
  _integrator.reset(new ZeroOrderHoldOSI());
  std11::static_pointer_cast<ZeroOrderHoldOSI>(_integrator)->setExtraAdditionalTerms(
      std11::shared_ptr<ControlZOHAdditionalTerms>(new ControlZOHAdditionalTerms()));
  _model->nonSmoothDynamicalSystem()->insertDynamicalSystem(_DS);

  // Add the necessary properties
  DynamicalSystemsGraph& DSG0 = *_model->nonSmoothDynamicalSystem()->topology()->dSG(0);
  DynamicalSystemsGraph::VDescriptor dsgVD = DSG0.descriptor(_DS);
  // Observer part
  DSG0.L[dsgVD] = _L;
  DSG0.e[dsgVD] = _e;

  // Was the original DynamicalSystem controlled ?
  if (isDSinDSG0 && originalDSG0.B.hasKey(originaldsgVD))
  {
    DSG0.B[dsgVD] = originalDSG0.B[originaldsgVD];
    assert(originalDSG0.u[originaldsgVD] && "A DynamicalSystem is controlled but its control input has not been initialized yet");
    DSG0.u[dsgVD] = originalDSG0.u[originaldsgVD];
  }

  // all necessary things for simulation
  _simulation.reset(new TimeStepping(_td, 0));
  _simulation->prepareIntegratorForDS(_integrator, _DS, _model, t0);
  _model->setSimulation(_simulation);
  _model->initialize();

  // initialize error
  *_y = _sensor->y();
}

void SlidingReducedOrderObserver::process()
{
  if (!_pass)
  {
    _pass = true;
    //update the estimate using the first value of y, such that C\hat{x}_0 = y_0
    const SiconosVector& y = _sensor->y();
    _e->zero();
    prod(*_C, *_xHat, *_e);
    *_e -= y;

    SiconosVector tmpV(_DS->n());
    SimpleMatrix tmpC(*_C);
    for (unsigned int i = 0; i < _e->size(); ++i)
      tmpV(i) = (*_e)(i);

    tmpC.SolveByLeastSquares(tmpV);
    *(_xHat) -= tmpV;
    *(_DS->x()) -= tmpV;
    _DS->swapInMemory();
  }
  else
  {
    // get measurement from sensor
    const SiconosVector& y = _sensor->y();
    // update the current measured value
    *_y = y;

////    prod(*_C, _DS->getx(), *_e);
//    *_e -= y;
//
//    SiconosVector tmpV(_DS->n());
//    SimpleMatrix tmpC(*_C);
//    for (unsigned int i = 0; i < _e->size(); ++i)
//      tmpV(i) = (*_e)(i);
//
//    tmpC.SolveByLeastSquares(tmpV);
//    *(_DS->x()) -= tmpV;
    // First pass, set _e to 0, integrate the system
    // and get the innovation term
    _e->zero();
    _simulation->computeOneStep();

    // e = C*xhat_{k+1} - y_{k+1}
    prod(*_C, _DS->getx(), *_e);
    *_e -= *_y;

    // Second pass, now we update the state
    //
    // But first we need to reset the state to the
    // previous value (at t_k)
    _DS->setX(*_DS->xMemory()->getSiconosVector(0));
    // integrate with the new innovation term
    _simulation->computeOneStep();

    // We can go one step forward
    _simulation->nextStep();

    *_xHat = _DS->getx();
  }
}

AUTO_REGISTER_OBSERVER(SLIDING_REDUCED_ORDER, SlidingReducedOrderObserver);
