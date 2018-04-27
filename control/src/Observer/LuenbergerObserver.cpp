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

#include "ModelingTools.hpp"
#include "SimulationTools.hpp"
#include "LuenbergerObserver.hpp"
#include "ControlSensor.hpp"
#include "ObserverFactory.hpp"
#include "ControlZOHAdditionalTerms.hpp"

//#define DEBUG_BEGIN_END_ONLY
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"

void LuenbergerObserver::initialize(const NonSmoothDynamicalSystem& nsds, const Simulation &s)
{
  DEBUG_BEGIN("void LuenbergerObserver::initialize(const NonSmoothDynamicalSystem& nsds, const Simulation &s)\n");
  if (!_C)
  {
    RuntimeException::selfThrow("LuenbergerObserver::initialize - you have to set C before initializing the Observer");
  }
  else
  {
    Observer::initialize(nsds, s);
  }
  bool isDSinDSG0 = true;
  DynamicalSystemsGraph& originalDSG0 = *nsds.topology()->dSG(0);
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
      RuntimeException::selfThrow("LuenbergerObserver is not yet implemented for system of type" + dsType);

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
  _DS->resetToInitialState();
  DEBUG_EXPR(_DS->display(););
  _e.reset(new SiconosVector(_C->size(0)));
  _y.reset(new SiconosVector(_C->size(0)));

  double t0 = nsds.t0();
  double h = s.currentTimeStep();
  double T = nsds.finalT() + h;
  _nsds.reset(new NonSmoothDynamicalSystem(t0, T));
  _integrator.reset(new ZeroOrderHoldOSI());
  
  std11::static_pointer_cast<ZeroOrderHoldOSI>(_integrator)->setExtraAdditionalTerms(
    std11::shared_ptr<ControlZOHAdditionalTerms>(new ControlZOHAdditionalTerms()));
  _nsds->insertDynamicalSystem(_DS);

  // Add the necessary properties
  DynamicalSystemsGraph& DSG0 = *_nsds->topology()->dSG(0);
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
  _simulation.reset(new TimeStepping(_nsds, _td, 0));
  _simulation->associate(_integrator, _DS);

  // initialize error
  *_y = _sensor->y();
  DEBUG_END("void LuenbergerObserver::initialize(const NonSmoothDynamicalSystem& nsds, const Simulation &s)\n");
}

void LuenbergerObserver::process()
{
  if (!_pass)
    _pass = true;
  else
  {
    prod(*_C, *_xHat, *_e); // e = C*xhat_k

    *_e -= *_y; // e -= y_k


    // get measurement from sensor
    const SiconosVector& y = _sensor->y();

    // TODO theta method on the error
    _simulation->computeOneStep();
//    _simulation->nextStep();

    // update the current measured value
    *_y = y;
    *_xHat = _DS->getx();
  }
}

AUTO_REGISTER_OBSERVER(LUENBERGER, LuenbergerObserver);
