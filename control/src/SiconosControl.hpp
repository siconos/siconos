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
/*! \file SiconosKernel.hpp
\brief Include files related to controlTools.
*/

#include "ControlTypeDef.hpp"

// Sensors - generic
#include "SensorFactory.hpp"
#include "SensorEvent.hpp"
#include "ControlSensor.hpp"

// Sensors - available
#include "LinearSensor.hpp"

// Actuator - generic
#include "ActuatorFactory.hpp"
#include "ActuatorEvent.hpp"
#include "CommonSMC.hpp"
// Actuator - available
#include "PID.hpp"
#include "LinearSMC.hpp"
#include "ExplicitLinearSMC.hpp"
#include "LinearSMCOT2.hpp"
#include "LinearSMCimproved.hpp"
#include "Twisting.hpp"

// Observer - generic
#include "Observer.hpp"

// Observer - available
#include "LuenbergerObserver.hpp"
#include "SlidingReducedOrderObserver.hpp"

// Misc
#include "ControlManager.hpp"
#include "MatrixIntegrator.hpp"

// sugar
#include "ControlZOHSimulation.hpp"
#include "ControlLsodarSimulation.hpp"

