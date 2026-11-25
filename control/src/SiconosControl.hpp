/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
/*! \file SiconosKernel.hpp
\brief Include files related to controlTools.
*/

// Sensors - generic
#include "ControlSensor.hpp"
#include "SensorEvent.hpp"

// Sensors - available
#include "LinearSensor.hpp"

// Actuator - generic
#include "ActuatorEvent.hpp"
#include "CommonSMC.hpp"
// Actuator - available
#include "ExplicitLinearSMC.hpp"
#include "LinearSMC.hpp"
#include "LinearSMCOT2.hpp"
#include "LinearSMCimproved.hpp"
#include "PID.hpp"
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
#include "ControlLsodarSimulation.hpp"
#include "ControlZOHSimulation.hpp"
