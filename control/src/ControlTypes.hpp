/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

/*! \file ControlTypes.hpp
  \brief Enum related to the control toolbox.  */

#ifndef ControlType_H
#define ControlType_H

namespace siconos::control {

/** Actuator types */

enum class ActuatorType {

  PID = 100,
  LINEAR_SMC = 101,
  EXPLICIT_LINEAR_SMC = 103,
  LINEAR_SMC_OT2 = 104,
  LINEAR_SMC_IMPROVED = 105,
  TWISTING = 106,
  REGULAR_TWISTING = 107,
  EXPLICIT_TWISTING = 108,
};

/** Sensor types */
enum class SensorType {
  LINEAR_SENSOR = 100,
};

/** Observer types */
enum class ObserverType {
  LUENBERGER = 100,
  SLIDING_REDUCED_ORDER = 101,
};

}  // namespace siconos::control

#endif
