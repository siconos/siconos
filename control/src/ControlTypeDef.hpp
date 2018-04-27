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

/*! \file ControlTypeDef.hpp
  \brief Typedef for control-related objects
  */

#ifndef ControlTypeDef_H
#define ControlTypeDef_H

/** Actuator types */
#define PID_                       100
#define LINEAR_SMC                 101
#define EXPLICIT_LINEAR_SMC        103
#define LINEAR_SMC_OT2             104
#define LINEAR_SMC_IMPROVED        105
#define TWISTING                   106


/** Sensor types */
#define LINEAR_SENSOR              100

/** Observer types */
#define LUENBERGER                 100
#define SLIDING_REDUCED_ORDER      101

/** Event types
  \warning You have also to update the
  description in Event.hpp
*/
#define SENSOR_EVENT               3
#define OBSERVER_EVENT             4
#define ACTUATOR_EVENT             5

/** Base type forward declaration */

#endif
