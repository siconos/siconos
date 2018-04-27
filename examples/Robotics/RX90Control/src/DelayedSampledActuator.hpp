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

/*! \file DelayedSampledActuator.h
  An example on how to implement a user-defined actuator.
*/

#ifndef DelayedSampledActuator_H
#define DelayedSampledActuator_H

#include "Actuator.hpp"

class SiconosMatrix;

/** \class DelayedSampledActuator
 *  \brief Specific Actuator used as an example on how to implement a user-defined actuator.
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.3.0.
 *  \date (Creation) february 09, 2007
 *
 *
 */
class DelayedSampledActuator : public Actuator
{
private:

  /** Default constructor
   */
  DelayedSampledActuator();

public:

  /** Constructor with a TimeDiscretisation.
   * \param a string, the type of the Actuator, which corresponds to the class type
   * \param a TimeDiscretisation*, (linked to a model).
   */
  DelayedSampledActuator(int, SP::TimeDiscretisation);

  /** Destructor
   */
  ~DelayedSampledActuator();

  /** initialize sensor data.
   */
  void initialize();

  /** capture data when the ActuatorEvent is processed ( for example set data[ActuatorEvent]=... )
   */
  void actuate();

  /** Encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Actuator*
   * \return a pointer on the derived type
   */
  static DelayedSampledActuator* convert(Actuator* s);
};

#endif
