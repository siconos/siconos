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

/*! \file SensorX.h
  A specific sensor, to capture q vector ("position") in Lagrangian systems.
  Used as an example on how to implement new user-sensors.
*/

#ifndef SensorX_H
#define SensorX_H

#include "Sensor.hpp"

class SiconosVector;

/** \class SensorX
 *  \brief Specific Sensor to get position
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.3.0.
 *  \date (Creation) february 01, 2007
 *
 * A specific sensor, to capture q vector ("position") in Lagrangian systems.
 * Used as an example on how to implement new user-sensors.
 *
 */
class SensorX : public Sensor
{
private:

  /** Default constructor
   */
  SensorX();

  /** A copy vector of X at the last capture time
   */
  SP::SiconosVector storedX;

public:

  /** Constructor with a TimeDiscretisation.
   * \param a string, the type of the Sensor, which corresponds to the class type
   * \param a TimeDiscretisation*, (linked to a model).
   */
  SensorX(int, SP::TimeDiscretisation);

  /** Destructor
   */
  ~SensorX();

  /** initialize sensor data.
   */
  void initialize();

  /** capture data when the SensorEvent is processed ( for example set data[SensorEvent]=... )
   */
  void capture();

  /** Encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Sensor*
   * \return a pointer on the derived type
   */
  static SensorX* convert(Sensor* s);
};

#endif
