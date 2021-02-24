/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

/*! \file RegularTwisting.hpp
  \brief twisting algorithm with a direct discretization
*/

#ifndef RegularTwisting_H
#define RegularTwisting_H

#include "Twisting.hpp"


/** Twisting Controller with a straightforward discretization */
class RegularTwisting : public Twisting
{
private:
  /** serialization hooks */
  ACCEPT_SERIALIZATION(RegularTwisting);


protected:
  /** default constructor */
  RegularTwisting() {};

public:

  /** Constructor for a nonlinear system or the ActuatorFactory
   * \param sensor the ControlSensor feeding the Actuator
   * \param hControl sampling period
   */
  RegularTwisting(SP::ControlSensor sensor);

  /** Constructor for the linear case
   * \param sensor the ControlSensor feeding the Actuator
   * \param gain control magnitude
   * \param beta twisting parameter
   * \param hControl sampling period
   */
  RegularTwisting(SP::ControlSensor sensor, double gain, double beta);

  /** destructor
   */
  virtual ~RegularTwisting();

};
#endif
