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

/*! \file ExplicitTwisting.hpp
  \brief twisting algorithm with an explicit discretization
*/

#ifndef ExplicitTwisting_H
#define ExplicitTwisting_H

#include "CommonSMC.hpp"


/** Twisting Controller with an explicit discretization */
class ExplicitTwisting : public CommonSMC
{
private:
  /** serialization hooks */
  ACCEPT_SERIALIZATION(ExplicitTwisting);


protected:
  /** default constructor */
  ExplicitTwisting() {};

public:

  /** Constructor for a nonlinear system or the ActuatorFactory
   * \param sensor the ControlSensor feeding the Actuator
   */
  ExplicitTwisting(SP::ControlSensor sensor);

  /** Constructor for the linear case
   * \param sensor the ControlSensor feeding the Actuator
   * \param gain control magnitude
   * \param beta twisting parameter
   */
  ExplicitTwisting(SP::ControlSensor sensor, double gain, double beta);

  /** destructor
   */
  virtual ~ExplicitTwisting();

  /** Compute the new control law at each event
   * Here we are using the following formula:
   */
  virtual void actuate();

  virtual void initialize(const NonSmoothDynamicalSystem & nsds, const Simulation& s);
};
#endif
