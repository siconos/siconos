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

/*! \file Twisting.hpp
  \brief twisting algorithm
*/

#ifndef Twisting_H
#define Twisting_H

#include "CommonSMC.hpp"


class Twisting : public CommonSMC
{
private:
  /** serialization hooks */
  ACCEPT_SERIALIZATION(Twisting);


protected:
  /** default constructor */
  Twisting() {};

public:
  /** Constructor for the ActuatorFactory
   * \param sensor the ControlSensor feeding the Actuator
   */
  Twisting(SP::ControlSensor sensor): CommonSMC(TWISTING, sensor) {};

  /** Constructor for a nonlinear system.
   * \param sensor the ControlSensor feeding the Actuator
   * \param hControl sampling period
   */
  Twisting(SP::ControlSensor sensor, double hControl);

  /** Constructor for the linear case
   * \param sensor the ControlSensor feeding the Actuator
   * \param gain control magnitude
   * \param beta twisting parameter
   * \param hControl sampling period
   */
  Twisting(SP::ControlSensor sensor, double gain, double beta, double hControl);

  /** destructor
   */
  virtual ~Twisting();

  /** Compute the new control law at each event
   * Here we are using the following formula:
   */
  virtual void actuate();

  /** set nonsmooth data: NormalConeNSL and AVI osnsp
   * \param hControl sampling period
   */
  void setNSdata(double hControl);

  virtual void initialize(const NonSmoothDynamicalSystem & nsds, const Simulation& s);
};
#endif
