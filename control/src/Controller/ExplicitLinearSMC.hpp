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

/*! \file ExplicitLinearSMC.hpp
  \brief General interface to define an actuator
  */

#ifndef ExplicitLinearSMC_H
#define ExplicitLinearSMC_H

#include "CommonSMC.hpp"

class ExplicitLinearSMC : public CommonSMC
{
private:
  /** default constructor */
  ExplicitLinearSMC() {};

  /** serialization hooks */
  ACCEPT_SERIALIZATION(ExplicitLinearSMC);

  /** \f$\sigma = Cx\f$ */
  SP::SiconosVector _sigma;

public:

  /** Constructor.
   * \param sensor the ControlSensor feeding the Actuator
   */
  ExplicitLinearSMC(SP::ControlSensor sensor);

  /** Constructor.with all data
   * \param sensor the ControlSensor feeding the Actuator
   * \param B the B matrix
   */
  ExplicitLinearSMC(SP::ControlSensor sensor, SP::SimpleMatrix B);

  /** destructor
  */
  virtual ~ExplicitLinearSMC();

  /** Initializer
   * \param m the Model of the Simulation
   */
  virtual void initialize(const NonSmoothDynamicalSystem& nsds, const Simulation &s);

  /** Compute the new control law at each event
   * Here we are using the following formula:
   * TODO
   */
  void actuate();

};
#endif
