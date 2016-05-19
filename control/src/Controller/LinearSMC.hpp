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

/*! \file LinearSMC.hpp
  \brief General interface to define an actuator
*/

#ifndef LinearSMC_H
#define LinearSMC_H

#include "CommonSMC.hpp"


class LinearSMC : public CommonSMC
{
private:
  /** serialization hooks */
  ACCEPT_SERIALIZATION(LinearSMC);


protected:
  /** default constructor */
  LinearSMC() {};


public:

  /** Constructor for the ActuatorFactory
   * \param sensor the ControlSensor feeding the Actuator
   * \param type do not set this yourself ! this is used in derived classes
   */
  LinearSMC(SP::ControlSensor sensor, unsigned int type = LINEAR_SMC);

  /** Constructor with all the data
   * \param sensor the ControlSensor feeding the Actuator
   * \param B the B matrix in the FirstOrderLinearR
   * \param D the D matrix in the FirstOrderLinearR
   * \param type do not set this yourself ! this is used in derived classes
   */
  LinearSMC(SP::ControlSensor sensor, SP::SimpleMatrix B, SP::SimpleMatrix D = std11::shared_ptr<SimpleMatrix>(), unsigned int type = LINEAR_SMC);

  /** destructor
   */
  virtual ~LinearSMC();

  /** Compute the new control law at each event
   * Here we are using the following formula:
   * TODO
   */
  virtual void actuate();

  /** Set the D matrix
  * \param D the new D matrix
  */
  inline void setDPtr(SP::SimpleMatrix D)
  {
    _D = D;
  };


};
#endif
