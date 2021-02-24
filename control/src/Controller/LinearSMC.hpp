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

/*! \file LinearSMC.hpp
 *
 *  \brief Linear sliding mode controller
 */

#ifndef LinearSMC_H
#define LinearSMC_H

#include "CommonSMC.hpp"


/** Linear sliding mode controller
 *
 *  This controller implements the following sliding mode control strategy:
 *  \f$\begin{equation} u = u_{\mathrm{eq}} + u_s\qquad\text{where}\begin{cases} u_{\mathrm{eq}} \text{is the equivalent control input}\\u_s = -sgn(\sigma)\end{cases},\end{equation}\f$
 *
 * where \f$ \sigma = Cx\f$ is the user-defined sliding variable.
 *
 *  */
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
   * \param B the matrix 
   * \param D the D matrix in the FirstOrderLinearR
   * \param type do not set this yourself ! this is used in derived classes
   */
  LinearSMC(SP::ControlSensor sensor, SP::SimpleMatrix B,
            SP::SimpleMatrix D = std::shared_ptr<SimpleMatrix>(),
            unsigned int type = LINEAR_SMC);

  /** destructor
   */
  virtual ~LinearSMC();

  /** Compute the new control law at each event
   * Here we are using the following formula:
   * \f$u = u_{\mathrm{eq}} + u_s\f$
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
