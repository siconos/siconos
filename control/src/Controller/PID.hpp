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

/*! \file PID.hpp
  \brief Proportional-Integral-Derivative Controller
*/

#ifndef PID_H
#define PID_H

#include "Actuator.hpp"
#include "SiconosAlgebraTypeDef.hpp"
#include <boost/circular_buffer.hpp>

class PID : public Actuator
{
private:
  /** default constructor */
  PID() {};

  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(PID);

  /** error vector */
  std11::shared_ptr<boost::circular_buffer<double> > _err;

  /** reference we are tracking */
  double _ref;

  double _curDeltaT;

  /** vector of gains */
  SP::SiconosVector _K;

public:

  /** Constructor.
   * \param sensor the ControlSensor feeding the Actuator
   * \param B the B matrix
   */
  PID(SP::ControlSensor sensor, SP::SimpleMatrix B = std11::shared_ptr<SimpleMatrix>());

  /** destructor
   */
  virtual ~PID();

  /** initialize actuator data.
   * \param m a SP::Model
   */
  virtual void initialize(const Model& m);

  /** Compute the new control law at each event
   * Here we are using the following formula:
   * \f$ u_k = u_{k-1} + c_1 e_k + c_2 e_{k-1} + c_3 e_{k-2} \f$ , where
   * \f{array} c_1 &= K_P - \frac{K_D}{\Delta t} + K_I \Delta t \\
   * c_2 &= -1 - \frac{2K_D}{\Delta t} \\
   * c_3 &= \frac{K_D}{\Delta t} \\
   * \f}
   */
  void actuate();

  /** Set K
   * \param K SP::SiconosVector \f$ [K_P, K_I, K_D] \f$
   */
  void setK(SP::SiconosVector K);

  /** Set the value of _ref to reference
   * \param reference the new value
   */
  void inline setRef(double reference)
  {
    _ref = reference;
  }

  /** Get the timestep from the TimeDiscretisation associated with this PID controller
  *  \param td the TimeDiscretisation for this Actuator
  */
  virtual void setTimeDiscretisation(const TimeDiscretisation& td);

  void setDeltaT(double deltaT)
  {
    _curDeltaT = deltaT;
  }
/** display the data of the Actuator on the standard output
   */
  virtual void display() const;

};
#endif
