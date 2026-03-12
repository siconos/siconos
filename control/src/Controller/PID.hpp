/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include <boost/circular_buffer_fwd.hpp>
#include "Actuator.hpp"

namespace siconos::control {

class PID : public Actuator {
 private:

  ACCEPT_SERIALIZATION(PID);

  /** error vector */
  std::shared_ptr<boost::circular_buffer<double> > _err{nullptr};

  /** reference we are tracking */
  double _ref{0.};

  double _curDeltaT{0.};

  /** vector of gains */
  std::shared_ptr<siconos::algebra::SiconosVector> _K{nullptr};

 public:
  /** Constructor.
   *
   *  \param sensor the ControlSensor feeding the Actuator
   *  \param B the B matrix
   */
  PID(std::shared_ptr<ControlSensor> sensor,
      std::shared_ptr<siconos::algebra::SiconosMatrix> B = nullptr);

  /** destructor
   */
  virtual ~PID() noexcept = default;

  /** initialize actuator data.
   *
   *  \param nsds a siconos::modeling::NonSmoothDynamicalSystem
   *  \param s the simulation
   */
  virtual void initialize(const siconos::modeling::NonSmoothDynamicalSystem& nsds,
                          const siconos::simulation::Simulation& s);

  /**
     Compute the new control law at each event
     Here we are using the following formula:
     \f$ u_k = u_{k-1} + c_1 e_k + c_2 e_{k-1} + c_3 e_{k-2} \f$ , where
     \f[
     c_1 &= K_P - \frac{K_D}{\Delta t} + K_I \Delta t \\
     c_2 &= -1 - \frac{2K_D}{\Delta t} \\
     c_3 &= \frac{K_D}{\Delta t}
     \f]
   */
  void actuate();

  /** Set K
   *
   *  \param K std::shared_ptr<siconos::algebra::SiconosVector> \f$ [K_P, K_I, K_D] \f$
   */
  void setK(std::shared_ptr<siconos::algebra::SiconosVector> K);

  /** Set the value of _ref to reference
   *
   *  \param reference the new value
   */
  void inline setRef(double reference) { _ref = reference; }

  /** Get the timestep from the TimeDiscretisation associated with this PID controller
   *
   *  \param td the TimeDiscretisation for this Actuator
   */
  virtual void setTimeDiscretisation(const siconos::simulation::TimeDiscretisation& td);

  void setDeltaT(double deltaT) { _curDeltaT = deltaT; }

  /** display the data of the Actuator on the standard output
   */
  virtual void display() const;
};
// Register the observer into the factory
static ActuatorRegistration<PID> reg_APID(ActuatorType::PID);

}  // namespace siconos::control
#endif
