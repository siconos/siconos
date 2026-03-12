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

/*! \file ExplicitLinearSMC.hpp
  \brief General interface to define an actuator
  */

#ifndef ExplicitLinearSMC_H
#define ExplicitLinearSMC_H

#include "CommonSMC.hpp"

namespace siconos::control {

class ExplicitLinearSMC : public CommonSMC {
 private:
  ACCEPT_SERIALIZATION(ExplicitLinearSMC);

  /** \f$ \sigma = Cx \f$ */
  std::shared_ptr<siconos::algebra::SiconosVector> _sigma{nullptr};

 public:
  /** Constructor
   *
   *  \param sensor the ControlSensor feeding the Actuator
   */
  ExplicitLinearSMC(std::shared_ptr<ControlSensor> sensor);

  /** Constructor.with all data
   *
   *  \param sensor the ControlSensor feeding the Actuator
   *  \param B the B matrix
   */
  ExplicitLinearSMC(std::shared_ptr<ControlSensor> sensor,
                    std::shared_ptr<siconos::algebra::SiconosMatrix> B);

  /** destructor
   */
  virtual ~ExplicitLinearSMC() noexcept = default;

  /** Initializer
   *
   *  \param nsds current nonsmooth dynamical system
   *  \param s current simulation setup
   */
  virtual void initialize(const siconos::modeling::NonSmoothDynamicalSystem& nsds,
                          const siconos::simulation::Simulation& s);

  /** Compute the new control law at each event
   *  Here we are using the following formula:
   *  TODO
   */
  void actuate();
};

// Register the observer into the factory
static ActuatorRegistration<ExplicitLinearSMC> reg_AELSMC(ActuatorType::ExplicitLinearSMC);

}  // namespace siconos::control
#endif
