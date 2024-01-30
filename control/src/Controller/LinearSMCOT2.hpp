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

/*! \file LinearSMCOT2.hpp
  \brief General interface to define a sliding mode controller with
  disturbance compensation. Reference: Su, W.C.; Drakunov, S.; Ozguner, U.
  An O(T2) boundary layer in sliding mode for sampled-data systems
  */

#ifndef LinearSMCOT2_H
#define LinearSMCOT2_H

#include "CommonSMC.hpp"

namespace siconos::modeling {
class FirstOrderLinearDS;
}
namespace siconos::integrators {
class LsodarOSI;
}

namespace siconos::simulation {
class EventDriven;
}

namespace siconos::control {

class LinearSMCOT2 : public CommonSMC {
 private:
  /** Current value of the state (\f$ x_k\f$)*/
  std::shared_ptr<siconos::algebra::SiconosVector> _X{nullptr};
  /** Predicted current value of the state (\f$ \hat{x}_k = \Phi x_{k-1} + \Gamma u_{k-1}\f$)*/
  std::shared_ptr<siconos::algebra::SiconosVector> _Xhat{nullptr};
  /** Next value of the state only with the influence of the dynamic \f$ \XPhi = \Phi x_k\f$*/
  std::shared_ptr<siconos::algebra::SiconosVector> _XPhi{nullptr};
  /** Model for the computation of _XPhi*/
  std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> _nsdsPhi{nullptr};
  /** DynamicalSystem for the computation of _XPhi*/
  std::shared_ptr<siconos::modeling::FirstOrderLinearDS> _DSPhi{nullptr};
  /** TimeDiscretisation for the computation of _XPhi*/
  std::shared_ptr<siconos::simulation::TimeDiscretisation> _tdPhi{nullptr};
  /** OneSteoIntegrator for the computation of _XPhi*/
  std::shared_ptr<siconos::integrators::LsodarOSI> _PhiOSI{nullptr};
  /** Simulation for the computation of _XPhi*/
  std::shared_ptr<siconos::simulation::EventDriven> _simulPhi{nullptr};
  /** Model for the computation of Xhat*/
  std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> _nsdsPred{nullptr};
  /** TimeDiscretisation for the computation of Xhat*/
  std::shared_ptr<siconos::simulation::TimeDiscretisation> _tdPred{nullptr};
  /** OneSteoIntegrator for the computation of Xhat*/
  std::shared_ptr<siconos::integrators::LsodarOSI> _PredOSI{nullptr};
  /** Simulation for the computation of Xhat*/
  std::shared_ptr<siconos::simulation::EventDriven> _simulPred{nullptr};
  /** DynamicalSystem for the computation of _Xhat*/
  std::shared_ptr<siconos::modeling::FirstOrderLinearDS> _DSPred{nullptr};
  /** Coefficient*/
  double _coeff{0.};

  ACCEPT_SERIALIZATION(LinearSMCOT2);

 public:
  /** Constructor
   *
   *  \param sensor the ControlSensor feeding the Actuator
   */
  LinearSMCOT2(std::shared_ptr<ControlSensor> sensor);

  /** destructor
   */
  virtual ~LinearSMCOT2() noexcept = default;

  /** initialize actuator data.
   *
   *  \param nsds current nonsmooth dynamical system
   *  \param s current simulation setup
   */
  void initialize(const siconos::modeling::NonSmoothDynamicalSystem& nsds,
                  const siconos::simulation::Simulation& s);

  /** Compute the new control law at each event
   *  Here we are using the following formula:
   *  TODO
   */
  void actuate();

  /** This is derived in child classes if they need to copy the TimeDiscretisation
   *  associated with this Sensor
   *
   *  \param td the TimeDiscretisation for this Sensor
   */
  virtual void setTimeDiscretisation(const siconos::simulation::TimeDiscretisation& td);
};
// Register the observer into the factory
static ActuatorRegistration<LinearSMCOT2> reg_LSMCOT2(ActuatorType::LinearSMCOT2);

}  // namespace siconos::control
#endif
