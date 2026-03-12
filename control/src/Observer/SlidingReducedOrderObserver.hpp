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

/*! \file SlidingReducedOrderObserver.hpp
  Discrete-time Sliding Observer
*/

#ifndef SlidingReducedOrderObserver_H
#define SlidingReducedOrderObserver_H

#include "Observer.hpp"
#include "SiconosMatrix.hpp"

namespace siconos::control {
class SlidingReducedOrderObserver : public Observer {
 private:
  ACCEPT_SERIALIZATION(SlidingReducedOrderObserver);

 protected:
  /** the vector defining the measurements (\f$ y = Cx \f$) */
  std::shared_ptr<siconos::algebra::SiconosMatrix> _C{nullptr};

  /** matrix multiplying the innovation term */
  std::shared_ptr<siconos::algebra::SiconosMatrix> _L{nullptr};

  double _theta{0.};

  // clumsy hack to do nothing the first time this Observer is called
  bool _pass{false};

 public:
  /** Constructor with the standard interface
   *
   *  \param sensor the std::shared_ptr<ControlSensor> that feed us with measurements
   *  \param xHat0 the initial guess for the state
   */
  SlidingReducedOrderObserver(std::shared_ptr<ControlSensor> sensor,
                              const siconos::algebra::SiconosVector& xHat0)
      : Observer(ObserverType::SlidingReducedOrder, sensor, xHat0) {}

  /** Constructor with all the data
   *
   *  \param sensor the sensor that feeds the Observer
   *  \param xHat0 the initial guess for the state
   *  \param C observation matrix
   *  \param L gain matrix
   */
  SlidingReducedOrderObserver(std::shared_ptr<ControlSensor> sensor,
                              const siconos::algebra::SiconosVector& xHat0,
                              std::shared_ptr<siconos::algebra::SiconosMatrix> C,
                              std::shared_ptr<siconos::algebra::SiconosMatrix> L)
      : Observer(ObserverType::SlidingReducedOrder, sensor, xHat0), _C(C), _L(L) {}

  ~SlidingReducedOrderObserver() noexcept = default;

  /** Update the estimate at each event
   */
  virtual void process();

  /** Initialization
   *
   *  \param nsds current nonsmooth dynamical system
   *  \param s current simulation setup
   */
  virtual void initialize(const siconos::modeling::NonSmoothDynamicalSystem& nsds,
                          const siconos::simulation::Simulation& s);

  /** Set the L matrix
   *
   *  \param L the new L matrix
   */
  inline void setLPtr(std::shared_ptr<siconos::algebra::SiconosMatrix> L) { _L = L; };

  /** Set the C matrix
   *
   *  \param C the new C matrix
   */
  inline void setCPtr(std::shared_ptr<siconos::algebra::SiconosMatrix> C) { _C = C; };
};

// Register the observer into the factory
static ObserverRegistration<SlidingReducedOrderObserver> reg_OSRO(
    ObserverType::SlidingReducedOrder);

}  // namespace siconos::control
#endif
