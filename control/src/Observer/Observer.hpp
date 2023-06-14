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

/*! \file Observer.hpp
  \brief General interface to define an Observer
*/

#ifndef Observer_H
#define Observer_H

#include <cassert>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <string>

#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "SiconosSerialization.hpp"

namespace siconos::modeling {
class DynamicalSystem;
class NonSmoothDynamicalSystem;

};  // namespace siconos::modeling

namespace siconos::integrators {
class OneStepIntegrator;
}
namespace siconos::simulation {
class TimeDiscretisation;
class TimeStepping;
class Simulation;

};  // namespace siconos::simulation

namespace siconos::control {

class ControlSensor;

/** Observer types */
enum class ObserverType {
  Luenberger,
  SlidingReducedOrder,
};

/**
 Observers Base Class

 Abstract class, interface to user-defined observers.

 An Observer is dedicated to estimate the state of a DynamicalSystem given
 its dynamics, inputs and a initial estimate of the state.
*/

class Observer {
 protected:
  ACCEPT_SERIALIZATION(Observer);

  /** type of the Observer */
  ObserverType _type{ObserverType::Luenberger};

  /** the DynamicalSystem used in the Observer */
  std::shared_ptr<siconos::modeling::DynamicalSystem> _DS{nullptr};

  /** The TimeDiscretisation */
  std::shared_ptr<siconos::simulation::TimeDiscretisation> _td{nullptr};

  /** the sensor that feed the observer */
  std::shared_ptr<ControlSensor> _sensor{nullptr};

  /** estimated state */
  std::shared_ptr<siconos::algebra::SiconosVector> _xHat{nullptr};

  /** The error \f$e=\hat{y}-y\f$ */
  std::shared_ptr<siconos::algebra::SiconosVector> _e{nullptr};

  /** The measurements from the sensor */
  std::shared_ptr<siconos::algebra::SiconosVector> _y{nullptr};

  /** id of the Observer */
  std::string _id{"none"};

  // /** Model for integration */
  std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> _nsds{nullptr};

  /** Simulation for integration */
  std::shared_ptr<siconos::simulation::TimeStepping> _simulation{nullptr};

  /** Integration for integration */
  std::shared_ptr<siconos::integrators::OneStepIntegrator> _integrator{nullptr};

  // Rule of five
  Observer() = delete;
  Observer(const Observer&) = delete;
  Observer(Observer&&) = delete;
  Observer& operator=(const Observer&) = delete;
  Observer& operator=(Observer&&) = delete;

 public:
  /** Constructor with a TimeDiscretisation.
   *
   *  \param type the type of the Observer, which corresponds to the class type
   *  \param sensor the std::shared_ptr<Sensor> to get the measurements
   *  \param xHat0 the initial guess for the state
   *  \param newId the id of the Observer
   */
  Observer(ObserverType type, std::shared_ptr<ControlSensor> sensor,
           const siconos::algebra::SiconosVector& xHat0, const std::string& newId = "none");

  /** Constructor with a TimeDiscretisation.
   *
   *  \param type the type of the Observer, which corresponds to the class type.
   *  \param sensor the std::shared_ptr<Sensor> to get the measurements
   *  \param xHat0 the initial guess for the state
   *  \param ds the std::shared_ptr<siconos::modeling::DynamicalSystem> used as a model for the
   * real DynamicalSystem \param newId the id of the Observer
   */
  Observer(ObserverType type, std::shared_ptr<ControlSensor> sensor,
           const siconos::algebra::SiconosVector& xHat0,
           std::shared_ptr<siconos::modeling::DynamicalSystem> ds,
           const std::string& newId = "none");

  /** destructor
   */
  virtual ~Observer() noexcept = default;

  /** set id of the Observer
   *
   *  \param newId the new id.
   */
  inline void setId(const std::string& newId) { _id = newId; };

  /** get id of the Observer
   *
   *  \return a string
   */
  inline const std::string getId() const { return _id; };

  /** get the type of the Observer (ie class name)
   *
   *  \return an integer
   */
  inline ObserverType getType() const { return _type; };

  /** This is derived in child classes if they need to copy the TimeDiscretisation
   *  associated with this Sensor
   *
   *  \param td the TimeDiscretisation for this Sensor
   */
  virtual void setTimeDiscretisation(const siconos::simulation::TimeDiscretisation& td);

  /** initialize observer data.
   *
   *  \param nsds current nonsmooth dynamical system
   *  \param s current simulation setup
   */
  virtual void initialize(const siconos::modeling::NonSmoothDynamicalSystem& nsds,
                          const siconos::simulation::Simulation& s);

  /** capture data when the ObserverEvent is processed
   */
  virtual void process() = 0;

  /** display the data of the Observer on the standard output
   */
  void display() const;

  /** get the error e
   *
   *  \return a pointer to e
   */
  inline std::shared_ptr<siconos::algebra::SiconosVector> e() { return _e; }

  /** \return the estimated state
   */
  inline std::shared_ptr<siconos::algebra::SiconosVector> xHat() { return _xHat; }

  /** Set the DynamicalSystem used in the Observer
   *
   *  \param ds the DynamicalSystem used in the Observer
   */
  inline void setDS(std::shared_ptr<siconos::modeling::DynamicalSystem> ds) { _DS = ds; }

  /** \return the Model used in the Observer
   */
  virtual std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> getInternalNSDS() const
  {
    return _nsds;
  };
};

/** A class to handle observers creation

  Requirements:
  - the Observer type must be known and register

  For a XXXXObserver, add in the file describing the XXXXObserver class:

  static siconos::simulation::ObserverRegistration<siconos::simulation::XXXObserver>
 reg(siconos::simulation::ObserverType::XXXX);

  See ObserverType enum for the available names.

  Usage:

  auto observer = ObserverFactory::instance()->create(time, ObserverType::XXXX)

*/
class ObserverFactory {
  // Signature of Observer constructor
  using ObserverCreator = std::function<std::shared_ptr<Observer>(
      (std::shared_ptr<ControlSensor>, const siconos::algebra::SiconosVector&))>;

  /** map to connect observer type and the function used to create them */
  std::map<ObserverType, ObserverCreator> m_factories;

 public:
  /** Factory function which creates and returns an observer

   *  \param sensor the std::shared_ptr<Sensor> to get the measurements
   *  \param xHat0 the initial guess for the state
      \param type type of the observer (must be a ObserverType (enum))
      \return a pointer to observer
  */
  std::shared_ptr<Observer> create(std::shared_ptr<ControlSensor> sensor,
                                   const siconos::algebra::SiconosVector& xHat0,
                                   ObserverType type)
  {
    assert(m_factories.contains(type) && "unknown Observer type");
    return m_factories[type](sensor, xHat0);
  }

  /** access to the (singleton) factory instance */
  static ObserverFactory* instance()
  {
    static ObserverFactory factory;
    return &factory;
  }

  void registerCreator(ObserverType newtype, ObserverCreator caller)
  {
    m_factories[newtype] = caller;
  }
};

template <class T>
class ObserverRegistration {
 public:
  ObserverRegistration(ObserverType newtype)
  {
    ObserverFactory::instance()->registerCreator(
        newtype, [](std::shared_ptr<ControlSensor> sensor,
                    const siconos::algebra::SiconosVector& xHat0) {
          return std::make_shared<T>(sensor, xHat0);
        });
  }
};

}  // namespace siconos::control
#endif
