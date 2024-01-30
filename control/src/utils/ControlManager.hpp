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

/*! \file ControlManager.hpp
  \brief Tools to provide control in a Simulation: Sensors, Observer and Actuators.
*/

#ifndef ControlManager_H
#define ControlManager_H

#include "SiconosVector.hpp"
#include <memory>
#include <set>

#include "SiconosSerialization.hpp"


namespace siconos::modeling {
class NonSmoothDynamicalSystem;
class DynamicalSystem;
}  // namespace siconos::modeling

namespace siconos::simulation {
class TimeDiscretisation;
class Simulation;
}  // namespace siconos::simulation

namespace siconos::control {

enum class ObserverType;
enum class SensorType;
enum class ActuatorType;
class Observer;
class ControlSensor;
class Actuator;
class Sensor;

/**
   ControlManager Class: tools to provide control in a Simulation (Sensors, Actuators,
   Observers)

   This class is used to handle all the sensors and actuators declared by the user and to
   schedule them into the simulation.

   A ControlManager has:
   - a list of Sensor
   - a list of Actuator
   - a link to an existing Simulation

   The usual way to define control over a system is as follows:
   - declare a ControlManager and associate it with a Simulation
   - add some sensors and actuators into the ControlManager
   - initialize the ControlManager (which will result in the recording of all actuators and
   sensors into the list of events processed during the simulation)
   - optionally add some new sensor/actuator at any time but with a specific function:
   addAndRecord(...). A call to this function results in the creation of a Sensor/Actuator and
   in the insertion of the corresponding event into the simulation eventsManager.
*/

class ControlManager {
 private:
  ACCEPT_SERIALIZATION(ControlManager);

 protected:
  /** A list of Sensors */
  std::set<std::shared_ptr<Sensor>> _allSensors;

  /** A list of Actuators */
  std::set<std::shared_ptr<Actuator>> _allActuators;

  /** A list of Observers */
  std::set<std::shared_ptr<Observer>> _allObservers;

  /** The simulation linked to this ControlManager */
  std::shared_ptr<siconos::simulation::Simulation> _sim{nullptr};

  // /** default constructor
  //  */
  // ControlManager(): _sim(std::shared_ptr<siconos::simulation::Simulation>()) {};

  // Rule of five
  ControlManager() = delete;
  ControlManager(const ControlManager&) = delete;
  ControlManager(ControlManager&&) = delete;
  ControlManager& operator=(const ControlManager&) = delete;
  ControlManager& operator=(ControlManager&&) = delete;

  /** Create associated Event and give the opportunity to get the TimeDiscretisation
   *
   *  \param s a Sensor
   *  \param td a TimeDiscretisation asociated with this Sensor
   */
  void linkSensorSimulation(std::shared_ptr<Sensor> s,
                            std::shared_ptr<siconos::simulation::TimeDiscretisation> td);

  /** Create associated Event and give the opportunity to get the TimeDiscretisation
   *
   *  \param act a Sensor
   *  \param td a TimeDiscretisation asociated with this Sensor
   */
  void linkActuatorSimulation(std::shared_ptr<Actuator> act,
                              std::shared_ptr<siconos::simulation::TimeDiscretisation> td);

  /** Create associated Event and give the opportunity to get the TimeDiscretisation
   *
   *  \param obs a Sensor
   *  \param td a TimeDiscretisation asociated with this Sensor
   */
  void linkObserverSimulation(std::shared_ptr<Observer> obs,
                              std::shared_ptr<siconos::simulation::TimeDiscretisation> td);

 public:
  /** Constructor with a Simulation, to which control will be applied.
   *
   *  \param sim the Simulation
   */
  ControlManager(std::shared_ptr<siconos::simulation::Simulation> sim);

  /** destructor
   */
  virtual ~ControlManager() noexcept = default;

  /** get the Simulation linked to this ControlManager
   *
   *  \return a std::shared_ptr<siconos::simulation::Simulation>
   */
  inline std::shared_ptr<siconos::simulation::Simulation> simulation() const { return _sim; };

  /**\return the list of Sensors associated to this manager.
   *
   */
  inline const std::set<std::shared_ptr<Sensor>> getSensors() const { return _allSensors; };

  /** \return the list of Actuators associated to this manager.
   *
   */
  inline const std::set<std::shared_ptr<Actuator>> getActuators() const
  {
    return _allActuators;
  };

  /** \return the list of Observers associated to this manager */
  inline const std::set<std::shared_ptr<Observer>> getObservers() const
  {
    return _allObservers;
  };

  /** To build and add a new Sensor in the Manager
   *
   *  \param name the type of the Sensor
   *  \param td the std::shared_ptr<siconos::simulation::TimeDiscretisation> of the Sensor
   *  \param ds the DynamicalSystem used in the Sensor
   *  \return a std::shared_ptr<Sensor> to the added Sensor
   */
  std::shared_ptr<Sensor> addSensor(
      SensorType type, std::shared_ptr<siconos::simulation::TimeDiscretisation> td,
      std::shared_ptr<siconos::modeling::DynamicalSystem> ds);

  /** To build, add, initialize a new Sensor in the Manager and record
   *  it in the simulation This function is only useful to add a new
   *  Sensor after the initialization of the manager else call
   *  addSensor()
   *
   *  \param name the type (int) of the Sensor
   *  \param td the std::shared_ptr<siconos::simulation::TimeDiscretisation> of the Sensor
   *  \param ds the DynamicalSystem used in the Sensor
   *  \param nsds the siconos::modeling::NonSmoothDynamicalSystem
   *  \return a std::shared_ptr<Sensor> to the added Sensor
   */
  std::shared_ptr<Sensor> addAndRecordSensor(
      SensorType type, std::shared_ptr<siconos::simulation::TimeDiscretisation> td,
      std::shared_ptr<siconos::modeling::DynamicalSystem> ds,
      const siconos::modeling::NonSmoothDynamicalSystem& nsds);

  /** Add an existing Sensor to the Manager
   *
   *  \param s a std::shared_ptr<Sensor> to the Sensor we want to add
   *  \param td the TimeDiscretisation used for the associated Event
   */
  void addSensorPtr(std::shared_ptr<Sensor> s,
                    std::shared_ptr<siconos::simulation::TimeDiscretisation> td);

  /** To add, initialize an existing Sensor in the manager and record
   *  it in the simulation This function is only useful to add a new
   *  Sensor after the initialization of the manager else call
   *  addSensor()
   *
   *  \param s a std::shared_ptr<Sensor> to the Sensor we want to add
   *  \param td the TimeDiscretisation used for the associated Event
   *  \param nsds current nonsmooth dynamical system
   */
  void addAndRecordSensorPtr(std::shared_ptr<Sensor> s,
                             std::shared_ptr<siconos::simulation::TimeDiscretisation> td,
                             const siconos::modeling::NonSmoothDynamicalSystem& nsds);

  /** To build and add a new Actuator in the Manager
   *
   *  \param name the type of the Actuator
   *  \param td the std::shared_ptr<siconos::simulation::TimeDiscretisation> of the Actuator
   *  \param sensor the ControlSensor used to feed the Actuator
   *  \return the added Actuator
   */
  std::shared_ptr<Actuator> addActuator(
      ActuatorType name, std::shared_ptr<siconos::simulation::TimeDiscretisation> td,
      std::shared_ptr<ControlSensor> sensor);

  /** To build, add, initialize a new Actuator in the manager and
   *  record it in the simulation This function is only useful to add a
   *  new Actuator after the initialization of the manager else call
   *  addActuator()
   *
   *  \param name the type of the Actuator
   *  \param t the std::shared_ptr<siconos::simulation::TimeDiscretisation> of the Actuator
   *  \param sensor the ControlSensor used to feed the Actuator
   *  \param nsds the siconos::modeling::NonSmoothDynamicalSystem
   *  \return a std::shared_ptr<Actuator> to the added Actuator
   */
  std::shared_ptr<Actuator> addAndRecordActuator(
      ActuatorType name, std::shared_ptr<siconos::simulation::TimeDiscretisation> t,
      std::shared_ptr<ControlSensor> sensor,
      const siconos::modeling::NonSmoothDynamicalSystem& nsds);

  /** Add an existing Actuator to the manager
   *
   *  \param act a std::shared_ptr<Actuator> to the Actuator we want to add
   *  \param td the TimeDiscretisation used for the associated Event
   */
  void addActuatorPtr(std::shared_ptr<Actuator> act,
                      std::shared_ptr<siconos::simulation::TimeDiscretisation> td);

  /** To add, initialize an existing Actuator in the manager and record
   *  it in the simulation This function is only useful to add a new
   *  Actuator after the initialization of the manager otherwise call
   *  addActuator()
   *
   *  \param act a std::shared_ptr<Actuator> to the Actuator we want to add
   *  \param td the TimeDiscretisation used for the associated Event
   *  \param nsds current nonsmooth dynamical system
   */
  void addAndRecordActuatorPtr(std::shared_ptr<Actuator> act,
                               std::shared_ptr<siconos::simulation::TimeDiscretisation> td,
                               const siconos::modeling::NonSmoothDynamicalSystem& nsds);

  /** To build and add a new Observer in the Manager
   *
   *  \param name the type of the Observer
   *  \param td the std::shared_ptr<siconos::simulation::TimeDiscretisation> of the Observer
   *  \param sensor the ControlSensor feeding the Observer
   *  \param xHat0 the initial guess for the state
   *  \return a SP::ACtuator to the added Observer
   */
  std::shared_ptr<Observer> addObserver(
      ObserverType name, std::shared_ptr<siconos::simulation::TimeDiscretisation> td,
      std::shared_ptr<ControlSensor> sensor, const siconos::algebra::SiconosVector& xHat0);

  /** To build, add, initialize a new Observer in the manager and
   *  record it in the simulation This function is only useful to add a
   *  new Observer after the initialization of the manager else call
   *  addObserver()
   *
   *  \param name the type of the Observer
   *  \param td the std::shared_ptr<siconos::simulation::TimeDiscretisation> of the Observer
   *  \param sensor the ControlSensor feeding the Observer
   *  \param xHat0 the initial guess for the state
   *  \param nsds current nonsmooth dynamical system
   *  \return the added Observer
   */
  std::shared_ptr<Observer> addAndRecordObserver(
      ObserverType name, std::shared_ptr<siconos::simulation::TimeDiscretisation> td,
      std::shared_ptr<ControlSensor> sensor, const siconos::algebra::SiconosVector& xHat0,
      const siconos::modeling::NonSmoothDynamicalSystem& nsds);

  /** Add an existing Observer to the manager
   *
   *  \param obs a std::shared_ptr<Observer> to the Observer we want to add
   *  \param td the TimeDiscretisation used for the associated Event
   */
  void addObserverPtr(std::shared_ptr<Observer> obs,
                      std::shared_ptr<siconos::simulation::TimeDiscretisation> td);

  /** To add, initialize an existing Observer in the manager and record
   *  it in the simulation This function is only useful to add a new
   *  Observer after the initialization of the manager otherwise call
   *  addObserver()
   *
   *  \param obs a std::shared_ptr<Observer> to the Observer we want to add
   *  \param td the TimeDiscretisation used for the associated Event
   *  \param nsds current nonsmooth dynamical system
   */
  void addAndRecordObserverPtr(std::shared_ptr<Observer> obs,
                               std::shared_ptr<siconos::simulation::TimeDiscretisation> td,
                               const siconos::modeling::NonSmoothDynamicalSystem& nsds);

  /** initialize all Sensors, Observers and Actuators.
   *
   *  \param nsds current nonsmooth dynamical system
   */
  void initialize(const siconos::modeling::NonSmoothDynamicalSystem& nsds);

  /** display the data of the ControlManager on the standard output
   */
  void display() const;
};
}  // namespace siconos::control
#endif
