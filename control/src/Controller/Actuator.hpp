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

/*! \file Actuator.hpp
  \brief General interface to define an actuator
*/

#ifndef Actuator_H
#define Actuator_H

#include <cassert>
#include <functional>
#include <map>
#include <memory>
#include <string>

#include "SiconosSerialization.hpp"

namespace siconos::algebra {
class SiconosVector;
class SimpleMatrix;
}  // namespace siconos::algebra

namespace siconos::modeling {
class NonSmoothDynamicalSystem;
}

namespace siconos::simulation {
class TimeDiscretisation;
class Simulation;
}  // namespace siconos::simulation

namespace siconos::control {

class ControlSensor;

enum class ActuatorType {

  PID,
  LinearSMC,
  ExplicitLinearSMC,
  LinearSMCOT2,
  LinearSMCimproved,
  Twisting,
  RegularTwisting,
  ExplicitTwisting
};

/**
   Actuators Base Class

   Abstract class, interface to user-defined actuators.

   An Actuator is dedicated to act on parameters of the Model
   (especially z param. in DynamicalSystem) according to some specific
   values recorded thanks to sensors. It gives an interface for User
   who can implement its own Actuator.  clearly define which data he
   needs to save.

   An Actuator handles a TimeDiscretisation, which defines the set of
   all instants where the Actuator must operate (i.e. each times
   where actuate() function will be called). An Event, inserted into
   the EventsManager of the Simulation, is linked to this
   TimeDiscretisation.

   Moreover, an Actuator is identified thanks to an id and a type (a
   number associated to the derived class type indeed).

   Construction

   To build an Actuator it is necessary to use the factory. Inputs are
   a number which identify the derived class type and a
   TimeDiscretisation:

   \code
   ActuatorFactory::instance()->create(sensor, type)
   \endcode

   The best way is to use the controlManager:
   \code
   // cm a ControlManager
   cm->addActuator(myType,t);
   // or if cm has already been initialized:
   cm->addAndRecordActuator(myType,t)
   \endcode
*/
class Actuator {
 private:
  ACCEPT_SERIALIZATION(Actuator);

 protected:
  /** type of the Actuator */
  ActuatorType _type{ActuatorType::PID};

  /** id of the Actuator */
  std::string _id{"none"};

  /** Control variable */
  std::shared_ptr<siconos::algebra::SiconosVector> _u{nullptr};

  /** B Matrix */
  std::shared_ptr<siconos::algebra::SimpleMatrix> _B{nullptr};

  /** name of the plugin for g (nonlinear affine in control system)*/
  std::string _plugingName;

  /** name of the plugin to compute \f$ \nabla_x g \f$ for the nonlinear case*/
  std::string _pluginJacgxName;

  /** ControlSensor feeding the Controller */
  std::shared_ptr<ControlSensor> _sensor{nullptr};

  // Rule of five
  Actuator() = delete;
  Actuator(const Actuator&) = delete;
  Actuator(Actuator&&) = delete;
  Actuator& operator=(const Actuator&) = delete;
  Actuator& operator=(Actuator&&) = delete;

 public:
  /** General Constructor
   *
   *  \param type the type of the Actuator, which corresponds to the class type
   *  \param sensor the ControlSensor feeding the Actuator
   */
  Actuator(ActuatorType type, std::shared_ptr<ControlSensor> sensor);

  /** General Constructor with dynamics affine in control
   *
   *  \param type the type of the Actuator, which corresponds to the class type
   *  \param sensor the ControlSensor feeding the Actuator
   */
  Actuator(ActuatorType type, std::shared_ptr<ControlSensor> sensor,
           std::shared_ptr<siconos::algebra::SimpleMatrix> B);

  /** destructor
   */
  virtual ~Actuator() noexcept = default;

  /** set id of the Actuator
   *
   *  \param newId the new id.
   */
  inline void setId(const std::string& newId) { _id = newId; };

  /** get id of the Actuator
   *
   *  \return a std::string
   */
  inline const std::string getId() const { return _id; };

  /** get the type of the Actuator (ie class name)
   *
   *  \return an integer
   */
  inline ActuatorType getType() const { return _type; };

  /** Get the control value
   *
   *  \return current control value u
   */
  inline const siconos::algebra::SiconosVector& u() const { return *_u; };

  /** Set the control size
   *
   *  \param size dimension of the control input u
   */
  void setSizeu(unsigned size);

  /** Set the B matrix
   *
   *  \param B the new B matrix
   */
  inline void setB(std::shared_ptr<siconos::algebra::SimpleMatrix> B) { _B = B; };

  /** Set the name of the plugin for computing g
   *
   *  \param g the name of the plugin to compute g
   */
  inline void setg(const std::string& g) { _plugingName = g; };

  /** add a Sensor in the actuator.
   *
   *  \param newSensor a Sensor that will be connected to the Actuator
   */
  void addSensorPtr(std::shared_ptr<ControlSensor> newSensor);

  /** This is derived in child classes if they need to copy the TimeDiscretisation
   *  associated with this Actuator
   *
   *  \param td the TimeDiscretisation for this Actuator
   */
  virtual void setTimeDiscretisation(const siconos::simulation::TimeDiscretisation& td){};

  /** initialize actuator data.
   *
   *  \param nsds the siconos::modeling::NonSmoothDynamicalSystem
   *  \param s the simulation
   */
  virtual void initialize(const siconos::modeling::NonSmoothDynamicalSystem& nsds,
                          const siconos::simulation::Simulation& s);

  /** capture data when the ActuatorEvent is processed
   */
  virtual void actuate() = 0;

  /** display the data of the Actuator on the standard output
   */
  virtual void display() const;

  /** get the NSDS used in the Controller, if there is one
   *
   *  \return "NULL" shared_ptr if there is no internal simulation, otherwise
   *  it returns the Model hoding the simulation
   */
  virtual std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> getInternalNSDS() const;
};

/** A class to handle actuators creation

  Requirements:
  - the Actuator type must be known and register

  For a XXXXActuator, add in the file describing the XXXXActuator class:

  static siconos::simulation::ActuatorRegistration<siconos::simulation::XXXActuator>
 reg(siconos::simulation::ActuatorType::XXXX);

  See ActuatorType enum for the available names.

  Usage:

  auto actuator = ActuatorFactory::instance()->create(sensor, ActuatorType::XXXX)

*/
class ActuatorFactory {
  // Signature of Actuator constructor
  using ActuatorCreator =
      std::function<std::shared_ptr<Actuator>(std::shared_ptr<ControlSensor>)>;

  /** map to connect actuator type and the function used to create them */
  std::map<ActuatorType, ActuatorCreator> m_factories;

 public:
  /** Factory function which creates and returns an actuator

      \param sensor the ControlSensor used by the Actuator
      \param type type of the actuator (must be a ActuatorType (enum))
      \return a pointer to actuator
  */
  std::shared_ptr<Actuator> create(std::shared_ptr<ControlSensor> sensor, ActuatorType type)
  {
    assert(m_factories.contains(type) && "unknown Actuator type");
    return m_factories[type](sensor);
  }

  /** access to the (singleton) factory instance */
  static ActuatorFactory* instance()
  {
    static ActuatorFactory factory;
    return &factory;
  }

  void registerCreator(ActuatorType newtype, ActuatorCreator caller)
  {
    m_factories[newtype] = caller;
  }
};

template <class T>
class ActuatorRegistration {
 public:
  ActuatorRegistration(ActuatorType newtype)
  {
    ActuatorFactory::instance()->registerCreator(
        newtype,
        [](std::shared_ptr<ControlSensor> sensor) { return std::make_shared<T>(sensor); });
  }
};

}  // namespace siconos::control
#endif
