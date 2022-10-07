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

/*! \file Sensor.hpp
  \brief General interface to define a sensor
*/

#ifndef Sensor_H
#define Sensor_H
#include <cassert>
#include <functional>
#include <map>
#include <memory>
#include <string>

#include "SiconosSerialization.hpp"

namespace siconos::algebra {
class SiconosVector;
class SimpleMatrix;
class SiconosMatrix;
}  // namespace siconos::algebra
namespace siconos::modeling {
class DynamicalSystem;
class NonSmoothDynamicalSystem;

};  // namespace siconos::modeling

namespace siconos::simulation {
class TimeDiscretisation;
};  // namespace siconos::simulation

namespace siconos::control {

/** Sensor types */
enum class SensorType {
  Linear,
};

/**
   Sensor Base Class

   Abstract class, interface to user-defined sensors.

   A Sensor is dedicated to data capture. It gives an interface for
   User who can implement its own Sensor to clearly define which data
   he needs to save.

   A Sensor handles a TimeDiscretisation, which defines the set of all
   instants where the sensor must operate (i.e. each times where
   capture() function will be called). An Event, inserted into the
   EventsManager of the Simulation, is linked to this
   TimeDiscretisation.

   Moreover, a Sensor is identified thanks to an id and a type (a
   number associated to the derived class type indeed).

   Construction

   To build a Sensor it is necessary to use the factory. Inputs are a
   number which identify the derived class type and a
   TimeDiscretisation:

   \code

   // Get the registry
   SensorFactory::Registry& regSensor(SensorFactory::Registry::get()) ;
   // Build a Sensor of type "myType" with t as a TimeDiscretisation.
   regSensor.instantiate(myType, t);

   \endcode

   The best way is to use the controlManager:
   \code
   // cm a ControlManager
   cm->addSensor(myType,t);
   // or if cm has already been initialized:
   cm->addAndRecordSensor(myType,t)
   \endcode


   The data are saved in a DataSet object named data, a map which
   associate to each Event another map.  This second map links a
   std::string, used to identify the data, and a SiconosVector.  As an
   example consider the case where you need to save the state vector x
   of a DynamicalSystem, then you can define a Data object, with
   "myDS_X" as an id and yourDS->x() as the SiconosVector. For
   myEvent being an Event where you need to save data, you get:
   (data[myEvent])["myDS_X"] =
   nonSmoothDynamicalSystem()->dynamicalSystem()->x()

   See users' manual for details on how to define its own Sensor.

 */
class Sensor {
 protected:
  ACCEPT_SERIALIZATION(Sensor);

  /** type of the Sensor */
  SensorType _type{SensorType::Linear};

  /** id of the Sensor */
  std::string _id{"none"};

  /** pointer to the DynamicalSystem we are measuring */
  std::shared_ptr<siconos::modeling::DynamicalSystem> _DS{nullptr};

  /** pointer to the state of the DynamicalSystem */
  std::shared_ptr<siconos::algebra::SiconosVector> _DSx{nullptr};

  // Rule of five
  Sensor() = delete;
  Sensor(const Sensor&) = delete;
  Sensor(Sensor&&) = delete;
  Sensor& operator=(const Sensor&) = delete;
  Sensor& operator=(Sensor&&) = delete;

 public:
  /** Constructor with a TimeDiscretisation.
   *
   *  \param type the type of the Sensor, which corresponds to the class type.
   *  \param ds the std::shared_ptr<siconos::modeling::DynamicalSystem> we observe.
   */
  Sensor(SensorType type, std::shared_ptr<siconos::modeling::DynamicalSystem> ds);

  /** destructor
   */
  virtual ~Sensor() noexcept = default;

  /** set id of the Sensor
   *
   *  \param newId the id of the Sensor
   */
  inline void setId(const std::string& newId) { _id = newId; };

  /** get id of the Sensor
   *
   *  \return a std::string
   */
  inline const std::string getId() const { return _id; };

  /**  \return the type of the Sensor
   *
   */
  inline SensorType getType() const { return _type; };

  /** get the DynamicalSystem linked to this Sensor
   *
   *  \return std::shared_ptr<siconos::modeling::DynamicalSystem>
   */
  inline std::shared_ptr<siconos::modeling::DynamicalSystem> getDS() const { return _DS; };

  /** This is derived in child classes if they need to copy the TimeDiscretisation
   *  associated with this Sensor
   *
   *  \param td the TimeDiscretisation for this Sensor
   */
  virtual void setTimeDiscretisation(const siconos::simulation::TimeDiscretisation& td){};

  /* get all the data saved for this sensor
   *  \return a DataSet
   */
  //  inline const DataSet getData() const
  //  {return _data;};

  /** initialize sensor data.
   *
   *  \param nsds the Model
   */
  virtual void initialize(const siconos::modeling::NonSmoothDynamicalSystem& nsds){};

  /** capture data when the SensorEvent is processed => set data[SensorEvent]=...
   */
  virtual void capture() = 0;

  /** display the data of the Sensor on the standard output
   */
  void display() const;
};

/** A class to handle sensors creation

  Requirements:
  - the Sensor type must be known and register

  For a XXXXSensor, add in the file describing the XXXXSensor class:

  static siconos::simulation::SensorRegistration<siconos::simulation::XXXSensor>
 reg(siconos::simulation::SensorType::XXXX);

  See SensorType enum for the available names.

  Usage:

  auto sensor = SensorFactory::instance()->create(time, SensorType::XXXX)

*/
class SensorFactory {
  // Signature of Sensor constructor
  using SensorCreator = std::function<std::shared_ptr<Sensor>(
      (std::shared_ptr<siconos::modeling::DynamicalSystem>))>;

  /** map to connect sensor type and the function used to create them */
  std::map<SensorType, SensorCreator> m_factories;

 public:
  /** Factory function which creates and returns an sensor

      \param  the dynamical system which is observed
      \param type type of the sensor (must be a SensorType (enum))
      \return a pointer to sensor
  */
  std::shared_ptr<Sensor> create(std::shared_ptr<siconos::modeling::DynamicalSystem> ds,
                                 SensorType type)
  {
    assert(m_factories.contains(type) && "unknown Sensor type");
    return m_factories[type](ds);
  }

  /** access to the (singleton) factory instance */
  static SensorFactory* instance()
  {
    static SensorFactory factory;
    return &factory;
  }

  void registerCreator(SensorType newtype, SensorCreator caller)
  {
    m_factories[newtype] = caller;
  }
};

template <class T>
class SensorRegistration {
 public:
  SensorRegistration(SensorType newtype)
  {
    SensorFactory::instance()->registerCreator(
        newtype, [](std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
          return std::make_shared<T>(ds);
        });
  }
};

}  // namespace siconos::control
#endif
