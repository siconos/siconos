/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

/*! \file SensorFactory.hpp
  Factory to generate user-defined sensors
*/

#ifndef SensorFactory_H
#define SensorFactory_H

#include<map>

#include "Sensor.hpp"
#include "TimeDiscretisation.hpp"

/** Namespace for Sensor factory related objects. */
namespace SensorFactory
{

/** A pointer to function, returning a SP::Sensor, built with its type (ie class name) and a SP::DynamicalSystem.*/
typedef SP::Sensor(*object_creator)(SP::DynamicalSystem) ;

/** The type of the factory map */
typedef std::map<int, object_creator> MapFactory;

/** An iterator through the MapFactory */
typedef MapFactory::iterator MapFactoryIt;

/** Template function to return a new object of type SubType
 * \param ds the DynamicalSystem for the Sensor
 * \return a Sensor
 */
template<class SubType> SP::Sensor factory(SP::DynamicalSystem ds)
{
  return std11::shared_ptr<SubType>(new SubType( ds));
}

/** Registry Class for sensors.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) February 01, 2007
 *
 * Sensor factory.
 * Use:
 *     SensorFactory::Registry& regSensor(SensorFactory::Registry::get()) ;
 *     SP::Sensor yourSensor = regSensor.instantiate(sensorType, timeD, myDS);
 * With sensorType a std::string, the name of the class of your Sensor (expl: "SensorPosition"), timeD a SP::TimeDiscretisation and
 * myDS a SP::DynamicalSystem.
 *
 */
class Registry
{

private :

  /** map that links a std::string, the type of the class, to a pointer
      to function, used to build the object. */
  MapFactory factory_map;

public :

  /** Access function to the Registry
   * \return a reference to the registry
   */
  static Registry& get() ;

  /** Add an object_creator into the factory_map
   * \param type the type of the added Sensor
   * \param creator object creator
   */
  void add(int type, object_creator creator);

  /** Function to instantiate a new Sensor
   * \param type the type of the Sensor we want to instantiate
   * \param ds a SP::DynamicalSystem that will be linked to this Sensor
   * \return a Sensor
   */
  SP::Sensor instantiate(int type, SP::DynamicalSystem ds);

} ;

/** Registration Class for sensors.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) February 01, 2007
 *
 * Class used for auto-registration of Sensor-type objects.
 *
 */
class Registration
{

public :

  /** To register some new object into the factory
   * \param type the type of the added Sensor
   * \param creator object creator
   */
  Registration(int type, object_creator creator) ;
} ;

#define AUTO_REGISTER_SENSOR(class_name, class_type) SensorFactory::Registration _registration_## class_type(class_name, &SensorFactory::factory<class_type>);
}
// end of namespace SensorFactory

#endif
