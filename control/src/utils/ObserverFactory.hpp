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

/*! \file ObserverFactory.hpp
  Factory to generate user-defined Observers
*/

#ifndef ObserverFactory_H
#define ObserverFactory_H

#include <map>

#include "Observer.hpp"
#include "TimeDiscretisation.hpp"

/** Namespace for Observer factory related objects. */
namespace ObserverFactory
{

/** A pointer to function, returning a SP::Observer, built with its type (ie class name) a ControlSensor
 * and a SiconosVector, the initial estimate*/
typedef SP::Observer(*object_creator)(SP::ControlSensor, const SiconosVector&) ;

/** The type of the factory map */
typedef std::map<unsigned int, object_creator> MapFactory;

/** An iterator through the MapFactory */
typedef MapFactory::iterator MapFactoryIt;

/** Template function to return a new object of type SubType
 * \param sensor ControlSensor used by this Observer
 * \param xHat0 initial state estimate
 * \return an Observer
 */
template<class SubType> SP::Observer factory(SP::ControlSensor sensor, const SiconosVector& xHat0)
{
  return std::shared_ptr<SubType>(new SubType(sensor, xHat0));
}

/** Registry Class for Observers.
 *
 * Observer factory.
 * Use:
 *     ObserverFactory::Registry& regObserver(ObserverFactory::Registry::get()) ;
 *     SP::Observer yourObserver = regObserver.instantiate(observerType, timeD, myDS);
 * With observerType a string, the name of the class of your Observer (expl: "ObserverPosition"), timeD a SP::TimeDiscretisation and
 * myDS a SP::DynamicalSystem.
 *
 */
class Registry
{

private :

  /** map that links a string, the type of the class, to a pointer
      to function, used to build the object. */
  MapFactory factory_map;

public :

  /** Access function to the Registry
   * \return reference to the registry
   */
  static Registry& get() ;

  /** Add an object_creator into the factory_map, factory_map[name] = object.
   * \param type the type of the added Observer
   * \param creator object creator
   */
  void add(unsigned int type, object_creator creator);

  /** Function to instantiate a new Observer
   * \param type the type of the Observer we want to instantiate
   * \param sensor the ControlSensor feeding this Observer
   * \param xHat0 the original estimate
   * \return a SP::Observer to the created Observer
   */
  SP::Observer instantiate(unsigned int type, SP::ControlSensor sensor, const SiconosVector& xHat0);

} ;

/** Registration Class for Observers.
 *
 * Class used for auto-registration of Observer-type objects.
 *
 */
class Registration
{

public :

  /** To register some new object into the factory
   * \param type the type of the added Observer
   * \param creator object creator
   */
  Registration(unsigned int type, object_creator creator) ;
} ;

#define AUTO_REGISTER_OBSERVER(class_name,class_type) ObserverFactory::Registration _registration_## class_type(class_name, &ObserverFactory::factory<class_type>);
}
// end of namespace ObserverFactory

#endif
