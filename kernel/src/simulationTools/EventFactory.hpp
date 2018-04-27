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

/*! \file EventFactory.hpp
\brief  Factory to generate user-defined Events
*/

#ifndef EventFactory_H
#define EventFactory_H

# include "Event.hpp"
# include <map>


/** Namespace for Events factory related objects. */
namespace EventFactory
{

/** A pointer to function, returning a pointer to Event, built with its type (ie class name) and a pointer to Model.*/
typedef SP::Event(*object_creator)(double, int);

/** The type of the factory map */
typedef std::map<int, object_creator> MapFactory;

/** An iterator through the MapFactory */
typedef MapFactory::iterator MapFactoryIt;

/** Template function to return a new object of type SubType
 * \param time time of the Event
 * \param type type of the Event
 * \return an Event
 */
template<class SubType> SP::Event factory(double time, int type)
{
  SP::Event e(new SubType(time, type));
  return e;
}

/** Registry Class for Events.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) February 01, 2007
 *
 */
class Registry
{

private :

  /** map that links a std::string, the type of the class, to a pointer to function, used to build the object. */
  MapFactory factory_map;

public :

  /** get access to the Registry
   * \return reference to the registry
   */
  static Registry& get() ;

  /** Add an object_creator into the factory_map, factory_map[name] = object.
   * \param type the type of the object added
   * \param creator object creator
   */
  void add(int type, object_creator creator);

  /**
   *  \param time time of Event
   *  \param type type of Event
   *  \return an Event
   */
  SP::Event instantiate(double time, int type);

} ;

/** Registration Class for Events.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) February 01, 2007
 *
 * Class used for auto-registration of Event-type objects.
 *
 */
class Registration
{

public :

  /** To register some new object into the factory
   * \param type the type of the object added
   * \param creator object creator
   */
  Registration(int type, object_creator creator);
} ;

}
// end of namespace EventFactory

#define AUTO_REGISTER_EVENT(class_name,class_type) EventFactory::Registration _registration_## class_type(class_name,&EventFactory::factory<class_type>);
#endif













