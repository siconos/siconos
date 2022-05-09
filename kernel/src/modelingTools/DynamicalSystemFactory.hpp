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

/*! \file DynamicalSystemFactory.hpp
\brief  Factory to generate user-defined DynamicalSystems
*/

#ifndef DynamicalSystemFactory_H
#define DynamicalSystemFactory_H

#include "SiconosPointers.hpp"
#include "DynamicalSystem.hpp"
#include <string>
#include <map>


/** Namespace for DynamicalSystem factory related objects. */
namespace DynamicalSystemFactory
{

/** A pointer to function, returning a pointer to DynamicalSystem,
    built with its type (related to the class name) and a vector of
    initial conditions. */
typedef SP::DynamicalSystem(*object_creator)(int, const SiconosVector&) ;

/** The type of the factory map */
typedef std::map<int, object_creator> MapFactory;

/** An iterator through the MapFactory */
typedef MapFactory::iterator MapFactoryIt;

/** Template function to return a new object of type SubType
 * \param name
 * \param x0
 * \return SP::DynamicalSystem
 */
template<class SubType> SP::DynamicalSystem factory(int name, const SiconosVector& x0)
{
  SP::DynamicalSystem res(new SubType(name, x0));
  return res;
}

/** Registry Class for DynamicalSystem.
 *
 * DynamicalSystem factory.
 * Use:

 *     DynamicalSystemFactory::Registry&
 *        regDynamicalSystem(DynamicalSystemFactory::Registry::get()) ;
 *
 *     SP::DynamicalSystem yourDynamicalSystem =
 *        regDynamicalSystem.instantiate(type, x0);
 */
class Registry
{

private :

  /** map that links a std::string, the type of the class, to a pointer
      to function, used to build the object. */
  MapFactory factory_map;

public :

  /** Access function to the Registry
   * \return static Registry&
   */
  static Registry& get() ;

  /** Add an object_creator into the factory_map, factory_map[name] = object.
   * \param i an int, the name of the object added
   * \param obj an object creator
   */
  void add(int i, object_creator obj);

  /** Function to instantiate a new DynamicalSystem
   * \param i an int, the name of the object added (type name!)
   * \param x0 the initial condition (SP)
   * \return  SP::DynamicalSystem
   */
  SP::DynamicalSystem instantiate(int i, const SiconosVector& x0);
} ;

/** Registration Class for sensors.
 *
 * Class used for auto-registration of DynamicalSystem-type objects.
 *
 */
class Registration
{
public :

  /** To register some new object into the factory
   * \param i an int, the name of the object to be registered
   * \param obj an object creator
   */
  Registration(int i, object_creator obj) ;
} ;

#define AUTO_REGISTER_DS(class_name,class_type) Registration _registration_## class_type(class_name,&factory<class_type>);
}
// end of namespace DynamicalSystemFactory

#endif
