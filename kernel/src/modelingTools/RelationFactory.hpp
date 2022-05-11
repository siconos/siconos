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

/*! \file RelationFactory.hpp
\brief  Factory to generate user-defined Relations
*/

#ifndef RelationFactory_H
#define RelationFactory_H

#include "SiconosPointers.hpp"
#include "Relation.hpp"
#include <string>
#include <map>

/** Namespace for Relation factory related objects. */
namespace RelationFactory
{

/** A pointer to function, returning a pointer to Relation, built with its type (related to the class name) */
typedef SP::Relation(*object_creator)(int) ;

/** The type of the factory map */
typedef std::map<int, object_creator> MapFactory;

/** An iterator through the MapFactory */
typedef MapFactory::iterator MapFactoryIt;

/** Template function to return a new object of type SubType
 * \param name an int, the name of the object added
 * \return SP::Relation
 */
template<class SubType> SP::Relation factory(int name)
{
  SP::Relation res(new SubType(name));
  return res;
}

/** Registry Class for sensors.
 *
 * Relation factory.
 * Use:
 *     RelationFactory::Registry& regRelation(RelationFactory::Registry::get()) ;
 *     SP::Relation yourRelation = regRelation.instantiate(type, x0);
 */
class Registry
{

private :

  /** map that links a std::string, the type of the class, to a pointer to function, used to build the object. */
  MapFactory factory_map;

public :

  /** Access function to the Registry 
   * \return Registry&
   */
  static Registry& get() ;

  /** Add an object_creator into the factory_map, factory_map[name] = object.
   * \param name an int, the name of the object added
   * \param obj an object creator
   */
  void add(int name, object_creator obj);

  /** Function to instantiate a new Relation
   * \param name an int, the name of the object added (type name!)
   * \return  SP::Relation
   */
  SP::Relation instantiate(int name);
} ;

/** Registration Class for relations.
 *
 * Class used for auto-registration of Relation-type objects.
 *
 */
class Registration
{

public :

  /** To register some new object into the factory
   * \param name an int, the name of the object to be registered
   * \param obj an object creator
   */
  Registration(int name, object_creator obj) ;
} ;

#define AUTO_REGISTER_RELATION(class_name,class_type) Registration _registration_## class_type(class_name,&factory<class_type>);
}
// end of namespace RelationFactory

#endif













