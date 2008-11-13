/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */

/*! \file RelationFactory.h
\brief  Factory to generate user-defined Relations
*/

#ifndef RelationFactory_H
#define RelationFactory_H

#include "SiconosPointers.hpp"
#include "Relation.h"
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

/** Template function to return a new object of type SubType*/
template<class SubType> SP::Relation factory(int name)
{
  SP::Relation res(new SubType(name));
  return res;
}

/** Registry Class for sensors.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) November 2008
 *
 * Relation factory.
 * Use:
 *     RelationFactory::Registry& regRelation(RelationFactory::Registry::get()) ;
 *     SP::Relation yourRelation = regRelation.instantiate(type, x0);
 */
class Registry
{

private :

  /** map that links a string, the type of the class, to a pointer to function, used to build the object. */
  MapFactory factory_map;

public :

  /** Access function to the Registry */
  static Registry& get() ;

  /** Add an object_creator into the factory_map, factory_map[name] = object.
   * \param an int, the name of the object added
   * \param an object creator
   */
  void add(int, object_creator);

  /** Function to instantiate a new Relation
   * \param an int, the name of the object added (type name!)
   */
  SP::Relation instantiate(int);
} ;

/** Registration Class for sensors.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) November 2008
 *
 * Class used for auto-registration of Relation-type objects.
 *
 */
class Registration
{

public :

  /** To register some new object into the factory
   * \param an int, the name of the object to be registered
   * \param an object creator
   */
  Registration(int, object_creator) ;
} ;

#define AUTO_REGISTER_RELATION(class_name,class_type) Registration _registration_## class_type(class_name,&factory<class_type>);
}
// end of namespace RelationFactory

#endif













