/* Siconos-Kernel version 1.2.0, Copyright INRIA 2005-2007.
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

/*! \file EventFactory.h
\brief  Factory to generate user-defined Events
*/

#ifndef EventFactory_H
#define EventFactory_H

# include "Event.h"
# include <string>
# include <map>

class Event;

/** Namespace for Events factory related objects. */
namespace EventFactory
{

/** A pointer to function, returning a pointer to Event, built with its type (ie class name) and a pointer to Model.*/
typedef Event* (*object_creator)(double, const std::string&);

/** The type of the factory map */
typedef std::map<const std::string, object_creator> MapFactory;

/** An iterator through the MapFactory */
typedef MapFactory::iterator MapFactoryIt;

/** Template function to return a new object of type SubType*/
template<class SubType> Event* factory(double time, const std::string& type)
{
  return new SubType(time, type);
}

/** Registry Class for Events.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
 *  \date (Creation) February 01, 2007
 *
 */
class Registry
{

private :

  /** map that links a string, the type of the class, to a pointer to function, used to build the object. */
  MapFactory factory_map;

public :

  /** */
  static Registry& get() ;

  /** Add an object_creator into the factory_map, factory_map[name] = object.
   * \param a string, the name of the object added
   * \param an object creator
   */
  void add(const std::string&, object_creator);

  /**
   *  \param a double, time of Event
   *  \param a string, type of Event
   */
  Event* instantiate(double, const std::string&);

} ;

/** Registration Class for Events.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
 *  \date (Creation) February 01, 2007
 *
 * Class used for auto-registration of Event-type objects.
 *
 */
class Registration
{

public :

  /** To register some new object into the factory
   * \param a string, the name of the object to be registered
   * \param an object creator
   */
  Registration(const std::string&, object_creator) ;
} ;

#define AUTO_REGISTER_EVENT(class_name,class_type) Registration _registration_## class_type(class_name,&factory<class_type>);
}
// end of namespace EventFactory

#endif













