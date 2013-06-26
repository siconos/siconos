/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
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

/** Template function to return a new object of type SubType*/
template<class SubType> SP::Observer factory(SP::ControlSensor sensor, const SiconosVector& xHat0)
{
  return std11::shared_ptr<SubType>(new SubType(sensor, xHat0));
}

/** Registry Class for Observers.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.6.0.
 *  \date (Creation) June 15, 2013
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

  /** Access function to the Registry */
  static Registry& get() ;

  /** Add an object_creator into the factory_map, factory_map[name] = object.
   * \param type the type of the added Observer
   * \param creator object creator
   */
  void add(unsigned int type, object_creator object);

  /** Function to instantiate a new Observer
   * \param type the type of the Observer we want to instantiate
   * \param t a SP::TimeDiscretisation.
   * \param sensor the ControlSensor feeding this Observer
   * \param xHat0 the original estimate
   * \return a SP::Observer to the created Observer
   */
  SP::Observer instantiate(unsigned int type, SP::ControlSensor sensor, const SiconosVector& xHat0);

} ;

/** Registration Class for Observers.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.6.0.
 *  \date (Creation) June 15, 2013
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
  Registration(unsigned int type, object_creator object) ;
} ;

#define AUTO_REGISTER_OBSERVER(class_name,class_type) ObserverFactory::Registration _registration_## class_type(class_name, &ObserverFactory::factory<class_type>);
}
// end of namespace ObserverFactory

#endif
