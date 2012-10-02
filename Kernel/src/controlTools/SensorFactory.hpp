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
typedef SP::Sensor(*object_creator)(SP::TimeDiscretisation, SP::DynamicalSystem) ;

/** The type of the factory map */
typedef std::map<int, object_creator> MapFactory;

/** An iterator through the MapFactory */
typedef MapFactory::iterator MapFactoryIt;

/** Template function to return a new object of type SubType*/
template<class SubType> SP::Sensor factory(SP::TimeDiscretisation t, SP::DynamicalSystem ds)
{
  return std11::shared_ptr<SubType>(new SubType(t, ds));
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
 * With sensorType a string, the name of the class of your Sensor (expl: "SensorPosition"), timeD a SP::TimeDiscretisation and
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
   * \param name the type of the added Sensor
   * \param creator object creator
   */
  void add(int name, object_creator object);

  /** Function to instantiate a new Sensor
   * \param name the type of the Sensor we want to instantiate
   * \param t a SP::TimeDiscretisation.
   * \param ds a SP::DynamicalSystem that will be linked to this Sensor
   * \param a SP::Sensor to the created Sensor
   */
  SP::Sensor instantiate(int name, SP::TimeDiscretisation t, SP::DynamicalSystem ds);

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
   * \param name the type of the added Sensor
   * \param creator object creator
   */
  Registration(int name, object_creator object) ;
} ;

#define AUTO_REGISTER_SENSOR(class_name,class_type) Registration _registration_## class_type(class_name,&factory<class_type>);
}
// end of namespace SensorFactory

#endif
