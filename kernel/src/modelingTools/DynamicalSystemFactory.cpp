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

#include "DynamicalSystemFactory.hpp"
#include "RuntimeException.hpp"

namespace DynamicalSystemFactory
{

Registry& Registry::get()
{
  static Registry instance;
  return instance;
}

void Registry::add(int name, object_creator creator)
{
  factory_map[name] = creator;
}

SP::DynamicalSystem Registry::instantiate(int name, const SiconosVector& x0)
{
  MapFactoryIt it = factory_map.find(name) ;

  if (it == factory_map.end())
    RuntimeException::selfThrow("Registry::instantiate (DynamicalSystemFactory) failed, no class named: " + name);

  return (it->second)(name, x0) ; // run our factory
}

Registration::Registration(int name, object_creator creator)
{
  // Add creator into the factory of DynamicalSystems
  Registry::get().add(name, creator) ;
}

}
