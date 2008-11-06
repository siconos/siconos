/* Siconos-Kernel version 1.2.0, Copyright INRIA 2005-2008.
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

#include "RelationFactory.h"
#include "RuntimeException.h"

namespace RelationFactory
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

SP::Relation Registry::instantiate(int name)
{
  MapFactoryIt it = factory_map.find(name) ;

  if (it == factory_map.end())
    RuntimeException::selfThrow("Registry::instantiate (RelationFactory) failed, no class named: " + name);

  return (it->second)(name) ;  // run our factory
}

Registration::Registration(int name, object_creator creator)
{
  // Add creator into the factory of Relations
  Registry::get().add(name, creator) ;
}

}








