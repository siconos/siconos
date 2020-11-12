/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

#include "RelationFactory.hpp"
#include "SiconosException.hpp"

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

  if(it == factory_map.end())
    THROW_EXCEPTION("Registry::instantiate (RelationFactory) failed, no class named: " + std::to_string(name));

  return (it->second)(name) ;  // run our factory
}

Registration::Registration(int name, object_creator creator)
{
  // Add creator into the factory of Relations
  Registry::get().add(name, creator) ;
}

}








