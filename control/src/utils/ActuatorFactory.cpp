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

#include "ActuatorFactory.hpp"
#include "SiconosException.hpp"



namespace ActuatorFactory
{

Registry& Registry::get()
{
  static Registry instance;
  return instance;
}

void Registry::add(unsigned int type, object_creator creator)
{
  factory_map[type] = creator;
}

SP::Actuator Registry::instantiate(unsigned int type, SP::ControlSensor sensor)
{
  MapFactoryIt it = factory_map.find(type);

  if(it == factory_map.end())
    THROW_EXCEPTION("Registry::instantiate (ActuatorFactory) \
        failed, no class numbered: " + std::to_string(type));

  // std::cout <<std::endl << "Factory instance for class" << name <<std::endl ; // for test purposes only
  return (it->second)(sensor) ;  // run our factory
}

Registration::Registration(unsigned int type, object_creator creator)
{
  //  std::cout <<std::endl << "Registration of " << name <<std::endl <<std::endl ;
  // Add creator into the factory of Actuators
  Registry::get().add(type, creator) ;
}

}
