/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

// \todo : create a work vector for all tmp vectors used in computeg, computeh ...

#include "LagrangianR.hpp"
#include "PluggedObject.hpp"
#include "SimpleMatrix.hpp"

#include <iostream>

void LagrangianR::_zeroPlugin()
{
  Relation::_zeroPlugin();
  _pluginJachq = std::make_shared<PluggedObject>();
}

void LagrangianR::display() const
{
  Relation::display();
  std::cout << " _jachq :" << "\n";
  if(_jachq)
    _jachq->display();
  std::cout << " _jachqDot :" << "\n";
  if(_jachqDot)
    _jachqDot->display();
  std::cout << " _jachlambda :" << "\n";
  if(_jachlambda)
    _jachlambda->display();
  else
    std::cout << " nullptr :" << "\n";

}
