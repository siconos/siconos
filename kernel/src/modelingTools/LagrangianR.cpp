/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include "Interaction.hpp"
#include "LagrangianDS.hpp"

#include <iostream>

void LagrangianR::initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM)
{
  // do nothing here, overload this if you need something done
}

void LagrangianR::initialize(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM)
{
  // Memory allocation for G[i], if required (depends on the chosen constructor).
  unsigned int sizeY = inter.getSizeOfY();
  unsigned int sizeDS = inter.getSizeOfDS();

  if (! _jachq)
    _jachq.reset(new SimpleMatrix(sizeY, sizeDS));
  initComponents(inter, DSlink, workV, workM);
}

void LagrangianR::_zeroPlugin()
{
  Relation::_zeroPlugin();
  _pluginJachq.reset(new PluggedObject());
}
void LagrangianR::display() const
{
  Relation::display();
  std::cout << " _jachq :" << std::endl;
  if (_jachq)
    _jachq->display();
  std::cout << " _jachqDot :" << std::endl;
  if (_jachqDot)
    _jachqDot->display();
  std::cout << " _jachlambda :" << std::endl;
  if (_jachlambda)
    _jachlambda->display();
  else
    std::cout << " NULL :" << std::endl;

}
