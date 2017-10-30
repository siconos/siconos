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
#include "LagrangianCompliantLinearTIR.hpp"
#include "Interaction.hpp"
//
#include "LagrangianDS.hpp"
#include "BlockVector.hpp"
#include "SimulationGraphs.hpp"

#include <iostream>

using namespace RELATION;

// Minimum data (C as pointer) constructor
LagrangianCompliantLinearTIR::LagrangianCompliantLinearTIR(SP::SimpleMatrix C, SP::SimpleMatrix D):
  LagrangianR(CompliantLinearTIR)
{
  _jachq = C;
  _jachlambda = D;
}

// Constructor from a complete set of data
LagrangianCompliantLinearTIR::LagrangianCompliantLinearTIR(SP::SimpleMatrix C, SP::SimpleMatrix D, SP::SimpleMatrix F, SP::SiconosVector e):
  LagrangianR(CompliantLinearTIR)
{
  _jachq = C;
  _jachlambda = D;
  _F = F;
  _e = e;
}

// Minimum data (C, e as pointers) constructor
LagrangianCompliantLinearTIR::LagrangianCompliantLinearTIR(SP::SimpleMatrix C, SP::SimpleMatrix D, SP::SiconosVector e):
  LagrangianR(CompliantLinearTIR)
{
  _jachq = C;
  _jachlambda = D;
  _e = e;
}

void LagrangianCompliantLinearTIR::initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM)
{
  unsigned int sizeY = inter.dimension();

  if (!(_jachq) || _jachq->size(1) !=  inter.getSizeOfDS() ||  _jachq->size(0) != sizeY)
    RuntimeException::selfThrow("LagrangianCompliantLinearTIR::initComponents inconsistent sizes between H matrix and the interaction.");

  if ((_jachlambda) && (_jachlambda->size(0) != sizeY || _jachlambda->size(1) != sizeY))
    RuntimeException::selfThrow("LagrangianCompliantLinearTIR::initComponents inconsistent sizes between D matrix and the interaction.");

  if ((_e) && _e->size() != sizeY)
    RuntimeException::selfThrow("LagrangianCompliantLinearTIR::initComponents inconsistent sizes between e vector and the dimension of the interaction.");

  unsigned int sizeZ = DSlink[LagrangianR::z]->size();
  if ((_F) && (
        _F->size(0) != sizeZ || _F->size(1) != sizeZ))
    RuntimeException::selfThrow("LagrangianCompliantLinearTIR::initComponents inconsistent sizes between F matrix and the interaction.");


}

void LagrangianCompliantLinearTIR::computeOutput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int derivativeNumber)
{
  // get y and lambda of the interaction
  SiconosVector& y = *inter.y(derivativeNumber);
  SiconosVector& lambda = *inter.lambda(derivativeNumber);
  VectorOfBlockVectors& DSlink = *interProp.DSlink;

  prod(*_jachq, *DSlink[LagrangianR::q0 + derivativeNumber], y);
  prod(*_jachlambda, lambda, y, false);

  if (derivativeNumber == 0)
  {
    if (_e)
      y += *_e;
    if (_F)
      prod(*_F, *DSlink[LagrangianR::z], y, false);
  }

}

void LagrangianCompliantLinearTIR::computeInput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level)
{
  // get lambda of the concerned interaction
  SiconosVector& lambda = *inter.lambda(level);
  VectorOfBlockVectors& DSlink = *interProp.DSlink;
  // computation of p = Ht lambda
  prod(lambda, *_jachq, *DSlink[LagrangianR::p0 + level], false);
}

void LagrangianCompliantLinearTIR::display() const
{
  LagrangianR::display();
  std::cout << "===== Lagrangian Linear Relation display ===== " <<std::endl;
  std::cout << " C: " <<std::endl;
  if (_jachq)
    _jachq->display();
  else
    std::cout << " -> NULL " <<std::endl;
  std::cout << " e: " <<std::endl;
  if (_e)
    _e->display();
  else
    std::cout << " -> NULL " <<std::endl;
  std::cout << " D: " <<std::endl;
  if (_jachlambda)
    _jachlambda->display();
  else
    std::cout << " -> NULL " <<std::endl;
  std::cout << " F: " <<std::endl;
  if (_F)
    _F->display();
  else
    std::cout << " -> NULL " <<std::endl;
  std::cout << "===================================== " <<std::endl;
}
