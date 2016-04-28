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
#include "LagrangianLinearTIR.hpp"
#include "Interaction.hpp"
//
#include "LagrangianDS.hpp"
#include "BlockVector.hpp"
#include "SimulationGraphs.hpp"

#include <iostream>

using namespace RELATION;

// Minimum data (C as pointer) constructor
LagrangianLinearTIR::LagrangianLinearTIR(SP::SimpleMatrix C):
  LagrangianR(LinearTIR)
{
  _jachq = C;
}

// Constructor from a complete set of data
LagrangianLinearTIR::LagrangianLinearTIR(SP::SimpleMatrix C, SP::SimpleMatrix D, SP::SimpleMatrix F, SP::SiconosVector e):
  LagrangianR(LinearTIR)
{
  _jachq = C;
  _jachlambda = D;
  _F = F;
  _e = e;
}

// Minimum data (C, e as pointers) constructor
LagrangianLinearTIR::LagrangianLinearTIR(SP::SimpleMatrix C, SP::SiconosVector e):
  LagrangianR(LinearTIR)
{
  _jachq = C;
  _e = e;
}

void LagrangianLinearTIR::initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM)
{
  unsigned int sizeY = inter.getSizeOfY();

  if (!(_jachq) || _jachq->size(1) !=  inter.getSizeOfDS() ||  _jachq->size(0) != sizeY)
    RuntimeException::selfThrow("LagrangianLinearTIR::initComponents inconsistent sizes between H matrix and the interaction.");

  if ((_jachlambda) && (_jachlambda->size(0) != sizeY || _jachlambda->size(1) != sizeY))
    RuntimeException::selfThrow("LagrangianLinearTIR::initComponents inconsistent sizes between D matrix and the interaction.");

  if ((_e) && _e->size() != sizeY)
    RuntimeException::selfThrow("LagrangianLinearTIR::initComponents inconsistent sizes between e vector and the dimension of the interaction.");

  unsigned int sizeZ = DSlink[LagrangianR::z]->size();
  if ((_F) && (
        _F->size(0) != sizeZ || _F->size(1) != sizeZ))
    RuntimeException::selfThrow("LagrangianLinearTIR::initComponents inconsistent sizes between F matrix and the interaction.");


}

void LagrangianLinearTIR::computeOutput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int derivativeNumber)
{
  // get y and lambda of the interaction
  SiconosVector& y = *inter.y(derivativeNumber);
  VectorOfBlockVectors& DSlink = *interProp.DSlink;

  prod(*_jachq, *DSlink[LagrangianR::q0 + derivativeNumber], y);

  if (derivativeNumber == 0)
  {
    if (_e)
      y += *_e;
    if (_F)
      prod(*_F, *DSlink[LagrangianR::z], y, false);
  }

  if (_jachlambda)
  {
    SiconosVector& lambda = *inter.lambda(derivativeNumber);
    prod(*_jachlambda, lambda, y, false);
  }


}

void LagrangianLinearTIR::computeInput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level)
{
  // get lambda of the concerned interaction
  SiconosVector& lambda = *inter.lambda(level);
  VectorOfBlockVectors& DSlink = *interProp.DSlink;
  // computation of p = Ht lambda
  prod(lambda, *_jachq, *DSlink[LagrangianR::p0 + level], false);
}

void LagrangianLinearTIR::display() const
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
