/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include "SimpleMatrix.hpp"
#include "FirstOrderLinearTIR.hpp"
#include "Interaction.hpp"
#include "BlockVector.hpp"
#include "SimulationGraphs.hpp"
#include <iostream>
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"

using namespace RELATION;


FirstOrderLinearTIR::FirstOrderLinearTIR():
  FirstOrderR(LinearTIR)
{
}

// Minimum data (C, B as pointers) constructor
FirstOrderLinearTIR::FirstOrderLinearTIR(SP::SimpleMatrix C, SP::SimpleMatrix B):
  FirstOrderR(LinearTIR)
{
  _C = C;
  _B = B;
}

// Constructor from a complete set of data
FirstOrderLinearTIR::FirstOrderLinearTIR(SP::SimpleMatrix C, SP::SimpleMatrix D, SP::SimpleMatrix F, SP::SiconosVector e, SP::SimpleMatrix B):
  FirstOrderR(LinearTIR)
{
  _C = C;
  _B = B;
  _D = D;
  _F = F;
  _e = e;
}

void FirstOrderLinearTIR::initialize(Interaction& inter)
{
  DEBUG_PRINT("FirstOrderLinearTIR::initialize(Interaction & inter)\n");

  FirstOrderR::initialize(inter); // ?

  if (!_C)
    RuntimeException::selfThrow("FirstOrderLinearTIR::initialize() C is null and is a required input.");
  if (!_B)
    RuntimeException::selfThrow("FirstOrderLinearTIR::initialize() B is null and is a required input.");

  checkSize(inter);
}

void FirstOrderLinearTIR::checkSize(Interaction& inter)
{
  DEBUG_PRINT("FirstOrderLinearTIR::checkSize(Interaction & inter)\n");
  DEBUG_PRINTF("_C->size(0) = %i,\t inter.dimension() = %i\n ",_C->size(0),inter.dimension() );
  DEBUG_PRINTF("_C->size(1) = %i,\t inter.getSizeOfDS() = %i\n ",_C->size(1),inter.getSizeOfDS() );

  assert((_C->size(0) == inter.dimension() && _C->size(1) == inter.getSizeOfDS()) && "FirstOrderLinearTIR::initialize , inconsistent size between C and Interaction.");

  assert((_B->size(1) == inter.dimension() && _B->size(0) ==  inter.getSizeOfDS()) && "FirstOrderLinearTIR::initialize , inconsistent size between B and interaction.");

  // C and B are the minimum inputs. The others may remain null.

  if (_D)
    assert((_D->size(0) == inter.dimension() || _D->size(1) == inter.dimension()) && "FirstOrderLinearTIR::initialize , inconsistent size between C and D.");


  if (_F)
    assert(((_F->size(0) != inter.dimension()) && (_F->size(1) != (inter.linkToDSVariables())[FirstOrderR::z]->size())) && "FirstOrderLinearTIR::initialize , inconsistent size between C and F.");
  if (_e)
    assert(_e->size() == inter.dimension() && "FirstOrderLinearTIR::initialize , inconsistent size between C and e.");

}

void FirstOrderLinearTIR::computeh(BlockVector& x, SiconosVector& lambda, BlockVector& z, SiconosVector& y)
{

  if (_C)
    prod(*_C, x, y, true);
  else
    y.zero();

  if (_D)
    prod(*_D, lambda, y, false);

  if (_e)
    y += *_e;

  if (_F)
    prod(*_F, z, y, false);

}

void FirstOrderLinearTIR::computeOutput(double time, Interaction& inter, unsigned int level)
{
  // We get y and lambda of the interaction (pointers)
  SiconosVector& y = *inter.y(level);
  SiconosVector& lambda = *inter.lambda(level);

  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  computeh(*DSlink[FirstOrderR::x], lambda, *DSlink[FirstOrderR::z], y);
}

void FirstOrderLinearTIR::computeg(SiconosVector& lambda, BlockVector& r)
{
  prod(*_B, lambda, r, false);
}

void FirstOrderLinearTIR::computeInput(double time, Interaction& inter, unsigned int level)
{
  DEBUG_BEGIN("FirstOrderLinearTIR::computeInput(double time, Interaction& inter, unsigned int level)\n")
  VectorOfBlockVectors& DSlink = inter.linkToDSVariables();
  DEBUG_EXPR(inter.lambda(level)->display(););
  DEBUG_EXPR(DSlink[FirstOrderR::r]->display(););
  computeg(*inter.lambda(level), *DSlink[FirstOrderR::r]);
  DEBUG_END("FirstOrderLinearTIR::computeInput(double time, Interaction& inter, unsigned int level)\n")
}

void FirstOrderLinearTIR::display() const
{
  std::cout << " ===== Linear Time Invariant relation display ===== " <<std::endl;
  std::cout << "| C " <<std::endl;
  if (_C) _C->display();
  else std::cout << "->NULL" <<std::endl;
  std::cout << "| D " <<std::endl;
  if (_D) _D->display();
  else std::cout << "->NULL" <<std::endl;
  std::cout << "| F " <<std::endl;
  if (_F) _F->display();
  else std::cout << "->NULL" <<std::endl;
  std::cout << "| e " <<std::endl;
  if (_e) _e->display();
  else std::cout << "->NULL" <<std::endl;
  std::cout << "| B " <<std::endl;
  if (_B) _B->display();
  else std::cout << "->NULL" <<std::endl;
  std::cout << " ================================================== " <<std::endl;
}
