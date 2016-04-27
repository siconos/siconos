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
#include "FirstOrderLinearTIR.hpp"
#include "Interaction.hpp"
#include "BlockVector.hpp"
#include "SimulationGraphs.hpp"

#include <iostream>

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
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

void FirstOrderLinearTIR::initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM)
{
  DEBUG_PRINT("FirstOrderLinearTIR::initialize(Interaction & inter)\n");
  // Note: do not call FirstOrderR::initialize to avoid jacobianH and jacobianG allocation.

  if (!_C)
    RuntimeException::selfThrow("FirstOrderLinearTIR::initialize() C is null and is a required input.");
  if (!_B)
    RuntimeException::selfThrow("FirstOrderLinearTIR::initialize() B is null and is a required input.");

  // Check if various operators sizes are consistent.
  // Reference: interaction.

  DEBUG_PRINTF("_C->size(0) = %i,\t inter.getSizeOfY() = %i\n ",_C->size(0),inter.getSizeOfY() );
  DEBUG_PRINTF("_C->size(1) = %i,\t inter.getSizeOfDS() = %i\n ",_C->size(1),inter.getSizeOfDS() );

  assert((_C->size(0) == inter.getSizeOfY() && _C->size(1) == inter.getSizeOfDS()) && "FirstOrderLinearTIR::initialize , inconsistent size between C and Interaction.");

  assert((_B->size(1) == inter.getSizeOfY() && _B->size(0) ==  inter.getSizeOfDS()) && "FirstOrderLinearTIR::initialize , inconsistent size between B and interaction.");

  // C and B are the minimum inputs. The others may remain null.

  if (_D)
    assert((_D->size(0) == inter.getSizeOfY() || _D->size(1) == inter.getSizeOfY()) && "FirstOrderLinearTIR::initialize , inconsistent size between C and D.");


  if (_F)
    assert(((_F->size(0) != inter.getSizeOfY()) && (_F->size(1) != DSlink[FirstOrderR::z]->size())) && "FirstOrderLinearTIR::initialize , inconsistent size between C and F.");
  if (_e)
    assert(_e->size() == inter.getSizeOfY() && "FirstOrderLinearTIR::initialize , inconsistent size between C and e.");
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

void FirstOrderLinearTIR::computeOutput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level)
{
  // We get y and lambda of the interaction (pointers)
  SiconosVector& y = *inter.y(0);
  SiconosVector& lambda = *inter.lambda(0);

  VectorOfBlockVectors& DSlink = *interProp.DSlink;
  computeh(*DSlink[FirstOrderR::x], lambda, *DSlink[FirstOrderR::z], y);
}

void FirstOrderLinearTIR::computeg(SiconosVector& lambda, BlockVector& r)
{
  prod(*_B, lambda, r, false);
}

void FirstOrderLinearTIR::computeInput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level)
{
  VectorOfBlockVectors& DSlink = *interProp.DSlink;
  computeg(*inter.lambda(level), *DSlink[FirstOrderR::r]);
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
