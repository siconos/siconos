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

void LagrangianR::zeroPlugin()
{
  Relation::zeroPlugin();
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
