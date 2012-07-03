/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
#include "RelationXML.hpp"
#include "Interaction.hpp"
#include "LagrangianDS.hpp"

using namespace std;

void LagrangianR::initComponents(Interaction& inter)
{
  unsigned int sizeY = inter.getSizeOfY();
  unsigned int sizeDS = inter.getSizeOfDS();

  // The initialization of Jach[0] depends on the way the Relation was built ie if the matrix
  // was read from xml or not
  if (! _jachq)
    _jachq.reset(new SimpleMatrix(sizeY, sizeDS));
  else
  {
    if (_jachq->size(0) == 0) // if the matrix dim are null
      _jachq->resize(sizeY, sizeDS);
    else
      assert((_jachq->size(1) == sizeDS && _jachq->size(0) == sizeY) &&
             "LagrangianR::initComponents inconsistent sizes between Jach[0] matrix and the interaction.");
  }
  // Added by Son Nguyen (8/12/2010)
  if (! _jachqDot)
    _jachqDot.reset(new SimpleMatrix(sizeY, sizeDS));
  else
  {
    if (_jachqDot->size(0) == 0) // if the matrix dimension are null
      _jachqDot->resize(sizeY, sizeDS);
    else
    {
      if ((_jachqDot->size(1) != sizeDS && _jachqDot->size(0) != sizeY))
        RuntimeException::selfThrow("LagrangianR::initComponents inconsistent sizes between Jach[1] matrix and the interaction.");
    }
  }
}

void LagrangianR::initialize(Interaction& inter)
{
  // Memory allocation for G[i], if required (depends on the chosen constructor).
  initComponents(inter);
}

void LagrangianR::computeh(const double time, Interaction& inter)
{
  RuntimeException::selfThrow("LagrangianR::computeh: not yet implemented (or useless) for Lagrangian relation of type " + _subType);
}

void LagrangianR::saveRelationToXML() const
{
  RuntimeException::selfThrow("LagrangianR::saveRelationToXML - not yet implemented.");
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
