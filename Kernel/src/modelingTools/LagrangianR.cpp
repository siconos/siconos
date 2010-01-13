/* Siconos-Kernel, Copyright INRIA 2005-2010.
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

void LagrangianR::initComponents()
{
  unsigned int sizeY = interaction()->getSizeOfY();
  unsigned int sizeDS = interaction()->getSizeOfDS();

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
             "LagrangianScleronomousR::initComponents inconsistent sizes between Jach[0] matrix and the interaction.");
  }

  _workX.reset(new SimpleVector(sizeDS));
  _workZ.reset(new SimpleVector(interaction()->getSizeZ()));
  _workY.reset(new SimpleVector(sizeY));
}

void LagrangianR::initialize(SP::Interaction inter)
{
  assert(inter && "FirstOrderR::initialize failed. No Interaction linked to the present relation.");
  _interaction = inter;

  // Memory allocation for G[i], if required (depends on the chosen constructor).
  initComponents();
  data.resize(sizeDataNames);

  DSIterator it;
  data[q0].reset(new BlockVector()); // displacement
  data[q1].reset(new BlockVector()); // velocity
  data[q2].reset(new BlockVector()); // acceleration
  data[z].reset(new BlockVector()); // z vector
  data[p0].reset(new BlockVector());
  data[p1].reset(new BlockVector());
  data[p2].reset(new BlockVector());
  SP::LagrangianDS lds;
  DS::TYPES type;
  for (it = interaction()->dynamicalSystemsBegin(); it != interaction()->dynamicalSystemsEnd(); ++it)
  {
    type = (*it)->getType();
    // check dynamical system type
    assert((type == DS::LLTIDS || type == DS::LNLDS) && "LagrangianR::initialize failed, not implemented for dynamical system of type: " + type);

    // convert vDS systems into LagrangianDS and put them in vLDS
    lds = boost::static_pointer_cast<LagrangianDS> (*it);
    // Put q/velocity/acceleration of each DS into a block. (Pointers links, no copy!!)
    data[q0]->insertPtr(lds->q());
    data[q1]->insertPtr(lds->velocity());
    data[q2]->insertPtr(lds->acceleration());
    data[p0]->insertPtr(lds->p(1));
    data[p1]->insertPtr(lds->p(1));
    data[p2]->insertPtr(lds->p(2));
    data[z]->insertPtr(lds->z());
  }
}


void LagrangianR::computeh(double)
{
  RuntimeException::selfThrow("LagrangianR::computeh: not yet implemented (or useless) for Lagrangian relation of type " + subType);
}

void LagrangianR::saveRelationToXML() const
{
  RuntimeException::selfThrow("LagrangianR1::saveRelationToXML - not yet implemented.");
}

void LagrangianR::display() const
{
  Relation::display();
}
