/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */

// \todo : create a work vector for all tmp vectors used in computeG, computeH ...

#include "LagrangianR.h"
#include "RelationXML.h"
#include "Interaction.h"
#include "LagrangianDS.h"

using namespace std;

void LagrangianR::initComponents()
{
  unsigned int sizeY = getInteractionPtr()->getSizeOfY();
  unsigned int sizeDS = getInteractionPtr()->getSizeOfDS();

  // The initialization of JacH[0] depends on the way the Relation was built ie if the matrix
  // was read from xml or not
  if (! JacQH)
    JacQH.reset(new SimpleMatrix(sizeY, sizeDS));
  else
  {
    if (JacQH->size(0) == 0) // if the matrix dim are null
      JacQH->resize(sizeY, sizeDS);
    else
      assert((JacQH->size(1) == sizeDS && JacQH->size(0) == sizeY) &&
             "LagrangianScleronomousR::initComponents inconsistent sizes between JacH[0] matrix and the interaction.");
  }

  workX.reset(new SimpleVector(sizeDS));
  workZ.reset(new SimpleVector(getInteractionPtr()->getSizeZ()));
  workY.reset(new SimpleVector(sizeY));
}

void LagrangianR::initialize(SP::Interaction inter)
{
  assert(inter && "FirstOrderR::initialize failed. No Interaction linked to the present relation.");
  interaction = inter;

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
  for (it = getInteractionPtr()->dynamicalSystemsBegin(); it != getInteractionPtr()->dynamicalSystemsEnd(); ++it)
  {
    type = (*it)->getType();
    // check dynamical system type
    assert((type == DS::LLTIDS || type == DS::LNLDS) && "LagrangianR::initialize failed, not implemented for dynamical system of type: " + type);

    // convert vDS systems into LagrangianDS and put them in vLDS
    lds = boost::static_pointer_cast<LagrangianDS> (*it);
    // Put q/velocity/acceleration of each DS into a block. (Pointers links, no copy!!)
    data[q0]->insertPtr(lds->getQPtr());
    data[q1]->insertPtr(lds->getVelocityPtr());
    data[q2]->insertPtr(lds->getAccelerationPtr());
    data[p0]->insertPtr(lds->getPPtr(1));
    data[p1]->insertPtr(lds->getPPtr(1));
    data[p2]->insertPtr(lds->getPPtr(2));
    data[z]->insertPtr(lds->getZPtr());
  }
}


void LagrangianR::computeH(double)
{
  RuntimeException::selfThrow("LagrangianR::computeH: not yet implemented (or useless) for Lagrangian relation of type " + subType);
}

//  void LagrangianR::computeJacXH(double)
// {
//   RuntimeException::selfThrow("FirstOrderR::computeJacobianXH, not (yet) implemented or forbidden for relations of type "+subType);
// }
//  void LagrangianR::computeJacLH(double)
// {
//   RuntimeException::selfThrow("FirstOrderR::computeJacobianLH, not (yet) implemented or forbidden for relations of type "+subType);
// }

void LagrangianR::saveRelationToXML() const
{
  RuntimeException::selfThrow("LagrangianR1::saveRelationToXML - not yet implemented.");
}

void LagrangianR::display() const
{
  Relation::display();
}
