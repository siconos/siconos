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
#include "LagrangianRXML.h"
#include "Interaction.h"
#include "LagrangianDS.h"

using namespace std;
using namespace RELATION;

// Default constructor
LagrangianR::LagrangianR():
  Relation(Lagrangian, NonLinearR), LagrangianRelationType(NonLinearR)
{}

// Basic constructor
LagrangianR::LagrangianR(RELATION::SUBTYPES lagType):
  Relation(Lagrangian, lagType), LagrangianRelationType(lagType)
{}

LagrangianR::LagrangianR(SP::RelationXML relxml, RELATION::SUBTYPES newSubType): Relation(relxml, Lagrangian, newSubType), LagrangianRelationType()
{}

void LagrangianR::readGInXML(SP::LagrangianRXML LRxml, unsigned int i)
{
  RELATION::PluginNames name;
  if (i == 0) name = G0;
  else if (i == 1)
    name = G1;
  else //if (i==2)
    name = G2;

  if (LRxml->isGPlugin(i))
  {
    pluginNames[name] = LRxml->getGPlugin(i);
    setComputeGFunction(SSL::getPluginName(pluginNames[name]), SSL::getPluginFunctionName(pluginNames[name]), i);

  }
  else
  {

    G[i].reset(new SimpleMatrix(LRxml->getGMatrix(i)));
    isPlugged[name] = false;
  }
}

LagrangianR::~LagrangianR()
{
  G.clear();
}

void LagrangianR::initialize()
{
  // Check if an Interaction is connected to the Relation.
  if (!interaction)
    RuntimeException::selfThrow("LagrangianR::initialize failed. No Interaction linked to the present relation.");

  // Memory allocation for G[i], if required (depends on the chosen constructor).
  initComponents();

  DSIterator it;
  data["q0"].reset(new BlockVector()); // displacement
  data["q1"].reset(new BlockVector()); // velocity
  data["q2"].reset(new BlockVector()); // acceleration
  data["z"].reset(new BlockVector()); // z vector
  data["p0"].reset(new BlockVector());
  data["p1"].reset(new BlockVector());
  data["p2"].reset(new BlockVector());
  SP::LagrangianDS lds;
  DS::TYPES type;
  for (it = interaction->dynamicalSystemsBegin(); it != interaction->dynamicalSystemsEnd(); ++it)
  {
    type = (*it)->getType();
    // check dynamical system type
    if (type != DS::LLTIDS && type != DS::LNLDS)
      RuntimeException::selfThrow("LagrangianR1::initialize failed, not implemented for dynamical system of type: " + type);

    // convert vDS systems into LagrangianDS and put them in vLDS
    lds = boost::static_pointer_cast<LagrangianDS> (*it);
    // Put q/velocity/acceleration of each DS into a block. (Pointers links, no copy!!)
    data["q0"]->insertPtr(lds->getQPtr());
    data["q1"]->insertPtr(lds->getVelocityPtr());
    data["q2"]->insertPtr(lds->getAccelerationPtr());
    data["p0"]->insertPtr(lds->getPPtr(1));
    data["p1"]->insertPtr(lds->getPPtr(1));
    data["p2"]->insertPtr(lds->getPPtr(2));
    data["z"]->insertPtr(lds->getZPtr());
  }
}

void LagrangianR::initComponents()
{
  unsigned int sizeY = interaction->getSizeOfY();
  unsigned int sizeQ = interaction->getSizeOfDS();

  if (! G[0])
    G[0].reset(new SimpleMatrix(sizeY, sizeQ));

  else // Check if dimension are consistents with interaction.
    if (sizeY != G[0]->size(0) || sizeQ != G[0]->size(1))
      RuntimeException::selfThrow("LagrangianR:: initComponents failed. Inconsistent sizes between Interaction and Relation matrices.");
  workX.reset(new SimpleVector(interaction->getSizeOfDS()));
  workZ.reset(new SimpleVector(interaction->getSizeZ()));
  workY.reset(new SimpleVector(sizeY));
}

void LagrangianR::setGVector(const VectorOfMatrices& newVector)
{
  unsigned int nG = G.size();
  if (newVector.size() != nG)
    RuntimeException::selfThrow("LagrangianR::setGVector(newG) failed. Inconsistent sizes between newG and the problem type. You might have forget to call setLagrangianRelationType before?");


  // If G[i] has been allocated before => delete
  RELATION::PluginNames name;
  G.clear();

  for (unsigned int i = 0; i < nG; i++)
  {
    if (i == 0)
      name = G0;
    else if (i == 1)
      name = G1;
    else
      name = G2;
    G[i] = newVector[i]; // Warning: links to pointers, no copy!!
    isPlugged[name] = false;
  }
}

void LagrangianR::setG(const SiconosMatrix& newValue, unsigned int index)
{
  if (index >= G.size())
    RuntimeException::selfThrow("LagrangianR:: setG(mat,index), index out of range. Maybe you do not call setLagrangianRelationType before?");

  RELATION::PluginNames name;
  if (index == 0)
    name = G0;
  else if (index == 1)
    name = G1;
  else
    name = G2;

  if (! G[index])
    G[index].reset(new SimpleMatrix(newValue));

  else
    *(G[index]) = newValue;

  isPlugged[name] = false;
}

void LagrangianR::setGPtr(SP::SiconosMatrix newPtr, unsigned int  index)
{
  if (index >= G.size())
    RuntimeException::selfThrow("LagrangianR:: setG(mat,index), index out of range. Maybe you do not call setLagrangianRelationType before?");
  RELATION::PluginNames name;
  if (index == 0)
    name = G0;
  else if (index == 1)
    name = G1;
  else
    name = G2;
  G[index] = newPtr;
  isPlugged[name] = false;
}

void LagrangianR::setComputeHFunction(const string& , const string&)
{
  RuntimeException::selfThrow("LagrangianR::setComputeHFunction: not yet implemented (or useless) for Lagrangian relation of type " + subType);
}

void LagrangianR::setComputeGFunction(const string& , const string& , unsigned int)
{
  RuntimeException::selfThrow("LagrangianR::setComputeGFunction: not yet implemented (or useless) for Lagrangian relation of type " + subType);
}

void LagrangianR::computeH(double)
{
  RuntimeException::selfThrow("LagrangianR::computeH: not yet implemented (or useless) for Lagrangian relation of type " + subType);
}

void LagrangianR::computeG(double, unsigned int)
{
  RuntimeException::selfThrow("LagrangianR::computeG: not yet implemented (or useless) for Lagrangian relation of type " + subType);
}

void LagrangianR::saveRelationToXML() const
{
  RuntimeException::selfThrow("LagrangianR1::saveRelationToXML - not yet implemented.");
}

void LagrangianR::display() const
{
  cout << "=====> Lagrangian Relation of type "  << LagrangianRelationType << endl;
  if (interaction) cout << "- Interaction id" << interaction->getId() << endl;
  else cout << "- Linked interaction -> NULL" << endl;
  RELATION::PluginList::const_iterator it;
  cout << "The following operators are linked to plug-in: " << endl;
  for (it = pluginNames.begin(); it != pluginNames.end(); ++it)
    cout << (*it).first << " plugged to:" << (*it).second << endl;
  cout << "===================================== " << endl;
}
