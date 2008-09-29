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

// Default constructor
LagrangianR::LagrangianR():
  Relation(Lagrangian, NonLinearR), LagrangianRelationType(NonLinearR)
{}

// Basic constructor
LagrangianR::LagrangianR(RELATIONSUBTYPES lagType):
  Relation(Lagrangian, lagType), LagrangianRelationType(lagType)
{}

LagrangianR::LagrangianR(RelationXMLSPtr relxml, RELATIONSUBTYPES newSubType): Relation(relxml, Lagrangian, newSubType), LagrangianRelationType()
{}

void LagrangianR::readGInXML(LagrangianRXMLSPtr LRxml, unsigned int i)
{
  string name = "G" + toString<unsigned int>(i);
  if (LRxml->isGPlugin(i))
  {
    pluginNames[name] = LRxml->getGPlugin(i);
    setComputeGFunction(cShared.getPluginName(pluginNames[name]), cShared.getPluginFunctionName(pluginNames[name]), i);

#ifndef WithSmartPtr
    isAllocatedIn[name] = false;
#endif
  }
  else
  {

#ifndef WithSmartPtr
    G[i] = new SimpleMatrix(LRxml->getGMatrix(i));
    isAllocatedIn[name] = true   ;
#else
    G[i].reset(new SimpleMatrix(LRxml->getGMatrix(i)));
#endif

    isPlugged[name] = false;
  }
}

LagrangianR::~LagrangianR()
{

#ifndef WithSmartPtr
  string name;
  for (unsigned int i = 0; i < G.size(); ++i)
  {
    name = "G" + toString<unsigned int>(i);
    if (isAllocatedIn[name]) delete G[i];
  }
#endif

  G.clear();
}

void LagrangianR::initialize()
{
  // Check if an Interaction is connected to the Relation.
  if (interaction == NULL)
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
  LagrangianDSSPtr lds;
  DSTYPES type;
  for (it = interaction->dynamicalSystemsBegin(); it != interaction->dynamicalSystemsEnd(); ++it)
  {
    type = (*it)->getType();
    // check dynamical system type
    if (type != LLTIDS && type != LNLDS)
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
  {

#ifndef WithSmartPtr
    G[0] = new SimpleMatrix(sizeY, sizeQ);
    isAllocatedIn["G0"] = true;
#else
    G[0].reset(new SimpleMatrix(sizeY, sizeQ));
#endif

  }
  else // Check if dimension are consistents with interaction.
    if (sizeY != G[0]->size(0) || sizeQ != G[0]->size(1))
      RuntimeException::selfThrow("LagrangianR:: initComponents failed. Inconsistent sizes between Interaction and Relation matrices.");

#ifndef WithSmartPtr
  workX = new SimpleVector(interaction->getSizeOfDS());
  workZ = new SimpleVector(interaction->getSizeZ());
  workY = new SimpleVector(sizeY);
#else
  workX.reset(new SimpleVector(interaction->getSizeOfDS()));
  workZ.reset(new SimpleVector(interaction->getSizeZ()));
  workY.reset(new SimpleVector(sizeY));
#endif
}

void LagrangianR::setGVector(const VectorOfMatrices& newVector)
{
  unsigned int nG = G.size();
  if (newVector.size() != nG)
    RuntimeException::selfThrow("LagrangianR::setGVector(newG) failed. Inconsistent sizes between newG and the problem type. You might have forget to call setLagrangianRelationType before?");


  // If G[i] has been allocated before => delete
  string name;

#ifndef WithSmartPtr
  for (unsigned int i = 0; i < nG; i++)
  {
    name = "G" + toString<unsigned int>(i);
    if (isAllocatedIn[name]) delete G[i];
    G[i] = NULL;
    isAllocatedIn[name] = false;
    isPlugged[name] = false;
  }
#endif

  G.clear();

  for (unsigned int i = 0; i < nG; i++)
  {
    G[i] = newVector[i]; // Warning: links to pointers, no copy!!
    name = "G" + toString<unsigned int>(i);

#ifndef WithSmartPtr
    isAllocatedIn[name] = false;
#endif

    isPlugged[name] = false;
  }
}

void LagrangianR::setG(const SiconosMatrix& newValue, unsigned int index)
{
  if (index >= G.size())
    RuntimeException::selfThrow("LagrangianR:: setG(mat,index), index out of range. Maybe you do not call setLagrangianRelationType before?");

  string name = "G" + toString<unsigned int>(index);
  if (! G[index])
  {

#ifndef WithSmartPtr
    G[index] =  new SimpleMatrix(newValue);
    isAllocatedIn[name] = true;
#else
    G[index].reset(new SimpleMatrix(newValue));
#endif

  }
  else
    *(G[index]) = newValue;

  isPlugged[name] = false;
}

void LagrangianR::setGPtr(SiconosMatrixSPtr newPtr, unsigned int  index)
{
  if (index >= G.size())
    RuntimeException::selfThrow("LagrangianR:: setG(mat,index), index out of range. Maybe you do not call setLagrangianRelationType before?");

  string name = "G" + toString<unsigned int>(index);

#ifndef WithSmartPtr
  if (isAllocatedIn[name]) delete G[index];
  isAllocatedIn[name] = false;
#endif

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
  if (interaction != NULL) cout << "- Interaction id" << interaction->getId() << endl;
  else cout << "- Linked interaction -> NULL" << endl;
  NamesConstIterator it;
  cout << "The following operators are linked to plug-in: " << endl;
  for (it = pluginNames.begin(); it != pluginNames.end(); ++it)
    cout << (*it).first << " plugged to:" << (*it).second << endl;
  cout << "===================================== " << endl;
}
