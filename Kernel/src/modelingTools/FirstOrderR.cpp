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
#include "FirstOrderR.h"
#include "FirstOrderRXML.h"
#include "Interaction.h"
#include "FirstOrderNonLinearDS.h"

using namespace std;
using namespace RELATION;

// Default constructor
FirstOrderR::FirstOrderR(): Relation(FirstOrder, NonLinearR), firstOrderType(NonLinearR)
{}
// Basic constructor
FirstOrderR::FirstOrderR(RELATION::SUBTYPES newType): Relation(FirstOrder, newType), firstOrderType(newType)
{}

// xml constructor
FirstOrderR::FirstOrderR(SP::RelationXML relxml, RELATION::SUBTYPES newType):
  Relation(relxml, FirstOrder, newType), firstOrderType(newType)
{}

FirstOrderR::~FirstOrderR()
{
  jacobianH.clear();
  jacobianG.clear();
}

void FirstOrderR::initDSLinks()
{
  // Get the DS concerned by the interaction of this relation
  DSIterator it;
  data["x"].reset(new BlockVector()); // displacements
  data["z"].reset(new BlockVector());
  data["r"].reset(new BlockVector());

  SP::FirstOrderNonLinearDS ds;
  for (it = interaction->dynamicalSystemsBegin(); it != interaction->dynamicalSystemsEnd(); ++it)
  {
    // Put x/r ... of each DS into a block. (Pointers links, no copy!!)
    ds = boost::static_pointer_cast<FirstOrderNonLinearDS> (*it);
    data["x"]->insertPtr(ds->getXPtr());
    data["z"]->insertPtr(ds->getZPtr());
    data["r"]->insertPtr(ds->getRPtr());
  }
}

void FirstOrderR::initialize()
{
  // Check if an Interaction is connected to the Relation.
  unsigned int sizeY = interaction->getSizeOfY();
  unsigned int sizeX = interaction->getSizeOfDS();
  unsigned int sizeZ = interaction->getSizeZ();
  if (!interaction)
    RuntimeException::selfThrow("FirstOrderR::initialize failed. No Interaction linked to the present relation.");

  // Update data member (links to DS variables)
  initDSLinks();
  // Initialize work vectors

  workX.reset(new SimpleVector(sizeX));
  workZ.reset(new SimpleVector(sizeZ));
  workY.reset(new SimpleVector(sizeY));
}

void FirstOrderR::setJacobianHVector(const VectorOfMatrices& newVector)
{
  unsigned int nJH = jacobianH.size();
  if (newVector.size() != nJH)
    RuntimeException::selfThrow("FirstOrderR::setJacobianHVector(newH) failed. Inconsistent sizes between newH and the problem type.");

  RELATION::PluginNames name;

  jacobianH.clear();

  for (unsigned int i = 0; i < nJH; i++)
  {
    if (i == 0)
      name = RELATION::jacobianH0;
    else if (i == 1)
      name = RELATION::jacobianH1;
    else
      name = RELATION::jacobianH2;
    jacobianH[i] = newVector[i]; // Warning: links to pointers, no copy!!
    isPlugged[name] = false;
  }
}

void FirstOrderR::setJacobianH(const SiconosMatrix& newValue, unsigned int index)
{
  if (index >= jacobianH.size())
    RuntimeException::selfThrow("FirstOrderR:: setJacobianH(mat,index), index out of range.");

  RELATION::PluginNames name;
  if (index == 0)
    name = RELATION::jacobianH0;
  else if (index == 1)
    name = RELATION::jacobianH1;
  else
    name = RELATION::jacobianH2;
  if (! jacobianH[index])
    jacobianH[index].reset(new SimpleMatrix(newValue));

  else
    *(jacobianH[index]) = newValue;

  isPlugged[name] = false;
}

void FirstOrderR::setJacobianHPtr(SP::SiconosMatrix newPtr, unsigned int  index)
{
  if (index >= jacobianH.size())
    RuntimeException::selfThrow("FirstOrderR:: setJacobianH(mat,index), index out of range.");

  RELATION::PluginNames name;
  if (index == 0)
    name = RELATION::jacobianH0;
  else if (index == 1)
    name = RELATION::jacobianH1;
  else
    name = RELATION::jacobianH2;
  jacobianH[index] = newPtr;
  isPlugged[name] = false;
}

void FirstOrderR::setJacobianGVector(const VectorOfMatrices& newVector)
{
  unsigned int nJG = jacobianG.size();
  if (newVector.size() != nJG)
    RuntimeException::selfThrow("FirstOrderR::setJacobianGVector(newG) failed. Inconsistent sizes between newG and the problem type.");
  RELATION::PluginNames name;
  jacobianG.clear();

  for (unsigned int i = 0; i < nJG; i++)
  {
    if (i == 0)
      name = RELATION::jacobianG0;
    else if (i == 1)
      name = RELATION::jacobianG1;
    else
      name = RELATION::jacobianG2;
    jacobianG[i] = newVector[i]; // Warning: links to pointers, no copy!!
    isPlugged[name] = false;
  }
}

void FirstOrderR::setJacobianG(const SiconosMatrix& newValue, unsigned int index)
{
  if (index >= jacobianG.size())
    RuntimeException::selfThrow("FirstOrderR:: setJacobianG(mat,index), index out of range.");

  RELATION::PluginNames name;
  if (index == 0)
    name = RELATION::jacobianG0;
  else if (index == 1)
    name = RELATION::jacobianG1;
  else
    name = RELATION::jacobianG2;
  if (! jacobianG[index])
    jacobianG[index].reset(new SimpleMatrix(newValue));

  else
    *(jacobianG[index]) = newValue;

  isPlugged[name] = false;
}

void FirstOrderR::setJacobianGPtr(SP::SiconosMatrix newPtr, unsigned int  index)
{
  if (index >= jacobianG.size())
    RuntimeException::selfThrow("FirstOrderR:: setJacobianG(mat,index), index out of range.");

  RELATION::PluginNames name;
  if (index == 0)
    name = RELATION::jacobianG0;
  else if (index == 1)
    name = RELATION::jacobianG1;
  else
    name = RELATION::jacobianG2;
  jacobianG[index] = newPtr;
  isPlugged[name] = false;
}

void FirstOrderR::setComputeHFunction(const string& pluginPath, const string& functionName)
{
  RuntimeException::selfThrow("FirstOrderR::setComputeHFunction, not (yet) implemented or forbidden for relations of type " + firstOrderType);
}

void FirstOrderR::setComputeJacobianHFunction(const string& pluginPath, const string& functionName, unsigned int index)
{
  RuntimeException::selfThrow("FirstOrderR::setComputeJacobianHFunction, not (yet) implemented or forbidden for relations of type " + firstOrderType);
  //   string name = "jacobianH"+toString<unsigned int>(index);
  //   // Warning: output[0] corresponds to h, thus use output[index+1]
  //   SSL::setFunction(&(output[index+1]), pluginPath, functionName);
  //   string plugin = pluginPath.substr(0, pluginPath.length()-3);
  //   pluginNames[name] = plugin + ":" + functionName;
  //   isPlugged[name] = true;
}

void FirstOrderR::setComputeGFunction(const string& pluginPath, const string& functionName)
{
  RuntimeException::selfThrow("FirstOrderR::setComputeGFunction, not (yet) implemented or forbidden for relations of type " + firstOrderType);
}

void FirstOrderR::setComputeJacobianGFunction(const string& pluginPath, const string& functionName, unsigned int index)
{
  RuntimeException::selfThrow("FirstOrderR::setComputeJacobianGFunction, not (yet) implemented or forbidden for relations of type " + firstOrderType);
  //   string name = "jacobianG"+toString<unsigned int>(index);
  //   SSL::setFunction(&(input[index+1]), pluginPath, functionName);
  //   string plugin = pluginPath.substr(0, pluginPath.length()-3);
  //   pluginNames[name] = plugin + ":" + functionName;
  //   isPlugged[name] = true;
}

void FirstOrderR::computeJacobianH(double, unsigned int)
{
  RuntimeException::selfThrow("FirstOrderR::computeJacobianH, not (yet) implemented or forbidden for relations of type " + firstOrderType);
}

void FirstOrderR::computeJacobianG(double, unsigned int)
{
  RuntimeException::selfThrow("FirstOrderR::computeJacobianG, not (yet) implemented or forbidden for relations of type " + firstOrderType);
}

void FirstOrderR::display() const
{
  cout << "=====> FirstOrderR of type " << firstOrderType << endl;
  if (interaction) cout << "- Interaction id" << interaction->getId() << endl;
  else cout << "- Linked interaction -> NULL" << endl;
  PluginList::const_iterator it;
  cout << "The following operators are linked to plug-in: " << endl;
  for (it = pluginNames.begin(); it != pluginNames.end(); ++it)
    cout << (*it).first << " plugged to:" << (*it).second << endl;
  cout << "===================================== " << endl;
}

