/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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

// void FirstOrderR::initAllocationFlags(bool in)
// {
// }

// void FirstOrderR::initPluginFlags(bool in)
// {
// }

// Default constructor
FirstOrderR::FirstOrderR(const string& newType): Relation("FirstOrder", newType), firstOrderType(newType)
{}

// xml constructor
FirstOrderR::FirstOrderR(RelationXML* relxml, const string& newType):
  Relation(relxml, "FirstOrder", newType), firstOrderType(newType)
{}

FirstOrderR::~FirstOrderR()
{
  string name;
  for (unsigned int i = 0; i < jacobianH.size(); ++i)
  {
    name = "jacobianH" + toString<unsigned int>(i);
    if (isAllocatedIn[name]) delete jacobianH[i];
  }
  jacobianH.clear();
  for (unsigned int i = 0; i < jacobianG.size(); ++i)
  {
    name = "jacobianG" + toString<unsigned int>(i);
    if (isAllocatedIn[name]) delete jacobianG[i];
  }
  jacobianG.clear();
}

void FirstOrderR::initDSLinks()
{
  // Get the DS concerned by the interaction of this relation
  DSIterator it;
  data["x"] = new BlockVector(); // displacements
  data["z"] = new BlockVector();
  data["r"] = new BlockVector();

  FirstOrderNonLinearDS* ds;
  for (it = interaction->dynamicalSystemsBegin(); it != interaction->dynamicalSystemsEnd(); ++it)
  {
    // Put x/r ... of each DS into a block. (Pointers links, no copy!!)
    ds = static_cast<FirstOrderNonLinearDS*>(*it);
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
  if (interaction == NULL)
    RuntimeException::selfThrow("FirstOrderR::initialize failed. No Interaction linked to the present relation.");

  // Update data member (links to DS variables)
  initDSLinks();
  // Initialize work vectors
  workX = new SimpleVector(sizeX);
  workZ = new SimpleVector(sizeZ);
  workY = new SimpleVector(sizeY);
}

void FirstOrderR::setJacobianHVector(const VectorOfMatrices& newVector)
{
  unsigned int nJH = jacobianH.size();
  if (newVector.size() != nJH)
    RuntimeException::selfThrow("FirstOrderR::setJacobianHVector(newH) failed. Inconsistent sizes between newH and the problem type.");

  // If jacobianH[i] has been allocated before => delete
  string name;
  for (unsigned int i = 0; i < nJH; i++)
  {
    name = "jacobianH" + toString<unsigned int>(i);
    if (isAllocatedIn[name]) delete jacobianH[i];
    jacobianH[i] = NULL;
    isAllocatedIn[name] = false;
    isPlugged[name] = false;
  }

  jacobianH.clear();

  for (unsigned int i = 0; i < nJH; i++)
  {
    jacobianH[i] = newVector[i]; // Warning: links to pointers, no copy!!
    name = "jacobianH" + toString<unsigned int>(i);
    isAllocatedIn[name] = false;
    isPlugged[name] = false;
  }
}

void FirstOrderR::setJacobianH(const SiconosMatrix& newValue, unsigned int index)
{
  if (index >= jacobianH.size())
    RuntimeException::selfThrow("FirstOrderR:: setJacobianH(mat,index), index out of range.");

  string name = "jacobianH" + toString<unsigned int>(index);
  if (jacobianH[index] == NULL)
  {
    jacobianH[index] =  new SimpleMatrix(newValue);
    isAllocatedIn[name] = true;
  }
  else
    *(jacobianH[index]) = newValue;

  isPlugged[name] = false;
}

void FirstOrderR::setJacobianHPtr(SiconosMatrix *newPtr, unsigned int  index)
{
  if (index >= jacobianH.size())
    RuntimeException::selfThrow("FirstOrderR:: setJacobianH(mat,index), index out of range.");

  string name = "jacobianH" + toString<unsigned int>(index);
  if (isAllocatedIn[name]) delete jacobianH[index];
  jacobianH[index] = newPtr;
  isAllocatedIn[name] = false;
  isPlugged[name] = false;
}

void FirstOrderR::setJacobianGVector(const VectorOfMatrices& newVector)
{
  unsigned int nJG = jacobianG.size();
  if (newVector.size() != nJG)
    RuntimeException::selfThrow("FirstOrderR::setJacobianGVector(newG) failed. Inconsistent sizes between newG and the problem type.");

  // If jacobianG[i] has been allocated before => delete
  string name;
  for (unsigned int i = 0; i < nJG; i++)
  {
    name = "jacobianG" + toString<unsigned int>(i);
    if (isAllocatedIn[name]) delete jacobianG[i];
    jacobianG[i] = NULL;
    isAllocatedIn[name] = false;
    isPlugged[name] = false;
  }

  jacobianG.clear();

  for (unsigned int i = 0; i < nJG; i++)
  {
    jacobianG[i] = newVector[i]; // Warning: links to pointers, no copy!!
    name = "jacobianG" + toString<unsigned int>(i);
    isAllocatedIn[name] = false;
    isPlugged[name] = false;
  }
}

void FirstOrderR::setJacobianG(const SiconosMatrix& newValue, unsigned int index)
{
  if (index >= jacobianG.size())
    RuntimeException::selfThrow("FirstOrderR:: setJacobianG(mat,index), index out of range.");

  string name = "jacobianG" + toString<unsigned int>(index);
  if (jacobianG[index] == NULL)
  {
    jacobianG[index] =  new SimpleMatrix(newValue);
    isAllocatedIn[name] = true;
  }
  else
    *(jacobianG[index]) = newValue;

  isPlugged[name] = false;
}

void FirstOrderR::setJacobianGPtr(SiconosMatrix *newPtr, unsigned int  index)
{
  if (index >= jacobianG.size())
    RuntimeException::selfThrow("FirstOrderR:: setJacobianG(mat,index), index out of range.");

  string name = "jacobianG" + toString<unsigned int>(index);
  if (isAllocatedIn[name]) delete jacobianG[index];
  jacobianG[index] = newPtr;
  isAllocatedIn[name] = false;
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
  //   cShared.setFunction(&(output[index+1]), pluginPath, functionName);
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
  //   cShared.setFunction(&(input[index+1]), pluginPath, functionName);
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
  if (interaction != NULL) cout << "- Interaction id" << interaction->getId() << endl;
  else cout << "- Linked interaction -> NULL" << endl;
  NamesConstIterator it;
  cout << "The following operators are linked to plug-in: " << endl;
  for (it = pluginNames.begin(); it != pluginNames.end(); ++it)
    cout << (*it).first << " plugged to:" << (*it).second << endl;
  cout << "===================================== " << endl;
}

