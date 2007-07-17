/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
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
#include "DynamicalSystemsSet.h"
#include "FirstOrderNonLinearDS.h"

using namespace std;

void FirstOrderR::initAllocationFlags(bool in)
{
  isAllocatedIn["jacobianH0"] = in;
  isAllocatedIn["jacobianH1"] = in;
  isAllocatedIn["jacobianG0"] = in;
}

void FirstOrderR::initPluginFlags(bool in)
{
  isPlugged["h"] = in;
  isPlugged["jacobianH0"] = in;
  isPlugged["jacobianH1"] = in;
  isPlugged["g"] = in;
  isPlugged["jacobianG0"] = in ;
}

// Default constructor (protected)
FirstOrderR::FirstOrderR(const string& newType): Relation("FirstOrder", newType), firstOrderType(newType)
{
  initAllocationFlags(false);
  initPluginFlags(false);
  jacobianH.resize(2, NULL);
  jacobianG.resize(1, NULL);
  input.resize(3, NULL);
  output.resize(3, NULL);
}

// xml constructor
FirstOrderR::FirstOrderR(RelationXML* relxml, const string& newType):
  Relation(relxml, "FirstOrder", newType), firstOrderType(newType)
{
  FirstOrderRXML * FORxml = static_cast<FirstOrderRXML*>(relationxml);

  initAllocationFlags(false);
  initPluginFlags(false);
  input.resize(3, NULL);
  output.resize(3, NULL);
  // Gradients
  jacobianH.resize(2, NULL);
  jacobianG.resize(1, NULL);
  // input g
  if (FORxml->hasG())
  {
    pluginNames["g"] = FORxml->getGPlugin();
    setComputeGFunction(cShared.getPluginName(pluginNames["g"]), cShared.getPluginFunctionName(pluginNames["g"]));
    // Gradients
    if (!FORxml->hasJacobianG())
      RuntimeException::selfThrow("FirstOrderR xml constructor failed. No input for gradient(s) of g function.");

    if (FORxml->isJacobianGPlugin())
    {
      string name = pluginNames["jacobianG0"] = FORxml->getJacobianGPlugin();
      setComputeJacobianGFunction(cShared.getPluginName(name), cShared.getPluginFunctionName(name));
    }
    else
    {
      jacobianG[0] = new SimpleMatrix(FORxml->getJacobianGMatrix());
      isAllocatedIn["jacobianG0"] = true   ;
    }
  }

  // output h
  if (FORxml->hasH())
  {
    pluginNames["h"] = FORxml->getHPlugin();
    setComputeHFunction(cShared.getPluginName(pluginNames["h"]), cShared.getPluginFunctionName(pluginNames["h"]));
    // Gradients
    if (!FORxml->hasJacobianH())
      RuntimeException::selfThrow("FirstOrderR xml constructor failed. No input for gradients of h function.");
    string name;
    for (unsigned int i = 0; i < jacobianH.size(); ++i)
    {
      name = "jacobianH" + toString<unsigned int>(i);
      if (FORxml->isJacobianHPlugin(i))
      {
        pluginNames[name] = FORxml->getJacobianHPlugin(i);
        setComputeJacobianHFunction(cShared.getPluginName(pluginNames[name]), cShared.getPluginFunctionName(pluginNames[name]), i);
      }
      else
      {
        jacobianH[i] = new SimpleMatrix(FORxml->getJacobianHMatrix(i));
        isAllocatedIn[name] = true   ;
      }
    }
  }
}

FirstOrderR::FirstOrderR(const string& computeOut, const string& computeIn):
  Relation("FirstOrder", "R"), firstOrderType("R")
{
  // Size vector of pointers to functions.
  input.resize(3, NULL);
  output.resize(3, NULL);
  // Connect input and output to plug-in
  setComputeHFunction(cShared.getPluginName(computeOut), cShared.getPluginFunctionName(computeOut));
  setComputeGFunction(cShared.getPluginName(computeIn), cShared.getPluginFunctionName(computeIn));

  // Nothing for the jacobians => a set is required. Add a constructor?
}

FirstOrderR::~FirstOrderR()
{
  input.clear();
  output.clear();
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
  data["xFree"] = new BlockVector(); // free displacements
  data["z"] = new BlockVector();
  data["r"] = new BlockVector();

  FirstOrderNonLinearDS* ds;
  for (it = interaction->dynamicalSystemsBegin(); it != interaction->dynamicalSystemsEnd(); ++it)
  {
    // Put x/r ... of each DS into a block. (Pointers links, no copy!!)
    ds = static_cast<FirstOrderNonLinearDS*>(*it);
    data["x"]->insertPtr(ds->getXPtr());
    data["xFree"]->insertPtr(ds->getXFreePtr());
    data["z"]->insertPtr(ds->getZPtr());
    data["r"]->insertPtr(ds->getRPtr());
  }
}

void FirstOrderR::initialize()
{
  // Check if an Interaction is connected to the Relation.
  if (interaction == NULL)
    RuntimeException::selfThrow("FirstOrderR::initialize failed. No Interaction linked to the present relation.");

  // Check and allocate (if necessary) jacobianH and G. Can not be done in constructors interaction is not known.
  string name;
  for (unsigned int index = 0; index < jacobianH.size(); ++index)
  {
    name = "jacobianH" + toString<unsigned int>(index);
    if (jacobianH[index] == NULL)
    {
      jacobianH[index] = new SimpleMatrix(interaction->getSizeOfY(), interaction->getSizeOfDS());
      isAllocatedIn[name] = true;
    }
  }

  for (unsigned int index = 0; index < jacobianG.size(); ++index)
  {
    name = "jacobianH" + toString<unsigned int>(index);
    if (jacobianG[index] == NULL)
    {
      jacobianG[index] = new SimpleMatrix(interaction->getSizeOfDS(), interaction->getSizeOfY());
      isAllocatedIn[name] = true;
    }
  }

  // Update data member (links to DS variables)
  initDSLinks();
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
  cShared.setFunction(&(output[0]), pluginPath, functionName);
  string plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames["h"] = plugin + ":" + functionName;
  isPlugged["h"] = true;
}

void FirstOrderR::setComputeJacobianHFunction(const string& pluginPath, const string& functionName, unsigned int index)
{
  string name = "jacobianH" + toString<unsigned int>(index);
  // Warning: output[0] corresponds to h, thus use output[index+1]
  cShared.setFunction(&(output[index + 1]), pluginPath, functionName);
  string plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames[name] = plugin + ":" + functionName;
  isPlugged[name] = true;
}

void FirstOrderR::setComputeGFunction(const string& pluginPath, const string& functionName)
{
  cShared.setFunction(&(input[0]), pluginPath, functionName);
  string plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames["g"] = plugin + ":" + functionName;
  isPlugged["g"] = true;
}

void FirstOrderR::setComputeJacobianGFunction(const string& pluginPath, const string& functionName, unsigned int index)
{
  string name = "jacobianG" + toString<unsigned int>(index);
  cShared.setFunction(&(input[index + 1]), pluginPath, functionName);
  string plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames[name] = plugin + ":" + functionName;
  isPlugged[name] = true;
}

void FirstOrderR::computeOutput(double time, unsigned int)
{
  // Note that the second argument remains unamed since it is not used: for first order systems, we always compute
  // y[0] (at the time).

  if (output[0] == NULL)
    RuntimeException::selfThrow("FirstOrderR::computeOutput() is not linked to a plugin function");

  SiconosVector *y = interaction->getYPtr(0);
  SiconosVector *lambda = interaction->getLambdaPtr(0);
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  SimpleVector * xCopy = new SimpleVector(*data["x"]);
  SimpleVector * zCopy = new SimpleVector(*data["z"]);
  SimpleVector * yCopy = new SimpleVector(*y);
  SimpleVector * lambdaCopy = new SimpleVector(*lambda);

  unsigned int sizeY = y->size();
  unsigned int sizeX = data["x"]->size();
  unsigned int sizeZ = data["z"]->size();

  (output[0])(sizeX, &(*xCopy)(0), time, sizeY, &(*lambdaCopy)(0), &(*yCopy)(0), sizeZ, &(*zCopy)(0));

  // Rebuilt y/z from Tmp
  *y = *yCopy;
  *data["z"] = *zCopy;

  delete zCopy;
  delete lambdaCopy;
  delete yCopy;
  delete xCopy;
}

void FirstOrderR::computeFreeOutput(double time, unsigned int)
{
  // Note that the second argument remains unamed since it is not used: for first order systems, we always compute
  // y[0] (at the time).

  if (output[0] == NULL)
    RuntimeException::selfThrow("FirstOrderR::computeOutput() is not linked to a plugin function");

  SiconosVector *y = interaction->getYPtr(0);
  SiconosVector *lambda = interaction->getLambdaPtr(0);
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  SimpleVector * xCopy = new SimpleVector(*data["xFree"]);
  SimpleVector * zCopy = new SimpleVector(*data["z"]);
  SimpleVector * yCopy = new SimpleVector(*y);
  SimpleVector * lambdaCopy = new SimpleVector(*lambda);

  unsigned int sizeY = y->size();
  unsigned int sizeX = data["x"]->size();
  unsigned int sizeZ = data["z"]->size();

  (output[0])(sizeX, &(*xCopy)(0), time, sizeY, &(*lambdaCopy)(0), &(*yCopy)(0), sizeZ, &(*zCopy)(0));

  // Rebuilt y/z from Tmp
  *y = *yCopy;
  *data["z"] = *zCopy;

  delete zCopy;
  delete lambdaCopy;
  delete yCopy;
  delete xCopy;
}

void FirstOrderR::computeInput(double time, unsigned int level)
{
  if (input[0] == NULL)
    RuntimeException::selfThrow("FirstOrderR::computeInput() is not linked to a plugin function");

  SiconosVector *lambda = interaction->getLambdaPtr(level);
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  SimpleVector * rCopy = new SimpleVector(*data["r"]);
  SimpleVector * zCopy = new SimpleVector(*data["z"]);
  SimpleVector * lambdaCopy = new SimpleVector(*lambda);

  unsigned int sizeY = lambda->size();
  unsigned int sizeZ = data["z"]->size();
  unsigned int sizeR = rCopy->size();

  (input[0])(sizeY, &(*lambdaCopy)(0), time, sizeR, &(*rCopy)(0), sizeZ, &(*zCopy)(0));

  *data["r"] = *rCopy;
  *data["z"] = *zCopy;

  delete rCopy;
  delete zCopy;
  delete lambdaCopy;
}

void FirstOrderR::computeJacobianH(double time, unsigned int i)
{
  if (output[i + 1] == NULL)
    RuntimeException::selfThrow("FirstOrderR::computeJacobianH() failed; not linked to a plug-in function.");

  SiconosVector *lambda = interaction->getLambdaPtr(0);
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  SimpleVector * xCopy = new SimpleVector(*data["x"]);
  SimpleVector * zCopy = new SimpleVector(*data["z"]);
  SimpleVector * lambdaCopy = new SimpleVector(*lambda);

  unsigned int sizeY = lambda->size();
  unsigned int sizeX = data["x"]->size();
  unsigned int sizeZ = data["z"]->size();

  (output[i + 1])(sizeX, &(*xCopy)(0), time, sizeY, &(*lambdaCopy)(0), &(*(jacobianH[i]))(0, 0), sizeZ, &(*zCopy)(0));

  // Rebuilt z from Tmp
  *data["z"] = *zCopy;

  delete zCopy;
  delete lambdaCopy;
  delete xCopy;
}

void FirstOrderR::computeJacobianG(double time, unsigned int i)
{
  // At the time, second parameter is not use: only one possible jacobian => i = 0 is the default value.
  if (input[i + 1] == NULL)
    RuntimeException::selfThrow("FirstOrderR::computeJacobianH() failed; not linked to a plug-in function.");

  SiconosVector *lambda = interaction->getLambdaPtr(0);
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  SimpleVector * zCopy = new SimpleVector(*data["z"]);
  SimpleVector * lambdaCopy = new SimpleVector(*lambda);

  unsigned int sizeY = lambda->size();
  unsigned int sizeX = data["x"]->size();
  unsigned int sizeZ = data["z"]->size();

  (input[i + 1])(sizeY, &(*lambdaCopy)(0), time, sizeX, &(*(jacobianG[i]))(0, 0), sizeZ, &(*zCopy)(0));

  // Rebuilt z from Tmp
  *data["z"] = *zCopy;

  delete zCopy;
  delete lambdaCopy;
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

FirstOrderR* FirstOrderR::convert(Relation *r)
{
  return dynamic_cast<FirstOrderR*>(r);
}

