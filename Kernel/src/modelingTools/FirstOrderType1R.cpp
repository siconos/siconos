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
#include "FirstOrderType1R.h"
#include "FirstOrderRXML.h"
#include "Interaction.h"
#include "FirstOrderNonLinearDS.h"

using namespace std;

void FirstOrderType1R::initAllocationFlags(bool in)
{
  isAllocatedIn["jacobianH0"] = in;
  isAllocatedIn["jacobianG0"] = in;
}

void FirstOrderType1R::initPluginFlags(bool in)
{
  isPlugged["h"] = in;
  isPlugged["jacobianH0"] = in;
  isPlugged["g"] = in;
  isPlugged["jacobianG0"] = in ;
}

// xml constructor
FirstOrderType1R::FirstOrderType1R(RelationXML* relxml):
  FirstOrderR(relxml, Type1R), output(NULL), jXOutput(NULL), input(NULL), jLInput(NULL)
{
  FirstOrderRXML * FORxml = static_cast<FirstOrderRXML*>(relationxml);

  initAllocationFlags(false);
  initPluginFlags(false);
  // Gradients
  jacobianH.resize(1, NULL);
  jacobianG.resize(1, NULL);
  // input g
  if (FORxml->hasG())
  {
    pluginNames["g"] = FORxml->getGPlugin();
    setComputeGFunction(cShared.getPluginName(pluginNames["g"]), cShared.getPluginFunctionName(pluginNames["g"]));
    // Gradients
    if (!FORxml->hasJacobianG())
      RuntimeException::selfThrow("FirstOrderType1R xml constructor failed. No input for gradient(s) of g function.");

    if (FORxml->isJacobianGPlugin(0))
    {
      pluginNames["jacobianG0"] = FORxml->getJacobianGPlugin(0);
      setComputeJacobianGFunction(cShared.getPluginName(pluginNames["jacobianG0"]), cShared.getPluginFunctionName(pluginNames["jacobianG0"]));
    }
    else
    {
      jacobianG[0] = new SimpleMatrix(FORxml->getJacobianGMatrix(0));
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
      RuntimeException::selfThrow("FirstOrderType1R xml constructor failed. No input for gradients of h function.");
    if (FORxml->isJacobianHPlugin(0))
    {
      pluginNames["jacobianH0"] = FORxml->getJacobianHPlugin(0);
      setComputeJacobianHFunction(cShared.getPluginName(pluginNames["jacobianH0"]), cShared.getPluginFunctionName(pluginNames["jacobianH0"]));
    }
    else
    {
      jacobianH[0] = new SimpleMatrix(FORxml->getJacobianHMatrix(0));
      isAllocatedIn["jacobianH0"] = true   ;
    }
  }
}

FirstOrderType1R::FirstOrderType1R(const string& computeOut, const string& computeIn):
  FirstOrderR(Type1R), output(NULL), jXOutput(NULL), input(NULL), jLInput(NULL)
{
  initAllocationFlags(false);
  initPluginFlags(false);
  // Size vector of pointers to functions.
  // Connect input and output to plug-in
  setComputeHFunction(cShared.getPluginName(computeOut), cShared.getPluginFunctionName(computeOut));
  setComputeGFunction(cShared.getPluginName(computeIn), cShared.getPluginFunctionName(computeIn));

  jacobianH.resize(1, NULL);
  jacobianG.resize(1, NULL);
  // The jacobians are not set, and thus considered as null matrices at this point.
}

FirstOrderType1R::FirstOrderType1R(const string& computeOut, const string& computeIn, const string& computeJX, const string& computeJL):
  FirstOrderR(Type1R), output(NULL), jXOutput(NULL), input(NULL), jLInput(NULL)
{
  initAllocationFlags(false);
  initPluginFlags(false);
  // Size vector of pointers to functions.
  jacobianH.resize(1, NULL);
  jacobianG.resize(1, NULL);
  // Connect input and output to plug-in
  setComputeHFunction(cShared.getPluginName(computeOut), cShared.getPluginFunctionName(computeOut));
  setComputeGFunction(cShared.getPluginName(computeIn), cShared.getPluginFunctionName(computeIn));
  setComputeJacobianHFunction(cShared.getPluginName(computeJX), cShared.getPluginFunctionName(computeJX));
  setComputeJacobianGFunction(cShared.getPluginName(computeJL), cShared.getPluginFunctionName(computeJL));
}

FirstOrderType1R::~FirstOrderType1R()
{
  input = NULL;
  output = NULL;
  jLInput = NULL;
  jXOutput = NULL;
}

void FirstOrderType1R::initialize()
{
  FirstOrderR::initialize();

  // if jacobianH0 is plugged and memory not allocated ...
  if (isPlugged["jacobianH0"] && jacobianH[0] == NULL)
  {
    jacobianH[0] = new SimpleMatrix(interaction->getSizeOfY(), interaction->getSizeOfDS());
    isAllocatedIn["jacobianH0"] = true;
  }

  // Same work for jacobianLambdaG
  if (isPlugged["jacobianG0"] && jacobianG[0] == NULL)
  {
    jacobianG[0] = new SimpleMatrix(interaction->getSizeOfDS(), interaction->getSizeOfY());
    isAllocatedIn["jacobianG0"] = true;
  }
}

void FirstOrderType1R::setComputeHFunction(const string& pluginPath, const string& functionName)
{
  cShared.setFunction(&output, pluginPath, functionName);
  string plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames["h"] = plugin + ":" + functionName;
  isPlugged["h"] = true;
}

void FirstOrderType1R::setComputeJacobianHFunction(const string& pluginPath, const string& functionName, unsigned int)
{
  string name = "jacobianH0";
  // Warning: output[0] corresponds to h, thus use output[index+1]
  cShared.setFunction(&jXOutput, pluginPath, functionName);
  string plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames[name] = plugin + ":" + functionName;
  isPlugged[name] = true;
}

void FirstOrderType1R::setComputeGFunction(const string& pluginPath, const string& functionName)
{
  cShared.setFunction(&input, pluginPath, functionName);
  string plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames["g"] = plugin + ":" + functionName;
  isPlugged["g"] = true;
}

void FirstOrderType1R::setComputeJacobianGFunction(const string& pluginPath, const string& functionName, unsigned int)
{
  string name = "jacobianG0";
  cShared.setFunction(&jLInput, pluginPath, functionName);
  string plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames[name] = plugin + ":" + functionName;
  isPlugged[name] = true;
}

void FirstOrderType1R::computeOutput(double, unsigned int)
{
  if (output == NULL)
    RuntimeException::selfThrow("FirstOrderType1R::computeOutput() is not linked to a plugin function");

  SiconosVector *y = interaction->getYPtr(0);
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.

  *workX = *data["x"];
  *workZ = *data["z"];
  *workY = *y;

  unsigned int sizeY = y->size();
  unsigned int sizeX = data["x"]->size();
  unsigned int sizeZ = data["z"]->size();

  output(sizeX, &(*workX)(0), sizeY, &(*workY)(0), sizeZ, &(*workZ)(0));

  // Rebuilt y/z from Tmp
  *y = *workY;
  *data["z"] = *workZ;
}

void FirstOrderType1R::computeInput(double, unsigned int level)
{
  if (input == NULL)
    RuntimeException::selfThrow("FirstOrderType1R::computeInput() is not linked to a plugin function");

  SiconosVector *lambda = interaction->getLambdaPtr(level);
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.

  *workX = *data["r"];
  *workZ = *data["z"];
  *workY = *lambda;

  unsigned int sizeY = lambda->size();
  unsigned int sizeZ = data["z"]->size();
  unsigned int sizeR = workX->size();

  input(sizeY, &(*workY)(0), sizeR, &(*workX)(0), sizeZ, &(*workZ)(0));

  *data["r"] = *workX;
  *data["z"] = *workZ;
}

void FirstOrderType1R::computeJacobianH(double, unsigned int)
{
  if (jXOutput == NULL)
    RuntimeException::selfThrow("FirstOrderType1R::computeJacobianH() failed; not linked to a plug-in function.");

  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  *workX = *data["x"];
  *workZ = *data["z"];

  unsigned int sizeY = interaction->getSizeOfY();
  unsigned int sizeX = data["x"]->size();
  unsigned int sizeZ = data["z"]->size();

  jXOutput(sizeX, &(*workX)(0), sizeY, &(*(jacobianH[1]))(0, 0), sizeZ, &(*workZ)(0));

  // Rebuilt z from Tmp
  *data["z"] = *workZ;
}

void FirstOrderType1R::computeJacobianG(double, unsigned int)
{
  if (jLInput == NULL)
    RuntimeException::selfThrow("FirstOrderType1R::computeJacobianH() failed; not linked to a plug-in function.");

  SiconosVector *lambda = interaction->getLambdaPtr(0);
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  *workZ = *data["z"];
  *workY = *lambda;

  unsigned int sizeY = lambda->size();
  unsigned int sizeX = data["x"]->size();
  unsigned int sizeZ = data["z"]->size();

  jLInput(sizeY, &(*workY)(0), sizeX, &(*(jacobianG[1]))(0, 0), sizeZ, &(*workZ)(0));

  // Rebuilt z from Tmp
  *data["z"] = *workZ;
}

FirstOrderType1R* FirstOrderType1R::convert(Relation *r)
{
  return dynamic_cast<FirstOrderType1R*>(r);
}

