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
using namespace RELATION;

void FirstOrderType1R::initPluginFlags(bool in)
{
  isPlugged[RELATION::h] = in;
  isPlugged[RELATION::jacobianH0] = in;
  isPlugged[RELATION::g] = in;
  isPlugged[RELATION::jacobianG0] = in ;
}

// xml constructor
FirstOrderType1R::FirstOrderType1R(SP::RelationXML relxml):
  FirstOrderR(relxml, Type1R)
{
  SP::FirstOrderRXML FORxml = boost::static_pointer_cast<FirstOrderRXML>(relationxml);

  initPluginFlags(false);
  // Gradients

  jacobianH.resize(1);
  jacobianG.resize(1);

  // input g
  if (FORxml->hasG())
  {
    pluginNames[RELATION::g] = FORxml->getGPlugin();
    setComputeGFunction(SSL::getPluginName(pluginNames[RELATION::g]), SSL::getPluginFunctionName(pluginNames[RELATION::g]));
    // Gradients
    if (!FORxml->hasJacobianG())
      RuntimeException::selfThrow("FirstOrderType1R xml constructor failed. No input for gradient(s) of g function.");

    if (FORxml->isJacobianGPlugin(0))
    {
      pluginNames[RELATION::jacobianG0] = FORxml->getJacobianGPlugin(0);
      setComputeJacobianGFunction(SSL::getPluginName(pluginNames[RELATION::jacobianG0]), SSL::getPluginFunctionName(pluginNames[RELATION::jacobianG0]));
    }
    else
      jacobianG[0].reset(new SimpleMatrix(FORxml->getJacobianGMatrix(0)));
  }

  // output h
  if (FORxml->hasH())
  {
    pluginNames[RELATION::h] = FORxml->getHPlugin();
    setComputeHFunction(SSL::getPluginName(pluginNames[RELATION::h]), SSL::getPluginFunctionName(pluginNames[RELATION::h]));
    // Gradients
    if (!FORxml->hasJacobianH())
      RuntimeException::selfThrow("FirstOrderType1R xml constructor failed. No input for gradients of h function.");
    if (FORxml->isJacobianHPlugin(0))
    {
      pluginNames[RELATION::jacobianH0] = FORxml->getJacobianHPlugin(0);
      setComputeJacobianHFunction(SSL::getPluginName(pluginNames[RELATION::jacobianH0]), SSL::getPluginFunctionName(pluginNames[RELATION::jacobianH0]));
    }
    else
      jacobianH[0].reset(new SimpleMatrix(FORxml->getJacobianHMatrix(0)));
  }
}

FirstOrderType1R::FirstOrderType1R(const string& computeOut, const string& computeIn):
  FirstOrderR(Type1R), output(NULL), jXOutput(NULL), input(NULL), jLInput(NULL)
{
  initPluginFlags(false);
  // Size vector of pointers to functions.
  // Connect input and output to plug-in
  setComputeHFunction(SSL::getPluginName(computeOut), SSL::getPluginFunctionName(computeOut));
  setComputeGFunction(SSL::getPluginName(computeIn), SSL::getPluginFunctionName(computeIn));

  jacobianH.resize(1);
  jacobianG.resize(1);

  // The jacobians are not set, and thus considered as null matrices at this point.
}

FirstOrderType1R::FirstOrderType1R(const string& computeOut, const string& computeIn, const string& computeJX, const string& computeJL):
  FirstOrderR(Type1R), output(NULL), jXOutput(NULL), input(NULL), jLInput(NULL)
{
  initPluginFlags(false);
  // Size vector of pointers to functions.

  jacobianH.resize(1);
  jacobianG.resize(1);

  // Connect input and output to plug-in
  setComputeHFunction(SSL::getPluginName(computeOut), SSL::getPluginFunctionName(computeOut));
  setComputeGFunction(SSL::getPluginName(computeIn), SSL::getPluginFunctionName(computeIn));
  setComputeJacobianHFunction(SSL::getPluginName(computeJX), SSL::getPluginFunctionName(computeJX));
  setComputeJacobianGFunction(SSL::getPluginName(computeJL), SSL::getPluginFunctionName(computeJL));
}

FirstOrderType1R::~FirstOrderType1R()
{}

void FirstOrderType1R::initialize()
{
  FirstOrderR::initialize();

  // if jacobianH0 is plugged and memory not allocated ...
  if (isPlugged[RELATION::jacobianH0] && ! jacobianH[0])
    jacobianH[0].reset(new SimpleMatrix(interaction->getSizeOfY(), interaction->getSizeOfDS()));
  // Same work for jacobianLambdaG
  if (isPlugged[RELATION::jacobianG0] && ! jacobianG[0])
    jacobianG[0].reset(new SimpleMatrix(interaction->getSizeOfDS(), interaction->getSizeOfY()));
}

void FirstOrderType1R::setComputeHFunction(const string& pluginPath, const string& functionName)
{
  isPlugged[RELATION::h] = Plugin::setFunction(&output, pluginPath, functionName, pluginNames[RELATION::h]);
}

void FirstOrderType1R::setComputeJacobianHFunction(const string& pluginPath, const string& functionName, unsigned int)
{
  isPlugged[RELATION::jacobianH0] = Plugin::setFunction(&jXOutput, pluginPath, functionName, pluginNames[RELATION::jacobianH0]);
}

void FirstOrderType1R::setComputeGFunction(const string& pluginPath, const string& functionName)
{
  isPlugged[RELATION::g] = Plugin::setFunction(&input, pluginPath, functionName, pluginNames[RELATION::g]);
}

void FirstOrderType1R::setComputeJacobianGFunction(const string& pluginPath, const string& functionName, unsigned int)
{
  isPlugged[RELATION::jacobianG0] = Plugin::setFunction(&jLInput, pluginPath, functionName, pluginNames[RELATION::jacobianG0]);
}

void FirstOrderType1R::computeOutput(double, unsigned int)
{
  if (!output)
    RuntimeException::selfThrow("FirstOrderType1R::computeOutput() is not linked to a plugin function");

  SP::SiconosVector y = interaction->getYPtr(0);
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
  if (!input)
    RuntimeException::selfThrow("FirstOrderType1R::computeInput() is not linked to a plugin function");

  SP::SiconosVector lambda = interaction->getLambdaPtr(level);
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
  if (!jXOutput)
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
  if (!jLInput)
    RuntimeException::selfThrow("FirstOrderType1R::computeJacobianH() failed; not linked to a plug-in function.");

  SP::SiconosVector lambda = interaction->getLambdaPtr(0);
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

