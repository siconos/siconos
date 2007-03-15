/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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

// Default constructor
FirstOrderR::FirstOrderR(const string& newType):
  Relation("FirstOrder", newType), firstOrderType(newType), computeOutputPtr(NULL), computeInputPtr(NULL)
{
  setComputeOutputFunction("DefaultPlugin.so", "computeOutput");
  isPlugged["output"] = false;
  setComputeInputFunction("DefaultPlugin.so", "computeInput");
  isPlugged["input"] = false;
}

// xml constructor
FirstOrderR::FirstOrderR(RelationXML* relxml, const string& newType):
  Relation(relxml, "FirstOrder", newType), firstOrderType(newType), computeOutputPtr(NULL), computeInputPtr(NULL)
{
  FirstOrderRXML * FORxml = static_cast<FirstOrderRXML*>(relationxml);
  string plugin;

  // computeInput
  if (FORxml->hasComputeInput())
  {
    plugin = (FORxml)->getComputeInputPlugin();
    setComputeInputFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
  }
  else
  {
    setComputeInputFunction("DefaultPlugin.so", "computeInput");
    isPlugged["input"] = false; //
  }

  // computeOutput
  if (FORxml->hasComputeOutput())
  {
    plugin = (FORxml)->getComputeOutputPlugin();
    setComputeOutputFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
  }
  else
  {
    setComputeOutputFunction("DefaultPlugin.so", "computeOutput");
    isPlugged["output"] = false; //
  }
}

FirstOrderR::FirstOrderR(const string& computeOut, const string& computeIn):
  Relation("FirstOrder", "R"), firstOrderType("R"), computeOutputPtr(NULL), computeInputPtr(NULL)
{
  pluginNames["output"] = computeOut;
  setComputeOutputFunction(cShared.getPluginName(pluginNames["output"]), cShared.getPluginFunctionName(pluginNames["output"]));
  pluginNames["input"] = computeIn;
  setComputeInputFunction(cShared.getPluginName(pluginNames["input"]), cShared.getPluginFunctionName(pluginNames["input"]));
}

FirstOrderR::~FirstOrderR()
{
  computeOutputPtr = NULL;
  computeInputPtr = NULL;
}

void FirstOrderR::initialize()
{
  // Check if an Interaction is connected to the Relation.
  if (interaction == NULL)
    RuntimeException::selfThrow("FirstOrderR::initialize failed. No Interaction linked to the present relation.");

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
    data["x"]->addPtr(ds->getXPtr());
    data["xFree"]->addPtr(ds->getXFreePtr());
    data["z"]->addPtr(ds->getZPtr());
    data["r"]->addPtr(ds->getRPtr());
  }

}

void FirstOrderR::setComputeOutputFunction(const string& pluginPath, const string& functionName)
{
  cShared.setFunction(&computeOutputPtr, pluginPath, functionName);
  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames["output"] = plugin + ":" + functionName;
  isPlugged["output"] = true;
}

void FirstOrderR::setComputeInputFunction(const string& pluginPath, const string& functionName)
{
  cShared.setFunction(&computeInputPtr, pluginPath, functionName);
  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames["input"] = plugin + ":" + functionName;
  isPlugged["input"] = true;
}

void FirstOrderR::computeOutput(double time, unsigned int)
{
  // Note that the second argument remains unamed since it is not used: for first order systems, we always compute
  // y[0] (at the time).

  if (computeOutputPtr == NULL)
    RuntimeException::selfThrow("computeOutput() is not linked to a plugin function");

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

  computeOutputPtr(sizeX, &(*xCopy)(0), time, sizeY, &(*lambdaCopy)(0), &(*yCopy)(0), sizeZ, &(*zCopy)(0));

  // Rebuilt lambda/y from Tmp
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

  if (computeOutputPtr == NULL)
    RuntimeException::selfThrow("computeOutput() is not linked to a plugin function");

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

  computeOutputPtr(sizeX, &(*xCopy)(0), time, sizeY, &(*lambdaCopy)(0), &(*yCopy)(0), sizeZ, &(*zCopy)(0));

  // Rebuilt lambda/y from Tmp
  *y = *yCopy;
  *data["z"] = *zCopy;

  delete zCopy;
  delete lambdaCopy;
  delete yCopy;
  delete xCopy;
}

void FirstOrderR::computeInput(double time, unsigned int level)
{
  if (computeInputPtr == NULL)
    RuntimeException::selfThrow("computeInput() is not linked to a plugin function");

  SiconosVector *lambda = interaction->getLambdaPtr(level);
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  SimpleVector * rCopy = new SimpleVector(*data["r"]);
  SimpleVector * zCopy = new SimpleVector(*data["z"]);
  SimpleVector * lambdaCopy = new SimpleVector(*lambda);

  unsigned int sizeY = lambda->size();
  unsigned int sizeZ = data["z"]->size();

  computeInputPtr(sizeY, &(*lambdaCopy)(0), time, &(*rCopy)(0), sizeZ, &(*zCopy)(0));

  *data["r"] = *rCopy;
  *data["z"] = *zCopy;

  delete rCopy;
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

