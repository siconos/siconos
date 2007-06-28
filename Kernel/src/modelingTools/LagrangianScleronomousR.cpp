/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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

#include "LagrangianScleronomousR.h"
#include "LagrangianRXML.h"
#include "Interaction.h"
#include "LagrangianDS.h"

using namespace std;

// Default constructor
LagrangianScleronomousR::LagrangianScleronomousR():
  LagrangianR("ScleronomousR"), hPtr(NULL), G0Ptr(NULL)
{
  isPlugged["G0"] = false ;
  isAllocatedIn["G0"] = false;
  G.resize(1, NULL);
}

// xml constructor
LagrangianScleronomousR::LagrangianScleronomousR(RelationXML* relxml): LagrangianR(relxml, "ScleronomousR"), hPtr(NULL), G0Ptr(NULL)
{
  LagrangianRXML * LRxml = static_cast<LagrangianRXML *>(relationxml);
  // h plug-in
  if (LRxml->hasH())
  {
    pluginNames["h"] = LRxml->getHPlugin();
    setComputeHFunction(cShared.getPluginName(pluginNames["h"]), cShared.getPluginFunctionName(pluginNames["h"]));
  }
  if (!LRxml->hasG())
    RuntimeException::selfThrow("LagrangianScleronomousR:: xml constructor failed, can not find a definition for G0.");
  G.resize(1, NULL);
  // Read G matrix or plug-in names.
  readGInXML(LRxml, 0);
}

// constructor from a set of data
LagrangianScleronomousR::LagrangianScleronomousR(const string& computeH, const std::string& computeG):
  LagrangianR("ScleronomousR"), hPtr(NULL), G0Ptr(NULL)
{
  setComputeHFunction(cShared.getPluginName(computeH), cShared.getPluginFunctionName(computeH));
  // Note that in this case, G is not allocated since we do not have its dimensions.
  // That will be done during initialize, with Interaction input.
  G.resize(1, NULL);
  string name = "G0";
  pluginNames[name] = computeG;
  setComputeGFunction(cShared.getPluginName(pluginNames[name]), cShared.getPluginFunctionName(pluginNames[name]), 0);
  isAllocatedIn[name] = false;
}

LagrangianScleronomousR::~LagrangianScleronomousR()
{
  hPtr = NULL;
  G0Ptr = NULL;
}

void LagrangianScleronomousR::setComputeHFunction(const string& pluginPath, const string& functionName)
{
  cShared.setFunction(&hPtr, pluginPath, functionName);
  string plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames["h"] = plugin + ":" + functionName;
  isPlugged["h"] = true;
}

void LagrangianScleronomousR::setComputeGFunction(const string& pluginPath, const string& functionName, unsigned int)
{
  cShared.setFunction(&G0Ptr, pluginPath, functionName);
  string plugin = pluginPath.substr(0, pluginPath.length() - 3);
  string name = "G0";
  pluginNames["G0"] = plugin + ":" + functionName;
  isPlugged["G0"] = true;
}

void LagrangianScleronomousR::computeH(double)
{
  // arg= time. Unused in this function but required for interface.
  if (isPlugged["h"])
  {
    // get vector y of the current interaction
    SiconosVector *y = interaction->getYPtr(0);

    // Warning: temporary method to have contiguous values in memory, copy of block to simple.
    SimpleVector * qCopy = new SimpleVector(*data["q0"]);
    SimpleVector * zCopy = new SimpleVector(*data["z"]);
    SimpleVector * yCopy = new SimpleVector(*y);

    unsigned int sizeQ = qCopy->size();
    unsigned int sizeY = y->size();
    unsigned int sizeZ = zCopy->size();

    if (hPtr == NULL)
      RuntimeException::selfThrow("LagrangianScleronomousR:computeH() failed, h is not linked to a plugin function");
    hPtr(sizeQ, &(*qCopy)(0) , sizeY, &(*yCopy)(0), sizeZ, &(*zCopy)(0));

    // Copy data that might have been changed in the plug-in call.
    *data["z"] = *zCopy;
    *y = *yCopy;
    delete qCopy;
    delete yCopy;
    delete zCopy;
  }
  // else nothing
}

void LagrangianScleronomousR::computeG(double, unsigned int)
{
  // First arg: time. Useless.
  // Last arg: index for G - Useless, always equal to 0 for this kind of relation.

  if (isPlugged["G0"])
  {
    // Warning: temporary method to have contiguous values in memory, copy of block to simple.
    SimpleVector * qCopy = new SimpleVector(*data["q0"]);
    SimpleVector * zCopy = new SimpleVector(*data["z"]);

    unsigned int sizeY = G[0]->size(0);
    unsigned int sizeQ = qCopy->size();
    unsigned int sizeZ = zCopy->size();

    if (G0Ptr == NULL)
      RuntimeException::selfThrow("computeG() is not linked to a plugin function");
    G0Ptr(sizeQ, &(*qCopy)(0), sizeY, &(*(G[0]))(0, 0), sizeZ, &(*zCopy)(0));

    // Copy data that might have been changed in the plug-in call.
    *data["z"] = *zCopy;
    delete qCopy;
    delete zCopy;
  }
  //  else nothing!
}

void LagrangianScleronomousR::computeHFree(double)
{
  if (isPlugged["h"])
  {
    // arg= time. Unused in this function but required for interface.

    // get vector y of the current interaction
    SiconosVector *y = interaction->getYPtr(0);

    // Warning: temporary method to have contiguous values in memory, copy of block to simple.
    SimpleVector * qCopy = new SimpleVector(*data["q0Free"]);
    SimpleVector * zCopy = new SimpleVector(*data["z"]);
    SimpleVector * yCopy = new SimpleVector(*y);

    unsigned int sizeQ = qCopy->size();
    unsigned int sizeY = y->size();
    unsigned int sizeZ = zCopy->size();

    if (hPtr == NULL)
      RuntimeException::selfThrow("LagrangianScleronomousR:computeH() is not linked to a plugin function");
    hPtr(sizeQ, &(*qCopy)(0) , sizeY, &(*yCopy)(0), sizeZ, &(*zCopy)(0));

    // Copy data that might have been changed in the plug-in call.
    *data["z"] = *zCopy;
    *y = *yCopy;
    delete qCopy;
    delete yCopy;
    delete zCopy;
  }
  //else nothing
}

void LagrangianScleronomousR::computeGFree(double, unsigned int)
{
  // First arg: time. Useless.
  // Last arg: index for G - Useless, always equal to 0 for this kind of relation.

  if (isPlugged["G0"])
  {
    unsigned int sizeY = interaction->getSizeOfY();
    unsigned int sizeQ = interaction->getSizeOfDS();

    // Warning: temporary method to have contiguous values in memory, copy of block to simple.
    SimpleVector * qCopy = new SimpleVector(*data["q0Free"]);
    SimpleVector * zCopy = new SimpleVector(*data["z"]);
    unsigned int sizeZ = zCopy->size();

    if (G0Ptr == NULL)
      RuntimeException::selfThrow("computeG() is not linked to a plugin function");
    G0Ptr(sizeQ, &(*qCopy)(0), sizeY, &(*(G[0]))(0, 0), sizeZ, &(*zCopy)(0));

    // Copy data that might have been changed in the plug-in call.
    *data["z"] = *zCopy;
    delete qCopy;
    delete zCopy;
  }
  //else nothing
}

void LagrangianScleronomousR::computeOutput(double time, unsigned int derivativeNumber)
{
  if (derivativeNumber == 0)
    computeH(time);
  else
  {
    computeG(time);
    SiconosVector *y = interaction->getYPtr(derivativeNumber) ;
    if (derivativeNumber == 1)
      *y = prod(*G[0], *data["q1"]);
    else if (derivativeNumber == 2)
      *y = prod(*G[0], *data["q2"]);
    else
      RuntimeException::selfThrow("LagrangianScleronomousR::computeOutput(t,index), index out of range");
  }
}

void LagrangianScleronomousR::computeFreeOutput(double time, unsigned int derivativeNumber)
{
  if (derivativeNumber == 0)
    computeHFree(time);

  else
  {
    computeGFree(time);
    SiconosVector *y = interaction->getYPtr(derivativeNumber) ;
    if (derivativeNumber == 1)
      *y = prod(*G[0], *data["q1Free"]);
    else if (derivativeNumber == 2)
      *y = prod(*G[0], *data["q2"]);
    else
      RuntimeException::selfThrow("LagrangianScleronomousR::computeOutput(t,index), index out of range");
  }
}

void LagrangianScleronomousR::computeInput(double time, unsigned int level)
{
  computeG(time, 0);
  string name = "p" + toString<unsigned int>(level);
  // get lambda of the concerned interaction
  SiconosVector *lambda = interaction->getLambdaPtr(level);
  SiconosMatrix * GT = new SimpleMatrix(*G[0]);
  GT->trans();
  *data[name] += prod(*GT, *lambda);
  delete GT;
  //gemv(CblasTrans, 1.0,*(G[0]), *lambda, 1.0, *data[name]); => not yet implemented for BlockVectors.
}

LagrangianScleronomousR* LagrangianScleronomousR::convert(Relation *r)
{
  return dynamic_cast<LagrangianScleronomousR*>(r);
}
