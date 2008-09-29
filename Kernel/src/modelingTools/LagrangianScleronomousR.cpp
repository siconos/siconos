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

#include "LagrangianScleronomousR.h"
#include "LagrangianRXML.h"
#include "Interaction.h"
#include "LagrangianDS.h"

using namespace std;

//default constructor
LagrangianScleronomousR::LagrangianScleronomousR() : LagrangianR(ScleronomousR), hPtr(NULL), G0Ptr(NULL)
{
  G.resize(1);
}


// xml constructor
LagrangianScleronomousR::LagrangianScleronomousR(RelationXMLSPtr relxml): LagrangianR(relxml, ScleronomousR)
{
  LagrangianRXMLSPtr LRxml = boost::static_pointer_cast<LagrangianRXML>(relationxml);
  // h plug-in
  if (LRxml->hasH())
  {
    pluginNames["h"] = LRxml->getHPlugin();
    setComputeHFunction(cShared.getPluginName(pluginNames["h"]), cShared.getPluginFunctionName(pluginNames["h"]));
  }
  if (!LRxml->hasG())
    RuntimeException::selfThrow("LagrangianScleronomousR:: xml constructor failed, can not find a definition for G0.");

#ifndef WithSmartPtr
  G.resize(1, NULL);
#else
  G.resize(1);
#endif

  // Read G matrix or plug-in names.
  readGInXML(LRxml, 0);
}

// constructor from a set of data
LagrangianScleronomousR::LagrangianScleronomousR(const string& computeH, const std::string& computeG):
  LagrangianR(ScleronomousR), hPtr(NULL), G0Ptr(NULL)
{
  setComputeHFunction(cShared.getPluginName(computeH), cShared.getPluginFunctionName(computeH));
  // Note that in this case, G is not allocated since we do not have its dimensions.
  // That will be done during initialize, with Interaction input.

#ifndef WithSmartPtr
  G.resize(1, NULL);
#else
  G.resize(1);
#endif

  string name = "G0";
  pluginNames[name] = computeG;
  setComputeGFunction(cShared.getPluginName(pluginNames[name]), cShared.getPluginFunctionName(pluginNames[name]), 0);

#ifndef WithSmartPtr
  isAllocatedIn[name] = false;
#endif

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
    SiconosVectorSPtr y = interaction->getYPtr(0);

    // Warning: temporary method to have contiguous values in memory, copy of block to simple.
    *workX = *data["q0"];
    *workZ = *data["z"];
    *workY = *y;

    unsigned int sizeQ = workX->size();
    unsigned int sizeY = y->size();
    unsigned int sizeZ = workZ->size();

    if (hPtr == NULL)
      RuntimeException::selfThrow("LagrangianScleronomousR:computeH() failed, h is not linked to a plugin function");
    hPtr(sizeQ, &(*workX)(0) , sizeY, &(*workY)(0), sizeZ, &(*workZ)(0));

    // Copy data that might have been changed in the plug-in call.
    *data["z"] = *workZ;
    *y = *workY;
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
    *workX = *data["q0"];
    *workZ = *data["z"];

    unsigned int sizeY = G[0]->size(0);
    unsigned int sizeQ = workX->size();
    unsigned int sizeZ = workZ->size();

    if (G0Ptr == NULL)
      RuntimeException::selfThrow("computeG() is not linked to a plugin function");
    G0Ptr(sizeQ, &(*workX)(0), sizeY, &(*(G[0]))(0, 0), sizeZ, &(*workZ)(0));

    // Copy data that might have been changed in the plug-in call.
    *data["z"] = *workZ;
  }
  //  else nothing!
}

void LagrangianScleronomousR::computeOutput(double time, unsigned int derivativeNumber)
{
  if (derivativeNumber == 0)
    computeH(time);
  else
  {
    computeG(time);
    SiconosVectorSPtr y = interaction->getYPtr(derivativeNumber) ;
    if (derivativeNumber == 1)
      prod(*G[0], *data["q1"], *y);
    else if (derivativeNumber == 2)
      prod(*G[0], *data["q2"], *y);
    else
      RuntimeException::selfThrow("LagrangianScleronomousR::computeOutput(t,index), index out of range");
  }
}

void LagrangianScleronomousR::computeInput(double time, unsigned int level)
{
  computeG(time, 0);
  string name = "p" + toString<unsigned int>(level);
  // get lambda of the concerned interaction
  SiconosVectorSPtr lambda = interaction->getLambdaPtr(level);
  // data[name] += trans(G) * lambda
  prod(*lambda, *G[0], *data[name], false);

  //   SP::SiconosMatrix  GT = new SimpleMatrix(*G[0]);
  //   GT->trans();
  //   *data[name] += prod(*GT, *lambda);
  //   delete GT;
  //gemv(CblasTrans, 1.0,*(G[0]), *lambda, 1.0, *data[name]); => not yet implemented for BlockVectors.
}

LagrangianScleronomousR* LagrangianScleronomousR::convert(Relation *r)
{
  return dynamic_cast<LagrangianScleronomousR*>(r);
}
