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

#include "LagrangianRheonomousR.h"
#include "LagrangianRXML.h"
#include "Interaction.h"
#include "LagrangianDS.h"

using namespace std;

// Default constructor
LagrangianRheonomousR::LagrangianRheonomousR(): LagrangianR(RheonomousR), hPtr(NULL), hDotPtr(NULL), G0Ptr(NULL)
{
  isPlugged["G0"] = false ;
  isPlugged["hDot"] = false ;

#ifndef WithSmartPtr
  isAllocatedIn["G0"] = false;
  isAllocatedIn["hDot"] = false;
  G.resize(1, NULL);
#else
  G.resize(1);
#endif

}

// xml constructor
LagrangianRheonomousR::LagrangianRheonomousR(RelationXMLSPtr relxml): LagrangianR(relxml, RheonomousR)
{
  LagrangianRXMLSPtr LRxml = boost::static_pointer_cast<LagrangianRXML>(relationxml);
  // h plug-in
  if (LRxml->hasH())
  {
    pluginNames["h"] = LRxml->getHPlugin();
    setComputeHFunction(cShared.getPluginName(pluginNames["h"]), cShared.getPluginFunctionName(pluginNames["h"]));
  }

  // Read hDot
  if (!LRxml->hasHDot())
    RuntimeException::selfThrow("LagrangianRheonomousR:: xml constructor failed, can not find a definition for hDot.");
  if (LRxml->isHDotPlugin())
  {
    pluginNames["hDot"] = LRxml->getHDotPlugin();
    setComputeHDotFunction(cShared.getPluginName(pluginNames["hDot"]), cShared.getPluginFunctionName(pluginNames["hDot"]));

  }
  else
  {

    hDot.reset(new SimpleVector(LRxml->getHDotVector()));


    isPlugged["hDot"] = false;
  }

  if (!LRxml->hasG())
    RuntimeException::selfThrow("LagrangianRheonomousR:: xml constructor failed, can not find a definition for G0.");

#ifndef WithSmartPtr
  G.resize(1, NULL);
#else
  G.resize(1);
#endif

  // Read G matrices or plug-in names.
  readGInXML(LRxml, 0);
}

// constructor from a set of data
LagrangianRheonomousR::LagrangianRheonomousR(const string& computeH, const string& computeHDot, const string& computeG):
  LagrangianR(RheonomousR), hPtr(NULL), hDotPtr(NULL), G0Ptr(NULL)
{
  // h
  setComputeHFunction(cShared.getPluginName(computeH), cShared.getPluginFunctionName(computeH));
  // hDot
  setComputeHDotFunction(cShared.getPluginName(computeHDot), cShared.getPluginFunctionName(computeHDot));

#ifndef WithSmartPtr
  isAllocatedIn["hDot"] = false;
  // G[0]

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

LagrangianRheonomousR::~LagrangianRheonomousR()
{
  hPtr = NULL;
  G0Ptr = NULL;
  hDotPtr = NULL;
}

void LagrangianRheonomousR::initComponents()
{
  LagrangianR::initComponents();

  unsigned int sizeY = interaction->getSizeOfY();
  // hDot
  if (hDot == NULL)
  {

#ifndef WithSmartPtr
    hDot = new SimpleVector(sizeY);
    isAllocatedIn["hDot"] = true;
#else
    hDot.reset(new SimpleVector(sizeY));
#endif

  }
  else if (sizeY != hDot->size())
    RuntimeException::selfThrow("LagrangianRheonomousR:: initComponents failed. Inconsistent sizes between Interaction and Relation matrices.");

}

void LagrangianRheonomousR::setHDot(const SiconosVector& newValue)
{
  if (hDot == NULL)
  {
#ifndef WithSmartPtr
    hDot =  new SimpleVector(newValue);
    isAllocatedIn["hDot"] = true;
#else
    hDot.reset(new SimpleVector(newValue));
#endif

  }
  else
    *(hDot) = newValue;

  isPlugged["hDot"] = false;
}

void LagrangianRheonomousR::setHDotPtr(SiconosVectorSPtr newPtr)
{

#ifndef WithSmartPtr
  if (isAllocatedIn["hDot"]) delete hDot;
  isAllocatedIn["hDot"] = false;
#endif

  hDot = newPtr;

  isPlugged["hDot"] = false;
}

void LagrangianRheonomousR::setComputeHFunction(const string& pluginPath, const string& functionName)
{
  cShared.setFunction(&hPtr, pluginPath, functionName);
  string plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames["h"] = plugin + ":" + functionName;
  isPlugged["h"] = true;
}

void LagrangianRheonomousR::setComputeHDotFunction(const string& pluginPath, const string& functionName)
{
  cShared.setFunction(&hDotPtr, pluginPath, functionName);
  string plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames["hDot"] = plugin + ":" + functionName;
  isPlugged["hDot"] = true;
}

void LagrangianRheonomousR::setComputeGFunction(const string& pluginPath, const string& functionName, unsigned int index)
{
  if (index != 0)
    RuntimeException::selfThrow("LagrangianRheonomousR:: setComputeGFunction(...,index), index out of range");

  cShared.setFunction(&G0Ptr, pluginPath, functionName);
  string plugin = pluginPath.substr(0, pluginPath.length() - 3);
  string name = "G" + toString<unsigned int>(index);
  pluginNames[name] = plugin + ":" + functionName;
  isPlugged[name] = true;
}

void LagrangianRheonomousR::computeH(double time)
{
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
      RuntimeException::selfThrow("LagrangianRheonomousR:computeH() failed, h is not linked to a plugin function");
    hPtr(sizeQ, &(*workX)(0), time, sizeY,  &(*workY)(0), sizeZ, &(*workZ)(0));

    // Copy data that might have been changed in the plug-in call.
    *data["z"] = *workZ;
    *y = *workY;
  }
  // else nothing
}

void LagrangianRheonomousR::computeHDot(double time)
{
  if (isPlugged["hDot"])
  {
    // Warning: temporary method to have contiguous values in memory, copy of block to simple.
    *workX = *data["q0"];
    *workZ = *data["z"];

    unsigned int sizeQ = workX->size();
    unsigned int sizeY = hDot->size();
    unsigned int sizeZ = workZ->size();

    if (hDotPtr == NULL)
      RuntimeException::selfThrow("LagrangianRheonomousR:computeHDot() failed, hDot is not linked to a plugin function");
    hDotPtr(sizeQ, &(*workX)(0), time, sizeY,  &(*hDot)(0), sizeZ, &(*workZ)(0));

    // Copy data that might have been changed in the plug-in call.
    *data["z"] = *workZ;
  }
  // else nothing
}

void LagrangianRheonomousR::computeG(double time, unsigned int)
{
  // Note that second input arg is useless.
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
    G0Ptr(sizeQ, &(*workX)(0), time, sizeY, &(*G[0])(0, 0), sizeZ, &(*workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *data["z"] = *workZ;
  }
  // else nothing.
}

void LagrangianRheonomousR::computeOutput(double time, unsigned int derivativeNumber)
{
  if (derivativeNumber == 0)
    computeH(time);
  else
  {
    SiconosVectorSPtr y = interaction->getYPtr(derivativeNumber);
    computeG(time);
    if (derivativeNumber == 1)
    {
      computeHDot(time); // \todo: save hDot directly into y[1] ?
      prod(*G[0], *data["q1"], *y);
      *y += *hDot;
    }
    else if (derivativeNumber == 2)
      prod(*G[0], *data["q2"], *y); // Approx: y[2] = Gq[2], other terms are neglected ...
    else
      RuntimeException::selfThrow("LagrangianRheonomousR::computeOutput(time,index), index out of range or not yet implemented.");
  }
}

void LagrangianRheonomousR::computeInput(double time, unsigned int level)
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
  //  gemv(CblasTrans, 1.0,*(G[0]), *lambda, 1.0, *data[name]);
}

LagrangianRheonomousR* LagrangianRheonomousR::convert(Relation *r)
{
  return dynamic_cast<LagrangianRheonomousR*>(r);
}

