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
using namespace RELATION;

// Default constructor
LagrangianRheonomousR::LagrangianRheonomousR(): LagrangianR(RheonomousR), hPtr(NULL), hDotPtr(NULL), G0Ptr(NULL)
{
  isPlugged[RELATION::G0] = false ;
  isPlugged[RELATION::hDot] = false ;
  G.resize(1);
}

// xml constructor
LagrangianRheonomousR::LagrangianRheonomousR(SP::RelationXML relxml): LagrangianR(relxml, RheonomousR)
{
  SP::LagrangianRXML LRxml = boost::static_pointer_cast<LagrangianRXML>(relationxml);
  // h plug-in
  if (LRxml->hasH())
  {
    pluginNames[RELATION::h] = LRxml->getHPlugin();
    setComputeHFunction(SSL::getPluginName(pluginNames[RELATION::h]), SSL::getPluginFunctionName(pluginNames[RELATION::h]));
  }

  // Read hDot
  if (!LRxml->hasHDot())
    RuntimeException::selfThrow("LagrangianRheonomousR:: xml constructor failed, can not find a definition for hDot.");
  if (LRxml->isHDotPlugin())
  {
    pluginNames[RELATION::hDot] = LRxml->getHDotPlugin();
    setComputeHDotFunction(SSL::getPluginName(pluginNames[RELATION::hDot]), SSL::getPluginFunctionName(pluginNames[RELATION::hDot]));

  }
  else
  {

    hDot.reset(new SimpleVector(LRxml->getHDotVector()));


    isPlugged[RELATION::hDot] = false;
  }

  if (!LRxml->hasG())
    RuntimeException::selfThrow("LagrangianRheonomousR:: xml constructor failed, can not find a definition for G0.");

  G.resize(1);
  // Read G matrices or plug-in names.
  readGInXML(LRxml, 0);
}

// constructor from a set of data
LagrangianRheonomousR::LagrangianRheonomousR(const string& computeH, const string& computeHDot, const string& computeG):
  LagrangianR(RheonomousR), hPtr(NULL), hDotPtr(NULL), G0Ptr(NULL)
{
  // h
  setComputeHFunction(SSL::getPluginName(computeH), SSL::getPluginFunctionName(computeH));
  // hDot
  setComputeHDotFunction(SSL::getPluginName(computeHDot), SSL::getPluginFunctionName(computeHDot));

  G.resize(1);
  pluginNames[RELATION::G0] = computeG;
  setComputeGFunction(SSL::getPluginName(pluginNames[RELATION::G0]), SSL::getPluginFunctionName(pluginNames[RELATION::G0]), 0);
}

LagrangianRheonomousR::~LagrangianRheonomousR()
{}

void LagrangianRheonomousR::initComponents()
{
  LagrangianR::initComponents();

  unsigned int sizeY = interaction->getSizeOfY();
  // hDot
  if (!hDot)
    hDot.reset(new SimpleVector(sizeY));
  else if (sizeY != hDot->size())
    RuntimeException::selfThrow("LagrangianRheonomousR:: initComponents failed. Inconsistent sizes between Interaction and Relation matrices.");

}

void LagrangianRheonomousR::setHDot(const SiconosVector& newValue)
{
  if (!hDot)
    hDot.reset(new SimpleVector(newValue));
  else
    *(hDot) = newValue;

  isPlugged[RELATION::hDot] = false;
}

void LagrangianRheonomousR::setHDotPtr(SP::SiconosVector newPtr)
{
  hDot = newPtr;
  isPlugged[RELATION::hDot] = false;
}

void LagrangianRheonomousR::setComputeHFunction(const string& pluginPath, const string& functionName)
{
  isPlugged[RELATION::h] = Plugin::setFunction(&hPtr, pluginPath, functionName, pluginNames[RELATION::h]);
}

void LagrangianRheonomousR::setComputeHDotFunction(const string& pluginPath, const string& functionName)
{
  isPlugged[RELATION::hDot] = Plugin::setFunction(&hDotPtr, pluginPath, functionName, pluginNames[RELATION::hDot]);
}

void LagrangianRheonomousR::setComputeGFunction(const string& pluginPath, const string& functionName, unsigned int)
{
  isPlugged[RELATION::G0] = Plugin::setFunction(&G0Ptr, pluginPath, functionName, pluginNames[RELATION::G0]);
}

void LagrangianRheonomousR::computeH(double time)
{
  if (isPlugged[RELATION::h])
  {
    // get vector y of the current interaction
    SP::SiconosVector y = interaction->getYPtr(0);

    // Warning: temporary method to have contiguous values in memory, copy of block to simple.
    *workX = *data["q0"];
    *workZ = *data["z"];
    *workY = *y;

    unsigned int sizeQ = workX->size();
    unsigned int sizeY = y->size();
    unsigned int sizeZ = workZ->size();

    if (!hPtr)
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
  if (isPlugged[RELATION::hDot])
  {
    // Warning: temporary method to have contiguous values in memory, copy of block to simple.
    *workX = *data["q0"];
    *workZ = *data["z"];

    unsigned int sizeQ = workX->size();
    unsigned int sizeY = hDot->size();
    unsigned int sizeZ = workZ->size();

    if (!hDotPtr)
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
  if (isPlugged[RELATION::G0])
  {
    // Warning: temporary method to have contiguous values in memory, copy of block to simple.
    *workX = *data["q0"];
    *workZ = *data["z"];

    unsigned int sizeY = G[0]->size(0);
    unsigned int sizeQ = workX->size();
    unsigned int sizeZ = workZ->size();
    if (!G0Ptr)
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
    SP::SiconosVector y = interaction->getYPtr(derivativeNumber);
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
  SP::SiconosVector lambda = interaction->getLambdaPtr(level);
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

