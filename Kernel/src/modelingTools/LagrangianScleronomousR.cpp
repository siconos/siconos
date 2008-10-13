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
using namespace RELATION;

//default constructor
LagrangianScleronomousR::LagrangianScleronomousR() : LagrangianR(ScleronomousR), hPtr(NULL), G0Ptr(NULL)
{
  G.resize(1);
}


// xml constructor
LagrangianScleronomousR::LagrangianScleronomousR(SP::RelationXML relxml): LagrangianR(relxml, ScleronomousR)
{
  SP::LagrangianRXML LRxml = boost::static_pointer_cast<LagrangianRXML>(relationxml);
  // h plug-in
  if (LRxml->hasH())
  {
    pluginNames[RELATION::h] = LRxml->getHPlugin();
    setComputeHFunction(SSL::getPluginName(pluginNames[RELATION::h]), SSL::getPluginFunctionName(pluginNames[RELATION::h]));
  }
  if (!LRxml->hasG())
    RuntimeException::selfThrow("LagrangianScleronomousR:: xml constructor failed, can not find a definition for G0.");
  G.resize(1);

  // Read G matrix or plug-in names.
  readGInXML(LRxml, 0);
}

// constructor from a set of data
LagrangianScleronomousR::LagrangianScleronomousR(const string& computeH, const std::string& computeG):
  LagrangianR(ScleronomousR), hPtr(NULL), G0Ptr(NULL)
{
  setComputeHFunction(SSL::getPluginName(computeH), SSL::getPluginFunctionName(computeH));
  // Note that in this case, G is not allocated since we do not have its dimensions.
  // That will be done during initialize, with Interaction input.

  G.resize(1);
  pluginNames[RELATION::G0] = computeG;
  setComputeGFunction(SSL::getPluginName(pluginNames[RELATION::G0]), SSL::getPluginFunctionName(pluginNames[RELATION::G0]), 0);
}

LagrangianScleronomousR::~LagrangianScleronomousR()
{}

void LagrangianScleronomousR::setComputeHFunction(const string& pluginPath, const string& functionName)
{
  isPlugged[RELATION::h] = Plugin::setFunction(&hPtr, pluginPath, functionName, pluginNames[RELATION::h]);
}

void LagrangianScleronomousR::setComputeGFunction(const string& pluginPath, const string& functionName, unsigned int)
{
  isPlugged[RELATION::G0] = Plugin::setFunction(&G0Ptr, pluginPath, functionName, pluginNames[RELATION::G0]);
}

void LagrangianScleronomousR::computeH(double)
{
  // arg= time. Unused in this function but required for interface.
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
    if (G[0])
    {
      computeG(time);
      SP::SiconosVector y = interaction->getYPtr(derivativeNumber) ;
      if (derivativeNumber == 1)
        prod(*G[0], *data["q1"], *y);
      else if (derivativeNumber == 2)
        prod(*G[0], *data["q2"], *y);
      else
        RuntimeException::selfThrow("LagrangianScleronomousR::computeOutput(t,index), index out of range");
    }
  }
}

void LagrangianScleronomousR::computeInput(double time, unsigned int level)
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
  //gemv(CblasTrans, 1.0,*(G[0]), *lambda, 1.0, *data[name]); => not yet implemented for BlockVectors.
}

LagrangianScleronomousR* LagrangianScleronomousR::convert(Relation *r)
{
  return dynamic_cast<LagrangianScleronomousR*>(r);
}
