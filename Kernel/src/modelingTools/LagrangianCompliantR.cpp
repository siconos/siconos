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

#include "LagrangianCompliantR.h"
#include "LagrangianRXML.h"
#include "Interaction.h"
#include "LagrangianDS.h"

using namespace std;
using namespace RELATION;

// Default constructor
LagrangianCompliantR::LagrangianCompliantR(): LagrangianR(CompliantR), hPtr(NULL), G0Ptr(NULL), G1Ptr(NULL)
{
  isPlugged[RELATION::G0] = false ;
  isPlugged[RELATION::G1] = false ;
  G.resize(2);
}

// xml constructor
LagrangianCompliantR::LagrangianCompliantR(SP::RelationXML relxml): LagrangianR(relxml, CompliantR), hPtr(NULL), G0Ptr(NULL), G1Ptr(NULL)
{
  SP::LagrangianRXML LRxml = boost::static_pointer_cast<LagrangianRXML>(relationxml);
  // h plug-in
  if (LRxml->hasH())
  {
    pluginNames[RELATION::h] = LRxml->getHPlugin();
    setComputeHFunction(SSL::getPluginName(pluginNames[RELATION::h]), SSL::getPluginFunctionName(pluginNames[RELATION::h]));
  }

  if (!LRxml->hasG())
    RuntimeException::selfThrow("LagrangianCompliantR:: xml constructor failed, can not find a definition for G0.");
  G.resize(2);

  // Read G matrices or plug-in names.
  readGInXML(LRxml, 0);
  readGInXML(LRxml, 1);
}

// constructor from a set of data
LagrangianCompliantR::LagrangianCompliantR(const string& computeH, const std::vector<string> & computeG): LagrangianR(CompliantR), hPtr(NULL), G0Ptr(NULL), G1Ptr(NULL)
{
  // h
  setComputeHFunction(SSL::getPluginName(computeH), SSL::getPluginFunctionName(computeH));
  // G
  unsigned int nG = 2;

  G.resize(nG);

  RELATION::PluginNames name;
  for (unsigned int i = 0; i < nG; ++i)
  {
    if (i == 0) name = G0;
    else if (i == 1) name = G1;
    else name = G2;
    pluginNames[name] = computeG[i];
    setComputeGFunction(SSL::getPluginName(pluginNames[name]), SSL::getPluginFunctionName(pluginNames[name]), i);
  }
}

LagrangianCompliantR::~LagrangianCompliantR()
{}

void LagrangianCompliantR::initComponents()
{
  LagrangianR::initComponents();
  unsigned int sizeY = interaction->getSizeOfY();
  workL.reset(new SimpleVector(sizeY));

  if (! G[1])
    G[1].reset(new SimpleMatrix(sizeY, sizeY));
  else if (sizeY != G[1]->size(0) || sizeY != G[1]->size(1))
    RuntimeException::selfThrow("LagrangianCompliantR:: initComponents failed. Inconsistent sizes between Interaction and Relation matrices.");
}

void LagrangianCompliantR::setComputeHFunction(const string& pluginPath, const string& functionName)
{
  isPlugged[RELATION::h] = Plugin::setFunction(&hPtr, pluginPath, functionName, pluginNames[RELATION::h]);
}

void LagrangianCompliantR::setComputeGFunction(const string& pluginPath, const string& functionName, unsigned int index)
{
  if (index == 0)
    isPlugged[RELATION::G0] = Plugin::setFunction(&G0Ptr, pluginPath, functionName, pluginNames[RELATION::G0]);
  else
    isPlugged[RELATION::G1] = Plugin::setFunction(&G1Ptr, pluginPath, functionName, pluginNames[RELATION::G1]);
}

void LagrangianCompliantR::computeH(double time)
{
  if (isPlugged[RELATION::h])
  {
    // get vector y of the current interaction
    SP::SiconosVector y = interaction->getYPtr(0);
    SP::SiconosVector lambda = interaction->getLambdaPtr(0);

    // Warning: temporary method to have contiguous values in memory, copy of block to simple.
    *workX = *data["q0"];
    *workZ = *data["z"];
    *workY = *y;
    *workL = *lambda;

    unsigned int sizeQ = workX->size();
    unsigned int sizeY = y->size();
    unsigned int sizeZ = workZ->size();

    if (!hPtr)
      RuntimeException::selfThrow("LagrangianCompliantR:computeH() failed, h is not linked to a plugin function");
    hPtr(sizeQ, &(*workX)(0), sizeY, &(*workL)(0), &(*workY)(0), sizeZ, &(*workZ)(0));

    // Copy data that might have been changed in the plug-in call.
    *data["z"] = *workZ;
    *y = *workY;
  }
}

void LagrangianCompliantR::computeG(double time, unsigned int  index)
{
  if (index >= G.size())
    RuntimeException::selfThrow("LagrangianCompliantR:: computeG(index), index out of range.");

  RELATION::PluginNames name;
  if (index == 0) name = G0;
  else if (index == 1) name = G1;
  else name = G2;

  if (isPlugged[name])
  {
    // Warning: temporary method to have contiguous values in memory, copy of block to simple.
    *workX = *data["q0"];
    *workZ = *data["z"];

    unsigned int sizeY = G[0]->size(0);
    unsigned int sizeQ = workX->size();
    unsigned int sizeZ = workZ->size();

    // get vector lambda of the current interaction
    *workL = *interaction->getLambdaPtr(0);
    if (index == 0)
    {
      if (!G0Ptr)
        RuntimeException::selfThrow("LagrangianCompliantR::computeG() is not linked to a plugin function");
      G0Ptr(sizeQ, &(*workX)(0), sizeY, &(*workL)(0), &(*G[0])(0, 0), sizeZ, &(*workZ)(0));
    }
    else if (index == 1)
    {
      if (!G1Ptr)
        RuntimeException::selfThrow("LagrangianCompliantR::computeG() is not linked to a plugin function");
      G1Ptr(sizeQ, &(*workX)(0), sizeY, &(*workL)(0), &(*G[1])(0, 0), sizeZ, &(*workZ)(0));
    }
    // Copy data that might have been changed in the plug-in call.
    *data["z"] = *workZ;
  }
}

void LagrangianCompliantR::computeOutput(double time, unsigned int derivativeNumber)
{
  if (derivativeNumber == 0)
    computeH(time);
  else
  {
    SP::SiconosVector y = interaction->getYPtr(derivativeNumber);
    computeG(time, 0);
    computeG(time, 1);
    if (derivativeNumber == 1)
    {
      SP::SiconosVector lambda = interaction->getLambdaPtr(derivativeNumber);
      // y = G0 q1 + G1 lambda
      prod(*G[0], *data["q1"], *y);
      prod(*G[1], *lambda, *y, false);
    }
    else if (derivativeNumber == 2)
      prod(*G[0], *data["q2"], *y); // Approx: y[2] = G0q[2], other terms are neglected ...
    else
      RuntimeException::selfThrow("LagrangianR2::computeOutput(time,index), index out of range or not yet implemented.");
  }
}

void LagrangianCompliantR::computeInput(const double time, const unsigned int level)
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
  //  gemv(CblasTrans, 1.0,*(G[0]), *lambda, 1.0, *data[name]); => not yet implemented for BlockVectors.
}

LagrangianCompliantR* LagrangianCompliantR::convert(Relation *r)
{
  LagrangianCompliantR* lnlr = dynamic_cast<LagrangianCompliantR*>(r);
  return lnlr;
}

