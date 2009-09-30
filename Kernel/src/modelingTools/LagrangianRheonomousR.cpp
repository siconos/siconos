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
#include "RelationXML.h"
#include "Interaction.h"
#include "LagrangianDS.h"

using namespace std;
using namespace RELATION;

// xml constructor
LagrangianRheonomousR::LagrangianRheonomousR(SP::RelationXML LRxml): LagrangianR(LRxml, RheonomousR)
{
  /*  // h plug-in
  if(!LRxml->hasH())
    RuntimeException::selfThrow("LagrangianRheonomousR:: xml constructor failed, can not find a definition for h.");
  hName = LRxml->getHPlugin();
  setComputeHFunction(SSL::getPluginName( hName ),SSL::getPluginFunctionName( hName ));

  // Read hDot
  if(!LRxml->hasHDot())
    RuntimeException::selfThrow("LagrangianRheonomousR:: xml constructor failed, can not find a definition for hDot.");
  if(LRxml->isHDotPlugin())
    hDot.reset(new PVT2(LRxml->getHDotPlugin()));
  else
    hDot.reset(new PVT2(LRxml->getHDotVector()));

  if(!LRxml->hasJacobianH())
    RuntimeException::selfThrow("LagrangianRheonomousR:: xml constructor failed, can not find a definition for G0.");
  JacH.resize(1);
  LRxml->readJacobianXML<PluggedMatrix, SP_PluggedMatrix>(JacH[0], LRxml, 0);*/
}

// constructor from a set of data
LagrangianRheonomousR::LagrangianRheonomousR(const string& computeH, const string& computeHDot, const string& strcomputeJacQH):
  LagrangianR(RheonomousR)
{
  // h
  setComputeHFunction(SSL::getPluginName(computeH), SSL::getPluginFunctionName(computeH));

  // hDot
  setComputeHDotFunction(SSL::getPluginName(computeHDot), SSL::getPluginFunctionName(computeHDot));
  Plugin::setFunction(&computeJacQHPtr, SSL::getPluginName(strcomputeJacQH), SSL::getPluginFunctionName(strcomputeJacQH));

  unsigned int sizeY = getInteractionPtr()->getSizeOfY();
  unsigned int sizeQ = workX->size();
  JacQH.reset(new SimpleMatrix(sizeY, sizeQ));
}

void LagrangianRheonomousR::initComponents()
{
  LagrangianR::initComponents();

  unsigned int sizeY = getInteractionPtr()->getSizeOfY();
  // hDot
  if (!hDot)
    hDot.reset(new SimpleVector(sizeY));
  else
    hDot->resize(sizeY);
}
void LagrangianRheonomousR::setComputeHFunction(const string& pluginPath, const string& functionName)
{
  Plugin::setFunction(&hPtr, pluginPath, functionName);
}

void LagrangianRheonomousR::setComputeHDotFunction(const string& pluginPath, const string& functionName)
{
  Plugin::setFunction(&hDotPtr, pluginPath, functionName);
}

void LagrangianRheonomousR::computeH(double time)
{
  if (hPtr)
  {
    // get vector y of the current interaction
    SP::SiconosVector y = getInteractionPtr()->getYPtr(0);

    // Warning: temporary method to have contiguous values in
    // memory, copy of block to simple.
    *workX = *data[q0];
    *workZ = *data[z];
    *workY = *y;

    unsigned int sizeQ = workX->size();
    unsigned int sizeY = y->size();
    unsigned int sizeZ = workZ->size();

    hPtr(sizeQ, &(*workX)(0), time, sizeY,  &(*workY)(0), sizeZ, &(*workZ)(0));

    // Copy data that might have been changed in the plug-in call.
    *data[z] = *workZ;
    *y = *workY;
  }
  // else nothing
}
void LagrangianRheonomousR::computeG(double time)
{
  assert(false && "LagrangianScleronomousR::computeG : G is computed in computeInput!\n");
}

void LagrangianRheonomousR::computeHDot(double time)
{
  if (hDotPtr)
  {
    // Warning: temporary method to have contiguous values in
    // memory, copy of block to simple.
    *workX = *data[q0];
    *workZ = *data[z];

    unsigned int sizeQ = workX->size();
    unsigned int sizeY = hDot->size();
    unsigned int sizeZ = workZ->size();

    hDotPtr(sizeQ, &(*workX)(0), time, sizeY,  &(*hDot)(0), sizeZ, &(*workZ)(0));

    // Copy data that might have been changed in the plug-in call.
    *data[z] = *workZ;
  }
  // else nothing
}

void LagrangianRheonomousR::computeJacQH(double time)
{
  // Note that second input arg is useless.
  if (computeJacQHPtr)
  {
    // Warning: temporary method to have contiguous values in
    // memory, copy of block to simple.
    *workX = *data[q0];
    *workZ = *data[z];

    unsigned int sizeY = JacQH->size(0);
    unsigned int sizeQ = workX->size();
    unsigned int sizeZ = workZ->size();
    ((FPtr4) computeJacQHPtr)(sizeQ, &(*workX)(0), time, sizeY, &(*JacQH)(0, 0), sizeZ, &(*workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *data[z] = *workZ;
  }
  // else nothing.
}

void LagrangianRheonomousR::computeOutput(double time, unsigned int derivativeNumber)
{
  if (derivativeNumber == 0)
    computeH(time);
  else
  {
    SP::SiconosVector y = getInteractionPtr()->getYPtr(derivativeNumber);
    computeJacQH(time);
    if (derivativeNumber == 1)
    {
      computeHDot(time); // \todo: save hDot directly into y[1] ?
      prod(*JacQH, *data[q1], *y);
      *y += *hDot;
    }
    else if (derivativeNumber == 2)
      prod(*JacQH, *data[q2], *y); // Approx:,  ...
    // \warning :  the computation of y[2] (in event-driven simulation for instance) is approximated by  y[2] = JacH[0]q[2]. For the moment, other terms are neglected (especially, partial derivatives with respect to time).
    else
      RuntimeException::selfThrow("LagrangianRheonomousR::computeOutput(time,index), index out of range or not yet implemented.");
  }
}

void LagrangianRheonomousR::computeInput(double time, unsigned int level)
{
  computeJacQH(time);
  // get lambda of the concerned interaction
  SP::SiconosVector lambda = getInteractionPtr()->getLambdaPtr(level);
  // data[name] += trans(G) * lambda
  prod(*lambda, *JacQH, *data[p0 + level], false);

  //   SP::SiconosMatrix  GT = new SimpleMatrix(*G[0]);
  //   GT->trans();
  //   *data[name] += prod(*GT, *lambda);
  //  gemv(CblasTrans, 1.0,*(G[0]), *lambda, 1.0, *data[name]);
}

LagrangianRheonomousR* LagrangianRheonomousR::convert(Relation *r)
{
  return dynamic_cast<LagrangianRheonomousR*>(r);
}

