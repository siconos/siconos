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
#include "LagrangianR.cpp"

using namespace std;
using namespace RELATION;

// xml constructor
LagrangianRheonomousR::LagrangianRheonomousR(SP::RelationXML LRxml): BaseClass(LRxml, RheonomousR)
{
  // h plug-in
  if (!LRxml->hasH())
    RuntimeException::selfThrow("LagrangianRheonomousR:: xml constructor failed, can not find a definition for h.");
  hName = LRxml->getHPlugin();
  setComputeHFunction(SSL::getPluginName(hName), SSL::getPluginFunctionName(hName));

  // Read hDot
  if (!LRxml->hasHDot())
    RuntimeException::selfThrow("LagrangianRheonomousR:: xml constructor failed, can not find a definition for hDot.");
  if (LRxml->isHDotPlugin())
    hDot.reset(new PVT2(LRxml->getHDotPlugin()));
  else
    hDot.reset(new PVT2(LRxml->getHDotVector()));

  if (!LRxml->hasJacobianH())
    RuntimeException::selfThrow("LagrangianRheonomousR:: xml constructor failed, can not find a definition for G0.");
  JacH.resize(1);
  LRxml->readJacobianXML<PluggedMatrix, SP_PluggedMatrix>(JacH[0], LRxml, 0);
}

// constructor from a set of data
LagrangianRheonomousR::LagrangianRheonomousR(const string& computeH, const string& computeHDot, const string& computeG):
  BaseClass(RheonomousR)
{
  // h
  setComputeHFunction(SSL::getPluginName(computeH), SSL::getPluginFunctionName(computeH));

  // hDot
  hDot.reset(new PVT2(computeHDot));

  JacH.resize(1);
  JacH[0].reset(new PluggedMatrix(computeG));
}

void LagrangianRheonomousR::initComponents()
{
  LagrangianR<FPtr4>::initComponents();

  unsigned int sizeY = interaction->getSizeOfY();
  // hDot
  if (!hDot)
    hDot.reset(new PVT2(sizeY));
  else
    hDot->resize(sizeY);
}

void LagrangianRheonomousR::setComputeHDotFunction(const string& pluginPath, const string& functionName)
{
  hDot->setComputeFunction(pluginPath, functionName);
}

void LagrangianRheonomousR::computeH(double time)
{
  if (hPlugged)
  {
    // get vector y of the current interaction
    SP::SiconosVector y = interaction->getYPtr(0);

    // Warning: temporary method to have contiguous values in memory, copy of block to simple.
    *workX = *data[q0];
    *workZ = *data[z];
    *workY = *y;

    unsigned int sizeQ = workX->size();
    unsigned int sizeY = y->size();
    unsigned int sizeZ = workZ->size();

    if (!hPtr)
      RuntimeException::selfThrow("LagrangianRheonomousR:computeH() failed, h is not linked to a plugin function");
    hPtr(sizeQ, &(*workX)(0), time, sizeY,  &(*workY)(0), sizeZ, &(*workZ)(0));

    // Copy data that might have been changed in the plug-in call.
    *data[z] = *workZ;
    *y = *workY;
  }
  // else nothing
}

void LagrangianRheonomousR::computeHDot(double time)
{
  if (hDot->isPlugged())
  {
    // Warning: temporary method to have contiguous values in memory, copy of block to simple.
    *workX = *data[q0];
    *workZ = *data[z];

    unsigned int sizeQ = workX->size();
    unsigned int sizeY = hDot->size();
    unsigned int sizeZ = workZ->size();

    if (!(hDot->fPtr))
      RuntimeException::selfThrow("LagrangianRheonomousR:computeHDot() failed, hDot is not linked to a plugin function");
    (hDot->fPtr)(sizeQ, &(*workX)(0), time, sizeY,  &(*hDot)(0), sizeZ, &(*workZ)(0));

    // Copy data that might have been changed in the plug-in call.
    *data[z] = *workZ;
  }
  // else nothing
}

void LagrangianRheonomousR::computeJacH(double time, unsigned int)
{
  // Note that second input arg is useless.
  assert(index == 0 && "LagrangianRheonomousR::computeJacH(index): index is out of range");
  if (JacH[0]->isPlugged())
  {
    // Warning: temporary method to have contiguous values in memory, copy of block to simple.
    *workX = *data[q0];
    *workZ = *data[z];

    unsigned int sizeY = JacH[0]->size(0);
    unsigned int sizeQ = workX->size();
    unsigned int sizeZ = workZ->size();
    if (!(JacH[0]->fPtr))
      RuntimeException::selfThrow("LagrangianRheonomousR::computeJacH(), JacH[0] is not linked to a plugin function");
    (JacH[0]->fPtr)(sizeQ, &(*workX)(0), time, sizeY, &(*JacH[0])(0, 0), sizeZ, &(*workZ)(0));
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
    SP::SiconosVector y = interaction->getYPtr(derivativeNumber);
    computeJacH(time, 0);
    if (derivativeNumber == 1)
    {
      computeHDot(time); // \todo: save hDot directly into y[1] ?
      prod(*JacH[0], *data[q1], *y);
      *y += *hDot;
    }
    else if (derivativeNumber == 2)
      prod(*JacH[0], *data[q2], *y); // Approx: y[2] = JacH[0]q[2], other terms are neglected ...
    else
      RuntimeException::selfThrow("LagrangianRheonomousR::computeOutput(time,index), index out of range or not yet implemented.");
  }
}

void LagrangianRheonomousR::computeInput(double time, unsigned int level)
{
  computeJacH(time, 0);
  // get lambda of the concerned interaction
  SP::SiconosVector lambda = interaction->getLambdaPtr(level);
  // data[name] += trans(G) * lambda
  prod(*lambda, *JacH[0], *data[p0 + level], false);
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

