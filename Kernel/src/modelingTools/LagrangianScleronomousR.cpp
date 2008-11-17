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
#include "RelationXML.h"
#include "Interaction.h"
#include "LagrangianDS.h"
#include "LagrangianR.cpp"

using namespace std;
using namespace RELATION;

// xml constructor
LagrangianScleronomousR::LagrangianScleronomousR(SP::RelationXML LRxml): LagrangianR<FPtr3>(LRxml, ScleronomousR)
{
  // h plug-in
  if (!LRxml->hasH())
    RuntimeException::selfThrow("LagrangianScleronomousR:: xml constructor failed, can not find a definition for h.");

  hName = LRxml->getHPlugin();
  setComputeHFunction(SSL::getPluginName(hName), SSL::getPluginFunctionName(hName));

  if (!LRxml->hasJacobianH())
    RuntimeException::selfThrow("LagrangianScleronomousR:: xml constructor failed, can not find a definition for JacH0.");
  JacH.resize(1);
  LRxml->readJacobianXML<PluggedMatrix, SP_PluggedMatrix>(JacH[0], LRxml, 0);
}

// constructor from a set of data
LagrangianScleronomousR::LagrangianScleronomousR(const string& computeH, const std::string& computeG):
  LagrangianR<FPtr3>(ScleronomousR)
{
  setComputeHFunction(SSL::getPluginName(computeH), SSL::getPluginFunctionName(computeH));

  // Warning: we cannot allocate memory for JacH[0] matrix since no interaction
  // is connected to the relation. This will be done during initialize.
  // We only set the name of the plugin-function and connect it to the user-defined function.
  JacH.resize(1);
  JacH[0].reset(new PluggedMatrix(computeG));
}

void LagrangianScleronomousR::computeH(double)
{
  // arg= time. Unused in this function but required for interface.
  if (hPlugged)
  {
    // get vector y of the current interaction
    SP::SiconosVector y = getInteractionPtr()->getYPtr(0);

    // Warning: temporary method to have contiguous values in memory, copy of block to simple.
    *workX = *data[q0];
    *workZ = *data[z];
    *workY = *y;

    unsigned int sizeQ = workX->size();
    unsigned int sizeY = y->size();
    unsigned int sizeZ = workZ->size();

    if (!hPtr)
      RuntimeException::selfThrow("LagrangianScleronomousR:computeH() failed, h is not linked to a plugin function");
    hPtr(sizeQ, &(*workX)(0) , sizeY, &(*workY)(0), sizeZ, &(*workZ)(0));

    // Copy data that might have been changed in the plug-in call.
    *data[z] = *workZ;
    *y = *workY;
  }
  // else nothing
}

void LagrangianScleronomousR::computeJacH(double, unsigned int index)
{
  // First arg: time. Useless.
  // Last arg: index for G - Useless, always equal to 0 for this kind of relation.

  //
  assert(index == 0 && "LagrangianScleronomousR::computeJacH(index): index is out of range");
  if (JacH[0]->isPlugged())
  {
    // Warning: temporary method to have contiguous values in memory, copy of block to simple.
    *workX = *data[q0];
    *workZ = *data[z];

    unsigned int sizeY = JacH[0]->size(0);
    unsigned int sizeQ = workX->size();
    unsigned int sizeZ = workZ->size();

    if (!(JacH[0]->fPtr))
      RuntimeException::selfThrow("LagrangianScleronomousR::computeJacH(), JacH[0] is not linked to a plugin function");
    (JacH[0]->fPtr)(sizeQ, &(*workX)(0), sizeY, &(*(JacH[0]))(0, 0), sizeZ, &(*workZ)(0));

    // Copy data that might have been changed in the plug-in call.
    *data[z] = *workZ;
  }
  //  else nothing!
}

void LagrangianScleronomousR::computeOutput(double time, unsigned int derivativeNumber)
{
  if (derivativeNumber == 0)
    computeH(time);
  else
  {
    computeJacH(time, 0);
    SP::SiconosVector y = getInteractionPtr()->getYPtr(derivativeNumber) ;
    // Approx: y[i] = jacH[0] q[i] , other terms are neglected.
    prod(*JacH[0], *data[q0 + derivativeNumber], *y);
  }
}

void LagrangianScleronomousR::computeInput(double time, unsigned int level)
{
  computeJacH(time, 0);
  // get lambda of the concerned interaction
  SP::SiconosVector lambda = getInteractionPtr()->getLambdaPtr(level);
  // data[name] += trans(G) * lambda
  prod(*lambda, *JacH[0], *data[p0 + level], false);

  //   SP::SiconosMatrix  GT = new SimpleMatrix(*G[0]);
  //   GT->trans();
  //   *data[name] += prod(*GT, *lambda);
  //gemv(CblasTrans, 1.0,*(G[0]), *lambda, 1.0, *data[name]); => not yet implemented for BlockVectors.
}

LagrangianScleronomousR* LagrangianScleronomousR::convert(Relation *r)
{
  return dynamic_cast<LagrangianScleronomousR*>(r);
}
