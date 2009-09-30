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

using namespace std;
using namespace RELATION;

// xml constructor
LagrangianScleronomousR::LagrangianScleronomousR(SP::RelationXML LRxml): LagrangianR(LRxml, ScleronomousR)
{
  /* // h plug-in
  if(!LRxml->hasH())
    RuntimeException::selfThrow("LagrangianScleronomousR:: xml constructor failed, can not find a definition for h.");

  hName = LRxml->getHPlugin();
  setComputeHFunction(SSL::getPluginName( hName ),SSL::getPluginFunctionName( hName ));

  if(!LRxml->hasJacobianH())
    RuntimeException::selfThrow("LagrangianScleronomousR:: xml constructor failed, can not find a definition for JacH0.");
  JacH.resize(1);
  LRxml->readJacobianXML<PluggedMatrix,SP_PluggedMatrix>(JacH[0], LRxml, 0);
  */
}

// constructor from a set of data
LagrangianScleronomousR::LagrangianScleronomousR(const string& computeH, const std::string& strcomputeJacQH):
  LagrangianR(ScleronomousR)
{
  Plugin::setFunction(&hPtr, SSL::getPluginName(computeH), SSL::getPluginFunctionName(computeH));
  Plugin::setFunction(&computeJacQHPtr, SSL::getPluginName(strcomputeJacQH), SSL::getPluginFunctionName(strcomputeJacQH));

  // Warning: we cannot allocate memory for JacH[0] matrix since no interaction
  // is connected to the relation. This will be done during initialize.
  // We only set the name of the plugin-function and connect it to the user-defined function.
}

void LagrangianScleronomousR::computeH(double)
{
  // arg= time. Unused in this function but required for interface.
  if (hPtr)
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

    hPtr(sizeQ, &(*workX)(0) , sizeY, &(*workY)(0), sizeZ, &(*workZ)(0));

    // Copy data that might have been changed in the plug-in call.
    *data[z] = *workZ;
    *y = *workY;
  }
  // else nothing
}
void LagrangianScleronomousR::computeG(double, unsigned int)
{
  assert(false && "LagrangianScleronomousR::computeG : G is computed in computeInput!\n");
}
void LagrangianScleronomousR::computeJacQH(double)
{
  // First arg: time. Useless.
  // Last arg: index for G - Useless, always equal to 0 for this kind of relation.

  //
  if (computeJacQHPtr)
  {
    // Warning: temporary method to have contiguous values in memory, copy of block to simple.
    *workX = *data[q0];
    *workZ = *data[z];

    unsigned int sizeY = JacQH->size(0);
    unsigned int sizeQ = workX->size();
    unsigned int sizeZ = workZ->size();

    (computeJacQHPtr)(sizeQ, &(*workX)(0), sizeY, &(*JacQH)(0, 0), sizeZ, &(*workZ)(0));

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
    computeJacQH(time);
    SP::SiconosVector y = getInteractionPtr()->getYPtr(derivativeNumber) ;
    if (derivativeNumber == 1)
      prod(*JacQH, *data[q1], *y);
    else if (derivativeNumber == 2)
      prod(*JacQH, *data[q2], *y);
    else
      RuntimeException::selfThrow("LagrangianScleronomousR::computeOutput(t,index), index out of range");
  }
}

void LagrangianScleronomousR::computeInput(double time, unsigned int level)
{
  computeJacQH(time);
  // get lambda of the concerned interaction
  SP::SiconosVector lambda = getInteractionPtr()->getLambdaPtr(level);
  // data[name] += trans(G) * lambda
  prod(*lambda, *JacQH, *data[p0 + level], false);

}

LagrangianScleronomousR* LagrangianScleronomousR::convert(Relation *r)
{
  return dynamic_cast<LagrangianScleronomousR*>(r);
}
