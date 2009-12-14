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

// \todo : create a work vector for all tmp vectors used in computeg, computeh ...

#include "LagrangianScleronomousR.hpp"
#include "RelationXML.hpp"
#include "Interaction.hpp"
#include "LagrangianDS.hpp"

using namespace std;
using namespace RELATION;

// xml constructor
LagrangianScleronomousR::LagrangianScleronomousR(SP::RelationXML LRxml): LagrangianR(LRxml, ScleronomousR)
{
  // h plug-in
  if (!LRxml->hasH())
    RuntimeException::selfThrow("LagrangianScleronomousR:: xml constructor failed, can not find a definition for h.");

  hName = LRxml->gethPlugin();
  setComputehFunction(SSL::getPluginName(hName), SSL::getPluginFunctionName(hName));

  if (!LRxml->hasJacobianH())
    RuntimeException::selfThrow("LagrangianScleronomousR:: xml constructor failed, can not find a definition for JacH0.");
  //  LRxml->readJacobianXML<PluggedMatrix,SP_PluggedMatrix>(JacH[0], LRxml, 0);
  if (LRxml->isJacobianHPlugin(0))
    pluginjqh->setComputeFunction(LRxml->getJacobianHPlugin(0));
  else
    Jacqh.reset(new SimpleMatrix(LRxml->getJacobianHMatrix(0)));

}

// constructor from a set of data
LagrangianScleronomousR::LagrangianScleronomousR(const string& computeh, const std::string& strcomputeJacqh):
  LagrangianR(ScleronomousR)
{
  setComputehFunction(SSL::getPluginName(computeh), SSL::getPluginFunctionName(computeh));

  pluginjqh.reset(new PluggedObject());
  pluginjqh->setComputeFunction(strcomputeJacqh);

  //  unsigned int sizeY = interaction()->getSizeOfY();
  //  unsigned int sizeQ = workX->size();
  //  Jacqh.reset(new SimpleMatrix(sizeY,sizeQ));

  // Warning: we cannot allocate memory for JacH[0] matrix since no interaction
  // is connected to the relation. This will be done during initialize.
  // We only set the name of the plugin-function and connect it to the user-defined function.
}

void LagrangianScleronomousR::computeh(double)
{
  // arg= time. Unused in this function but required for interface.
  if (_pluginh->fPtr)
  {
    // get vector y of the current interaction
    SP::SiconosVector y = interaction()->y(0);

    // Warning: temporary method to have contiguous values in memory, copy of block to simple.
    *_workX = *data[q0];
    *_workZ = *data[z];
    *_workY = *y;

    unsigned int sizeQ = _workX->size();
    unsigned int sizeY = y->size();
    unsigned int sizeZ = _workZ->size();

    ((FPtr3)(_pluginh->fPtr))(sizeQ, &(*_workX)(0) , sizeY, &(*_workY)(0), sizeZ, &(*_workZ)(0));

    // Copy data that might have been changed in the plug-in call.
    *data[z] = *_workZ;
    *y = *_workY;
  }
  // else nothing
}
void LagrangianScleronomousR::computeg(double, unsigned int)
{
  assert(false && "LagrangianScleronomousR::computeg : G is computed in computeInput!\n");
}
void LagrangianScleronomousR::computeJacqh(double)
{
  // First arg: time. Useless.
  // Last arg: index for G - Useless, always equal to 0 for this kind of relation.

  //
  if (pluginjqh->fPtr)
  {
    // Warning: temporary method to have contiguous values in memory, copy of block to simple.
    *_workX = *data[q0];
    *_workZ = *data[z];

    unsigned int sizeY = Jacqh->size(0);
    unsigned int sizeQ = _workX->size();
    unsigned int sizeZ = _workZ->size();

    ((FPtr3)(pluginjqh->fPtr))(sizeQ, &(*_workX)(0), sizeY, &(*Jacqh)(0, 0), sizeZ, &(*_workZ)(0));

    // Copy data that might have been changed in the plug-in call.
    *data[z] = *_workZ;
  }
  //  else nothing!
}

void LagrangianScleronomousR::computeOutput(double time, unsigned int derivativeNumber)
{
  if (derivativeNumber == 0)
    computeh(time);
  else
  {
    computeJacqh(time);
    SP::SiconosVector y = interaction()->y(derivativeNumber) ;
    if (derivativeNumber == 1)
      prod(*Jacqh, *data[q1], *y);
    else if (derivativeNumber == 2)
      prod(*Jacqh, *data[q2], *y);
    else
      RuntimeException::selfThrow("LagrangianScleronomousR::computeOutput(t,index), index out of range");
  }
}

void LagrangianScleronomousR::computeInput(double time, unsigned int level)
{
  computeJacqh(time);
  // get lambda of the concerned interaction
  SP::SiconosVector lambda = interaction()->lambda(level);
  // data[name] += trans(G) * lambda
  prod(*lambda, *Jacqh, *data[p0 + level], false);

}

LagrangianScleronomousR* LagrangianScleronomousR::convert(Relation *r)
{
  return dynamic_cast<LagrangianScleronomousR*>(r);
}
