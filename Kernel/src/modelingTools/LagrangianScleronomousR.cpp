/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
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

  setComputehFunction(SSL::getPluginName(LRxml->gethPlugin()), SSL::getPluginFunctionName(LRxml->gethPlugin()));

  if (!LRxml->hasJacobianH())
    RuntimeException::selfThrow("LagrangianScleronomousR:: xml constructor failed, can not find a definition for Jach0.");
  //  LRxml->readJacobianXML<PluggedMatrix,SP_PluggedMatrix>(Jach[0], LRxml, 0);
  _pluginjqh.reset(new PluggedObject());
  if (LRxml->isJacobianHPlugin(0))
  {
    _pluginjqh->setComputeFunction(LRxml->getJacobianHPlugin(0));
  }
  else
    _jachq.reset(new SimpleMatrix(LRxml->getJacobianHMatrix(0)));

}

// constructor from a set of data
LagrangianScleronomousR::LagrangianScleronomousR(const string& computeh, const std::string& strcomputeJachq):
  LagrangianR(ScleronomousR)
{
  setComputehFunction(SSL::getPluginName(computeh), SSL::getPluginFunctionName(computeh));

  _pluginjqh.reset(new PluggedObject());
  _pluginjqh->setComputeFunction(strcomputeJachq);

  //  unsigned int sizeY = interaction()->getSizeOfY();
  //  unsigned int sizeQ = workX->size();
  //  _jachq.reset(new SimpleMatrix(sizeY,sizeQ));

  // Warning: we cannot allocate memory for Jach[0] matrix since no interaction
  // is connected to the relation. This will be done during initialize.
  // We only set the name of the plugin-function and connect it to the user-defined function.
}
// constructor from a data used for EventDriven scheme
LagrangianScleronomousR::LagrangianScleronomousR(const std::string& computeh, const std::string& strcomputeJachq, const std::string& computeJachqdot):
  LagrangianR(ScleronomousR)
{
  setComputehFunction(SSL::getPluginName(computeh), SSL::getPluginFunctionName(computeh));

  _pluginjqh.reset(new PluggedObject());
  _pluginjqh->setComputeFunction(strcomputeJachq);

  _pluginjqhdot.reset(new PluggedObject());
  _pluginjqhdot->setComputeFunction(computeJachqdot);
}
void LagrangianScleronomousR::computeh(double)
{
  if (_pluginh)
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
  }
  // else nothing
}

void LagrangianScleronomousR::computeJachq(double)
{
  // First arg: time. Useless.
  // Last arg: index for G - Useless, always equal to 0 for this kind of relation.

  //
  if (_pluginjqh)
  {
    if (_pluginjqh->fPtr)
    {
      // Warning: temporary method to have contiguous values in memory, copy of block to simple.
      *_workX = *data[q0];
      *_workZ = *data[z];

      unsigned int sizeY = _jachq->size(0);
      unsigned int sizeQ = _workX->size();
      unsigned int sizeZ = _workZ->size();

      ((FPtr3)(_pluginjqh->fPtr))(sizeQ, &(*_workX)(0), sizeY, &(*_jachq)(0, 0), sizeZ, &(*_workZ)(0));

      // Copy data that might have been changed in the plug-in call.
      *data[z] = *_workZ;
    }
  }
  //  else nothing!
}

void LagrangianScleronomousR::computeJachqDot(double time)
{
  if (_pluginjqhdot)
  {
    if (_pluginjqhdot->fPtr)
    {
      // Warning: temporary method to have contiguous values in memory, copy of block to simple.
      *_workX = *data[q0];
      *_workXdot = *data[q1];
      *_workZ = *data[z];

      unsigned int sizeS = _jachqDot->size(0);
      unsigned int sizeQ = _workX->size();
      unsigned int sizeQdot = _workXdot->size();
      unsigned int sizeZ = _workZ->size();

      ((FPtr5bis)(_pluginjqhdot->fPtr))(sizeQ, &(*_workX)(0), sizeQdot, &(*_workXdot)(0) , sizeS, &(*_jachqDot)(0, 0), sizeZ, &(*_workZ)(0));

      // Copy data that might have been changed in the plug-in call.
      *data[z] = *_workZ;
    }
  }
}
void  LagrangianScleronomousR::computeNonLinearH2dot(double time)
{
  // Compute the H Jacobian dot
  LagrangianScleronomousR::computeJachqDot(time);
  //
  _NLh2dot.reset(new SiconosVector(_jachqDot->size(0)));
  prod(*_jachqDot, *_workXdot, *_NLh2dot);
}

void LagrangianScleronomousR::computeOutput(double time, unsigned int derivativeNumber)
{
  if (derivativeNumber == 0)
    computeh(time);
  else
  {
    computeJachq(time);
    SP::SiconosVector y = interaction()->y(derivativeNumber) ;
    if (derivativeNumber == 1)
      prod(*_jachq, *data[q1], *y);
    else if (derivativeNumber == 2)
      prod(*_jachq, *data[q2], *y);
    else
      RuntimeException::selfThrow("LagrangianScleronomousR::computeOutput(t,index), index out of range");
  }
}

void LagrangianScleronomousR::computeInput(double time, unsigned int level)
{
  computeJachq(time);
  // get lambda of the concerned interaction
  SP::SiconosVector lambda = interaction()->lambda(level);
  // data[name] += trans(G) * lambda
  prod(*lambda, *_jachq, *data[p0 + level], false);

}
const std::string LagrangianScleronomousR::getJachqName() const
{
  if (_pluginjqh->fPtr)
    return _pluginjqh->getPluginName();
  return "unamed";

}
LagrangianScleronomousR* LagrangianScleronomousR::convert(Relation *r)
{
  return dynamic_cast<LagrangianScleronomousR*>(r);
}
