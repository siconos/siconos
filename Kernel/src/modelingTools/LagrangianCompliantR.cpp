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

#include "LagrangianCompliantR.hpp"
#include "RelationXML.hpp"
#include "Interaction.hpp"
#include "LagrangianDS.hpp"
#include "Plugin.hpp"

using namespace std;
using namespace RELATION;

// xml constructor
LagrangianCompliantR::LagrangianCompliantR(SP::RelationXML LRxml): LagrangianR(LRxml, CompliantR)
{
  zeroPlugin();
  if (!LRxml->hasH())
    RuntimeException::selfThrow("LagrangianCompliantR:: xml constructor failed, can not find a definition for h.");

  setComputehFunction(SSL::getPluginName(LRxml->gethPlugin()), SSL::getPluginFunctionName(LRxml->gethPlugin()));

  if (!LRxml->hasJacobianH())
    RuntimeException::selfThrow("LagrangianCompliantR:: xml constructor failed, can not find a definition for Jach0.");
  //   Jach.resize(2);
  if (LRxml->isJacobianHPlugin(0))
    _pluginJachq->setComputeFunction(SSL::getPluginName(LRxml->getJacobianHPlugin(0)), SSL::getPluginFunctionName(LRxml->getJacobianHPlugin(0)));
  if (LRxml->isJacobianHPlugin(1))
    _pluginJachlambda->setComputeFunction(SSL::getPluginName(LRxml->getJacobianHPlugin(1)), SSL::getPluginFunctionName(LRxml->getJacobianHPlugin(1)));
}

// constructor from a set of data
LagrangianCompliantR::LagrangianCompliantR(const string& computeh, const std::vector<string> & computeg): LagrangianR(CompliantR)
{
  zeroPlugin();
  setComputehFunction(SSL::getPluginName(computeh), SSL::getPluginFunctionName(computeh));
  _pluginJachq->setComputeFunction(SSL::getPluginName(computeg[0]), SSL::getPluginFunctionName(computeg[0]));
  _pluginJachlambda->setComputeFunction(SSL::getPluginName(computeg[1]), SSL::getPluginFunctionName(computeg[1]));
}

void LagrangianCompliantR::zeroPlugin()
{
  _pluginJachq.reset(new PluggedObject());
  _pluginJachlambda.reset(new PluggedObject());
}

const std::string LagrangianCompliantR::getJachlambdaName() const
{
  if (_pluginJachlambda->fPtr)
    return _pluginJachlambda->getPluginName();
  return "unamed";

}
const std::string LagrangianCompliantR::getJachqName() const
{
  if (_pluginJachq->fPtr)
    return _pluginJachq->getPluginName();
  return "unamed";
}


void LagrangianCompliantR::initComponents(Interaction& inter)
{
  LagrangianR::initComponents(inter);
  unsigned int sizeY = inter.getSizeOfY();

  if (! _jachlambda)
    _jachlambda.reset(new SimpleMatrix(sizeY, sizeY));
  else
    _jachlambda->resize(sizeY, sizeY);
}

void LagrangianCompliantR::computeh(const double time, Interaction& inter)
{
  if (_pluginh->fPtr)
  {
    // get vector y of the current interaction
    SiconosVector& y = *inter.y(0);
    SiconosVector& lambda = *inter.lambda(0);

    // Warning: temporary method to have contiguous values in memory, copy of block to simple.
    SiconosVector workQ = *inter.data(q0);
    SiconosVector workZ = *inter.data(z);

    ((FPtr2)(_pluginh->fPtr))(workQ.size(), &(workQ)(0), y.size(), &(lambda)(0), &(y)(0), workZ.size(), &(workZ)(0));

    // Copy data that might have been changed in the plug-in call.
    *inter.data(z) = workZ;
  }
}

void LagrangianCompliantR::computeJachq(const double time, Interaction& inter)
{

  if (_pluginJachq->fPtr)
  {
    // Warning: temporary method to have contiguous values in memory, copy of block to simple.
    SiconosVector workQ = *inter.data(q0);
    SiconosVector workZ = *inter.data(z);
    SiconosVector lambda = *inter.lambda(0);

    // get vector lambda of the current interaction
    ((FPtr2)(_pluginJachq->fPtr))(workQ.size(), &(workQ)(0), lambda.size(), &(lambda)(0), &(*_jachq)(0, 0), workZ.size(), &(workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *inter.data(z) = workZ;
  }
}
void LagrangianCompliantR::computeJachlambda(const double time, Interaction& inter)
{

  if (_pluginJachlambda->fPtr)
  {
    SiconosVector workQ = *inter.data(q0);
    SiconosVector workZ = *inter.data(z);
    SiconosVector& lambda = *inter.lambda(0);

    // get vector lambda of the current interaction
    ((FPtr2)_pluginJachlambda->fPtr)(workQ.size(), &(workQ)(0), lambda.size(), &(lambda)(0), &(*_jachlambda)(0, 0), workZ.size(), &(workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *inter.data(z) = workZ;
  }
}

void LagrangianCompliantR::computeOutput(const double time, Interaction& inter, unsigned int derivativeNumber)
{
  if (derivativeNumber == 0)
    computeh(time, inter);
  else
  {
    SiconosVector& y = *inter.y(derivativeNumber);
    computeJachq(time, inter);
    computeJachlambda(time, inter);
    if (derivativeNumber == 1)
    {
      SiconosVector& lambda = *inter.lambda(derivativeNumber);
      // y = Jach[0] q1 + Jach[1] lambda
      prod(*_jachq, *inter.data(q1), y);
      prod(*_jachlambda, lambda, y, false);
    }
    else if (derivativeNumber == 2)
      prod(*_jachq, *inter.data(q2), y); // Approx: y[2] = Jach[0]q[2], other terms are neglected ...
    else
      RuntimeException::selfThrow("LagrangianCompliantR::computeOutput(time,index), index out of range or not yet implemented.");
  }
}

void LagrangianCompliantR::computeInput(const double time, Interaction& inter, const unsigned int level)
{
  computeJachq(time, inter);
  // get lambda of the concerned interaction
  SiconosVector& lambda = *inter.lambda(level);

  // data[name] += trans(G) * lambda
  prod(lambda, *_jachq, *inter.data(p0 + level), false);
}

LagrangianCompliantR* LagrangianCompliantR::convert(Relation *r)
{
  LagrangianCompliantR* lnlr = dynamic_cast<LagrangianCompliantR*>(r);
  return lnlr;
}

