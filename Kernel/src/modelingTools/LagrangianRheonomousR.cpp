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

#include "LagrangianRheonomousR.hpp"
#include "RelationXML.hpp"
#include "Interaction.hpp"
#include "LagrangianDS.hpp"

using namespace std;
using namespace RELATION;

// xml constructor
LagrangianRheonomousR::LagrangianRheonomousR(SP::RelationXML LRxml): LagrangianR(LRxml, RheonomousR)
{
  zeroPlugin();
  // h plug-in
  if (!LRxml->hasH())
    RuntimeException::selfThrow("LagrangianRheonomousR:: xml constructor failed, can not find a definition for h.");
  setComputehFunction(SSL::getPluginName(LRxml->gethPlugin()), SSL::getPluginFunctionName(LRxml->gethPlugin()));

  // Read hDot
  if (!LRxml->hasHDot())
    RuntimeException::selfThrow("LagrangianRheonomousR:: xml constructor failed, can not find a definition for hDot.");
  if (LRxml->isHDotPlugin())
  {
    //    hDot.reset(new PVT2(LRxml->gethDotPlugin()));
    _pluginhDot->setComputeFunction(SSL::getPluginName(LRxml->gethDotPlugin()), SSL::getPluginFunctionName(LRxml->gethDotPlugin()));
  }
  else
    _hDot.reset(new SiconosVector(LRxml->gethDotVector()));

  if (!LRxml->hasJacobianH())
    RuntimeException::selfThrow("LagrangianRheonomousR:: xml constructor failed, can not find a definition for G0.");
  if (LRxml->isJacobianHPlugin(0))
    _pluginJachq->setComputeFunction(LRxml->getJacobianHPlugin(0));
  else
    _jachq.reset(new SimpleMatrix(LRxml->getJacobianHMatrix(0)));
}

// constructor from a set of data
LagrangianRheonomousR::LagrangianRheonomousR(const string& computeh, const string& computehDot, const string& strcomputeJachq):
  LagrangianR(RheonomousR)
{
  zeroPlugin();
  // h
  setComputehFunction(SSL::getPluginName(computeh), SSL::getPluginFunctionName(computeh));

  // hDot
  setComputehDotFunction(SSL::getPluginName(computehDot), SSL::getPluginFunctionName(computehDot));
  _pluginJachq->setComputeFunction(strcomputeJachq);

}

void LagrangianRheonomousR::initComponents(Interaction& inter)
{
  LagrangianR::initComponents(inter);

  unsigned int sizeY = inter.getSizeOfY();
  // hDot
  if (!_hDot)
    _hDot.reset(new SiconosVector(sizeY));
  else
    _hDot->resize(sizeY);
  if (_pluginJachq->fPtr && !_jachq)
  {
    unsigned int sizeY = inter.getSizeOfY();
    unsigned int sizeQ = inter.data(q0)->size();
    _jachq.reset(new SimpleMatrix(sizeY, sizeQ));
  }
}
// void LagrangianRheonomousR::setComputehFunction(const string& pluginPath, const string& functionName){
//   Plugin::setFunction(&hPtr, pluginPath, functionName);
// }

void LagrangianRheonomousR::setComputehDotFunction(const string& pluginPath, const string& functionName)
{
  _pluginhDot->setComputeFunction(pluginPath, functionName);
}
void LagrangianRheonomousR::zeroPlugin()
{
  LagrangianR::zeroPlugin();
  _pluginhDot.reset(new PluggedObject());
}
const std::string LagrangianRheonomousR::getJachqName() const
{
  if (_pluginJachq->fPtr)
    return _pluginJachq->getPluginName();
  return "unamed";
}
const std::string LagrangianRheonomousR::gethDotName() const
{

  if (_pluginhDot->fPtr)
    return _pluginhDot->getPluginName();
  return "unamed";
}
void LagrangianRheonomousR::computeh(const double time, Interaction& inter)
{
  if (_pluginh->fPtr)
  {
    // get vector y of the current interaction
    SiconosVector& y = *inter.y(0);

    // Warning: temporary method to have contiguous values in
    // memory, copy of block to simple.
    SiconosVector workQ = *inter.data(q0);
    SiconosVector workZ = *inter.data(z);

    ((FPtr4)(_pluginh->fPtr))(workQ.size(), &(workQ)(0), time, y.size(),  &(y)(0), workZ.size(), &(workZ)(0));

    // Copy data that might have been changed in the plug-in call.
    *inter.data(z) = workZ;
  }
  // else nothing
}

void LagrangianRheonomousR::computehDot(const double time, Interaction& inter)
{
  if (_pluginhDot->fPtr)
  {
    // Warning: temporary method to have contiguous values in
    // memory, copy of block to simple.
    SiconosVector workQ = *inter.data(q0);
    SiconosVector workZ = *inter.data(z);

    ((FPtr4)(_pluginhDot->fPtr))(workQ.size(), &(workQ)(0), time, _hDot->size(),  &(*_hDot)(0), workZ.size(), &(workZ)(0));

    // Copy data that might have been changed in the plug-in call.
    *inter.data(z) = workZ;
  }
  // else nothing
}

void LagrangianRheonomousR::computeJachq(const double time, Interaction& inter)
{
  // Note that second input arg is useless.
  if (_pluginJachq->fPtr)
  {
    // Warning: temporary method to have contiguous values in
    // memory, copy of block to simple.
    SiconosVector workQ = *inter.data(q0);
    SiconosVector workZ = *inter.data(z);

    ((FPtr4)(_pluginJachq->fPtr))(workQ.size(), &(workQ)(0), time, _jachq->size(0), &(*_jachq)(0, 0), workZ.size(), &(workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *inter.data(z) = workZ;
  }
  // else nothing.
}

void LagrangianRheonomousR::computeOutput(const double time, Interaction& inter, unsigned int derivativeNumber)
{
  if (derivativeNumber == 0)
    computeh(time, inter);
  else
  {
    SiconosVector& y = *inter.y(derivativeNumber);
    computeJachq(time, inter);
    if (derivativeNumber == 1)
    {
      // Computation of the partial derivative w.r.t time of h(q,t)
      computehDot(time, inter);
      // Computation of the partial derivative w.r.t q of h(q,t) : \nabla_q h(q,t) \dot q
      prod(*_jachq, *inter.data(q1), y);
      // Sum of the terms
      y += *_hDot;
    }
    else if (derivativeNumber == 2)
      prod(*_jachq, *inter.data(q2), y); // Approx:,  ...
    // \warning : the computation of y[2] (in event-driven
    // simulation for instance) is approximated by y[2] =
    // Jach[0]q[2]. For the moment, other terms are neglected
    // (especially, partial derivatives with respect to time).
    else
      RuntimeException::selfThrow("LagrangianRheonomousR::computeOutput(time,index), index >2  not yet implemented.");
  }
}

void LagrangianRheonomousR::computeInput(const double time, Interaction& inter, unsigned int level)
{
  computeJachq(time, inter);
  // get lambda of the concerned interaction
  SiconosVector& lambda = *inter.lambda(level);
  // data[name] += trans(G) * lambda
  prod(lambda, *_jachq, *inter.data(p0 + level), false);
}

LagrangianRheonomousR* LagrangianRheonomousR::convert(Relation *r)
{
  return dynamic_cast<LagrangianRheonomousR*>(r);
}

