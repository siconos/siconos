/* Siconos-Kernel, Copyright INRIA 2005-2012.
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

//#define DEBUG_MESSAGES
//#define DEBUG_STDOUT
#include "debug.h"



using namespace std;
using namespace RELATION;

// xml constructor
LagrangianScleronomousR::LagrangianScleronomousR(SP::RelationXML LRxml): LagrangianR(LRxml, ScleronomousR)
{
  zeroPlugin();
  // h plug-in
  if (!LRxml->hasH())
    RuntimeException::selfThrow("LagrangianScleronomousR:: xml constructor failed, can not find a definition for h.");

  setComputehFunction(SSLH::getPluginName(LRxml->gethPlugin()), SSLH::getPluginFunctionName(LRxml->gethPlugin()));

  if (!LRxml->hasJacobianH())
    RuntimeException::selfThrow("LagrangianScleronomousR:: xml constructor failed, can not find a definition for Jach0.");
  //  LRxml->readJacobianXML<PluggedMatrix,SP_PluggedMatrix>(Jach[0], LRxml, 0);
  if (LRxml->isJacobianHPlugin(0))
  {
    _pluginJachq->setComputeFunction(LRxml->getJacobianHPlugin(0));
  }
  else
    _jachq.reset(new SimpleMatrix(LRxml->getJacobianHMatrix(0)));

}

// constructor from a set of data
LagrangianScleronomousR::LagrangianScleronomousR(const string& computeh, const std::string& strcomputeJachq):
  LagrangianR(ScleronomousR)
{
  zeroPlugin();
  setComputehFunction(SSLH::getPluginName(computeh), SSLH::getPluginFunctionName(computeh));

  _pluginJachq->setComputeFunction(strcomputeJachq);

  //  unsigned int sizeY = inter.getSizeOfY();
  //  unsigned int sizeQ = workX->size();
  //  _jachq.reset(new SimpleMatrix(sizeY,sizeQ));

  // Warning: we cannot allocate memory for Jach[0] matrix since no interaction
  // is connected to the relation. This will be done during initialize.
  // We only set the name of the plugin-function and connect it to the user-defined function.
}
// constructor from a data used for EventDriven scheme
LagrangianScleronomousR::LagrangianScleronomousR(const std::string& computeh, const std::string& strcomputeJachq, const std::string& computeDotJachq):
  LagrangianR(ScleronomousR)
{
  zeroPlugin();
  setComputehFunction(SSLH::getPluginName(computeh), SSLH::getPluginFunctionName(computeh));

  _pluginJachq->setComputeFunction(strcomputeJachq);

  _plugindotjacqh->setComputeFunction(computeDotJachq);
}



void LagrangianScleronomousR::zeroPlugin()
{
  LagrangianR::zeroPlugin();
  _pluginJachq.reset(new PluggedObject());
  _plugindotjacqh.reset(new PluggedObject());
}

void LagrangianScleronomousR::computeh(const double time, Interaction& inter)
{

  DEBUG_PRINT("LagrangianScleronomousR::computeh(const double time, Interaction& inter)\n");
  computeh(inter, inter.data(q0),inter.data(z));

  // if (_pluginh)
  // {
  //   // arg= time. Unused in this function but required for interface.
  //   if (_pluginh->fPtr)
  //   {
  //     // get vector y of the current interaction
  //     SiconosVector& y = *inter.y(0);

  //     // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  //     SiconosVector workQ = *inter.data(q0);
  //     SiconosVector workZ = *inter.data(z);

  //     ((FPtr3)(_pluginh->fPtr))(workQ.size(), &(workQ(0)) , y.size(), &(y(0)), workZ.size(), &(workZ(0)));

  //     // Copy data that might have been changed in the plug-in call.
  //     *inter.data(z) = workZ;
  //   }
  // }
  // else nothing
}

void LagrangianScleronomousR::computeh(Interaction& inter, SP::BlockVector q, SP::BlockVector z )
{
  DEBUG_PRINT(" LagrangianScleronomousR::computeh(Interaction& inter, SP::BlockVector q, SP::BlockVector z)\n");
  if (_pluginh)
  {
    // arg= time. Unused in this function but required for interface.
    if (_pluginh->fPtr)
    {
      // get vector y of the current interaction
      SiconosVector& y = *inter.y(0);

      // Warning: temporary method to have contiguous values in memory, copy of block to simple.
      SiconosVector workQ = *q;
      SiconosVector workZ = *z;

      ((FPtr3)(_pluginh->fPtr))(workQ.size(), &(workQ(0)) , y.size(), &(y(0)), workZ.size(), &(workZ(0)));

      // Copy data that might have been changed in the plug-in call.
      *z = workZ;

      DEBUG_EXPR(q->display());
      DEBUG_EXPR(z->display());
      DEBUG_EXPR(y.display());




    }
  }
  // else nothing
}
void LagrangianScleronomousR::computeJachq(const double time, Interaction& inter)
{
  if (_pluginJachq)
  {
    if (_pluginJachq->fPtr)
    {
      // Warning: temporary method to have contiguous values in memory, copy of block to simple.
      SiconosVector workQ = *inter.data(q0);
      SiconosVector workZ = *inter.data(z);

      // get vector lambda of the current interaction
      ((FPtr3)(_pluginJachq->fPtr))(workQ.size(), &(workQ)(0), _jachq->size(0), &(*_jachq)(0, 0), workZ.size(), &(workZ)(0));
      // Copy data that might have been changed in the plug-in call.
      *inter.data(z) = workZ;
    }
  }
}

void LagrangianScleronomousR::computeDotJachq(const double time, Interaction& inter)
{
  if (_plugindotjacqh)
  {
    if (_plugindotjacqh->fPtr)
    {
      // Warning: temporary method to have contiguous values in memory, copy of block to simple.
      SiconosVector workQ = *inter.data(q0);
      SiconosVector workZ = *inter.data(z);
      SiconosVector workQdot = *inter.data(q1);


      // get vector _jachqDo of the current interaction
      if (! _jachqDot)
      {
        unsigned int sizeY = inter.getSizeOfY();
        unsigned int sizeDS = inter.getSizeOfDS();
        _dotjachq.reset(new SimpleMatrix(sizeY, sizeDS));
      }
      // else
      // {
      //   if (_dotjachq->size(0) == 0) // if the matrix dimension are null
      //     _dotjachq->resize(sizeY, sizeDS);
      //   else
      //   {
      //     if ((_dotjachq->size(1) != sizeDS && _dotjachq->size(0) != sizeY))
      //       RuntimeException::selfThrow("LagrangianR::initComponents inconsistent sizes between Jach[1] matrix and the interaction.");
      //   }
      // }

      ((FPtr2)(_plugindotjacqh->fPtr))(workQ.size(), &(workQ)(0), workQdot.size(), &(workQdot)(0), &(*_dotjachq)(0, 0), workZ.size(), &(workZ)(0));
      // Copy data that might have been changed in the plug-in call.
      *inter.data(z) = workZ;
    }
  }
}

void  LagrangianScleronomousR::computedotjacqhXqdot(const double time, Interaction& inter)
{
  DEBUG_PRINT("LagrangianScleronomousR::computeNonLinearH2dot starts");
  // Compute the H Jacobian dot
  LagrangianScleronomousR::computeDotJachq(time, inter);
  _dotjacqhXqdot.reset(new SiconosVector(_dotjachq->size(0)));
  SiconosVector workQdot = *inter.data(q1);
  DEBUG_EXPR(workQdot.display(););
  DEBUG_EXPR(_dotjachq->display(););
  prod(*_dotjachq, workQdot, *_dotjacqhXqdot);
  DEBUG_PRINT("LagrangianScleronomousR::computeNonLinearH2dot ends");
}

void LagrangianScleronomousR::computeOutput(const double time, Interaction& inter, unsigned int derivativeNumber)
{

  DEBUG_PRINTF("LagrangianScleronomousR::computeOutput(const double time, Interaction& inter, unsigned int derivativeNumber) with time = %f and derivativeNumber = %i\n", time, derivativeNumber);

  if (derivativeNumber == 0)
    computeh(time, inter);
  else
  {
    computeJachq(time, inter);

    SiconosVector& y = *inter.y(derivativeNumber);
    if (derivativeNumber == 1)
      prod(*_jachq, *inter.data(q1), y);
    else if (derivativeNumber == 2)
    {
      computeDotJachq(time, inter);
      prod(*_jachq, *inter.data(q2), y);
      prod(*_dotjachq, *inter.data(q1), y, false);
    }
    else
      RuntimeException::selfThrow("LagrangianScleronomousR::computeOutput(t,index), index out of range");
  }
}

void LagrangianScleronomousR::computeInput(const double time, Interaction& inter, unsigned int level)
{
  DEBUG_PRINT("void LagrangianScleronomousR::computeInput(const double time, Interaction& inter, unsigned int level)\n");
  DEBUG_PRINTF("level = %i\n", level);

  computeJachq(time, inter);
  // get lambda of the concerned interaction
  SiconosVector& lambda = *inter.lambda(level);
  // data[name] += trans(G) * lambda
  prod(lambda, *_jachq, *inter.data(p0 + level), false);
  DEBUG_EXPR(inter.data(p0 + level)->display(););
}
const std::string LagrangianScleronomousR::getJachqName() const
{
  if (_pluginJachq->fPtr)
    return _pluginJachq->getPluginName();
  return "unamed";

}
LagrangianScleronomousR* LagrangianScleronomousR::convert(Relation *r)
{
  return dynamic_cast<LagrangianScleronomousR*>(r);
}
