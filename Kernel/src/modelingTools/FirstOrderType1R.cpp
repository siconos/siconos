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
#include "FirstOrderType1R.hpp"
#include "RelationXML.hpp"
#include "Interaction.hpp"
#include "FirstOrderNonLinearDS.hpp"

using namespace std;

// xml constructor
FirstOrderType1R::FirstOrderType1R(SP::RelationXML FORxml):
  FirstOrderR(FORxml, RELATION::Type1R)
{
  // input g
  if (FORxml->hasG())
  {
    setComputegFunction(SSL::getPluginName(FORxml->getgPlugin()), SSL::getPluginFunctionName(FORxml->getgPlugin()));
    // Gradients
    if (!FORxml->hasJacobianG())
      RuntimeException::selfThrow("FirstOrderType1R xml constructor failed. No input for gradient(s) of g function.");

    if (FORxml->isJacobianGPlugin(0))
    {
      //Jacg[0].reset(new PluggedMatrix(FORxml->getJacobianGPlugin(0)));
      setComputeJacglambdaFunction(SSL::getPluginName(FORxml->getgPlugin()), SSL::getPluginFunctionName(FORxml->getgPlugin()));
    }
    else
    {
      _jacglambda.reset(new SimpleMatrix(FORxml->getJacobianGMatrix(0)));
    }
  }

  // output h
  if (FORxml->hasH())
  {
    setComputehFunction(SSL::getPluginName(FORxml->gethPlugin()), SSL::getPluginFunctionName(FORxml->gethPlugin()));
    // Gradients
    if (!FORxml->hasJacobianH())
      RuntimeException::selfThrow("FirstOrderType1R xml constructor failed. No input for gradients of h function.");
    if (FORxml->isJacobianHPlugin(0))
    {
      setComputeJachxFunction(SSL::getPluginName(FORxml->getgPlugin()), SSL::getPluginFunctionName(FORxml->getgPlugin()));
      //Jach[0].reset(new PluggedMatrix(FORxml->getJacobianHPlugin(0)));
    }
    else
    {
      _jachx.reset(new SimpleMatrix(FORxml->getJacobianHMatrix(0)));
    }
  }
}

FirstOrderType1R::FirstOrderType1R(const string& computeOut, const string& computeIn):
  FirstOrderR(RELATION::Type1R)
{
  // Size vector of pointers to functions.
  // Connect input and output to plug-in
  setComputehFunction(SSL::getPluginName(computeOut), SSL::getPluginFunctionName(computeOut));
  setComputegFunction(SSL::getPluginName(computeIn), SSL::getPluginFunctionName(computeIn));
  // The jacobians are not set, and thus considered as null matrices at this point.
}

FirstOrderType1R::FirstOrderType1R(const string& computeOut, const string& computeIn, const string& computeJX, const string& computeJL):
  FirstOrderR(RELATION::Type1R)
{
  // Size vector of pointers to functions.
  // Connect input and output to plug-in
  setComputehFunction(SSL::getPluginName(computeOut), SSL::getPluginFunctionName(computeOut));
  setComputegFunction(SSL::getPluginName(computeIn), SSL::getPluginFunctionName(computeIn));
  setComputeJachxFunction(SSL::getPluginName(computeJX), SSL::getPluginFunctionName(computeJX));
  setComputeJacglambdaFunction(SSL::getPluginName(computeJL), SSL::getPluginFunctionName(computeJL));
}

void FirstOrderType1R::initialize(Interaction& inter)
{
  FirstOrderR::initialize(inter);

  // Check if an Interaction is connected to the Relation.
  unsigned int sizeY = inter.getSizeOfY();
  unsigned int sizeDS = inter.getSizeOfDS();

  // The initialization of each component depends on the way the Relation was built ie if the matrix/vector
  // was read from xml or not
  if (_jachx->size(0) == 0) // if the matrix dim are null
    _jachx->resize(sizeY, sizeDS);
  else
    assert((_jachx->size(1) == sizeDS && _jachx->size(0) == sizeY) &&
           "FirstOrderType1R::initialize inconsistent sizes between _jach[0] matrix and the interaction.");

  // Same work for jacobianLambdaG
  if (_jacglambda->size(0) == 0) // if the matrix dim are null
    _jacglambda->resize(sizeDS, sizeY);
  else
    assert((_jacglambda->size(0) == sizeDS && _jacglambda->size(1) == sizeY) &&
           "FirstOrderType1R::initialize inconsistent sizes between _jacg[0] matrix and the interaction.");
}

void FirstOrderType1R::computeh(const double time, Interaction& inter)
{
  computeOutput(time, inter, 0);
}

void FirstOrderType1R::computeg(const double time, Interaction& inter)
{
  computeInput(time, inter, 0);
}

void FirstOrderType1R::computeOutput(const double time, Interaction& inter, unsigned int)
{
  assert(_pluginh && "FirstOrderType1R::computeOutput() is not linked to a plugin function");

  SiconosVector& y = *inter.y(0);
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.

  SiconosVector workX = *inter.data(x);
  SiconosVector workZ = *inter.data(z);

  ((Type1Ptr)(_pluginh->fPtr))(workX.size(), &(workX)(0), y.size(), &(y)(0), workZ.size(), &(workZ)(0));

  // Rebuilt z from Tmp
  *inter.data(z) = workZ;
}

void FirstOrderType1R::computeInput(const double time, Interaction& inter, unsigned int level)
{
  assert(_pluging && "FirstOrderType1R::computeInput() is not linked to a plugin function");

  SiconosVector& lambda = *inter.lambda(level);
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.

  SiconosVector workR = *inter.data(r);
  SiconosVector workZ = *inter.data(z);

  ((Type1Ptr)(_pluging->fPtr))(lambda.size(), &(lambda)(0), workR.size(), &(workR)(0), workZ.size(), &(workZ)(0));

  *inter.data(r) = workR;
  *inter.data(z) = workZ;
}

void FirstOrderType1R::computeJachx(const double time, Interaction& inter)
{
  //
  assert(_pluginJachx && "FirstOrderType1R::computeJacobianH() failed; not linked to a plug-in function.");

  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  SiconosVector workX = *inter.data(x);
  SiconosVector workZ = *inter.data(z);

  unsigned int sizeY = inter.getSizeOfY();

  ((Type1Ptr)(_pluginJachx->fPtr))(workX.size(), &(workX)(0), sizeY, &(*(_jachx))(0, 0), workZ.size(), &(workZ)(0));

  // Rebuilt z from Tmp
  *inter.data(z) = workZ;
}

void FirstOrderType1R::computeJacglambda(const double time, Interaction& inter)
{
  assert(_pluginJacLg && "FirstOrderType1R::computeJacobiang() failed; not linked to a plug-in function.");

  SiconosVector& lambda = *inter.lambda(0);
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  SiconosVector workZ = *inter.data(z);

  unsigned int sizeX = inter.data(x)->size();

  ((Type1Ptr)(_pluginJacLg->fPtr))(lambda.size(), &(lambda)(0), sizeX, &(*(_jacglambda))(0, 0), workZ.size(), &(workZ)(0));

  // Rebuilt z from Tmp
  *inter.data(z) = workZ;
}

FirstOrderType1R* FirstOrderType1R::convert(Relation *r)
{
  return dynamic_cast<FirstOrderType1R*>(r);
}

