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
#include "FirstOrderType1R.hpp"
#include "Interaction.hpp"
#include "FirstOrderNonLinearDS.hpp"

#include "BlockVector.hpp"


FirstOrderType1R::FirstOrderType1R(const std::string& computeOut, const std::string& computeIn):
  FirstOrderR(RELATION::Type1R)
{
  // Size vector of pointers to functions.
  // Connect input and output to plug-in
  setComputehFunction(SSLH::getPluginName(computeOut), SSLH::getPluginFunctionName(computeOut));
  setComputegFunction(SSLH::getPluginName(computeIn), SSLH::getPluginFunctionName(computeIn));
  // The jacobians are not set, and thus considered as null matrices at this point.
}

FirstOrderType1R::FirstOrderType1R(const std::string& computeOut, const std::string& computeIn, const std::string& computeJX, const std::string& computeJL):
  FirstOrderR(RELATION::Type1R)
{
  // Size vector of pointers to functions.
  // Connect input and output to plug-in
  setComputehFunction(SSLH::getPluginName(computeOut), SSLH::getPluginFunctionName(computeOut));
  setComputegFunction(SSLH::getPluginName(computeIn), SSLH::getPluginFunctionName(computeIn));
  setComputeJachxFunction(SSLH::getPluginName(computeJX), SSLH::getPluginFunctionName(computeJX));
  setComputeJacglambdaFunction(SSLH::getPluginName(computeJL), SSLH::getPluginFunctionName(computeJL));
}

void FirstOrderType1R::initialize(Interaction& inter)
{
  FirstOrderR::initialize(inter);

  // Check if an Interaction is connected to the Relation.
  unsigned int sizeY = inter.getSizeOfY();
  unsigned int sizeDS = inter.getSizeOfDS();
  unsigned int sizeZ = inter.data(z)->size();

  if (!_jachx)
    _jachx.reset(new SimpleMatrix(sizeY, sizeDS));
  if (!_jachz)
    _jachz.reset(new SimpleMatrix(sizeY, sizeZ));
  if (!_jacglambda)
    _jacglambda.reset(new SimpleMatrix(sizeDS, sizeY));

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

void FirstOrderType1R::computeh(double time, Interaction& inter)
{
  computeOutput(time, inter, 0);
}

void FirstOrderType1R::computeg(double time, Interaction& inter)
{
  computeInput(time, inter, 0);
}

void FirstOrderType1R::computeOutput(double time, Interaction& inter, unsigned int)
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

void FirstOrderType1R::computeInput(double time, Interaction& inter, unsigned int level)
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

void FirstOrderType1R::computeJachx(double time, Interaction& inter)
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

void FirstOrderType1R::computeJachz(double time, Interaction& inter)
{
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  SiconosVector workX = *inter.data(x);
  SiconosVector workZ = *inter.data(z);

  unsigned int sizeZ = inter.data(z)->size();

  if (_pluginJachz && _pluginJachz->fPtr)
    ((Type1Ptr)(_pluginJachz->fPtr))(workX.size(), &(workX)(0), sizeZ, &(*(_jachz))(0, 0), workZ.size(), &(workZ)(0));

  // Rebuilt z from Tmp
  *inter.data(z) = workZ;
}

void FirstOrderType1R::computeJacglambda(double time, Interaction& inter)
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

