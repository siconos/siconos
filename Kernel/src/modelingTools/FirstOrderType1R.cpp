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
      //  Jacg[0].reset(new PluggedMatrix(FORxml->getJacobianGPlugin(0)));
      setComputeJacglambdaFunction(SSL::getPluginName(FORxml->getgPlugin()), SSL::getPluginFunctionName(FORxml->getgPlugin()));
    }
    else
    {
      Jacglambda.reset(new SimpleMatrix(FORxml->getJacobianGMatrix(0)));
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
      //  Jach[0].reset(new PluggedMatrix(FORxml->getJacobianHPlugin(0)));
    }
    else
    {
      Jachx.reset(new SimpleMatrix(FORxml->getJacobianHMatrix(0)));
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

void FirstOrderType1R::initialize(SP::Interaction inter)
{
  FirstOrderR::initialize(inter);

  // Check if an Interaction is connected to the Relation.
  unsigned int sizeY = interaction()->getSizeOfY();
  unsigned int sizeDS = interaction()->getSizeOfDS();
  unsigned int sizeZ = interaction()->getSizez();
  if (!interaction())
    RuntimeException::selfThrow("FirstOrderR::initialize failed. No Interaction linked to the present relation.");

  // Update data member (links to DS variables)
  initDSLinks();
  // Initialize work vectors

  _workX.reset(new SiconosVector());
  _workZ.reset(new SiconosVector());
  _workY.reset(new SiconosVector(sizeY));

  // The initialization of each component depends on the way the Relation was built ie if the matrix/vector
  // was read from xml or not
  if (Jachx->size(0) == 0) // if the matrix dim are null
    Jachx->resize(sizeY, sizeDS);
  else
    assert((Jachx->size(1) == sizeDS && Jachx->size(0) == sizeY) &&
           "FirstOrderType1R::initialize inconsistent sizes between Jach[0] matrix and the interaction.");

  // Same work for jacobianLambdaG
  if (Jacglambda->size(0) == 0) // if the matrix dim are null
    Jacglambda->resize(sizeDS, sizeY);
  else
    assert((Jacglambda->size(0) == sizeDS && Jacglambda->size(1) == sizeY) &&
           "FirstOrderType1R::initialize inconsistent sizes between Jacg[0] matrix and the interaction.");
}

void FirstOrderType1R::computeh(double t)
{
  computeOutput(t, 0);
}

void FirstOrderType1R::computeg(double t)
{
  computeInput(t, 0);
}

void FirstOrderType1R::computeOutput(double t, unsigned int)
{
  assert(_pluginh && "FirstOrderType1R::computeOutput() is not linked to a plugin function");

  SP::SiconosVector y = interaction()->y(0);
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.

  *_workX = *data[x];
  *_workZ = *data[z];
  *_workY = *y;

  unsigned int sizeY = y->size();
  unsigned int sizeX = data[x]->size();
  unsigned int sizeZ = data[z]->size(); //

  ((Type1Ptr)(_pluginh->fPtr))(sizeX, &(*_workX)(0), sizeY, &(*_workY)(0), sizeZ, &(*_workZ)(0));

  // Rebuilt y/z from Tmp
  *y = *_workY;
  *data[z] = *_workZ;
}

void FirstOrderType1R::computeInput(double t, unsigned int level)
{
  assert(_pluging && "FirstOrderType1R::computeInput() is not linked to a plugin function");

  SP::SiconosVector lambda = interaction()->lambda(level);
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.

  *_workX = *data[r];
  *_workZ = *data[z];
  *_workY = *lambda;

  unsigned int sizeY = lambda->size();
  unsigned int sizeZ = data[z]->size();
  unsigned int sizeR = _workX->size();


  ((Type1Ptr)(_pluging->fPtr))(sizeY, &(*_workY)(0), sizeR, &(*_workX)(0), sizeZ, &(*_workZ)(0));

  *data[r] = *_workX;
  *data[z] = *_workZ;
}

void FirstOrderType1R::computeJachx(double)
{
  //
  assert(_pluginJachx && "FirstOrderType1R::computeJacobianH() failed; not linked to a plug-in function.");

  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  *_workX = *data[x];
  *_workZ = *data[z];

  unsigned int sizeY = interaction()->getSizeOfY();
  unsigned int sizeX = data[x]->size();
  unsigned int sizeZ = data[z]->size();

  ((Type1Ptr)(_pluginJachx->fPtr))(sizeX, &(*_workX)(0), sizeY, &(*(Jachx))(0, 0), sizeZ, &(*_workZ)(0));

  // Rebuilt z from Tmp
  *data[z] = *_workZ;
}

void FirstOrderType1R::computeJacglambda(double)
{
  assert(_pluginJacLg && "FirstOrderType1R::computeJacobiang() failed; not linked to a plug-in function.");

  SP::SiconosVector lambda = interaction()->lambda(0);
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  *_workZ = *data[z];
  *_workY = *lambda;

  unsigned int sizeY = lambda->size();
  unsigned int sizeX = data[x]->size();
  unsigned int sizeZ = data[z]->size();

  ((Type1Ptr)(_pluginJacLg->fPtr))(sizeY, &(*_workY)(0), sizeX, &(*(Jacglambda))(0, 0), sizeZ, &(*_workZ)(0));

  // Rebuilt z from Tmp
  *data[z] = *_workZ;
}

FirstOrderType1R* FirstOrderType1R::convert(Relation *r)
{
  return dynamic_cast<FirstOrderType1R*>(r);
}

