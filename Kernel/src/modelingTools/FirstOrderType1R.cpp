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
#include "FirstOrderType1R.h"
#include "RelationXML.h"
#include "Interaction.h"
#include "FirstOrderNonLinearDS.h"

using namespace std;

// xml constructor
FirstOrderType1R::FirstOrderType1R(SP::RelationXML FORxml):
  FirstOrderR(FORxml, RELATION::Type1R)
{
  // input g
  if (FORxml->hasG())
  {
    gName = FORxml->getGPlugin();
    setComputeGFunction(SSL::getPluginName(gName), SSL::getPluginFunctionName(gName));
    // Gradients
    if (!FORxml->hasJacobianG())
      RuntimeException::selfThrow("FirstOrderType1R xml constructor failed. No input for gradient(s) of g function.");

    if (FORxml->isJacobianGPlugin(0))
    {
      //  JacG[0].reset(new PluggedMatrix(FORxml->getJacobianGPlugin(0)));
      setComputeJacLGFunction(SSL::getPluginName(gName), SSL::getPluginFunctionName(gName));
    }
    else
    {
      JacLG.reset(new SimpleMatrix(FORxml->getJacobianGMatrix(0)));
    }
  }

  // output h
  if (FORxml->hasH())
  {
    hName = FORxml->getHPlugin();
    setComputeHFunction(SSL::getPluginName(hName), SSL::getPluginFunctionName(hName));
    // Gradients
    if (!FORxml->hasJacobianH())
      RuntimeException::selfThrow("FirstOrderType1R xml constructor failed. No input for gradients of h function.");
    if (FORxml->isJacobianHPlugin(0))
    {
      setComputeJacXHFunction(SSL::getPluginName(gName), SSL::getPluginFunctionName(gName));
      //  JacH[0].reset(new PluggedMatrix(FORxml->getJacobianHPlugin(0)));
    }
    else
    {
      JacXH.reset(new SimpleMatrix(FORxml->getJacobianHMatrix(0)));
    }
  }
}

FirstOrderType1R::FirstOrderType1R(const string& computeOut, const string& computeIn):
  FirstOrderR(RELATION::Type1R)
{
  // Size vector of pointers to functions.
  // Connect input and output to plug-in
  setComputeHFunction(SSL::getPluginName(computeOut), SSL::getPluginFunctionName(computeOut));
  setComputeGFunction(SSL::getPluginName(computeIn), SSL::getPluginFunctionName(computeIn));
  // The jacobians are not set, and thus considered as null matrices at this point.
}

FirstOrderType1R::FirstOrderType1R(const string& computeOut, const string& computeIn, const string& computeJX, const string& computeJL):
  FirstOrderR(RELATION::Type1R)
{
  // Size vector of pointers to functions.
  // Connect input and output to plug-in
  setComputeHFunction(SSL::getPluginName(computeOut), SSL::getPluginFunctionName(computeOut));
  setComputeGFunction(SSL::getPluginName(computeIn), SSL::getPluginFunctionName(computeIn));
  setComputeJacXHFunction(SSL::getPluginName(computeJX), SSL::getPluginFunctionName(computeJX));
  setComputeJacLGFunction(SSL::getPluginName(computeJL), SSL::getPluginFunctionName(computeJL));
}

void FirstOrderType1R::initialize(SP::Interaction inter)
{
  FirstOrderR::initialize(inter);

  // Check if an Interaction is connected to the Relation.
  unsigned int sizeY = interaction()->getSizeOfY();
  unsigned int sizeDS = interaction()->getSizeOfDS();
  unsigned int sizeZ = interaction()->getSizeZ();
  if (!interaction())
    RuntimeException::selfThrow("FirstOrderR::initialize failed. No Interaction linked to the present relation.");

  // Update data member (links to DS variables)
  initDSLinks();
  // Initialize work vectors

  _workX.reset(new SimpleVector(sizeDS));
  _workZ.reset(new SimpleVector(sizeZ));
  _workY.reset(new SimpleVector(sizeY));

  // The initialization of each component depends on the way the Relation was built ie if the matrix/vector
  // was read from xml or not
  if (JacXH->size(0) == 0) // if the matrix dim are null
    JacXH->resize(sizeY, sizeDS);
  else
    assert((JacXH->size(1) == sizeDS && JacXH->size(0) == sizeY) &&
           "FirstOrderType1R::initialize inconsistent sizes between JacH[0] matrix and the interaction.");

  // Same work for jacobianLambdaG
  if (JacLG->size(0) == 0) // if the matrix dim are null
    JacLG->resize(sizeDS, sizeY);
  else
    assert((JacLG->size(0) == sizeDS && JacLG->size(1) == sizeY) &&
           "FirstOrderType1R::initialize inconsistent sizes between JacG[0] matrix and the interaction.");
}

void FirstOrderType1R::computeH(double t)
{
  computeOutput(t, 0);
}

void FirstOrderType1R::computeG(double t)
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

void FirstOrderType1R::computeJacXH(double)
{
  //
  assert(index == 0 && "FirstOrderType1R::computeJacobianH(index): index is out of range");
  assert(_plunginJacxh && "FirstOrderType1R::computeJacobianH() failed; not linked to a plug-in function.");

  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  *_workX = *data[x];
  *_workZ = *data[z];

  unsigned int sizeY = interaction()->getSizeOfY();
  unsigned int sizeX = data[x]->size();
  unsigned int sizeZ = data[z]->size();

  ((Type1Ptr)(_plunginJacxh->fPtr))(sizeX, &(*_workX)(0), sizeY, &(*(JacXH))(0, 0), sizeZ, &(*_workZ)(0));

  // Rebuilt z from Tmp
  *data[z] = *_workZ;
}

void FirstOrderType1R::computeJacLG(double)
{
  assert(index == 0 && "FirstOrderType1R::computeJacobianG(index): index is out of range");
  assert(_pluginJacLg && "FirstOrderType1R::computeJacobianG() failed; not linked to a plug-in function.");

  SP::SiconosVector lambda = interaction()->lambda(0);
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  *_workZ = *data[z];
  *_workY = *lambda;

  unsigned int sizeY = lambda->size();
  unsigned int sizeX = data[x]->size();
  unsigned int sizeZ = data[z]->size();

  ((Type1Ptr)(_pluginJacLg->fPtr))(sizeY, &(*_workY)(0), sizeX, &(*(JacLG))(0, 0), sizeZ, &(*_workZ)(0));

  // Rebuilt z from Tmp
  *data[z] = *_workZ;
}

FirstOrderType1R* FirstOrderType1R::convert(Relation *r)
{
  return dynamic_cast<FirstOrderType1R*>(r);
}

