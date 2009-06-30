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
#include "FirstOrderType2R.h"
#include "RelationXML.h"
#include "Interaction.h"
#include "FirstOrderNonLinearDS.h"
#include "FirstOrderR.cpp"

using namespace std;

FirstOrderType2R::FirstOrderType2R():
  BaseClass(RELATION::Type2R)
{}
// xml constructor
FirstOrderType2R::FirstOrderType2R(SP::RelationXML FORxml):
  BaseClass(FORxml, RELATION::Type2R)
{
  JacH.resize(1);
  JacG.resize(1);
  // input g
  if (FORxml->hasG())
  {
    gName = FORxml->getGPlugin();
    setComputeGFunction(SSL::getPluginName(gName), SSL::getPluginFunctionName(gName));
    // Gradients
    if (!FORxml->hasJacobianG())
      RuntimeException::selfThrow("FirstOrderType2R xml constructor failed. No input for gradient(s) of g function.");

    if (FORxml->isJacobianGPlugin(0))
      JacG[0].reset(new PluggedMatrix(FORxml->getJacobianGPlugin(0)));
    else
      JacG[0].reset(new PluggedMatrix(FORxml->getJacobianGMatrix(0)));
  }

  // output h
  if (FORxml->hasH())
  {
    hName = FORxml->getHPlugin();
    setComputeHFunction(SSL::getPluginName(hName), SSL::getPluginFunctionName(hName));
    // Gradients
    if (!FORxml->hasJacobianH())
      RuntimeException::selfThrow("FirstOrderType2R xml constructor failed. No input for gradients of h function.");
    if (FORxml->isJacobianHPlugin(0))
      JacH[0].reset(new PluggedMatrix(FORxml->getJacobianHPlugin(0)));
    else
      JacH[0].reset(new PluggedMatrix(FORxml->getJacobianHMatrix(0)));
  }
}

FirstOrderType2R::FirstOrderType2R(const string& computeOut, const string& computeIn):
  BaseClass(RELATION::Type2R)
{
  // Size vector of pointers to functions.
  // Connect input and output to plug-in
  setComputeHFunction(SSL::getPluginName(computeOut), SSL::getPluginFunctionName(computeOut));
  setComputeGFunction(SSL::getPluginName(computeIn), SSL::getPluginFunctionName(computeIn));
  // The jacobians are not set, and thus considered as null matrices at this point.
  JacG.resize(1);
  JacH.resize(1);
  JacH[0].reset(new PluggedMatrix());
  JacG[0].reset(new PluggedMatrix());
}

FirstOrderType2R::FirstOrderType2R(const string& computeOut, const string& computeIn, const string& computeJX, const string& computeJL):
  BaseClass(RELATION::Type2R)
{
  // Size vector of pointers to functions.
  // Connect input and output to plug-in
  setComputeHFunction(SSL::getPluginName(computeOut), SSL::getPluginFunctionName(computeOut));
  setComputeGFunction(SSL::getPluginName(computeIn), SSL::getPluginFunctionName(computeIn));
  JacG.resize(1);
  JacH.resize(1);
  JacH[0].reset(new PluggedMatrix(computeJX));
  JacG[0].reset(new PluggedMatrix(computeJL));
}

void FirstOrderType2R::initialize(SP::Interaction inter)
{
  BaseClass::initialize(inter);

  // Check if an Interaction is connected to the Relation.
  unsigned int sizeY = getInteractionPtr()->getSizeOfY();
  unsigned int sizeDS = getInteractionPtr()->getSizeOfDS();
  unsigned int sizeZ = getInteractionPtr()->getSizeZ();
  if (!getInteractionPtr())
    RuntimeException::selfThrow("FirstOrderR::initialize failed. No Interaction linked to the present relation.");

  // Update data member (links to DS variables)
  initDSLinks();
  // Initialize work vectors

  workX.reset(new SimpleVector(sizeDS));
  workZ.reset(new SimpleVector(sizeZ));
  workY.reset(new SimpleVector(sizeY));

  // The initialization of each component depends on the way the Relation was built ie if the matrix/vector
  // was read from xml or not
  if (JacH[0]->size(0) == 0) // if the matrix dim are null
    JacH[0]->resize(sizeY, sizeDS);
  else
    assert((JacH[0]->size(1) == sizeDS && JacH[0]->size(0) == sizeY) &&
           "FirstOrderType2R::initialize inconsistent sizes between JacH[0] matrix and the interaction.");

  // Same work for jacobianLambdaG
  if (JacG[0]->size(0) == 0) // if the matrix dim are null
    JacG[0]->resize(sizeDS, sizeY);
  else
    assert((JacG[0]->size(0) == sizeDS && JacG[0]->size(1) == sizeY) &&
           "FirstOrderType2R::initialize inconsistent sizes between JacG[0] matrix and the interaction.");
}

void FirstOrderType2R::computeH(double t)
{
  computeOutput(t, 0);
}

void FirstOrderType2R::computeG(double t)
{
  computeInput(t, 0);
}

void FirstOrderType2R::computeOutput(double t, unsigned int)
{
  computeH(t);
}

void FirstOrderType2R::computeInput(double t, unsigned int level)
{

  /**compute the newr */
  SP::SiconosVector lambda = getInteractionPtr()->getLambdaPtr(level);
  *workX = *data[g_alpha];
  *workL = *lambda;
  *workL -= *(getInteractionPtr()->getLambdaOldPtr(level));
  //  cout<<"FirstOrderType2R::computeInput : diff lambda"<<endl;
  //  getInteractionPtr()->getLambdaOldPtr(level)->display();
  //  lambda->display();
  //  workL->display();
  //  cout<<"FirstOrderType2R::computeInput : g_alpha"<<endl;
  //  workX->display();
  prod(*(getB()), *workL, *workX, false);
  //  cout<<"FirstOrderType2R::computeInput : result g_alpha - B*diffL"<<endl;
  //  workX->display();
  *data[r] += *workX;
  //  cout<<"FirstOrderType2R::computeInput : r alpha+1"<<endl;
  //  data[r]->display();

  /*compute the new g_alpha*/
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  *workX = *data[x];
  *workZ = *data[z];


  computeG(t);
  //  cout<<"next g_alpha"<<endl;
  //  data[g_alpha]->display();
  /*  unsigned int sizeL = lambda->size();
  unsigned int sizeZ = data[z]->size();
  unsigned int sizeR = workX->size();

  input(sizeL, &(*workL)(0), sizeR, &(*workX)(0),sizeR, &(*workR)(0), sizeZ, &(*workZ)(0));


  *data[g_alpha] = *workR;
  *data[z] = *workZ;
  */

}

void FirstOrderType2R::computeJacH(double, unsigned int index)
{
  //
  /*assert(index==0&&"FirstOrderType2R::computeJacobianH(index): index is out of range");
  assert(JacH[0]->fPtr&&"FirstOrderType2R::computeJacobianH() failed; not linked to a plug-in function.");

  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  *workX = *data[x];
  *workZ = *data[z];

  unsigned int sizeY = getInteractionPtr()->getSizeOfY();
  unsigned int sizeX = data[x]->size();
  unsigned int sizeZ = data[z]->size();

  (JacH[0]->fPtr)(sizeX, &(*workX)(0), sizeY, &(*(JacH[0]))(0,0), sizeZ, &(*workZ)(0));

  // Rebuilt z from Tmp
  *data[z] = *workZ;*/
}

void FirstOrderType2R::computeJacG(double, unsigned int index)
{
  /*assert(index==0&&"FirstOrderType2R::computeJacobianG(index): index is out of range");
  assert(JacG[0]->fPtr&&"FirstOrderType2R::computeJacobianG() failed; not linked to a plug-in function.");

  SP::SiconosVector lambda = getInteractionPtr()->getLambdaPtr(0);
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  *workZ = *data[z];
  *workY = *lambda;

  unsigned int sizeY = lambda->size();
  unsigned int sizeX = data[x]->size();
  unsigned int sizeZ = data[z]->size();

  (JacG[0]->fPtr)(sizeY, &(*workY)(0), sizeX, &(*(JacG[0]))(0,0), sizeZ, &(*workZ)(0));

  // Rebuilt z from Tmp
  *data[z] = *workZ;*/
}

void FirstOrderType2R::preparNewtonIteration()
{

  /** compute the comtribution in xp, for the next iteration*/
  /** */
  SP::SiconosVector lambda = getInteractionPtr()->getLambdaPtr(0);
  *workL = *lambda;

  //  cout<<"FirstOrderType2R::preparNewtonIteration, lambda: \n";
  //  workL->display();

  scal(-1.0, *workL, *workL);
  prod(*(getB()), *workL, *workX, true);

  //  cout<<"FirstOrderType2R::preparNewtonIteration, -B*lambda: \n";
  //  workX->display();

  //  cout<<"FirstOrderType2R::preparNewtonIteration, g_alpha: \n";
  //  data[g_alpha]->display();

  *workX += *data[g_alpha];


  *data[ds_xp] += *workX;
  //  cout<<"FirstOrderType2R::preparNewtonIteration, g_alpha -B*lambda : \n";
  //  workX->display();
}



FirstOrderType2R* FirstOrderType2R::convert(Relation *r)
{
  return dynamic_cast<FirstOrderType2R*>(r);
}

