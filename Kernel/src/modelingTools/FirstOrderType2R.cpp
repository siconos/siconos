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

using namespace std;

FirstOrderType2R::FirstOrderType2R():
  FirstOrderR(RELATION::Type2R)
{}
// xml constructor
FirstOrderType2R::FirstOrderType2R(SP::RelationXML FORxml):
  FirstOrderR(FORxml, RELATION::Type2R)
{
  RuntimeException::selfThrow("FirstOrderR::FirstOrderType2R xml constructor not implemented.");
}

FirstOrderType2R::FirstOrderType2R(const string& computeOut, const string& computeIn):
  FirstOrderR(RELATION::Type2R)
{
  // Size vector of pointers to functions.
  // Connect input and output to plug-in
  setComputeHFunction(SSL::getPluginName(computeOut), SSL::getPluginFunctionName(computeOut));
  setComputeGFunction(SSL::getPluginName(computeIn), SSL::getPluginFunctionName(computeIn));
  // The jacobians are not set, and thus considered as null matrices at this point.
}

FirstOrderType2R::FirstOrderType2R(const string& computeOut, const string& computeIn, const string& computeJX, const string& computeJL):
  FirstOrderR(RELATION::Type2R)
{
  // Size vector of pointers to functions.
  // Connect input and output to plug-in
  setComputeHFunction(SSL::getPluginName(computeOut), SSL::getPluginFunctionName(computeOut));
  setComputeGFunction(SSL::getPluginName(computeIn), SSL::getPluginFunctionName(computeIn));

  setComputeJacXHFunction(SSL::getPluginName(computeJX), SSL::getPluginFunctionName(computeJX));
  setComputeJacLGFunction(SSL::getPluginName(computeJL), SSL::getPluginFunctionName(computeJL));
}

void FirstOrderType2R::initialize(SP::Interaction inter)
{
  FirstOrderR::initialize(inter);

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
  if (!JacXH)
    JacXH.reset(new SimpleMatrix(sizeY, sizeDS));
  if (!JacLH)
    JacLH.reset(new SimpleMatrix(sizeY, sizeY));
  if (!JacXG)
    JacXG.reset(new SimpleMatrix(sizeDS, sizeDS));
  if (!JacLG)
    JacLG.reset(new SimpleMatrix(sizeDS, sizeY));


  assert((JacXH->size(1) == sizeDS && JacXH->size(0) == sizeY) &&
         "FirstOrderType2R::initialize inconsistent sizes between JacH[0] matrix and the interaction.");

  assert((JacLG->size(0) == sizeDS && JacLG->size(1) == sizeY) &&
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
  prod(*getBPtr(), *workL, *workX, false);
  //  cout<<"FirstOrderType2R::computeInput : result g_alpha - B*diffL"<<endl;
  //  workX->display();
  *data[r] += *workX;

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

void FirstOrderType2R::computeJacLH(double)
{
  RuntimeException::selfThrow("FirstOrderType2R::computeJacLH must be overload.");
}
void FirstOrderType2R::computeJacXH(double)
{
  RuntimeException::selfThrow("FirstOrderType2R::computeJacXH must be overload.");
}

void FirstOrderType2R::computeJacLG(double)
{
  RuntimeException::selfThrow("FirstOrderType2R::computeJacLG must be overload.");
}
void FirstOrderType2R::computeJacXG(double)
{
  RuntimeException::selfThrow("FirstOrderType2R::computeJacXG must be overload.");
}
void FirstOrderType2R::computeJacG(double t)
{
  computeJacLG(t);
  computeJacXG(t);
}

void FirstOrderType2R::preparNewtonIteration()
{

  /** compute the comtribution in xp, for the next iteration*/
  /** */
  SP::SiconosVector lambda = getInteractionPtr()->getLambdaPtr(0);
  *workL = *lambda;

  //     cout<<"FirstOrderType2R::preparNewtonIteration, lambda: \n";
  //     workL->display();

  scal(-1.0, *workL, *workL);
  prod(*(getBPtr()), *workL, *workX, true);

  //      cout<<"FirstOrderType2R::preparNewtonIteration, -B*lambda: \n";
  //      workX->display();

  //      cout<<"FirstOrderType2R::preparNewtonIteration, g_alpha: \n";
  //      data[g_alpha]->display();

  *workX += *data[g_alpha];


  *data[ds_xp] += *workX;
  //     cout<<"FirstOrderType2R::preparNewtonIteration,xp= g_alpha -B*lambda : \n";
  //     workX->display();
}



FirstOrderType2R* FirstOrderType2R::convert(Relation *r)
{
  return dynamic_cast<FirstOrderType2R*>(r);
}

