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
#include "FirstOrderType2R.hpp"
#include "RelationXML.hpp"
#include "Interaction.hpp"
#include "FirstOrderNonLinearDS.hpp"

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
  setComputehFunction(SSL::getPluginName(computeOut), SSL::getPluginFunctionName(computeOut));
  setComputegFunction(SSL::getPluginName(computeIn), SSL::getPluginFunctionName(computeIn));
  // The jacobians are not set, and thus considered as null matrices at this point.
}

FirstOrderType2R::FirstOrderType2R(const string& computeOut, const string& computeIn, const string& computeJX, const string& computeJL):
  FirstOrderR(RELATION::Type2R)
{
  // Size vector of pointers to functions.
  // Connect input and output to plug-in
  setComputehFunction(SSL::getPluginName(computeOut), SSL::getPluginFunctionName(computeOut));
  setComputegFunction(SSL::getPluginName(computeIn), SSL::getPluginFunctionName(computeIn));

  setComputeJachxFunction(SSL::getPluginName(computeJX), SSL::getPluginFunctionName(computeJX));
  setComputeJacglambdaFunction(SSL::getPluginName(computeJL), SSL::getPluginFunctionName(computeJL));
}

void FirstOrderType2R::initialize(SP::Interaction inter)
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

  _workX.reset(new SiconosVector(sizeDS));
  _workZ.reset(new SiconosVector(sizeZ));
  _workY.reset(new SiconosVector(sizeY));

  // The initialization of each component depends on the way the Relation was built ie if the matrix/vector
  // was read from xml or not
  if (!Jachx)
    Jachx.reset(new SimpleMatrix(sizeY, sizeDS));
  if (!_jachlambda)
    _jachlambda.reset(new SimpleMatrix(sizeY, sizeY));
  if (!jacgx)
    jacgx.reset(new SimpleMatrix(sizeDS, sizeDS));
  if (!Jacglambda)
    Jacglambda.reset(new SimpleMatrix(sizeDS, sizeY));


  assert((Jachx->size(1) == sizeDS && Jachx->size(0) == sizeY) &&
         "FirstOrderType2R::initialize inconsistent sizes between Jach[0] matrix and the interaction.");

  assert((Jacglambda->size(0) == sizeDS && Jacglambda->size(1) == sizeY) &&
         "FirstOrderType2R::initialize inconsistent sizes between Jacg[0] matrix and the interaction.");
}

void FirstOrderType2R::computeh(double t)
{
  computeOutput(t, 0);
}

void FirstOrderType2R::computeg(double t)
{
  computeInput(t, 0);
}

void FirstOrderType2R::computeOutput(double t, unsigned int)
{
  computeh(t);
}

void FirstOrderType2R::computeInput(double t, unsigned int level)
{

  /**compute the newr */
  SP::SiconosVector lambda = interaction()->lambda(level);
  *_workX = *data[g_alpha];
  *_workL = *lambda;
  *_workL -= *(interaction()->lambdaOld(level));
  //  cout<<"FirstOrderType2R::computeInput : diff lambda"<<endl;
  //  interaction()->lambdaOld(level)->display();
  //  lambda->display();
  //  _workL->display();
  //  cout<<"FirstOrderType2R::computeInput : g_alpha"<<endl;
  //  _workX->display();
  prod(*B(), *_workL, *_workX, false);
  //  cout<<"FirstOrderType2R::computeInput : result g_alpha - B*diffL"<<endl;
  //  _workX->display();

  SP::BlockVector tmp(new BlockVector(*data[r]));
  *tmp = *_workX;
  *data[r] += *tmp;

  /*compute the new g_alpha*/
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  *_workX = *data[x];
  *_workZ = *data[z];


  computeg(t);
  //  cout<<"next g_alpha"<<endl;
  //  data[g_alpha]->display();
  /*  unsigned int sizeL = lambda->size();
  unsigned int sizeZ = data[z]->size();
  unsigned int sizeR = _workX->size();

  input(sizeL, &(*_workL)(0), sizeR, &(*_workX)(0),sizeR, &(*__workR)(0), sizeZ, &(*_workZ)(0));


  *data[g_alpha] = *_workR;
  *data[z] = *_workZ;
  */

}

void FirstOrderType2R::computeJachlambda(double)
{
  RuntimeException::selfThrow("FirstOrderType2R::computeJachlambda must be overload.");
}
void FirstOrderType2R::computeJachx(double)
{
  RuntimeException::selfThrow("FirstOrderType2R::computeJachx must be overload.");
}

void FirstOrderType2R::computeJacglambda(double)
{
  RuntimeException::selfThrow("FirstOrderType2R::computeJacglambda must be overload.");
}
void FirstOrderType2R::computeJacgx(double)
{
  RuntimeException::selfThrow("FirstOrderType2R::computejacgx must be overload.");
}
void FirstOrderType2R::computeJacg(double t)
{
  computeJacglambda(t);
  computeJacgx(t);
}

void FirstOrderType2R::preparNewtonIteration()
{

  /** compute the comtribution in xp, for the next iteration*/
  /** */
  SP::SiconosVector lambda = interaction()->lambda(0);
  *_workL = *lambda;

  //     cout<<"FirstOrderType2R::preparNewtonIteration, lambda: \n";
  //     _workL->display();

  scal(-1.0, *_workL, *_workL);
  prod(*(B()), *_workL, *_workX, true);

  //      cout<<"FirstOrderType2R::preparNewtonIteration, -B*lambda: \n";
  //      _workX->display();

  //      cout<<"FirstOrderType2R::preparNewtonIteration, g_alpha: \n";
  //      data[g_alpha]->display();

  SP::SiconosVector tmp(new SiconosVector(*_workX));
  *tmp = *data[g_alpha];
  *_workX += *tmp;

  SP::BlockVector ntmp(new BlockVector(*data[ds_xp]));
  *ntmp = *_workX;
  *data[ds_xp] += *ntmp;
  //     cout<<"FirstOrderType2R::preparNewtonIteration,xp= g_alpha -B*lambda : \n";
  //     _workX->display();
}



FirstOrderType2R* FirstOrderType2R::convert(Relation *r)
{
  return dynamic_cast<FirstOrderType2R*>(r);
}

