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
  setComputehFunction(SSLH::getPluginName(computeOut), SSLH::getPluginFunctionName(computeOut));
  setComputegFunction(SSLH::getPluginName(computeIn), SSLH::getPluginFunctionName(computeIn));
  // The jacobians are not set, and thus considered as null matrices at this point.
}

FirstOrderType2R::FirstOrderType2R(const string& computeOut, const string& computeIn, const string& computeJX, const string& computeJL):
  FirstOrderR(RELATION::Type2R)
{
  // Size vector of pointers to functions.
  // Connect input and output to plug-in
  setComputehFunction(SSLH::getPluginName(computeOut), SSLH::getPluginFunctionName(computeOut));
  setComputegFunction(SSLH::getPluginName(computeIn), SSLH::getPluginFunctionName(computeIn));

  setComputeJachxFunction(SSLH::getPluginName(computeJX), SSLH::getPluginFunctionName(computeJX));
  setComputeJacglambdaFunction(SSLH::getPluginName(computeJL), SSLH::getPluginFunctionName(computeJL));
}

void FirstOrderType2R::initialize(Interaction& inter)
{
  FirstOrderR::initialize(inter);

  // Check if an Interaction is connected to the Relation.
  unsigned int sizeY = inter.getSizeOfY();
  unsigned int sizeDS = inter.getSizeOfDS();

  // The initialization of each component depends on the way the Relation was built ie if the matrix/vector
  // was read from xml or not
  if (!_jachx)
    _jachx.reset(new SimpleMatrix(sizeY, sizeDS));
  if (!_jachlambda)
    _jachlambda.reset(new SimpleMatrix(sizeY, sizeY));
  if (!_jacgx)
    _jacgx.reset(new SimpleMatrix(sizeDS, sizeDS));
  if (!_jacglambda)
    _jacglambda.reset(new SimpleMatrix(sizeDS, sizeY));


  assert((_jachx->size(1) == sizeDS && _jachx->size(0) == sizeY) &&
         "FirstOrderType2R::initialize inconsistent sizes between _jach[0] matrix and the interaction.");

  assert((_jacglambda->size(0) == sizeDS && _jacglambda->size(1) == sizeY) &&
         "FirstOrderType2R::initialize inconsistent sizes between _jacg[0] matrix and the interaction.");
}

void FirstOrderType2R::computeh(const double time, Interaction& inter)
{
  computeOutput(time, inter, 0);
}

void FirstOrderType2R::computeg(const double time, Interaction& inter)
{
  /* this function calls computeInput, which in return calls this function.
  * Broken logic, put a selfthrow so people may fix it */
  RuntimeException::selfThrow("The logic in this function is broken. Please fix it before using it");
  computeInput(time, inter, 0);
}

void FirstOrderType2R::computeOutput(const double time, Interaction& inter, unsigned int)
{
  computeh(time, inter);
}

void FirstOrderType2R::computeInput(const double time, Interaction& inter, unsigned int level)
{

  /**compute the newr */

  SiconosVector workL = *inter.lambda(level);
  workL -= *(inter.lambdaOld(level));

  //  cout<<"FirstOrderType2R::computeInput : diff lambda"<<endl;
  //  inter.lambdaOld(level)->display();
  //  lambda->display();
  //  _workL->display();
  //  cout<<"FirstOrderType2R::computeInput : g_alpha"<<endl;
  //  _workX->display();
  prod(*B(), workL, *inter.data(g_alpha), false);
  //  cout<<"FirstOrderType2R::computeInput : result g_alpha - B*diffL"<<endl;
  //  _workX->display();

  *inter.data(r) += *inter.data(g_alpha);

  //compute the new g_alpha
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  computeg(time, inter);


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
void FirstOrderType2R::preparNewtonIteration(Interaction& inter)
{

  /* compute the comtribution in xp, for the next iteration */
  SiconosVector& lambda = *inter.lambda(0);
  SiconosVector tmpV = SiconosVector(inter.data(ds_xp)->size());

  prod(*(B()), lambda, tmpV, true);

  *inter.data(ds_xp) -= tmpV;

  *inter.data(ds_xp) += *inter.data(g_alpha);
}

void FirstOrderType2R::computeJachlambda(const double time, Interaction& inter)
{
  RuntimeException::selfThrow("FirstOrderType2R::computeJachlambda must be overload.");
}
void FirstOrderType2R::computeJachx(const double time, Interaction& inter)
{
  RuntimeException::selfThrow("FirstOrderType2R::computeJachx must be overload.");
}

void FirstOrderType2R::computeJacglambda(const double time, Interaction& inter)
{
  RuntimeException::selfThrow("FirstOrderType2R::computeJacglambda must be overload.");
}
void FirstOrderType2R::computeJacgx(const double time, Interaction& inter)
{
  RuntimeException::selfThrow("FirstOrderType2R::computejacgx must be overload.");
}
void FirstOrderType2R::computeJacg(const double time, Interaction& inter)
{
  computeJacglambda(time, inter);
  computeJacgx(time, inter);
}


FirstOrderType2R* FirstOrderType2R::convert(Relation *r)
{
  return dynamic_cast<FirstOrderType2R*>(r);
}

