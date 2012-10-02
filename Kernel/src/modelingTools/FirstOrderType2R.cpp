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
#include "FirstOrderType2R.hpp"
#include "RelationXML.hpp"
#include "Interaction.hpp"
#include "FirstOrderNonLinearDS.hpp"

using namespace std;
#define DEBUG_STDOUT
//#define DEBUG_MESSAGES 1

#include <debug.h>


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
  /* this function calls computeOutput, which in return calls this function.
   * Broken logic, put a selfthrow so people may fix it */
  RuntimeException::selfThrow("FirstOrderType2R::computeh. The logic in this function is broken. Please fix it before using it");
  computeOutput(time, inter, 0);
}

void FirstOrderType2R::computeg(const double time, Interaction& inter)
{
  /* this function calls computeInput, which in return calls this function.
  * Broken logic, put a selfthrow so people may fix it */
  RuntimeException::selfThrow("FirstOrderType2R::computeg . The logic in this function is broken. Please fix it before using it");
  computeInput(time, inter, 0);
}

void FirstOrderType2R::computeOutput(const double time, Interaction& inter, unsigned int level)
{
  DEBUG_PRINT("FirstOrderType2R::computeOutput \n");
  // compute the new y  obtained by linearisation (see DevNotes)
  // y_{alpha+1}_{k+1} = h(x_{k+1}^{alpha},lambda_{k+1}^{alpha},t_k+1)
  //                     + C_{k+1}^alpha ( x_{k+1}^{alpha+1}- x_{k+1}^{alpha} )
  //                     + D_{k+1}^alpha ( lambda_{k+1}^{alpha+1} - lambda_{k+1}^{alpha} )
  // or equivalently
  // y_{alpha+1}_{k+1} = y_{alpha}_{k+1} - ResiduY_{k+1}^{alpha}
  //                     + C_{k+1}^alpha ( x_{k+1}^{alpha+1}- x_{k+1}^{alpha} )
  //                     + D_{k+1}^alpha ( lambda_{k+1}^{alpha+1} - lambda_{k+1}^{alpha} )


  DEBUG_PRINT("FirstOrderType2R::computeOutput : y(level) before \n");
#ifdef DEBUG_MESSAGES
  inter.y(level)->display();
#endif
  //SiconosVector yOld = *inter.yOld(0); // Retrieve  y_{alpha}_{k+1}
  DEBUG_PRINT("FirstOrderType2R::computeOutput : yOld(level) \n");
#ifdef DEBUG_MESSAGES
  inter.yOld(level)->display();
#endif
  *inter.y(level) = *inter.yOld(level);

  DEBUG_PRINT("FirstOrderType2R::computeOutput : ResiduY() \n");
#ifdef DEBUG_MESSAGES
  inter.residuY()->display();
#endif
  *inter.y(level) -= *inter.residuY();

  DEBUG_PRINT("FirstOrderType2R::computeOutput : y(level) \n");
#ifdef DEBUG_MESSAGES
  inter.y(level)->display();
#endif
  SiconosVector Deltax = *inter.data(deltax);
  //  deltax -= *(inter.data(xold));
  DEBUG_PRINT("FirstOrderType2R::computeOutput : deltax \n");
#ifdef DEBUG_MESSAGES
  Deltax.display();
#endif

  prod(*C(), Deltax, *inter.y(level), false);


  SiconosVector deltalambda = *inter.lambda(level);
  deltalambda -= *(inter.lambdaOld(level));

  DEBUG_PRINT("FirstOrderType2R::computeOutput : deltalambda \n");
#ifdef DEBUG_MESSAGES
  deltalambda.display();
#endif
  prod(*D(), deltalambda, *inter.y(level), false);

  DEBUG_PRINT("FirstOrderType2R::computeOutput : new linearized y \n");
#ifdef DEBUG_MESSAGES
  inter.y(level)->display();
#endif

  computeh(time, inter);





}

void FirstOrderType2R::computeInput(const double time, Interaction& inter, unsigned int level)
{
  DEBUG_PRINT("FirstOrderType2R::computeInput \n");
  // compute the new r  obtained by linearisation
  // r_{alpha+1}_{k+1} = g(x_{k+1}^{alpha},lambda_{k+1}^{alpha},t_k+1)
  //                     + B_{k+1}^alpha ( lambda_{k+1}^{alpha+1}- lambda_{k+1}^{alpha} )


  SiconosVector lambda = *inter.lambda(level);
  lambda -= *(inter.lambdaOld(level));

  //  cout<<"FirstOrderType2R::computeInput : diff lambda"<<endl;
  //  inter.lambdaOld(level)->display();
  //  lambda->display();
  //  _lambda->display();
  //  cout<<"FirstOrderType2R::computeInput : g_alpha"<<endl;
  //  _workX->display();
  prod(*B(), lambda, *inter.data(g_alpha), false);
  //  cout<<"FirstOrderType2R::computeInput : result g_alpha - B*diffL"<<endl;
  //  _workX->display();

  *inter.data(r) += *inter.data(g_alpha);

  //compute the new g_alpha
  // Warning: temporary method to have contiguous values in memory, copy of block to simple.
  computeg(time, inter);


  //  cout<<"next g_alpha"<<endl;
  //  data[g_alpha]->display();


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

