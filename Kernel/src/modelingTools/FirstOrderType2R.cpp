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
#include "Interaction.hpp"
#include "FirstOrderNonLinearDS.hpp"

#include "BlockVector.hpp"

#define DEBUG_STDOUT
#define DEBUG_MESSAGES 1

#include <debug.h>


FirstOrderType2R::FirstOrderType2R():
  FirstOrderR(RELATION::Type2R)
{}

FirstOrderType2R::FirstOrderType2R(const std::string& pluginh, const std::string& pluging):
  FirstOrderR(RELATION::Type2R)
{
  // Size vector of pointers to functions.
  // Connect input and output to plug-in
  setComputehFunction(SSLH::getPluginName(pluginh), SSLH::getPluginFunctionName(pluginh));
  setComputegFunction(SSLH::getPluginName(pluging), SSLH::getPluginFunctionName(pluging));
  // The jacobians are not set, and thus considered as null matrices at this point.
}

FirstOrderType2R::FirstOrderType2R(const std::string& pluginh, const std::string& pluging, const std::string& pluginJacobianhx, const std::string& pluginJacobianglambda):
  FirstOrderR(RELATION::Type2R)
{
  // Size vector of pointers to functions.
  // Connect input and output to plug-in
  setComputehFunction(SSLH::getPluginName(pluginh), SSLH::getPluginFunctionName(pluginh));
  setComputegFunction(SSLH::getPluginName(pluging), SSLH::getPluginFunctionName(pluging));

  setComputeJachxFunction(SSLH::getPluginName(pluginJacobianhx), SSLH::getPluginFunctionName(pluginJacobianhx));
  setComputeJacglambdaFunction(SSLH::getPluginName(pluginJacobianglambda), SSLH::getPluginFunctionName(pluginJacobianglambda));
}

void FirstOrderType2R::initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM)
{

  // Check if an Interaction is connected to the Relation.
  unsigned int sizeY = inter.getSizeOfY();
  unsigned int sizeDS = inter.getSizeOfDS();

  workV.resize(FirstOrderRVec::workVecSize);
//  workV[FirstOrderRVec::z].reset(new SiconosVector(sizeZ));
  workV[FirstOrderRVec::x].reset(new SiconosVector(sizeDS));
  workV[FirstOrderRVec::r].reset(new SiconosVector(sizeDS));

  workM.resize(FirstOrderRMat::workMatSize);
  if (!_C)
    workM[FirstOrderRMat::C].reset(new SimpleMatrix(sizeY, sizeDS));
  if (!_D)
    workM[FirstOrderRMat::D].reset(new SimpleMatrix(sizeY, sizeY));
//  if (!_jacgx)
//  {
//    workM[FirstOrderRMat::K].reset(new SimpleMatrix(sizeDS, sizeDS));
    // TODO add this back to workV of the DS -> needed for X partial NS
//  }
  if (!_B)
    workM[FirstOrderRMat::B].reset(new SimpleMatrix(sizeDS, sizeY));


  assert((_C->size(1) == sizeDS && _C->size(0) == sizeY) &&
         "FirstOrderType2R::initComponents inconsistent sizes between _jach[0] matrix and the interaction.");

  assert((_D->size(0) == sizeDS && _D->size(1) == sizeY) &&
         "FirstOrderType2R::initComponents inconsistent sizes between _jacg[0] matrix and the interaction.");
}

void FirstOrderType2R::computeh(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& y)
{
  ((Type2PtrH)(_pluginh->fPtr))(x.size(), x.getArray(), lambda.size(), lambda.getArray(), y.size(), y.getArray());
}

void FirstOrderType2R::computeg(double time, SiconosVector& lambda, SiconosVector& r)
{
  ((Type2PtrG)(_pluging->fPtr))(lambda.size(), lambda.getArray(), r.size(), r.getArray());
}

void FirstOrderType2R::computeg(double time, SiconosVector& lambda, SiconosVector& r);

void FirstOrderType2R::computeOutput(double time, Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM, SiconosMatrix& osnsM, unsigned int level)
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
  SiconosVector& y = *inter.y(level);
  DEBUG_EXPR(y.display());

  if (_D)
    prod(*_D, *(inter.lambdaOld(level)), y, true);
  else
    prod(*workM[FirstOrderRMat::D], *(inter.lambdaOld(level)), y, true);

  y *= -1.0;
  //SiconosVector yOld = *inter.yOld(0); // Retrieve  y_{alpha}_{k+1}
  DEBUG_PRINT("FirstOrderType2R::computeOutput : yOld(level) \n");
  DEBUG_EXPR(inter.yOld(level)->display());

  y += *inter.yOld(level);

  DEBUG_PRINT("FirstOrderType2R::computeOutput : ResiduY() \n");
  DEBUG_EXPR(inter.residuY()->display());

  y -= *inter.residuY();
  DEBUG_PRINT("FirstOrderType2R::computeOutput : y(level) \n");
  DEBUG_EXPR(y.display());

  BlockVector& deltax = *DSlink[FirstOrderRDS::deltax];
  //  deltax -= *(DSlink[FirstOrderRDS::xold)];
  DEBUG_PRINT("FirstOrderType2R::computeOutput : deltax \n");
  DEBUG_EXPR(deltax.display());

  if (_C)
    prod(*_C, deltax, y, false);
  else
    prod(*workM[FirstOrderRMat::C], deltax, y, false);


  prod(osnsM, *inter.lambda(level), y, false);
  DEBUG_PRINT("FirstOrderType2R::computeOutput : new linearized y \n");
  DEBUG_EXPR(y.display());

  SiconosVector& x = *workV[FirstOrderRVec::x];
  x = *DSlink[FirstOrderRDS::x];

  computeh(time, x, *inter.lambda(level), *inter.Halpha());

}

void FirstOrderType2R::computeInput(double time, Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM, SiconosMatrix& osnsM, unsigned int level)
{
  DEBUG_PRINT("FirstOrderType2R::computeInput \n");
  // compute the new r  obtained by linearisation
  // r_{alpha+1}_{k+1} = g(lambda_{k+1}^{alpha},t_k+1)
  //                     + B_{k+1}^alpha ( lambda_{k+1}^{alpha+1}- lambda_{k+1}^{alpha} )


  SiconosVector lambda = *inter.lambda(level);
  lambda -= *(inter.lambdaOld(level));

  //  std::cout<<"FirstOrderType2R::computeInput : diff lambda"<<endl;
  //  inter.lambdaOld(level)->display();
  //  lambda->display();
  //  _lambda->display();
  //  std::cout<<"FirstOrderType2R::computeInput : g_alpha"<<endl;
  //  _workX->display();
  if (_B)
    prod(*_B, lambda, *workV[FirstOrderRVec::g_alpha], false);
  else
    prod(*workM[FirstOrderRMat::B], lambda, *workV[FirstOrderRVec::g_alpha], false);
  //  std::cout<<"FirstOrderType2R::computeInput : result g_alpha - B*diffL"<<endl;
  //  _workX->display();

  *DSlink[FirstOrderRDS::r] += *workV[FirstOrderRVec::g_alpha];

  //compute the new g_alpha
  computeg(time, *inter.lambda(level), *workV[FirstOrderRVec::g_alpha]);



}

void FirstOrderType2R::preparNewtonIteration(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM)
{

  /* compute the contribution in xPartialNS for the next iteration */
  DEBUG_PRINT("FirstOrderType2R::preparNewtonIteration\n");
  if (_B)
    prod(*_B, *inter.lambda(0), *workV[FirstOrderRVec::x], true);
  else
    prod(*workM[FirstOrderRMat::B], *inter.lambda(0), *workV[FirstOrderRVec::x], true);

  *DSlink[FirstOrderRDS::xPartialNS] = *workV[FirstOrderRVec::g_alpha];
  *DSlink[FirstOrderRDS::xPartialNS] -= *workV[FirstOrderRVec::x];
}

void FirstOrderType2R::computeJachlambda(double time, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& D)
{
  RuntimeException::selfThrow("FirstOrderType2R::computeJachlambda must be overload.");
}
void FirstOrderType2R::computeJachx(double time, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& C)
{
  RuntimeException::selfThrow("FirstOrderType2R::computeJachx must be overload.");
}

void FirstOrderType2R::computeJach(double time, Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM)
{
  if (!_C)
  {
    SiconosVector& x = *workV[FirstOrderRVec::x];
    x = *DSlink[FirstOrderRDS::x];
    computeJachx(time, x, *inter.lambda(0), *workM[FirstOrderRMat::C]);
  }
  if (!_D)
  {
    SiconosVector& x = *workV[FirstOrderRVec::x];
    x = *DSlink[FirstOrderRDS::x];
    computeJachlambda(time, x, *inter.lambda(0), *workM[FirstOrderRMat::D]);
  }
}

void FirstOrderType2R::computeJacglambda(double time, SiconosVector& lambda, SimpleMatrix& B)
{
  RuntimeException::selfThrow("FirstOrderType2R::computeJacglambda must be overload.");
}

void FirstOrderType2R::computeJacg(double time, Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM)
{
  if (!_B)
  {
    computeJacglambda(time, *inter.lambda(0), *workM[FirstOrderRMat::B]);
  }
}
