/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#include "FirstOrderType2R.hpp"
#include "Interaction.hpp"
#include "FirstOrderNonLinearDS.hpp"

#include "BlockVector.hpp"
#include "SimulationGraphs.hpp"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES 1

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


//  workV[FirstOrderR::vec_z].reset(new SiconosVector(sizeZ));
  workV[FirstOrderR::vec_x].reset(new SiconosVector(sizeDS));
  workV[FirstOrderR::vec_r].reset(new SiconosVector(sizeDS));
  workV[FirstOrderR::h_alpha].reset(new SiconosVector(sizeY));
  workV[FirstOrderR::g_alpha].reset(new SiconosVector(sizeDS));

  if (!_C)
    workM[FirstOrderR::mat_C].reset(new SimpleMatrix(sizeY, sizeDS));
  if (!_D)
    workM[FirstOrderR::mat_D].reset(new SimpleMatrix(sizeY, sizeY));
//  if (!_jacgx)
//  {
//    workM[FirstOrderR::mat_K].reset(new SimpleMatrix(sizeDS, sizeDS));
    // TODO add this back to workV of the DS -> needed for X partial NS
//  }
  if (!_B)
    workM[FirstOrderR::mat_B].reset(new SimpleMatrix(sizeDS, sizeY));


}

void FirstOrderType2R::computeh(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& y)
{
  ((Type2PtrH)(_pluginh->fPtr))(x.size(), x.getArray(), lambda.size(), lambda.getArray(), y.size(), y.getArray());
}

void FirstOrderType2R::computeg(double time, SiconosVector& lambda, SiconosVector& r)
{
  ((Type2PtrG)(_pluging->fPtr))(lambda.size(), lambda.getArray(), r.size(), r.getArray());
}


void FirstOrderType2R::computeOutput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level)
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
  VectorOfBlockVectors& DSlink = *interProp.DSlink;
  VectorOfVectors& workV = *interProp.workVectors;
  VectorOfSMatrices& workM = *interProp.workMatrices;
  SiconosMatrix& osnsM = *interProp.block;

  if (_D)
    prod(*_D, *(inter.lambdaOld(level)), y, true);
  else
    prod(*workM[FirstOrderR::mat_D], *(inter.lambdaOld(level)), y, true);

  y *= -1.0;
  //SiconosVector yOld = *inter.yOld(0); // Retrieve  y_{alpha}_{k+1}
  DEBUG_PRINT("FirstOrderType2R::computeOutput : yOld(level) \n");
  DEBUG_EXPR(inter.yOld(level)->display());

  y += *inter.yOld(level);

  DEBUG_PRINT("FirstOrderType2R::computeOutput : ResiduY() \n");
  SiconosVector& residuY = *workV[FirstOrderR::vec_residuY];
  DEBUG_EXPR(residuY().display());

  y -= residuY;
  DEBUG_PRINT("FirstOrderType2R::computeOutput : y(level) \n");
  DEBUG_EXPR(y.display());

  BlockVector& deltax = *DSlink[FirstOrderR::deltax];
  //  deltax -= *(DSlink[FirstOrderR::xold)];
  DEBUG_PRINT("FirstOrderType2R::computeOutput : deltax \n");
  DEBUG_EXPR(deltax.display());

  if (_C)
    prod(*_C, deltax, y, false);
  else
    prod(*workM[FirstOrderR::mat_C], deltax, y, false);


  DEBUG_PRINT("FirstOrderType2R::computeOutput : y before osnsM\n");
  DEBUG_EXPR(y.display());
  prod(osnsM, *inter.lambda(level), y, false);
  DEBUG_PRINT("FirstOrderType2R::computeOutput : new linearized y \n");
  DEBUG_EXPR(y.display());

  SiconosVector& x = *workV[FirstOrderR::vec_x];
  x = *DSlink[FirstOrderR::x];

  SiconosVector& hAlpha= *workV[FirstOrderR::h_alpha];
  computeh(time, x, *inter.lambda(level), hAlpha);
  DEBUG_PRINT("FirstOrderType2R::computeOutput : new Halpha \n");
  DEBUG_EXPR(hAlpha()->display());

}

void FirstOrderType2R::computeInput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level)
{
  DEBUG_PRINT("FirstOrderType2R::computeInput \n");
  // compute the new r  obtained by linearisation
  // r_{alpha+1}_{k+1} = g(lambda_{k+1}^{alpha},t_k+1)
  //                     + B_{k+1}^alpha ( lambda_{k+1}^{alpha+1}- lambda_{k+1}^{alpha} )

  VectorOfBlockVectors& DSlink = *interProp.DSlink;
  VectorOfVectors& workV = *interProp.workVectors;
  VectorOfSMatrices& workM = *interProp.workMatrices;

  SiconosVector lambda = *inter.lambda(level);
  lambda -= *(inter.lambdaOld(level));

  //  std::cout<<"FirstOrderType2R::computeInput : diff lambda"<<endl;
  //  inter.lambdaOld(level)->display();
  //  lambda->display();
  //  _lambda->display();
  //  std::cout<<"FirstOrderType2R::computeInput : g_alpha"<<endl;
  //  _workX->display();
  if (_B)
    prod(*_B, lambda, *workV[FirstOrderR::g_alpha], false);
  else
    prod(*workM[FirstOrderR::mat_B], lambda, *workV[FirstOrderR::g_alpha], false);
  //  std::cout<<"FirstOrderType2R::computeInput : result g_alpha - B*diffL"<<endl;
  //  _workX->display();

  *DSlink[FirstOrderR::r] += *workV[FirstOrderR::g_alpha];

  //compute the new g_alpha
  computeg(time, *inter.lambda(level), *workV[FirstOrderR::g_alpha]);



}

void FirstOrderType2R::prepareNewtonIteration(Interaction& inter, InteractionProperties& interProp)
{

  /* compute the contribution in xPartialNS for the next iteration */
  VectorOfBlockVectors& DSlink = *interProp.DSlink;
  VectorOfVectors& workV = *interProp.workVectors;
  VectorOfSMatrices& workM = *interProp.workMatrices;
  DEBUG_PRINT("FirstOrderType2R::preparNewtonIteration\n");
  if (_B)
    prod(*_B, *inter.lambda(0), *workV[FirstOrderR::vec_x], true);
  else
    prod(*workM[FirstOrderR::mat_B], *inter.lambda(0), *workV[FirstOrderR::vec_x], true);

  *DSlink[FirstOrderR::xPartialNS] = *workV[FirstOrderR::g_alpha];
  *DSlink[FirstOrderR::xPartialNS] -= *workV[FirstOrderR::vec_x];
}

void FirstOrderType2R::computeJachlambda(double time, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& D)
{
  RuntimeException::selfThrow("FirstOrderType2R::computeJachlambda must be overload.");
}
void FirstOrderType2R::computeJachx(double time, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& C)
{
  RuntimeException::selfThrow("FirstOrderType2R::computeJachx must be overload.");
}

void FirstOrderType2R::computeJach(double time, Interaction& inter, InteractionProperties& interProp)
{
  VectorOfBlockVectors& DSlink = *interProp.DSlink;
  VectorOfVectors& workV = *interProp.workVectors;
  VectorOfSMatrices& workM = *interProp.workMatrices;
  if (!_C)
  {
    SiconosVector& x = *workV[FirstOrderR::vec_x];
    x = *DSlink[FirstOrderR::x];
    computeJachx(time, x, *inter.lambda(0), *workM[FirstOrderR::mat_C]);
  }
  if (!_D)
  {
    SiconosVector& x = *workV[FirstOrderR::vec_x];
    x = *DSlink[FirstOrderR::x];
    computeJachlambda(time, x, *inter.lambda(0), *workM[FirstOrderR::mat_D]);
  }
}

void FirstOrderType2R::computeJacglambda(double time, SiconosVector& lambda, SimpleMatrix& B)
{
  RuntimeException::selfThrow("FirstOrderType2R::computeJacglambda must be overload.");
}

void FirstOrderType2R::computeJacg(double time, Interaction& inter, InteractionProperties& interProp)
{
  if (!_B)
  {
    VectorOfSMatrices& workM = *interProp.workMatrices;
    computeJacglambda(time, *inter.lambda(0), *workM[FirstOrderR::mat_B]);
  }
}
