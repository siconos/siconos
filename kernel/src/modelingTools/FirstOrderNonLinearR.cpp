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
#include "FirstOrderNonLinearR.hpp"
#include "Interaction.hpp"
#include "FirstOrderNonLinearDS.hpp"

#include "BlockVector.hpp"
#include "SimulationGraphs.hpp"

#include "debug.h"

typedef void (*FONLR_h)(double, unsigned, double*, unsigned, double*, double*, unsigned, double*);
typedef FONLR_h FONLR_g;
typedef FONLR_h FONLR_C;
typedef FONLR_h FONLR_B;
typedef FONLR_h FONLR_K;
typedef FONLR_h FONLR_D;



void FirstOrderNonLinearR::initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM)
{
  unsigned int sizeY = inter.getSizeOfY();
  unsigned int sizeDS = inter.getSizeOfDS();
  unsigned int sizeZ = DSlink[FirstOrderR::z]->size();
  workV.resize(FirstOrderR::workVecSize);
  workM.resize(FirstOrderR::mat_workMatSize);

  workV[FirstOrderR::vec_r].reset(new SiconosVector(sizeDS));
  workV[FirstOrderR::vec_x].reset(new SiconosVector(sizeDS));
  workV[FirstOrderR::vec_z].reset(new SiconosVector(sizeZ));
  
  workV[FirstOrderR::h_alpha].reset(new SiconosVector(sizeY));
  workV[FirstOrderR::g_alpha].reset(new SiconosVector(sizeDS));

  workM[FirstOrderR::mat_C].reset(new SimpleMatrix(sizeY, sizeDS));
  workM[FirstOrderR::mat_D].reset(new SimpleMatrix(sizeY, sizeY));
  workM[FirstOrderR::mat_B].reset(new SimpleMatrix(sizeDS, sizeY));
  workM[FirstOrderR::mat_K].reset(new SimpleMatrix(sizeDS, sizeDS));

  workM[FirstOrderR::mat_Khat].reset(new SimpleMatrix(sizeDS, sizeY));


}

void FirstOrderNonLinearR::computeh(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SiconosVector& y)
{
  if (_pluginh)
  {
    ((FONLR_h)_pluginh->fPtr)(time, x.size(), x.getArray(), lambda.size(), lambda.getArray(), y.getArray(), z.size(), z.getArray());
  }
  else
  {
    RuntimeException::selfThrow("FirstOrderNonLinearR::computeh - no plugin detected, you should provide one or derive this class and implement this function");
  }
}

void FirstOrderNonLinearR::computeg(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SiconosVector& r)
{
  if (_pluging)
  {
    ((FONLR_g)_pluging->fPtr)(time, x.size(), x.getArray(), lambda.size(), lambda.getArray(), r.getArray(), z.size(), z.getArray());
  }
  else
  {
    RuntimeException::selfThrow("FirstOrderNonLinearR::computeg - no plugin detected, you should provide one or derive this class and implement this function");
  }
}

void FirstOrderNonLinearR::computeJachx(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SimpleMatrix& C)
{
  if (_pluginJachx)
  {
    ((FONLR_C)_pluginJachx->fPtr)(time, x.size(), x.getArray(), lambda.size(), lambda.getArray(), C.getArray(), z.size(), z.getArray());
  }
  else
    RuntimeException::selfThrow("FirstOrderNonLinearR::computeJachx, you need to derive this function in order to use it");
}

void FirstOrderNonLinearR::computeJachlambda(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SimpleMatrix& D)
{
  if (_pluginJachlambda)
  {
    ((FONLR_D)_pluginJachlambda->fPtr)(time, x.size(), x.getArray(), lambda.size(), lambda.getArray(), D.getArray(), z.size(), z.getArray());
  }
  else
    RuntimeException::selfThrow("FirstOrderNonLinearR::computeJachlambda, you need to either provide a matrix D or derive this function in order to use it");
}

void FirstOrderNonLinearR::computeJacglambda(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SimpleMatrix& B)
{
  if (_pluginJacglambda)
  {
    ((FONLR_B)_pluginJacglambda->fPtr)(time, x.size(), x.getArray(), lambda.size(), lambda.getArray(), B.getArray(), z.size(), z.getArray());
  }
  else
    RuntimeException::selfThrow("FirstOrderNonLinearR::computeJacglambda, you need to either provide a matrix B or derive this function in order to use it");
}

void FirstOrderNonLinearR::computeJacgx(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SimpleMatrix& K)
{
  if (_pluginJacgx)
  {
    ((FONLR_K)_pluginJacgx->fPtr)(time, x.size(), x.getArray(), lambda.size(), lambda.getArray(), K.getArray(), z.size(), z.getArray());
  }
  else
    RuntimeException::selfThrow("FirstOrderNonLinearR::computeJacgx, you need to either provide a matrix K or derive this function in order to use it");
}




void FirstOrderNonLinearR::computeOutput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level)
{
  DEBUG_PRINT("FirstOrderNonLinearR::computeOutput \n");
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
  DEBUG_PRINT("FirstOrderNonLinearR::computeOutput : yOld(level) \n");
  DEBUG_EXPR(inter.yOld(level)->display());

  y += *inter.yOld(level);

  DEBUG_PRINT("FirstOrderNonLinearR::computeOutput : ResiduY() \n");
  SiconosVector& residuY = *workV[FirstOrderR::vec_residuY];
  DEBUG_EXPR(residuY().display());

  y -= residuY;
  
  DEBUG_PRINT("FirstOrderNonLinearR::computeOutput : y(level) \n");
  DEBUG_EXPR(y.display());

  BlockVector& deltax = *DSlink[FirstOrderR::deltax];
  DEBUG_PRINT("FirstOrderNonLinearR::computeOutput : deltax \n");
  DEBUG_EXPR(deltax.display());

  if (_C)
    prod(*_C, deltax, y, false);
  else
    prod(*workM[FirstOrderR::mat_C], deltax, y, false);

  // osnsM = h * C * W^-1 * B + D
  prod(osnsM, *inter.lambda(level), y, false);
  DEBUG_PRINT("FirstOrderNonLinearR::computeOutput : new linearized y \n");
  DEBUG_EXPR(y.display());

  SiconosVector& x = *workV[FirstOrderR::vec_x];
  x = *DSlink[FirstOrderR::x];
  SiconosVector& z = *workV[FirstOrderR::vec_z];
  z = *DSlink[FirstOrderR::z];

  SiconosVector hAlpha =  *workV[FirstOrderR::h_alpha];
  computeh(time, x, *inter.lambda(level), z, hAlpha);

  *DSlink[FirstOrderR::z] = z;
}

void FirstOrderNonLinearR::computeInput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level)
{
  DEBUG_PRINT("FirstOrderNonLinearR::computeInput \n");
  // compute the new r  obtained by linearisation
  // r_{alpha+1}_{k+1} = g(lambda_{k+1}^{alpha},t_k+1)
  //                     + B_{k+1}^alpha ( lambda_{k+1}^{alpha+1}- lambda_{k+1}^{alpha} )


  SiconosVector lambda = *inter.lambda(level);
  lambda -= *(inter.lambdaOld(level));

  VectorOfBlockVectors& DSlink = *interProp.DSlink;
  VectorOfVectors& workV = *interProp.workVectors;
  VectorOfSMatrices& workM = *interProp.workMatrices;

  SiconosVector& g_alpha = *workV[FirstOrderR::g_alpha];

  if (_B)
    prod(*_B, lambda, g_alpha, false);
  else
    prod(*workM[FirstOrderR::mat_B], lambda, g_alpha, false);

  BlockVector& deltax = *DSlink[FirstOrderR::deltax];
  DEBUG_PRINT("FirstOrderNonLinearR::computeInput : deltax \n");
  DEBUG_EXPR(deltax.display());

  if (_K)
    prod(*_K, deltax, g_alpha, false);
  else
    prod(*workM[FirstOrderR::mat_K], deltax, g_alpha, false);


  // Khat = h * K * W^-1 * B
  prod(*workM[FirstOrderR::mat_Khat], *inter.lambda(level), g_alpha, false);

  *DSlink[FirstOrderR::r] += g_alpha;

  //compute the new g_alpha
  SiconosVector& x = *workV[FirstOrderR::vec_x];
  x = *DSlink[FirstOrderR::x];
  SiconosVector& z = *workV[FirstOrderR::vec_z];
  z = *DSlink[FirstOrderR::z];
  computeg(time, x, *inter.lambda(level), z, g_alpha);
  *DSlink[FirstOrderR::z] = z;
}

void FirstOrderNonLinearR::prepareNewtonIteration(Interaction& inter, InteractionProperties& interProp)
{

  /* compute the contribution in xPartialNS for the next iteration */
  DEBUG_PRINT("FirstOrderNonLinearR::prepareNewtonIteration\n");

  VectorOfBlockVectors& DSlink = *interProp.DSlink;
  VectorOfVectors& workV = *interProp.workVectors;
  VectorOfSMatrices& workM = *interProp.workMatrices;

  BlockVector& xPartialNS = *DSlink[FirstOrderR::xPartialNS];
  xPartialNS = *workV[FirstOrderR::g_alpha];

  if (_B)
    prod(*_B, *inter.lambda(0), *workV[FirstOrderR::vec_x], true);
  else
    prod(*workM[FirstOrderR::mat_B], *inter.lambda(0), *workV[FirstOrderR::vec_x], true);

  xPartialNS -= *workV[FirstOrderR::vec_x];


}

void FirstOrderNonLinearR::computeJach(double time, Interaction& inter, InteractionProperties& interProp)
{
  VectorOfBlockVectors& DSlink = *interProp.DSlink;
  VectorOfVectors& workV = *interProp.workVectors;
  VectorOfSMatrices& workM = *interProp.workMatrices;

  SiconosVector& x = *workV[FirstOrderR::vec_x];
  x = *DSlink[FirstOrderR::x];
  SiconosVector& z = *workV[FirstOrderR::vec_z];
  z =  *DSlink[FirstOrderR::z];
  SiconosVector& lambda = *inter.lambda(0);

  if (!_C)
  {
    computeJachx(time, x, lambda, z, *workM[FirstOrderR::mat_C]);
  }

  if (!_D)
  {
    computeJachlambda(time, x, lambda, z, *workM[FirstOrderR::mat_D]);
  }
  *DSlink[FirstOrderR::z] = z;
}

void FirstOrderNonLinearR::computeJacg(double time, Interaction& inter, InteractionProperties& interProp)
{
  VectorOfBlockVectors& DSlink = *interProp.DSlink;
  VectorOfVectors& workV = *interProp.workVectors;
  VectorOfSMatrices& workM = *interProp.workMatrices;
  SiconosVector& x = *workV[FirstOrderR::vec_x];
  x = *DSlink[FirstOrderR::x];
  SiconosVector& z = *workV[FirstOrderR::vec_z];
  z =  *DSlink[FirstOrderR::z];
  SiconosVector& lambda = *inter.lambda(0);
  if (!_B)
  {
    computeJacglambda(time, x, lambda, z, *workM[FirstOrderR::mat_B]);
  }
  if (!_K)
  {
    computeJacgx(time, x, lambda, z, *workM[FirstOrderR::mat_K]);
  }
  *DSlink[FirstOrderR::z] = z;
}
