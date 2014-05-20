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
#include "FirstOrderNonLinearR.hpp"
#include "Interaction.hpp"
#include "FirstOrderNonLinearDS.hpp"

#include "BlockVector.hpp"
#include "SimulationTypeDef.hpp"

#include "debug.h"


void FirstOrderNonLinearR::initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM)
{
  unsigned int sizeY = inter.getSizeOfY();
  unsigned int sizeDS = inter.getSizeOfDS();
  workV.resize(FirstOrderR::workVecSize);
  workM.resize(FirstOrderR::mat_workMatSize);

  workV[FirstOrderR::vec_r].reset(new SiconosVector(sizeDS));
  workV[FirstOrderR::vec_x].reset(new SiconosVector(sizeDS));
  workV[FirstOrderR::g_alpha].reset(new SiconosVector(sizeDS));

  workM[FirstOrderR::mat_C].reset(new SimpleMatrix(sizeY, sizeDS));
  workM[FirstOrderR::mat_D].reset(new SimpleMatrix(sizeY, sizeY));
  workM[FirstOrderR::mat_B].reset(new SimpleMatrix(sizeDS, sizeY));
  workM[FirstOrderR::mat_K].reset(new SimpleMatrix(sizeDS, sizeDS));

  workM[FirstOrderR::mat_Khat].reset(new SimpleMatrix(sizeDS, sizeY));


}

void FirstOrderNonLinearR::computeJachx(double time, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& C)
{
  RuntimeException::selfThrow("FirstOrderNonLinearR::computeJachx, you need to derive this function in order to use it");
}
void FirstOrderNonLinearR::computeJachz(double time, Interaction& inter)
{
  RuntimeException::selfThrow("FirstOrderNonLinearR::computeJachx, you need to derive this function in order to use it");
}
void FirstOrderNonLinearR::computeJachlambda(double time, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& D)
{
  RuntimeException::selfThrow("FirstOrderNonLinearR::computeJachlambda, you need to derive this function in order to use it");
}

void FirstOrderNonLinearR::computeJacglambda(double time, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& B)
{
  RuntimeException::selfThrow("FirstOrderNonLinearR::computeJacglambda, you need to derive this function in order to use it");
}
void FirstOrderNonLinearR::computeJacgx(double time, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& K)
{
  RuntimeException::selfThrow("FirstOrderNonLinearR::computeJacglambda, you need to derive this function in order to use it");
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
  DEBUG_EXPR(inter.residuY()->display());

  y -= *inter.residuY();
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

  computeh(time, x, *inter.lambda(level), *inter.Halpha());

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
  computeg(time, x, *inter.lambda(level), g_alpha);
}

void FirstOrderNonLinearR::prepareNewtonIteration(Interaction& inter, InteractionProperties& interProp)
{

  /* compute the contribution in xPartialNS for the next iteration */
  DEBUG_PRINT("FirstOrderNonLinearR::preparNewtonIteration\n");

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
  SiconosVector& lambda = *inter.lambda(0);

  if (!_C)
  {
    computeJachx(time, x, lambda, *workM[FirstOrderR::mat_C]);
  }

  if (!_D)
  {
    computeJachlambda(time, x, lambda, *workM[FirstOrderR::mat_D]);
  }
}

void FirstOrderNonLinearR::computeJacg(double time, Interaction& inter, InteractionProperties& interProp)
{
  VectorOfBlockVectors& DSlink = *interProp.DSlink;
  VectorOfVectors& workV = *interProp.workVectors;
  VectorOfSMatrices& workM = *interProp.workMatrices;
  SiconosVector& x = *workV[FirstOrderR::vec_x];
  x = *DSlink[FirstOrderR::x];
  SiconosVector& lambda = *inter.lambda(0);
  if (!_B)
  {
    computeJacglambda(time, x, lambda, *workM[FirstOrderR::mat_B]);
  }
  if (!_K)
  {
    computeJacgx(time, x, lambda, *workM[FirstOrderR::mat_K]);
  }
}
