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

#include "debug.h"


void FirstOrderNonLinearR::initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM)
{
  unsigned int sizeY = inter.getSizeOfY();
  unsigned int sizeDS = inter.getSizeOfDS();
  workV.resize(FirstOrderRVec::workVecSize);
  workM.resize(FirstOrderRMat::workMatSize);

  workV[FirstOrderRVec::r].reset(new SiconosVector(sizeDS));
  workV[FirstOrderRVec::x].reset(new SiconosVector(sizeDS));
  workV[FirstOrderRVec::g_alpha].reset(new SiconosVector(sizeDS));

  workM[FirstOrderRMat::C].reset(new SimpleMatrix(sizeY, sizeDS));
  workM[FirstOrderRMat::D].reset(new SimpleMatrix(sizeY, sizeY));
  workM[FirstOrderRMat::B].reset(new SimpleMatrix(sizeDS, sizeY));
  workM[FirstOrderRMat::K].reset(new SimpleMatrix(sizeDS, sizeDS));
  workM[FirstOrderRMat::Khat].reset(new SimpleMatrix(sizeDS, sizeY));


}

void FirstOrderNonLinearR::computeJachx(double time, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& C)
{
  //RuntimeException::selfThrow("FirstOrderNonLinearR::computeJachx, not (yet) implemented or forbidden for relations of type "+subType);
}
void FirstOrderNonLinearR::computeJachz(double time, Interaction& inter)
{
  //RuntimeException::selfThrow("FirstOrderNonLinearR::computeJachx, not (yet) implemented or forbidden for relations of type "+subType);
}
void FirstOrderNonLinearR::computeJachlambda(double time, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& D)
{
  //RuntimeException::selfThrow("FirstOrderNonLinearR::computeJachlambda, not (yet) implemented or forbidden for relations of type "+subType);
}

void FirstOrderNonLinearR::computeJacglambda(double time, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& B)
{
  //RuntimeException::selfThrow("FirstOrderNonLinearR::computeJacglambda, not (yet) implemented or forbidden for relations of type "+subType);
}
void FirstOrderNonLinearR::computeJacgx(double time, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& K)
{
  //RuntimeException::selfThrow("FirstOrderNonLinearR::computeJacglambda, not (yet) implemented or forbidden for relations of type "+subType);
}




void FirstOrderNonLinearR::computeOutput(double time, Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM, SiconosMatrix& osnsM, unsigned int level)
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

  if (_D)
    prod(*_D, *(inter.lambdaOld(level)), y, true);
  else
    prod(*workM[FirstOrderRMat::D], *(inter.lambdaOld(level)), y, true);

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

  BlockVector& deltax = *DSlink[FirstOrderRDS::deltax];
  DEBUG_PRINT("FirstOrderNonLinearR::computeOutput : deltax \n");
  DEBUG_EXPR(deltax.display());

  if (_C)
    prod(*_C, deltax, y, false);
  else
    prod(*workM[FirstOrderRMat::C], deltax, y, false);

  // osnsM = h * C * W^-1 * B + D
  prod(osnsM, *inter.lambda(level), y, false);
  DEBUG_PRINT("FirstOrderNonLinearR::computeOutput : new linearized y \n");
  DEBUG_EXPR(y.display());

  SiconosVector& x = *workV[FirstOrderRVec::x];
  x = *DSlink[FirstOrderRDS::x];

  computeh(time, x, *inter.lambda(level), *inter.Halpha());

}

void FirstOrderNonLinearR::computeInput(double time, Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM, SiconosMatrix& osnsM, unsigned int level)
{
  DEBUG_PRINT("FirstOrderNonLinearR::computeInput \n");
  // compute the new r  obtained by linearisation
  // r_{alpha+1}_{k+1} = g(lambda_{k+1}^{alpha},t_k+1)
  //                     + B_{k+1}^alpha ( lambda_{k+1}^{alpha+1}- lambda_{k+1}^{alpha} )


  SiconosVector lambda = *inter.lambda(level);
  lambda -= *(inter.lambdaOld(level));

  SiconosVector& g_alpha = *workV[FirstOrderRVec::g_alpha];

  if (_B)
    prod(*_B, lambda, g_alpha, false);
  else
    prod(*workM[FirstOrderRMat::B], lambda, g_alpha, false);

  BlockVector& deltax = *DSlink[FirstOrderRDS::deltax];
  DEBUG_PRINT("FirstOrderNonLinearR::computeInput : deltax \n");
  DEBUG_EXPR(deltax.display());

  if (_K)
    prod(*_K, deltax, g_alpha, false);
  else
    prod(*workM[FirstOrderRMat::K], deltax, g_alpha, false);

  // Khat = h * K * W^-1 * B
  prod(*workM[FirstOrderRMat::Khat], *inter.lambda(level), g_alpha, false);

  *DSlink[FirstOrderRDS::r] += g_alpha;

  //compute the new g_alpha
  SiconosVector& x = *workV[FirstOrderRVec::x];
  x = *DSlink[FirstOrderRDS::x];
  computeg(time, x, *inter.lambda(level), g_alpha);
}

void FirstOrderNonLinearR::preparNewtonIteration(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM)
{

  /* compute the contribution in xPartialNS for the next iteration */
  DEBUG_PRINT("FirstOrderNonLinearR::preparNewtonIteration\n");

  BlockVector& xPartialNS = *DSlink[FirstOrderRDS::xPartialNS];
  xPartialNS = *workV[FirstOrderRVec::g_alpha];

  if (_B)
    prod(*_B, *inter.lambda(0), *workV[FirstOrderRVec::x], true);
  else
    prod(*workM[FirstOrderRMat::B], *inter.lambda(0), *workV[FirstOrderRVec::x], true);

  xPartialNS -= *workV[FirstOrderRVec::x];

}


