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

// \todo : create a work vector for all tmp vectors used in computeG, computeH ...

#include "LagrangianCompliantR.h"
#include "RelationXML.h"
#include "Interaction.h"
#include "LagrangianDS.h"
#include "LagrangianR.cpp"

using namespace std;
using namespace RELATION;

// xml constructor
LagrangianCompliantR::LagrangianCompliantR(SP::RelationXML LRxml): BaseClass(LRxml, CompliantR)
{
  if (!LRxml->hasH())
    RuntimeException::selfThrow("LagrangianCompliantR:: xml constructor failed, can not find a definition for h.");

  hName = LRxml->getHPlugin();
  setComputeHFunction(SSL::getPluginName(hName), SSL::getPluginFunctionName(hName));

  if (!LRxml->hasJacobianH())
    RuntimeException::selfThrow("LagrangianCompliantR:: xml constructor failed, can not find a definition for JacH0.");
  JacH.resize(2);
  LRxml->readJacobianXML<PluggedMatrix, SP_PluggedMatrix>(JacH[0], LRxml, 0);
  LRxml->readJacobianXML<PluggedMatrix, SP_PluggedMatrix>(JacH[1], LRxml, 1);
}

// constructor from a set of data
LagrangianCompliantR::LagrangianCompliantR(const string& computeH, const std::vector<string> & computeG): BaseClass(CompliantR)
{
  setComputeHFunction(SSL::getPluginName(computeH), SSL::getPluginFunctionName(computeH));

  // Warning: we cannot allocate memory for JacH[0] matrix since no interaction
  // is connected to the relation. This will be done during initialize.
  // We only set the name of the plugin-function and connect it to the user-defined function.
  JacH.resize(2);
  assert(computeG.size() == 2 && "LagrangianCompliantR(string,string,vector<string> computeJH) error, computeJH size must be 2. ");
  JacH[0].reset(new PluggedMatrix(computeG[0]));
  JacH[1].reset(new PluggedMatrix(computeG[1]));
}

void LagrangianCompliantR::initComponents()
{
  BaseClass::initComponents();
  unsigned int sizeY = interaction->getSizeOfY();
  workL.reset(new SimpleVector(sizeY));

  if (! JacH[1])
    JacH[1].reset(new PluggedMatrix(sizeY, sizeY));
  else
    JacH[1]->resize(sizeY, sizeY);
}

void LagrangianCompliantR::computeH(double time)
{
  if (hPlugged)
  {
    // get vector y of the current interaction
    SP::SiconosVector y = interaction->getYPtr(0);
    SP::SiconosVector lambda = interaction->getLambdaPtr(0);

    // Warning: temporary method to have contiguous values in memory, copy of block to simple.
    *workX = *data[q0];
    *workZ = *data[z];
    *workY = *y;
    *workL = *lambda;

    unsigned int sizeQ = workX->size();
    unsigned int sizeY = y->size();
    unsigned int sizeZ = workZ->size();

    if (!hPtr)
      RuntimeException::selfThrow("LagrangianCompliantR:computeH() failed, h is not linked to a plugin function");
    hPtr(sizeQ, &(*workX)(0), sizeY, &(*workL)(0), &(*workY)(0), sizeZ, &(*workZ)(0));

    // Copy data that might have been changed in the plug-in call.
    *data[z] = *workZ;
    *y = *workY;
  }
}

void LagrangianCompliantR::computeJacH(double time, unsigned int  index)
{
  assert(index < JacH.size() && "LagrangianCompliantR:: computeJacH(index), index out of range.");

  if (JacH[index]->isPlugged())
  {
    // Warning: temporary method to have contiguous values in memory, copy of block to simple.
    *workX = *data[q0];
    *workZ = *data[z];

    unsigned int sizeY = JacH[0]->size(0);
    unsigned int sizeQ = workX->size();
    unsigned int sizeZ = workZ->size();

    // get vector lambda of the current interaction
    *workL = *interaction->getLambdaPtr(0);
    if (!(JacH[index]->fPtr))
      RuntimeException::selfThrow("LagrangianCompliantR::computeJacH(i) JacH[i] is not linked to a plugin function");
    (JacH[index]->fPtr)(sizeQ, &(*workX)(0), sizeY, &(*workL)(0), &(*JacH[index])(0, 0), sizeZ, &(*workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *data[z] = *workZ;
  }
}

void LagrangianCompliantR::computeOutput(double time, unsigned int derivativeNumber)
{
  if (derivativeNumber == 0)
    computeH(time);
  else
  {
    SP::SiconosVector y = interaction->getYPtr(derivativeNumber);
    computeJacH(time, 0);
    computeJacH(time, 1);
    if (derivativeNumber == 1)
    {
      SP::SiconosVector lambda = interaction->getLambdaPtr(derivativeNumber);
      // y = JacH[0] q1 + JacH[1] lambda
      prod(*JacH[0], *data[q1], *y);
      prod(*JacH[1], *lambda, *y, false);
    }
    else if (derivativeNumber == 2)
      prod(*JacH[0], *data[q2], *y); // Approx: y[2] = JacH[0]q[2], other terms are neglected ...
    else
      RuntimeException::selfThrow("LagrangianCompliantR::computeOutput(time,index), index out of range or not yet implemented.");
  }
}

void LagrangianCompliantR::computeInput(const double time, const unsigned int level)
{
  computeJacH(time, 0);
  // get lambda of the concerned interaction
  SP::SiconosVector lambda = interaction->getLambdaPtr(level);

  // data[name] += trans(G) * lambda
  prod(*lambda, *JacH[0], *data[p0 + level], false);

  //   SP::SiconosMatrix  GT = new SimpleMatrix(*G[0]);
  //   GT->trans();
  //   *data[name] += prod(*GT, *lambda);
  //  gemv(CblasTrans, 1.0,*(G[0]), *lambda, 1.0, *data[name]); => not yet implemented for BlockVectors.
}

LagrangianCompliantR* LagrangianCompliantR::convert(Relation *r)
{
  LagrangianCompliantR* lnlr = dynamic_cast<LagrangianCompliantR*>(r);
  return lnlr;
}

