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
#include "FirstOrderLinearR.hpp"
#include "Interaction.hpp"
#include "PluggedObject.hpp"
#include "SimpleMatrix.hpp"
#include "BlockVector.hpp"
#include "SimulationGraphs.hpp"

#include <iostream>

using namespace RELATION;

FirstOrderLinearR::FirstOrderLinearR():
  FirstOrderR(LinearR)
{
  ;
}

// Constructor with C and B plug-in names
FirstOrderLinearR::FirstOrderLinearR(const std::string& Cname, const std::string& Bname):
  FirstOrderR(LinearR)
{
  // Warning: we cannot allocate memory for C/D matrix since no interaction
  // is connected to the relation. This will be done during initialize.
  // We only set the name of the plugin-function and connect it to the user-defined function.
  _pluginJachx->setComputeFunction(Cname);
  _pluginJacglambda->setComputeFunction(Bname);
}

// Constructor from a complete set of data (plugin)
FirstOrderLinearR::FirstOrderLinearR(const std::string& Cname, const std::string& Dname, const std::string& Fname, const std::string& Ename, const std::string& Bname): FirstOrderR(LinearR)
{
  _pluginJachx->setComputeFunction(Cname);
  _pluginJachlambda->setComputeFunction(Dname);
  _pluginJacglambda->setComputeFunction(Bname);
  _pluginf->setComputeFunction(Fname);
  _plugine->setComputeFunction(Ename);
}

// Minimum data (C, B as pointers) constructor
FirstOrderLinearR::FirstOrderLinearR(SP::SimpleMatrix C, SP::SimpleMatrix B):
  FirstOrderR(LinearR)
{
  _C = C;
  _B = B;
}

// // Constructor from a complete set of data
FirstOrderLinearR::FirstOrderLinearR(SP::SimpleMatrix C, SP::SimpleMatrix D, SP::SimpleMatrix F, SP::SiconosVector E, SP::SimpleMatrix B):
  FirstOrderR(LinearR)
{
  _C = C;
  _B = B;
  _D = D;
  _F = F;
  _e = E;
}

void FirstOrderLinearR::initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM)
{
  // Note: do not call FirstOrderR::initialize to avoid jacobianH and jacobianG allocation.

  // get interesting size
  unsigned int sizeY = inter.getSizeOfY();
  unsigned int sizeX = inter.getSizeOfDS();
  unsigned int sizeZ = DSlink[FirstOrderR::z]->size();

  // Update workV (copy of DS variables)
  workV.resize(FirstOrderR::workVecSize);
  workV[FirstOrderR::vec_z].reset(new SiconosVector(sizeZ));

  workM.resize(FirstOrderR::mat_workMatSize);

  if (!_C && _pluginJachx->fPtr)
    workM[FirstOrderR::mat_C].reset(new SimpleMatrix(sizeY, sizeX));
  if (!_D && _pluginJachlambda->fPtr)
    workM[FirstOrderR::mat_D].reset(new SimpleMatrix(sizeY, sizeY));
  if (!_B && _pluginJacglambda->fPtr)
    workM[FirstOrderR::mat_B].reset(new SimpleMatrix(sizeX, sizeY));
  if (!_F && _pluginf->fPtr)
  {
    workM[FirstOrderR::mat_F].reset(new SimpleMatrix(sizeY, sizeZ));
  }
  if (!_e && _plugine->fPtr)
  {
    workV[FirstOrderR::e].reset(new SiconosVector(sizeY));
  }

  // Check if various operators sizes are consistent.
  // Reference: interaction.

  if (_C)
  {
    if (_C->size(0) == 0)
      _C->resize(sizeX, sizeY);
    else
      assert((_C->size(0) == sizeY && _C->size(1) == sizeX) && "FirstOrderLinearR::initialize , inconsistent size between C and Interaction.");
  }
  if (_B)
  {
    if (_B->size(0) == 0)
    _B->resize(sizeY, sizeX);
  else
    assert((_B->size(1) == sizeY && _B->size(0) == sizeX) && "FirstOrderLinearR::initialize , inconsistent size between B and interaction.");
  }
  // C and B are the minimum inputs. The others may remain null.

  if (_D)
  {
    if (_D->size(0) == 0)
      _D->resize(sizeY, sizeY);
    else
      assert((_D->size(0) == sizeY || _D->size(1) == sizeY) && "FirstOrderLinearR::initialize , inconsistent size between C and D.");
  }

  if (_F)
  {
    if (_F->size(0) == 0)
      _F->resize(sizeY, sizeZ);
    else
      assert(((_F->size(0) == sizeY) && (_F->size(1) == sizeZ)) && "FirstOrderLinearR::initialize , inconsistent size between C and F.");
  }

  if (_e)
  {
    if (_e->size() == 0)
      _e->resize(sizeY);
    else
      assert(_e->size() == sizeY && "FirstOrderLinearR::initialize , inconsistent size between C and e.");
  }
}

void FirstOrderLinearR::computeC(double time, SiconosVector& z, SimpleMatrix& C)
{
  if (_pluginJachx->fPtr)
  {
    ((FOMatPtr1)(_pluginJachx->fPtr))(time, C.size(0), C.size(1), &(C)(0, 0), z.size(), &(z)(0));
  }
}

void FirstOrderLinearR::computeD(double time, SiconosVector& z, SimpleMatrix& D)
{
  if (_pluginJachlambda->fPtr)
  {
    ((FOMatPtr1)(_pluginJachlambda->fPtr))(time, D.size(0), D.size(1), &(D)(0, 0), z.size(), &(z)(0));
  }
}

void FirstOrderLinearR::computeF(double time, SiconosVector& z, SimpleMatrix& F)
{
  if (_pluginf->fPtr)
  {
    ((FOMatPtr1)(_pluginf->fPtr))(time, F.size(0), F.size(1), &(F)(0, 0), z.size(), &(z)(0));
  }
}

void FirstOrderLinearR::computee(double time, SiconosVector& z, SiconosVector& e)
{

  if (_plugine->fPtr)
  {
    ((FOVecPtr) _plugine->fPtr)(time, e.size(), &(e)(0), z.size(), &(z)(0));
  }
}

void FirstOrderLinearR::computeB(double time, SiconosVector& z, SimpleMatrix& B)
{
  if (_pluginJacglambda->fPtr)
  {
    ((FOMatPtr1) _pluginJacglambda->fPtr)(time, B.size(0), B.size(1), &(B)(0, 0), z.size(), &(z)(0));
  }
}

void FirstOrderLinearR::computeh(double time, VectorOfVectors& workV, VectorOfSMatrices& workM, BlockVector& x, SiconosVector& lambda, SiconosVector& z, SiconosVector& y)
{

  y.zero();

  if (_pluginJachx->fPtr)
  {
    SimpleMatrix& C = *workM[FirstOrderR::mat_C];
    computeC(time, z, C);
    prod(C, x, y, false);
  }
  if (_pluginJachlambda->fPtr)
  {
    SimpleMatrix& D = *workM[FirstOrderR::mat_D];
    computeD(time, z, D);
    prod(D, lambda, y, false);
  }
  if (_pluginf->fPtr)
  {
    SimpleMatrix& F = *workM[FirstOrderR::mat_F];
    computeF(time, z, F);
    prod(F, z, y, false);
  }
  if (_plugine->fPtr)
  {
    SiconosVector& e = *workV[FirstOrderR::e];
    computee(time, z, e);
    y+= e;
  }

  if (_C)
    prod(*_C, x, y, false);

  if (_D)
    prod(*_D, lambda, y, false);

  if (_e)
    y += *_e;

  if (_F)
    prod(*_F, z, y, false);

}

void FirstOrderLinearR::computeOutput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level)
{
  VectorOfBlockVectors& DSlink = *interProp.DSlink;
  VectorOfVectors& workV = *interProp.workVectors;
  VectorOfSMatrices& workM = *interProp.workMatrices;
  SiconosVector& z = *workV[FirstOrderR::vec_z];
  z = *DSlink[FirstOrderR::z];
  // We get y and lambda of the interaction (pointers)
  SiconosVector& y = *inter.y(0);
  SiconosVector& lambda = *inter.lambda(0);

  computeh(time, workV, workM, *DSlink[FirstOrderR::x], lambda, z, y);

  *DSlink[FirstOrderR::z] = z;
}

void FirstOrderLinearR::computeg(double time, VectorOfSMatrices& workM, SiconosVector& lambda, SiconosVector& z, BlockVector& r)
{
  if (_B)
  {
    prod(*_B, lambda, r, false);
  }
  else
  {
    SimpleMatrix& B = *workM[FirstOrderR::mat_B];
    computeB(time, z, B);
    prod(B, lambda, r, false);
  }
}


void FirstOrderLinearR::computeInput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level)
{
  SiconosVector& lambda = *inter.lambda(level);
  VectorOfBlockVectors& DSlink = *interProp.DSlink;
  VectorOfVectors& workV = *interProp.workVectors;
  VectorOfSMatrices& workM = *interProp.workMatrices;
  SiconosVector& z = *workV[FirstOrderR::z];
  z = *DSlink[FirstOrderR::z];

  computeg(time, workM, lambda, z, *DSlink[FirstOrderR::r]);

  *DSlink[FirstOrderR::z] = z;
}

void FirstOrderLinearR::display() const
{
  std::cout << " ===== Linear Time Invariant relation display ===== " <<std::endl;
  std::cout << "| C " <<std::endl;
  if (_C) _C->display();
  else std::cout << "->NULL" <<std::endl;
  std::cout << "| D " <<std::endl;
  if (_D) _D->display();
  else std::cout << "->NULL" <<std::endl;
  std::cout << "| F " <<std::endl;
  if (_F) _F->display();
  else std::cout << "->NULL" <<std::endl;
  std::cout << "| e " <<std::endl;
  if (_e) _e->display();
  else std::cout << "->NULL" <<std::endl;
  std::cout << "| B " <<std::endl;
  if (_B) _B->display();
  else std::cout << "->NULL" <<std::endl;
  std::cout << " ================================================== " <<std::endl;
}
