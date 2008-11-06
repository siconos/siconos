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
#include "FirstOrderLinearR.h"
#include "LinearRXML.h"
#include "Interaction.h"
#include "FirstOrderR.cpp"

using namespace std;
using namespace RELATION;

// xml constructor
FirstOrderLinearR::FirstOrderLinearR(SP::RelationXML relxml):
  FirstOrderR<MatrixFunctionOfTime>(relxml, LinearR)
{
  SP::LinearRXML folrXML = boost::static_pointer_cast<LinearRXML>(relationxml);
  // get matrices values. All are optional.

  JacH.reserve(2);
  JacG.resize(1);
  string plugin;
  if (folrXML->hasC())
  {
    JacH.resize(1);
    if (folrXML->isCPlugin())
      JacH[0].reset(new PluggedMatrix(folrXML->getCPlugin()));
    else
      JacH[0].reset(new PluggedMatrix(folrXML->getC()));
  }
  else
    RuntimeException::selfThrow("FirstOrderLinearR:: xml constructor failed, can not find a definition for C.");

  if (folrXML->hasD())
  {
    JacH.resize(2);
    if (folrXML->isDPlugin())
      JacH[1].reset(new PluggedMatrix(folrXML->getDPlugin()));
    else
      JacH[1].reset(new PluggedMatrix(folrXML->getD()));
  }

  if (folrXML->hasF())
  {
    if (folrXML->isFPlugin())
      F.reset(new PluggedMatrix(folrXML->getFPlugin()));
    else
      F.reset(new PluggedMatrix(folrXML->getF()));
  }

  if (folrXML->hasE())
  {
    if (folrXML->isEPlugin())
      e.reset(new PVTime(folrXML->getEPlugin()));
    else
      e.reset(new PVTime(folrXML->getE()));
  }

  if (folrXML->hasB())
  {
    if (folrXML->isBPlugin())
      JacG[0].reset(new PluggedMatrix(folrXML->getBPlugin()));
    else
      JacG[0].reset(new PluggedMatrix(folrXML->getB()));
  }
  else
    RuntimeException::selfThrow("FirstOrderLinearR:: xml constructor failed, can not find a definition for B.");
}

// Constructor with C and B plug-in names
FirstOrderLinearR::FirstOrderLinearR(const string& CName, const string& BName):
  FirstOrderR<MatrixFunctionOfTime>(LinearR)
{
  // Warning: we cannot allocate memory for C/D matrix since no interaction
  // is connected to the relation. This will be done during initialize.
  // We only set the name of the plugin-function and connect it to the user-defined function.
  JacH.reserve(2);
  JacG.resize(1);
  JacH.resize(1);
  JacH[0].reset(new PluggedMatrix(CName));
  JacG[0].reset(new PluggedMatrix(BName));
}

// Constructor from a complete set of data (plugin)
FirstOrderLinearR::FirstOrderLinearR(const string& CName, const string& DName, const string& FName, const string& EName, const string& BName): FirstOrderR<MatrixFunctionOfTime>(LinearR)
{
  JacG.resize(1);
  JacH.resize(2);
  JacH[0].reset(new PluggedMatrix(CName));
  JacH[1].reset(new PluggedMatrix(DName));
  F.reset(new PluggedMatrix(FName));
  e.reset(new PVTime(EName));
  JacG[0].reset(new PluggedMatrix(BName));
}

// Minimum data (C, B as pointers) constructor
FirstOrderLinearR::FirstOrderLinearR(SP_PluggedMatrix newC, SP_PluggedMatrix newB):
  FirstOrderR<MatrixFunctionOfTime>(LinearR)
{
  JacH.push_back(newC);
  JacG.push_back(newB);
}

// Constructor from a complete set of data
FirstOrderLinearR::FirstOrderLinearR(SP_PluggedMatrix newC, SP_PluggedMatrix newD, SP_PluggedMatrix newF, SP::PVTime newE, SP_PluggedMatrix newB):
  FirstOrderR<MatrixFunctionOfTime>(LinearR)
{
  JacH.push_back(newC);
  JacG.push_back(newB);
  JacH.push_back(newD);
  F = newF;
  e = newE;
}

// Minimum data (C, B as matrices) constructor
FirstOrderLinearR::FirstOrderLinearR(const SiconosMatrix& newC, const SiconosMatrix& newB):
  FirstOrderR<MatrixFunctionOfTime>(LinearR)
{
  JacH.reserve(2);
  JacG.resize(1);
  JacH.resize(1);
  JacH[0].reset(new PluggedMatrix(newC));
  JacG[0].reset(new PluggedMatrix(newB));
}

// Constructor from a complete set of data (matrices)
FirstOrderLinearR::FirstOrderLinearR(const SiconosMatrix& newC, const SiconosMatrix& newD, const SiconosMatrix& newF, const SiconosVector& newE, const SiconosMatrix& newB):
  FirstOrderR<MatrixFunctionOfTime>(LinearR)
{
  JacG.resize(1);
  JacH.resize(2);
  JacH[0].reset(new PluggedMatrix(newC));
  JacH[1].reset(new PluggedMatrix(newD));
  F.reset(new PluggedMatrix(newF));
  e.reset(new PVTime(newE));
  JacG[0].reset(new PluggedMatrix(newB));
}

void FirstOrderLinearR::initialize(SP::Interaction inter)
{
  assert(inter && "FirstOrderLinearR::initialize failed. No Interaction linked to the present relation.");
  interaction = inter;

  // Note: do not call FirstOrderR::initialize to avoid jacobianH and jacobianG allocation.

  // Update data member (links to DS variables)
  initDSLinks();
  if (!JacH[0])
    RuntimeException::selfThrow("FirstOrderLinearR::initialize() C is null and is a required input.");
  if (!JacG[0])
    RuntimeException::selfThrow("FirstOrderLinearR::initialize() B is null and is a required input.");

  // Check if various operators sizes are consistent.
  // Reference: interaction.
  unsigned int sizeY = interaction->getSizeOfY();
  unsigned int sizeX = interaction->getSizeOfDS();
  unsigned int sizeZ = interaction->getSizeZ();

  // The initialization of each matrix/vector depends on the way the Relation was built ie if the matrix/vector
  // was read from xml or not
  if (JacH[0]->size(0) == 0)
    JacH[0]->resize(sizeX, sizeY);
  else
    assert((JacH[0]->size(0) == sizeY && JacH[0]->size(1) == sizeX) && "FirstOrderLinearR::initialize , inconsistent size between C and Interaction.");

  if (JacG[0]->size(0) == 0)
    JacG[0]->resize(sizeY, sizeX);
  else
    assert((JacG[0]->size(1) == sizeY && JacG[0]->size(0) == sizeX) && "FirstOrderLinearR::initialize , inconsistent size between B and interaction.");

  // C and B are the minimum inputs. The others may remain null.

  if (JacH[1])
  {
    if (JacH[1]->size(0) == 0)
      JacH[1]->resize(sizeY, sizeY);
    else
      assert((JacH[1]->size(0) == sizeY || JacH[1]->size(1) == sizeY) && "FirstOrderLinearR::initialize , inconsistent size between C and D.");
  }

  if (F)
  {
    if (F->size(0) == 0)
      F->resize(sizeY, sizeZ);
    else
      assert(((F->size(0) != sizeY) && (F->size(1) != sizeZ)) && "FirstOrderLinearR::initialize , inconsistent size between C and F.");
  }

  if (e)
  {
    if (e->size() == 0)
      e->resize(sizeY);
    else
      assert(e->size() == sizeY && "FirstOrderLinearR::initialize , inconsistent size between C and e.");
  }

  workZ.reset(new SimpleVector(sizeZ));
}

// setters

void FirstOrderLinearR::setComputeCFunction(const string& pluginPath, const string& functionName)
{
  setComputeJacobianHFunction(pluginPath, functionName, 0);
}

void FirstOrderLinearR::setComputeDFunction(const string& pluginPath, const string& functionName)
{
  if (JacH.size() < 2) JacH.resize(2);
  if (!JacH[1])
    JacH[1].reset(new PluggedMatrix()); // Required since D may not have been created by constructor
  setComputeJacobianHFunction(pluginPath, functionName, 1);
}

void FirstOrderLinearR::setComputeFFunction(const string& pluginPath, const string& functionName)
{
  if (!F)
    F.reset(new PluggedMatrix()); // Required since F may not have been created by constructor
  F->setComputeFunction(pluginPath, functionName);
}

void FirstOrderLinearR::setComputeEFunction(VectorFunctionOfTime ptrFunct)
{
  if (!e)
    e.reset(new PVTime()); // Required since e may not have been created by constructor
  e->setComputeFunction(ptrFunct);
}

void FirstOrderLinearR::setComputeEFunction(const string& pluginPath, const string& functionName)
{
  if (!e)
    e.reset(new PVTime()); // Required since e may not have been created by constructor
  e->setComputeFunction(pluginPath, functionName);
}

void FirstOrderLinearR::setComputeBFunction(const string& pluginPath, const string& functionName)
{
  setComputeJacobianGFunction(pluginPath, functionName, 0);
}

void FirstOrderLinearR::computeC(const double time)
{
  if (JacH[0])
  {
    if (JacH[0]->isPlugged())
    {
      if (!JacH[0]->fPtr)
        RuntimeException::selfThrow("FirstOrderLinearR::computeC() is not linked to a plugin function");
      unsigned int sizeY = interaction->getSizeOfY();
      unsigned int sizeX = interaction->getSizeOfDS();
      unsigned int sizeZ = interaction->getSizeZ();
      *workZ = *data[z];
      (JacH[0]->fPtr)(time, sizeY, sizeX, &(*JacH[0])(0, 0), sizeZ, &(*workZ)(0));
      // Copy data that might have been changed in the plug-in call.
      *data[z] = *workZ;
    }
    // else nothing
  }
}

void FirstOrderLinearR::computeD(const double time)
{
  if (JacH[1])
  {
    if (JacH[1]->isPlugged())
    {
      if (!JacH[1]->fPtr)
        RuntimeException::selfThrow("FirstOrderLinearR::computeD() is not linked to a plugin function");
      unsigned int sizeY = interaction->getSizeOfY();
      unsigned int sizeZ = interaction->getSizeZ();
      *workZ = *data[z];
      JacH[1]->fPtr(time, sizeY, sizeY, &(*JacH[1])(0, 0), sizeZ, &(*workZ)(0));
      // Copy data that might have been changed in the plug-in call.
      *data[z] = *workZ;
    }
    // else nothing
  }
}

void FirstOrderLinearR::computeF(const double time)
{
  if (F)
  {
    if (F->isPlugged())
    {
      if (!F->fPtr)
        RuntimeException::selfThrow("FirstOrderLinearR::computeF() is not linked to a plugin function");
      unsigned int sizeY = interaction->getSizeOfY();
      unsigned int sizeZ = interaction->getSizeZ();
      *workZ = *data[z];
      (F->fPtr)(time, sizeY, sizeZ, &(*F)(0, 0), sizeZ, &(*workZ)(0));
      // Copy data that might have been changed in the plug-in call.
      *data[z] = *workZ;
    }
    // else nothing
  }
}

void FirstOrderLinearR::computeE(const double time)
{
  if (e)
  {
    if (e->isPlugged())
    {
      if (!e->fPtr)
        RuntimeException::selfThrow("FirstOrderLinearR::computeE() is not linked to a plugin function");
      unsigned int sizeY = interaction->getSizeOfY();
      unsigned int sizeZ = interaction->getSizeZ();
      *workZ = *data[z];
      (e->fPtr)(time, sizeY, &(*e)(0), sizeZ, &(*workZ)(0));
      // Copy data that might have been changed in the plug-in call.
      *data[z] = *workZ;
    }
    // else nothing
  }
}

void FirstOrderLinearR::computeB(const double time)
{
  if (JacG[0]->isPlugged())
  {
    if (!JacG[0]->fPtr)
      RuntimeException::selfThrow("FirstOrderLinearR::computeB() is not linked to a plugin function");
    unsigned int sizeY = interaction->getSizeOfY();
    unsigned int sizeX = interaction->getSizeOfDS();
    unsigned int sizeZ = interaction->getSizeZ();
    *workZ = *data[z];
    JacG[0]->fPtr(time, sizeX, sizeY, &(*JacG[0])(0, 0), sizeZ, &(*workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *data[z] = *workZ;
  }
  // else nothing
}

void FirstOrderLinearR::computeH(double time)
{
  computeOutput(time, 0);
}

void FirstOrderLinearR::computeG(double time)
{
  computeInput(time, 0);
}

void FirstOrderLinearR::computeJacH(double time, unsigned int index)
{
  if (index == 0)
    computeC(time);
  else if (index == 1)
    computeD(time);
  else
    RuntimeException::selfThrow("FirstOrderLinearR::computeJacH(t,index), index is out of range.");
}

void FirstOrderLinearR::computeJacG(double time, unsigned int index)
{
  if (index == 0)
    computeB(time);
  else
    RuntimeException::selfThrow("FirstOrderLinearR::computeJacG(t,index), index is out of range.");
}

void FirstOrderLinearR::computeOutput(double time, unsigned int)
{
  computeC(time);
  computeD(time);
  computeF(time);
  computeE(time);

  // Note that the second argument remains unamed since it is not used: for first order systems, we always compute
  // y[0]

  // We get y and lambda of the interaction (pointers)
  SP::SiconosVector y = interaction->getYPtr(0);
  SP::SiconosVector lambda = interaction->getLambdaPtr(0);

  // compute y
  if (JacH[0])
    prod(*JacH[0], *data[x], *y);
  else
    y->zero();

  if (JacH[1])
    prod(*JacH[1], *lambda, *y, false);

  if (e)
    *y += *e;

  if (F)
    prod(*F, *data[z], *y, false);
}

void FirstOrderLinearR::computeInput(double time, unsigned int level)
{
  computeB(time);

  // We get lambda of the interaction (pointers)
  SP::SiconosVector lambda = interaction->getLambdaPtr(level);
  prod(*JacG[0], *lambda, *data[r], false);
}

void FirstOrderLinearR::display() const
{
  cout << " ===== Linear Time Invariant relation display ===== " << endl;
  cout << "| C " << endl;
  if (JacH[0]) JacH[0]->display();
  else cout << "->NULL" << endl;
  cout << "| D " << endl;
  if (JacH[1]) JacH[1]->display();
  else cout << "->NULL" << endl;
  cout << "| F " << endl;
  if (F) F->display();
  else cout << "->NULL" << endl;
  cout << "| e " << endl;
  if (e) e->display();
  else cout << "->NULL" << endl;
  cout << "| B " << endl;
  if (JacG[0]) JacG[0]->display();
  else cout << "->NULL" << endl;
  cout << " ================================================== " << endl;
}

void FirstOrderLinearR::saveRelationToXML() const
{
  //   if(!relationxml)
  RuntimeException::selfThrow("FirstOrderLinearR::saveRelationToXML, no yet implemented.");

  //   SP::LinearRXML folrXML = (boost::static_pointer_cast<LinearRXML>(relationxml));
  //   folrXML->setC( *C );
  //   folrXML->setD( *D );
  //   folrXML->setF( *F );
  //   folrXML->setE( *e );
  //   folrXML->setB( *B );
}

FirstOrderLinearR* FirstOrderLinearR::convert(Relation *r)
{
  return dynamic_cast<FirstOrderLinearR*>(r);
}

