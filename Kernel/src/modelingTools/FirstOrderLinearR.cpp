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
#include "FirstOrderLinearRXML.h"
#include "Interaction.h"

using namespace std;
using namespace RELATION;

void FirstOrderLinearR::initPluginFlags(bool in)
{
  isPlugged[RELATION::C] = in;
  isPlugged[RELATION::D] = in;
  isPlugged[RELATION::F] = in;
  isPlugged[RELATION::e] = in;
  isPlugged[RELATION::B] = in ;
}

// Default (private) constructor
FirstOrderLinearR::FirstOrderLinearR(): FirstOrderR(LinearR)
{
  initPluginFlags(false);
}


// xml constructor
FirstOrderLinearR::FirstOrderLinearR(SP::RelationXML relxml):
  FirstOrderR(relxml, LinearR)
{
  SP::FirstOrderLinearRXML folrXML = boost::static_pointer_cast<FirstOrderLinearRXML>(relationxml);
  // get matrices values. All are optional.

  initPluginFlags(false);

  if (folrXML->hasH())
    RuntimeException::selfThrow(" FirstOrderLinearR xml constructor failed. Too many inputs: you can not give h or its jacobian.");

  string plugin;
  if (folrXML->hasC())
  {
    if (folrXML->isCPlugin())
    {
      plugin = folrXML->getCPlugin();
      setComputeCFunction(SSL::getPluginName(plugin), SSL::getPluginFunctionName(plugin));
    }
    else
      C.reset(new SimpleMatrix(folrXML->getC()));
  }

  if (folrXML->hasD())
  {
    if (folrXML->isDPlugin())
    {
      plugin = folrXML->getDPlugin();
      setComputeDFunction(SSL::getPluginName(plugin), SSL::getPluginFunctionName(plugin));
    }
    else
      D.reset(new SimpleMatrix(folrXML->getD()));
  }

  if (folrXML->hasF())
  {
    if (folrXML->isFPlugin())
    {
      plugin = folrXML->getFPlugin();
      setComputeFFunction(SSL::getPluginName(plugin), SSL::getPluginFunctionName(plugin));
    }
    else
      F.reset(new SimpleMatrix(folrXML->getF()));
  }

  if (folrXML->hasE())
  {
    if (folrXML->isEPlugin())
    {
      plugin = folrXML->getEPlugin();
      setComputeEFunction(SSL::getPluginName(plugin), SSL::getPluginFunctionName(plugin));
    }
    else
      e.reset(new SimpleVector(folrXML->getE()));
  }

  if (folrXML->hasB())
  {
    if (folrXML->hasG())
      RuntimeException::selfThrow(" FirstOrderLinearR xml constructor failed. Too many inputs: you can not give B and g or its jacobian.");

    if (folrXML->isBPlugin())
    {
      plugin = folrXML->getBPlugin();
      setComputeBFunction(SSL::getPluginName(plugin), SSL::getPluginFunctionName(plugin));
    }
    else
      B.reset(new SimpleMatrix(folrXML->getB()));
  }
}

// Constructor with C and B plug-in names
FirstOrderLinearR::FirstOrderLinearR(const string& CName, const string& BName):
  FirstOrderR(LinearR)
{
  initPluginFlags(false);

  setComputeCFunction(SSL::getPluginName(CName), SSL::getPluginFunctionName(CName));
  setComputeBFunction(SSL::getPluginName(BName), SSL::getPluginFunctionName(BName));
}

// Constructor with e plug-in name
FirstOrderLinearR::FirstOrderLinearR(const string& EName):
  FirstOrderR(LinearR)
{
  initPluginFlags(false);
  setComputeEFunction(SSL::getPluginName(EName), SSL::getPluginFunctionName(EName));
}

// Constructor from a complete set of data (plugin)
FirstOrderLinearR::FirstOrderLinearR(const string& CName, const string& DName, const string& FName, const string& EName, const string& BName):
  FirstOrderR(LinearR)
{
  initPluginFlags(true);

  setComputeCFunction(SSL::getPluginName(CName), SSL::getPluginFunctionName(CName));
  setComputeDFunction(SSL::getPluginName(DName), SSL::getPluginFunctionName(DName));
  setComputeFFunction(SSL::getPluginName(FName), SSL::getPluginFunctionName(FName));
  setComputeEFunction(SSL::getPluginName(EName), SSL::getPluginFunctionName(EName));
  setComputeBFunction(SSL::getPluginName(BName), SSL::getPluginFunctionName(BName));
}

// Minimum data (C, B as pointers) constructor
FirstOrderLinearR::FirstOrderLinearR(SP::SiconosMatrix newC, SP::SiconosMatrix newB):
  FirstOrderR(LinearR)
{
  initPluginFlags(false);

  C = newC;

}

// Constructor from a complete set of data
FirstOrderLinearR::FirstOrderLinearR(SP::SiconosMatrix newC, SP::SiconosMatrix newD, SP::SiconosMatrix newF, SP::SiconosVector newE, SP::SiconosMatrix newB):
  FirstOrderR(LinearR)
{
  initPluginFlags(false);

  C = newC;
  D = newD;
  F = newF;
  e = newE;
  B = newB;
}

FirstOrderLinearR::~FirstOrderLinearR()
{}

void FirstOrderLinearR::initialize()
{
  // Note: do not call FirstOrderR::initialize to avoid jacobianH and jacobianG allocation.

  // Check if an Interaction is connected to the Relation.
  if (!interaction)
    RuntimeException::selfThrow("FirstOrderR::initialize failed. No Interaction linked to the present relation.");

  // Update data member (links to DS variables)
  initDSLinks();

  // Check if various operators sizes are consistent.
  // Reference: interaction.
  unsigned int sizeY = interaction->getSizeOfY();
  unsigned int sizeX = interaction->getSizeOfDS();
  unsigned int sizeZ = interaction->getSizeZ();

  if (C)
  {
    if (C->size(0) != sizeY || C->size(1) != sizeX)
      RuntimeException::selfThrow("FirstOrderLinearR::initialize , inconsistent size between C and Interaction.");
  }
  if (B)
  {
    if (B->size(0) != sizeX || B->size(1) != sizeY)
      RuntimeException::selfThrow("FirstOrderLinearR::initialize , inconsistent size between C and B.");
  }

  if (D)
  {
    if (D->size(0) != sizeY || D->size(1) != sizeY)
      RuntimeException::selfThrow("FirstOrderLinearR::initialize , inconsistent size between C and D.");
  }

  if (F && (F->size(0) != sizeY))
    RuntimeException::selfThrow("FirstOrderLinearR::initialize , inconsistent size between C and F.");

  if (e  && e->size() != sizeY)
    RuntimeException::selfThrow("FirstOrderLinearR::initialize , inconsistent size between C and e.");

  // Memory allocation if required (ie if plugged and not allocated: must be done here since in constructors, interaction is not knonw).
  if (! C && isPlugged[RELATION::C])
  {
    C.reset(new SimpleMatrix(sizeY, sizeX));
  }

  if (! D && isPlugged[RELATION::D])
  {
    D.reset(new SimpleMatrix(sizeY, sizeY));
  }

  if (! F && isPlugged[RELATION::F])
  {
    F.reset(new SimpleMatrix(sizeY, sizeZ));
  }

  if (! e && isPlugged[RELATION::e])
  {
    e.reset(new SimpleVector(sizeY));
  }

  if (! B && isPlugged[RELATION::B])
  {
    B.reset(new SimpleMatrix(sizeX, sizeY));
  }
  workZ.reset(new SimpleVector(sizeZ));
}

// setters

void FirstOrderLinearR::setC(const SiconosMatrix& newValue)
{
  if (! C)
    C.reset(new SimpleMatrix(newValue));
  else
  {
    if ((newValue.size(0) == C->size(0)) && (newValue.size(1) == C->size(1)))
      *C = newValue;
    else
      RuntimeException::selfThrow("FirstOrderLinearR - setC: inconsistent dimensions with problem size for input matrix C.");
  }

  isPlugged[RELATION::C] = false;

}

void FirstOrderLinearR::setCPtr(SP::SiconosMatrix newPtr)
{
  C = newPtr;
  isPlugged[RELATION::C] = false;
}

void FirstOrderLinearR::setComputeCFunction(const string& pluginPath, const string& functionName)
{
  isPlugged[RELATION::C] = Plugin::setFunction(&CPtr, pluginPath, functionName, pluginNames[RELATION::C]);
}

void FirstOrderLinearR::setD(const SiconosMatrix& newValue)
{
  if (! D)
    D.reset(new SimpleMatrix(newValue));
  else
  {
    if ((newValue.size(0) == D->size(0)) && (newValue.size(1) == D->size(1)))
      *D = newValue;
    else
      RuntimeException::selfThrow("FirstOrderLinearR - setD: inconsistent dimensions with problem size for input matrix D.");
  }
  isPlugged[RELATION::D] = false;
}

void FirstOrderLinearR::setDPtr(SP::SiconosMatrix newPtr)
{
  D = newPtr;
  isPlugged[RELATION::D] = false;
}

void FirstOrderLinearR::setComputeDFunction(const string& pluginPath, const string& functionName)
{
  isPlugged[RELATION::D] = Plugin::setFunction(&DPtr, pluginPath, functionName, pluginNames[RELATION::D]);
}

void FirstOrderLinearR::setF(const SiconosMatrix& newValue)
{
  if (! F)
    F.reset(new SimpleMatrix(newValue));
  else
  {
    if ((newValue.size(0) == F->size(0)) && (newValue.size(1) == F->size(1)))
      *F = newValue;
    else
      RuntimeException::selfThrow("FirstOrderLinearR - setF: inconsistent dimensions with problem size for input matrix F.");
  }
  isPlugged[RELATION::F] = false;
}

void FirstOrderLinearR::setFPtr(SP::SiconosMatrix newPtr)
{
  F = newPtr;
  isPlugged[RELATION::F] = false;
}

void FirstOrderLinearR::setComputeFFunction(const string& pluginPath, const string& functionName)
{
  isPlugged[RELATION::F] = Plugin::setFunction(&FPtr, pluginPath, functionName, pluginNames[RELATION::F]);
}

void FirstOrderLinearR::setE(const SiconosVector& newValue)
{
  if (! e)
    e.reset(new SimpleVector(newValue));
  else
  {
    if (newValue.size() == e->size())
      *e = newValue;
    else
      RuntimeException::selfThrow("FirstOrderLinearR - setE: inconsistent dimensions with problem size for input vector e.");
  }
  isPlugged[RELATION::e] = false;
}

void FirstOrderLinearR::setEPtr(SP::SiconosVector newPtr)
{
  e = newPtr;
  isPlugged[RELATION::e] = false;
}
void FirstOrderLinearR::setComputeEFunction(FOVecPtr ptrFunct)
{
  ePtr = ptrFunct;
  isPlugged[RELATION::e] = true;
}

void FirstOrderLinearR::setComputeEFunction(const string& pluginPath, const string& functionName)
{
  isPlugged[RELATION::e] = Plugin::setFunction(&ePtr, pluginPath, functionName, pluginNames[RELATION::e]);
}

void FirstOrderLinearR::setB(const SiconosMatrix& newValue)
{
  if (! B)
    B.reset(new SimpleMatrix(newValue));
  else
  {
    if ((newValue.size(0) == B->size(0)) && (newValue.size(1) == B->size(1)))
      *B = newValue;
    else
      RuntimeException::selfThrow("FirstOrderLinearR - setB: inconsistent dimensions with problem size for input matrix B.");
  }
  isPlugged[RELATION::B] = false;
}

void FirstOrderLinearR::setBPtr(SP::SiconosMatrix newPtr)
{
  B = newPtr;
  isPlugged[RELATION::B] = false;
}

void FirstOrderLinearR::setComputeBFunction(const string& pluginPath, const string& functionName)
{
  isPlugged[RELATION::B] = Plugin::setFunction(&BPtr, pluginPath, functionName, pluginNames[RELATION::B]);
}

void FirstOrderLinearR::computeOutput(double time, unsigned int)
{
  if (isPlugged[RELATION::C])
    computeC(time);
  if (isPlugged[RELATION::D])
    computeD(time);
  if (isPlugged[RELATION::F])
    computeF(time);
  if (isPlugged[RELATION::e])
    computeE(time);

  // Note that the second argument remains unamed since it is not used: for first order systems, we always compute
  // y[0]

  // We get y and lambda of the interaction (pointers)
  SP::SiconosVector y = interaction->getYPtr(0);
  SP::SiconosVector lambda = interaction->getLambdaPtr(0);

  // compute y
  if (C)
    prod(*C, *data["x"], *y);
  else
    y->zero();

  if (D)
    prod(*D, *lambda, *y, false);

  if (e)
    *y += *e;

  if (F)
    prod(*F, *data["z"], *y, false);
}

void FirstOrderLinearR::computeInput(double time, unsigned int level)
{
  if (isPlugged[RELATION::B])
    computeB(time);

  // We get lambda of the interaction (pointers)
  SP::SiconosVector lambda = interaction->getLambdaPtr(level);
  prod(*B, *lambda, *data["r"], false);
}

void FirstOrderLinearR::computeC(const double time)
{
  if (isPlugged[RELATION::C])
  {
    if (!CPtr)
      RuntimeException::selfThrow("computeC() is not linked to a plugin function");
    unsigned int sizeY = interaction->getSizeOfY();
    unsigned int sizeX = interaction->getSizeOfDS();
    unsigned int sizeZ = interaction->getSizeZ();
    *workZ = *data["z"];
    CPtr(time, sizeY, sizeX, &(*C)(0, 0), sizeZ, &(*workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *data["z"] = *workZ;
  }
  // else nothing
}

void FirstOrderLinearR::computeD(const double time)
{
  if (isPlugged[RELATION::D])
  {
    if (!DPtr)
      RuntimeException::selfThrow("computeD() is not linked to a plugin function");
    unsigned int sizeY = interaction->getSizeOfY();
    unsigned int sizeZ = interaction->getSizeZ();
    *workZ = *data["z"];
    DPtr(time, sizeY, &(*D)(0, 0), sizeZ, &(*workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *data["z"] = *workZ;
  }
  // else nothing
}

void FirstOrderLinearR::computeF(const double time)
{
  if (isPlugged[RELATION::F])
  {
    if (!FPtr)
      RuntimeException::selfThrow("computeF() is not linked to a plugin function");
    unsigned int sizeY = interaction->getSizeOfY();
    unsigned int sizeZ = interaction->getSizeZ();
    *workZ = *data["z"];
    FPtr(time, sizeY, &(*F)(0, 0), sizeZ, &(*workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *data["z"] = *workZ;
  }
  // else nothing
}

void FirstOrderLinearR::computeE(const double time)
{
  if (isPlugged[RELATION::e])
  {
    if (!ePtr)
      RuntimeException::selfThrow("computeE() is not linked to a plugin function");
    unsigned int sizeY = interaction->getSizeOfY();
    unsigned int sizeZ = interaction->getSizeZ();
    *workZ = *data["z"];
    ePtr(time, sizeY, &(*e)(0), sizeZ, &(*workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *data["z"] = *workZ;
  }
  // else nothing
}

void FirstOrderLinearR::computeB(const double time)
{
  if (isPlugged[RELATION::B])
  {
    if (!BPtr)
      RuntimeException::selfThrow("computeB() is not linked to a plugin function");
    unsigned int sizeY = interaction->getSizeOfY();
    unsigned int sizeX = interaction->getSizeOfDS();
    unsigned int sizeZ = interaction->getSizeZ();
    *workZ = *data["z"];
    BPtr(time, sizeX, sizeY, &(*B)(0, 0), sizeZ, &(*workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *data["z"] = *workZ;
  }
  // else nothing
}

void FirstOrderLinearR::display() const
{
  cout << " ===== Linear Time Invariant relation display ===== " << endl;
  cout << "| C " << endl;
  if (C) C->display();
  else cout << "->NULL" << endl;
  cout << "| D " << endl;
  if (D) D->display();
  else cout << "->NULL" << endl;
  cout << "| F " << endl;
  if (F) F->display();
  else cout << "->NULL" << endl;
  cout << "| e " << endl;
  if (e) e->display();
  else cout << "->NULL" << endl;
  cout << "| B " << endl;
  if (B) B->display();
  else cout << "->NULL" << endl;
  cout << " ================================================== " << endl;
}

void FirstOrderLinearR::saveRelationToXML() const
{
  if (!relationxml)
    RuntimeException::selfThrow("FirstOrderLinearR::saveRelationToXML, no xml object found.");

  SP::FirstOrderLinearRXML folrXML = (boost::static_pointer_cast<FirstOrderLinearRXML>(relationxml));
  folrXML->setC(*C);
  folrXML->setD(*D);
  folrXML->setF(*F);
  folrXML->setE(*e);
  folrXML->setB(*B);
}

FirstOrderLinearR* FirstOrderLinearR::convert(Relation *r)
{
  return dynamic_cast<FirstOrderLinearR*>(r);
}

