/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
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

void FirstOrderLinearR::initAllocationFlags(bool in)
{
  isAllocatedIn["C"] = in;
  isAllocatedIn["D"] = in;
  isAllocatedIn["F"] = in;
  isAllocatedIn["e"] = in;
  isAllocatedIn["B"] = in;
}

void FirstOrderLinearR::initPluginFlags(bool in)
{
  isPlugged["C"] = in;
  isPlugged["D"] = in;
  isPlugged["F"] = in;
  isPlugged["e"] = in;
  isPlugged["B"] = in ;
}

// Default (private) constructor
FirstOrderLinearR::FirstOrderLinearR(): FirstOrderR("LinearR"), C(NULL), D(NULL), F(NULL), e(NULL), B(NULL),
  CPtr(NULL), DPtr(NULL), FPtr(NULL), ePtr(NULL), BPtr(NULL)
{
  initPluginFlags(false);
  initAllocationFlags(false);
}

// xml constructor
FirstOrderLinearR::FirstOrderLinearR(RelationXML* relxml):
  FirstOrderR(relxml, "LinearR"), C(NULL), D(NULL), F(NULL), e(NULL), B(NULL),
  CPtr(NULL), DPtr(NULL), FPtr(NULL), ePtr(NULL), BPtr(NULL)
{
  FirstOrderLinearRXML * folrXML = static_cast<FirstOrderLinearRXML *>(relationxml);
  // get matrices values. All are optional.

  initAllocationFlags(false);
  initPluginFlags(false);

  string plugin;
  if (folrXML->hasC())
  {
    if (folrXML->hasH())
      RuntimeException::selfThrow(" FirstOrderLinearR xml constructor failed. Too many inputs: you can not give C and h or its jacobian.");

    if (folrXML->isCPlugin())
    {
      plugin = folrXML->getCPlugin();
      setComputeCFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else
    {
      C = new SimpleMatrix(folrXML->getC());
      isAllocatedIn["C"] = true;
    }
  }

  if (folrXML->hasD())
  {
    if (folrXML->hasH())
      RuntimeException::selfThrow(" FirstOrderLinearR xml constructor failed. Too many inputs: you can not give D and h or its jacobian.");

    if (folrXML->isDPlugin())
    {
      plugin = folrXML->getDPlugin();
      setComputeDFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else
    {
      D = new SimpleMatrix(folrXML->getD());
      isAllocatedIn["D"] = true;
    }
  }

  if (folrXML->hasF())
  {
    if (folrXML->isFPlugin())
    {
      plugin = folrXML->getFPlugin();
      setComputeFFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else
    {
      F = new SimpleMatrix(folrXML->getF());
      isAllocatedIn["F"] = true;
    }
  }

  if (folrXML->hasE())
  {
    if (folrXML->isEPlugin())
    {
      plugin = folrXML->getEPlugin();
      setComputeEFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else
    {
      e = new SimpleVector(folrXML->getE());
      isAllocatedIn["e"] = true;
    }
  }

  if (folrXML->hasB())
  {
    if (folrXML->hasG())
      RuntimeException::selfThrow(" FirstOrderLinearR xml constructor failed. Too many inputs: you can not give B and g or its jacobian.");

    if (folrXML->isBPlugin())
    {
      plugin = folrXML->getBPlugin();
      setComputeBFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else
    {
      B = new SimpleMatrix(folrXML->getB());
      isAllocatedIn["B"] = true;
    }
  }
}

// Constructor with C and B plug-in names
FirstOrderLinearR::FirstOrderLinearR(const string& CName, const string& BName):
  FirstOrderR("LinearR"), C(NULL), D(NULL), F(NULL), e(NULL), B(NULL),
  CPtr(NULL), DPtr(NULL), FPtr(NULL), ePtr(NULL), BPtr(NULL)
{
  initPluginFlags(false);
  initAllocationFlags(false);
  setComputeCFunction(cShared.getPluginName(CName), cShared.getPluginFunctionName(CName));
  setComputeBFunction(cShared.getPluginName(BName), cShared.getPluginFunctionName(BName));
}

// Constructor with e plug-in name
FirstOrderLinearR::FirstOrderLinearR(const string& EName):
  FirstOrderR("LinearR"), C(NULL), D(NULL), F(NULL), e(NULL), B(NULL),
  CPtr(NULL), DPtr(NULL), FPtr(NULL), ePtr(NULL), BPtr(NULL)
{
  initPluginFlags(false);
  initAllocationFlags(false);
  setComputeEFunction(cShared.getPluginName(EName), cShared.getPluginFunctionName(EName));
}

// Constructor from a complete set of data (plugin)
FirstOrderLinearR::FirstOrderLinearR(const string& CName, const string& DName, const string& FName, const string& EName, const string& BName):
  FirstOrderR("LinearR"), C(NULL), D(NULL), F(NULL), e(NULL), B(NULL),
  CPtr(NULL), DPtr(NULL), FPtr(NULL), ePtr(NULL), BPtr(NULL)
{
  initPluginFlags(true);
  initAllocationFlags(false);
  setComputeCFunction(cShared.getPluginName(CName), cShared.getPluginFunctionName(CName));
  setComputeDFunction(cShared.getPluginName(DName), cShared.getPluginFunctionName(DName));
  setComputeFFunction(cShared.getPluginName(FName), cShared.getPluginFunctionName(FName));
  setComputeEFunction(cShared.getPluginName(EName), cShared.getPluginFunctionName(EName));
  setComputeBFunction(cShared.getPluginName(BName), cShared.getPluginFunctionName(BName));
}

// Minimum data (C, B as pointers) constructor
FirstOrderLinearR::FirstOrderLinearR(SiconosMatrix * newC, SiconosMatrix * newB):
  FirstOrderR("LinearR"), C(newC), D(NULL), F(NULL), e(NULL), B(newB),
  CPtr(NULL), DPtr(NULL), FPtr(NULL), ePtr(NULL), BPtr(NULL)
{
  initPluginFlags(false);
  initAllocationFlags(false);
}

// Constructor from a complete set of data
FirstOrderLinearR::FirstOrderLinearR(SiconosMatrix* newC, SiconosMatrix* newD, SiconosMatrix* newF, SiconosVector* newE, SiconosMatrix* newB):
  FirstOrderR("LinearR"), C(newC), D(newD), F(newF), e(newE), B(newB),
  CPtr(NULL), DPtr(NULL), FPtr(NULL), ePtr(NULL), BPtr(NULL)
{
  initPluginFlags(false);
  initAllocationFlags(false);
}

FirstOrderLinearR::~FirstOrderLinearR()
{
  if (isAllocatedIn["C"]) delete C;
  C = NULL;
  if (isAllocatedIn["D"]) delete D;
  D = NULL;
  if (isAllocatedIn["F"]) delete F;
  F = NULL;
  if (isAllocatedIn["e"]) delete e;
  e = NULL;
  if (isAllocatedIn["B"]) delete B;
  B = NULL;
  CPtr = NULL;
  DPtr = NULL;
  FPtr = NULL;
  ePtr = NULL;
  BPtr = NULL;
}

void FirstOrderLinearR::initialize()
{
  // Note: do not call FirstOrderR::initialize to avoid jacobianH and jacobianG allocation.

  // Check if an Interaction is connected to the Relation.
  if (interaction == NULL)
    RuntimeException::selfThrow("FirstOrderR::initialize failed. No Interaction linked to the present relation.");

  // Update data member (links to DS variables)
  initDSLinks();

  // Check if various operators sizes are consistent.
  // Reference: interaction.
  unsigned int sizeY = interaction->getSizeOfY();
  unsigned int sizeX = interaction->getSizeOfDS();
  unsigned int sizeZ = interaction->getSizeZ();

  if (C != NULL)
  {
    if (C->size(0) != sizeY || C->size(1) != sizeX)
      RuntimeException::selfThrow("FirstOrderLinearR::initialize , inconsistent size between C and Interaction.");
    if (jacobianH[0] != NULL)
      RuntimeException::selfThrow("FirstOrderLinearR::initialize failed. Either C and jacobianH are defined: possible conflict.");
    jacobianH[0] = C;
  }
  if (B != NULL)
  {
    if (B->size(0) != sizeX || B->size(1) != sizeY)
      RuntimeException::selfThrow("FirstOrderLinearR::initialize , inconsistent size between C and B.");
    if (jacobianG[0] != NULL)
      RuntimeException::selfThrow("FirstOrderLinearR::initialize failed. Either B and jacobianG are defined: possible conflict.");
    jacobianG[0] = B;
  }

  if (D != NULL)
  {
    if (D->size(0) != sizeY || D->size(1) != sizeY)
      RuntimeException::selfThrow("FirstOrderLinearR::initialize , inconsistent size between C and D.");
    if (jacobianH[1] != NULL)
      RuntimeException::selfThrow("FirstOrderLinearR::initialize failed. Either D and jacobianH are defined: possible conflict.");
    jacobianH[1] = D;
  }

  if (F != NULL && (F->size(0) != sizeY))
    RuntimeException::selfThrow("FirstOrderLinearR::initialize , inconsistent size between C and F.");

  if (e != NULL && e->size() != sizeY)
    RuntimeException::selfThrow("FirstOrderLinearR::initialize , inconsistent size between C and e.");

  // Memory allocation if required (ie if plugged and not allocated: must be done here since in constructors, interaction is not knonw).
  if (C == NULL && isPlugged["C"])
  {
    C = new SimpleMatrix(sizeY, sizeX);
    isAllocatedIn["C"] = true ;
  }

  if (D == NULL && isPlugged["D"])
  {
    D = new SimpleMatrix(sizeY, sizeY);
    isAllocatedIn["D"] = true ;
  }

  if (F == NULL && isPlugged["F"])
  {
    F = new SimpleMatrix(sizeY, sizeZ);
    isAllocatedIn["F"] = true ;
  }

  if (e == NULL && isPlugged["e"])
  {
    e = new SimpleVector(sizeY);
    isAllocatedIn["e"] = true ;
  }

  if (B == NULL && isPlugged["B"])
  {
    B = new SimpleMatrix(sizeX, sizeY);
    isAllocatedIn["B"] = true ;
  }
}

// setters

void FirstOrderLinearR::setC(const SiconosMatrix& newValue)
{
  if (C == NULL)
  {
    C = new SimpleMatrix(newValue);
    isAllocatedIn["C"] = true;
  }
  else
  {
    if ((newValue.size(0) == C->size(0)) && (newValue.size(1) == C->size(1)))
      *C = newValue;
    else
      RuntimeException::selfThrow("FirstOrderLinearR - setC: inconsistent dimensions with problem size for input matrix C.");
  }

  isPlugged["C"] = false;

}

void FirstOrderLinearR::setCPtr(SiconosMatrix *newPtr)
{
  if (isAllocatedIn["C"]) delete C;
  C = newPtr;
  isAllocatedIn["C"] = false;
  isPlugged["C"] = false;
}

void FirstOrderLinearR::setComputeCFunction(const string& pluginPath, const string& functionName)
{
  CPtr = NULL;
  cShared.setFunction(&CPtr, pluginPath, functionName);

  string plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames["C"] = plugin + ":" + functionName;
  isPlugged["C"] = true;
}

void FirstOrderLinearR::setD(const SiconosMatrix& newValue)
{
  if (D == NULL)
  {
    D = new SimpleMatrix(newValue);
    isAllocatedIn["D"] = true;
  }
  else
  {
    if ((newValue.size(0) == D->size(0)) && (newValue.size(1) == D->size(1)))
      *D = newValue;
    else
      RuntimeException::selfThrow("FirstOrderLinearR - setD: inconsistent dimensions with problem size for input matrix D.");
  }
  isPlugged["D"] = false;
}

void FirstOrderLinearR::setDPtr(SiconosMatrix *newPtr)
{
  if (isAllocatedIn["D"]) delete D;
  D = newPtr;
  isAllocatedIn["D"] = false;
  isPlugged["D"] = false;
}

void FirstOrderLinearR::setComputeDFunction(const string& pluginPath, const string& functionName)
{
  DPtr = NULL;
  cShared.setFunction(&DPtr, pluginPath, functionName);

  string plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames["D"] = plugin + ":" + functionName;
  isPlugged["D"] = true;
}

void FirstOrderLinearR::setF(const SiconosMatrix& newValue)
{
  if (F == NULL)
  {
    F = new SimpleMatrix(newValue);
    isAllocatedIn["F"] = true;
  }
  else
  {
    if ((newValue.size(0) == F->size(0)) && (newValue.size(1) == F->size(1)))
      *F = newValue;
    else
      RuntimeException::selfThrow("FirstOrderLinearR - setF: inconsistent dimensions with problem size for input matrix F.");
  }
  isPlugged["F"] = false;
}

void FirstOrderLinearR::setFPtr(SiconosMatrix *newPtr)
{
  if (isAllocatedIn["F"]) delete F;
  F = newPtr;
  isAllocatedIn["F"] = false;
  isPlugged["F"] = false;
}

void FirstOrderLinearR::setComputeFFunction(const string& pluginPath, const string& functionName)
{
  FPtr = NULL;
  cShared.setFunction(&FPtr, pluginPath, functionName);

  string plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames["F"] = plugin + ":" + functionName;
  isPlugged["F"] = true;
}

void FirstOrderLinearR::setE(const SiconosVector& newValue)
{
  if (e == NULL)
  {
    e = new SimpleVector(newValue);
    isAllocatedIn["e"] = true;
  }
  else
  {
    if (newValue.size() == e->size())
      *e = newValue;
    else
      RuntimeException::selfThrow("FirstOrderLinearR - setE: inconsistent dimensions with problem size for input vector e.");
  }
  isPlugged["e"] = false;
}

void FirstOrderLinearR::setEPtr(SiconosVector* newPtr)
{
  if (isAllocatedIn["e"]) delete e;
  e = newPtr;
  isAllocatedIn["e"] = false;
  isPlugged["e"] = false;
}

void FirstOrderLinearR::setComputeEFunction(const string& pluginPath, const string& functionName)
{
  ePtr = NULL;
  cShared.setFunction(&ePtr, pluginPath, functionName);

  string plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames["e"] = plugin + ":" + functionName;
  isPlugged["e"] = true;
}

void FirstOrderLinearR::setB(const SiconosMatrix& newValue)
{
  if (B == NULL)
  {
    B = new SimpleMatrix(newValue);
    isAllocatedIn["B"] = true;
  }
  else
  {
    if ((newValue.size(0) == B->size(0)) && (newValue.size(1) == B->size(1)))
      *B = newValue;
    else
      RuntimeException::selfThrow("FirstOrderLinearR - setB: inconsistent dimensions with problem size for input matrix B.");
  }
  isPlugged["B"] = false;
}

void FirstOrderLinearR::setBPtr(SiconosMatrix *newPtr)
{
  if (isAllocatedIn["B"]) delete B;
  B = newPtr;
  isAllocatedIn["B"] = false;
  isPlugged["B"] = false;
}

void FirstOrderLinearR::setComputeBFunction(const string& pluginPath, const string& functionName)
{

  BPtr = NULL;
  cShared.setFunction(&BPtr, pluginPath, functionName);

  string plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames["B"] = plugin + ":" + functionName;
  isPlugged["B"] = true;
}

void FirstOrderLinearR::computeOutput(double time, unsigned int)
{
  if (isPlugged["C"])
    computeC(time);
  if (isPlugged["D"])
    computeD(time);
  if (isPlugged["F"])
    computeF(time);
  if (isPlugged["e"])
    computeE(time);

  // Note that the second argument remains unamed since it is not used: for first order systems, we always compute
  // y[0]

  // We get y and lambda of the interaction (pointers)
  SiconosVector *y = interaction->getYPtr(0);
  SiconosVector *lambda = interaction->getLambdaPtr(0);

  // compute y
  if (C != NULL)
    prod(*C, *data["x"], *y);
  else
    y->zero();

  if (D != NULL)
    prod(*D, *lambda, *y, false);

  if (e != NULL)
    *y += *e;

  if (F != NULL)
    prod(*F, *data["z"], *y, false);
}

void FirstOrderLinearR::computeFreeOutput(double time, unsigned int)
{
  if (isPlugged["C"])
    computeC(time);
  if (isPlugged["D"])
    computeD(time);
  if (isPlugged["F"])
    computeF(time);
  if (isPlugged["e"])
    computeE(time);
  // Note that the second argument remains unamed since it is not used: for first order systems, we always compute
  // y[0]

  SiconosVector *yFree = interaction->getYPtr(0);
  // warning : yFree is saved in y !!

  // compute yFree
  if (C != NULL)
    prod(*C, *data["xFree"], *yFree);
  else
    yFree->zero();

  if (e != NULL)
    *yFree += *e;

  if (F != NULL)
    prod(*F, *data["z"], *yFree, false);
}

void FirstOrderLinearR::computeInput(double time, unsigned int level)
{
  if (isPlugged["B"])
    computeB(time);

  // We get lambda of the interaction (pointers)
  SiconosVector *lambda = interaction->getLambdaPtr(level);
  prod(*B, *lambda, *data["r"], false);
}

void FirstOrderLinearR::computeC(const double time)
{
  if (isPlugged["C"])
  {
    if (CPtr == NULL)
      RuntimeException::selfThrow("computeC() is not linked to a plugin function");
    unsigned int sizeY = interaction->getSizeOfY();
    unsigned int sizeX = interaction->getSizeOfDS();
    unsigned int sizeZ = interaction->getSizeZ();
    SimpleVector * zCopy = new SimpleVector(*data["z"]);
    CPtr(time, sizeY, sizeX, &(*C)(0, 0), sizeZ, &(*zCopy)(0));
    // Copy data that might have been changed in the plug-in call.
    *data["z"] = *zCopy;
    delete zCopy;
  }
  // else nothing
}

void FirstOrderLinearR::computeD(const double time)
{
  if (isPlugged["D"])
  {
    if (DPtr == NULL)
      RuntimeException::selfThrow("computeD() is not linked to a plugin function");
    unsigned int sizeY = interaction->getSizeOfY();
    unsigned int sizeZ = interaction->getSizeZ();
    SimpleVector * zCopy = new SimpleVector(*data["z"]);
    DPtr(time, sizeY, &(*D)(0, 0), sizeZ, &(*zCopy)(0));
    // Copy data that might have been changed in the plug-in call.
    *data["z"] = *zCopy;
    delete zCopy;
  }
  // else nothing
}

void FirstOrderLinearR::computeF(const double time)
{
  if (isPlugged["F"])
  {
    if (FPtr == NULL)
      RuntimeException::selfThrow("computeF() is not linked to a plugin function");
    unsigned int sizeY = interaction->getSizeOfY();
    unsigned int sizeZ = interaction->getSizeZ();
    SimpleVector * zCopy = new SimpleVector(*data["z"]);
    FPtr(time, sizeY, &(*F)(0, 0), sizeZ, &(*zCopy)(0));
    // Copy data that might have been changed in the plug-in call.
    *data["z"] = *zCopy;
    delete zCopy;
  }
  // else nothing
}

void FirstOrderLinearR::computeE(const double time)
{
  if (isPlugged["e"])
  {
    if (ePtr == NULL)
      RuntimeException::selfThrow("computeE() is not linked to a plugin function");
    unsigned int sizeY = interaction->getSizeOfY();
    unsigned int sizeZ = interaction->getSizeZ();
    SimpleVector * zCopy = new SimpleVector(*data["z"]);
    ePtr(time, sizeY, &(*e)(0), sizeZ, &(*zCopy)(0));
    // Copy data that might have been changed in the plug-in call.
    *data["z"] = *zCopy;
    delete zCopy;
  }
  // else nothing
}

void FirstOrderLinearR::computeB(const double time)
{
  if (isPlugged["B"])
  {
    if (BPtr == NULL)
      RuntimeException::selfThrow("computeB() is not linked to a plugin function");
    unsigned int sizeY = interaction->getSizeOfY();
    unsigned int sizeX = interaction->getSizeOfDS();
    unsigned int sizeZ = interaction->getSizeZ();
    SimpleVector * zCopy = new SimpleVector(*data["z"]);
    BPtr(time, sizeX, sizeY, &(*B)(0, 0), sizeZ, &(*zCopy)(0));
    // Copy data that might have been changed in the plug-in call.
    *data["z"] = *zCopy;
    delete zCopy;
  }
  // else nothing
}

void FirstOrderLinearR::computeJacobianH(double time, unsigned int i)
{
  if (i == 0)
    computeC(time);
  else if (i == 1)
    computeD(time);
  else
    RuntimeException::selfThrow("FirstOrderLinearR::computeJacobianH() failed; index out of range.");
}

void FirstOrderLinearR::computeJacobianG(double time, unsigned int i)
{
  computeB(time);
}

void FirstOrderLinearR::display() const
{
  cout << " ===== Linear Time Invariant relation display ===== " << endl;
  cout << "| C " << endl;
  if (C != NULL) C->display();
  else cout << "->NULL" << endl;
  cout << "| D " << endl;
  if (D != NULL) D->display();
  else cout << "->NULL" << endl;
  cout << "| F " << endl;
  if (F != NULL) F->display();
  else cout << "->NULL" << endl;
  cout << "| e " << endl;
  if (e != NULL) e->display();
  else cout << "->NULL" << endl;
  cout << "| B " << endl;
  if (B != NULL) B->display();
  else cout << "->NULL" << endl;
  cout << " ================================================== " << endl;
}

void FirstOrderLinearR::saveRelationToXML() const
{
  if (relationxml == NULL)
    RuntimeException::selfThrow("FirstOrderLinearR::saveRelationToXML, no xml object found.");

  FirstOrderLinearRXML * folrXML = (static_cast<FirstOrderLinearRXML*>(relationxml));
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

