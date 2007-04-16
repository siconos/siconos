/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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
  FirstOrderLinearRXML * lTIRxml = static_cast<FirstOrderLinearRXML *>(relationxml);
  // get matrices values. All are optional.

  initAllocationFlags(false);
  initPluginFlags(false);

  string plugin;
  if (lTIRxml->hasC())
  {
    if (lTIRxml->hasH())
      RuntimeException::selfThrow(" FirstOrderLinearR xml constructor failed. Too many inputs: you can not give C and h or its jacobian.");

    if (lTIRxml->isCPlugin())
    {
      plugin = lTIRxml->getCPlugin();
      setComputeCFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else
    {
      C = new SimpleMatrix(lTIRxml->getC());
      isAllocatedIn["C"] = true;
    }
  }

  if (lTIRxml->hasD())
  {
    if (lTIRxml->hasH())
      RuntimeException::selfThrow(" FirstOrderLinearR xml constructor failed. Too many inputs: you can not give D and h or its jacobian.");

    if (lTIRxml->isDPlugin())
    {
      plugin = lTIRxml->getDPlugin();
      setComputeDFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else
    {
      D = new SimpleMatrix(lTIRxml->getD());
      isAllocatedIn["D"] = true;
    }
  }

  if (lTIRxml->hasF())
  {
    if (lTIRxml->isFPlugin())
    {
      plugin = lTIRxml->getFPlugin();
      setComputeFFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else
    {
      F = new SimpleMatrix(lTIRxml->getF());
      isAllocatedIn["F"] = true;
    }
  }

  if (lTIRxml->hasE())
  {
    if (lTIRxml->isEPlugin())
    {
      plugin = lTIRxml->getEPlugin();
      setComputeEFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else
    {
      e = new SimpleVector(lTIRxml->getE());
      isAllocatedIn["e"] = true;
    }
  }

  if (lTIRxml->hasB())
  {
    if (lTIRxml->hasG())
      RuntimeException::selfThrow(" FirstOrderLinearR xml constructor failed. Too many inputs: you can not give B and g or its jacobian.");

    if (lTIRxml->isBPlugin())
    {
      plugin = lTIRxml->getBPlugin();
      setComputeBFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else
    {
      B = new SimpleMatrix(lTIRxml->getB());
      isAllocatedIn["B"] = true;
    }
  }
}

// Minimum data (C, B) constructor
FirstOrderLinearR::FirstOrderLinearR(const SiconosMatrix& newC, const SiconosMatrix& newB):
  FirstOrderR("LinearR"), C(NULL), D(NULL), F(NULL), e(NULL), B(NULL),
  CPtr(NULL), DPtr(NULL), FPtr(NULL), ePtr(NULL), BPtr(NULL)
{
  C = new SimpleMatrix(newC);
  B = new SimpleMatrix(newB);
  isAllocatedIn["C"] = true;
  isAllocatedIn["B"] = true;
  isAllocatedIn["D"] = false;
  isAllocatedIn["F"] = false;
  isAllocatedIn["e"] = false;
  initPluginFlags(false);
}

// Constructor from a complete set of data
FirstOrderLinearR::FirstOrderLinearR(const SiconosMatrix& newC, const SiconosMatrix& newD,
                                     const SiconosMatrix& newF, const SimpleVector& newE,
                                     const SiconosMatrix& newB):
  FirstOrderR("LinearR"), C(NULL), D(NULL), F(NULL), e(NULL), B(NULL),
  CPtr(NULL), DPtr(NULL), FPtr(NULL), ePtr(NULL), BPtr(NULL)
{
  C = new SimpleMatrix(newC);
  isAllocatedIn["C"] = true;

  D = new SimpleMatrix(newD);
  isAllocatedIn["D"] = true;

  F = new SimpleMatrix(newF);
  isAllocatedIn["F"] = true;

  e = new SimpleVector(newE);
  isAllocatedIn["e"] = true;

  B = new SimpleMatrix(newB);
  isAllocatedIn["B"] = true;
  initPluginFlags(false);
}

// Minimum data (C, B as pointers) constructor
FirstOrderLinearR::FirstOrderLinearR(SiconosMatrix * newC, SiconosMatrix * newB):
  FirstOrderR("LinearR"), C(NULL), D(NULL), F(NULL), e(NULL), B(NULL),
  CPtr(NULL), DPtr(NULL), FPtr(NULL), ePtr(NULL), BPtr(NULL)
{
  initPluginFlags(false);
  initAllocationFlags(false);
}

// Constructor from a complete set of data
FirstOrderLinearR::FirstOrderLinearR(SiconosMatrix* newC, SiconosMatrix* newD, SiconosMatrix* newF, SimpleVector* newE, SiconosMatrix* newB):
  FirstOrderR("LinearR"), C(NULL), D(NULL), F(NULL), e(NULL), B(NULL),
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
}

void FirstOrderLinearR::initialize()
{
  FirstOrderR::initialize();

  // Check if various operators sizes are consistent.
  // Reference: interaction.
  unsigned int sizeY = interaction->getInteractionSize();
  unsigned int sizeX = interaction->getSizeOfDS();

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
  if (C == NULL)
  {
    unsigned int sizeY = interaction->getInteractionSize();
    unsigned int sizeX = interaction->getSizeOfDS();
    C = new SimpleMatrix(sizeY, sizeX);
    isAllocatedIn["C"] = true ;
  }

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
  if (D == NULL)
  {
    unsigned int sizeY = interaction->getInteractionSize();
    D = new SimpleMatrix(sizeY, sizeY);
    isAllocatedIn["D"] = true ;
  }

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
  if (F == NULL)
  {
    unsigned int sizeY = interaction->getInteractionSize();
    unsigned int sizeZ = interaction->getSizeZ();
    F = new SimpleMatrix(sizeY, sizeZ);
    isAllocatedIn["F"] = true ;
  }

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
  if (e == NULL)
  {
    unsigned int sizeY = interaction->getInteractionSize();
    e = new SimpleVector(sizeY);
    isAllocatedIn["e"] = true ;
  }

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
  if (B == NULL)
  {
    unsigned int sizeY = interaction->getInteractionSize();
    unsigned int sizeX = interaction->getSizeOfDS();
    B = new SimpleMatrix(sizeX, sizeY);
    isAllocatedIn["B"] = true ;
  }

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
    *y = prod(*C, *data["x"]);

  if (D != NULL)
    *y += prod(*D, *lambda);

  if (e != NULL)
    *y += *e;

  if (F != NULL)
    *y += prod(*F, *data["z"]);
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
  *yFree = prod(*C, *data["xFree"]);

  if (e != NULL)
    *yFree += *e;

  if (F != NULL)
    *yFree += prod(*F, *data["z"]);
}

void FirstOrderLinearR::computeInput(double time, unsigned int level)
{
  if (isPlugged["B"])
    computeB(time);

  // We get lambda of the interaction (pointers)
  SiconosVector *lambda = interaction->getLambdaPtr(level);
  *data["r"] += prod(*B, *lambda);
}

void FirstOrderLinearR::computeC(const double time)
{
  if (isPlugged["C"])
  {
    if (CPtr == NULL)
      RuntimeException::selfThrow("computeC() is not linked to a plugin function");
    unsigned int sizeY = interaction->getInteractionSize();
    unsigned int sizeX = interaction->getSizeOfDS();
    unsigned int sizeZ = interaction->getSizeZ();
    SimpleVector * zCopy = new SimpleVector(*data["z"]);
    CPtr(time, sizeY, sizeX, &(*C)(0, 0), sizeZ, &(*zCopy)(0));
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
    unsigned int sizeY = interaction->getInteractionSize();
    unsigned int sizeZ = interaction->getSizeZ();
    SimpleVector * zCopy = new SimpleVector(*data["z"]);
    DPtr(time, sizeY, &(*D)(0, 0), sizeZ, &(*zCopy)(0));
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
    unsigned int sizeY = interaction->getInteractionSize();
    unsigned int sizeZ = interaction->getSizeZ();
    SimpleVector * zCopy = new SimpleVector(*data["z"]);
    FPtr(time, sizeY, &(*F)(0, 0), sizeZ, &(*zCopy)(0));
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
    unsigned int sizeY = interaction->getInteractionSize();
    unsigned int sizeZ = interaction->getSizeZ();
    SimpleVector * zCopy = new SimpleVector(*data["z"]);
    ePtr(time, sizeY, &(*e)(0), sizeZ, &(*zCopy)(0));
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
    unsigned int sizeY = interaction->getInteractionSize();
    unsigned int sizeX = interaction->getSizeOfDS();
    unsigned int sizeZ = interaction->getSizeZ();
    SimpleVector * zCopy = new SimpleVector(*data["z"]);
    BPtr(time, sizeX, sizeY, &(*B)(0, 0), sizeZ, &(*zCopy)(0));
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

  FirstOrderLinearRXML * lTIRxml = (static_cast<FirstOrderLinearRXML*>(relationxml));
  lTIRxml->setC(*C);
  lTIRxml->setD(*D);
  lTIRxml->setF(*F);
  lTIRxml->setE(*e);
  lTIRxml->setB(*B);
}

FirstOrderLinearR* FirstOrderLinearR::convert(Relation *r)
{
  return dynamic_cast<FirstOrderLinearR*>(r);
}

