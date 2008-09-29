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

void FirstOrderLinearR::initPluginFlags(bool in)
{
  isPlugged["C"] = in;
  isPlugged["D"] = in;
  isPlugged["F"] = in;
  isPlugged["e"] = in;
  isPlugged["B"] = in ;
}

// Default (private) constructor
FirstOrderLinearR::FirstOrderLinearR(): FirstOrderR(LinearR)
{
  initPluginFlags(false);
}


// xml constructor
FirstOrderLinearR::FirstOrderLinearR(RelationXMLSPtr relxml):
  FirstOrderR(relxml, LinearR)
{
  FirstOrderLinearRXMLSPtr folrXML = boost::static_pointer_cast<FirstOrderLinearRXML>(relationxml);
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
      setComputeCFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else
    {

#ifndef WithSmartPtr
      C = new SimpleMatrix(folrXML->getC());
      isAllocatedIn["C"] = true;
#else
      C.reset(new SimpleMatrix(folrXML->getC()));
#endif

    }
  }

  if (folrXML->hasD())
  {
    if (folrXML->isDPlugin())
    {
      plugin = folrXML->getDPlugin();
      setComputeDFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else
    {

#ifndef WithSmartPtr
      D = new SimpleMatrix(folrXML->getD());
      isAllocatedIn["D"] = true;
#else
      D.reset(new SimpleMatrix(folrXML->getD()));
#endif
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

#ifndef WithSmartPtr
      F = new SimpleMatrix(folrXML->getF());
      isAllocatedIn["F"] = true;
#else
      F.reset(new SimpleMatrix(folrXML->getF()));
#endif

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

#ifndef WithSmartPtr
      e = new SimpleVector(folrXML->getE());
      isAllocatedIn["e"] = true;
#else
      e.reset(new SimpleVector(folrXML->getE()));
#endif

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

#ifndef WithSmartPtr
      B = new SimpleMatrix(folrXML->getB());
      isAllocatedIn["B"] = true;
#else
      B.reset(new SimpleMatrix(folrXML->getB()));
#endif
    }
  }
}

// Constructor with C and B plug-in names
FirstOrderLinearR::FirstOrderLinearR(const string& CName, const string& BName):
  FirstOrderR(LinearR)
{
  initPluginFlags(false);

  setComputeCFunction(cShared.getPluginName(CName), cShared.getPluginFunctionName(CName));
  setComputeBFunction(cShared.getPluginName(BName), cShared.getPluginFunctionName(BName));
}

// Constructor with e plug-in name
FirstOrderLinearR::FirstOrderLinearR(const string& EName):
  FirstOrderR(LinearR)
{
  initPluginFlags(false);

  setComputeEFunction(cShared.getPluginName(EName), cShared.getPluginFunctionName(EName));
}

// Constructor from a complete set of data (plugin)
FirstOrderLinearR::FirstOrderLinearR(const string& CName, const string& DName, const string& FName, const string& EName, const string& BName):
  FirstOrderR(LinearR)
{
  initPluginFlags(true);

  setComputeCFunction(cShared.getPluginName(CName), cShared.getPluginFunctionName(CName));
  setComputeDFunction(cShared.getPluginName(DName), cShared.getPluginFunctionName(DName));
  setComputeFFunction(cShared.getPluginName(FName), cShared.getPluginFunctionName(FName));
  setComputeEFunction(cShared.getPluginName(EName), cShared.getPluginFunctionName(EName));
  setComputeBFunction(cShared.getPluginName(BName), cShared.getPluginFunctionName(BName));
}

// Minimum data (C, B as pointers) constructor
FirstOrderLinearR::FirstOrderLinearR(SiconosMatrixSPtr newC, SiconosMatrixSPtr newB):
  FirstOrderR(LinearR)
{
  initPluginFlags(false);

  C = newC;

}

// Constructor from a complete set of data
FirstOrderLinearR::FirstOrderLinearR(SiconosMatrixSPtr newC, SiconosMatrixSPtr newD, SiconosMatrixSPtr newF, SiconosVectorSPtr newE, SiconosMatrixSPtr newB):
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
{

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
#ifndef WithSmartPtr
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
  workZ = new SimpleVector(sizeZ);
#else
  if (! C && isPlugged["C"])
  {
    C.reset(new SimpleMatrix(sizeY, sizeX));
  }

  if (! D && isPlugged["D"])
  {
    D.reset(new SimpleMatrix(sizeY, sizeY));
  }

  if (! F && isPlugged["F"])
  {
    F.reset(new SimpleMatrix(sizeY, sizeZ));
  }

  if (! e && isPlugged["e"])
  {
    e.reset(new SimpleVector(sizeY));
  }

  if (! B && isPlugged["B"])
  {
    B.reset(new SimpleMatrix(sizeX, sizeY));
  }
  workZ.reset(new SimpleVector(sizeZ));
#endif
}

// setters

void FirstOrderLinearR::setC(const SiconosMatrix& newValue)
{
  if (! C)
  {
#ifndef WithSmartPtr
    C = new SimpleMatrix(newValue);
    isAllocatedIn["C"] = true;
#else
    C.reset(new SimpleMatrix(newValue));
#endif

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

void FirstOrderLinearR::setCPtr(SiconosMatrixSPtr newPtr)
{

#ifndef WithSmartPtr
  if (isAllocatedIn["C"]) delete C;
  isAllocatedIn["C"] = false;
#endif

  C = newPtr;
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
  if (! D)
  {

#ifndef WithSmartPtr
    D = new SimpleMatrix(newValue);
    isAllocatedIn["D"] = true;
#else
    D.reset(new SimpleMatrix(newValue));
#endif

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

void FirstOrderLinearR::setDPtr(SiconosMatrixSPtr newPtr)
{

#ifndef WithSmartPtr
  if (isAllocatedIn["D"]) delete D;
  isAllocatedIn["D"] = false;
#endif

  D = newPtr;
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
  if (! F)
  {

#ifndef WithSmartPtr
    F = new SimpleMatrix(newValue);
    isAllocatedIn["F"] = true;
#else
    F.reset(new SimpleMatrix(newValue));
#endif

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

void FirstOrderLinearR::setFPtr(SiconosMatrixSPtr newPtr)
{

#ifndef WithSmartPtr
  if (isAllocatedIn["F"]) delete F;
  isAllocatedIn["F"] = false;
#endif

  F = newPtr;
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
  if (! e)
  {

#ifndef WithSmartPtr
    e = new SimpleVector(newValue);
    isAllocatedIn["e"] = true;
#else
    e.reset(new SimpleVector(newValue));
#endif

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

void FirstOrderLinearR::setEPtr(SiconosVectorSPtr newPtr)
{

#ifndef WithSmartPtr
  if (isAllocatedIn["e"]) delete e;
  isAllocatedIn["e"] = false;
#endif

  e = newPtr;
  isPlugged["e"] = false;
}
void FirstOrderLinearR::setComputeEFunction(FOVecPtr ptrFunct)
{
  ePtr = ptrFunct;
  isPlugged["e"] = true;
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
  if (! B)
  {

#ifndef WithSmartPtr
    B = new SimpleMatrix(newValue);
    isAllocatedIn["B"] = true;
#else
    B.reset(new SimpleMatrix(newValue));
#endif

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

void FirstOrderLinearR::setBPtr(SiconosMatrixSPtr newPtr)
{

#ifndef WithSmartPtr
  if (isAllocatedIn["B"]) delete B;
  isAllocatedIn["B"] = false;
#endif

  B = newPtr;
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
  SiconosVectorSPtr y = interaction->getYPtr(0);
  SiconosVectorSPtr lambda = interaction->getLambdaPtr(0);

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
  if (isPlugged["B"])
    computeB(time);

  // We get lambda of the interaction (pointers)
  SiconosVectorSPtr lambda = interaction->getLambdaPtr(level);
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
    *workZ = *data["z"];
    CPtr(time, sizeY, sizeX, &(*C)(0, 0), sizeZ, &(*workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *data["z"] = *workZ;
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
    *workZ = *data["z"];
    DPtr(time, sizeY, &(*D)(0, 0), sizeZ, &(*workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *data["z"] = *workZ;
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
    *workZ = *data["z"];
    FPtr(time, sizeY, &(*F)(0, 0), sizeZ, &(*workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *data["z"] = *workZ;
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
    *workZ = *data["z"];
    ePtr(time, sizeY, &(*e)(0), sizeZ, &(*workZ)(0));
    // Copy data that might have been changed in the plug-in call.
    *data["z"] = *workZ;
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
  if (relationxml == NULL)
    RuntimeException::selfThrow("FirstOrderLinearR::saveRelationToXML, no xml object found.");

  FirstOrderLinearRXMLSPtr folrXML = (boost::static_pointer_cast<FirstOrderLinearRXML>(relationxml));
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

