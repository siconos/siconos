/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
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
#include "LagrangianLinearTIDS.h"
using namespace std;

void LagrangianLinearTIDS::connectToDS()
{
  // This function set properly links with members of class DynamicalSystem and LagrangianDS and allocate memory if necessary.
  // Warning: computation of jacobianX or FInt is done in initialize, not here, since this function is called at the
  // end of constructor, and that operators are likely to be updated with set functions or else.

  if (K != NULL || C != NULL)
  {
    // FInt
    fInt = new SimpleVector(ndof);
    isAllocatedIn["fInt"] =  true;

    // jacobianXF
    workMatrix["zero-matrix"] = new SimpleMatrix(ndof, ndof);
    workMatrix["Id-matrix"] = new SimpleMatrix(ndof, ndof);
    workMatrix["jacob-block10"] = new SimpleMatrix(ndof, ndof);
    workMatrix["jacob-block11"] = new SimpleMatrix(ndof, ndof);
    workMatrix["Id-matrix"]->eye();

    jacobianXF = new BlockMatrix(workMatrix["zero-matrix"], workMatrix["Id-matrix"], workMatrix["jacob-block10"], workMatrix["jacob-block11"]);
    isAllocatedIn["jacobianXF"] = true;
  }

  if (K != NULL)
    jacobianQFInt = K;

  if (C != NULL)
    jacobianVelocityFInt = C;

  // Remark: all other operators are NULL pointers (NNL, jacobian ...). See LagrangianTIDS description in User Manual for details.

}
void LagrangianLinearTIDS::initAllocationFlags(const bool& in)
{
  isAllocatedIn["C"] = in; // By default, neither K nor C are allocated in, whatever the value of in is.
  isAllocatedIn["K"] = in;
}

void LagrangianLinearTIDS::initPluginFlags(const bool& val)
{
  // For LagrangianLinearTIDS, the only authorized plug-in is for fExt. All other plug-in flags are set to false.
  isPlugin["mass"] = false;
  isPlugin["fInt"] = false;
  isPlugin["NNL"] = false;
  isPlugin["jacobianQFInt"] = false;
  isPlugin["jacobianVelocityFInt"] = false;
  isPlugin["jacobianQNNL"] = false;
  isPlugin["jacobianVelocityNNL"] = false;
}

// --- Default (private) constructor ---
LagrangianLinearTIDS::LagrangianLinearTIDS(): LagrangianDS(), K(NULL), C(NULL)
{
  initAllocationFlags(false);
  DSType = LTIDS;
}

// --- Constructor from an xml file (newNsds is optional) ---
LagrangianLinearTIDS::LagrangianLinearTIDS(DynamicalSystemXML * dsXML,  NonSmoothDynamicalSystem* newNsds):
  LagrangianDS(dsXML, newNsds), K(NULL), C(NULL)
{
  initAllocationFlags(true);
  initPluginFlags(false);

  DSType = LTIDS;
  if (dsXML != NULL)
  {
    LagrangianLinearTIDSXML * lltidsxml = (static_cast <LagrangianLinearTIDSXML*>(dsxml));
    if (lltidsxml->hasK())
      K = new SimpleMatrix(lltidsxml->getK());
    else isAllocatedIn["K"] = false;
    if (lltidsxml->hasC())
      C = new SimpleMatrix(lltidsxml->getC());
    else isAllocatedIn["C"] = false;
  }
  else RuntimeException::selfThrow("LagrangianLinearTIDS::LagrangianLinearTIDS - DynamicalSystemXML paramater must not be NULL");

  connectToDS();// set connections with FInt and jacobianX

  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// --- Constructor from a set of data - Mass, K and C ---
LagrangianLinearTIDS::LagrangianLinearTIDS(const int& newNumber, const unsigned int& newNdof,
    const SimpleVector& newQ0, const SimpleVector& newVelocity0,
    const SiconosMatrix& newMass, const SiconosMatrix& newK, const SiconosMatrix& newC):
  LagrangianDS(newNumber, newNdof, newQ0, newVelocity0, newMass), K(NULL), C(NULL)
{
  if (newK.size(0) != ndof || newK.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - constructor from data, inconsistent size between K and ndof");

  if (newC.size(0) != ndof || newC.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - constructor from data, inconsistent size between C and ndof");

  initPluginFlags(false);

  K = new SimpleMatrix(newK);
  C = new SimpleMatrix(newC);
  connectToDS(); // set connections with FInt and jacobianX

  DSType = LTIDS;
  initAllocationFlags(true);
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// --- Constructor from a set of data - Mass, no K and no C ---
LagrangianLinearTIDS::LagrangianLinearTIDS(const int& newNumber, const unsigned int& newNdof,
    const SimpleVector& newQ0, const SimpleVector& newVelocity0,
    const SiconosMatrix& newMass):
  LagrangianDS(newNumber, newNdof, newQ0, newVelocity0, newMass), K(NULL), C(NULL)
{
  initPluginFlags(false);

  DSType = LTIDS;
  initAllocationFlags(false);
  // no connections with FInt or jacobianX since K and C are NULL.

  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// Copy constructor
LagrangianLinearTIDS::LagrangianLinearTIDS(const DynamicalSystem & newDS):
  LagrangianDS(), K(NULL), C(NULL) // warning: call default constructor of LagrangianDS !!
{
  if (newDS.getType() != LTIDS)
    RuntimeException::selfThrow("LagrangianLinearTIDS - copy constructor: try to copy into a LagrangianLinearTIDS a DS of type: " + newDS.getType());

  DSType = LTIDS;
  id = "copy";
  stepsInMemory = newDS.getStepsInMemory();

  // === convert newDS to lagrangianLinearTIDS and keeps const options ===
  const LagrangianLinearTIDS * ltids = static_cast<const LagrangianLinearTIDS*>(&newDS);

  isPlugin = ltids ->getIsPlugin();
  initPluginFlags(false); // set all plug-in to false but not fExt!
  setParameters(newDS.getParameters());   // Copy !!

  // === copy LagrangianDS up-class members ===
  LagrangianDS::initAllocationFlags(); // default

  ndof = ltids->getNdof();
  mass = new SimpleMatrix(ltids->getMass());
  q0 = new SimpleVector(ltids->getQ0());
  q = new SimpleVector(ltids->getQ());
  qFree = new SimpleVector(ltids->getQFree());
  if (ltids->getQMemoryPtr() != NULL)
    qMemory = new SiconosMemory(ltids->getQMemory());
  else isAllocatedIn["qMemory"] = false;
  velocity0 = new SimpleVector(ltids->getVelocity0());
  velocity  = new SimpleVector(ltids->getVelocity0());
  velocityFree = new SimpleVector(ltids->getVelocityFree());
  if (ltids->getVelocityMemoryPtr() != NULL)
    velocityMemory = new SiconosMemory(ltids->getVelocityMemory());
  else isAllocatedIn["velocityMemory"] = false;
  p = new SimpleVector(ltids->getP());
  if (ltids->getFExtPtr() != NULL)
  {
    fExt = new SimpleVector(ltids->getFExt());
    isAllocatedIn["fExt"] = true;
  }

  if (isPlugin["fExt"])
  {
    string functionName, pluginPath;
    fExtFunctionName = ltids -> getFExtFunctionName();
    functionName = cShared.getPluginFunctionName(fExtFunctionName);
    pluginPath  = cShared.getPluginName(fExtFunctionName);
    setComputeFExtFunction(pluginPath, functionName);
  }

  // === copy LagrangianLinearTIDS members ===
  initAllocationFlags(false);
  if (ltids->getKPtr() != NULL)
  {
    K = new SimpleMatrix(ltids->getK());
    isAllocatedIn["K"] = true;
  }

  if (ltids->getCPtr() != NULL)
  {
    C = new SimpleMatrix(ltids->getC());
    isAllocatedIn["C"] = true;
  }

  connectToDS(); // set connections with FInt and jacobianX

  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

LagrangianLinearTIDS::~LagrangianLinearTIDS()
{
  if (isAllocatedIn["K"]) delete K;
  K = NULL;
  if (isAllocatedIn["C"]) delete C;
  C = NULL;
}


bool LagrangianLinearTIDS::checkDynamicalSystem()
{
  bool output = true;
  // ndof
  if (ndof == 0)
  {
    RuntimeException::selfThrow("LagrangianLinearTIDS::checkDynamicalSystem - number of degrees of freedom is equal to 0.");
    output = false;
  }

  // q0 and velocity0 != NULL
  if (q0 == NULL || velocity0 == NULL)
  {
    RuntimeException::selfThrow("LagrangianLinearTIDS::checkDynamicalSystem - initial conditions are badly set.");
    output = false;
  }

  // Mass
  if (mass == NULL)
  {
    RuntimeException::selfThrow("LagrangianLinearTIDS::checkDynamicalSystem - Mass is not set.");
    output = false;
  }

  // fInt, NNL and their Jacobian
  if (isPlugin["fInt"] || isPlugin["jacobianQFInt"] || isPlugin["jacobianVelocityFInt"])
    // ie if fInt is defined and not constant => its Jacobian must be defined (but not necessarily plugged)
  {
    RuntimeException::selfThrow("LagrangianLinearTIDS::checkDynamicalSystem - Fint and/or their Jacobian are plugged, which should not be.");
    output = false;
  }

  if (isPlugin["NNL"] || isPlugin["jacobianQNNL"] || isPlugin["jacobianVelocityNNL"])
    // ie if fInt is defined and not constant => its Jacobian must be defined (but not necessarily plugged)
  {
    RuntimeException::selfThrow("LagrangianLinearTIDS::checkDynamicalSystem -NNl and/or their Jacobian are plugged, which should not be.");
    output = false;
  }

  if (isPlugin["f"] || isPlugin["jacobianXF"])
  {
    RuntimeException::selfThrow("LagrangianLinearTIDS::checkDynamicalSystem - f and/or its Jacobian can not be plugged for a Lagrangian system.");
    output = false;
  }
  return output;
}

void LagrangianLinearTIDS::initialize(const double& time, const unsigned int& sizeOfMemory)
{
  // set q and velocity to q0 and velocity0
  *q = *q0;
  *velocity = *velocity0;

  // reset r and free vectors
  qFree->zero();
  velocityFree->zero();
  p->zero();

  // Initialize memory vectors
  initMemory(sizeOfMemory);

  bool flag = false;

  // === fExt, fInt, rhs and jacobianXRhs ===

  // right-hand side
  SiconosVector* vField = rhs->getVectorPtr(1); // Pointer link!
  vField ->zero();

  // jacobianX:
  if (K != NULL || C != NULL)
  {
    workMatrix["jacob-block10"]->zero(); // = -M^-1 K
    workMatrix["jacob-block11"]->zero(); // = -M^-1 C
  }

  // Compute M-1 if necessary
  if (fExt != NULL || K != NULL || C != NULL)
  {
    workMatrix["inverseOfMass"] = new SimpleMatrix(*mass);
    workMatrix["inverseOfMass"]->PLUFactorizationInPlace();
    workMatrix["inverseOfMass"]->PLUInverseInPlace();
  }

  // fExt and update rhs
  if (fExt != NULL)
  {
    computeFExt(time);
    *vField += *fExt;
    flag = true;
  }

  // fInt and update jacobianXF
  if (K != NULL)
  {
    *fInt += *K **q;
    *workMatrix["jacob-block10"] -= *workMatrix["inverseOfMass"] * *K;
    flag = true;
  }
  if (C != NULL)
  {
    *fInt += *C ** velocity;
    *workMatrix["jacob-block11"] -= *workMatrix["inverseOfMass"] * *C;
    flag = true;
  }

  if (flag)
  {
    if (fInt != NULL)
      *vField -= *fInt;
    *vField = *(workMatrix["inverseOfMass"])* *vField;
  }

  // todo: use linearSolve to avoid inversion ? Or save M-1 to avoid re-computation. See this when "M" will be added in DS or LDS.
  // \todo: control terms handling
}


void LagrangianLinearTIDS::setK(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != ndof || newValue.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - setK: inconsistent input matrix size ");

  if (K == NULL)
  {
    K = new SimpleMatrix(newValue);
    isAllocatedIn["K"] = true;
  }
  else
    *K = newValue;
}

void LagrangianLinearTIDS::setKPtr(SiconosMatrix *newPtr)
{
  if (newPtr->size(0) != ndof || newPtr->size(1) != ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - setKPtr: inconsistent input matrix size ");

  if (isAllocatedIn["K"]) delete K;
  K = newPtr;
  isAllocatedIn["K"] = false;
}

void LagrangianLinearTIDS::setC(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != ndof || newValue.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - setC: inconsistent input matrix size ");

  if (C == NULL)
  {
    C = new SimpleMatrix(newValue);
    isAllocatedIn["C"] = true;
  }
  else
    *C = newValue;
}

void LagrangianLinearTIDS::setCPtr(SiconosMatrix *newPtr)
{
  if (newPtr->size(0) != ndof || newPtr->size(1) != ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - setCPtr: inconsistent input matrix size ");

  if (isAllocatedIn["C"]) delete C;
  C = newPtr;
  isAllocatedIn["C"] = false;
}

void LagrangianLinearTIDS::display() const
{
  LagrangianDS::display();
  cout << "===== Lagrangian Linear Time Invariant System display ===== " << endl;
  cout << "- Stiffness Matrix K : " << endl;
  if (K != NULL) K->display();
  else cout << "-> NULL" << endl;
  cout << "- Viscosity Matrix C : " << endl;
  if (C != NULL) C->display();
  else cout << "-> NULL" << endl;
  cout << "=========================================================== " << endl;
}

void LagrangianLinearTIDS::computeRhs(const double& time, const bool &)
{
  // second argument is useless but present because of top-class function overloading.

  // note that xDot(0) = velocity with pointer link must already be set.
  SiconosVector* vField = rhs->getVectorPtr(1); // Pointer link!
  vField->zero();
  if (fExt != NULL || K != NULL || C != NULL) // else rhs = 0
  {
    // Compute M-1 if necessary - Only in the case where initialize has not been called before.
    map<string, SiconosMatrix*>::iterator it = workMatrix.find("inverseOfMass");
    if (it == workMatrix.end())
    {
      workMatrix["inverseOfMass"] = new SimpleMatrix(*mass);
      workMatrix["inverseOfMass"]->PLUFactorizationInPlace();
      workMatrix["inverseOfMass"]->PLUInverseInPlace();
    }
    bool flag = false;
    if (fExt != NULL)
    {
      *vField += *fExt;
      flag = true;
    }
    if (K != NULL)
    {
      *vField -= *K**q;
      flag = true;
    }
    if (C != NULL)
    {
      *vField -= *C**velocity;
      flag = true;
    }
    if (flag)
      *vField = *(workMatrix["inverseOfMass"])**vField;
  }
}

void LagrangianLinearTIDS::computeJacobianXRhs(const double& time, const bool &)
{
  // second argument is useless but present because of top-class function overloading.
  if (K != NULL || C != NULL) // else jacobianX = 0
  {
    // Compute M-1 if required - Only in the case where initialize has not been called before.
    map<string, SiconosMatrix*>::iterator it = workMatrix.find("inverseOfMass");
    if (it == workMatrix.end())
    {
      workMatrix["inverseOfMass"] = new SimpleMatrix(*mass);
      workMatrix["inverseOfMass"]->PLUFactorizationInPlace();
      workMatrix["inverseOfMass"]->PLUInverseInPlace();
    }

    if (K != NULL)
    {
      workMatrix["jacob-block10"]->zero(); // = -M^-1 K
      *workMatrix["jacob-block10"] -= *workMatrix["inverseOfMass"] * *K;
    }
    if (C != NULL)
    {
      workMatrix["jacob-block11"]->zero(); // = -M^-1 C
      *workMatrix["jacob-block11"] -= *workMatrix["inverseOfMass"] * *C;
    }
  }
}

void LagrangianLinearTIDS::saveDSToXML()
{
  // --- general DS data---
  saveDSDataToXML();
  // --- other data  ---
  if (dsxml != NULL)
  {
    LagrangianDSXML* lgptr = static_cast <LagrangianDSXML*>(dsxml);
    lgptr->setNdof(ndof);
    lgptr->setMMatrix(*mass);
    lgptr->setQ(*q);
    lgptr->setQ0(*q0);
    lgptr->setQMemory(*qMemory);
    lgptr->setVelocity(*velocity);
    lgptr->setVelocity0(*velocity0);
    lgptr->setVelocityMemory(*velocityMemory);

    // FExt
    if (lgptr->hasFext())
    {
      if (!lgptr->isFextPlugin())
      {
        lgptr->setFextVector(*fExt);
      }
    }
    else
    {
      lgptr->setFextPlugin(fExtFunctionName);
    }
    (static_cast <LagrangianLinearTIDSXML*>(dsxml))->setK(*K);
    (static_cast <LagrangianLinearTIDSXML*>(dsxml))->setC(*C);
  }
  else RuntimeException::selfThrow("LagrangianLinearTIDS::saveDSToXML - object DynamicalSystemXML does not exist");
}

LagrangianLinearTIDS* LagrangianLinearTIDS::convert(DynamicalSystem* ds)
{
  LagrangianLinearTIDS* ltids = dynamic_cast<LagrangianLinearTIDS*>(ds);
  return ltids;
}

