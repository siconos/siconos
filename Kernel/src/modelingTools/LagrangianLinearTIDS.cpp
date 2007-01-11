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
#include "LagrangianLinearTIDS.h"
using namespace std;

void LagrangianLinearTIDS::connectToDS()
{
  // This function set properly links with members of class DynamicalSystem and LagrangianDS and allocate memory if necessary.
  // Warning: computation of jacobianX or FInt is done in initialize, not here, since this function is called at the
  // end of constructor, and that operators are likely to be updated with set functions or else.

  n = 2 * ndof;

  // x and related vectors
  x[0] = new BlockVector(q[0], q[1]);
  isAllocatedIn["x"] = true;
  x0 = new BlockVector(q0, velocity0);
  isAllocatedIn["x0"] = true;
  xFree = new BlockVector(qFree[0], qFree[1]);
  isAllocatedIn["xFree"] = true;

  // r
  r = new BlockVector(NULL, p[2]);
  isAllocatedIn["r"] = true;
  r->zero();

  // rhs // add conditional allocation ??
  x[1] = new BlockVector(q[1], q[2]);
  isAllocatedIn["rhs"] = true;

  // jacobianXRhs

  isPlugin["f"] = false;
  isPlugin["jacobianXF"] = false;

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

    jacobianXRhs = new BlockMatrix(workMatrix["zero-matrix"], workMatrix["Id-matrix"], workMatrix["jacob-block10"], workMatrix["jacob-block11"]);
    isAllocatedIn["jacobianXF"] = true;
  }

  if (K != NULL)
    jacobianFInt[0] = K;

  if (C != NULL)
    jacobianFInt[0] = C;

  // Remark: all other operators are NULL pointers (NNL, jacobian ...). See LagrangianTIDS description in User Manual for details.

}
void LagrangianLinearTIDS::initAllocationFlags(const bool in)
{
  isAllocatedIn["C"] = in; // By default, neither K nor C are allocated in, whatever the value of in is.
  isAllocatedIn["K"] = in;
}

void LagrangianLinearTIDS::initPluginFlags(const bool val)
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
  LagrangianLinearTIDSXML * lltidsxml = (static_cast <LagrangianLinearTIDSXML*>(dsxml));
  if (lltidsxml->hasK())
    K = new SimpleMatrix(lltidsxml->getK());
  else isAllocatedIn["K"] = false;
  if (lltidsxml->hasC())
    C = new SimpleMatrix(lltidsxml->getC());
  else isAllocatedIn["C"] = false;
}

// --- Constructor from a set of data - Mass, K and C ---
LagrangianLinearTIDS::LagrangianLinearTIDS(const int newNumber, const SimpleVector& newQ0, const SimpleVector& newVelocity0,
    const SiconosMatrix& newMass, const SiconosMatrix& newK, const SiconosMatrix& newC):
  LagrangianDS(newNumber, newQ0, newVelocity0, newMass), K(NULL), C(NULL)
{
  if (newK.size(0) != ndof || newK.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - constructor from data, inconsistent size between K and ndof");

  if (newC.size(0) != ndof || newC.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - constructor from data, inconsistent size between C and ndof");

  initPluginFlags(false);

  K = new SimpleMatrix(newK);
  C = new SimpleMatrix(newC);

  DSType = LTIDS;
  initAllocationFlags(true);
}

// --- Constructor from a set of data - Mass, no K and no C ---
LagrangianLinearTIDS::LagrangianLinearTIDS(const int newNumber, const SimpleVector& newQ0, const SimpleVector& newVelocity0,
    const SiconosMatrix& newMass):
  LagrangianDS(newNumber, newQ0, newVelocity0, newMass), K(NULL), C(NULL)
{
  initPluginFlags(false);

  DSType = LTIDS;
  initAllocationFlags(false);
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
  q[0] = new SimpleVector(ltids->getQ());
  if (ltids->getQMemoryPtr() != NULL)
    qMemory = new SiconosMemory(ltids->getQMemory());
  else isAllocatedIn["qMemory"] = false;
  velocity0 = new SimpleVector(ltids->getVelocity0());
  q[1]  = new SimpleVector(ltids->getVelocity0());
  if (ltids->getVelocityMemoryPtr() != NULL)
    velocityMemory = new SiconosMemory(ltids->getVelocityMemory());
  else isAllocatedIn["velocityMemory"] = false;

  p.resize(3, NULL);

  for (unsigned int i = 0; i < 3; ++i)
  {
    if (ltids->getPPtr(i) != NULL)
    {
      p[i] = new SimpleVector(ltids->getP(i));
      string stringValue;
      stringstream sstr;
      sstr << i;
      sstr >> stringValue;
      stringValue = "p" + stringValue;
      isAllocatedIn[stringValue] = true;
    }
  }

  if (ltids->getFExtPtr() != NULL)
  {
    fExt = new SimpleVector(ltids->getFExt());
    isAllocatedIn["fExt"] = true;
  }

  if (isPlugin["fExt"])
  {
    string functionName, pluginPath;
    pluginNames["fExt"] = ltids -> getFunctionName("fExt");
    functionName = cShared.getPluginFunctionName(pluginNames["fExt"]);
    pluginPath  = cShared.getPluginName(pluginNames["fExt"]);
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

void LagrangianLinearTIDS::initialize(const string& simulationType, double time, unsigned int sizeOfMemory)
{
  initFreeVectors(simulationType);

  initP(simulationType);

  // Set variables of top-class DynamicalSystem
  connectToDS(); // note that connection can not be done during constructor call, since user can complete the ds after (add plugin or anything else).
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;

  // Warning: ->zero for free vectors must be done before q/q[1] initialization, because for EventDriven, qFree[0]/vFree = q/v (pointers equality)

  // set q and q[1] to q0 and velocity0
  *q[0] = *q0;
  *q[1] = *velocity0;

  // Initialize memory vectors
  initMemory(sizeOfMemory);

  bool flag = false;

  // === fExt, fInt, rhs and jacobianXRhs ===

  // right-hand side. vField corresponds to the second derivative of q.
  SiconosVector* vField = x[1]->getVectorPtr(1); // Pointer link!
  vField ->zero();

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
    *fInt += prod(*K, *q[0]);
    *workMatrix["jacob-block10"] -= *workMatrix["inverseOfMass"] * *K;
    flag = true;
  }
  if (C != NULL)
  {
    *fInt += prod(*C, *q[1]);
    *workMatrix["jacob-block11"] -= *workMatrix["inverseOfMass"] * *C;
    flag = true;
  }

  if (flag) // flag == true means that among FExt and Fint, one at least is not Null.
  {
    if (fInt != NULL)
      *vField -= *fInt;
    *vField = prod(*(workMatrix["inverseOfMass"]), *vField);
  }

  // todo: use linearSolve to avoid inversion ? Or save M-1 to avoid re-computation. See this when "M" will be added in DS or LDS.
  // \todo: control terms handling
}

void LagrangianLinearTIDS::update(const double time)
{
  if (fExt != NULL)
    computeFExt(time);
  computeRhs(time);
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

void LagrangianLinearTIDS::computeRhs(const double time, const bool)
{
  // second argument is useless but present because of top-class function overloading.
  // note that rhs(0) = q[1] with pointer link must already be set.
  SiconosVector* vField = x[1]->getVectorPtr(1); // Pointer link!

  *vField = *p[2]; // Warning: p update is done in Interactions/Relations

  // Compute M-1 if necessary - Only in the case where initialize has not been called before.
  map<string, SiconosMatrix*>::iterator it = workMatrix.find("inverseOfMass");
  if (it == workMatrix.end())
  {
    workMatrix["inverseOfMass"] = new SimpleMatrix(*mass);
    workMatrix["inverseOfMass"]->PLUFactorizationInPlace();
    workMatrix["inverseOfMass"]->PLUInverseInPlace();
  }

  if (fExt != NULL || K != NULL || C != NULL) // else rhs = p
  {
    if (fExt != NULL)
      *vField += *fExt; // This supposes that fExt is up to date!!

    if (K != NULL)
      *vField -= prod(*K, *q[0]);

    if (C != NULL)
      *vField -= prod(*C, *q[1]);
  }

  *vField = prod(*(workMatrix["inverseOfMass"]), *vField);
}

void LagrangianLinearTIDS::computeJacobianXRhs(const double time, const bool)
{
  // Nothing to be done since jacobianXRhs is not time-dependent.
  // But this function is required, since it is called from Lsodar (if not present, the one of LagrangianDS will be called)
}

void LagrangianLinearTIDS::saveDSToXML()
{
  // --- general DS data---
  saveDSDataToXML();
  // --- other data  ---
  if (dsxml != NULL)
  {
    LagrangianDSXML* lgptr = static_cast <LagrangianDSXML*>(dsxml);
    lgptr->setMMatrix(*mass);
    lgptr->setQ(*q[0]);
    lgptr->setQ0(*q0);
    lgptr->setQMemory(*qMemory);
    lgptr->setVelocity(*q[1]);
    lgptr->setVelocity0(*velocity0);
    lgptr->setVelocityMemory(*velocityMemory);

    // FExt
    if (lgptr->hasFExt())
    {
      if (!lgptr->isFExtPlugin())
      {
        lgptr->setFExtVector(*fExt);
      }
    }
    else
    {
      lgptr->setFExtPlugin(pluginNames["fExt"]);
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

