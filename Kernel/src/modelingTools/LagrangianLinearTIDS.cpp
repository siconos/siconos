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

void LagrangianLinearTIDS::initAllocationFlags(const bool& val)
{
  isAllocatedIn["C"] = val;
  isAllocatedIn["K"] = val;
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
  DSType = LTIDS;
  if (dsXML != NULL)
  {
    LagrangianLinearTIDSXML * lltidsxml = (static_cast <LagrangianLinearTIDSXML*>(dsxml));
    if (lltidsxml->hasK())
      K = new SimpleMatrix(lltidsxml->getK());
    if (lltidsxml->hasC())
      C = new SimpleMatrix(lltidsxml->getC());
  }
  else RuntimeException::selfThrow("LagrangianLinearTIDS::LagrangianLinearTIDS - DynamicalSystemXML paramater must not be NULL");
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// --- Constructor from a set of data - Mass (from a matrix), K and C ---
LagrangianLinearTIDS::LagrangianLinearTIDS(const int& newNumber, const unsigned int& newNdof,
    const SimpleVector& newQ0, const SimpleVector& newVelocity0,
    const SiconosMatrix& newMass, const SiconosMatrix& newK, const SiconosMatrix& newC):
  LagrangianDS(newNumber, newNdof, newQ0, newVelocity0, newMass), K(NULL), C(NULL)
{
  if (newK.size(0) != ndof || newK.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - constructor from data, inconsistent size between K and ndof");

  if (newC.size(0) != ndof || newC.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - constructor from data, inconsistent size between C and ndof");

  K = new SimpleMatrix(newK);
  C = new SimpleMatrix(newC);
  DSType = LTIDS;
  initAllocationFlags(true);
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// --- Constructor from a set of data - Mass (from plugin), K and C ---
LagrangianLinearTIDS::LagrangianLinearTIDS(const int& newNumber, const unsigned int& newNdof,
    const SimpleVector& newQ0, const SimpleVector& newVelocity0,
    const std::string& massName, const SiconosMatrix& newK, const SiconosMatrix& newC):
  LagrangianDS(newNumber, newNdof, newQ0, newVelocity0, massName), K(NULL), C(NULL)
{
  if (newK.size(0) != ndof || newK.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - constructor from data, inconsistent size between K and ndof");

  if (newC.size(0) != ndof || newC.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - constructor from data, inconsistent size between C and ndof");

  K = new SimpleMatrix(newK);
  C = new SimpleMatrix(newC);
  DSType = LTIDS;
  initAllocationFlags(true);
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// --- Constructor from a set of data - Mass (from a matrix), no K and no C ---
LagrangianLinearTIDS::LagrangianLinearTIDS(const int& newNumber, const unsigned int& newNdof,
    const SimpleVector& newQ0, const SimpleVector& newVelocity0,
    const SiconosMatrix& newMass):
  LagrangianDS(newNumber, newNdof, newQ0, newVelocity0, newMass), K(NULL), C(NULL)
{
  DSType = LTIDS;
  initAllocationFlags(false);
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}

// --- Constructor from a set of data - Mass (from plugin), no K and no C ---
LagrangianLinearTIDS::LagrangianLinearTIDS(const int& newNumber, const unsigned int& newNdof,
    const SimpleVector& newQ0, const SimpleVector& newVelocity0,
    const std::string& massName):
  LagrangianDS(newNumber, newNdof, newQ0, newVelocity0, massName), K(NULL), C(NULL)
{
  DSType = LTIDS;
  initAllocationFlags(false);
  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
}



// Copy constructor
LagrangianLinearTIDS::LagrangianLinearTIDS(const DynamicalSystem & newDS):
  LagrangianDS(newDS), K(NULL), C(NULL)
{
  if (newDS.getType() != LTIDS)
    RuntimeException::selfThrow("LagrangianLinearTIDS - copy constructor: try to copy into a LagrangianLinearTIDS a DS of type: " + newDS.getType());

  DSType = LTIDS;

  // convert newDS to lagrangianLinearTIDS by keeping const options
  const LagrangianLinearTIDS * ltids = static_cast<const LagrangianLinearTIDS*>(&newDS);

  if (ltids->getKPtr() != NULL)
  {
    K = new SimpleMatrix(ltids->getK());
    isAllocatedIn["K"] = true;
  }
  else isAllocatedIn["K"] = false;

  if (ltids->getCPtr() != NULL)
  {
    C = new SimpleMatrix(ltids->getC());
    isAllocatedIn["C"] = true;
  }
  else isAllocatedIn["C"] = false;
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

  if (isPlugin["vectorField"] || isPlugin["jacobianX"])
  {
    RuntimeException::selfThrow("LagrangianLinearTIDS::checkDynamicalSystem - vectorField and/or ist Jacobian can not be plugged for a Lagrangian system.");
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

  // FInt
  * fInt = *C ** velocity + *K **q;

  computeFExt(time);

  // vectorField
  // note that xDot(0) = velocity, with pointer link, must already be set.
  SiconosVector* tmp = xDot->getVectorPtr(1); // Pointer link!
  // Compute M-1
  workMatrix["inverseOfMass"] = new SimpleMatrix(*mass);
  workMatrix["inverseOfMass"]->PLUFactorizationInPlace();
  workMatrix["inverseOfMass"]->PLUInverseInPlace();

  *tmp = *(workMatrix["inverseOfMass"]) * (*fExt - *fInt);
  // todo: use linearSolve to avoid inversion ? Or save M-1 to avoid re-computation. See this when "M" will be added in DS or LDS.

  // jacobianX
  SiconosMatrix * tmp2;
  tmp2 = jacobianX->getBlockPtr(1, 0); // Pointer link!
  // !!! jacobian of M according to q is not take into account at the time !!!
  *tmp2 = -1 * *workMatrix["inverseOfMass"] * *K;
  tmp2 = jacobianX->getBlockPtr(1, 1); // Pointer link!
  *tmp2 = -1 * *workMatrix["inverseOfMass"] * *C;

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

void LagrangianLinearTIDS::computeVectorField(const double& time)
{
  // note that xDot(0) = velocity with pointer link must already be set.

  SiconosVector* tmp = xDot->getVectorPtr(1); // Pointer link!
  // Compute M-1
  map<string, SiconosMatrix*>::iterator it = workMatrix.find("inverseOfMass");
  if (it == workMatrix.end()) // if it is the first call to computeVectorField or computeJacobianX
  {
    workMatrix["inverseOfMass"] = new SimpleMatrix(*mass);
    workMatrix["inverseOfMass"]->PLUFactorizationInPlace();
    workMatrix["inverseOfMass"]->PLUInverseInPlace();
  }

  *tmp = *(workMatrix["inverseOfMass"]) * (*fExt - *C ** velocity - *K **q);;
  // todo: use linearSolve to avoid inversion ? Or save M-1 to avoid re-computation. See this when "M" will be added in DS or LDS.
}

void LagrangianLinearTIDS::computeJacobianX(const double& time)
{
  SiconosMatrix* tmp = jacobianX->getBlockPtr(1, 0); // Pointer link!
  // Compute M-1 if required
  map<string, SiconosMatrix*>::iterator it = workMatrix.find("inverseOfMass");
  if (it == workMatrix.end()) // if it is the first call to computeVectorField or computeJacobianX
  {
    workMatrix["inverseOfMass"] = new SimpleMatrix(*mass);
    workMatrix["inverseOfMass"]->PLUFactorizationInPlace();
    workMatrix["inverseOfMass"]->PLUInverseInPlace();
  }
  *tmp = -1 * *workMatrix["inverseOfMass"] * *K;
  tmp = jacobianX->getBlockPtr(1, 1); // Pointer link!
  *tmp = -1 * *workMatrix["inverseOfMass"] * *C;
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

