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

// --- Default (private) constructor ---
LagrangianLinearTIDS::LagrangianLinearTIDS(): LagrangianDS(), K(NULL), C(NULL), isKAllocatedIn(false), isCAllocatedIn(false)
{
  DSType = LTIDS;
}

// --- Constructor from an xml file (newNsds is optional) ---
LagrangianLinearTIDS::LagrangianLinearTIDS(DynamicalSystemXML * dsXML,  NonSmoothDynamicalSystem* newNsds):
  LagrangianDS(dsXML, newNsds), K(NULL), C(NULL),
  isKAllocatedIn(true), isCAllocatedIn(true)
{
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
}

// --- Constructor from a set of data - Mass (from a matrix), K and C ---
LagrangianLinearTIDS::LagrangianLinearTIDS(const int& newNumber, const unsigned int& newNdof,
    const SimpleVector& newQ0, const SimpleVector& newVelocity0,
    const SiconosMatrix& newMass, const SiconosMatrix& newK, const SiconosMatrix& newC):
  LagrangianDS(newNumber, newNdof, newQ0, newVelocity0, newMass), K(NULL), C(NULL), isKAllocatedIn(true), isCAllocatedIn(true)
{
  IN("LagrangianLinearTIDS::LagrangianLinearTIDS -  Constructor from a minimum set of data\n");
  if (newK.size(0) != ndof || newK.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - constructor from data, inconsistent size between K and ndof");

  if (newC.size(0) != ndof || newC.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - constructor from data, inconsistent size between C and ndof");

  K = new SimpleMatrix(newK);
  C = new SimpleMatrix(newC);
  DSType = LTIDS;
  OUT("LagrangianLinearTIDS::LagrangianLinearTIDS - Constructor from a minimum set of data\n");
}

// --- Constructor from a set of data - Mass (from plugin), K and C ---
LagrangianLinearTIDS::LagrangianLinearTIDS(const int& newNumber, const unsigned int& newNdof,
    const SimpleVector& newQ0, const SimpleVector& newVelocity0,
    const std::string& massName, const SiconosMatrix& newK, const SiconosMatrix& newC):
  LagrangianDS(newNumber, newNdof, newQ0, newVelocity0, massName), K(NULL), C(NULL), isKAllocatedIn(true), isCAllocatedIn(true)
{
  IN("LagrangianLinearTIDS::LagrangianLinearTIDS -  Constructor from a minimum set of data\n");

  if (newK.size(0) != ndof || newK.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - constructor from data, inconsistent size between K and ndof");

  if (newC.size(0) != ndof || newC.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - constructor from data, inconsistent size between C and ndof");

  K = new SimpleMatrix(newK);
  C = new SimpleMatrix(newC);
  DSType = LTIDS;
  OUT("LagrangianLinearTIDS::LagrangianLinearTIDS - Constructor from a minimum set of data\n");
}

// --- Constructor from a set of data - Mass (from a matrix), no K and no C ---
LagrangianLinearTIDS::LagrangianLinearTIDS(const int& newNumber, const unsigned int& newNdof,
    const SimpleVector& newQ0, const SimpleVector& newVelocity0,
    const SiconosMatrix& newMass):
  LagrangianDS(newNumber, newNdof, newQ0, newVelocity0, newMass), K(NULL), C(NULL), isKAllocatedIn(true), isCAllocatedIn(true)
{
  IN("LagrangianLinearTIDS::LagrangianLinearTIDS -  Constructor from a minimum set of data\n");
  DSType = LTIDS;
  OUT("LagrangianLinearTIDS::LagrangianLinearTIDS - Constructor from a minimum set of data\n");
}

// --- Constructor from a set of data - Mass (from plugin), no K and no C ---
LagrangianLinearTIDS::LagrangianLinearTIDS(const int& newNumber, const unsigned int& newNdof,
    const SimpleVector& newQ0, const SimpleVector& newVelocity0,
    const std::string& massName):
  LagrangianDS(newNumber, newNdof, newQ0, newVelocity0, massName), K(NULL), C(NULL), isKAllocatedIn(true), isCAllocatedIn(true)
{
  IN("LagrangianLinearTIDS::LagrangianLinearTIDS -  Constructor from a minimum set of data\n");
  DSType = LTIDS;
  OUT("LagrangianLinearTIDS::LagrangianLinearTIDS - Constructor from a minimum set of data\n");
}



// Copy constructor
LagrangianLinearTIDS::LagrangianLinearTIDS(const DynamicalSystem & newDS):
  LagrangianDS(newDS), K(NULL), C(NULL), isKAllocatedIn(false), isCAllocatedIn(false)
{
  if (newDS.getType() != LTIDS)
    RuntimeException::selfThrow("LagrangianLinearTIDS - copy constructor: try to copy into a LagrangianLinearTIDS a DS of type: " + newDS.getType());

  DSType = LTIDS;

  // convert newDS to lagrangianLinearTIDS by keeping const options
  const LagrangianLinearTIDS * ltids = static_cast<const LagrangianLinearTIDS*>(&newDS);

  if (ltids->getKPtr() != NULL)
  {
    K = new SimpleMatrix(ltids->getK());
    isKAllocatedIn = true;
  }

  if (ltids->getCPtr() != NULL)
  {
    C = new SimpleMatrix(ltids->getC());
    isCAllocatedIn = true;
  }
}

LagrangianLinearTIDS::~LagrangianLinearTIDS()
{
  IN("LagrangianLinearTIDS::~LagrangianLinearTIDS()\n");
  if (isKAllocatedIn) delete K;
  K = NULL;
  if (isCAllocatedIn) delete C;
  C = NULL;
  OUT("LagrangianLinearTIDS::~LagrangianLinearTIDS()\n");
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

  // compute values for mass, FInt etc ...
  isConstant["mass"] = true;

  // FInt
  * fInt = *C ** velocity + *K **q;

  SiconosVector * param;

  if (isLDSPlugin[2]) // FExt
  {
    if (computeFExtPtr == NULL)
      RuntimeException::selfThrow("LagrangianLinearTIDS initialize, computeFExt() is not linked to a plugin function");
    param = parametersList[2];
    computeFExtPtr(ndof, &time, &(*fExt)(0), &(*param)(0));
  }

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

  isDSup = true;

  // \todo: control terms handling
}


void LagrangianLinearTIDS::setK(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != ndof || newValue.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - setK: inconsistent input matrix size ");

  if (K == NULL)
  {
    K = new SimpleMatrix(newValue);
    isKAllocatedIn = true;
  }
  else
    *K = newValue;
}

void LagrangianLinearTIDS::setKPtr(SiconosMatrix *newPtr)
{
  if (newPtr->size(0) != ndof || newPtr->size(1) != ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - setKPtr: inconsistent input matrix size ");

  if (isKAllocatedIn) delete K;
  K = newPtr;
  isKAllocatedIn = false;
}

void LagrangianLinearTIDS::setC(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != ndof || newValue.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - setC: inconsistent input matrix size ");

  if (C == NULL)
  {
    C = new SimpleMatrix(newValue);
    isCAllocatedIn = true;
  }
  else
    *C = newValue;
}

void LagrangianLinearTIDS::setCPtr(SiconosMatrix *newPtr)
{
  if (newPtr->size(0) != ndof || newPtr->size(1) != ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - setCPtr: inconsistent input matrix size ");

  if (isCAllocatedIn) delete C;
  C = newPtr;
  isCAllocatedIn = false;
}

void LagrangianLinearTIDS::display() const
{
  IN("LagrangianLinearTIDS::display\n");

  LagrangianDS::display();
  cout << "===== Lagrangian Linear Time Invariant System display ===== " << endl;
  cout << "- Stiffness Matrix K : " << endl;
  if (K != NULL) K->display();
  else cout << "-> NULL" << endl;
  cout << "- Viscosity Matrix C : " << endl;
  if (C != NULL) C->display();
  else cout << "-> NULL" << endl;
  cout << "=========================================================== " << endl;
  OUT("LagrangianLinearTIDS::display\n");
}

void LagrangianLinearTIDS::computeVectorField(const double& time)
{
  // note that xDot(0) = velocity with pointer link must already be set.

  SiconosVector* tmp = xDot->getVectorPtr(1); // Pointer link!
  // Compute M-1
  map<string, SiconosMatrix*>::iterator it = workMatrix.find("inverseOfMass");
  if (it == workMatrix.end()) // if it is the first call to computeVectorField
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
  IN("LagrangianLinearTIDS::saveDSToXML\n");

  // --- general DS data---
  saveDSDataToXML();
  // --- other data  ---
  if (dsxml != NULL)
  {
    LagrangianDSXML* lgptr = static_cast <LagrangianDSXML*>(dsxml);
    lgptr->setNdof(ndof);
    lgptr->setMMatrix(mass);
    lgptr->setQ(q);
    lgptr->setQ0(q0);
    lgptr->setQMemory(qMemory);
    lgptr->setVelocity(velocity);
    lgptr->setVelocity0(velocity0);
    lgptr->setVelocityMemory(velocityMemory);

    // FExt
    if (lgptr->hasFext())
    {
      if (!lgptr->isFextPlugin())
      {
        lgptr->setFextVector(fExt);
      }
    }
    else
    {
      lgptr->setFextPlugin(fExtFunctionName);
    }
    (static_cast <LagrangianLinearTIDSXML*>(dsxml))->setK(K);
    (static_cast <LagrangianLinearTIDSXML*>(dsxml))->setC(C);
  }
  else RuntimeException::selfThrow("LagrangianLinearTIDS::saveDSToXML - object DynamicalSystemXML does not exist");
  OUT("LagrangianLinearTIDS::saveDSToXML\n");
}

LagrangianLinearTIDS* LagrangianLinearTIDS::convert(DynamicalSystem* ds)
{
  LagrangianLinearTIDS* ltids = dynamic_cast<LagrangianLinearTIDS*>(ds);
  return ltids;
}

