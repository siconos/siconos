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
#include "LagrangianLinearTIDS.h"
#include "LagrangianLinearTIDSXML.h"
#include "BlockMatrix.h"

using namespace std;
using namespace DS;

// --- Default (private) constructor ---
LagrangianLinearTIDS::LagrangianLinearTIDS(): LagrangianDS()
{
  DSType = LLTIDS;
}

// --- Constructor from an xml file---
LagrangianLinearTIDS::LagrangianLinearTIDS(SP::DynamicalSystemXML dsxml): LagrangianDS(dsxml)
{
  SP::LagrangianLinearTIDSXML lltidsxml = boost::static_pointer_cast <LagrangianLinearTIDSXML>(dsxml);

  // If Fint or NNL is given: ignored.
  if (lltidsxml->hasFInt() ||  lltidsxml->hasNNL())
    cout << "!!!!! Warning : LagrangianLinearTIDS: xml constructor, Fint or NNL input will be ignored in xml file." << endl;

  // K and C

  if (lltidsxml->hasK())
  {
    K.reset(new SimpleMatrix(lltidsxml->getK()));
  }
  if (lltidsxml->hasC())
  {
    C.reset(new SimpleMatrix(lltidsxml->getC()));
  }
}

// --- Constructor from a set of data - Mass, K and C ---
LagrangianLinearTIDS::LagrangianLinearTIDS(SP::SimpleVector newQ0, SP::SimpleVector newVelocity0,
    SP::SiconosMatrix newMass,  SP::SiconosMatrix newK, SP::SiconosMatrix newC):
  LagrangianDS(newQ0, newVelocity0, newMass)
{
  assert((newK->size(0) == ndof && newK->size(1) == ndof) &&
         "LagrangianLinearTIDS - constructor from data, inconsistent size between K and ndof");

  assert((newC->size(0) == ndof && newC->size(1) == ndof) &&
         "LagrangianLinearTIDS - constructor from data, inconsistent size between C and ndof");

  K = newK;
  C = newC;

  DSType = LLTIDS;
}

// --- Constructor from a set of data - Mass, no K and no C ---
LagrangianLinearTIDS::LagrangianLinearTIDS(SP::SimpleVector newQ0, SP::SimpleVector newVelocity0,
    SP::SiconosMatrix newMass):
  LagrangianDS(newQ0, newVelocity0, newMass)
{
  DSType = LLTIDS;
}
LagrangianLinearTIDS::LagrangianLinearTIDS(const SimpleVector& newQ0, const SimpleVector& newVelocity0, const SiconosMatrix& newMass):
  LagrangianDS(createSPtrSiconosVector((SimpleVector&)newQ0), createSPtrSiconosVector((SimpleVector&)newVelocity0), createSPtrSiconosMatrix((SimpleMatrix&)newMass))
{
  DSType = LLTIDS;
}
LagrangianLinearTIDS::~LagrangianLinearTIDS()
{}


bool LagrangianLinearTIDS::checkDynamicalSystem()
{
  bool output = true;
  // ndof
  if (ndof == 0)
  {
    RuntimeException::selfThrow("LagrangianLinearTIDS::checkDynamicalSystem - number of degrees of freedom is equal to 0.");
    output = false;
  }

  // q0 and velocity0
  if (! q0 || ! velocity0)
  {
    RuntimeException::selfThrow("LagrangianLinearTIDS::checkDynamicalSystem - initial conditions are badly set.");
    output = false;
  }

  // Mass
  if (! mass)
  {
    RuntimeException::selfThrow("LagrangianLinearTIDS::checkDynamicalSystem - Mass is not set.");
    output = false;
  }

  if (!output) cout << "LagrangianLinearTIDS Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
  return output;
}

void LagrangianLinearTIDS::initRhs(double time)
{
  workMatrix.resize(sizeWorkMat);
  // Copy of Mass into workMatrix for LU-factorization.
  workMatrix[invMass].reset(new SimpleMatrix(*mass));
  // compute x[1] (and thus fExt if required)
  computeRhs(time);
  workMatrix[zeroMatrix].reset(new SimpleMatrix(ndof, ndof, ZERO));
  workMatrix[idMatrix].reset(new SimpleMatrix(ndof, ndof, IDENTITY));

  // jacobianXRhs
  if (K)
  {
    //  bloc10 of jacobianX is solution of Mass*Bloc10 = K
    workMatrix[jacobianXBloc10].reset(new SimpleMatrix(-1 * *K));
    workMatrix[invMass]->PLUForwardBackwardInPlace(*workMatrix[jacobianXBloc10]);
  }
  else
    workMatrix[jacobianXBloc10] = workMatrix[zeroMatrix] ;

  if (C)
  {
    //  bloc11 of jacobianX is solution of Mass*Bloc11 = C
    workMatrix[jacobianXBloc11].reset(new SimpleMatrix(-1 * *C));
    workMatrix[invMass]->PLUForwardBackwardInPlace(*workMatrix[jacobianXBloc11]);
  }
  else
    workMatrix[jacobianXBloc11] = workMatrix[zeroMatrix] ;

  jacobianXRhs.reset(new BlockMatrix(workMatrix[zeroMatrix], workMatrix[idMatrix], workMatrix[jacobianXBloc10], workMatrix[jacobianXBloc11]));
}

void LagrangianLinearTIDS::initialize(const string& simulationType, double time, unsigned int sizeOfMemory)
{
  // Memory allocation for p[0], p[1], p[2].
  initP(simulationType);

  // set q and q[1] to q0 and velocity0, initialize acceleration.
  *q[0] = *q0;
  *q[1] = *velocity0;

  // If z has not been set, we initialize it with a null vector of size 1, since z is required in plug-in functions call.
  if (! z)
    z.reset(new SimpleVector(1));

  // Set variables of top-class DynamicalSystem
  connectToDS(); // note that connection can not be done during constructor call, since user can complete the ds after (add plugin or anything else).

  checkDynamicalSystem();

  // Initialize memory vectors
  initMemory(sizeOfMemory);

  // rhs and its jacobian
  initRhs(time);
}

void LagrangianLinearTIDS::setK(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != ndof || newValue.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - setK: inconsistent input matrix size ");

  if (!K)
    K.reset(new SimpleMatrix(newValue));
  else
    *K = newValue;
}

void LagrangianLinearTIDS::setKPtr(SP::SiconosMatrix newPtr)
{
  if (newPtr->size(0) != ndof || newPtr->size(1) != ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - setKPtr: inconsistent input matrix size ");
  K = newPtr;
}


void LagrangianLinearTIDS::setCPtr(SP::SiconosMatrix newPtr)
{
  if (newPtr->size(0) != ndof || newPtr->size(1) != ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - setCPtr: inconsistent input matrix size ");

  C = newPtr;
}

void LagrangianLinearTIDS::display() const
{
  LagrangianDS::display();
  cout << "===== Lagrangian Linear Time Invariant System display ===== " << endl;
  cout << "- Stiffness Matrix K : " << endl;
  if (K) K->display();
  else cout << "-> NULL" << endl;
  cout << "- Viscosity Matrix C : " << endl;
  if (C) C->display();
  else cout << "-> NULL" << endl;
  cout << "=========================================================== " << endl;
}

void LagrangianLinearTIDS::computeRhs(double time, bool)
{
  // second argument is useless but present because of top-class function overloading.

  *q[2] = *p[2]; // Warning: p update is done in Interactions/Relations

  if (fExt)
  {
    computeFExt(time);
    *q[2] += *fExt; // This supposes that fExt is up to date!!
  }

  if (K)
    *q[2] -= prod(*K, *q[0]);

  if (C)
    *q[2] -= prod(*C, *q[1]);

  // Then we search for q[2], such as Mass*q[2] = fExt - Cq[1] - Kq[0] + p.
  workMatrix[invMass]->PLUForwardBackwardInPlace(*q[2]);
}

void LagrangianLinearTIDS::computeJacobianXRhs(double time, bool)
{
  // Nothing to be done since jacobianXRhs is constant and filled during initialize.
  // But this function is required, since it is called from Lsodar (if not present, the one of LagrangianDS will be called)
}

void LagrangianLinearTIDS::saveSpecificDataToXML()
{
  assert(dsxml &&
         "LagrangianLinearTIDS::saveDSToXML - object DynamicalSystemXML does not exist");

  /*  SP::LagrangianDSXML lgptr = boost::static_pointer_cast <LagrangianDSXML>(dsxml);
  lgptr->setMassMatrix( *mass );
  lgptr->setQ( *q[0] );
  lgptr->setQ0( *q0 );
  lgptr->setQMemory( *qMemory );
  lgptr->setVelocity( *q[1] );
  lgptr->setVelocity0( *velocity0 );
  lgptr->setVelocityMemory( *velocityMemory );

  // FExt
  if( lgptr->hasFExt() )
    {
      if( !lgptr->isFExtPlugin())
  {
    lgptr->setFExtVector( *fExt );
  }
    }
  else
    {
      lgptr->setFExtPlugin(fExt->getPluginName());
    }
  (boost::static_pointer_cast <LagrangianLinearTIDSXML>(dsxml))->setK( *K );
  (boost::static_pointer_cast <LagrangianLinearTIDSXML>(dsxml))->setC( *C );*/
}

LagrangianLinearTIDS* LagrangianLinearTIDS::convert(DynamicalSystem* ds)
{
  LagrangianLinearTIDS* ltids = dynamic_cast<LagrangianLinearTIDS*>(ds);
  return ltids;
}

