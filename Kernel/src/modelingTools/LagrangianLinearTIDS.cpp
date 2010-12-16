/* Siconos-Kernel, Copyright INRIA 2005-2010.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/
#include "LagrangianLinearTIDS.hpp"
#include "LagrangianLinearTIDSXML.hpp"
#include "BlockMatrix.hpp"

using namespace std;

// --- Default (private) constructor ---
LagrangianLinearTIDS::LagrangianLinearTIDS(): LagrangianDS()
{
}

// --- Constructor from an xml file---
LagrangianLinearTIDS::LagrangianLinearTIDS(SP::DynamicalSystemXML dsxml): LagrangianDS(dsxml)
{
  SP::LagrangianLinearTIDSXML lltidsxml = boost::static_pointer_cast <LagrangianLinearTIDSXML>(dsxml);

  // If FInt or NNL is given: ignored.
  if (lltidsxml->hasFInt() ||  lltidsxml->hasNNL())
    cout << "!!!!! Warning : LagrangianLinearTIDS: xml constructor, FInt or NNL input will be ignored in xml file." << endl;

  // K and C

  if (lltidsxml->hasK())
  {
    _K.reset(new SimpleMatrix(lltidsxml->getK()));
  }
  if (lltidsxml->hasC())
  {
    _C.reset(new SimpleMatrix(lltidsxml->getC()));
  }

  if (lltidsxml->hasFExt())
  {
    string plugin = lltidsxml->getFExtPlugin();
    setComputeFExtFunction(SSL::getPluginName(plugin), SSL::getPluginFunctionName(plugin));
    _fExt.reset(new SimpleVector(_ndof));
  }
}

// --- Constructor from a set of data - _Mass, K and C ---
LagrangianLinearTIDS::LagrangianLinearTIDS(SP::SimpleVector newQ0, SP::SimpleVector newVelocity0,
    SP::SiconosMatrix newMass,  SP::SiconosMatrix newK, SP::SiconosMatrix newC):
  LagrangianDS(newQ0, newVelocity0, newMass)
{
  assert((newK->size(0) == _ndof && newK->size(1) == _ndof) &&
         "LagrangianLinearTIDS - constructor from data, inconsistent size between K and _ndof");

  assert((newC->size(0) == _ndof && newC->size(1) == _ndof) &&
         "LagrangianLinearTIDS - constructor from data, inconsistent size between C and ndof");

  _K = newK;
  _C = newC;

}

// --- Constructor from a set of data - Mass, no K and no C ---
LagrangianLinearTIDS::LagrangianLinearTIDS(SP::SimpleVector newQ0, SP::SimpleVector newVelocity0,
    SP::SiconosMatrix newMass):
  LagrangianDS(newQ0, newVelocity0, newMass)
{
}
LagrangianLinearTIDS::LagrangianLinearTIDS(const SimpleVector& newQ0, const SimpleVector& newVelocity0, const SiconosMatrix& newMass):
  LagrangianDS(createSPtrSiconosVector((SimpleVector&)newQ0), createSPtrSiconosVector((SimpleVector&)newVelocity0), createSPtrSiconosMatrix((SimpleMatrix&)newMass))
{
}
LagrangianLinearTIDS::~LagrangianLinearTIDS()
{}

bool LagrangianLinearTIDS::checkDynamicalSystem()
{
  bool output = true;
  // _ndof
  if (_ndof == 0)
  {
    RuntimeException::selfThrow("LagrangianLinearTIDS::checkDynamicalSystem - number of degrees of freedom is equal to 0.");
    output = false;
  }

  // q0 and velocity0
  if (! _q0 || ! _velocity0)
  {
    RuntimeException::selfThrow("LagrangianLinearTIDS::checkDynamicalSystem - initial conditions are badly set.");
    output = false;
  }

  // Mass
  if (! _mass)
  {
    RuntimeException::selfThrow("LagrangianLinearTIDS::checkDynamicalSystem - Mass is not set.");
    output = false;
  }

  if (!output) cout << "LagrangianLinearTIDS Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
  return output;
}

void LagrangianLinearTIDS::initRhs(double time)
{
  _workMatrix.resize(sizeWorkMat);
  // Copy of Mass into _workMatrix for LU-factorization.
  _workMatrix[invMass].reset(new SimpleMatrix(*_mass));

  // compute x[1] (and thus _fExt if required)
  computeRhs(time);
  _workMatrix[zeroMatrix].reset(new SimpleMatrix(_ndof, _ndof, Siconos::ZERO));
  _workMatrix[idMatrix].reset(new SimpleMatrix(_ndof, _ndof, Siconos::IDENTITY));

  // jacobianRhsx
  if (_K)
  {
    //  bloc10 of jacobianX is solution of Mass*Bloc10 = K
    _workMatrix[jacobianXBloc10].reset(new SimpleMatrix(-1 * *_K));
    _workMatrix[invMass]->PLUForwardBackwardInPlace(*_workMatrix[jacobianXBloc10]);
  }
  else
    _workMatrix[jacobianXBloc10] = _workMatrix[zeroMatrix] ;

  if (_C)
  {
    //  bloc11 of jacobianX is solution of Mass*Bloc11 = C
    _workMatrix[jacobianXBloc11].reset(new SimpleMatrix(-1 * *_C));
    _workMatrix[invMass]->PLUForwardBackwardInPlace(*_workMatrix[jacobianXBloc11]);
  }
  else
    _workMatrix[jacobianXBloc11] = _workMatrix[zeroMatrix] ;

  _jacxRhs.reset(new BlockMatrix(_workMatrix[zeroMatrix], _workMatrix[idMatrix],
                                 _workMatrix[jacobianXBloc10], _workMatrix[jacobianXBloc11]));
}

void LagrangianLinearTIDS::initialize(const string& simulationType,
                                      double time, unsigned int sizeOfMemory)
{
  // Memory allocation for p[0], p[1], p[2].
  initP(simulationType);

  // set q and _q[1] to _q0 and velocity0, initialize acceleration.
  *_q[0] = *_q0;
  *_q[1] = *_velocity0;

  if (_boundaryConditions)
  {
    _reactionToBoundaryConditions.reset(new SimpleVector(_boundaryConditions->velocityIndices()->size()));
  }

  if (!_workFree)
    _workFree.reset(new SimpleVector(getDim()));

  // If z has not been set, we initialize it with a null vector of
  // size 1, since z is required in plug-in functions call.
  if (! _z)
    _z.reset(new SimpleVector(1));

  // Set variables of top-class DynamicalSystem
  connectToDS(); // note that connection can not be done during
  // constructor call, since user can complete the ds
  // after (add plugin or anything else).

  checkDynamicalSystem();

  // Initialize memory vectors
  initMemory(sizeOfMemory);

  // rhs and its jacobian
  //initRhs(time);
}

void LagrangianLinearTIDS::setK(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != _ndof || newValue.size(1) != _ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - setK: inconsistent input matrix size ");

  if (!_K)
    _K.reset(new SimpleMatrix(newValue));
  else
    *_K = newValue;
}

void LagrangianLinearTIDS::setKPtr(SP::SiconosMatrix newPtr)
{
  if (newPtr->size(0) != _ndof || newPtr->size(1) != _ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - setKPtr: inconsistent input matrix size ");
  _K = newPtr;
}


void LagrangianLinearTIDS::setCPtr(SP::SiconosMatrix newPtr)
{
  if (newPtr->size(0) != _ndof || newPtr->size(1) != _ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - setCPtr: inconsistent input matrix size ");

  _C = newPtr;
}

void LagrangianLinearTIDS::display() const
{
  LagrangianDS::display();
  cout << "===== Lagrangian Linear Time Invariant System display ===== " << endl;
  cout << "- Stiffness Matrix K : " << endl;
  if (_K) _K->display();
  else cout << "-> NULL" << endl;
  cout << "- Viscosity Matrix C : " << endl;
  if (_C) _C->display();
  else cout << "-> NULL" << endl;
  cout << "=========================================================== " << endl;
}

void LagrangianLinearTIDS::computeRhs(double time, bool)
{
  // second argument is useless but present because of top-class function overloading.
  *_q[2] = *_p[2]; // Warning: p update is done in Interactions/Relations
  //   std::cout << "LagrangianTIDS :: computeRhs " << std::endl ;
  //   std::cout << " p[2] " << std::endl ;
  //  _p[2]->display();
  if (_fExt)
  {
    computeFExt(time);
    *_q[2] += *_fExt; // This supposes that _fExt is up to date!!
  }
  if (_K)
    *_q[2] -= prod(*_K, *_q[0]);

  if (_C)
    *_q[2] -= prod(*_C, *_q[1]);

  // Then we search for _q[2], such as Mass*_q[2] = _fExt - C_q[1] - K_q[0] + p.
  _workMatrix[invMass]->PLUForwardBackwardInPlace(*_q[2]);
  _workFree->zero();
  if (_fExt)
  {
    computeFExt(time);
    *_workFree += *_fExt; // This supposes that _fExt is up to date!!
  }

  if (_K)
    *_workFree -= prod(*_K, *_q[0]);

  if (_C)
    *_workFree -= prod(*_C, *_q[1]);

  // Then we search for _workFree, such as Mass*_workfree = _fExt - C_q[1] - K_q[0] .
  _workMatrix[invMass]->PLUForwardBackwardInPlace(*_workFree);


  //   std::cout << "LagrangianTIDS :: computeRhs " << std::endl ;
  //   std::cout << " q[2] " << std::endl ;
  //   _q[2]->display();
  //   std::cout << " _workFree " << std::endl ;
  //   _workFree->display();

}

void LagrangianLinearTIDS::computeJacobianRhsx(double time, bool)
{
  // Nothing to be done since jacobianRhsx is constant and filled
  // during initialize.  But this function is required, since it is
  // called from Lsodar (if not present, the one of LagrangianDS will
  // be called)
}

void LagrangianLinearTIDS::saveSpecificDataToXML()
{
  assert(_dsxml &&
         "LagrangianLinearTIDS::saveDSToXML - object DynamicalSystemXML does not exist");

  SP::LagrangianDSXML lgptr = boost::static_pointer_cast <LagrangianDSXML>(_dsxml);
  lgptr->setMassMatrix(*_mass);
  lgptr->setQ(*_q[0]);
  lgptr->setQ0(*_q0);
  lgptr->setQMemory(*_qMemory);
  lgptr->setVelocity(*_q[1]);
  lgptr->setVelocity0(*_velocity0);
  lgptr->setVelocityMemory(*_velocityMemory);

  // FExt
  if (lgptr->hasFExt())
  {
    if (!lgptr->isFExtPlugin())
    {
      lgptr->setFExtVector(*_fExt);
    }
  }
  else
  {
    lgptr->setFExtPlugin(_pluginFExt->getPluginName());
  }
  (boost::static_pointer_cast <LagrangianLinearTIDSXML>(_dsxml))->setK(*_K);
  (boost::static_pointer_cast <LagrangianLinearTIDSXML>(_dsxml))->setC(*_C);
}

LagrangianLinearTIDS* LagrangianLinearTIDS::convert(DynamicalSystem* ds)
{
  LagrangianLinearTIDS* ltids = dynamic_cast<LagrangianLinearTIDS*>(ds);
  return ltids;
}

