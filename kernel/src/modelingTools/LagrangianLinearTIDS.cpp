/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#include "LagrangianLinearTIDS.hpp"
#include "BlockMatrix.hpp"
#include "debug.h"

#include <iostream>

// --- Constructor from a set of data - _Mass, K and C ---
LagrangianLinearTIDS::LagrangianLinearTIDS(SP::SiconosVector newQ0, SP::SiconosVector newVelocity0,
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
LagrangianLinearTIDS::LagrangianLinearTIDS(SP::SiconosVector newQ0, SP::SiconosVector newVelocity0,
    SP::SiconosMatrix newMass):
  LagrangianDS(newQ0, newVelocity0, newMass)
{
}
//LagrangianLinearTIDS::LagrangianLinearTIDS(const SiconosVector& newQ0, const SiconosVector& newVelocity0, const SiconosMatrix& newMass):
//  LagrangianDS(createSPtrSiconosVector((SiconosVector&)newQ0), createSPtrSiconosVector((SiconosVector&)newVelocity0), createSPtrSiconosMatrix((SimpleMatrix&)newMass)){
//}
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

  if (!output) std::cout << "LagrangianLinearTIDS Warning: your dynamical system seems to be uncomplete (check = false)" <<std::endl;
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

void LagrangianLinearTIDS::initialize(double time, unsigned int sizeOfMemory)
{

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
  std::cout << "===== Lagrangian Linear Time Invariant System display ===== " <<std::endl;
  std::cout << "- Mass Matrix M : " <<std::endl;
  if (_mass) _mass->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "- Stiffness Matrix K : " <<std::endl;
  if (_K) _K->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "- Viscosity Matrix C : " <<std::endl;
  if (_C) _C->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "=========================================================== " <<std::endl;
}

void LagrangianLinearTIDS::computeRhs(double time, bool)
{

  computeForces(time, _q[0], _q[1]);

  // second argument is useless but present because of top-class function overloading.
  *_q[2] = *_p[2]; // Warning: p update is done in Interactions/Relations
  //    std::cout << "LagrangianTIDS :: computeRhs " << std::endl ;
  //    std::cout << " p[2] " << std::endl ;
  //  _p[2]->display();
  *_q[2] += *_forces;
  // Then we search for _q[2], such as Mass*_q[2] = _fExt - C_q[1] - K_q[0] + p.
  _workMatrix[invMass]->PLUForwardBackwardInPlace(*_q[2]);

  _workspace[free]->zero();
  computeForces(time, _q[0], _q[1]);
  *_workspace[free] = *_forces;
  // Then we search for _workspace[free], such as Mass*_workfree = _fExt - C_q[1] - K_q[0] .
  _workMatrix[invMass]->PLUForwardBackwardInPlace(*_workspace[free]);

  //    std::cout << "LagrangianTIDS :: computeRhs " << std::endl ;
  //    std::cout << " q[2] " << std::endl ;
  //   _q[2]->display();
  //    std::cout << " _workspace[free] " << std::endl ;
  //   _workspace[free]->display();

}
void LagrangianLinearTIDS::computeForces(double time)
{
  computeForces(time, _q[0], _q[1]);
}

void LagrangianLinearTIDS::computeForces(double time, SP::SiconosVector q2, SP::SiconosVector v2)
{

  DEBUG_PRINT("LagrangianLinearTIDS::computeForces(double time, SP::SiconosVector q2, SP::SiconosVector v2) \n");
  // Warning: an operator (fInt ...) may be set (ie allocated and not NULL) but not plugged, that's why two steps are required here.
  if (!_forces)
  {
    _forces.reset(new SiconosVector(_ndof));
  }
  // 1 - Computes the required forces
  //_forces->zero();

  if (_fExt)
  {
    computeFExt(time);
    *_forces = *_fExt; // This supposes that _fExt is up to date!!
  }
  if (_K)
    *_forces -= prod(*_K, *q2);
  if (_C)
    *_forces -= prod(*_C, *v2);
}

void LagrangianLinearTIDS::computeJacobianRhsx(double time, bool)
{
  // Nothing to be done since jacobianRhsx is constant and filled
  // during initialize.  But this function is required, since it is
  // called from LsodarOSI (if not present, the one of LagrangianDS will
  // be called)
}
