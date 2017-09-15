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
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"

#include <iostream>

// --- Constructor from a initial conditions and matrix-operators
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

void LagrangianLinearTIDS::initRhs(double time)
{
  // _rhsMatrices.resize(numberOfRhsMatrices);
  // // Copy of Mass into _workMatrix for LU-factorization.
  // _inverseMass.reset(new SimpleMatrix(*_mass));

  // // compute x[1] (and thus _fExt if required)
  // computeRhs(time);

  LagrangianDS::initRhs(time);

  // jacobianRhsx
  if (_K)
  {
    //  bloc10 of jacobianX is solution of Mass*Bloc10 = K
    if(!_rhsMatrices[jacobianXBloc10])
      _rhsMatrices[jacobianXBloc10].reset(new SimpleMatrix(-1 * *_K));
    _inverseMass->PLUForwardBackwardInPlace(*_rhsMatrices[jacobianXBloc10]);
  }
  else
    _rhsMatrices[jacobianXBloc10] = _rhsMatrices[zeroMatrix] ;

  if (_C)
  {
    //  bloc11 of jacobianX is solution of Mass*Bloc11 = C
    if(!_rhsMatrices[jacobianXBloc11])
      _rhsMatrices[jacobianXBloc11].reset(new SimpleMatrix(-1 * *_C));
    _inverseMass->PLUForwardBackwardInPlace(*_rhsMatrices[jacobianXBloc11]);
  }
  else
    _rhsMatrices[jacobianXBloc11] = _rhsMatrices[zeroMatrix] ;

  if(_C || _K)
    _jacxRhs.reset(new BlockMatrix(_rhsMatrices[zeroMatrix], _rhsMatrices[idMatrix],
                                   _rhsMatrices[jacobianXBloc10], _rhsMatrices[jacobianXBloc11]));

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

void LagrangianLinearTIDS::setC(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != _ndof || newValue.size(1) != _ndof)
    RuntimeException::selfThrow("LagrangianLinearTIDS - setC: inconsistent input matrix size ");

  if (!_C)
    _C.reset(new SimpleMatrix(newValue));
  else
    *_C = newValue;
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

void LagrangianLinearTIDS::computeForces(double time, SP::SiconosVector q2, SP::SiconosVector v2)
{
  DEBUG_PRINT("LagrangianLinearTIDS::computeForces(double time, SP::SiconosVector q2, SP::SiconosVector v2) \n");

  if (!_forces)
  {
    _forces.reset(new SiconosVector(_ndof));
  }
  else
    _forces->zero();

  if (_fExt)
  {
    computeFExt(time);
    *_forces += *_fExt;
  }
  if (_K)
    *_forces -= prod(*_K, *q2);
  if (_C)
    *_forces -= prod(*_C, *v2);
}

