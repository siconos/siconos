/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#include "CableDS.hpp"

#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"

siconos::mechanics::fem::CableDS::CableDS(
    std::shared_ptr<siconos::algebra::SiconosVector> q0, std::shared_ptr<siconos::algebra::SiconosVector> velocity0,
    std::shared_ptr<siconos::algebra::SiconosMatrix> mass, ExternalForcesFunction fext)
    : LagrangianDS(q0, velocity0, mass), computefext_{fext}
{
  // Constructor with initial state and mass.
  // We assume that q0, v0 and mass are computed by the "cable model" based on mesh and other
  // inputs Maybe these things are to be computed inside the DS? To be discussed ...
  // The call to LagrangianDS base constructor  ensures a proper allocation of memories for q0,
  // v0 Mass is just a pointer link. Mass alloc : to be done in cable model if mass is a shared
  // pointer input. _mass = std::make_shared<siconos::algebra::SimpleMatrix>(_ndof, _ndof,
  // UBLAS_TYPE::SPARSE);
  // We can deal with variable mass later.
  _hasConstantMass = true;

  // What may happen in cable model: (see examples in testCableDS)
  // Case 1: no fext
  // cableDS{q0, v0, mass}

  // Case 2: constant fext
  // cableDS{q0, v0, mass}
  // call to setFExtPtr()
  // cable

  if (fext) {
    _hasConstantFExt = false; // Indeed, this is the default for SecondOrderDS
    _fExt = std::make_shared<siconos::algebra::SiconosVector>(_ndof);
  }
  else
    _hasConstantFExt = true;
  // In that case _fExt = nullptr
  // setFextPtr is to be called later by cable model to set a constant fext.

  // _ndof is given by the size of q0 during SecondOrderDS build
  _forces = std::make_shared<siconos::algebra::SiconosVector>(_ndof);

  // We will use _jacobianqForces and _jacobianvForces to save tangent stiffness and damping
  // matrices.
  // Those are attributes of LagrangianDS class.

  _jacobianqForces =
      std::make_shared<siconos::algebra::SimpleMatrix>(_ndof, _ndof, siconos::algebra::UblasType::SPARSE);

  _jacobianqDotForces =
      std::make_shared<siconos::algebra::SimpleMatrix>(_ndof, _ndof, siconos::algebra::UblasType::SPARSE);
}

void siconos::mechanics::fem::CableDS::computeFExt(double time)
{
  assert(_fExt);
  assert(computefext_);
  // Call the std::function attribute that must be connected to some external function
  computefext_(time, _fExt);
  // ...
}

void siconos::mechanics::fem::CableDS::computeForces(double time,
                                                     std::shared_ptr<siconos::algebra::SiconosVector> q,
                                                     std::shared_ptr<siconos::algebra::SiconosVector> velocity)
{
  assert(_forces);
  // if (!_forces) {
  //   _forces = std::make_shared<siconos::algebra::SiconosVector>(_ndof);
  // } // --> done during constructor call.
  // else
  _forces->zero();

  // Update internal and external forces
  // Here you must:
  //  - compute internal forces and save them somewhere or directly in _forces
  //  - assume that jacobian are uptodate or maybe recompute them
  // *_forces = ...
  //  - compute external forces in fExt (or just get them if they are constant)
  if (!_hasConstantFExt)
    computeFExt(time);
  if (_fExt)
    *_forces += *_fExt;
}

// \f$ \nabla_q F \f$
void siconos::mechanics::fem::CableDS::computeJacobianqForces(double time)
{
  // Call a local routine which compute tangent stiffness and update a local operator

  tangentStiffnessMatrix();
}

// \f$ \nabla_v F \f$
void siconos::mechanics::fem::CableDS::computeJacobianvForces(double time)
{

  // Call a local routine to compute damping matrix and update a local operator
  dampingMatrix();
}

// Probably not needed since mass will be constant. Called  by the integrator at each time step
// to override mass operator.
void siconos::mechanics::fem::CableDS::computeMass(std::shared_ptr<siconos::algebra::SiconosVector> position)
{

  // _mass = ...
}

void siconos::mechanics::fem::CableDS::tangentStiffnessMatrix(/** ...*/)
{
  // must update jacobianqForces
}

void siconos::mechanics::fem::CableDS::dampingMatrix(/** ...*/)
{

  // must update jacobianqDotForces
}
