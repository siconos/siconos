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
#include "LagrangianDS.hpp"

#include <iostream>

#include "BlockMatrix.hpp"
#include "BlockVector.hpp"
#include "PluggedObject.hpp"  // for getPluginfunctionname ...
#include "SiconosConst.hpp"
#include "SiconosMatrixVectorOp.hpp"  // for matrix-vector prod
#include "SiconosVector.hpp"
#include "SiconosVectorOp.hpp"  // for inner_prod
#include "SiconosVisitor.hpp"
#include "SimpleMatrix.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES

#include "siconos_debug.h"

// Build from initial state only
siconos::modeling::LagrangianDS::LagrangianDS(
    std::shared_ptr<siconos::algebra::SiconosVector> q0,
    std::shared_ptr<siconos::algebra::SiconosVector> v0)
    : SecondOrderDS(2 * q0->size(), v0->size()), _hasConstantFExt(true) {
  assert(_ndof > 0 && "lagrangian dynamical system dimension should be greater than 0.");

  // Set initial conditions
  _q0 = q0;
  _velocity0 = v0;

  // -- Memory allocation for vector and matrix members --
  _q[0] = std::make_shared<siconos::algebra::SiconosVector>(*_q0);
  _q[1] = std::make_shared<siconos::algebra::SiconosVector>(*_velocity0);

  /** \todo lazy Memory allocation */
  _p[1] = std::make_shared<siconos::algebra::SiconosVector>(_ndof);

  _zeroPlugin();
}

// From initial state and constant mass matrix, \f$ M\ddot q = p \f$
siconos::modeling::LagrangianDS::LagrangianDS(
    std::shared_ptr<siconos::algebra::SiconosVector> q0,
    std::shared_ptr<siconos::algebra::SiconosVector> v0,
    std::shared_ptr<siconos::algebra::SiconosMatrix> newMass)
    : LagrangianDS(q0, v0) {
  _mass = newMass;
}

void siconos::modeling::LagrangianDS::allocateMass() {
  if (!_mass) {
    _mass = std::make_shared<siconos::algebra::SimpleMatrix>(_ndof, _ndof);
  }
}

// From a set of data - Mass loaded from a plugin
// This constructor leads to the minimum Lagrangian System form: \f$ M(q)\ddot q = p \f$
siconos::modeling::LagrangianDS::LagrangianDS(
    std::shared_ptr<siconos::algebra::SiconosVector> q0,
    std::shared_ptr<siconos::algebra::SiconosVector> v0, const std::string& massName)
    : LagrangianDS(q0, v0) {
  _hasConstantMass = false;
  // Mass
  allocateMass();
  setComputeMassFunction(siconos::plugins::getPluginName(massName),
                         siconos::plugins::getPluginFunctionName(massName));
}

void siconos::modeling::LagrangianDS::_zeroPlugin() {
  _pluginMass = std::make_shared<siconos::plugins::PluggedObject>();
  _pluginFInt = std::make_shared<siconos::plugins::PluggedObject>();
  _pluginFExt = std::make_shared<siconos::plugins::PluggedObject>();
  _pluginFGyr = std::make_shared<siconos::plugins::PluggedObject>();
  _pluginJacqFInt = std::make_shared<siconos::plugins::PluggedObject>();
  _pluginJacqDotFInt = std::make_shared<siconos::plugins::PluggedObject>();
  _pluginJacqFGyr = std::make_shared<siconos::plugins::PluggedObject>();
  _pluginJacqDotFGyr = std::make_shared<siconos::plugins::PluggedObject>();
}

void siconos::modeling::LagrangianDS::initializeNonSmoothInput(unsigned int level) {
  if (!_p[level]) _p[level] = std::make_shared<siconos::algebra::SiconosVector>(_ndof);
}

void siconos::modeling::LagrangianDS::resetToInitialState() {
  if (_q0) {
    *(_q[0]) = *_q0;
  } else
    THROW_EXCEPTION(
        "siconos::modeling::LagrangianDS::resetToInitialState - initial position _q0 is null");
  if (_velocity0) {
    *(_q[1]) = *_velocity0;
  } else
    THROW_EXCEPTION(
        "siconos::modeling::LagrangianDS::resetToInitialState - initial velocity _velocity0 "
        "is null");
}

void siconos::modeling::LagrangianDS::init_generalized_coordinates(unsigned int level) {
  assert(level > 1);
  if (!_q[level]) _q[level] = std::make_shared<siconos::algebra::SiconosVector>(_ndof);
}

void siconos::modeling::LagrangianDS::init_inverse_mass() {
  if (_mass && !_inverseMass) {
    computeMass();
    _inverseMass = std::make_shared<siconos::algebra::SimpleMatrix>(*_mass);
  }
}

void siconos::modeling::LagrangianDS::update_inverse_mass() {
  if (_mass && _inverseMass && !_hasConstantMass) {
    computeMass();
    *_inverseMass = *_mass;
  }
}

void siconos::modeling::LagrangianDS::init_forces() {
  // Allocate memory for forces and its jacobians.void
  // siconos::modeling::LagrangianDS::init_forces() Needed only for integrators with
  // first-order formulation.
  if (_fInt || _fExt || _fGyr) {
    if (!_forces) _forces = std::make_shared<siconos::algebra::SiconosVector>(_ndof);
  }

  if (_fInt || _fGyr) {
    if (!_jacobianqForces)
      _jacobianqForces = std::make_shared<siconos::algebra::SimpleMatrix>(_ndof, _ndof);
    if (!_jacobianqDotForces)
      _jacobianqDotForces = std::make_shared<siconos::algebra::SimpleMatrix>(_ndof, _ndof);
  }
}

void siconos::modeling::LagrangianDS::initRhs(double time) {
  DEBUG_BEGIN("siconos::modeling::LagrangianDS::initRhs(double time)\n");
  // dim
  _n = 2 * _ndof;

  // All links between DS and LagrangianDS class members are pointer links, which means
  // that no useless memory is allocated when connection is established.
  // One exception: zero and identity matrices, used to filled in M and jacobianfx.

  // Initial conditions and state

  // WARNING : this function is supposed to be called
  // by the OneStepIntegrator, and maybe several times for the same DS
  // if the system is involved in more than one interaction. So, we must check
  // if p2 and q2 already exist to be sure that DSlink won't be lost.

  _x0 = std::make_shared<siconos::algebra::SiconosVector>(*_q0, *_velocity0);

  _x[0] = std::make_shared<siconos::algebra::SiconosVector>(*_q[0], *_q[1]);

  if (!_q[2]) _q[2] = std::make_shared<siconos::algebra::SiconosVector>(_ndof);

  _x[1] = std::make_shared<siconos::algebra::SiconosVector>(*_q[1], *_q[2]);

  // Everything concerning rhs and its jacobian is handled in initRhs and computeXXX related
  // functions.

  _rhsMatrices.resize(numberOfRhsMatrices_);

  if (!_p[2]) _p[2] = std::make_shared<siconos::algebra::SiconosVector>(_ndof);

  init_forces();
  init_inverse_mass();

  computeRhs(time);

  bool flag1 = false, flag2 = false;
  if (_jacobianqForces) {
    // Solve MjacobianX(1,0) = jacobianFL[0]
    computeJacobianqForces(time);

    _rhsMatrices[jacobianXBloc10_] =
        std::make_shared<siconos::algebra::SimpleMatrix>(*_jacobianqForces);
    _inverseMass->Solve(*_rhsMatrices[jacobianXBloc10_]);
    flag1 = true;
  }

  if (_jacobianqDotForces) {
    // Solve MjacobianX(1,1) = jacobianFL[1]
    computeJacobianqDotForces(time);
    _rhsMatrices[jacobianXBloc11_] =
        std::make_shared<siconos::algebra::SimpleMatrix>(*_jacobianqDotForces);
    _inverseMass->Solve(*_rhsMatrices[jacobianXBloc11_]);
    flag2 = true;
  }

  if (!_rhsMatrices[zeroMatrix_])
    _rhsMatrices[zeroMatrix_] = std::make_shared<siconos::algebra::SimpleMatrix>(
        _ndof, _ndof, siconos::algebra::UblasType::ZERO);
  if (!_rhsMatrices[idMatrix_])
    _rhsMatrices[idMatrix_] = std::make_shared<siconos::algebra::SimpleMatrix>(
        _ndof, _ndof, siconos::algebra::UblasType::IDENTITY);

  if (flag1 && flag2)
    _jacxRhs = std::make_shared<siconos::algebra::BlockMatrix>(
        _rhsMatrices[zeroMatrix_], _rhsMatrices[idMatrix_], _rhsMatrices[jacobianXBloc10_],
        _rhsMatrices[jacobianXBloc11_]);
  else if (flag1)  // flag2 = false
    _jacxRhs = std::make_shared<siconos::algebra::BlockMatrix>(
        _rhsMatrices[zeroMatrix_], _rhsMatrices[idMatrix_], _rhsMatrices[jacobianXBloc10_],
        _rhsMatrices[zeroMatrix_]);
  else if (flag2)  // flag1 = false
    _jacxRhs = std::make_shared<siconos::algebra::BlockMatrix>(
        _rhsMatrices[zeroMatrix_], _rhsMatrices[idMatrix_], _rhsMatrices[zeroMatrix_],
        _rhsMatrices[jacobianXBloc11_]);
  else
    _jacxRhs = std::make_shared<siconos::algebra::BlockMatrix>(
        _rhsMatrices[zeroMatrix_], _rhsMatrices[idMatrix_], _rhsMatrices[zeroMatrix_],
        _rhsMatrices[zeroMatrix_]);
  DEBUG_EXPR(display(););
  DEBUG_END("siconos::modeling::LagrangianDS::initRhs(double time)\n");
}

// --- GETTERS/SETTERS ---

void siconos::modeling::LagrangianDS::setQ(const siconos::algebra::SiconosVector& newValue) {
  if (newValue.size() != _ndof)
    THROW_EXCEPTION("LagrangianDS - setQ: inconsistent input vector size ");

  if (!_q[0])
    _q[0] = std::make_shared<siconos::algebra::SiconosVector>(newValue);
  else
    *_q[0] = newValue;
}

void siconos::modeling::LagrangianDS::setQPtr(
    std::shared_ptr<siconos::algebra::SiconosVector> newPtr) {
  if (newPtr->size() != _ndof)
    THROW_EXCEPTION("LagrangianDS - setQPtr: inconsistent input vector size ");
  _q[0] = newPtr;
}

void siconos::modeling::LagrangianDS::setQ0(const siconos::algebra::SiconosVector& newValue) {
  if (newValue.size() != _ndof)
    THROW_EXCEPTION("LagrangianDS - setQ0: inconsistent input vector size ");

  if (!_q0)
    _q0 = std::make_shared<siconos::algebra::SiconosVector>(newValue);
  else
    *_q0 = newValue;
}

void siconos::modeling::LagrangianDS::setQ0Ptr(
    std::shared_ptr<siconos::algebra::SiconosVector> newPtr) {
  if (newPtr->size() != _ndof)
    THROW_EXCEPTION("LagrangianDS - setQ0Ptr: inconsistent input vector size ");
  _q0 = newPtr;
}

void siconos::modeling::LagrangianDS::setVelocity0(
    const siconos::algebra::SiconosVector& newValue) {
  if (newValue.size() != _ndof)
    THROW_EXCEPTION("LagrangianDS - setVelocity0: inconsistent input vector size ");

  if (!_velocity0)
    _velocity0 = std::make_shared<siconos::algebra::SiconosVector>(newValue);
  else
    *_velocity0 = newValue;
}

void siconos::modeling::LagrangianDS::setVelocity(
    const siconos::algebra::SiconosVector& newValue) {
  if (newValue.size() != _ndof)
    THROW_EXCEPTION("LagrangianDS - setVelocity: inconsistent input vector size ");

  if (!_q[1])
    _q[1] = std::make_shared<siconos::algebra::SiconosVector>(newValue);
  else
    *_q[1] = newValue;
}

void siconos::modeling::LagrangianDS::setVelocityPtr(
    std::shared_ptr<siconos::algebra::SiconosVector> newPtr) {
  if (newPtr->size() != _ndof)
    THROW_EXCEPTION("LagrangianDS - setVelocityPtr: inconsistent input vector size ");
  _q[1] = newPtr;
}

void siconos::modeling::LagrangianDS::setVelocity0Ptr(
    std::shared_ptr<siconos::algebra::SiconosVector> newPtr) {
  if (newPtr->size() != _ndof)
    THROW_EXCEPTION("LagrangianDS - setVelocity0Ptr: inconsistent input vector size ");
  _velocity0 = newPtr;
}

void siconos::modeling::LagrangianDS::computeMass() {
  DEBUG_BEGIN("siconos::modeling::LagrangianDS::computeMass()\n");
  DEBUG_EXPR(_q[0]->display());
  computeMass(_q[0]);
  DEBUG_EXPR(_mass->display());
  DEBUG_END("siconos::modeling::LagrangianDS::computeMass()\n");
}

void siconos::modeling::LagrangianDS::computeMass(
    std::shared_ptr<siconos::algebra::SiconosVector> position) {
  if (_mass && !_hasConstantMass && _pluginMass->fPtr) {
    ((siconos::plugins::FPtr7)_pluginMass->fPtr)(_ndof, &(*position)(0), &(*_mass)(0, 0),
                                                 _z->size(), &(*_z)(0));
    _mass->resetFactorizationFlags();
  }
}

void siconos::modeling::LagrangianDS::computeFInt(double time) {
  if (_fInt && _pluginFInt->fPtr)
    ((siconos::plugins::FPtr6)_pluginFInt->fPtr)(time, _ndof, &(*_q[0])(0), &(*_q[1])(0),
                                                 &(*_fInt)(0), _z->size(), &(*_z)(0));
}
void siconos::modeling::LagrangianDS::computeFInt(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> position,
    std::shared_ptr<siconos::algebra::SiconosVector> velocity) {
  if (_fInt && _pluginFInt->fPtr)
    ((siconos::plugins::FPtr6)_pluginFInt->fPtr)(time, _ndof, &(*position)(0), &(*velocity)(0),
                                                 &(*_fInt)(0), _z->size(), &(*_z)(0));
}

void siconos::modeling::LagrangianDS::computeFExt(double time) {
  if (!_hasConstantFExt) {
    if (_fExt && _pluginFExt->fPtr)
      ((siconos::plugins::VectorFunctionOfTime)_pluginFExt->fPtr)(time, _ndof, &(*_fExt)(0),
                                                                  _z->size(), &(*_z)(0));
  }
}
void siconos::modeling::LagrangianDS::computeFGyr() {
  if (_fGyr && _pluginFGyr->fPtr)
    ((siconos::plugins::FPtr5)_pluginFGyr->fPtr)(_ndof, &(*_q[0])(0), &(*_q[1])(0),
                                                 &(*_fGyr)(0), _z->size(), &(*_z)(0));
}

void siconos::modeling::LagrangianDS::computeFGyr(
    std::shared_ptr<siconos::algebra::SiconosVector> position,
    std::shared_ptr<siconos::algebra::SiconosVector> velocity) {
  if (_fGyr && _pluginFGyr->fPtr)
    ((siconos::plugins::FPtr5)_pluginFGyr->fPtr)(_ndof, &(*position)(0), &(*velocity)(0),
                                                 &(*_fGyr)(0), _z->size(), &(*_z)(0));
}

void siconos::modeling::LagrangianDS::computeJacobianFIntq(double time) {
  DEBUG_BEGIN("siconos::modeling::LagrangianDS::computeJacobianFIntq()\n");
  DEBUG_EXPR(_q[0]->display());
  DEBUG_EXPR(_q[1]->display());
  if (_jacobianFIntq && _pluginJacqFInt->fPtr)
    ((siconos::plugins::FPtr6)_pluginJacqFInt->fPtr)(time, _ndof, &(*_q[0])(0), &(*_q[1])(0),
                                                     &(*_jacobianFIntq)(0, 0), _z->size(),
                                                     &(*_z)(0));
  DEBUG_EXPR(if (_jacobianFIntq) _jacobianFIntq->display(););
  DEBUG_END("siconos::modeling::LagrangianDS::computeJacobianFIntq()\n");
}
void siconos::modeling::LagrangianDS::computeJacobianFIntqDot(double time) {
  DEBUG_BEGIN("siconos::modeling::LagrangianDS::computeJacobianFIntqDot()\n");
  DEBUG_EXPR(_q[0]->display());
  DEBUG_EXPR(_q[1]->display());
  DEBUG_EXPR(_z->display());
  if (_jacobianFIntqDot && _pluginJacqDotFInt->fPtr)
    ((siconos::plugins::FPtr6)_pluginJacqDotFInt->fPtr)(
        time, _ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_jacobianFIntqDot)(0, 0), _z->size(),
        &(*_z)(0));
  DEBUG_EXPR(if (_jacobianFIntqDot) _jacobianFIntqDot->display(););
  DEBUG_END("siconos::modeling::LagrangianDS::computeJacobianFIntqDot()\n");
}

void siconos::modeling::LagrangianDS::computeJacobianFIntq(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> position,
    std::shared_ptr<siconos::algebra::SiconosVector> velocity) {
  DEBUG_BEGIN("siconos::modeling::LagrangianDS::computeJacobianFIntq()\n");
  DEBUG_EXPR(position->display());
  DEBUG_EXPR(velocity->display());
  if (_jacobianFIntq && _pluginJacqFInt->fPtr)
    ((siconos::plugins::FPtr6)_pluginJacqFInt->fPtr)(time, _ndof, &(*position)(0),
                                                     &(*velocity)(0), &(*_jacobianFIntq)(0, 0),
                                                     _z->size(), &(*_z)(0));
  DEBUG_EXPR(if (_jacobianFIntq) _jacobianFIntq->display(););
  DEBUG_END("siconos::modeling::LagrangianDS::computeJacobianFIntq()\n");
}
void siconos::modeling::LagrangianDS::computeJacobianFIntqDot(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> position,
    std::shared_ptr<siconos::algebra::SiconosVector> velocity) {
  if (_jacobianFIntqDot && _pluginJacqDotFInt->fPtr)
    ((siconos::plugins::FPtr6)_pluginJacqDotFInt->fPtr)(
        time, _ndof, &(*position)(0), &(*velocity)(0), &(*_jacobianFIntqDot)(0, 0), _z->size(),
        &(*_z)(0));
}

void siconos::modeling::LagrangianDS::computeJacobianFGyrq() {
  if (_pluginJacqFGyr->fPtr)
    ((siconos::plugins::FPtr5)_pluginJacqFGyr->fPtr)(
        _ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_jacobianFGyrq)(0, 0), _z->size(), &(*_z)(0));
}
void siconos::modeling::LagrangianDS::computeJacobianFGyrqDot() {
  if (_jacobianFGyrqDot && _pluginJacqDotFGyr->fPtr)
    ((siconos::plugins::FPtr5)_pluginJacqDotFGyr->fPtr)(
        _ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_jacobianFGyrqDot)(0, 0), _z->size(), &(*_z)(0));
}

void siconos::modeling::LagrangianDS::computeJacobianFGyrq(
    std::shared_ptr<siconos::algebra::SiconosVector> position,
    std::shared_ptr<siconos::algebra::SiconosVector> velocity) {
  if (_jacobianFGyrq && _pluginJacqFGyr->fPtr)
    ((siconos::plugins::FPtr5)_pluginJacqFGyr->fPtr)(_ndof, &(*position)(0), &(*velocity)(0),
                                                     &(*_jacobianFGyrq)(0, 0), _z->size(),
                                                     &(*_z)(0));
}

void siconos::modeling::LagrangianDS::computeJacobianFGyrqDot(
    std::shared_ptr<siconos::algebra::SiconosVector> position,
    std::shared_ptr<siconos::algebra::SiconosVector> velocity) {
  if (_jacobianFGyrqDot && _pluginJacqDotFGyr->fPtr)
    ((siconos::plugins::FPtr5)_pluginJacqDotFGyr->fPtr)(
        _ndof, &(*position)(0), &(*velocity)(0), &(*_jacobianFGyrqDot)(0, 0), _z->size(),
        &(*_z)(0));
}

void siconos::modeling::LagrangianDS::computeRhs(double time) {
  DEBUG_BEGIN("siconos::modeling::LagrangianDS::computeRhs(double time)");
  *_q[2] = *(_p[2]);  // Warning: r/p update is done in Interactions/Relations

  // if(_forces)
  //   {
  computeForces(time, _q[0], _q[1]);
  *_q[2] += *_forces;
  DEBUG_EXPR(_forces->display(););
  // #  }

  // Computes q[2] = inv(mass)*(fL+p) by solving Mq[2]=fL+p.
  // -- Case 1: if mass is constant, then a copy of imass is LU-factorized during
  // initialization and saved into _inverseMass
  // -- Case 2: mass is not constant, we copy it into _inverseMass
  // Then we proceed with PLUForwardBackward.
  // mass and inv(mass) computation
  if (_mass && !_hasConstantMass)  // if it is necessary to re-compute mass, FInt ..., ie if
                                   // they have not been compiled during the present time step
  {
    computeMass();
    *_inverseMass = *_mass;
  }

  //  if(mass->isPlugged()) : mass may be not plugged in LagrangianDS children
  if (_inverseMass) _inverseMass->Solve(*_q[2]);

  _x[1]->setBlock(0, *_q[1]);
  _x[1]->setBlock(_ndof, *_q[2]);
  DEBUG_END("siconos::modeling::LagrangianDS::computeRhs(double time)");
}

void siconos::modeling::LagrangianDS::computeJacobianRhsx(double time) {
  if (!_hasConstantMass) computeMass();

  //  if(mass->isPlugged()) : mass may b not plugged in LagrangianDS children

  if (_jacobianqForces || _jacobianqDotForces) {
    if (!_hasConstantMass)  // else inverseMass is already uptodate
      *_inverseMass = *_mass;
  }

  if (_jacobianqForces) {
    /** \warning the Jacobian of the inverse of the mass matrix
     * w.r.t q is not taken into account */

    std::shared_ptr<siconos::algebra::SiconosMatrix> bloc10 = _jacxRhs->block(1, 0);
    computeJacobianqForces(time);
    *bloc10 = *_jacobianqForces;
    _inverseMass->Solve(*bloc10);
  }

  if (_jacobianqDotForces) {
    std::shared_ptr<siconos::algebra::SiconosMatrix> bloc11 = _jacxRhs->block(1, 1);
    computeJacobianqDotForces(time);
    *bloc11 = *_jacobianqDotForces;
    _inverseMass->Solve(*bloc11);
  }
}

void siconos::modeling::LagrangianDS::computeForces(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> position,
    std::shared_ptr<siconos::algebra::SiconosVector> velocity) {
  if (!_forces) {
    _forces = std::make_shared<siconos::algebra::SiconosVector>(_ndof);
  } else
    _forces->zero();

  // 1 - Computes the required function
  computeFInt(time, position, velocity);
  computeFExt(time);
  computeFGyr(position, velocity);

  // seems ok.
  // if (_forces.use_count() == 1)
  // {
  //   //if not that means that fL is already (pointer-)connected with
  //   // either fInt, FGyr OR fExt.
  //_forces->zero();

  if (_fInt) *_forces -= *_fInt;

  if (_fExt) *_forces += *_fExt;

  if (_fGyr) *_forces -= *_fGyr;
  // }
}

void siconos::modeling::LagrangianDS::computeJacobianqForces(double time) {
  if (_jacobianqForces) {
    computeJacobianFIntq(time);
    computeJacobianFGyrq();

    // not true!
    // if( jacobianFL[i].use_count() == 1 )
    {
      // if not that means that jacobianFL[i] is already (pointer-)connected with
      //  either jacobianFInt or jacobianFGyr
      _jacobianqForces->zero();
      if (_jacobianFIntq) *_jacobianqForces -= *_jacobianFIntq;
      if (_jacobianFGyrq) *_jacobianqForces -= *_jacobianFGyrq;
    }
  }
  // else nothing.
}
void siconos::modeling::LagrangianDS::computeJacobianvForces(double time) {
  if (_jacobianqDotForces) {
    computeJacobianFIntqDot(time);
    computeJacobianFGyrqDot();

    // not true!
    // if( jacobianFL[i].use_count() == 1 )
    {
      // if not that means that jacobianFL[i] is already (pointer-)connected with
      //  either jacobianFInt or jacobianFGyr
      _jacobianqDotForces->zero();
      if (_jacobianFIntqDot) *_jacobianqDotForces -= *_jacobianFIntqDot;
      if (_jacobianFGyrqDot) *_jacobianqDotForces -= *_jacobianFGyrqDot;
    }
  }
  // else nothing.
}
// void siconos::modeling::LagrangianDS::computeJacobianZFL( double time){
//    THROW_EXCEPTION("siconos::modeling::LagrangianDS::computeJacobianZFL - not implemented");
// }

void siconos::modeling::LagrangianDS::display(bool brief) const {
  std::cout << "=====> Lagrangian System display (number: " << _number << ")." << std::endl;
  std::cout << "- _ndof : " << _ndof << std::endl;
  std::cout << "- q " << std::endl;
  if (_q[0])
    _q[0]->display();
  else
    std::cout << "-> nullptr" << std::endl;
  std::cout << "- q0 " << std::endl;
  if (_q0) _q0->display();
  std::cout << "- velocity " << std::endl;
  if (_q[1])
    _q[1]->display();
  else
    std::cout << "-> nullptr" << std::endl;
  std::cout << "- acceleration " << std::endl;
  if (_q[2])
    _q[2]->display();
  else
    std::cout << "-> nullptr" << std::endl;
  std::cout << "- v0 " << std::endl;
  if (_velocity0)
    _velocity0->display();
  else
    std::cout << "-> nullptr" << std::endl;
  std::cout << "- p[0] " << std::endl;
  if (_p[0])
    _p[0]->display();
  else
    std::cout << "-> nullptr" << std::endl;
  std::cout << "- p[1] " << std::endl;
  if (_p[1])
    _p[1]->display();
  else
    std::cout << "-> nullptr" << std::endl;
  std::cout << "- p[2] " << std::endl;
  if (_p[2])
    _p[2]->display();
  else
    std::cout << "-> nullptr" << std::endl;

  if (!brief) {
    std::cout << "- Mass " << std::endl;
    if (_mass)
      _mass->display();
    else
      std::cout << "-> nullptr" << std::endl;

    std::cout << "- Forces " << std::endl;
    if (_forces)
      _forces->display();
    else
      std::cout << "-> nullptr" << std::endl;
    std::cout << "- FInt " << std::endl;
    if (_fInt)
      _fInt->display();
    else
      std::cout << "-> nullptr" << std::endl;

    std::cout << "- jacobianqForces " << std::endl;
    if (_jacobianqForces)
      _jacobianqForces->display();
    else
      std::cout << "-> nullptr" << std::endl;
    std::cout << "- jacobianFIntq " << std::endl;
    if (_jacobianFIntq)
      _jacobianFIntq->display();
    else
      std::cout << "-> nullptr" << std::endl;

    std::cout << "- jacobianqDotForces " << std::endl;
    if (_jacobianqDotForces)
      _jacobianqDotForces->display();
    else
      std::cout << "-> nullptr" << std::endl;
  }

  std::cout << "===================================== " << std::endl;
}

// --- Functions for memory handling ---
void siconos::modeling::LagrangianDS::initMemory(unsigned int steps) {
  DEBUG_PRINTF(
      "siconos::modeling::LagrangianDS::initMemory(unsigned int steps) with steps = %i\n",
      steps);
  if (steps == 0)
    std::cout << "Warning : LagragianDS::initMemory with size equal to zero" << std::endl;
  else {
    _qMemory.setMemorySize(steps, _ndof);
    _velocityMemory.setMemorySize(steps, _ndof);
    _forcesMemory.setMemorySize(steps, _ndof);
    _pMemory.resize(3);

    // TODO : initMemory in graph + f(OSI/level)
    for (unsigned int level = 0; level < 3; ++level) {
      if (_pMemory[level].size() == 0) _pMemory[level].setMemorySize(steps, _ndof);
    }

    // swapInMemory();
  }
}

void siconos::modeling::LagrangianDS::swapInMemory() {
  _qMemory.swap(*_q[0]);
  _velocityMemory.swap(*_q[1]);
  if (_forces) _forcesMemory.swap(*_forces);

  // initialization of the reaction force due to the non smooth law
  // note: these are a no-op if either memory or vector is null
  _pMemory[0].swap(_p[0]);
  _pMemory[1].swap(_p[1]);
  _pMemory[2].swap(_p[2]);
  _xMemory.swap(_x[0]);
}

void siconos::modeling::LagrangianDS::resetAllNonSmoothParts() {
  if (_p[0]) _p[0]->zero();
  if (_p[1]) _p[1]->zero();
  if (_p[2]) _p[2]->zero();
}

void siconos::modeling::LagrangianDS::resetNonSmoothPart(unsigned int level) {
  if (level < siconos::internal::LEVELMAX)
    if (_p[level]) _p[level]->zero();
}

void siconos::modeling::LagrangianDS::computePostImpactVelocity() {
  // When this function is call, q[1] is supposed to be pre-impact velocity.
  // We solve M(v+ - v-) = p - The result is saved in(place of) p[1].
  DEBUG_BEGIN("siconos::modeling::LagrangianDS::computePostImpactV()\n");
  siconos::algebra::SiconosVector tmp(*_p[1]);
  if (_inverseMass) _inverseMass->Solve(tmp);
  *_q[1] += tmp;  // v+ = v- + p
  DEBUG_BEGIN("siconos::modeling::LagrangianDS::computePostImpactV() END \n");
}
void siconos::modeling::LagrangianDS::allocateFExt() {
  if (!_fExt) _fExt = std::make_shared<siconos::algebra::SiconosVector>(_ndof);
}
void siconos::modeling::LagrangianDS::allocateFInt() {
  if (!_fInt) _fInt = std::make_shared<siconos::algebra::SiconosVector>(_ndof);
}

void siconos::modeling::LagrangianDS::setComputeMassFunction(const std::string& pluginPath,
                                                             const std::string& functionName) {
  _pluginMass->setComputeFunction(pluginPath, functionName);
  if (!_mass) _mass = std::make_shared<siconos::algebra::SimpleMatrix>(_ndof, _ndof);
  _hasConstantMass = false;
}

void siconos::modeling::LagrangianDS::setComputeMassFunction(siconos::plugins::FPtr7 fct) {
  _pluginMass->setComputeFunction((void*)fct);
  if (!_mass) _mass = std::make_shared<siconos::algebra::SimpleMatrix>(_ndof, _ndof);
  _hasConstantMass = false;
}

void siconos::modeling::LagrangianDS::setComputeFIntFunction(const std::string& pluginPath,
                                                             const std::string& functionName) {
  _pluginFInt->setComputeFunction(pluginPath, functionName);
  allocateFInt();
  //    Plugin::setFunction(&computeFIntPtr, pluginPath,functionName);
}

void siconos::modeling::LagrangianDS::setComputeFIntFunction(siconos::plugins::FPtr6 fct) {
  _pluginFInt->setComputeFunction((void*)fct);
  allocateFInt();
  //    computeFIntPtr = fct;
}

void siconos::modeling::LagrangianDS::setComputeFExtFunction(const std::string& pluginPath,
                                                             const std::string& functionName) {
  _pluginFExt->setComputeFunction(pluginPath, functionName);
  if (!_fExt) _fExt = std::make_shared<siconos::algebra::SiconosVector>(_ndof);
  _hasConstantFExt = false;
}

/** set a specified function to compute fExt
 *
 *  \param fct a pointer on the plugin function
 */
void siconos::modeling::LagrangianDS::setComputeFExtFunction(
    siconos::plugins::VectorFunctionOfTime fct) {
  _pluginFExt->setComputeFunction((void*)fct);
  if (!_fExt) _fExt = std::make_shared<siconos::algebra::SiconosVector>(_ndof);
  //   computeFExtPtr = fct ;
  _hasConstantFExt = false;
}

void siconos::modeling::LagrangianDS::setComputeFGyrFunction(const std::string& pluginPath,
                                                             const std::string& functionName) {
  _pluginFGyr->setComputeFunction(pluginPath, functionName);
  if (!_fGyr) _fGyr = std::make_shared<siconos::algebra::SiconosVector>(_ndof);
  init_forces();
}

void siconos::modeling::LagrangianDS::setComputeFGyrFunction(siconos::plugins::FPtr5 fct) {
  _pluginFGyr->setComputeFunction((void*)fct);
  if (!_fGyr) _fGyr = std::make_shared<siconos::algebra::SiconosVector>(_ndof);
  init_forces();
}

void siconos::modeling::LagrangianDS::allocateJacobianFIntq() {
  if (!_jacobianFIntq)
    _jacobianFIntq = std::make_shared<siconos::algebra::SimpleMatrix>(_ndof, _ndof);
}

void siconos::modeling::LagrangianDS::setComputeJacobianFIntqFunction(
    const std::string& pluginPath, const std::string& functionName) {
  _pluginJacqFInt->setComputeFunction(pluginPath, functionName);
  allocateJacobianFIntq();
  init_forces();
}

void siconos::modeling::LagrangianDS::allocateJacobianFIntqDot() {
  if (!_jacobianFIntqDot)
    _jacobianFIntqDot = std::make_shared<siconos::algebra::SimpleMatrix>(_ndof, _ndof);
}

void siconos::modeling::LagrangianDS::setComputeJacobianFIntqDotFunction(
    const std::string& pluginPath, const std::string& functionName) {
  _pluginJacqDotFInt->setComputeFunction(pluginPath, functionName);
  allocateJacobianFIntqDot();
  init_forces();
}

void siconos::modeling::LagrangianDS::setComputeJacobianFIntqFunction(
    siconos::plugins::FPtr6 fct) {
  _pluginJacqFInt->setComputeFunction((void*)fct);
  allocateJacobianFIntq();
  init_forces();
}

void siconos::modeling::LagrangianDS::setComputeJacobianFIntqDotFunction(
    siconos::plugins::FPtr6 fct) {
  _pluginJacqDotFInt->setComputeFunction((void*)fct);
  allocateJacobianFIntqDot();
  init_forces();
}

void siconos::modeling::LagrangianDS::setComputeJacobianFGyrqFunction(
    const std::string& pluginPath, const std::string& functionName) {
  _pluginJacqFGyr->setComputeFunction(pluginPath, functionName);
  if (!_jacobianFGyrq)
    _jacobianFGyrq = std::make_shared<siconos::algebra::SimpleMatrix>(_ndof, _ndof);
  init_forces();
}

void siconos::modeling::LagrangianDS::setComputeJacobianFGyrqDotFunction(
    const std::string& pluginPath, const std::string& functionName) {
  _pluginJacqDotFGyr->setComputeFunction(pluginPath, functionName);
  if (!_jacobianFGyrqDot)
    _jacobianFGyrqDot = std::make_shared<siconos::algebra::SimpleMatrix>(_ndof, _ndof);
  init_forces();
}

void siconos::modeling::LagrangianDS::setComputeJacobianFGyrqFunction(
    siconos::plugins::FPtr5 fct) {
  _pluginJacqFGyr->setComputeFunction((void*)fct);
  if (!_jacobianFGyrq)
    _jacobianFGyrq = std::make_shared<siconos::algebra::SimpleMatrix>(_ndof, _ndof);
  init_forces();
}

void siconos::modeling::LagrangianDS::setComputeJacobianFGyrqDotFunction(
    siconos::plugins::FPtr5 fct) {
  _pluginJacqDotFGyr->setComputeFunction((void*)fct);
  if (!_jacobianFGyrqDot)
    _jacobianFGyrqDot = std::make_shared<siconos::algebra::SimpleMatrix>(_ndof, _ndof);
}

double siconos::modeling::LagrangianDS::computeKineticEnergy() {
  DEBUG_BEGIN("NewtonEulerDS::computeKineticEnergy()\n");
  auto velo = velocity();
  assert(velo);
  DEBUG_EXPR(velo->display());

  auto tmp = std::make_shared<siconos::algebra::SiconosVector>(*velo);
  if (_mass) siconos::algebra::prod(*_mass, *velo, *tmp, true);

  double K = 0.5 * siconos::algebra::inner_prod(*tmp, *velo);

  DEBUG_PRINTF("Kinetic Energy = %e\n", K);
  DEBUG_END("siconos::modeling::LagrangianDS::computeKineticEnergy()\n");
  return K;
}

void siconos::modeling::LagrangianDS::acceptSP(
    std::shared_ptr<siconos::internal::SiconosVisitor> tourist) const {
  tourist->visit(*this);
}
