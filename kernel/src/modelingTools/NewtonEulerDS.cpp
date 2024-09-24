/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#include "NewtonEulerDS.hpp"

#include <boost/math/quaternion.hpp>
#include <iostream>

#include "BlockMatrix.hpp"
#include "RotationQuaternion.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixOp.hpp"        // setBlock
#include "SiconosMatrixVectorOp.hpp"  // mat-vec prod
#include "SiconosVector.hpp"
#include "SiconosVectorOp.hpp"  // inner_prod, setBlock
#include "SiconosVisitor.hpp"
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

namespace siconos::modeling {
static void computeJacobianConvectedVectorInBodyFrame(
    double q0, double q1, double q2, double q3,
    std::shared_ptr<SecondOrderDS::Matrix> jacobian,
    std::shared_ptr<siconos::algebra::SiconosVector> v) {
  /* This routine compute the jacobian with respect to p of R^T(p)v */
  jacobian->zero();

  double v0 = v->getValue(0);
  double v1 = v->getValue(1);
  double v2 = v->getValue(2);

  jacobian->setValue(0, 3, q0 * v0 + q3 * v1 - q2 * v2);
  jacobian->setValue(0, 4, q1 * v0 + q2 * v1 + q3 * v2);
  jacobian->setValue(0, 5, -q2 * v0 + q1 * v1 - q0 * v2);
  jacobian->setValue(0, 6, -q3 * v0 + q0 * v1 + q1 * v2);

  jacobian->setValue(1, 3, -q3 * v0 + q0 * v1 + q1 * v2);
  jacobian->setValue(1, 4, q2 * v0 - q1 * v1 + q0 * v2);
  jacobian->setValue(1, 5, q1 * v0 + q2 * v1 + q3 * v2);
  jacobian->setValue(1, 6, -q0 * v0 - q3 * v1 + q2 * v2);

  jacobian->setValue(2, 3, q2 * v0 - q1 * v1 + q0 * v2);
  jacobian->setValue(2, 4, q3 * v0 - q0 * v1 - q1 * v2);
  jacobian->setValue(2, 5, q0 * v0 + q3 * v1 - q2 * v2);
  jacobian->setValue(2, 6, q1 * v0 + q2 * v1 + q3 * v2);

  *jacobian *= 2.0;
}

/** This function has been added to avoid Swig director to wrap _MExt into numpy.array
 * when we call  siconos::modeling::NewtonEulerDS::computeMExt(double time,
 * std::shared_ptr<siconos::algebra::SiconosVector> q,
 * std::shared_ptr<siconos::algebra::SiconosVector> mExt) that calls in turn computeMExt(time,
 * q, _mExt);
 */
static void computeMExt_internal(
    double time, bool hasConstantMExt, unsigned int qDim,
    std::shared_ptr<siconos::algebra::MapVectorType> q0,
    std::shared_ptr<siconos::plugins::PluggedObject> pluginMExt,
    std::shared_ptr<siconos::algebra::SiconosVector> mExt_attributes,
    std::shared_ptr<siconos::algebra::SiconosVector> mExt) {
  /* if the pointer has been set to an external vector
   * after setting the plugin, we do not call the plugin */
  if (hasConstantMExt) {
    if (mExt != mExt_attributes) *mExt = *mExt_attributes;
  } else if (pluginMExt->fPtr)
    ((FExt_NE)pluginMExt->fPtr)(time, &(*mExt)(0), qDim,
                                &(*q0)(0));  // parameter z are assumed to be equal to q0
}

}  // namespace siconos::modeling

void siconos::modeling::computeT(std::shared_ptr<siconos::algebra::SiconosVector> q,
                                 std::shared_ptr<SecondOrderDS::Matrix> T) {
  DEBUG_BEGIN(
      "computeT(std::shared_ptr<siconos::algebra::SiconosVector> q, "
      "std::shared_ptr<Matrix> T)\n")
  //  std::cout <<"\n NewtonEulerDS::computeT(std::shared_ptr<siconos::algebra::SiconosVector>
  //  q)\n  " <<std::endl;
  double q0 = q->getValue(3) / 2.0;
  double q1 = q->getValue(4) / 2.0;
  double q2 = q->getValue(5) / 2.0;
  double q3 = q->getValue(6) / 2.0;
  T->setValue(3, 3, -q1);
  T->setValue(3, 4, -q2);
  T->setValue(3, 5, -q3);
  T->setValue(4, 3, q0);
  T->setValue(4, 4, -q3);
  T->setValue(4, 5, q2);
  T->setValue(5, 3, q3);
  T->setValue(5, 4, q0);
  T->setValue(5, 5, -q1);
  T->setValue(6, 3, -q2);
  T->setValue(6, 4, q1);
  T->setValue(6, 5, q0);
  DEBUG_END(
      "computeT(std::shared_ptr<siconos::algebra::SiconosVector> q, "
      "std::shared_ptr<Matrix> T)\n")
}

siconos::modeling::NewtonEulerDS::NewtonEulerDS(
    Eigen::Ref<siconos::algebra::SiconosVector> initial_position,
    Eigen::Ref<siconos::algebra::SiconosVector> initial_twist, double mass,
    Eigen::Ref<siconos::algebra::SiconosMatrix> inertia)
    : SecondOrderDS(13, 6), _scalarMass{mass} {
  DEBUG_BEGIN(
      "siconos::modeling::NewtonEulerDS::NewtonEulerDS(std::shared_ptr<siconos::algebra::"
      "SiconosVector> Q0, "
      "std::shared_ptr<siconos::algebra::SiconosVector> Twist0,double  mass, "
      "std::shared_ptr<Matrix> inertialMatrix)\n");

  _zeroPlugin();

  // Initial conditions
  // Q0 : contains the center of mass coordinate, and the quaternion initial. (dim(Q0)=7)
  // Twist0 : contains the initial velocity of center of mass and the omega initial.
  // (dim(VTwist0)=6)
  _q0 = std::make_shared<siconos::algebra::MapVectorType>(initial_position.data(),
                                                          initial_position.size());
  _twist0 = std::make_shared<siconos::algebra::MapVectorType>(initial_twist.data(),
                                                              initial_twist.size());

  // Current state
  _q = std::make_shared<siconos::algebra::SiconosVector>(*_q0);
  _twist = std::make_shared<siconos::algebra::SiconosVector>(*_twist0);
  _dotq = std::make_shared<siconos::algebra::SiconosVector>(_qDim);
  _dotq->setZero();

  /** \todo lazy Memory allocation */
  _p.resize(3);
  _p[1] = std::make_shared<siconos::algebra::SiconosVector>(ndof_);  // Needed in NewtonEulerR

  //  -- Mass --
  //   Remind that inertial matrix is accessed with inertialMatrix() method, a view on
  //   the bloxk
  mass_internal_storage_ = std::make_unique<std::vector<double>>(ndof_ * ndof_);
  mass_view_ = std::make_shared<MapType>(mass_internal_storage_->data(), ndof_, ndof_);
  mass_view_->setZero();
  (*mass_view_)(0, 0) = _scalarMass;
  (*mass_view_)(1, 1) = _scalarMass;
  (*mass_view_)(2, 2) = _scalarMass;
  mass_view_->block<3, 3>(3, 3) = inertia;

  _T = std::make_shared<Matrix>(_qDim, ndof_);

  // Allocate memory for forces and its jacobians.
  // Needed only for integrators with first-order formulation.
  if (!_wrench) _wrench = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
  if (!_mGyr) {
    _mGyr = std::make_shared<siconos::algebra::SiconosVector>(3);
    _mGyr->setZero();
  }

  // The follwing jacobian are always allocated since always have
  // Gyroscopical forces that have non linear forces
  // This should be remove if the integration is explicit or _nullifyMGyr(false) is set to
  // true?
  _jacobianMGyrtwist = std::make_shared<Matrix>(3, ndof_);
  _jacobianWrenchTwist = std::make_shared<Matrix>(ndof_, ndof_);

  _T->setZero();
  (*_T)(0, 0) = 1.0;
  (*_T)(1, 1) = 1.0;
  (*_T)(2, 2) = 1.0;
  siconos::modeling::computeT(_q, _T);

  DEBUG_END(
      "siconos::modeling::NewtonEulerDS::NewtonEulerDS(std::shared_ptr<siconos::algebra::"
      "SiconosVector> Q0, "
      "std::shared_ptr<siconos::algebra::SiconosVector> Twist0,double  mass, "
      "std::shared_ptr<Matrix> inertialMatrix)\n");
}

// void siconos::modeling::NewtonEulerDS::updateMassMatrix() {
//   // _mass->zero();
//   // _mass->setValue(0, 0, _scalarMass);
//   // _mass->setValue(1, 1, _scalarMass);
//   // _mass->setValue(2, 2, _scalarMass);

//   mass_view_->setZero();
//   (*mass_view_)(0, 0) = _scalarMass;
//   (*mass_view_)(1, 1) = _scalarMass;
//   (*mass_view_)(2, 2) = _scalarMass;
//   (*mass_view_)(3, 3) = (*inertialMatrix_view_)(0, 0);
//   (*mass_view_)(4, 4) = (*inertialMatrix_view_)(1, 1);
//   (*mass_view_)(5, 5) = (*inertialMatrix_view_)(2, 2);
// }

void siconos::modeling::NewtonEulerDS::_zeroPlugin() {
  _pluginMExt = std::make_shared<siconos::plugins::PluggedObject>();
  _pluginFInt = std::make_shared<siconos::plugins::PluggedObject>();
  _pluginMInt = std::make_shared<siconos::plugins::PluggedObject>();
  _pluginJacqFInt = std::make_shared<siconos::plugins::PluggedObject>();
  _pluginJactwistFInt = std::make_shared<siconos::plugins::PluggedObject>();
  _pluginJacqMInt = std::make_shared<siconos::plugins::PluggedObject>();
  _pluginJactwistMInt = std::make_shared<siconos::plugins::PluggedObject>();
}

void siconos::modeling::NewtonEulerDS::setInertia(double ix, double iy, double iz) {
  (*mass_view_)(3, 3) = ix;
  (*mass_view_)(4, 4) = iy;
  (*mass_view_)(5, 5) = iz;
}

void siconos::modeling::NewtonEulerDS::initializeNonSmoothInput(unsigned int level) {
  DEBUG_PRINTF(
      "siconos::modeling::NewtonEulerDS::initializeNonSmoothInput(unsigned int level) for "
      "level = %i\n",
      level);

  if (!_p[level]) {
    if (level == 0) {
      _p[level] = std::make_shared<siconos::algebra::SiconosVector>(_qDim);
    } else
      _p[level] = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
  }

#ifdef DEBUG_MESSAGES
  DEBUG_PRINT("display() after initialization");
  display();
#endif
}

void siconos::modeling::NewtonEulerDS::initRhs(double time) {
  DEBUG_BEGIN("siconos::modeling::NewtonEulerDS::initRhs(double time)\n");
  // dim
  _n = _qDim + 6;

  _x0 = std::make_shared<siconos::algebra::SiconosVector>(_q0->size() + _twist0->size());
  *_x0 << *_q0, *_twist0;

  _x[0] = std::make_shared<siconos::algebra::SiconosVector>(_q->size() + _twist->size());
  *(_x[0]) << *_q, *_twist;

  if (!_acceleration) _acceleration = std::make_shared<siconos::algebra::SiconosVector>(6);

  // Compute _dotq
  computeT();
  siconos::algebra::prod(*_T, *_twist, *_dotq, true);
  _x[1] =
      std::make_shared<siconos::algebra::SiconosVector>(_dotq->size() + _acceleration->size());
  *(_x[1]) << *_dotq, *_acceleration;

  // Nothing to do for the initialization of the wrench

  // Everything concerning rhs and its jacobian is handled in initRhs and computeXXX
  // related functions.
  _rhsMatrices.resize(numberOfRhsMatrices_);

  if (!_p[2]) _p[2] = std::make_shared<siconos::algebra::SiconosVector>(6);

  init_inverse_mass();

  computeRhs(time);

  /** \warning the derivative of T w.r.t to q is neglected */
  _rhsMatrices[jacobianXBloc00_] = std::make_shared<Matrix>(_qDim, _qDim);

  _rhsMatrices[jacobianXBloc01_] = std::make_shared<Matrix>(*_T);
  bool flag1 = false, flag2 = false;
  if (_jacobianWrenchq) {
    // Solve MjacobianX(1,0) = jacobianFL[0]
    computeJacobianqForces(time);

    _rhsMatrices[jacobianXBloc10_] = std::make_shared<Matrix>(*_jacobianWrenchq);
    algebra::solveInPlace(*_inverseMass, *_rhsMatrices[jacobianXBloc10_]);
    flag1 = true;
  }

  if (_jacobianWrenchTwist) {
    // Solve MjacobianX(1,1) = jacobianFL[1]
    computeJacobianvForces(time);
    _rhsMatrices[jacobianXBloc11_] = std::make_shared<Matrix>(*_jacobianWrenchTwist);
    algebra::solveInPlace(*_inverseMass, *_rhsMatrices[jacobianXBloc11_]);
    flag2 = true;
  }

  if (!_rhsMatrices[zeroMatrix_]) {
    _rhsMatrices[zeroMatrix_] = std::make_shared<Matrix>(6, 6);
    _rhsMatrices[zeroMatrix_]->setZero();
  }

  if (!_rhsMatrices[zeroMatrixqDim_]) {
    _rhsMatrices[zeroMatrixqDim_] = std::make_shared<Matrix>(6, _qDim);
    _rhsMatrices[zeroMatrixqDim_]->setZero();
  }

  if (flag1 && flag2)
    _jacxRhs = std::make_shared<siconos::algebra::BlockMatrix>(
        _rhsMatrices[jacobianXBloc00_], _rhsMatrices[jacobianXBloc01_],
        _rhsMatrices[jacobianXBloc10_], _rhsMatrices[jacobianXBloc11_]);
  else if (flag1)  // flag2 = false
    _jacxRhs = std::make_shared<siconos::algebra::BlockMatrix>(
        _rhsMatrices[jacobianXBloc00_], _rhsMatrices[jacobianXBloc01_],
        _rhsMatrices[jacobianXBloc10_], _rhsMatrices[zeroMatrix_]);
  else if (flag2)  // flag1 = false
    _jacxRhs = std::make_shared<siconos::algebra::BlockMatrix>(
        _rhsMatrices[jacobianXBloc00_], _rhsMatrices[jacobianXBloc01_],
        _rhsMatrices[zeroMatrixqDim_], _rhsMatrices[jacobianXBloc11_]);
  else
    _jacxRhs = std::make_shared<siconos::algebra::BlockMatrix>(
        _rhsMatrices[jacobianXBloc00_], _rhsMatrices[jacobianXBloc01_],
        _rhsMatrices[zeroMatrixqDim_], _rhsMatrices[zeroMatrix_]);
  DEBUG_EXPR(display(););
  DEBUG_END("siconos::modeling::NewtonEulerDS::initRhs(double time)\n");
}

void siconos::modeling::NewtonEulerDS::setQ(const siconos::algebra::SiconosVector& newValue) {
  if (newValue.size() != _qDim)
    THROW_EXCEPTION("NewtonEulerDS - setQ: inconsistent input vector size ");

  if (!_q)
    _q = std::make_shared<siconos::algebra::SiconosVector>(newValue);
  else
    *_q = newValue;
}

void siconos::modeling::NewtonEulerDS::setQPtr(
    std::shared_ptr<siconos::algebra::SiconosVector> newPtr) {
  if (newPtr->size() != _qDim)
    THROW_EXCEPTION("NewtonEulerDS - setQPtr: inconsistent input vector size ");
  _q = newPtr;
}

void siconos::modeling::NewtonEulerDS::setVelocity(
    const siconos::algebra::SiconosVector& newValue) {
  if (newValue.size() != ndof_)
    THROW_EXCEPTION("NewtonEulerDS - setVelocity: inconsistent input vector size ");

  if (!_twist) {
    if (!twist0_internal_storage) {
      twist0_internal_storage = std::make_unique<std::vector<double>>(
          newValue.data(), newValue.data() + newValue.size());
      _twist0 = std::make_shared<siconos::algebra::MapVectorType>(
          twist0_internal_storage->data(), newValue.size());
    } else
      std::copy(newValue.data(), newValue.data() + newValue.size(),
                twist0_internal_storage->begin());
  } else {
    *_twist = newValue;
  }
}

void siconos::modeling::NewtonEulerDS::setVelocityPtr(
    std::shared_ptr<siconos::algebra::SiconosVector> newPtr) {
  if (newPtr->size() != ndof_)
    THROW_EXCEPTION("NewtonEulerDS - setVelocityPtr: inconsistent input vector size ");
  _twist = newPtr;
}

void siconos::modeling::NewtonEulerDS::resetToInitialState() {
  // set q and q[1] to q0 and Twist0
  if (_q0) {
    *_q = *_q0;
  } else
    THROW_EXCEPTION(
        "siconos::modeling::NewtonEulerDS::resetToInitialState - initial position _q0 is "
        "null");

  if (_twist0) {
    *_twist = *_twist0;
  } else
    THROW_EXCEPTION(
        "siconos::modeling::NewtonEulerDS::resetToInitialState - initial twist _twist0 is "
        "null");
}

void siconos::modeling::NewtonEulerDS::init_inverse_mass() {
  if (mass_view_ && !_inverseMass) {
    // updateMassMatrix();
    _inverseMass = std::make_shared<siconos::algebra::SiconosMatrix>(*mass_view_);  // copy!
  }
}

void siconos::modeling::NewtonEulerDS::update_inverse_mass(double time) {
  if (mass_view_ && _inverseMass) {
    // updateMassMatrix();
    *_inverseMass = *mass_view_;  // copy!
  }
}

void siconos::modeling::NewtonEulerDS::setConstantFext(
    Eigen::Ref<siconos::algebra::SiconosVector> newValue) {
  /**  Must:

   - create the Map (view onto memory handled by newValue) for Fext_
   - set the corresponding booleans
   - reset internal storage (should already be null but who knows ...)
   */

  fext_internal_storage_ = nullptr;

  fext_view_ =
      std::make_shared<siconos::algebra::MapVectorType>(newValue.data(), newValue.size());
  hasFext_ = true;
  hasConstantFext_ = true;
  computefext_ = nullptr;
}

void siconos::modeling::NewtonEulerDS::setComputeFextFunction(
    ExternalForcesFunction fext_func) {
  // Ensure that memory is properly allocated for fext_
  if (!fext_internal_storage_) {
    fext_internal_storage_ = std::make_unique<std::vector<double>>(ndof_);
  }
  fext_view_ =
      std::make_shared<siconos::algebra::MapVectorType>(fext_internal_storage_->data(), ndof_);

  hasFext_ = true;
  hasConstantFext_ = false;
  computefext_ = fext_func;
}

void siconos::modeling::NewtonEulerDS::computeFext(double time) {
  if (computefext_)
    // in that case, internal_storage must have been allocated by
    // setCompute... call
    computefext_(time, *fext_view_);
}

void siconos::modeling::NewtonEulerDS::computeMExt(double time) {
  DEBUG_BEGIN("N3ewtonEulerDS::computeMExt(double time)\n");
  computeMExt_internal(time, _hasConstantMExt, _qDim, _q0, _pluginMExt, _mExt, _mExt);
  DEBUG_END("siconos::modeling::NewtonEulerDS::computeMExt(double time)\n");
}

void siconos::modeling::NewtonEulerDS::computeMExt(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> mExt) {
  DEBUG_BEGIN(
      "siconos::modeling::NewtonEulerDS::computeMExt(double time, "
      "std::shared_ptr<siconos::algebra::SiconosVector> mExt)\n");
  computeMExt_internal(time, _hasConstantMExt, _qDim, _q0, _pluginMExt, _mExt, mExt);
  DEBUG_END(
      "siconos::modeling::NewtonEulerDS::computeMExt(double time, "
      "std::shared_ptr<siconos::algebra::SiconosVector> mExt)\n");
}

void siconos::modeling::NewtonEulerDS::computeJacobianMExtqExpressedInInertialFrameByFD(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> q) {
  DEBUG_BEGIN(
      "siconos::modeling::NewtonEulerDS::computeJacobianMExtqExpressedInInertialFrameByFD("
      "...)"
      "\n");

  /* The computation of Jacobian of R^T mExt is somehow very rough since the pertubation
   * that we apply to q  that gives qeps does not provide a unit quaternion. The rotation
   * is computed assuming that the quaternion is unit (see quaternionRotate(double q0,
   * double q1, double q2, double q3, std::shared_ptr<siconos::algebra::SiconosVector> v)).
   */

  auto mExt = std::make_shared<siconos::algebra::SiconosVector>(3);
  computeMExt(time, mExt);
  if (_isMextExpressedInInertialFrame) siconos::geometry::changeFrameAbsToBody(*q, *mExt);
  DEBUG_EXPR(q->display());
  DEBUG_EXPR(mExt->display(););

  double mExt0 = mExt->getValue(0);
  double mExt1 = mExt->getValue(1);
  double mExt2 = mExt->getValue(2);

  auto qeps = std::make_shared<siconos::algebra::SiconosVector>(*q);
  _jacobianMExtq->zero();
  (*qeps)(3) += _epsilonFD;
  for (int j = 3; j < 7; j++) {
    computeMExt(time, mExt);
    if (_isMextExpressedInInertialFrame) siconos::geometry::changeFrameAbsToBody(*qeps, *mExt);
    DEBUG_EXPR(mExt->display(););
    _jacobianMExtq->setValue(0, j, (mExt->getValue(0) - mExt0) / _epsilonFD);
    _jacobianMExtq->setValue(1, j, (mExt->getValue(1) - mExt1) / _epsilonFD);
    _jacobianMExtq->setValue(2, j, (mExt->getValue(2) - mExt2) / _epsilonFD);
    (*qeps)(j) -= _epsilonFD;
    if (j < 6) (*qeps)(j + 1) += _epsilonFD;
  }
  DEBUG_EXPR(_jacobianMExtq->display(););
  DEBUG_END(
      "siconos::modeling::NewtonEulerDS::computeJacobianMExtqExpressedInInertialFrameByFD("
      "...)"
      "\n");
}

void siconos::modeling::NewtonEulerDS::computeJacobianMExtqExpressedInInertialFrame(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> q) {
  DEBUG_BEGIN(
      "siconos::modeling::NewtonEulerDS::computeJacobianMExtqExpressedInInertialFrame(...)"
      "\n");
  bool isMextExpressedInInertialFrame_save = _isMextExpressedInInertialFrame;
  _isMextExpressedInInertialFrame = false;
  auto mExt = std::make_shared<siconos::algebra::SiconosVector>(3);
  computeMExt(time, mExt);
  if (_isMextExpressedInInertialFrame) siconos::geometry::changeFrameAbsToBody(*q, *mExt);
  DEBUG_EXPR(q->display());
  DEBUG_EXPR(mExt->display());

  _isMextExpressedInInertialFrame = isMextExpressedInInertialFrame_save;

  double q0 = q->getValue(3);
  double q1 = q->getValue(4);
  double q2 = q->getValue(5);
  double q3 = q->getValue(6);

  computeJacobianConvectedVectorInBodyFrame(q0, q1, q2, q3, _jacobianMExtq, mExt);

  DEBUG_EXPR(_jacobianMExtq->display());

  // std::shared_ptr<Matrix> jacobianMExtqtmp (new
  // SiconosMatrix(*_jacobianMExtq));
  // computeJacobianMExtqExpressedInInertialFrameByFD(time, q);

  // std::cout << "#################  " << (*jacobianMExtqtmp- *_jacobianMExtq).normInf()
  // << std::endl; assert((*jacobianMExtqtmp- *_jacobianMExtq).normInf()< 1e-10);

  // DEBUG_EXPR(_jacobianMExtq->display(););
  DEBUG_END(
      "siconos::modeling::NewtonEulerDS::computeJacobianMExtqExpressedInInertialFrame(...)"
      "\n");
}
void siconos::modeling::NewtonEulerDS::computeFInt(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SiconosVector> v) {
  computeFInt(time, q, v, _fInt);
}

void siconos::modeling::NewtonEulerDS::computeFInt(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SiconosVector> v,
    std::shared_ptr<siconos::algebra::SiconosVector> fInt) {
  if (_pluginFInt->fPtr)
    ((FInt_NE)_pluginFInt->fPtr)(time, &(*q)(0), &(*v)(0), &(*fInt)(0), _qDim,
                                 &(*_q0)(0));  // parameter z are assumed to be equal to q0
}

void siconos::modeling::NewtonEulerDS::computeMInt(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SiconosVector> v) {
  DEBUG_BEGIN(
      "siconos::modeling::NewtonEulerDS::computeMInt(double time, "
      "std::shared_ptr<siconos::algebra::SiconosVector> q, "
      "std::shared_ptr<siconos::algebra::SiconosVector> v)\n");
  computeMInt(time, q, v, _mInt);
  DEBUG_END(
      "siconos::modeling::NewtonEulerDS::computeMInt(double time, "
      "std::shared_ptr<siconos::algebra::SiconosVector> q, "
      "std::shared_ptr<siconos::algebra::SiconosVector> v)\n");
}

void siconos::modeling::NewtonEulerDS::computeMInt(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SiconosVector> v,
    std::shared_ptr<siconos::algebra::SiconosVector> mInt) {
  DEBUG_BEGIN(
      "siconos::modeling::NewtonEulerDS::computeMInt(double time, "
      "std::shared_ptr<siconos::algebra::SiconosVector> q, "
      "std::shared_ptr<siconos::algebra::SiconosVector> v, "
      "std::shared_ptr<siconos::algebra::SiconosVector> mInt)\n");
  if (_pluginMInt->fPtr)
    ((FInt_NE)_pluginMInt->fPtr)(time, &(*q)(0), &(*v)(0), &(*mInt)(0), _qDim,
                                 &(*_q0)(0));  // parameter z are assumed to be equal to q0
  DEBUG_END(
      "siconos::modeling::NewtonEulerDS::computeMInt(double time, "
      "std::shared_ptr<siconos::algebra::SiconosVector> q, "
      "std::shared_ptr<siconos::algebra::SiconosVector> v, "
      "std::shared_ptr<siconos::algebra::SiconosVector> mInt)\n");
}

void siconos::modeling::NewtonEulerDS::computeJacobianFIntq(double time) {
  computeJacobianFIntq(time, _q, _twist);
}
void siconos::modeling::NewtonEulerDS::computeJacobianFIntv(double time) {
  computeJacobianFIntv(time, _q, _twist);
}

void siconos::modeling::NewtonEulerDS::computeJacobianFIntq(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SiconosVector> twist) {
  DEBUG_PRINT("siconos::modeling::NewtonEulerDS::computeJacobianFIntq(...) starts");
  if (_pluginJacqFInt->fPtr)
    ((FInt_NE)_pluginJacqFInt->fPtr)(time, &(*q)(0), &(*twist)(0), &(*_jacobianFIntq)(0, 0),
                                     _qDim, &(*_q0)(0));
  else if (_computeJacobianFIntqByFD)
    computeJacobianFIntqByFD(time, q, twist);
  DEBUG_EXPR(_jacobianFIntq->display(););
  DEBUG_END("siconos::modeling::NewtonEulerDS::computeJacobianFIntq(...)");
}

void siconos::modeling::NewtonEulerDS::computeJacobianFIntqByFD(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SiconosVector> twist) {
  DEBUG_BEGIN("siconos::modeling::NewtonEulerDS::computeJacobianFIntqByFD(...)\n");
  auto fInt = std::make_shared<siconos::algebra::SiconosVector>(3);
  computeFInt(time, q, twist, fInt);

  double fInt0 = fInt->getValue(0);
  double fInt1 = fInt->getValue(1);
  double fInt2 = fInt->getValue(2);

  auto qeps = std::make_shared<siconos::algebra::SiconosVector>(*q);
  _jacobianFIntq->zero();
  (*qeps)(0) += _epsilonFD;
  for (int j = 0; j < 7; j++) {
    computeFInt(time, qeps, twist, fInt);
    _jacobianFIntq->setValue(0, j, (fInt->getValue(0) - fInt0) / _epsilonFD);
    _jacobianFIntq->setValue(1, j, (fInt->getValue(1) - fInt1) / _epsilonFD);
    _jacobianFIntq->setValue(2, j, (fInt->getValue(2) - fInt2) / _epsilonFD);
    (*qeps)(j) -= _epsilonFD;
    if (j < 6) (*qeps)(j + 1) += _epsilonFD;
  }
  DEBUG_END("siconos::modeling::NewtonEulerDS::computeJacobianFIntqByFD(...)\n");
}

void siconos::modeling::NewtonEulerDS::computeJacobianFIntv(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SiconosVector> twist) {
  if (_pluginJactwistFInt->fPtr)
    ((FInt_NE)_pluginJactwistFInt->fPtr)(time, &(*q)(0), &(*twist)(0),
                                         &(*_jacobianFInttwist)(0, 0), _qDim, &(*_q0)(0));
  else if (_computeJacobianFInttwistByFD)
    computeJacobianFIntvByFD(time, q, twist);
}

void siconos::modeling::NewtonEulerDS::computeJacobianFIntvByFD(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SiconosVector> twist) {
  DEBUG_BEGIN("siconos::modeling::NewtonEulerDS::computeJacobianFIntvByFD(...)\n");
  auto fInt = std::make_shared<siconos::algebra::SiconosVector>(3);
  computeFInt(time, q, twist, fInt);

  double fInt0 = fInt->getValue(0);
  double fInt1 = fInt->getValue(1);
  double fInt2 = fInt->getValue(2);

  auto veps = std::make_shared<siconos::algebra::SiconosVector>(*twist);
  _jacobianFInttwist->zero();

  (*veps)(0) += _epsilonFD;
  for (int j = 0; j < 6; j++) {
    computeFInt(time, q, veps, fInt);
    _jacobianFInttwist->setValue(0, j, (fInt->getValue(0) - fInt0) / _epsilonFD);
    _jacobianFInttwist->setValue(1, j, (fInt->getValue(1) - fInt1) / _epsilonFD);
    _jacobianFInttwist->setValue(2, j, (fInt->getValue(2) - fInt2) / _epsilonFD);
    (*veps)(j) -= _epsilonFD;
    if (j < 5) (*veps)(j + 1) += _epsilonFD;
  }

  DEBUG_END("siconos::modeling::NewtonEulerDS::computeJacobianFIntvByFD(...)\n");
}
void siconos::modeling::NewtonEulerDS::computeJacobianMGyrtwistByFD(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SiconosVector> twist) {
  DEBUG_BEGIN("siconos::modeling::NewtonEulerDS::computeJacobianMGyrvByFD(...)\n");

  siconos::algebra::SiconosVector3 mGyr;
  computeMGyr(inertia_view(), *twist, mGyr);

  double mGyr0 = mGyr(0);
  double mGyr1 = mGyr(1);
  double mGyr2 = mGyr(2);

  siconos::algebra::SiconosVector veps{*twist};
  _jacobianMGyrtwist->zero();

  veps(0) += _epsilonFD;
  for (int j = 0; j < 6; j++) {
    computeMGyr(inertia_view(), veps, mGyr);
    _jacobianMGyrtwist->setValue(3, j, (mGyr(0) - mGyr0) / _epsilonFD);
    _jacobianMGyrtwist->setValue(4, j, (mGyr(1) - mGyr1) / _epsilonFD);
    _jacobianMGyrtwist->setValue(5, j, (mGyr(2) - mGyr2) / _epsilonFD);
    veps(j) -= _epsilonFD;
    if (j < 5) veps(j + 1) += _epsilonFD;
  }
  DEBUG_EXPR(_jacobianMGyrtwist->display());
  DEBUG_END("siconos::modeling::NewtonEulerDS::computeJacobianMGyrvByFD(...)\n");
}
void siconos::modeling::NewtonEulerDS::computeJacobianMIntq(double time) {
  computeJacobianMIntq(time, _q, _twist);
}
void siconos::modeling::NewtonEulerDS::computeJacobianMIntv(double time) {
  computeJacobianMIntv(time, _q, _twist);
}

void siconos::modeling::NewtonEulerDS::computeJacobianMIntq(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SiconosVector> twist) {
  DEBUG_PRINT("siconos::modeling::NewtonEulerDS::computeJacobianMIntq(...) starts");
  if (_pluginJacqMInt->fPtr)
    ((FInt_NE)_pluginJacqMInt->fPtr)(time, &(*q)(0), &(*twist)(0), &(*_jacobianMIntq)(0, 0),
                                     _qDim, &(*_q0)(0));
  else if (_computeJacobianMIntqByFD)
    computeJacobianMIntqByFD(time, q, twist);
  DEBUG_EXPR(_jacobianMIntq->display());
  DEBUG_PRINT("siconos::modeling::NewtonEulerDS::computeJacobianMIntq(...) ends");
}

void siconos::modeling::NewtonEulerDS::computeJacobianMIntqByFD(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SiconosVector> twist) {
  DEBUG_PRINT("siconos::modeling::NewtonEulerDS::computeJacobianMIntqByFD(...) starts\n");

  auto mInt = std::make_shared<siconos::algebra::SiconosVector>(3);
  computeMInt(time, q, twist, mInt);
  double mInt0 = mInt->getValue(0);
  double mInt1 = mInt->getValue(1);
  double mInt2 = mInt->getValue(2);

  auto qeps = std::make_shared<siconos::algebra::SiconosVector>(*q);

  (*qeps)(0) += _epsilonFD;
  for (int j = 0; j < 7; j++) {
    computeMInt(time, qeps, twist, mInt);
    _jacobianMIntq->setValue(0, j, (mInt->getValue(0) - mInt0) / _epsilonFD);
    _jacobianMIntq->setValue(1, j, (mInt->getValue(1) - mInt1) / _epsilonFD);
    _jacobianMIntq->setValue(2, j, (mInt->getValue(2) - mInt2) / _epsilonFD);
    (*qeps)(j) -= _epsilonFD;
    if (j < 6) (*qeps)(j + 1) += _epsilonFD;
  }
  DEBUG_PRINT("siconos::modeling::NewtonEulerDS::computeJacobianMIntqByFD(...) ends\n");
}

void siconos::modeling::NewtonEulerDS::computeJacobianMIntv(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SiconosVector> twist) {
  if (_pluginJactwistMInt->fPtr)
    ((FInt_NE)_pluginJactwistMInt->fPtr)(time, &(*q)(0), &(*twist)(0),
                                         &(*_jacobianMInttwist)(0, 0), _qDim, &(*_q0)(0));
  else if (_computeJacobianMInttwistByFD)
    computeJacobianMIntvByFD(time, q, twist);
}

void siconos::modeling::NewtonEulerDS::computeJacobianMIntvByFD(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SiconosVector> twist) {
  DEBUG_PRINT("siconos::modeling::NewtonEulerDS::computeJacobianMIntvByFD(...) starts\n");

  auto mInt = std::make_shared<siconos::algebra::SiconosVector>(3);
  computeMInt(time, q, twist, mInt);
  double mInt0 = mInt->getValue(0);
  double mInt1 = mInt->getValue(1);
  double mInt2 = mInt->getValue(2);

  auto veps = std::make_shared<siconos::algebra::SiconosVector>(*twist);

  (*veps)(0) += _epsilonFD;
  for (int j = 0; j < 6; j++) {
    computeMInt(time, q, veps, mInt);
    _jacobianMInttwist->setValue(0, j, (mInt->getValue(0) - mInt0) / _epsilonFD);
    _jacobianMInttwist->setValue(1, j, (mInt->getValue(1) - mInt1) / _epsilonFD);
    _jacobianMInttwist->setValue(2, j, (mInt->getValue(2) - mInt2) / _epsilonFD);
    (*veps)(j) -= _epsilonFD;
    if (j < 5) (*veps)(j + 1) += _epsilonFD;
  }
  DEBUG_PRINT("siconos::modeling::NewtonEulerDS::computeJacobianMIntvByFD(...) ends\n");
}

void siconos::modeling::NewtonEulerDS::computeRhs(double time) {
  DEBUG_BEGIN("siconos::modeling::NewtonEulerDS::computeRhs(double time)");
  *_acceleration = *(_p[2]);  // Warning: r/p update is done in Interactions/Relations

  computeForces(time, _q, _twist);
  *_acceleration += *_wrench;
  DEBUG_EXPR(_wrench->display(););

  if (_inverseMass) algebra::solveInPlace(*_inverseMass, *_acceleration);

  // Compute _dotq
  computeT();
  siconos::algebra::prod(*_T, *_twist, *_dotq, true);

  _x[1]->setBlock(0, *_dotq);
  _x[1]->setBlock(_qDim, *_acceleration);
}

void siconos::modeling::NewtonEulerDS::computeJacobianRhsx(double time) {
  if (_jacobianWrenchq) {
    std::shared_ptr<Matrix> bloc10 = _jacxRhs->block(1, 0);
    computeJacobianqForces(time);
    *bloc10 = *_jacobianWrenchq;
    algebra::solveInPlace(*_inverseMass, *bloc10);
  }
  if (_jacobianWrenchTwist) {
    std::shared_ptr<Matrix> bloc11 = _jacxRhs->block(1, 1);
    computeJacobianvForces(time);
    *bloc11 = *_jacobianWrenchTwist;
    algebra::solveInPlace(*_inverseMass, *bloc11);
  }
}

void siconos::modeling::NewtonEulerDS::computeForces(double time) {
  computeForces(time, _q, _twist);
}

void siconos::modeling::NewtonEulerDS::computeMGyr(
    Eigen::Ref<siconos::algebra::SiconosMatrix> inertia_matrix,
    Eigen::Ref<siconos::algebra::SiconosVector> twist,
    Eigen::Ref<siconos::algebra::SiconosVector> mGyr) {
  DEBUG_EXPR(twist->display());
  siconos::algebra::SiconosVector3 omega;
  omega << twist(3), twist(4), twist(5);
  siconos::algebra::SiconosVector3 iomega;
  iomega = inertia_matrix * omega;
  mGyr = omega.cross(iomega);
}

void siconos::modeling::NewtonEulerDS::computeForces(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SiconosVector> twist) {
  DEBUG_BEGIN(
      "siconos::modeling::NewtonEulerDS::computeForces(double time, "
      "std::shared_ptr<siconos::algebra::SiconosVector> q, "
      "std::shared_ptr<siconos::algebra::SiconosVector> twist)\n")

  if (_wrench) {
    _wrench->zero();

    // External wrench

    if (hasFext_) {
      computeFext(time);
      assert(!std::isnan(fext_view_->sum()));
      _wrench->setBlock(0, *fext_view_);
    }
    if (_mExt) {
      computeMExt(time);
      assert(!std::isnan(_mExt->vector_sum()));
      if (_isMextExpressedInInertialFrame) {
        auto mExt = std::make_shared<siconos::algebra::SiconosVector>(*_mExt);
        siconos::geometry::changeFrameAbsToBody(*q, *mExt);
        _wrench->setBlock(3, *mExt);
      } else
        _wrench->setBlock(3, *_mExt);
    }

    // Internal wrench

    if (_fInt) {
      computeFInt(time, q, twist);
      assert(!std::isnan(_fInt->vector_sum()));
      _wrench->setValue(0, _wrench->getValue(0) - _fInt->getValue(0));
      _wrench->setValue(1, _wrench->getValue(1) - _fInt->getValue(1));
      _wrench->setValue(2, _wrench->getValue(2) - _fInt->getValue(2));
    }

    if (_mInt) {
      computeMInt(time, q, twist);
      assert(!std::isnan(_mInt->vector_sum()));
      _wrench->setValue(3, _wrench->getValue(3) - _mInt->getValue(0));
      _wrench->setValue(4, _wrench->getValue(4) - _mInt->getValue(1));
      _wrench->setValue(5, _wrench->getValue(5) - _mInt->getValue(2));
    }

    // Gyroscopical effect
    if (!_nullifyMGyr) {
      computeMGyr(inertia_view(), *twist, *_mGyr);
      assert(!std::isnan(_mGyr->vector_sum()));
      _wrench->setValue(3, _wrench->getValue(3) - _mGyr->getValue(0));
      _wrench->setValue(4, _wrench->getValue(4) - _mGyr->getValue(1));
      _wrench->setValue(5, _wrench->getValue(5) - _mGyr->getValue(2));
    }
    DEBUG_EXPR(_wrench->display());
    DEBUG_END(
        "siconos::modeling::NewtonEulerDS::computeForces(double time, "
        "std::shared_ptr<siconos::algebra::SiconosVector> q, "
        "std::shared_ptr<siconos::algebra::SiconosVector> twist)\n")
  } else {
    THROW_EXCEPTION("siconos::modeling::NewtonEulerDS::computeForces _wrench is null");
  }
  // else nothing.
}

void siconos::modeling::NewtonEulerDS::computeJacobianqForces(double time) {
  DEBUG_BEGIN("siconos::modeling::NewtonEulerDS::computeJacobianqWrench(double time) \n");
  if (_jacobianWrenchq) {
    _jacobianWrenchq->zero();
    if (_jacobianFIntq) {
      computeJacobianFIntq(time);
      std::dynamic_pointer_cast<Matrix>(_jacobianWrenchq)
          ->setBlock(0, 0, -1.0 * *_jacobianFIntq);
    }
    if (_jacobianMIntq) {
      computeJacobianMIntq(time);
    }
    if (_isMextExpressedInInertialFrame && _mExt) {
      computeJacobianMExtqExpressedInInertialFrame(time, _q);
      std::dynamic_pointer_cast<Matrix>(_jacobianWrenchq)
          ->setBlock(3, 0, 1.0 * *_jacobianMExtq);
    }
    DEBUG_EXPR(_jacobianWrenchq->display(););
  } else {
    THROW_EXCEPTION(
        "siconos::modeling::NewtonEulerDS::computeJacobianqForces _jacobianWrenchq is "
        "null");
  }
  // else nothing.
  DEBUG_END("siconos::modeling::NewtonEulerDS::computeJacobianqForces(double time) \n");
}

void siconos::modeling::NewtonEulerDS::computeJacobianvForces(double time) {
  DEBUG_BEGIN("siconos::modeling::NewtonEulerDS::computeJacobiantwistForces(double time) \n");
  if (_jacobianWrenchTwist) {
    _jacobianWrenchTwist->zero();
    if (_jacobianFInttwist) {
      computeJacobianFIntv(time);
      std::dynamic_pointer_cast<Matrix>(_jacobianWrenchTwist)
          ->setBlock(0, 0, -1.0 * *_jacobianFInttwist);
    }
    if (_jacobianMInttwist) {
      computeJacobianMIntv(time);
      std::dynamic_pointer_cast<Matrix>(_jacobianWrenchTwist)
          ->setBlock(3, 0, -1.0 * *_jacobianMInttwist);
    }
    if (!_nullifyMGyr) {
      if (_jacobianMGyrtwist) {
        // computeJacobianMGyrtwistByFD(time,_q,_twist);
        computeJacobianMGyrtwist(time);
        std::dynamic_pointer_cast<Matrix>(_jacobianWrenchTwist)
            ->setBlock(3, 0, -1.0 * *_jacobianMGyrtwist);
      }
    }
  }
  // else nothing.
  DEBUG_END("siconos::modeling::NewtonEulerDS::computeJacobiantwistForces(double time) \n");
}

void siconos::modeling::NewtonEulerDS::computeJacobianMGyrtwist(double time) {
  DEBUG_BEGIN("siconos::modeling::NewtonEulerDS::computeJacobianMGyrtwist(double time) \n");
  if (_jacobianMGyrtwist) {
    // Omega /\ I \Omega:
    _jacobianMGyrtwist->setZero();

    siconos::algebra::SiconosVector3 omega;
    omega << (*_twist)(3), (*_twist)(4), (*_twist)(5);
    siconos::algebra::SiconosVector3 iomega;
    iomega = inertia_view() * omega;
    siconos::algebra::SiconosVector3 ei;
    siconos::algebra::SiconosVector3 iei;
    siconos::algebra::SiconosVector3 ei_iomega;
    siconos::algebra::SiconosVector3 omega_iei;

    /*See equation of DevNotes.pdf, equation with label eq:NE_nablaFL1*/
    for (int i = 0; i < 3; i++) {
      ei(i) = 1.0;
      iei = inertia_view() * ei;
      omega_iei = omega.cross(iei);
      ei_iomega = ei.cross(iomega);
      for (int j = 0; j < 3; j++)
        _jacobianMGyrtwist->setValue(j, 3 + i, ei_iomega(j) + omega_iei(j));
    }
    // Check if Jacobian is valid. Warning to the transpose operation in
    // _jacobianMGyrtwist->setValue(3 + j, 3 + i, ei_Iomega.getValue(j) +
    // omega_Iei.getValue(j));
  }
  // else nothing.
  DEBUG_EXPR(_jacobianMGyrtwist->display());
  // _jacobianMGyrtwist->display();
  // std::shared_ptr<Matrix> jacobianMGyrtmp (new
  // SiconosMatrix(*_jacobianMGyrtwist)); computeJacobianMGyrtwistByFD(time, _q, _twist);
  // jacobianMGyrtmp->display();
  // std::cout << "#################  " << (*jacobianMGyrtmp -
  // *_jacobianMGyrtwist).normInf()
  // << std::endl; assert((*jacobianMGyrtmp - *_jacobianMGyrtwist).normInf()< 1e-10);
  DEBUG_END("siconos::modeling::NewtonEulerDS::computeJacobianMGyrtwist(double time) \n");
}

void siconos::modeling::NewtonEulerDS::display(bool brief) const {
  std::cout << "=====> NewtonEuler System display (number: " << _number << ")." << std::endl;
  std::cout << "- q " << std::endl;
  if (_q)
    _q->display();
  else
    std::cout << "-> nullptr" << std::endl;
  std::cout << "- q0 " << std::endl;
  if (_q0) _q0->display();
  std::cout << "- twist " << std::endl;
  if (_twist)
    _twist->display();
  else
    std::cout << "-> nullptr" << std::endl;
  std::cout << "- twist0 " << std::endl;
  if (_twist0)
    _twist0->display();
  else
    std::cout << "-> nullptr" << std::endl;
  std::cout << "- dotq " << std::endl;
  if (_dotq)
    _dotq->display();
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
  std::cout << "mass :" << _scalarMass << std::endl;
  std::cout << "Inertia :" << inertia_view() << "\n";
  std::cout << "===================================== " << std::endl;
}

void siconos::modeling::NewtonEulerDS::setIsMextExpressedInInertialFrame(bool value) {
  _isMextExpressedInInertialFrame = value;
  if (!_jacobianMExtq) _jacobianMExtq = std::make_shared<Matrix>(3, _qDim);
  if (!_jacobianWrenchq) _jacobianWrenchq = std::make_shared<Matrix>(ndof_, _qDim);
}

// --- Functions for memory handling ---
void siconos::modeling::NewtonEulerDS::initMemory(unsigned int steps) {
  DynamicalSystem::initMemory(steps);

  if (steps == 0)
    std::cout
        << "Warning : siconos::modeling::NewtonEulerDS::initMemory with size equal to zero"
        << std::endl;
  else {
    _qMemory.setMemorySize(steps, _qDim);
    _twistMemory.setMemorySize(steps, ndof_);
    _forcesMemory.setMemorySize(steps, ndof_);
    _dotqMemory.setMemorySize(steps, _qDim);
    //    swapInMemory(); Useless, done in osi->initializeWorkVectorsForDS
  }
}

void siconos::modeling::NewtonEulerDS::swapInMemory() {
  //  _xMemory->swap(_x[0]);
  _qMemory.swap(*_q);
  _twistMemory.swap(*_twist);
  _dotqMemory.swap(*_dotq);
  _forcesMemory.swap(*_wrench);
}

void siconos::modeling::NewtonEulerDS::resetAllNonSmoothParts() {
  if (_p[1])
    _p[1]->zero();
  else
    _p[1] = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
}
void siconos::modeling::NewtonEulerDS::resetNonSmoothPart(unsigned int level) {
  if (_p[level]) _p[level]->zero();
}

void siconos::modeling::NewtonEulerDS::computeT() { siconos::modeling::computeT(_q, _T); }

void siconos::modeling::NewtonEulerDS::computeTdot() {
  if (!_Tdot) {
    _Tdot = std::make_shared<Matrix>(_qDim, ndof_);
    _Tdot->zero();
  }

  siconos::modeling::computeT(_dotq, _Tdot);
}

void siconos::modeling::NewtonEulerDS::normalizeq() { siconos::geometry::normalizeq(*_q); }

void siconos::modeling::NewtonEulerDS::setComputeMExtFunction(
    const std::string& pluginPath, const std::string& functionName) {
  _pluginMExt->setComputeFunction(pluginPath, functionName);
  if (!_mExt) _mExt = std::make_shared<siconos::algebra::SiconosVector>(3, 0);
  _hasConstantMExt = false;
}

void siconos::modeling::NewtonEulerDS::setComputeMExtFunction(FExt_NE fct) {
  _pluginMExt->setComputeFunction((void*)fct);
  if (!_mExt) _mExt = std::make_shared<siconos::algebra::SiconosVector>(3, 0);
  _hasConstantMExt = false;
}

void siconos::modeling::NewtonEulerDS::setComputeFIntFunction(
    const std::string& pluginPath, const std::string& functionName) {
  _pluginFInt->setComputeFunction(pluginPath, functionName);
  if (!_fInt) _fInt = std::make_shared<siconos::algebra::SiconosVector>(3, 0);
}

void siconos::modeling::NewtonEulerDS::setComputeMIntFunction(
    const std::string& pluginPath, const std::string& functionName) {
  _pluginMInt->setComputeFunction(pluginPath, functionName);
  if (!_mInt) _mInt = std::make_shared<siconos::algebra::SiconosVector>(3, 0);
}

void siconos::modeling::NewtonEulerDS::setComputeFIntFunction(FInt_NE fct) {
  _pluginFInt->setComputeFunction((void*)fct);
  if (!_fInt) _fInt = std::make_shared<siconos::algebra::SiconosVector>(3, 0);
}

void siconos::modeling::NewtonEulerDS::setComputeMIntFunction(FInt_NE fct) {
  _pluginMInt->setComputeFunction((void*)fct);
  if (!_mInt) _mInt = std::make_shared<siconos::algebra::SiconosVector>(3, 0);
}

void siconos::modeling::NewtonEulerDS::setComputeJacobianFIntqFunction(
    const std::string& pluginPath, const std::string& functionName) {
  //    Plugin::setFunction(&computeJacobianFIntqPtr, pluginPath,functionName);
  _pluginJacqFInt->setComputeFunction(pluginPath, functionName);
  if (!_jacobianFIntq) _jacobianFIntq = std::make_shared<Matrix>(3, _qDim);
  if (!_jacobianWrenchq) _jacobianWrenchq = std::make_shared<Matrix>(ndof_, _qDim);
  _computeJacobianFIntqByFD = false;
}
void siconos::modeling::NewtonEulerDS::setComputeJacobianFIntvFunction(
    const std::string& pluginPath, const std::string& functionName) {
  //    Plugin::setFunction(&computeJacobianFIntvPtr, pluginPath,functionName);
  _pluginJactwistFInt->setComputeFunction(pluginPath, functionName);
  if (!_jacobianFInttwist) _jacobianFInttwist = std::make_shared<Matrix>(3, ndof_);
  _computeJacobianFInttwistByFD = false;
}
void siconos::modeling::NewtonEulerDS::setComputeJacobianFIntqFunction(FInt_NE fct) {
  _pluginJacqFInt->setComputeFunction((void*)fct);
  if (!_jacobianFIntq) _jacobianFIntq = std::make_shared<Matrix>(3, _qDim);
  if (!_jacobianWrenchq) _jacobianWrenchq = std::make_shared<Matrix>(ndof_, _qDim);
  _computeJacobianFIntqByFD = false;
}
void siconos::modeling::NewtonEulerDS::setComputeJacobianFIntvFunction(FInt_NE fct) {
  _pluginJactwistFInt->setComputeFunction((void*)fct);
  if (!_jacobianFInttwist) _jacobianFInttwist = std::make_shared<Matrix>(3, ndof_);
  if (!_jacobianWrenchTwist) _jacobianWrenchTwist = std::make_shared<Matrix>(ndof_, ndof_);
  _computeJacobianFInttwistByFD = false;
}

void siconos::modeling::NewtonEulerDS::setComputeJacobianMIntqFunction(
    const std::string& pluginPath, const std::string& functionName) {
  _pluginJacqMInt->setComputeFunction(pluginPath, functionName);
  if (!_jacobianMIntq) _jacobianMIntq = std::make_shared<Matrix>(3, _qDim);
  if (!_jacobianWrenchq) _jacobianWrenchq = std::make_shared<Matrix>(ndof_, _qDim);
  _computeJacobianMIntqByFD = false;
}
void siconos::modeling::NewtonEulerDS::setComputeJacobianMIntvFunction(
    const std::string& pluginPath, const std::string& functionName) {
  _pluginJactwistMInt->setComputeFunction(pluginPath, functionName);
  if (!_jacobianMInttwist) _jacobianMInttwist = std::make_shared<Matrix>(3, ndof_);
  if (!_jacobianWrenchTwist) _jacobianWrenchTwist = std::make_shared<Matrix>(ndof_, ndof_);
  _computeJacobianMInttwistByFD = false;
}
void siconos::modeling::NewtonEulerDS::setComputeJacobianMIntqFunction(FInt_NE fct) {
  _pluginJacqMInt->setComputeFunction((void*)fct);
  if (!_jacobianMIntq) _jacobianMIntq = std::make_shared<Matrix>(3, _qDim);
  if (!_jacobianWrenchq) _jacobianWrenchq = std::make_shared<Matrix>(ndof_, _qDim);
  _computeJacobianMIntqByFD = false;
}
void siconos::modeling::NewtonEulerDS::setComputeJacobianMIntvFunction(FInt_NE fct) {
  _pluginJactwistMInt->setComputeFunction((void*)fct);
  if (!_jacobianMInttwist) _jacobianMInttwist = std::make_shared<Matrix>(3, ndof_);
  if (!_jacobianWrenchTwist) _jacobianWrenchTwist = std::make_shared<Matrix>(ndof_, ndof_);
  _computeJacobianMInttwistByFD = false;
}

double siconos::modeling::NewtonEulerDS::computeKineticEnergy() {
  DEBUG_BEGIN("siconos::modeling::NewtonEulerDS::computeKineticEnergy()\n");
  assert(_twist);
  assert(mass_view_);
  DEBUG_EXPR(_twist->display());
  DEBUG_EXPR(_mass->display());

  auto tmp = *mass_view_ * *_twist;

  double K = 0.5 * tmp.dot(*_twist);

  DEBUG_PRINTF("Kinetic Energy = %e\n", K);
  DEBUG_END("siconos::modeling::NewtonEulerDS::computeKineticEnergy()\n");
  return K;
}

std::shared_ptr<siconos::algebra::SiconosVector>
siconos::modeling::NewtonEulerDS::linearVelocity(bool absoluteRef) const {
  // Short-cut: return the _twist 6-vector without modification, first
  // 3 components are the expected linear velocity.
  if (absoluteRef) return _twist;

  auto v = std::make_shared<siconos::algebra::SiconosVector>(3);
  linearVelocity(absoluteRef, *v);
  return v;
}

void siconos::modeling::NewtonEulerDS::linearVelocity(
    bool absoluteRef, siconos::algebra::SiconosVector& v) const {
  v(0) = (*_twist)(0);
  v(1) = (*_twist)(1);
  v(2) = (*_twist)(2);

  /* See _twist: linear velocity is in absolute frame */
  if (!absoluteRef) siconos::geometry::changeFrameAbsToBody(*_q, v);
}

std::shared_ptr<siconos::algebra::SiconosVector>
siconos::modeling::NewtonEulerDS::angularVelocity(bool absoluteRef) const {
  auto w = std::make_shared<siconos::algebra::SiconosVector>(3);
  angularVelocity(absoluteRef, *w);
  return w;
}

void siconos::modeling::NewtonEulerDS::angularVelocity(
    bool absoluteRef, siconos::algebra::SiconosVector& w) const {
  w(0) = (*_twist)(3);
  w(1) = (*_twist)(4);
  w(2) = (*_twist)(5);

  /* See _twist: angular velocity is in relative frame */
  if (absoluteRef) siconos::geometry::changeFrameBodyToAbs(*_q, w);
}

void siconos::modeling::computeExtForceAtPos(
    std::shared_ptr<siconos::algebra::SiconosVector> q, bool isMextExpressedInInertialFrame,
    std::shared_ptr<siconos::algebra::SiconosVector> force, bool forceAbsRef,
    std::shared_ptr<siconos::algebra::SiconosVector> pos, bool posAbsRef,
    Eigen::Ref<siconos::algebra::MapVectorType> fext,
    std::shared_ptr<siconos::algebra::SiconosVector> mExt, bool accumulate) {
  assert(fext.size() == 3);
  assert(!!force && force->size() == 3);
  if (pos) assert(!!mExt && mExt->size() == 3);

  siconos::algebra::SiconosVector abs_frc(*force), local_frc(*force);

  if (forceAbsRef) {
    if (pos) siconos::geometry::changeFrameAbsToBody(*q, local_frc);
  } else
    siconos::geometry::changeFrameBodyToAbs(*q, abs_frc);

  if (pos) {
    assert(!!mExt && mExt->size() >= 3);
    siconos::algebra::SiconosVector moment(3);
    if (posAbsRef) {
      siconos::algebra::SiconosVector local_pos(*pos);
      local_pos(0) -= (*q)(0);
      local_pos(1) -= (*q)(1);
      local_pos(2) -= (*q)(2);
      siconos::geometry::changeFrameAbsToBody(*q, local_pos);
      siconos::algebra::cross_product(local_pos, local_frc, moment);
    } else {
      siconos::algebra::cross_product(*pos, local_frc, moment);
    }

    if (isMextExpressedInInertialFrame) siconos::geometry::changeFrameBodyToAbs(*q, moment);

    if (accumulate)
      *mExt = *mExt + moment;
    else
      *mExt = moment;
  }

  if (accumulate)
    fext += fext + abs_frc;
  else
    fext = abs_frc;
}

void siconos::modeling::NewtonEulerDS::addExtForceAtPos(
    std::shared_ptr<siconos::algebra::SiconosVector> force, bool forceAbsRef,
    std::shared_ptr<siconos::algebra::SiconosVector> pos, bool posAbsRef) {
  assert(fext_view_->size() == 3);
  assert(!!force && force->size() == 3);
  if (pos) assert(!!_mExt && _mExt->size() == 3);

  computeExtForceAtPos(_q, _isMextExpressedInInertialFrame, force, forceAbsRef, pos, posAbsRef,
                       *fext_view_, _mExt, true);
}

void siconos::modeling::NewtonEulerDS::acceptSP(
    std::shared_ptr<siconos::internal::SiconosVisitor> tourist) const {
  tourist->visit(*this);
}
