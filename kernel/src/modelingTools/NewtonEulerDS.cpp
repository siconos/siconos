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
#include "BlockVector.hpp"
#include "RotationQuaternion.hpp"
#include "SiconosAlgebraProd.hpp"
#include "SiconosMatrixSetBlock.hpp"

// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

static void computeJacobianConvectedVectorInBodyFrame(double q0, double q1, double q2,
                                                      double q3, SP::SimpleMatrix jacobian,
                                                      SP::SiconosVector v) {
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

void computeT(SP::SiconosVector q, SP::SimpleMatrix T) {
  DEBUG_BEGIN("computeT(SP::SiconosVector q, SP::SimpleMatrix T)\n")
  //  std::cout <<"\n NewtonEulerDS::computeT(SP::SiconosVector q)\n  " <<std::endl;
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
  DEBUG_END("computeT(SP::SiconosVector q, SP::SimpleMatrix T)\n")
}

// From a set of data; Mass filled-in directly from a siconosMatrix -
// This constructor leads to the minimum NewtonEuler System form: \f$ M\ddot q = p \f$
/*
Q0 : contains the center of mass coordinate, and the quaternion initial. (dim(Q0)=7)
Twist0 : contains the initial velocity of center of mass and the omega initial.
(dim(VTwist0)=6)
*/
NewtonEulerDS::NewtonEulerDS()
    : SecondOrderDS(13, 6),
      _hasConstantFExt(false),
      _hasConstantMExt(false),
      _isMextExpressedInInertialFrame(false),
      _nullifyMGyr(false),
      _computeJacobianFIntqByFD(true),
      _computeJacobianFInttwistByFD(true),
      _computeJacobianMIntqByFD(true),
      _computeJacobianMInttwistByFD(true),
      _epsilonFD(sqrt(std::numeric_limits<double>::epsilon())) {
  /* common code for constructors
   * would be better to use delagation of constructors in c++11
   */
  _init();
}

NewtonEulerDS::NewtonEulerDS(SP::SiconosVector Q0, SP::SiconosVector Twist0, double mass,
                             SP::SiconosMatrix inertialMatrix)
    : SecondOrderDS(13, 6),
      _hasConstantFExt(false),
      _hasConstantMExt(false),
      _isMextExpressedInInertialFrame(false),
      _nullifyMGyr(false),
      _computeJacobianFIntqByFD(true),
      _computeJacobianFInttwistByFD(true),
      _computeJacobianMIntqByFD(true),
      _computeJacobianMInttwistByFD(true),
      _epsilonFD(sqrt(std::numeric_limits<double>::epsilon()))

{
  DEBUG_BEGIN(
      "NewtonEulerDS::NewtonEulerDS(SP::SiconosVector Q0, SP::SiconosVector Twist0,double  "
      "mass, SP::SiconosMatrix inertialMatrix)\n");

  /* common code for constructors
   * would be better to use delegation of constructors in c++11
   */
  _init();

  // Initial conditions
  _q0 = Q0;
  _twist0 = Twist0;
  resetToInitialState();

  _scalarMass = mass;
  if (inertialMatrix) _I = inertialMatrix;
  updateMassMatrix();

  _T->zero();
  _T->setValue(0, 0, 1.0);
  _T->setValue(1, 1, 1.0);
  _T->setValue(2, 2, 1.0);
  computeT();

  DEBUG_END(
      "NewtonEulerDS::NewtonEulerDS(SP::SiconosVector Q0, SP::SiconosVector Twist0,double  "
      "mass, SP::SiconosMatrix inertialMatrix)\n");
}

void NewtonEulerDS::init_forces() {
  // Allocate memory for forces and its jacobians.
  // Needed only for integrators with first-order formulation.

  if (!_wrench) _wrench.reset(new SiconosVector(_ndof));
  if (!_mGyr) _mGyr.reset(new SiconosVector(3, 0.0));

  /** The follwing jacobian are always allocated since we have always
   * Gyroscopical forces that has non linear forces
   * This should be remove if the integration is explicit or _nullifyMGyr(false) is set to true
   * ?
   */
  _jacobianMGyrtwist.reset(new SimpleMatrix(3, _ndof));
  _jacobianWrenchTwist.reset(new SimpleMatrix(_ndof, _ndof));
}

void NewtonEulerDS::_init() {
  // --- NEWTONEULER INHERITED CLASS MEMBERS ---
  // -- Memory allocation for vector and matrix members --

  _qDim = 7;
  _ndof = 6;
  _n = _qDim + _ndof;

  _zeroPlugin();

  // Current state
  _q.reset(new SiconosVector(_qDim));
  _twist.reset(new SiconosVector(_ndof));
  _dotq.reset(new SiconosVector(_qDim));

  /** \todo lazy Memory allocation */
  _p.resize(3);
  _p[1].reset(new SiconosVector(_ndof));  // Needed in NewtonEulerR

  _mass.reset(new SimpleMatrix(_ndof, _ndof));
  _mass->zero();
  _T.reset(new SimpleMatrix(_qDim, _ndof));

  _scalarMass = 1.;
  _I.reset(new SimpleMatrix(3, 3));
  _I->eye();
  updateMassMatrix();

  init_forces();

  // We initialize _z with a null vector of size 1, since z is required in plug-in functions
  // call.
  _z.reset(new SiconosVector(1));
}

void NewtonEulerDS::updateMassMatrix() {
  // _mass->zero();
  // _mass->setValue(0, 0, _scalarMass);
  // _mass->setValue(1, 1, _scalarMass);
  // _mass->setValue(2, 2, _scalarMass);

  _mass->eye();
  *_mass *= _scalarMass;

  Index dimIndex(2);
  dimIndex[0] = 3;
  dimIndex[1] = 3;
  Index startIndex(4);
  startIndex[0] = 0;
  startIndex[1] = 0;
  startIndex[2] = 3;
  startIndex[3] = 3;
  setBlock(_I, _mass, dimIndex, startIndex);
}

void NewtonEulerDS::_zeroPlugin() {
  _pluginFExt.reset(new PluggedObject());
  _pluginMExt.reset(new PluggedObject());
  _pluginFInt.reset(new PluggedObject());
  _pluginMInt.reset(new PluggedObject());
  _pluginJacqFInt.reset(new PluggedObject());
  _pluginJactwistFInt.reset(new PluggedObject());
  _pluginJacqMInt.reset(new PluggedObject());
  _pluginJactwistMInt.reset(new PluggedObject());
}

// Destructor
NewtonEulerDS::~NewtonEulerDS() {}

void NewtonEulerDS::setInertia(double ix, double iy, double iz) {
  _I->zero();

  (*_I)(0, 0) = ix;
  (*_I)(1, 1) = iy;
  (*_I)(2, 2) = iz;

  updateMassMatrix();
}

void NewtonEulerDS::initializeNonSmoothInput(unsigned int level) {
  DEBUG_PRINTF("NewtonEulerDS::initializeNonSmoothInput(unsigned int level) for level = %i\n",
               level);

  if (!_p[level]) {
    if (level == 0) {
      _p[level].reset(new SiconosVector(_qDim));
    } else
      _p[level].reset(new SiconosVector(_ndof));
  }

#ifdef DEBUG_MESSAGES
  DEBUG_PRINT("display() after initialization");
  display();
#endif
}

void NewtonEulerDS::initRhs(double time) {
  DEBUG_BEGIN("NewtonEulerDS::initRhs(double time)\n");
  // dim
  _n = _qDim + 6;

  _x0.reset(new SiconosVector(*_q0, *_twist0));

  _x[0].reset(new SiconosVector(*_q, *_twist));

  if (!_acceleration) _acceleration.reset(new SiconosVector(6));

  // Compute _dotq
  computeT();
  prod(*_T, *_twist, *_dotq, true);
  _x[1].reset(new SiconosVector(*_dotq, *_acceleration));

  // Nothing to do for the initialization of the wrench

  // Everything concerning rhs and its jacobian is handled in initRhs and computeXXX related
  // functions.
  _rhsMatrices.resize(numberOfRhsMatrices);

  if (!_p[2]) _p[2].reset(new SiconosVector(6));

  init_inverse_mass();

  computeRhs(time);

  /** \warning the derivative of T w.r.t to q is neglected */
  _rhsMatrices[jacobianXBloc00].reset(new SimpleMatrix(_qDim, _qDim, siconos::ZERO));

  _rhsMatrices[jacobianXBloc01].reset(new SimpleMatrix(*_T));
  bool flag1 = false, flag2 = false;
  if (_jacobianWrenchq) {
    // Solve MjacobianX(1,0) = jacobianFL[0]
    computeJacobianqForces(time);

    _rhsMatrices[jacobianXBloc10].reset(new SimpleMatrix(*_jacobianWrenchq));
    _inverseMass->Solve(*_rhsMatrices[jacobianXBloc10]);
    flag1 = true;
  }

  if (_jacobianWrenchTwist) {
    // Solve MjacobianX(1,1) = jacobianFL[1]
    computeJacobianvForces(time);
    _rhsMatrices[jacobianXBloc11].reset(new SimpleMatrix(*_jacobianWrenchTwist));
    _inverseMass->Solve(*_rhsMatrices[jacobianXBloc11]);
    flag2 = true;
  }

  if (!_rhsMatrices[zeroMatrix])
    _rhsMatrices[zeroMatrix].reset(new SimpleMatrix(6, 6, siconos::ZERO));

  if (!_rhsMatrices[zeroMatrixqDim])
    _rhsMatrices[zeroMatrixqDim].reset(new SimpleMatrix(6, _qDim, siconos::ZERO));

  if (flag1 && flag2)
    _jacxRhs.reset(
        new BlockMatrix(_rhsMatrices[jacobianXBloc00], _rhsMatrices[jacobianXBloc01],
                        _rhsMatrices[jacobianXBloc10], _rhsMatrices[jacobianXBloc11]));
  else if (flag1)  // flag2 = false
    _jacxRhs.reset(new BlockMatrix(_rhsMatrices[jacobianXBloc00],
                                   _rhsMatrices[jacobianXBloc01],
                                   _rhsMatrices[jacobianXBloc10], _rhsMatrices[zeroMatrix]));
  else if (flag2)  // flag1 = false
    _jacxRhs.reset(new BlockMatrix(_rhsMatrices[jacobianXBloc00],
                                   _rhsMatrices[jacobianXBloc01], _rhsMatrices[zeroMatrixqDim],
                                   _rhsMatrices[jacobianXBloc11]));
  else
    _jacxRhs.reset(new BlockMatrix(_rhsMatrices[jacobianXBloc00],
                                   _rhsMatrices[jacobianXBloc01], _rhsMatrices[zeroMatrixqDim],
                                   _rhsMatrices[zeroMatrix]));
  DEBUG_EXPR(display(););
  DEBUG_END("NewtonEulerDS::initRhs(double time)\n");
}

void NewtonEulerDS::setQ(const SiconosVector& newValue) {
  if (newValue.size() != _qDim)
    THROW_EXCEPTION("NewtonEulerDS - setQ: inconsistent input vector size ");

  if (!_q)
    _q.reset(new SiconosVector(newValue));
  else
    *_q = newValue;
}

void NewtonEulerDS::setQPtr(SP::SiconosVector newPtr) {
  if (newPtr->size() != _qDim)
    THROW_EXCEPTION("NewtonEulerDS - setQPtr: inconsistent input vector size ");
  _q = newPtr;
}
void NewtonEulerDS::setQ0(const SiconosVector& newValue) {
  if (newValue.size() != _qDim)
    THROW_EXCEPTION("NewtonEulerDS - setQ0: inconsistent input vector size ");

  if (!_q0)
    _q0.reset(new SiconosVector(newValue));
  else
    *_q0 = newValue;
}

void NewtonEulerDS::setQ0Ptr(SP::SiconosVector newPtr) {
  if (newPtr->size() != _qDim)
    THROW_EXCEPTION("NewtonEulerDS - setQ0Ptr: inconsistent input vector size ");
  _q0 = newPtr;
}
void NewtonEulerDS::setVelocity(const SiconosVector& newValue) {
  if (newValue.size() != _ndof)
    THROW_EXCEPTION("NewtonEulerDS - setVelocity: inconsistent input vector size ");

  if (!_twist)
    _twist0.reset(new SiconosVector(newValue));
  else
    *_twist = newValue;
}

void NewtonEulerDS::setVelocityPtr(SP::SiconosVector newPtr) {
  if (newPtr->size() != _ndof)
    THROW_EXCEPTION("NewtonEulerDS - setVelocityPtr: inconsistent input vector size ");
  _twist = newPtr;
}

void NewtonEulerDS::setVelocity0(const SiconosVector& newValue) {
  if (newValue.size() != _ndof)
    THROW_EXCEPTION("NewtonEulerDS - setVelocity0: inconsistent input vector size ");

  if (!_twist0)
    _twist0.reset(new SiconosVector(newValue));
  else
    *_twist0 = newValue;
}

void NewtonEulerDS::setVelocity0Ptr(SP::SiconosVector newPtr) {
  if (newPtr->size() != _ndof)
    THROW_EXCEPTION("NewtonEulerDS - setVelocity0Ptr: inconsistent input vector size ");
  _twist0 = newPtr;
}

void NewtonEulerDS::resetToInitialState() {
  // set q and q[1] to q0 and Twist0
  if (_q0) {
    *_q = *_q0;
  } else
    THROW_EXCEPTION("NewtonEulerDS::resetToInitialState - initial position _q0 is null");

  if (_twist0) {
    *_twist = *_twist0;
  } else
    THROW_EXCEPTION("NewtonEulerDS::resetToInitialState - initial twist _twist0 is null");
}

void NewtonEulerDS::init_inverse_mass() {
  if (_mass && !_inverseMass) {
    updateMassMatrix();
    _inverseMass.reset(new SimpleMatrix(*_mass));
  }
}

void NewtonEulerDS::update_inverse_mass() {
  if (_mass && _inverseMass) {
    updateMassMatrix();
    *_inverseMass = *_mass;
  }
}

void NewtonEulerDS::computeFExt(double time) {
  // computeFExt(time, _fExt);

  if (_pluginFExt->fPtr) {
    ((FExt_NE)_pluginFExt->fPtr)(time, &(*_fExt)(0), _qDim,
                                 &(*_q0)(0));  // parameter z are assumed to be equal to q0
  }
}

void NewtonEulerDS::computeFExt(double time, SP::SiconosVector fExt) {
  /* if the pointer has been set to an external vector
   * after setting the plugin, we do not call the plugin */
  if (_hasConstantFExt) {
    if (fExt != _fExt) *fExt = *_fExt;
  } else {
    if (_pluginFExt->fPtr) {
      ((FExt_NE)_pluginFExt->fPtr)(time, &(*fExt)(0), _qDim,
                                   &(*_q0)(0));  // parameter z are assumed to be equal to q0
    }
  }
}

/** This function has been added to avoid Swig director to wrap _MExt into numpy.array
 * when we call  NewtonEulerDS::computeMExt(double time, SP::SiconosVector q, SP::SiconosVector
 * mExt) that calls in turn computeMExt(time, q, _mExt);
 */
static void computeMExt_internal(double time, bool hasConstantMExt, unsigned int qDim,
                                 SP::SiconosVector q0, SP::PluggedObject pluginMExt,
                                 SP::SiconosVector mExt_attributes, SP::SiconosVector mExt) {
  /* if the pointer has been set to an external vector
   * after setting the plugin, we do not call the plugin */
  if (hasConstantMExt) {
    if (mExt != mExt_attributes) *mExt = *mExt_attributes;
  } else if (pluginMExt->fPtr)
    ((FExt_NE)pluginMExt->fPtr)(time, &(*mExt)(0), qDim,
                                &(*q0)(0));  // parameter z are assumed to be equal to q0
}

void NewtonEulerDS::computeMExt(double time) {
  DEBUG_BEGIN("N3ewtonEulerDS::computeMExt(double time)\n");
  computeMExt_internal(time, _hasConstantMExt, _qDim, _q0, _pluginMExt, _mExt, _mExt);
  DEBUG_END("NewtonEulerDS::computeMExt(double time)\n");
}

void NewtonEulerDS::computeMExt(double time, SP::SiconosVector mExt) {
  DEBUG_BEGIN("NewtonEulerDS::computeMExt(double time, SP::SiconosVector mExt)\n");
  computeMExt_internal(time, _hasConstantMExt, _qDim, _q0, _pluginMExt, _mExt, mExt);
  DEBUG_END("NewtonEulerDS::computeMExt(double time, SP::SiconosVector mExt)\n");
}

void NewtonEulerDS::computeJacobianMExtqExpressedInInertialFrameByFD(double time,
                                                                     SP::SiconosVector q) {
  DEBUG_BEGIN("NewtonEulerDS::computeJacobianMExtqExpressedInInertialFrameByFD(...)\n");

  /* The computation of Jacobian of R^T mExt is somehow very rough since the pertubation
   * that we apply to q  that gives qeps does not provide a unit quaternion. The rotation
   * is computed assuming that the quaternion is unit (see quaternionRotate(double q0, double
   * q1, double q2, double q3, SP::SiconosVector v)).
   */

  SP::SiconosVector mExt(new SiconosVector(3));
  computeMExt(time, mExt);
  if (_isMextExpressedInInertialFrame) ::changeFrameAbsToBody(q, mExt);
  DEBUG_EXPR(q->display());
  DEBUG_EXPR(mExt->display(););

  double mExt0 = mExt->getValue(0);
  double mExt1 = mExt->getValue(1);
  double mExt2 = mExt->getValue(2);

  SP::SiconosVector qeps(new SiconosVector(*q));
  _jacobianMExtq->zero();
  (*qeps)(3) += _epsilonFD;
  for (int j = 3; j < 7; j++) {
    computeMExt(time, mExt);
    if (_isMextExpressedInInertialFrame) ::changeFrameAbsToBody(qeps, mExt);
    DEBUG_EXPR(mExt->display(););
    _jacobianMExtq->setValue(0, j, (mExt->getValue(0) - mExt0) / _epsilonFD);
    _jacobianMExtq->setValue(1, j, (mExt->getValue(1) - mExt1) / _epsilonFD);
    _jacobianMExtq->setValue(2, j, (mExt->getValue(2) - mExt2) / _epsilonFD);
    (*qeps)(j) -= _epsilonFD;
    if (j < 6) (*qeps)(j + 1) += _epsilonFD;
  }
  DEBUG_EXPR(_jacobianMExtq->display(););
  DEBUG_END("NewtonEulerDS::computeJacobianMExtqExpressedInInertialFrameByFD(...)\n");
}

void NewtonEulerDS::computeJacobianMExtqExpressedInInertialFrame(double time,
                                                                 SP::SiconosVector q) {
  DEBUG_BEGIN("NewtonEulerDS::computeJacobianMExtqExpressedInInertialFrame(...)\n");
  bool isMextExpressedInInertialFrame_save = _isMextExpressedInInertialFrame;
  _isMextExpressedInInertialFrame = false;
  SP::SiconosVector mExt(new SiconosVector(3));
  computeMExt(time, mExt);
  if (_isMextExpressedInInertialFrame) ::changeFrameAbsToBody(q, mExt);
  DEBUG_EXPR(q->display());
  DEBUG_EXPR(mExt->display());

  _isMextExpressedInInertialFrame = isMextExpressedInInertialFrame_save;

  double q0 = q->getValue(3);
  double q1 = q->getValue(4);
  double q2 = q->getValue(5);
  double q3 = q->getValue(6);

  computeJacobianConvectedVectorInBodyFrame(q0, q1, q2, q3, _jacobianMExtq, mExt);

  DEBUG_EXPR(_jacobianMExtq->display());

  // SP::SimpleMatrix jacobianMExtqtmp (new SimpleMatrix(*_jacobianMExtq));
  // computeJacobianMExtqExpressedInInertialFrameByFD(time, q);

  // std::cout << "#################  " << (*jacobianMExtqtmp- *_jacobianMExtq).normInf() <<
  // std::endl; assert((*jacobianMExtqtmp- *_jacobianMExtq).normInf()< 1e-10);

  // DEBUG_EXPR(_jacobianMExtq->display(););
  DEBUG_END("NewtonEulerDS::computeJacobianMExtqExpressedInInertialFrame(...)\n");
}
void NewtonEulerDS::computeFInt(double time, SP::SiconosVector q, SP::SiconosVector v) {
  computeFInt(time, q, v, _fInt);
}

void NewtonEulerDS::computeFInt(double time, SP::SiconosVector q, SP::SiconosVector v,
                                SP::SiconosVector fInt) {
  if (_pluginFInt->fPtr)
    ((FInt_NE)_pluginFInt->fPtr)(time, &(*q)(0), &(*v)(0), &(*fInt)(0), _qDim,
                                 &(*_q0)(0));  // parameter z are assumed to be equal to q0
}

void NewtonEulerDS::computeMInt(double time, SP::SiconosVector q, SP::SiconosVector v) {
  DEBUG_BEGIN(
      "NewtonEulerDS::computeMInt(double time, SP::SiconosVector q, SP::SiconosVector v)\n");
  computeMInt(time, q, v, _mInt);
  DEBUG_END(
      "NewtonEulerDS::computeMInt(double time, SP::SiconosVector q, SP::SiconosVector v)\n");
}

void NewtonEulerDS::computeMInt(double time, SP::SiconosVector q, SP::SiconosVector v,
                                SP::SiconosVector mInt) {
  DEBUG_BEGIN(
      "NewtonEulerDS::computeMInt(double time, SP::SiconosVector q, SP::SiconosVector v, "
      "SP::SiconosVector mInt)\n");
  if (_pluginMInt->fPtr)
    ((FInt_NE)_pluginMInt->fPtr)(time, &(*q)(0), &(*v)(0), &(*mInt)(0), _qDim,
                                 &(*_q0)(0));  // parameter z are assumed to be equal to q0
  DEBUG_END(
      "NewtonEulerDS::computeMInt(double time, SP::SiconosVector q, SP::SiconosVector v, "
      "SP::SiconosVector mInt)\n");
}

void NewtonEulerDS::computeJacobianFIntq(double time) {
  computeJacobianFIntq(time, _q, _twist);
}
void NewtonEulerDS::computeJacobianFIntv(double time) {
  computeJacobianFIntv(time, _q, _twist);
}

void NewtonEulerDS::computeJacobianFIntq(double time, SP::SiconosVector q,
                                         SP::SiconosVector twist) {
  DEBUG_PRINT("NewtonEulerDS::computeJacobianFIntq(...) starts");
  if (_pluginJacqFInt->fPtr)
    ((FInt_NE)_pluginJacqFInt->fPtr)(time, &(*q)(0), &(*twist)(0), &(*_jacobianFIntq)(0, 0),
                                     _qDim, &(*_q0)(0));
  else if (_computeJacobianFIntqByFD)
    computeJacobianFIntqByFD(time, q, twist);
  DEBUG_EXPR(_jacobianFIntq->display(););
  DEBUG_END("NewtonEulerDS::computeJacobianFIntq(...)");
}

void NewtonEulerDS::computeJacobianFIntqByFD(double time, SP::SiconosVector q,
                                             SP::SiconosVector twist) {
  DEBUG_BEGIN("NewtonEulerDS::computeJacobianFIntqByFD(...)\n");
  SP::SiconosVector fInt(new SiconosVector(3));
  computeFInt(time, q, twist, fInt);

  double fInt0 = fInt->getValue(0);
  double fInt1 = fInt->getValue(1);
  double fInt2 = fInt->getValue(2);

  SP::SiconosVector qeps(new SiconosVector(*q));
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
  DEBUG_END("NewtonEulerDS::computeJacobianFIntqByFD(...)\n");
}

void NewtonEulerDS::computeJacobianFIntv(double time, SP::SiconosVector q,
                                         SP::SiconosVector twist) {
  if (_pluginJactwistFInt->fPtr)
    ((FInt_NE)_pluginJactwistFInt->fPtr)(time, &(*q)(0), &(*twist)(0),
                                         &(*_jacobianFInttwist)(0, 0), _qDim, &(*_q0)(0));
  else if (_computeJacobianFInttwistByFD)
    computeJacobianFIntvByFD(time, q, twist);
}

void NewtonEulerDS::computeJacobianFIntvByFD(double time, SP::SiconosVector q,
                                             SP::SiconosVector twist) {
  DEBUG_BEGIN("NewtonEulerDS::computeJacobianFIntvByFD(...)\n");
  SP::SiconosVector fInt(new SiconosVector(3));
  computeFInt(time, q, twist, fInt);

  double fInt0 = fInt->getValue(0);
  double fInt1 = fInt->getValue(1);
  double fInt2 = fInt->getValue(2);

  SP::SiconosVector veps(new SiconosVector(*twist));
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

  DEBUG_END("NewtonEulerDS::computeJacobianFIntvByFD(...)\n");
}
void NewtonEulerDS::computeJacobianMGyrtwistByFD(double time, SP::SiconosVector q,
                                                 SP::SiconosVector twist) {
  DEBUG_BEGIN("NewtonEulerDS::computeJacobianMGyrvByFD(...)\n");
  SP::SiconosVector mGyr(new SiconosVector(3));
  computeMGyr(twist, mGyr);

  double mGyr0 = mGyr->getValue(0);
  double mGyr1 = mGyr->getValue(1);
  double mGyr2 = mGyr->getValue(2);

  SP::SiconosVector veps(new SiconosVector(*twist));
  _jacobianMGyrtwist->zero();

  (*veps)(0) += _epsilonFD;
  for (int j = 0; j < 6; j++) {
    computeMGyr(veps, mGyr);
    _jacobianMGyrtwist->setValue(3, j, (mGyr->getValue(0) - mGyr0) / _epsilonFD);
    _jacobianMGyrtwist->setValue(4, j, (mGyr->getValue(1) - mGyr1) / _epsilonFD);
    _jacobianMGyrtwist->setValue(5, j, (mGyr->getValue(2) - mGyr2) / _epsilonFD);
    (*veps)(j) -= _epsilonFD;
    if (j < 5) (*veps)(j + 1) += _epsilonFD;
  }
  DEBUG_EXPR(_jacobianMGyrtwist->display());
  DEBUG_END("NewtonEulerDS::computeJacobianMGyrvByFD(...)\n");
}
void NewtonEulerDS::computeJacobianMIntq(double time) {
  computeJacobianMIntq(time, _q, _twist);
}
void NewtonEulerDS::computeJacobianMIntv(double time) {
  computeJacobianMIntv(time, _q, _twist);
}

void NewtonEulerDS::computeJacobianMIntq(double time, SP::SiconosVector q,
                                         SP::SiconosVector twist) {
  DEBUG_PRINT("NewtonEulerDS::computeJacobianMIntq(...) starts");
  if (_pluginJacqMInt->fPtr)
    ((FInt_NE)_pluginJacqMInt->fPtr)(time, &(*q)(0), &(*twist)(0), &(*_jacobianMIntq)(0, 0),
                                     _qDim, &(*_q0)(0));
  else if (_computeJacobianMIntqByFD)
    computeJacobianMIntqByFD(time, q, twist);
  DEBUG_EXPR(_jacobianMIntq->display());
  DEBUG_PRINT("NewtonEulerDS::computeJacobianMIntq(...) ends");
}

void NewtonEulerDS::computeJacobianMIntqByFD(double time, SP::SiconosVector q,
                                             SP::SiconosVector twist) {
  DEBUG_PRINT("NewtonEulerDS::computeJacobianMIntqByFD(...) starts\n");

  SP::SiconosVector mInt(new SiconosVector(3));
  computeMInt(time, q, twist, mInt);
  double mInt0 = mInt->getValue(0);
  double mInt1 = mInt->getValue(1);
  double mInt2 = mInt->getValue(2);

  SP::SiconosVector qeps(new SiconosVector(*q));

  (*qeps)(0) += _epsilonFD;
  for (int j = 0; j < 7; j++) {
    computeMInt(time, qeps, twist, mInt);
    _jacobianMIntq->setValue(0, j, (mInt->getValue(0) - mInt0) / _epsilonFD);
    _jacobianMIntq->setValue(1, j, (mInt->getValue(1) - mInt1) / _epsilonFD);
    _jacobianMIntq->setValue(2, j, (mInt->getValue(2) - mInt2) / _epsilonFD);
    (*qeps)(j) -= _epsilonFD;
    if (j < 6) (*qeps)(j + 1) += _epsilonFD;
  }
  DEBUG_PRINT("NewtonEulerDS::computeJacobianMIntqByFD(...) ends\n");
}

void NewtonEulerDS::computeJacobianMIntv(double time, SP::SiconosVector q,
                                         SP::SiconosVector twist) {
  if (_pluginJactwistMInt->fPtr)
    ((FInt_NE)_pluginJactwistMInt->fPtr)(time, &(*q)(0), &(*twist)(0),
                                         &(*_jacobianMInttwist)(0, 0), _qDim, &(*_q0)(0));
  else if (_computeJacobianMInttwistByFD)
    computeJacobianMIntvByFD(time, q, twist);
}

void NewtonEulerDS::computeJacobianMIntvByFD(double time, SP::SiconosVector q,
                                             SP::SiconosVector twist) {
  DEBUG_PRINT("NewtonEulerDS::computeJacobianMIntvByFD(...) starts\n");

  SP::SiconosVector mInt(new SiconosVector(3));
  computeMInt(time, q, twist, mInt);
  double mInt0 = mInt->getValue(0);
  double mInt1 = mInt->getValue(1);
  double mInt2 = mInt->getValue(2);

  SP::SiconosVector veps(new SiconosVector(*twist));

  (*veps)(0) += _epsilonFD;
  for (int j = 0; j < 6; j++) {
    computeMInt(time, q, veps, mInt);
    _jacobianMInttwist->setValue(0, j, (mInt->getValue(0) - mInt0) / _epsilonFD);
    _jacobianMInttwist->setValue(1, j, (mInt->getValue(1) - mInt1) / _epsilonFD);
    _jacobianMInttwist->setValue(2, j, (mInt->getValue(2) - mInt2) / _epsilonFD);
    (*veps)(j) -= _epsilonFD;
    if (j < 5) (*veps)(j + 1) += _epsilonFD;
  }
  DEBUG_PRINT("NewtonEulerDS::computeJacobianMIntvByFD(...) ends\n");
}

void NewtonEulerDS::computeRhs(double time) {
  DEBUG_BEGIN("NewtonEulerDS::computeRhs(double time)");
  *_acceleration = *(_p[2]);  // Warning: r/p update is done in Interactions/Relations

  computeForces(time, _q, _twist);
  *_acceleration += *_wrench;
  DEBUG_EXPR(_wrench->display(););

  if (_inverseMass) _inverseMass->Solve(*_acceleration);

  // Compute _dotq
  computeT();
  prod(*_T, *_twist, *_dotq, true);

  _x[1]->setBlock(0, *_dotq);
  _x[1]->setBlock(_qDim, *_acceleration);
}

void NewtonEulerDS::computeJacobianRhsx(double time) {
  if (_jacobianWrenchq) {
    SP::SiconosMatrix bloc10 = _jacxRhs->block(1, 0);
    computeJacobianqForces(time);
    *bloc10 = *_jacobianWrenchq;
    _inverseMass->Solve(*bloc10);
  }
  if (_jacobianWrenchTwist) {
    SP::SiconosMatrix bloc11 = _jacxRhs->block(1, 1);
    computeJacobianvForces(time);
    *bloc11 = *_jacobianWrenchTwist;
    _inverseMass->Solve(*bloc11);
  }
}

void NewtonEulerDS::computeForces(double time) { computeForces(time, _q, _twist); }

/** This function has been added to avoid Swig director to wrap _mGyr into numpy.array
 * when we call  NewtonEulerDS::computeMGyr(SP::SiconosVector twist) that calls in turn
 * computeMGyr(twist, _mGyr)
 */
static void computeMGyr_internal(SP::SiconosMatrix I, SP::SiconosVector twist,
                                 SP::SiconosVector mGyr) {
  if (I) {
    DEBUG_EXPR(I->display());
    DEBUG_EXPR(twist->display());
    SiconosVector omega(3);
    SiconosVector iomega(3);
    omega.setValue(0, twist->getValue(3));
    omega.setValue(1, twist->getValue(4));
    omega.setValue(2, twist->getValue(5));
    prod(*I, omega, iomega, true);
    cross_product(omega, iomega, *mGyr);
  }
}
void NewtonEulerDS::computeMGyr(SP::SiconosVector twist, SP::SiconosVector mGyr) {
  // computation of \Omega times I \Omega (MGyr is in the l.h.s of the equation of motion)
  DEBUG_BEGIN("NewtonEulerDS::computeMGyr(SP::SiconosVector twist, SP::SiconosVector mGyr)\n");

  ::computeMGyr_internal(_I, twist, mGyr);

  DEBUG_END("NewtonEulerDS::computeMGyr(SP::SiconosVector twist, SP::SiconosVector mGyr)\n");
}
void NewtonEulerDS::computeMGyr(SP::SiconosVector twist) {
  /*computation of \Omega times I \Omega*/
  // DEBUG_BEGIN("NewtonEulerDS::computeMGyr(SP::SiconosVector twist)\n");
  ::computeMGyr_internal(_I, twist, _mGyr);
  // DEBUG_END("NewtonEulerDS::computeMGyr(SP::SiconosVector twist)\n");
}

void NewtonEulerDS::computeForces(double time, SP::SiconosVector q, SP::SiconosVector twist) {
  DEBUG_BEGIN(
      "NewtonEulerDS::computeForces(double time, SP::SiconosVector q, SP::SiconosVector "
      "twist)\n")

  if (_wrench) {
    _wrench->zero();

    // External wrench

    if (_fExt) {
      computeFExt(time);
      assert(!std::isnan(_fExt->vector_sum()));
      _wrench->setBlock(0, *_fExt);
    }
    if (_mExt) {
      computeMExt(time);
      assert(!std::isnan(_mExt->vector_sum()));
      if (_isMextExpressedInInertialFrame) {
        SP::SiconosVector mExt(std::make_shared<SiconosVector>(*_mExt));
        ::changeFrameAbsToBody(q, mExt);
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
      computeMGyr(twist);
      assert(!std::isnan(_mGyr->vector_sum()));
      _wrench->setValue(3, _wrench->getValue(3) - _mGyr->getValue(0));
      _wrench->setValue(4, _wrench->getValue(4) - _mGyr->getValue(1));
      _wrench->setValue(5, _wrench->getValue(5) - _mGyr->getValue(2));
    }
    DEBUG_EXPR(_wrench->display());
    DEBUG_END(
        "NewtonEulerDS::computeForces(double time, SP::SiconosVector q, SP::SiconosVector "
        "twist)\n")

  } else {
    THROW_EXCEPTION("NewtonEulerDS::computeForces _wrench is null");
  }
  // else nothing.
}

void NewtonEulerDS::computeJacobianqForces(double time) {
  DEBUG_BEGIN("NewtonEulerDS::computeJacobianqWrench(double time) \n");
  if (_jacobianWrenchq) {
    _jacobianWrenchq->zero();
    if (_jacobianFIntq) {
      computeJacobianFIntq(time);
      _jacobianWrenchq->setBlock(0, 0, -1.0 * *_jacobianFIntq);
    }
    if (_jacobianMIntq) {
      computeJacobianMIntq(time);
      _jacobianWrenchq->setBlock(3, 0, -1.0 * *_jacobianMIntq);
    }
    if (_isMextExpressedInInertialFrame && _mExt) {
      computeJacobianMExtqExpressedInInertialFrame(time, _q);
      _jacobianWrenchq->setBlock(3, 0, 1.0 * *_jacobianMExtq);
    }
    DEBUG_EXPR(_jacobianWrenchq->display(););
  } else {
    THROW_EXCEPTION("NewtonEulerDS::computeJacobianqForces _jacobianWrenchq is null");
  }
  // else nothing.
  DEBUG_END("NewtonEulerDS::computeJacobianqForces(double time) \n");
}

void NewtonEulerDS::computeJacobianvForces(double time) {
  DEBUG_BEGIN("NewtonEulerDS::computeJacobiantwistForces(double time) \n");
  if (_jacobianWrenchTwist) {
    _jacobianWrenchTwist->zero();
    if (_jacobianFInttwist) {
      computeJacobianFIntv(time);
      _jacobianWrenchTwist->setBlock(0, 0, -1.0 * *_jacobianFInttwist);
    }
    if (_jacobianMInttwist) {
      computeJacobianMIntv(time);
      _jacobianWrenchTwist->setBlock(3, 0, -1.0 * *_jacobianMInttwist);
    }
    if (!_nullifyMGyr) {
      if (_jacobianMGyrtwist) {
        // computeJacobianMGyrtwistByFD(time,_q,_twist);
        computeJacobianMGyrtwist(time);
        _jacobianWrenchTwist->setBlock(3, 0, -1.0 * *_jacobianMGyrtwist);
      }
    }
  }
  // else nothing.
  DEBUG_END("NewtonEulerDS::computeJacobiantwistForces(double time) \n");
}

void NewtonEulerDS::computeJacobianMGyrtwist(double time) {
  DEBUG_BEGIN("NewtonEulerDS::computeJacobianMGyrtwist(double time) \n");
  if (_jacobianMGyrtwist) {
    // Omega /\ I \Omega:
    _jacobianMGyrtwist->zero();
    SiconosVector omega(3);
    omega.setValue(0, _twist->getValue(3));
    omega.setValue(1, _twist->getValue(4));
    omega.setValue(2, _twist->getValue(5));
    SiconosVector Iomega(3);
    prod(*_I, omega, Iomega, true);
    SiconosVector ei(3);
    SiconosVector Iei(3);
    SiconosVector ei_Iomega(3);
    SiconosVector omega_Iei(3);

    /*See equation of DevNotes.pdf, equation with label eq:NE_nablaFL1*/
    for (int i = 0; i < 3; i++) {
      ei.zero();
      ei.setValue(i, 1.0);
      prod(*_I, ei, Iei, true);
      cross_product(omega, Iei, omega_Iei);
      cross_product(ei, Iomega, ei_Iomega);
      for (int j = 0; j < 3; j++)
        _jacobianMGyrtwist->setValue(j, 3 + i, ei_Iomega.getValue(j) + omega_Iei.getValue(j));
    }
    // Check if Jacobian is valid. Warning to the transpose operation in
    // _jacobianMGyrtwist->setValue(3 + j, 3 + i, ei_Iomega.getValue(j) +
    // omega_Iei.getValue(j));
  }
  // else nothing.
  DEBUG_EXPR(_jacobianMGyrtwist->display());
  // _jacobianMGyrtwist->display();
  // SP::SimpleMatrix jacobianMGyrtmp (new SimpleMatrix(*_jacobianMGyrtwist));
  // computeJacobianMGyrtwistByFD(time, _q, _twist);
  // jacobianMGyrtmp->display();
  // std::cout << "#################  " << (*jacobianMGyrtmp - *_jacobianMGyrtwist).normInf()
  // << std::endl; assert((*jacobianMGyrtmp - *_jacobianMGyrtwist).normInf()< 1e-10);
  DEBUG_END("NewtonEulerDS::computeJacobianMGyrtwist(double time) \n");
}

void NewtonEulerDS::display(bool brief) const {
  std::cout << "=====> NewtonEuler System display (number: " << _number << ")." << std::endl;
  std::cout << "- _ndof : " << _ndof << std::endl;
  std::cout << "- _qDim : " << _qDim << std::endl;
  std::cout << "- _n : " << _n << std::endl;
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
  std::cout << "Inertia :" << std::endl;
  if (_I)
    _I->display();
  else
    std::cout << "-> nullptr" << std::endl;
  std::cout << "===================================== " << std::endl;
}

// --- Functions for memory handling ---
void NewtonEulerDS::initMemory(unsigned int steps) {
  DynamicalSystem::initMemory(steps);

  if (steps == 0)
    std::cout << "Warning : NewtonEulerDS::initMemory with size equal to zero" << std::endl;
  else {
    _qMemory.setMemorySize(steps, _qDim);
    _twistMemory.setMemorySize(steps, _ndof);
    _forcesMemory.setMemorySize(steps, _ndof);
    _dotqMemory.setMemorySize(steps, _qDim);
    //    swapInMemory(); Useless, done in osi->initializeWorkVectorsForDS
  }
}

void NewtonEulerDS::swapInMemory() {
  //  _xMemory->swap(_x[0]);
  _qMemory.swap(*_q);
  _twistMemory.swap(*_twist);
  _dotqMemory.swap(*_dotq);
  _forcesMemory.swap(*_wrench);
}

void NewtonEulerDS::resetAllNonSmoothParts() {
  if (_p[1])
    _p[1]->zero();
  else
    _p[1].reset(new SiconosVector(_ndof));
}
void NewtonEulerDS::resetNonSmoothPart(unsigned int level) {
  if (_p[level]) _p[level]->zero();
}

void NewtonEulerDS::computeT() { ::computeT(_q, _T); }

void NewtonEulerDS::computeTdot() {
  if (!_Tdot) {
    _Tdot.reset(new SimpleMatrix(_qDim, _ndof));
    _Tdot->zero();
  }

  ::computeT(_dotq, _Tdot);
}

void NewtonEulerDS::normalizeq() { ::normalizeq(_q); }

void NewtonEulerDS::setComputeJacobianFIntqFunction(const std::string& pluginPath,
                                                    const std::string& functionName) {
  //    Plugin::setFunction(&computeJacobianFIntqPtr, pluginPath,functionName);
  _pluginJacqFInt->setComputeFunction(pluginPath, functionName);
  if (!_jacobianFIntq) _jacobianFIntq.reset(new SimpleMatrix(3, _qDim));
  if (!_jacobianWrenchq) _jacobianWrenchq.reset(new SimpleMatrix(_ndof, _qDim));
  _computeJacobianFIntqByFD = false;
}
void NewtonEulerDS::setComputeJacobianFIntvFunction(const std::string& pluginPath,
                                                    const std::string& functionName) {
  //    Plugin::setFunction(&computeJacobianFIntvPtr, pluginPath,functionName);
  _pluginJactwistFInt->setComputeFunction(pluginPath, functionName);
  if (!_jacobianFInttwist) _jacobianFInttwist.reset(new SimpleMatrix(3, _ndof));
  _computeJacobianFInttwistByFD = false;
}
void NewtonEulerDS::setComputeJacobianFIntqFunction(FInt_NE fct) {
  _pluginJacqFInt->setComputeFunction((void*)fct);
  if (!_jacobianFIntq) _jacobianFIntq.reset(new SimpleMatrix(3, _qDim));
  if (!_jacobianWrenchq) _jacobianWrenchq.reset(new SimpleMatrix(_ndof, _qDim));
  _computeJacobianFIntqByFD = false;
}
void NewtonEulerDS::setComputeJacobianFIntvFunction(FInt_NE fct) {
  _pluginJactwistFInt->setComputeFunction((void*)fct);
  if (!_jacobianFInttwist) _jacobianFInttwist.reset(new SimpleMatrix(3, _ndof));
  if (!_jacobianWrenchTwist) _jacobianWrenchTwist.reset(new SimpleMatrix(_ndof, _ndof));
  _computeJacobianFInttwistByFD = false;
}

void NewtonEulerDS::setComputeJacobianMIntqFunction(const std::string& pluginPath,
                                                    const std::string& functionName) {
  _pluginJacqMInt->setComputeFunction(pluginPath, functionName);
  if (!_jacobianMIntq) _jacobianMIntq.reset(new SimpleMatrix(3, _qDim));
  if (!_jacobianWrenchq) _jacobianWrenchq.reset(new SimpleMatrix(_ndof, _qDim));
  _computeJacobianMIntqByFD = false;
}
void NewtonEulerDS::setComputeJacobianMIntvFunction(const std::string& pluginPath,
                                                    const std::string& functionName) {
  _pluginJactwistMInt->setComputeFunction(pluginPath, functionName);
  if (!_jacobianMInttwist) _jacobianMInttwist.reset(new SimpleMatrix(3, _ndof));
  if (!_jacobianWrenchTwist) _jacobianWrenchTwist.reset(new SimpleMatrix(_ndof, _ndof));
  _computeJacobianMInttwistByFD = false;
}
void NewtonEulerDS::setComputeJacobianMIntqFunction(FInt_NE fct) {
  _pluginJacqMInt->setComputeFunction((void*)fct);
  if (!_jacobianMIntq) _jacobianMIntq.reset(new SimpleMatrix(3, _qDim));
  if (!_jacobianWrenchq) _jacobianWrenchq.reset(new SimpleMatrix(_ndof, _qDim));
  _computeJacobianMIntqByFD = false;
}
void NewtonEulerDS::setComputeJacobianMIntvFunction(FInt_NE fct) {
  _pluginJactwistMInt->setComputeFunction((void*)fct);
  if (!_jacobianMInttwist) _jacobianMInttwist.reset(new SimpleMatrix(3, _ndof));
  if (!_jacobianWrenchTwist) _jacobianWrenchTwist.reset(new SimpleMatrix(_ndof, _ndof));
  _computeJacobianMInttwistByFD = false;
}

double NewtonEulerDS::computeKineticEnergy() {
  DEBUG_BEGIN("NewtonEulerDS::computeKineticEnergy()\n");
  assert(_twist);
  assert(_mass);
  DEBUG_EXPR(_twist->display());
  DEBUG_EXPR(_mass->display());

  SiconosVector tmp(6);
  prod(*_mass, *_twist, tmp, true);
  double K = 0.5 * inner_prod(tmp, *_twist);

  DEBUG_PRINTF("Kinetic Energy = %e\n", K);
  DEBUG_END("NewtonEulerDS::computeKineticEnergy()\n");
  return K;
}

SP::SiconosVector NewtonEulerDS::linearVelocity(bool absoluteRef) const {
  // Short-cut: return the _twist 6-vector without modification, first
  // 3 components are the expected linear velocity.
  if (absoluteRef) return _twist;

  SP::SiconosVector v(std::make_shared<SiconosVector>(3));
  linearVelocity(absoluteRef, *v);
  return v;
}

void NewtonEulerDS::linearVelocity(bool absoluteRef, SiconosVector& v) const {
  v(0) = (*_twist)(0);
  v(1) = (*_twist)(1);
  v(2) = (*_twist)(2);

  /* See _twist: linear velocity is in absolute frame */
  if (!absoluteRef) changeFrameAbsToBody(*_q, v);
}

SP::SiconosVector NewtonEulerDS::angularVelocity(bool absoluteRef) const {
  SP::SiconosVector w(std::make_shared<SiconosVector>(3));
  angularVelocity(absoluteRef, *w);
  return w;
}

void NewtonEulerDS::angularVelocity(bool absoluteRef, SiconosVector& w) const {
  w(0) = (*_twist)(3);
  w(1) = (*_twist)(4);
  w(2) = (*_twist)(5);

  /* See _twist: angular velocity is in relative frame */
  if (absoluteRef) changeFrameBodyToAbs(*_q, w);
}

void computeExtForceAtPos(SP::SiconosVector q, bool isMextExpressedInInertialFrame,
                          SP::SiconosVector force, bool forceAbsRef, SP::SiconosVector pos,
                          bool posAbsRef, SP::SiconosVector fExt, SP::SiconosVector mExt,
                          bool accumulate) {
  assert(!!fExt && fExt->size() == 3);
  assert(!!force && force->size() == 3);
  if (pos) assert(!!mExt && mExt->size() == 3);

  SiconosVector abs_frc(*force), local_frc(*force);

  if (forceAbsRef) {
    if (pos) changeFrameAbsToBody(*q, local_frc);
  } else
    changeFrameBodyToAbs(*q, abs_frc);

  if (pos) {
    assert(!!mExt && mExt->size() >= 3);
    SiconosVector moment(3);
    if (posAbsRef) {
      SiconosVector local_pos(*pos);
      local_pos(0) -= (*q)(0);
      local_pos(1) -= (*q)(1);
      local_pos(2) -= (*q)(2);
      changeFrameAbsToBody(*q, local_pos);
      cross_product(local_pos, local_frc, moment);
    } else {
      cross_product(*pos, local_frc, moment);
    }

    if (isMextExpressedInInertialFrame) changeFrameBodyToAbs(*q, moment);

    if (accumulate)
      *mExt = *mExt + moment;
    else
      *mExt = moment;
  }

  if (accumulate)
    *fExt += *fExt + abs_frc;
  else
    *fExt = abs_frc;
}

void NewtonEulerDS::addExtForceAtPos(SP::SiconosVector force, bool forceAbsRef,
                                     SP::SiconosVector pos, bool posAbsRef) {
  assert(!!_fExt && _fExt->size() == 3);
  assert(!!force && force->size() == 3);
  if (pos) assert(!!_mExt && _mExt->size() == 3);

  computeExtForceAtPos(_q, _isMextExpressedInInertialFrame, force, forceAbsRef, pos, posAbsRef,
                       _fExt, _mExt, true);
}
