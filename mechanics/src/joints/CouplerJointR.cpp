/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
/*! \file CouplerJointRR.cpp
 */

#include "CouplerJointR.hpp"

#include <Interaction.hpp>

#include "BlockVector.hpp"
#include "NewtonEulerDS.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
// #include <boost/math/quaternion.hpp>
// #include <cfloat>
// #include <iostream>

// #define DEBUG_BEGIN_END_ONLY
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

/** Initialize a coupler joint. For use with EqualityConditionNSL to
 * bind two degrees of freedom into one by a ratio and offset. */
siconos::joints::CouplerJointR::CouplerJointR(
    std::shared_ptr<NewtonEulerJointR> joint1, siconos::algebra::Index dof1,
    std::shared_ptr<NewtonEulerJointR> joint2, siconos::algebra::Index dof2, double ratio,
    std::optional<siconos::algebra::SiconosVector> ref1, siconos::algebra::Index ref1_index,
    std::optional<siconos::algebra::SiconosVector> ref2, siconos::algebra::Index ref2_index)
    : NewtonEulerJointR{},
      _joint1(joint1),
      _joint2(joint2),
      // _ref1(ref1),
      // _ref2(ref2),
      _dof1(dof1),
      _dof2(dof2),
      _ref1_index(ref1_index),
      _ref2_index(ref2_index),
      _ratio(ratio),
      _offset(0.0) {
  if (ref1)
    _ref1 = std::make_shared<siconos::algebra::SiconosVector>(*ref1);  // ref1.value()
  else
    _ref1 = nullptr;

  if (ref2)
    _ref2 = std::make_shared<siconos::algebra::SiconosVector>(*ref2);
  else
    _ref2 = nullptr;
  assert(_dof1 < _joint1->numberOfDoF());
  assert(_dof2 < _joint2->numberOfDoF());
}

/* A constructor taking a DS exists just because it's hard to pass
 * ds->q() through Python without it automatically converting to numpy
 * array and back, which messes up the shared_ptr reference.
 *
 * NOTE that using q() as the reference is not quite right, in fact it
 * should be using the reference ds's temporary work vector in order
 * to perform correctly in the Newton loop. (TODO) */
siconos::joints::CouplerJointR::CouplerJointR(
    std::shared_ptr<NewtonEulerJointR> joint1, siconos::algebra::Index dof1,
    std::shared_ptr<NewtonEulerJointR> joint2, siconos::algebra::Index dof2, double ratio,
    std::shared_ptr<siconos::modeling::NewtonEulerDS> refds1,
    siconos::algebra::Index ref1_index,
    std::shared_ptr<siconos::modeling::NewtonEulerDS> refds2,
    siconos::algebra::Index ref2_index)
    : NewtonEulerJointR(),
      _joint1(joint1),
      _joint2(joint2),
      _ref1(refds1->q()),
      _ref2(refds2->q()),
      _dof1(dof1),
      _dof2(dof2),
      _ref1_index(ref1_index),
      _ref2_index(ref2_index),
      _ratio(ratio) {
  assert(_dof1 < _joint1->numberOfDoF());
  assert(_dof2 < _joint2->numberOfDoF());
}

void siconos::joints::CouplerJointR::setReferences(
    std::shared_ptr<NewtonEulerJointR> joint1, siconos::algebra::Index dof1,
    std::shared_ptr<NewtonEulerJointR> joint2, siconos::algebra::Index dof2,
    std::shared_ptr<siconos::algebra::SiconosVector> ref1, siconos::algebra::Index ref1_index,
    std::shared_ptr<siconos::algebra::SiconosVector> ref2,
    siconos::algebra::Index ref2_index) {
  _joint1 = joint1;
  _joint2 = joint2;
  _dof1 = dof1;
  _dof2 = dof2;
  _ref1 = ref1;
  _ref2 = ref2;
  _ref1_index = ref1_index;
  _ref2_index = ref2_index;
}

void siconos::joints::CouplerJointR::setReferences(
    std::shared_ptr<NewtonEulerJointR> joint1, siconos::algebra::Index dof1,
    std::shared_ptr<NewtonEulerJointR> joint2, siconos::algebra::Index dof2,
    std::shared_ptr<siconos::modeling::NewtonEulerDS> refds1,
    siconos::algebra::Index ref1_index,
    std::shared_ptr<siconos::modeling::NewtonEulerDS> refds2,
    siconos::algebra::Index ref2_index) {
  _joint1 = joint1;
  _joint2 = joint2;
  _dof1 = dof1;
  _dof2 = dof2;
  _ref1 = refds1->q();
  _ref2 = refds2->q();
  _ref1_index = ref1_index;
  _ref2_index = ref2_index;
}

void siconos::joints::CouplerJointR::setRatio(double ratio) { _ratio = ratio; }

void siconos::joints::CouplerJointR::resolveVectors(
    const siconos::algebra::SiconosVector* q1, const siconos::algebra::SiconosVector* q2,
    const siconos::algebra::SiconosVector*& v1_1, const siconos::algebra::SiconosVector*& v1_2,
    const siconos::algebra::SiconosVector*& v2_1,
    const siconos::algebra::SiconosVector*& v2_2) const {
  // Case 1: no references, q2 present → coupling q1/q2
  if (!_ref1 && !_ref2 && q2) {
    v1_1 = q1;
    v1_2 = q2;
    v2_1 = q1;
    v2_2 = q2;
  }
  // Case 2: no q2, maybe ref1
  else if (!q2) {
    if (_ref1) {
      v1_1 = (_ref1_index == 0) ? _ref1.get() : q1;
      v1_2 = (_ref1_index == 0) ? q1 : _ref1.get();
      v2_1 = v1_1;
      v2_2 = v1_2;
    } else {
      v1_1 = q1;
      v1_2 = nullptr;
      v2_1 = q1;
      v2_2 = nullptr;
    }
  }
  // Case 3: q2 present, potentially ref1 and ref2
  else {
    // v1
    if (_ref1) {
      v1_1 = (_ref1_index == 0) ? _ref1.get() : q1;
      v1_2 = (_ref1_index == 0) ? q1 : _ref1.get();
    } else {
      v1_1 = q1;
      v1_2 = q2;
    }

    // v2
    if (_ref2) {
      v2_1 = (_ref2_index == 0) ? _ref2.get() : q2;
      v2_2 = (_ref2_index == 0) ? q2 : _ref2.get();
    } else {
      v2_1 = q1;
      v2_2 = q2;
    }
  }
}

void siconos::joints::CouplerJointR::setBasePositions(
    const Eigen::Ref<const siconos::algebra::SiconosVector>& q1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>& q2) {
  // TODO: time=0?

  // Get current positions of the implicated degrees of freedom
  siconos::algebra::SiconosVector y1(1), y2(1);

  // Make local copies to ensure stable pointers for resolveVectors
  // (Eigen::Ref parameters can create temporaries)
  siconos::algebra::SiconosVector q1_copy(q1);
  std::optional<siconos::algebra::SiconosVector> q2_copy;
  if (q2) {
    q2_copy = siconos::algebra::SiconosVector(*q2);
  }

  const siconos::algebra::SiconosVector* v1_1 = nullptr;
  const siconos::algebra::SiconosVector* v1_2 = nullptr;
  const siconos::algebra::SiconosVector* v2_1 = nullptr;
  const siconos::algebra::SiconosVector* v2_2 = nullptr;
  resolveVectors(&q1_copy, q2_copy ? &(*q2_copy) : nullptr, v1_1, v1_2, v2_1, v2_2);
  if (v1_2)
    _joint1->computehDoF(*v1_1, *v1_2, y1, _dof1);
  else
    _joint1->computehDoF(*v1_1, std::nullopt, y1, _dof1);
  if (v2_2)
    _joint2->computehDoF(*v2_1, *v2_2, y2, _dof2);
  else
    _joint2->computehDoF(*v2_1, std::nullopt, y2, _dof2);

  // Compute initial offset between the DoFs
  _offset = y1(0) * _ratio - y2(0);
}

void siconos::joints::CouplerJointR::computeh(
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>& q2,
    Eigen::Ref<siconos::algebra::SiconosVector> y) {
  siconos::algebra::SiconosVector y1(y), y2(y);

  // Make local copies to ensure stable pointers for resolveVectors
  // TEMP TO BE REVIEWED
  siconos::algebra::SiconosVector q1_copy(q1);
  std::optional<siconos::algebra::SiconosVector> q2_copy;
  if (q2) {
    q2_copy = siconos::algebra::SiconosVector(*q2);
  }

  // Get current positions of the implicated degrees of freedom
  const siconos::algebra::SiconosVector* v1_1 = nullptr;
  const siconos::algebra::SiconosVector* v1_2 = nullptr;
  const siconos::algebra::SiconosVector* v2_1 = nullptr;
  const siconos::algebra::SiconosVector* v2_2 = nullptr;

  if (q2) {
    resolveVectors(&q1_copy, &(*q2_copy), v1_1, v1_2, v2_1, v2_2);
  } else {
    resolveVectors(&q1_copy, nullptr, v1_1, v1_2, v2_1, v2_2);
  }

  // Compute hDoF for both joints
  if (v1_2)
    _joint1->computehDoF(*v1_1, *v1_2, y1, _dof1);
  else
    _joint1->computehDoF(*v1_1, std::nullopt, y1, _dof1);
  if (v2_2)
    _joint2->computehDoF(*v2_1, *v2_2, y2, _dof2);
  else
    _joint2->computehDoF(*v2_1, std::nullopt, y2, _dof2);

  // Constraint is the linear relation between them
  y(0) = y2(0) - y1(0) * _ratio + _offset;
}

void siconos::joints::CouplerJointR::computeH_NE_(double time,
                                                  siconos::modeling::Interaction& inter,
                                                  const siconos::algebra::BlockVector& q0) {
  auto jachq1 = std::make_shared<siconos::algebra::SiconosMatrix>(1, q0.size());
  auto jachq2 = std::make_shared<siconos::algebra::SiconosMatrix>(1, q0.size());

  // Get jacobians for the implicated degrees of freedom
  const siconos::algebra::SiconosVector* v1_1 = nullptr;
  const siconos::algebra::SiconosVector* v1_2 = nullptr;
  const siconos::algebra::SiconosVector* v2_1 = nullptr;
  const siconos::algebra::SiconosVector* v2_2 = nullptr;

  resolveVectors(q0.vector(0).get(), q0.numberOfBlocks() > 1 ? q0.vector(1).get() : nullptr,
                 v1_1, v1_2, v2_1, v2_2);

  if (v1_2) {
    _joint1->computeJachqDoF(inter, *v1_1, *v1_2, *jachq1, _dof1);
  } else {
    _joint1->computeJachqDoF(inter, *v1_1, std::nullopt, *jachq1, _dof1);
  }

  if (v2_2) {
    _joint2->computeJachqDoF(inter, *v2_1, *v2_2, *jachq2, _dof2);
  } else {
    _joint2->computeJachqDoF(inter, *v2_1, std::nullopt, *jachq2, _dof2);
  }

  // Constraint is the linear relation between them
  for (siconos::algebra::Index i = 0; i < 1; i++)
    for (siconos::algebra::Index j = 0; j < H_NE_view_->cols(); j++)
      H_NE_view_->setValue(i, j, (*jachq2)(i, j) - (*jachq1)(i, j) * _ratio);
}

// void siconos::joints::CouplerJointR::computeH_NE_(double time,
//                                                   siconos::modeling::Interaction& inter,
//                                                   const siconos::algebra::BlockVector& q0) {
//   const auto& q1 = *q0.vector(0);
//   std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>> q2 =
//       (q0.numberOfBlocks() > 1)
//           ? std::make_optional(
//                 Eigen::Ref<const siconos::algebra::SiconosVector>(*q0.vector(1)))
//           : std::nullopt;

//   const siconos::algebra::SiconosVector* v1_1 = nullptr;
//   const siconos::algebra::SiconosVector* v1_2 = nullptr;
//   const siconos::algebra::SiconosVector* v2_1 = nullptr;
//   const siconos::algebra::SiconosVector* v2_2 = nullptr;

//   resolveVectors(q1, q2, v1_1, v1_2, v2_1, v2_2);
//   auto jachq1 = std::make_shared<siconos::algebra::SiconosMatrix>(1, q0.size());
//   auto jachq2 = std::make_shared<siconos::algebra::SiconosMatrix>(1, q0.size());

//   // Get jacobians for the implicated degrees of freedom
//   // Compute the jacobian for the required range of axes
//   _joint1->computeJachqDoF(
//       inter, *v1_1,
//       v1_2 ? std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>(*v1_2)
//            : std::nullopt,
//       *jachq1, _dof1);

//   _joint2->computeJachqDoF(
//       inter, *v2_1,
//       v2_2 ? std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>(*v2_2)
//            : std::nullopt,
//       *jachq2, _dof2);

//   // Constraint is the linear relation between them
//   for (unsigned int i = 0; i < 1; i++)
//     for (unsigned int j = 0; j < H_NE_view_->cols(); j++)
//       H_NE_view_->setValue(i, j, (*jachq2)(i, j) - (*jachq1)(i, j) * _ratio);
// }

// siconos::algebra::SiconosVector siconos::joints::CouplerJointR::normalDoF(
//     const siconos::algebra::BlockVector& q0, int axis, bool absoluteRef) {
//   assert("siconos::joints::CouplerJointR::_normalDoF is not defined.");

//   // We could define the normal of this constraint as simply the cross
//   // product of the normals of the two coupled DoFs, as calculated
//   // below, however this is not really sensible as the constraint is
//   // defined in terms of a line through a manifold, and would need to
//   // be projected down to world coordinates somehow here.  For
//   // example, what is the "direction of the constraint" in world
//   // coordinates for a coupling between rotational and lineawr motion?

//   /*
//   siconos::algebra::SiconosVector n1(3), n2(3);
//   siconos::algebra::BlockVector q01, q02;
//   makeBlockVectors(q0.getAllVect()[0],
//                    q0.getAllVect().size()>1 ? q0.getAllVect()[1] :
//   std::shared_ptr<siconos::algebra::SiconosVector>(), q01, q02); _joint1->normalDoF(n1,
//   q01, _dof1, absoluteRef); _joint2->normalDoF(n2, q02, _dof2, absoluteRef);
//   cross_product(n1, n2, ans);
//   */
// }

void siconos::joints::CouplerJointR::computehDoF(
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>& q2,
    Eigen::Ref<siconos::algebra::SiconosVector> y, unsigned int axis) {
  // The DoF of the constraint is the same as the constraint itself
  assert(axis == 0);
  computeh(q1, q2, y);
}

void siconos::joints::CouplerJointR::computeJachqDoF(
    siconos::modeling::Interaction& inter,
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>& q2,
    Eigen::Ref<siconos::algebra::SiconosMatrix> jachq, unsigned int axis) {
  // The Jacobian of the DoF of the constraint is the same as the
  // Jacobian of the constraint itself. (Same as computeJacobianhOver_q(), but
  // don't store result in member object.)
  assert(axis == 0);

  // Make local copies to ensure stable pointers for resolveVectors
  siconos::algebra::SiconosVector q1_copy(q1);
  std::optional<siconos::algebra::SiconosVector> q2_copy;
  if (q2) {
    q2_copy = siconos::algebra::SiconosVector(*q2);
  }

  const siconos::algebra::SiconosVector* v1_1 = nullptr;
  const siconos::algebra::SiconosVector* v1_2 = nullptr;
  const siconos::algebra::SiconosVector* v2_1 = nullptr;
  const siconos::algebra::SiconosVector* v2_2 = nullptr;
  resolveVectors(&q1_copy, q2_copy ? &(*q2_copy) : nullptr, v1_1, v1_2, v2_1, v2_2);

  auto jachq1 = std::make_shared<siconos::algebra::SiconosMatrix>(
      1, q1_copy.size() + (q2_copy ? q2_copy->size() : 0));
  auto jachq2 = std::make_shared<siconos::algebra::SiconosMatrix>(
      1, q1_copy.size() + (q2_copy ? q2_copy->size() : 0));

  if (v1_2) {
    _joint1->computeJachqDoF(inter, *v1_1, *v1_2, *jachq1, _dof1);
  } else {
    _joint1->computeJachqDoF(inter, *v1_1, std::nullopt, *jachq1, _dof1);
  }

  if (v2_2) {
    _joint2->computeJachqDoF(inter, *v2_1, *v2_2, *jachq2, _dof2);
  } else {
    _joint2->computeJachqDoF(inter, *v2_1, std::nullopt, *jachq2, _dof2);
  }

  // Constraint is the linear relation between them
  for (siconos::algebra::Index i = 0; i < 1; i++)
    for (siconos::algebra::Index j = 0; j < jachq.cols(); j++)
      jachq(i, j) = (*jachq2)(i, j) - (*jachq1)(i, j) * _ratio;
}
