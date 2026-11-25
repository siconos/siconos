// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "RigidBodyDS.hpp"

#include "RotationQuaternion.hpp"  // for quaternionFromTwistVector and compositionLawLieGroup
#include "SiconosContactor.hpp"
#include "SiconosVector.hpp"
#include "StorageTools.hpp"

// siconos::collision::RigidBodyDS::RigidBodyDS(
//     Eigen::Ref<siconos::algebra::SiconosVector7> position,
//     Eigen::Ref<siconos::algebra::SiconosVector6> twist, double mass,
//     Eigen::Ref<siconos::algebra::SiconosMatrix33> inertia, siconos::algebra::AliasTag)
//     : siconos::modeling::NewtonEulerDS{position, twist, mass, inertia,
//                                        siconos::algebra::alias_t},
//       contactors_(std::make_shared<siconos::collision::SiconosContactorSet>()) {}

siconos::collision::RigidBodyDS::RigidBodyDS(const siconos::algebra::SiconosVector7& position,
                                             const siconos::algebra::SiconosVector6& twist,
                                             double mass,
                                             const siconos::algebra::SiconosMatrix33& inertia)
    : siconos::modeling::NewtonEulerDS{position, twist, mass, inertia,
                                       siconos::algebra::copy_t},
      contactors_(std::make_shared<siconos::collision::SiconosContactorSet>()) {}

void siconos::collision::RigidBodyDS::compute_extrapolated_position(
    double extrapolationCoefficient) {
  // we compute an extrapolation of the position
  if (!_qExtrapolated) {
    _qExtrapolated = std::make_shared<siconos::algebra::SiconosVector>(state_q_->size());
  }

  // auto velocityIncrement =
  // std::make_shared<siconos::algebra::SiconosVector>(twist_->size()); all the following is
  // useless since velocityIncrement = 0 ...
  //     _qExtrapolated->setZero() will do the job
  // (*_qExtrapolated)(0) = (*velocityIncrement)(0);
  // (*_qExtrapolated)(1) = (*velocityIncrement)(1);
  // (*_qExtrapolated)(2) = (*velocityIncrement)(2);
  // siconos::geometry::quaternionFromTwistVector(*velocityIncrement, *_qExtrapolated);
  _qExtrapolated->setZero();
  const siconos::algebra::SiconosVector7& qold = qMemory().getSiconosVector(0);
  siconos::geometry::compositionLawLieGroup(qold, *_qExtrapolated);
}
