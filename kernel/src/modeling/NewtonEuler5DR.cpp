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

#include "NewtonEuler5DR.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "RotationQuaternion.hpp"  // siconos::geometry::computeRotationMatrix
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

/*
See devNotes.pdf for details. A detailed documentation is available in DevNotes.pdf: chapter
'NewtonEulerR: computation of \nabla q H'. Subsection 'Case RFC3D: using the local frame local
velocities'
*/
void siconos::modeling::NewtonEuler5DR::initialize(Interaction& inter) {
  DEBUG_BEGIN(
      "siconos::modeling::NewtonEuler5DR::siconos::modeling::NewtonEuler5DR::initialize("
      "Interaction& inter)\n");
  auto qSize = 7 * (inter.getSizeOfDS() / 6);
  /*keep only the distance.*/

  H_NE_internal_storage_.resize(5, qSize);
  H_NE_view_ =
      std::make_shared<siconos::algebra::MapType>(H_NE_internal_storage_.data(), 5, qSize);
  H_NE_view_->setZero();
  NewtonEulerR::initialize(inter);

  //  _isContact=1;
  DEBUG_END(
      "siconos::modeling::NewtonEuler5DR::siconos::modeling::NewtonEuler5DR::initialize("
      "Interaction& inter)\n");
}

void siconos::modeling::NewtonEuler5DR::RFC3DcomputeJachqTFromContacts(
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1) {
  DEBUG_BEGIN("siconos::modeling::NewtonEuler5DR::RFC3DcomputeJachqTFromContacts()\n");
  DEBUG_PRINT("contact normal:\n");
  DEBUG_EXPR(siconos::algebra::print(nc_););
  DEBUG_PRINTF("nc_.norm() -1.0 = %e\n", nc_.norm() - 1.0);
  DEBUG_PRINT("contact point :\n");
  DEBUG_EXPR(siconos::algebra::print(contactPoint1_););
  DEBUG_PRINT("center of mass :\n");
  DEBUG_EXPR(siconos::algebra::print(q1););

  assert(nc_.norm() > 0.0 && std::abs(nc_.norm() - 1.0) < 1e-6 &&
         "siconos::modeling::NewtonEuler5DR::RFC3DcomputeJachqTFromContacts. Normal vector "
         "not consistent ");

  // 1 - Construction of the local contact frame from the normal vector
  siconos::algebra::SiconosVector3 t1, t2;
  bool res = siconos::geometry::orthoBaseFromVector(nc_, t1, t2);
  if (!res) {
    THROW_EXCEPTION(
        "siconos::modeling::NewtonEuler5DR::RFC3DcomputeJachqTFromContacts. Problem in "
        "calling orthoBaseFromVector");
  }
  // 2 - Construction of the rotation matrix from the absolute frame to the local contact frame
  rotationAbsoluteToContactFrame_.row(0) = nc_;
  rotationAbsoluteToContactFrame_.row(1) = t1;
  rotationAbsoluteToContactFrame_.row(2) = t2;
  DEBUG_PRINT("_rotationAbsoluteToContactFrame:\n");
  DEBUG_EXPR(siconos::algebra::print(*_rotationAbsoluteToContactFrame););

  // 3 - Construction of the lever arm matrix in  the absolute frame
  const auto v = q1.head<3>() - contactPoint1_.head<3>();
  NPG_buffer_ << 0.0, -v.z(), v.y(), v.z(), 0.0, -v.x(), -v.y(), v.x(), 0.0;

  DEBUG_PRINT("lever arm skew matrix :\n");
  DEBUG_EXPR(siconos::algebra::print(NPG_buffer_););

  /* The Jacobian matrix (H) is given by the product
   * H = _rotationAbsoluteToContactFrame
   * for the translation part and
   * H = _rotationAbsoluteToContactFrame * leverArmMatrix * _rotationBodyToAbsoluteFrame
   * for the rotation part and
   */

  // 4 - Compute the rotation matrix from the body-fixed frame to the absolute frame
  siconos::geometry::computeRotationMatrix(q1, rotationBodyToAbsoluteFrame_);
  DEBUG_EXPR(siconos::algebra::print(*_rotationBodyToAbsoluteFrame););

  // 5 - compose the body lever arm matrix with the rotation matrix, rotate the resulting
  // matrix in the contact frame and fill the Jacobian

  H_NE_prod_T_->block(0, 0, 3, 3) = rotationAbsoluteToContactFrame_;
  H_NE_prod_T_->block(0, 3, 3, 3) =
      rotationAbsoluteToContactFrame_ * NPG_buffer_ * rotationBodyToAbsoluteFrame_;
  H_NE_prod_T_->block(3, 3, 2, 3) =
      (rotationAbsoluteToContactFrame_ * rotationBodyToAbsoluteFrame_).bottomRows(2);

  DEBUG_EXPR(siconos::algebra::print(*jacobianhOver_q_T););
  DEBUG_END(
      "siconos::modeling::NewtonEuler5DR::RFC3DcomputeJachqTFromContacts(std::shared_ptr<"
      "siconos::algebra::SiconosVector> q1)\n");
}

void siconos::modeling::NewtonEuler5DR::RFC3DcomputeJachqTFromContacts(
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q2) {
  DEBUG_BEGIN("siconos::modeling::NewtonEuler5DR::RFC3DcomputeJachqTFromContacts()\n");
  RFC3DcomputeJachqTFromContacts(q1);
  const auto v = q2.head<3>() - contactPoint1_.head<3>();
  NPG_buffer_ << 0.0, -v.z(), v.y(), v.z(), 0.0, -v.x(), -v.y(), v.x(), 0.0;

  siconos::geometry::computeRotationMatrix(q2, rotationBodyToAbsoluteFrame_);
  H_NE_prod_T_->block(0, 6, 3, 3) = -rotationAbsoluteToContactFrame_;
  H_NE_prod_T_->block(0, 9, 3, 3) =
      rotationAbsoluteToContactFrame_ * NPG_buffer_ * rotationBodyToAbsoluteFrame_;
  H_NE_prod_T_->block(3, 3, 2, 3) =
      -(rotationAbsoluteToContactFrame_ * rotationBodyToAbsoluteFrame_).bottomRows(2);

  DEBUG_EXPR(siconos::algebra::print(*jacobianhOver_q_T););

  DEBUG_END(
      "siconos::modeling::NewtonEuler5DR::RFC3DcomputeJachqTFromContacts(std::shared_ptr<"
      "siconos::algebra::SiconosVector> q1, std::shared_ptr<siconos::algebra::SiconosVector> "
      "q2)\n");
}

void siconos::modeling::NewtonEuler5DR::computeH_NE_prod_T(
    const Interaction& inter, const siconos::algebra::BlockVector& q0) {
  DEBUG_BEGIN(
      "siconos::modeling::NewtonEuler5DR::computeH_NE_prod_T(Interaction& inter,  "
      "std::shared_ptr<siconos::algebra::BlockVector> q0)\n");
  if (q0.numberOfBlocks() > 1) {
    RFC3DcomputeJachqTFromContacts(*q0.vector(0), *q0.vector(1));
  } else {
    RFC3DcomputeJachqTFromContacts(*q0.vector(0));
  }
  DEBUG_END(
      "siconos::modeling::NewtonEuler5DR::computeH_NE_prod_T(Interaction& inter,  "
      "std::shared_ptr<siconos::algebra::BlockVector> q0)\n");
}
