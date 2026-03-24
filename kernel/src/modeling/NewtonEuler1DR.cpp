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

#include "NewtonEuler1DR.hpp"

#include <boost/math/quaternion.hpp>

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "RotationQuaternion.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

// #define NERI_DEBUG

// #define NEFC3D_DEBUG
//  #define DEBUG_NOCOLOR
//  #define DEBUG_STDOUT
//  #define DEBUG_MESSAGES
#include "siconos_debug.h"
/*
See devNotes.pdf for details. A detailed documentation is available in DevNotes.pdf: chapter
'NewtonEulerR: computation of \nabla q H'. Subsection 'Case FC3D: using the local frame local
velocities'
*/
void siconos::modeling::NewtonEuler1DR::NIcomputeJachqTFromContacts(
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1) {
#ifdef NEFC3D_DEBUG
  printf("contact normal:\n");
  siconos::algebra::print(nc_);
  printf("point de contact :\n");
  siconos::algebra::print(contactPoint1_);
  printf("center of masse :\n");
  siconos::algebra::print(q1);
#endif
  rotationAbsoluteToContactFrame_ << nc_.x(), nc_.y(), nc_.z();

  const auto v = q1.head<3>() - contactPoint1_.head<3>();
  NPG_buffer_ << 0.0, -v.z(), v.y(), v.z(), 0.0, -v.x(), -v.y(), v.x(), 0.0;

  siconos::geometry::computeRotationMatrix(q1, rotationBodyToAbsoluteFrame_);
  H_NE_prod_T_->row(0).segment(0, 3) = rotationAbsoluteToContactFrame_;
  H_NE_prod_T_->row(0).segment(3, 3) =
      rotationAbsoluteToContactFrame_ * NPG_buffer_ * rotationBodyToAbsoluteFrame_;

#ifdef NEFC3D_DEBUG
  printf("NewtonEuler1DR jhqt\n");
  siconos::algebra::print(*jacobianhOver_q_T);
#endif
}

void siconos::modeling::NewtonEuler1DR::NIcomputeJachqTFromContacts(
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q2) {
  NIcomputeJachqTFromContacts(q1);

  const auto v = q2.head<3>() - contactPoint1_.head<3>();
  NPG_buffer_ << 0.0, -v.z(), v.y(), v.z(), 0.0, -v.x(), -v.y(), v.x(), 0.0;

  siconos::geometry::computeRotationMatrix(q2, rotationBodyToAbsoluteFrame_);

  H_NE_prod_T_->row(0).segment(6, 3) = rotationAbsoluteToContactFrame_;
  H_NE_prod_T_->row(0).segment(9, 3) =
      -rotationAbsoluteToContactFrame_ * NPG_buffer_ * rotationBodyToAbsoluteFrame_;
}

void siconos::modeling::NewtonEuler1DR::initialize(Interaction& inter) {
  // proj_with_q  jacobianhOver_q_Proj =
  // std::make_shared<siconos::algebra::SiconosMatrix>(jacobianhOver_q_->rows(),jacobianhOver_q_->cols()));
  auto qSize = 7;
  if (inter.has2Bodies()) qSize *= 2;

  // H_NE_internal_storage_ = std::make_unique<std::vector<double>>(qSize);
  H_NE_internal_storage_.resize(1, qSize);
  H_NE_view_ =
      std::make_shared<siconos::algebra::MapType>(H_NE_internal_storage_.data(), 1, qSize);
  H_NE_view_->setZero();
  NewtonEulerR::initialize(inter);
  //  _isContact=1;
}

void siconos::modeling::NewtonEuler1DR::computeH_NE_(double time,
                                                     siconos::modeling::Interaction& inter,
                                                     const siconos::algebra::BlockVector& q0) {
  DEBUG_BEGIN(
      "siconos::modeling::NewtonEuler1DR::computeJacobianhOver_q(double time, Interaction& "
      "inter, "
      "std::shared_ptr<siconos::algebra::BlockVector> q0 ) \n");
  DEBUG_PRINTF("with time =  %f\n", time);
  DEBUG_PRINTF("with inter =  %p\n", &inter);

  H_NE_view_->row(0).leftCols(3) = nc_;
  if (inter.has2Bodies()) {
    H_NE_view_->row(0).segment(7, 3) = -nc_;
  }

  for (unsigned int iDS = 0; iDS < q0.numberOfBlocks(); iDS++) {
    auto q = q0.vector(iDS);
    double sign = 1.0;
    DEBUG_PRINTF("siconos::modeling::NewtonEuler1DR::computeJachq : ds%d->q :", iDS);
    DEBUG_EXPR_WE(siconos::algebra::print(*q));

    boost::math::quaternion<double> quatGP;
    if (iDS == 0) {
      boost::math::quaternion<double> quatAux(0, contactPoint1_(0) - (*q)(0),
                                              contactPoint1_(1) - (*q)(1),
                                              contactPoint1_(2) - (*q)(2));
      quatGP = quatAux;
    } else {
      sign = -1.0;
      // cout<<"siconos::modeling::NewtonEuler1DR::computeJachq sign is -1 \n";
      boost::math::quaternion<double> quatAux(0, contactPoint2_(0) - (*q)(0),
                                              contactPoint2_(1) - (*q)(1),
                                              contactPoint2_(2) - (*q)(2));
      quatGP = quatAux;
    }
    DEBUG_PRINTF("siconos::modeling::NewtonEuler1DR::computeJachq :GP :%lf, %lf, %lf\n",
                 quatGP.R_component_2(), quatGP.R_component_3(), quatGP.R_component_4());
    DEBUG_PRINTF("siconos::modeling::NewtonEuler1DR::computeJachq :Q :%e,%e, %e, %e\n",
                 (*q)(3), (*q)(4), (*q)(5), (*q)(6));
    boost::math::quaternion<double> quatQ((*q)(3), (*q)(4), (*q)(5), (*q)(6));
    boost::math::quaternion<double> quatcQ((*q)(3), -(*q)(4), -(*q)(5), -(*q)(6));
    boost::math::quaternion<double> quatBuff;
    boost::math::quaternion<double> _2qiquatGP;
    _2qiquatGP = quatGP;
    _2qiquatGP *= 2 * ((*q)(3));
    quatBuff = (quatGP * quatQ) + (quatcQ * quatGP) - _2qiquatGP;

    DEBUG_PRINTF("siconos::modeling::NewtonEuler1DR::computeJachq :quattBuuf : %e,%e,%e \n",
                 quatBuff.R_component_2(), quatBuff.R_component_3(), quatBuff.R_component_4());

    H_NE_view_->setValue(
        0, 7 * iDS + 3,
        sign * (quatBuff.R_component_2() * nc_(0) + quatBuff.R_component_3() * nc_(1) +
                quatBuff.R_component_4() * nc_(2)));
    // cout<<"WARNING NewtonEuler1DR set jachq \n";
    // jacobianhOver_q_->setValue(0,7*iDS+3,0);
    for (unsigned int i = 1; i < 4; i++) {
      boost::math::quaternion<double> quatei(0, (i == 1) ? 1 : 0, (i == 2) ? 1 : 0,
                                             (i == 3) ? 1 : 0);
      _2qiquatGP = quatGP;
      _2qiquatGP *= 2 * ((*q)(3 + i));
      quatBuff = quatei * quatcQ * quatGP - quatGP * quatQ * quatei - _2qiquatGP;
      H_NE_view_->setValue(
          0, 7 * iDS + 3 + i,
          sign * (quatBuff.R_component_2() * nc_(0) + quatBuff.R_component_3() * nc_(1) +
                  quatBuff.R_component_4() * nc_(2)));
    }
  }

  DEBUG_EXPR(siconos::algebra::print(*jacobianhOver_q_););
  DEBUG_END(
      "siconos::modeling::NewtonEuler1DR::computeJacobianhOver_q(double time, Interaction& "
      "inter, "
      "std::shared_ptr<siconos::algebra::BlockVector> q0 \n");
}

void siconos::modeling::NewtonEuler1DR::computeH_NE_prod_T(
    const Interaction& inter, const siconos::algebra::BlockVector& q0) {
  DEBUG_BEGIN(
      "siconos::modeling::NewtonEuler1DR::computeH_NE_prod_T(Interaction& inter, "
      "std::shared_ptr<siconos::algebra::BlockVector> q0 \n")

  if (q0.numberOfBlocks() > 1) {
    NIcomputeJachqTFromContacts(*q0.vector(0), *q0.vector(1));
  } else {
    NIcomputeJachqTFromContacts(*q0.vector(0));
  }

  DEBUG_END(
      "siconos::modeling::NewtonEuler1DR::computeH_NE_prod_T(Interaction& inter, "
      "std::shared_ptr<siconos::algebra::BlockVector> q0) \n");
}

double siconos::modeling::NewtonEuler1DR::distance() const {
  siconos::algebra::SiconosVector3 dpc = contactPoint2_ - contactPoint1_;
  return dpc.norm() * (nc_.dot(dpc) >= 0 ? -1 : 1);
}

void siconos::modeling::NewtonEuler1DR::computehFromRelativeContactPoints(
    double time, const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>& q2,
    Eigen::Ref<siconos::algebra::SiconosVector> y) {
  // Contact points and normal are stored as relative to q1 and q2, if
  // no q2 then pc2 and normal are absolute.

  // Update pc1 based on q0 and relPc1
  boost::math::quaternion<double> qq1(q1(3), q1(4), q1(5), q1(6));
  boost::math::quaternion<double> qpc1(0, relPc1_(0), relPc1_(1), relPc1_(2));

  // apply q1 rotation and add
  qpc1 = qq1 * qpc1 / qq1;
  contactPoint1_(0) = qpc1.R_component_2() + q1(0);
  contactPoint1_(1) = qpc1.R_component_3() + q1(1);
  contactPoint1_(2) = qpc1.R_component_4() + q1(2);

  if (q2) {
    // Update pc2 based on q0 and relPc2
    boost::math::quaternion<double> qq2((*q2)(3), (*q2)(4), (*q2)(5), (*q2)(6));
    boost::math::quaternion<double> qpc2(0, relPc2_(0), relPc2_(1), relPc2_(2));

    // apply q2 rotation and add
    qpc2 = qq2 * qpc2 / qq2;
    contactPoint2_(0) = qpc2.R_component_2() + (*q2)(0);
    contactPoint2_(1) = qpc2.R_component_3() + (*q2)(1);
    contactPoint2_(2) = qpc2.R_component_4() + (*q2)(2);

    // same for normal
    boost::math::quaternion<double> qnc(0, relNc_(0), relNc_(1), relNc_(2));
    qnc = qq2 * qnc / qq2;
    nc_ << qnc.R_component_2(), qnc.R_component_3(), qnc.R_component_4();
    NewtonEulerR::computeh(q1, q2, y);
  } else {
    contactPoint2_ = relPc2_;
    nc_ = relNc_;
    NewtonEulerR::computeh(q1, std::nullopt, y);
  }
}
