/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include "NewtonEulerNormalVelocityR.hpp"

#include <boost/math/quaternion.hpp>

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "RotationQuaternion.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

//  #define DEBUG_NOCOLOR
//  #define DEBUG_STDOUT
//  #define DEBUG_MESSAGES
#include "siconos_debug.h"

namespace {
/** \return the skew-symmetric matrix \f$ [a]_\times \f$ such that
 *  \f$ [a]_\times b = a \times b \f$ */
siconos::algebra::SiconosMatrix33 skew(const Eigen::Ref<const Eigen::Vector3d>& a) {
  siconos::algebra::SiconosMatrix33 m;
  m << 0.0, -a.z(), a.y(), a.z(), 0.0, -a.x(), -a.y(), a.x(), 0.0;
  return m;
}
}  // namespace

/* See the class documentation and DevNotes, chapter 'NewtonEulerR: computation of
 * \nabla_q H', subsection 'Case FC3D: using the local frame local velocities'. */
void siconos::modeling::NewtonEulerNormalVelocityR::computeContactFrame_(bool withTangents) {
  if (!withTangents) {
    contactFrame_.row(0) = nc_;
    return;
  }
  assert(std::abs(nc_.norm() - 1.0) < 1e-6 &&
         "NewtonEulerNormalVelocityR::computeContactFrame_: the normal vector is not a unit "
         "vector");
  siconos::algebra::SiconosVector3 n = nc_, t1, t2;
  if (!siconos::geometry::orthoBaseFromVector(n, t1, t2))
    THROW_EXCEPTION(
        "NewtonEulerNormalVelocityR::computeContactFrame_: cannot build the contact frame "
        "from a zero "
        "normal vector");
  contactFrame_.row(0) = n;
  contactFrame_.row(1) = t1;
  contactFrame_.row(2) = t2;
}

void siconos::modeling::NewtonEulerNormalVelocityR::computeContactVelocityJacobian_(
    const siconos::algebra::BlockVector& q0, Eigen::Index nRows) {
  const auto Rc = contactFrame_.topRows(nRows);
  siconos::algebra::SiconosMatrix33 R;
  for (Eigen::Index body = 0; body < static_cast<Eigen::Index>(q0.numberOfBlocks()); ++body) {
    const auto& q = *q0.vector(body);
    // body 1 contributes +R_c [I, [G1 - P]x R1], body 2 contributes -R_c [I, [G2 - P]x R2]
    const double sign = (body == 0) ? 1.0 : -1.0;
    siconos::geometry::computeRotationMatrix(q, R);
    contactMatrix_->block(0, 6 * body, nRows, 3) = sign * Rc;
    contactMatrix_->block(0, 6 * body + 3, nRows, 3) =
        sign * Rc * skew(q.head<3>() - contactPoint1_) * R;
  }
  DEBUG_EXPR(siconos::algebra::print(*contactMatrix_););
}

void siconos::modeling::NewtonEulerNormalVelocityR::initialize(Interaction& inter) {
  auto qSize = 7;
  if (inter.has2Bodies()) qSize *= 2;

  constraintJacobian_storage_.resize(1, qSize);
  constraintJacobian_view_ = std::make_shared<siconos::algebra::MapType>(
      constraintJacobian_storage_.data(), 1, qSize);
  constraintJacobian_view_->setZero();
  NewtonEulerR::initialize(inter);
  //  _isContact=1;
}

void siconos::modeling::NewtonEulerNormalVelocityR::computeConstraintJacobian(
    double time, siconos::modeling::Interaction& inter,
    const siconos::algebra::BlockVector& q0) {
  DEBUG_BEGIN(
      "siconos::modeling::NewtonEulerNormalVelocityR::computeContactMatrix(double time, "
      "Interaction& "
      "inter, "
      "std::shared_ptr<siconos::algebra::BlockVector> q0 ) \n");
  DEBUG_PRINTF("with time =  %f\n", time);
  DEBUG_PRINTF("with inter =  %p\n", &inter);

  constraintJacobian_view_->row(0).leftCols(3) = nc_;
  if (inter.has2Bodies()) {
    constraintJacobian_view_->row(0).segment(7, 3) = -nc_;
  }

  for (unsigned int iDS = 0; iDS < q0.numberOfBlocks(); iDS++) {
    auto q = q0.vector(iDS);
    double sign = 1.0;
    DEBUG_PRINTF("siconos::modeling::NewtonEulerNormalVelocityR::computeJachq : ds%d->q :",
                 iDS);
    DEBUG_EXPR_WE(siconos::algebra::print(*q));

    boost::math::quaternion<double> quatGP;
    if (iDS == 0) {
      boost::math::quaternion<double> quatAux(0, contactPoint1_(0) - (*q)(0),
                                              contactPoint1_(1) - (*q)(1),
                                              contactPoint1_(2) - (*q)(2));
      quatGP = quatAux;
    } else {
      sign = -1.0;
      // cout<<"siconos::modeling::NewtonEulerNormalVelocityR::computeJachq sign is -1 \n";
      boost::math::quaternion<double> quatAux(0, contactPoint2_(0) - (*q)(0),
                                              contactPoint2_(1) - (*q)(1),
                                              contactPoint2_(2) - (*q)(2));
      quatGP = quatAux;
    }
    DEBUG_PRINTF(
        "siconos::modeling::NewtonEulerNormalVelocityR::computeJachq :GP :%lf, %lf, %lf\n",
        quatGP.R_component_2(), quatGP.R_component_3(), quatGP.R_component_4());
    DEBUG_PRINTF(
        "siconos::modeling::NewtonEulerNormalVelocityR::computeJachq :Q :%e,%e, %e, %e\n",
        (*q)(3), (*q)(4), (*q)(5), (*q)(6));
    boost::math::quaternion<double> quatQ((*q)(3), (*q)(4), (*q)(5), (*q)(6));
    boost::math::quaternion<double> quatcQ((*q)(3), -(*q)(4), -(*q)(5), -(*q)(6));
    boost::math::quaternion<double> quatBuff;
    boost::math::quaternion<double> _2qiquatGP;
    _2qiquatGP = quatGP;
    _2qiquatGP *= 2 * ((*q)(3));
    quatBuff = (quatGP * quatQ) + (quatcQ * quatGP) - _2qiquatGP;

    DEBUG_PRINTF(
        "siconos::modeling::NewtonEulerNormalVelocityR::computeJachq :quattBuuf : %e,%e,%e \n",
        quatBuff.R_component_2(), quatBuff.R_component_3(), quatBuff.R_component_4());

    constraintJacobian_view_->setValue(
        0, 7 * iDS + 3,
        sign * (quatBuff.R_component_2() * nc_(0) + quatBuff.R_component_3() * nc_(1) +
                quatBuff.R_component_4() * nc_(2)));
    // cout<<"WARNING NewtonEulerNormalVelocityR set G \n";
    // jacobianhOver_q_->setValue(0,7*iDS+3,0);
    for (unsigned int i = 1; i < 4; i++) {
      boost::math::quaternion<double> quatei(0, (i == 1) ? 1 : 0, (i == 2) ? 1 : 0,
                                             (i == 3) ? 1 : 0);
      _2qiquatGP = quatGP;
      _2qiquatGP *= 2 * ((*q)(3 + i));
      quatBuff = quatei * quatcQ * quatGP - quatGP * quatQ * quatei - _2qiquatGP;
      constraintJacobian_view_->setValue(
          0, 7 * iDS + 3 + i,
          sign * (quatBuff.R_component_2() * nc_(0) + quatBuff.R_component_3() * nc_(1) +
                  quatBuff.R_component_4() * nc_(2)));
    }
  }

  DEBUG_EXPR(siconos::algebra::print(*constraintJacobian_view_););
  DEBUG_END(
      "siconos::modeling::NewtonEulerNormalVelocityR::computeContactMatrix(double time, "
      "Interaction& "
      "inter, "
      "std::shared_ptr<siconos::algebra::BlockVector> q0 \n");
}

void siconos::modeling::NewtonEulerNormalVelocityR::computeContactMatrix(
    const Interaction& inter, const siconos::algebra::BlockVector& q0) {
  DEBUG_BEGIN(
      "siconos::modeling::NewtonEulerNormalVelocityR::computeContactMatrix(Interaction& "
      "inter, "
      "std::shared_ptr<siconos::algebra::BlockVector> q0 \n")

  computeContactFrame_(false);
  computeContactVelocityJacobian_(q0, 1);

  DEBUG_END(
      "siconos::modeling::NewtonEulerNormalVelocityR::computeContactMatrix(Interaction& "
      "inter, "
      "std::shared_ptr<siconos::algebra::BlockVector> q0) \n");
}

double siconos::modeling::NewtonEulerNormalVelocityR::distance() const {
  siconos::algebra::SiconosVector3 dpc = contactPoint2_ - contactPoint1_;
  return dpc.norm() * (nc_.dot(dpc) >= 0 ? -1 : 1);
}

void siconos::modeling::NewtonEulerNormalVelocityR::computehFromRelativeContactPoints(
    double time, const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2,
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

void siconos::modeling::NewtonEulerNormalVelocityR::
    computeContactPointsFromRelativeContactPoints(
        const siconos::algebra::BlockVector& q0, const siconos::algebra::SiconosVector3& r_pc1,
        const siconos::algebra::SiconosVector3& r_pc2,
        const siconos::algebra::SiconosVector3& r_nc,
        const siconos::algebra::SiconosVector3& r_t1,
        const siconos::algebra::SiconosVector3& r_t2, siconos::algebra::SiconosVector3& pc1,
        siconos::algebra::SiconosVector3& pc2, siconos::algebra::SiconosVector3& nc,
        siconos::algebra::SiconosVector3& t1, siconos::algebra::SiconosVector3& t2) {
  // Contact points and normal are stored as relative to q1 and q2, if
  // no q2 then pc2 and normal are absolute.

  // Update pc1 based on q0 and r_pc1
  const auto& q1 = *q0.vector(0);
  boost::math::quaternion<double> qq1(q1(3), q1(4), q1(5), q1(6));
  boost::math::quaternion<double> qpc1(0, r_pc1(0), r_pc1(1), r_pc1(2));

  // apply q1 rotation and add
  qpc1 = qq1 * qpc1 / qq1;
  pc1(0) = qpc1.R_component_2() + q1(0);
  pc1(1) = qpc1.R_component_3() + q1(1);
  pc1(2) = qpc1.R_component_4() + q1(2);

  if (q0.numberOfBlocks() > 1) {
    // Update pc2 based on q0 and r_pc2
    const auto& q2 = *q0.vector(1);
    boost::math::quaternion<double> qq2(q2(3), q2(4), q2(5), q2(6));
    boost::math::quaternion<double> qpc2(0, r_pc2(0), r_pc2(1), r_pc2(2));

    // apply q2 rotation and add
    qpc2 = qq2 * qpc2 / qq2;
    pc2(0) = qpc2.R_component_2() + q2(0);
    pc2(1) = qpc2.R_component_3() + q2(1);
    pc2(2) = qpc2.R_component_4() + q2(2);

    // same for normal
    boost::math::quaternion<double> qnc(0, r_nc(0), r_nc(1), r_nc(2));
    qnc = qq2 * qnc / qq2;
    nc(0) = qnc.R_component_2();
    nc(1) = qnc.R_component_3();
    nc(2) = qnc.R_component_4();

    // tangent 1
    boost::math::quaternion<double> qt1(0, r_t1(0), r_t1(1), r_t1(2));
    qt1 = qq2 * qt1 / qq2;
    t1(0) = qt1.R_component_2();
    t1(1) = qt1.R_component_3();
    t1(2) = qt1.R_component_4();

    // tangent 2
    boost::math::quaternion<double> qt2(0, r_t2(0), r_t2(1), r_t2(2));
    qt2 = qq2 * qt2 / qq2;
    t2(0) = qt2.R_component_2();
    t2(1) = qt2.R_component_3();
    t2(2) = qt2.R_component_4();
  } else {
    // Single body case - use relative coordinates as absolute
    pc2 = r_pc2;
    nc = r_nc;
    t1 = r_t1;
    t2 = r_t2;
  }
}
