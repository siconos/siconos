/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#include "NewtonEulerRelationsTest.hpp"

#include <memory>
#include <random>
#include <sstream>

#include "RotationQuaternion.hpp"
#include "SiconosKernel.hpp"

using siconos::algebra::BlockVector;
using siconos::algebra::SiconosMatrix;
using siconos::algebra::SiconosMatrix33;
using siconos::algebra::SiconosVector;
using siconos::algebra::SiconosVector3;
using siconos::algebra::SiconosVector6;
using siconos::algebra::SiconosVector7;
using namespace siconos::modeling;

CPPUNIT_TEST_SUITE_REGISTRATION(NewtonEulerRelationsTest);

void NewtonEulerRelationsTest::setUp() {}

void NewtonEulerRelationsTest::tearDown() {}

namespace {

constexpr int nTrials = 20;
constexpr double tol = 1e-12;

template <class Vec>
Vec randomVector(std::mt19937& gen) {
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  Vec v;
  for (Eigen::Index i = 0; i < v.size(); ++i) v(i) = dist(gen);
  return v;
}

SiconosVector7 randomPose(std::mt19937& gen) {
  SiconosVector7 q = randomVector<SiconosVector7>(gen);
  q.tail<4>().normalize();
  return q;
}

/** Builds the relation and the interaction for one or two random bodies, computes
 *  H and returns || H * twists - y_ref ||.
 */
template <class Rel>
double jacobianError(std::shared_ptr<NonSmoothLaw> nslaw, bool twoBodies, std::mt19937& gen) {
  const auto ySize = nslaw->size();
  auto rel = std::make_shared<Rel>();
  Interaction inter(nslaw, rel);
  inter.setHas2Bodies(twoBodies);
  inter.setDSSizes(twoBodies ? 12 : 6);

  const SiconosVector7 q1 = randomPose(gen);
  const SiconosVector7 q2 = randomPose(gen);
  const SiconosMatrix33 inertia = SiconosMatrix33::Identity();
  NewtonEulerDS ds1(q1, SiconosVector6::Zero(), 1.0, inertia, siconos::algebra::copy_t);
  NewtonEulerDS ds2(q2, SiconosVector6::Zero(), 1.0, inertia, siconos::algebra::copy_t);
  inter.initialize_read_access_to_ds_variables(ds1, twoBodies ? ds2 : ds1);

  // Contact points (in the body frames) and normal (in the frame of body 2, absolute if
  // there is only one body)
  rel->setRelPc1(randomVector<SiconosVector3>(gen));
  rel->setRelPc2(randomVector<SiconosVector3>(gen));
  rel->setRelNc(randomVector<SiconosVector3>(gen).normalized());
  SiconosVector y(ySize);
  if (twoBodies) {
    rel->computehFromRelativeContactPoints(0., q1, q2, y);
  } else {
    rel->computehFromRelativeContactPoints(0., q1, std::nullopt, y);
  }

  BlockVector q0;
  q0.insertPtr(std::make_shared<SiconosVector>(q1));
  if (twoBodies) q0.insertPtr(std::make_shared<SiconosVector>(q2));
  rel->computeContactMatrix(inter, q0);
  const SiconosMatrix H = rel->contactMatrix();

  // Twists: translational velocity (absolute frame), angular velocity (body frame)
  const SiconosVector6 twist1 = randomVector<SiconosVector6>(gen);
  const SiconosVector6 twist2 = randomVector<SiconosVector6>(gen);
  SiconosVector twists(twoBodies ? 12 : 6);
  twists.head<6>() = twist1;
  if (twoBodies) twists.tail<6>() = twist2;
  const SiconosVector yH = H * twists;

  // Reference: relative velocity at pc1 in the contact frame
  SiconosMatrix33 R1, R2;
  siconos::geometry::computeRotationMatrix(q1, R1);
  siconos::geometry::computeRotationMatrix(q2, R2);
  const SiconosVector3 pc = rel->pc1();
  SiconosVector3 n = rel->nc(), t1, t2;
  CPPUNIT_ASSERT(siconos::geometry::orthoBaseFromVector(n, t1, t2));

  const SiconosVector3 omega1 = R1 * twist1.tail<3>();
  const SiconosVector3 omega2 = R2 * twist2.tail<3>();
  SiconosVector3 u = twist1.head<3>() + omega1.cross(pc - q1.head<3>());
  SiconosVector3 omegaRel = omega1;
  if (twoBodies) {
    u -= twist2.head<3>() + omega2.cross(pc - q2.head<3>());
    omegaRel -= omega2;
  }

  SiconosVector yRef(ySize);
  yRef(0) = n.dot(u);
  if (ySize >= 3) {
    yRef(1) = t1.dot(u);
    yRef(2) = t2.dot(u);
  }
  if (ySize == 5) {
    yRef(3) = t1.dot(omegaRel);
    yRef(4) = t2.dot(omegaRel);
  }
  return (yH - yRef).norm();
}

/** T(q) for one body */
siconos::algebra::SiconosMatrix76 computeT(const SiconosVector7& q) {
  siconos::algebra::SiconosMatrix76 T = siconos::algebra::SiconosMatrix76::Zero();
  T.topLeftCorner<3, 3>().setIdentity();
  siconos::modeling::newton_euler::computeT(q, T);
  return T;
}

/** Builds a contact with pc1 - pc2 = g0 n, checks that the first row of J_h (G) is the
 *  gradient of the gap (central finite differences) and that J_h T v is equal to the first
 *  row of H v. Returns the largest of the two relative errors.
 */
template <class Rel>
double gapGradientError(std::shared_ptr<NonSmoothLaw> nslaw, bool twoBodies,
                        std::mt19937& gen) {
  const auto ySize = nslaw->size();
  auto rel = std::make_shared<Rel>();
  Interaction inter(nslaw, rel);
  inter.setHas2Bodies(twoBodies);
  inter.setDSSizes(twoBodies ? 12 : 6);

  SiconosVector7 q1 = randomPose(gen);
  SiconosVector7 q2 = randomPose(gen);
  const SiconosMatrix33 inertia = SiconosMatrix33::Identity();
  NewtonEulerDS ds1(q1, SiconosVector6::Zero(), 1.0, inertia, siconos::algebra::copy_t);
  NewtonEulerDS ds2(q2, SiconosVector6::Zero(), 1.0, inertia, siconos::algebra::copy_t);
  inter.initialize_read_access_to_ds_variables(ds1, twoBodies ? ds2 : ds1);

  // Geometry in the inertial frame: P1 = P2 + g0 n
  SiconosMatrix33 R1, R2;
  siconos::geometry::computeRotationMatrix(q1, R1);
  siconos::geometry::computeRotationMatrix(q2, R2);
  const SiconosVector3 n = randomVector<SiconosVector3>(gen).normalized();
  const SiconosVector3 P2 = randomVector<SiconosVector3>(gen);
  const SiconosVector3 P1 = P2 + 0.3 * n;
  rel->setRelPc1(R1.transpose() * (P1 - q1.head<3>()));
  if (twoBodies) {
    rel->setRelPc2(R2.transpose() * (P2 - q2.head<3>()));
    rel->setRelNc(R2.transpose() * n);
  } else {
    rel->setRelPc2(P2);
    rel->setRelNc(n);
  }

  SiconosVector y(ySize);
  auto gap = [&](const SiconosVector7& a, const SiconosVector7& b) {
    if (twoBodies)
      rel->computehFromRelativeContactPoints(0., a, b, y);
    else
      rel->computehFromRelativeContactPoints(0., a, std::nullopt, y);
    return rel->distance();
  };

  // Gradient of the gap by central finite differences
  const int nq = twoBodies ? 14 : 7;
  SiconosVector gradFD(nq);
  const double eps = 1e-6;
  for (int k = 0; k < nq; ++k) {
    SiconosVector7 a = q1, b = q2;
    SiconosVector7& qk = (k < 7) ? a : b;
    qk(k % 7) += eps;
    const double gp = gap(a, b);
    qk(k % 7) -= 2 * eps;
    const double gm = gap(a, b);
    gradFD(k) = (gp - gm) / (2 * eps);
  }

  // J_h and H at q
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.3, gap(q1, q2), 1e-12);
  rel->computeJach(0., inter);
  const SiconosVector Jh = rel->constraintJacobian().row(0).transpose();
  const SiconosMatrix H = rel->contactMatrix();
  const double errGrad = (Jh - gradFD).norm() / gradFD.norm();

  // J_h T v == (H v)_0
  const SiconosVector6 v1 = randomVector<SiconosVector6>(gen);
  const SiconosVector6 v2 = randomVector<SiconosVector6>(gen);
  SiconosVector v(twoBodies ? 12 : 6), qdot(nq);
  v.head<6>() = v1;
  qdot.head<7>() = computeT(q1) * v1;
  if (twoBodies) {
    v.tail<6>() = v2;
    qdot.tail<7>() = computeT(q2) * v2;
  }
  const double Hv = H.row(0).dot(v);
  const double errVel = std::abs(Jh.dot(qdot) - Hv) / std::max(1.0, std::abs(Hv));
  return std::max(errGrad, errVel);
}

template <class Rel>
void checkGapGradient(const std::string& name, std::shared_ptr<NonSmoothLaw> nslaw) {
  std::mt19937 gen(7);
  for (bool twoBodies : {false, true}) {
    for (int trial = 0; trial < nTrials; ++trial) {
      const double err = gapGradientError<Rel>(nslaw, twoBodies, gen);
      std::ostringstream msg;
      msg << name << " (" << (twoBodies ? "two bodies" : "one body") << ", trial " << trial
          << "): J_h vs gradient of the gap, relative error = " << err;
      CPPUNIT_ASSERT_MESSAGE(msg.str(), err < 1e-7);
    }
  }
}

template <class Rel>
void checkJacobian(const std::string& name, std::shared_ptr<NonSmoothLaw> nslaw) {
  std::mt19937 gen(42);
  for (bool twoBodies : {false, true}) {
    for (int trial = 0; trial < nTrials; ++trial) {
      const double err = jacobianError<Rel>(nslaw, twoBodies, gen);
      std::ostringstream msg;
      msg << name << " (" << (twoBodies ? "two bodies" : "one body") << ", trial " << trial
          << "): || H * twist - u_ref || = " << err;
      CPPUNIT_ASSERT_MESSAGE(msg.str(), err < tol);
    }
  }
}

}  // namespace

void NewtonEulerRelationsTest::testNewtonEuler1DRJacobian() {
  checkJacobian<NewtonEulerNormalVelocityR>("NewtonEulerNormalVelocityR",
                                            std::make_shared<NewtonImpactNSL>(0.));
}

void NewtonEulerRelationsTest::testNewtonEuler3DRJacobian() {
  checkJacobian<NewtonEulerVelocityR>(
      "NewtonEulerVelocityR", std::make_shared<NewtonImpactFrictionNSL>(0., 0., 0.5, 3));
}

void NewtonEulerRelationsTest::testNewtonEulerVelocityAngularVelocityRJacobian() {
  checkJacobian<NewtonEulerVelocityAngularVelocityR>(
      "NewtonEulerVelocityAngularVelocityR",
      std::make_shared<NewtonImpactRollingFrictionNSL>(0., 0., 0.5, 0.1, 5));
}

void NewtonEulerRelationsTest::testNewtonEuler1DRGapGradient() {
  checkGapGradient<NewtonEulerNormalVelocityR>("NewtonEulerNormalVelocityR",
                                               std::make_shared<NewtonImpactNSL>(0.));
}

void NewtonEulerRelationsTest::testNewtonEuler3DRGapGradient() {
  checkGapGradient<NewtonEulerVelocityR>(
      "NewtonEulerVelocityR", std::make_shared<NewtonImpactFrictionNSL>(0., 0., 0.5, 3));
}

void NewtonEulerRelationsTest::testNewtonEulerVelocityAngularVelocityRGapGradient() {
  checkGapGradient<NewtonEulerVelocityAngularVelocityR>(
      "NewtonEulerVelocityAngularVelocityR",
      std::make_shared<NewtonImpactRollingFrictionNSL>(0., 0., 0.5, 0.1, 5));
}

void NewtonEulerRelationsTest::testNewtonEulerRSecondOrderTerms() {
  std::mt19937 gen(3);
  for (bool twoBodies : {false, true}) {
    const int m = 2;
    const int nq = twoBodies ? 14 : 7;
    auto rel = std::make_shared<NewtonEulerR>();
    Interaction inter(std::make_shared<EqualityConditionNSL>(m), rel);
    inter.setHas2Bodies(twoBodies);
    inter.setDSSizes(twoBodies ? 12 : 6);

    const SiconosVector7 q1 = randomPose(gen), q2 = randomPose(gen);
    const SiconosVector6 v1 = randomVector<SiconosVector6>(gen);
    const SiconosVector6 v2 = randomVector<SiconosVector6>(gen);
    const SiconosMatrix33 inertia = SiconosMatrix33::Identity();
    auto ds1 = std::make_shared<NewtonEulerDS>(q1, v1, 1.0, inertia, siconos::algebra::copy_t);
    auto ds2 = std::make_shared<NewtonEulerDS>(q2, v2, 1.0, inertia, siconos::algebra::copy_t);
    *ds1->dotq() = computeT(q1) * v1;
    *ds2->dotq() = computeT(q2) * v2;

    SiconosMatrix J(m, nq), Jdot(m, nq);
    for (int i = 0; i < m; ++i)
      for (int j = 0; j < nq; ++j) {
        J(i, j) = std::uniform_real_distribution<double>(-1., 1.)(gen);
        Jdot(i, j) = std::uniform_real_distribution<double>(-1., 1.)(gen);
      }
    rel->setConstantConstraintJacobian(J);
    inter.initialize_read_access_to_ds_variables(*ds1, twoBodies ? *ds2 : *ds1);

    // Set after initialize: H_NE_dot must be allocated at first use.
    rel->setComputeConstraintJacobianDotFunction(
        [&Jdot](const BlockVector&, const BlockVector&,
                Eigen::Ref<siconos::algebra::MapType> r) { r = Jdot; });
    rel->computeContactMatrixDotV(0., inter, ds1, twoBodies ? ds2 : ds1);

    // Reference: Jdot qdot + d/dt(J T(q)) v, T being linear in q (central differences are
    // exact up to round-off)
    SiconosVector qdot(nq), v(twoBodies ? 12 : 6);
    qdot.head<7>() = *ds1->dotq();
    v.head<6>() = v1;
    if (twoBodies) {
      qdot.tail<7>() = *ds2->dotq();
      v.tail<6>() = v2;
    }
    const double eps = 1e-3;
    SiconosVector ref = Jdot * qdot;
    for (int body = 0; body < (twoBodies ? 2 : 1); ++body) {
      const SiconosVector7 q = (body == 0) ? q1 : q2;
      const SiconosVector7 dq = qdot.segment<7>(7 * body);
      const siconos::algebra::SiconosMatrix76 Tdot =
          (computeT(q + eps * dq) - computeT(q - eps * dq)) / (2 * eps);
      ref += J.middleCols(7 * body, 7) * (Tdot * v.segment<6>(6 * body));
    }
    const double err = (SiconosVector(rel->contactMatrixDotV()) - ref).norm();
    std::ostringstream msg;
    msg << "NewtonEulerR second-order terms (" << (twoBodies ? "two bodies" : "one body")
        << "): error = " << err;
    CPPUNIT_ASSERT_MESSAGE(msg.str(), err < 1e-12);
  }
}
