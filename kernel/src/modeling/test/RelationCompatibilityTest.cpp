/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
// This test uses the deprecated names on purpose.
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

#include "RelationCompatibilityTest.hpp"

#include <memory>
#include <type_traits>

#include "SiconosKernel.hpp"

using siconos::algebra::BlockVector;
using siconos::algebra::SiconosMatrix;
using siconos::algebra::SiconosVector;
using namespace siconos::modeling;

CPPUNIT_TEST_SUITE_REGISTRATION(RelationCompatibilityTest);

namespace {

/** q = (1, 2, 3) as a block vector (BlockVector is neither copyable nor movable) */
void fillQ(BlockVector& q) {
  auto v = std::make_shared<SiconosVector>(3);
  *v << 1., 2., 3.;
  q.insertPtr(v);
}

/** allocate a 1 x 3 H in a Lagrangian relation (normally done by initialize) */
template <class R>
struct WithH : R {
  WithH() {
    this->contactMatrix_storage_ = std::make_unique<std::vector<double>>(3, 0.);
    this->contactMatrix_view_ = std::make_shared<siconos::algebra::MapType>(
        this->contactMatrix_storage_->data(), 1, 3);
  }
};

// --- Scleronomic ---
/** a user class written with the former API: overrides computeJacobianhOver_q (without
 *  `override`). (The former member names, e.g. jacobianhOver_q_view_, have been removed.) */
struct OldScleronomousR : WithH<LagrangianScleronomousR> {
  int calls = 0;
  void computeJacobianhOver_q(const BlockVector& q) {
    ++calls;
    (*contactMatrix_view_)(0, 0) = q(0);
  }
};
/** a class written with the new API */
struct NewScleronomousR : WithH<LagrangianScleronomousR> {
  int calls = 0;
  void computeContactMatrix(const BlockVector& q) override {
    ++calls;
    (*contactMatrix_view_)(0, 1) = q(1);
  }
};

// --- Rheonomic ---
struct OldRheonomousR : WithH<LagrangianRheonomousR> {
  int callsH = 0, callsDt = 0;
  OldRheonomousR() { dhdt_ = std::make_shared<SiconosVector>(1); }
  void computeJacobianhOver_q(const BlockVector& q, double time) {
    ++callsH;
    (*contactMatrix_view_)(0, 2) = time;
  }
  void computehdot(const BlockVector& q, double time) {
    ++callsDt;
    (*dhdt_)(0) = 2 * time;
  }
};

// --- Compliant ---
struct OldCompliantR : WithH<LagrangianCompliantR> {
  int callsH = 0, callsD = 0;
  OldCompliantR() { D_ = std::make_shared<SiconosMatrix>(1, 1); }
  void computeJacobianhOver_q(const BlockVector& q,
                              const Eigen::Ref<const SiconosVector>& lambda) {
    ++callsH;
  }
  void computeJacobianhOver_lambda(const BlockVector& q,
                                   const Eigen::Ref<const SiconosVector>& lambda) {
    ++callsD;
    (*D_)(0, 0) = lambda(0);
  }
};

// --- Newton-Euler ---
/** allocate G (1 x 7) and H (1 x 6), normally done by initialize */
struct WithGH : NewtonEulerR {
  WithGH() {
    constraintJacobian_storage_.resize(1, 7);
    constraintJacobian_storage_.setZero();
    constraintJacobian_view_ =
        std::make_shared<siconos::algebra::MapType>(constraintJacobian_storage_.data(), 1, 7);
    contactMatrix_ = std::make_shared<SiconosMatrix>(1, 6);
    contactMatrix_->setZero();
    T_buffer_.setZero();  // constant part of T, as in NewtonEulerR::initialize
    T_buffer_.topLeftCorner<3, 3>().setIdentity();
  }
  /** computeConstraintJacobian is protected: expose it */
  void callComputeG(Interaction& inter, const BlockVector& q) {
    computeConstraintJacobian(0., inter, q);
  }
};
struct OldNewtonEulerR : WithGH {
  int callsG = 0, callsH = 0;
  void computeH_NE_(double, Interaction&, const BlockVector&) { ++callsG; }
  void computeH_NE_prod_T(const Interaction&, const BlockVector&) {
    ++callsH;
    (*contactMatrix_)(0, 0) = 5.;
  }
};
struct NewNewtonEulerR : WithGH {
  int callsH = 0;
  void computeContactMatrix(const Interaction&, const BlockVector&) override { ++callsH; }
};

}  // namespace

void RelationCompatibilityTest::testLagrangianScleronomousR() {
  BlockVector q;
  fillQ(q);

  // former override, called through the new name (the one used by the kernel)
  OldScleronomousR oldR;
  oldR.computeContactMatrix(q);
  CPPUNIT_ASSERT_EQUAL(1, oldR.calls);
  CPPUNIT_ASSERT_EQUAL(1., oldR.contactMatrix()(0, 0));
  oldR.computeJacobianhOver_q(q);  // and through the former name, directly
  CPPUNIT_ASSERT_EQUAL(2, oldR.calls);

  // new override, reached through the former name
  NewScleronomousR newR;
  static_cast<LagrangianScleronomousR&>(newR).computeJacobianhOver_q(q);
  CPPUNIT_ASSERT_EQUAL(1, newR.calls);
  CPPUNIT_ASSERT_EQUAL(2., newR.contactMatrix()(0, 1));

  // no override: the user function is called once, through either name
  WithH<LagrangianScleronomousR> plain;
  int fcalls = 0;
  plain.setComputeJacobianhOver_qFunction(
      [&fcalls](const BlockVector& q, Eigen::Ref<siconos::algebra::MapType> H) {
        ++fcalls;
        H(0, 2) = q(2);
      });
  plain.computeContactMatrix(q);
  CPPUNIT_ASSERT_EQUAL(1, fcalls);
  plain.computeJacobianhOver_q(q);
  CPPUNIT_ASSERT_EQUAL(2, fcalls);
  CPPUNIT_ASSERT_EQUAL(3., plain.jacobianhOver_q()(0, 2));
  CPPUNIT_ASSERT(plain.jacobianhOver_q() == plain.contactMatrix());
  CPPUNIT_ASSERT(plain.constraintJacobian() == plain.contactMatrix());
}

void RelationCompatibilityTest::testLagrangianRheonomousR() {
  BlockVector q;
  fillQ(q);
  OldRheonomousR rel;
  rel.computeContactMatrix(q, 0.5);
  CPPUNIT_ASSERT_EQUAL(1, rel.callsH);
  CPPUNIT_ASSERT_EQUAL(0.5, rel.contactMatrix()(0, 2));
  rel.computedhdt(q, 0.5);
  CPPUNIT_ASSERT_EQUAL(1, rel.callsDt);
  CPPUNIT_ASSERT_EQUAL(1., rel.dhdt()(0));
  CPPUNIT_ASSERT_EQUAL(1., rel.hdot()(0));
}

void RelationCompatibilityTest::testLagrangianCompliantR() {
  BlockVector q;
  fillQ(q);
  SiconosVector lambda(1);
  lambda << 4.;
  OldCompliantR rel;
  rel.computeContactMatrix(q, lambda);
  rel.computeD(q, lambda);
  CPPUNIT_ASSERT_EQUAL(1, rel.callsH);
  CPPUNIT_ASSERT_EQUAL(1, rel.callsD);
  CPPUNIT_ASSERT_EQUAL(4., rel.D()(0, 0));
  // the Relation interface returns D
  CPPUNIT_ASSERT_EQUAL(4., static_cast<const Relation&>(rel).jacobianhOver_lambda()(0, 0));
}

void RelationCompatibilityTest::testNewtonEulerR() {
  BlockVector q;
  auto q1 = std::make_shared<SiconosVector>(7);
  *q1 << 0., 0., 0., 1., 0., 0., 0.;
  q.insertPtr(q1);

  auto oldR = std::make_shared<OldNewtonEulerR>();
  Interaction inter(std::make_shared<NewtonImpactNSL>(0.), oldR);
  oldR->callComputeG(inter, q);
  oldR->computeContactMatrix(inter, q);
  CPPUNIT_ASSERT_EQUAL(1, oldR->callsG);
  CPPUNIT_ASSERT_EQUAL(1, oldR->callsH);
  CPPUNIT_ASSERT_EQUAL(5., oldR->contactMatrix()(0, 0));
  CPPUNIT_ASSERT(oldR->H_NE_prod_T() == oldR->contactMatrix());
  CPPUNIT_ASSERT(oldR->H_NE() == oldR->constraintJacobian());

  auto newR = std::make_shared<NewNewtonEulerR>();
  Interaction inter2(std::make_shared<NewtonImpactNSL>(0.), newR);
  static_cast<NewtonEulerR&>(*newR).computeH_NE_prod_T(inter2, q);
  CPPUNIT_ASSERT_EQUAL(1, newR->callsH);

  // no override: H = G T
  auto plain = std::make_shared<WithGH>();
  Interaction inter3(std::make_shared<NewtonImpactNSL>(0.), plain);
  SiconosMatrix G(1, 7);
  G << 1., 2., 3., 0., 0., 0., 0.;
  plain->setConstantH_NE(G);
  plain->computeContactMatrix(inter3, q);
  CPPUNIT_ASSERT_EQUAL(1., plain->contactMatrix()(0, 0));
  CPPUNIT_ASSERT_EQUAL(3., plain->contactMatrix()(0, 2));
}

void RelationCompatibilityTest::testClassAliases() {
  static_assert(std::is_same_v<NewtonEuler1DR, NewtonEulerNormalVelocityR>);
  static_assert(std::is_same_v<NewtonEuler3DR, NewtonEulerVelocityR>);
  static_assert(std::is_same_v<Lagrangian2d1DR, Lagrangian2dNormalVelocityR>);
  static_assert(std::is_same_v<Lagrangian2d2DR, Lagrangian2dVelocityR>);
  static_assert(std::is_same_v<Lagrangian2d3DR, Lagrangian2dVelocityAngularVelocityR>);
  auto rel = std::make_shared<NewtonEuler3DR>();
  CPPUNIT_ASSERT(rel != nullptr);
}
