// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "SiconosVisitorTest.hpp"

#include "ComplementarityConditionNSL.hpp"
#include "FremondImpactFrictionNSL.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "Question.hpp"

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(SiconosVisitorTest);

void SiconosVisitorTest::setUp() {}

void SiconosVisitorTest::tearDown() {}

struct ForMu : public siconos::modeling::nonsmooth_laws::Question<double> {
  void visit(const siconos::modeling::NewtonImpactFrictionNSL& nsl) override {
    answer = nsl.mu();
  }
  void visit(const siconos::modeling::FremondImpactFrictionNSL& nsl) override {
    answer = nsl.mu();
  }

  void visit(const siconos::modeling::ComplementarityConditionNSL& nsl) override {
    answer = 0.;
  }
};

void SiconosVisitorTest::testVisitNSL() {
  auto e = 0.3;
  auto mu = 0.4;
  auto nslaw0 = std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(e, e, mu, 3);
  auto nslaw1 = std::make_shared<siconos::modeling::ComplementarityConditionNSL>(4);
  auto nslaw3 = std::make_shared<siconos::modeling::FremondImpactFrictionNSL>(e, e, mu, 3);

  auto muVal = siconos::modeling::nonsmooth_laws::ask<ForMu>(*nslaw0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - visit a nonsmooth law:", muVal == mu, true);
  auto muVal2 = siconos::modeling::nonsmooth_laws::ask<ForMu>(*nslaw3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - visit a nonsmooth law:", muVal2 == mu, true);
  auto muVal3 = siconos::modeling::nonsmooth_laws::ask<ForMu>(*nslaw3);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("test - visit a nonsmooth law:", muVal2 == mu, true);
}
