// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef __EulerMoreauTest__
#define __EulerMoreauTest__

#include <cppunit/extensions/HelperMacros.h>

#include "EulerMoreauOSI.hpp"
#include "FirstOrderLinearDS.hpp"
#include "FirstOrderLinearTIR.hpp"
#include "Interaction.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Relay.hpp"
#include "RelayNSL.hpp"
#include "TimeDiscretisation.hpp"
#include "TimeStepping.hpp"

class EulerMoreauTest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(EulerMoreauTest);

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(EulerMoreauTest);

  // tests to be done ...

  CPPUNIT_TEST(testCstGradTIDS);
  CPPUNIT_TEST(testCstGradDS);
  CPPUNIT_TEST(testCstGradNLDS);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test
  void init(bool initDS);
  void testCstGradTIDS();
  void testCstGradDS();
  void testCstGradNLDS();
  // Members

  unsigned int _n{2};
  double _h{0.001};
  double _t0{0.};
  double _T{3.};
  double _tol{1.e-6};
  std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> _model{nullptr};
  std::shared_ptr<siconos::simulation::TimeStepping> _sim{nullptr};
  std::shared_ptr<siconos::modeling::DynamicalSystem> _DS{nullptr};
  std::shared_ptr<siconos::simulation::TimeDiscretisation> _TD{nullptr};
  std::shared_ptr<siconos::integrators::EulerMoreauOSI> _EulerMoreau{nullptr};
  std::shared_ptr<siconos::algebra::SiconosMatrix> _A{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _b{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _x0{nullptr};

 public:
  EulerMoreauTest() = default;
  void setUp();
  void tearDown();
};

#endif
