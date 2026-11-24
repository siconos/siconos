/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __AVITest__
#define __AVITest__

#include <cppunit/extensions/HelperMacros.h>

#include "FirstOrderLinearDS.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "TimeDiscretisation.hpp"
#include "TimeStepping.hpp"

class AVITest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(AVITest);

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(AVITest);

  // tests to be done ...

#ifdef HAS_EXTREME_POINT_ALGO
  CPPUNIT_TEST(testAVI);
#endif

  CPPUNIT_TEST_SUITE_END();

  void testAVI();

  unsigned int _n{2};
  double _h{0.1};
  double _t0{0.};
  double _T{10.};
  double _tol{1.e-12};
  double _theta{0.5};
  std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> _nsds{nullptr};
  std::shared_ptr<siconos::simulation::TimeStepping> _sim{nullptr};
  std::shared_ptr<siconos::modeling::FirstOrderLinearDS> _DS{nullptr};
  std::shared_ptr<siconos::simulation::TimeDiscretisation> _TD{nullptr};
  std::shared_ptr<siconos::integrators::OneStepIntegrator> _osi{nullptr};
  std::shared_ptr<siconos::algebra::SiconosMatrix> _A{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _b{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _x0{nullptr};

 public:
  AVITest() = default;
  void setUp();
  void tearDown();
};

#endif
