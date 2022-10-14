/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#ifndef __ZOHTest__
#define __ZOHTest__

#include <cppunit/extensions/HelperMacros.h>

#include "FirstOrderLinearTIDS.hpp"
#include "FirstOrderLinearTIR.hpp"
#include "Interaction.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Relay.hpp"
#include "RelayNSL.hpp"
#include "TimeDiscretisation.hpp"
#include "TimeStepping.hpp"
#include "ZeroOrderHoldOSI.hpp"

class ZOHTest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(ZOHTest);

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(ZOHTest);

  // tests to be done ...

  CPPUNIT_TEST(testMatrixExp0);
  CPPUNIT_TEST(testMatrixExp1);
  CPPUNIT_TEST(testMatrixIntegration1);
  CPPUNIT_TEST(testMatrixIntegration2);
  CPPUNIT_TEST(testMatrixIntegration3);
  CPPUNIT_TEST(testMatrixIntegration4);

  CPPUNIT_TEST_SUITE_END();

  // \todo exception test
  void init();
  void testMatrixExp0();
  void testMatrixExp1();
  void testMatrixIntegration1();
  void testMatrixIntegration2();
  void testMatrixIntegration3();
  void testMatrixIntegration4();
  // Members

  unsigned int _n{2};
  double _h{0.1};
  double _t0{0.};
  double _T{10.};
  double _tol{1.e-12};
  std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> _model{nullptr};
  std::shared_ptr<siconos::simulation::TimeStepping> _sim{nullptr};
  std::shared_ptr<siconos::modeling::FirstOrderLinearDS> _DS{nullptr};
  std::shared_ptr<siconos::simulation::TimeDiscretisation> _TD{nullptr};
  std::shared_ptr<siconos::integrators::ZeroOrderHoldOSI> _ZOH;
  std::shared_ptr<siconos::algebra::SiconosMatrix> _A{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _b{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _x0{nullptr};

 public:
  ZOHTest() = default;
  void setUp();
  void tearDown();
};

#endif
