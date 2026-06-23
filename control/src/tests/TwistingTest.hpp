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
#ifndef __TwistingTest__
#define __TwistingTest__

#include <cppunit/extensions/HelperMacros.h>

#include <FirstOrderLinearDS.hpp>

#include "ExplicitTwisting.hpp"
#include "LinearSensor.hpp"

class TwistingTest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(TwistingTest);

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(TwistingTest);

  // tests to be done ...

#ifdef HAS_EXTREME_POINT_ALGO
  CPPUNIT_TEST(test_Twisting_ZOH);
  CPPUNIT_TEST(test_Twisting_Lsodar);
  CPPUNIT_TEST(test_RegularTwisting_ZOH);
  CPPUNIT_TEST(test_RegularTwisting_Lsodar);
#endif

  CPPUNIT_TEST(test_ExplicitTwisting_ZOH);
  CPPUNIT_TEST(test_ExplicitTwisting_Lsodar);

  CPPUNIT_TEST_SUITE_END();

  void initTwisting();
  void initRegularTwisting();
  void initExplicitTwisting();
  void test_ExplicitTwisting_ZOH();
  void test_ExplicitTwisting_Lsodar();

#ifdef HAS_EXTREME_POINT_ALGO
  void test_Twisting_ZOH();
  void test_Twisting_Lsodar();
  void test_RegularTwisting_ZOH();
  void test_RegularTwisting_Lsodar();
#endif
  // Members

  size_t _n{2};
  double _h{0.05};
  double _t0{0.};
  double _T{100};
  double _tol{7.5e-11};
  double _beta{0.2};
  std::shared_ptr<siconos::modeling::FirstOrderLinearDS> _DS{nullptr};
  std::shared_ptr<siconos::algebra::SiconosMatrix> _A{nullptr};
  std::shared_ptr<siconos::algebra::SiconosMatrix> _B{nullptr};
  std::shared_ptr<siconos::algebra::SiconosMatrix> _C{nullptr};
  std::shared_ptr<siconos::algebra::SiconosMatrix> _Csurface{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _b{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _x0{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _K{nullptr};
  std::shared_ptr<siconos::control::LinearSensor> _sensor{nullptr};
#ifdef HAS_EXTREME_POINT_ALGO
  std::shared_ptr<siconos::control::Twisting> _itw;
  std::shared_ptr<siconos::control::RegularTwisting> _reg_itw;
#endif
  std::shared_ptr<siconos::control::ExplicitTwisting> _expl_tw;

 public:
  TwistingTest() = default;
  void setUp();
  void tearDown();
};

#endif
