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
#ifndef __SMCTest__
#define __SMCTest__

#include <cppunit/extensions/HelperMacros.h>

#include "ExplicitLinearSMC.hpp"
#include "FirstOrderLinearDS.hpp"
#include "LinearSMC.hpp"
#include "LinearSensor.hpp"
#include "SiconosConfig.h"
#include "Twisting.hpp"

class SMCTest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(SMCTest);

  // Name of the tests suite
  CPPUNIT_TEST_SUITE(SMCTest);

  // tests to be done ...

  CPPUNIT_TEST(test_iSMC_ZOH);
  CPPUNIT_TEST(test_iSMC_Lsodar);
  CPPUNIT_TEST(test_eSMC_ZOH);
  CPPUNIT_TEST(test_eSMC_Lsodar);
#ifdef HAS_EXTREME_POINT_ALGO
  CPPUNIT_TEST(test_itw_ZOH);
  CPPUNIT_TEST(test_itw_Lsodar);
#endif

  CPPUNIT_TEST_SUITE_END();

  void init();
  void init2();
  void initTwisting();
  void test_iSMC_ZOH();
  void test_iSMC_Lsodar();
  void test_eSMC_ZOH();
  void test_eSMC_Lsodar();

#ifdef HAS_EXTREME_POINT_ALGO
  void test_itw_ZOH();
  void test_itw_Lsodar();
#endif
  // Members

  unsigned int _n;
  double _h;
  double _t0;
  double _T;
  double _tol;
  double _beta;
  double _xFinal;
  std::shared_ptr<siconos::modeling::FirstOrderLinearDS> _DS;
  std::shared_ptr<siconos::algebra::SiconosMatrix> _A;
  std::shared_ptr<siconos::algebra::SiconosMatrix> _B;
  std::shared_ptr<siconos::algebra::SiconosMatrix> _C;
  std::shared_ptr<siconos::algebra::SiconosMatrix> _Csurface;
  std::shared_ptr<siconos::algebra::SiconosVector> _b;
  std::shared_ptr<siconos::algebra::SiconosVector> _x0;
  std::shared_ptr<siconos::algebra::SiconosVector> _K;
  std::shared_ptr<siconos::control::LinearSensor> _sensor;
  std::shared_ptr<siconos::control::LinearSMC> _iSMC;
  std::shared_ptr<siconos::control::ExplicitLinearSMC> _eSMC;
#ifdef HAS_EXTREME_POINT_ALGO
  std::shared_ptr<siconos::control::Twisting> _itw;
#endif

 public:
  SMCTest() : _n(2), _h(0.05), _t0(0.0), _T(100.0), _tol(1.5e-8), _beta(0.1) {}
  void setUp();
  void tearDown(){};
};

#endif
