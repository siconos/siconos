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
#ifndef __TwistingTest__
#define __TwistingTest__

#include <cppunit/extensions/HelperMacros.h>
#include "SiconosFwd.hpp"
#include "SiconosControlFwd.hpp"
#include <FirstOrderLinearTIDS.hpp>

#include <SiconosConfig.h>

class TwistingTest : public CppUnit::TestFixture
{

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

  unsigned int _n;
  double _h;
  double _t0;
  double _T;
  double _tol;
  double _beta;
  double _xFinal;
  SP::FirstOrderLinearTIDS _DS;
  SP::SiconosMatrix _A;
  SP::SimpleMatrix _B;
  SP::SimpleMatrix _C;
  SP::SimpleMatrix _Csurface;
  SP::SiconosVector _b;
  SP::SiconosVector _x0;
  SP::SiconosVector _K;
  SP::LinearSensor _sensor;
#ifdef HAS_EXTREME_POINT_ALGO
  SP::Twisting _itw;
  SP::RegularTwisting _reg_itw;
#endif
  SP::ExplicitTwisting _expl_tw;


public:

  TwistingTest(): _n(2), _h(0.05), _t0(0.0), _T(100.0), _tol(7.5e-11), _beta(0.2) {}
  void setUp();
  void tearDown();

};

#endif




