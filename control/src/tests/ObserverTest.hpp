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
#ifndef __ObserverTest__
#define __ObserverTest__

#include <cppunit/extensions/HelperMacros.h>
#include <SiconosFwd.hpp>
#include "SiconosControlFwd.hpp"
#include <FirstOrderLinearTIDS.hpp>

class ObserverTest : public CppUnit::TestFixture
{

private:
  ACCEPT_SERIALIZATION(ObserverTest);


  // Name of the tests suite
  CPPUNIT_TEST_SUITE(ObserverTest);

  // tests to be done ...

  CPPUNIT_TEST(test_SMO_ZOH);
  CPPUNIT_TEST(test_SMO_Lsodar);
  CPPUNIT_TEST(test_Luenberger_ZOH);
  CPPUNIT_TEST(test_Luenberger_Lsodar);

  CPPUNIT_TEST_SUITE_END();

  void init();
  void init2();
  void test_SMO_ZOH();
  void test_SMO_Lsodar();
  void test_Luenberger_ZOH();
  void test_Luenberger_Lsodar();
  // Members

  unsigned int _n;
  double _h;
  double _t0;
  double _T;
  double _tol;
  double _xFinal;
  SP::FirstOrderLinearTIDS _DS;
  SP::SiconosMatrix _A;
  SP::SimpleMatrix _B;
  SP::SimpleMatrix _C;
  SP::SimpleMatrix _Csurface;
  SP::SiconosVector _b;
  SP::SiconosVector _x0;
  SP::SiconosVector _xHat0;
  SP::SiconosVector _K;
  SP::SimpleMatrix _L;
  SP::LinearSensor _sensor;
  SP::PID _pid;


public:

  ObserverTest(): _n(2), _h(0.05), _t0(0.0), _T(100.0), _tol(7e-11) {}
  void setUp();
  void tearDown();

};

#endif




