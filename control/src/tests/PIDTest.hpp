/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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
#ifndef __PIDTest__
#define __PIDTest__

#include <cppunit/extensions/HelperMacros.h>
#include <SiconosFwd.hpp>
#include "SiconosControlFwd.hpp"
#include <FirstOrderLinearTIDS.hpp>

class PIDTest : public CppUnit::TestFixture
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(PIDTest);


  // Name of the tests suite
  CPPUNIT_TEST_SUITE(PIDTest);

  // tests to be done ...

  CPPUNIT_TEST(testPIDZOH);
  CPPUNIT_TEST(testPIDLsodar);

  CPPUNIT_TEST_SUITE_END();

  void init();
  void testPIDZOH();
  void testPIDLsodar();
  // Members

  unsigned int _n;
  double _h;
  double _t0;
  double _T;
  double _tol;
  double _xFinal;
  SP::FirstOrderLinearTIDS _DS;
  SP::SiconosMatrix _A;
  SP::SiconosVector _b;
  SP::SiconosVector _x0;
  SP::SiconosVector _K;
  SP::LinearSensor _sensor;
  SP::PID _PIDcontroller;



public:

  PIDTest(): _n(2), _h(0.05), _t0(0.0), _T(100.0), _tol(5e-12) {}
  void setUp();
  void tearDown();

};

#endif




