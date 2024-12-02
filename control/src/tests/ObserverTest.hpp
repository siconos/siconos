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
#ifndef __ObserverTest__
#define __ObserverTest__

#include <cppunit/extensions/HelperMacros.h>
#include <FirstOrderLinearTIDS.hpp>
#include "LinearSensor.hpp"
#include "PID.hpp"


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
  std::shared_ptr<siconos::modeling::FirstOrderLinearTIDS> _DS;
  std::shared_ptr<siconos::algebra::SiconosMatrix> _A{nullptr};
  std::shared_ptr<siconos::algebra::SimpleMatrix> _B{nullptr};
  std::shared_ptr<siconos::algebra::SimpleMatrix> _C{nullptr};
  std::shared_ptr<siconos::algebra::SimpleMatrix> _Csurface{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _b{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _x0{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _xHat0{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _K{nullptr};
  std::shared_ptr<siconos::algebra::SimpleMatrix> _L{nullptr};
  std::shared_ptr<siconos::control::LinearSensor> _sensor{nullptr};
  std::shared_ptr<siconos::control::PID> _pid{nullptr};


public:

  ObserverTest(): _n(2), _h(0.05), _t0(0.0), _T(100.0), _tol(7e-11) {}
  void setUp();
  void tearDown();

};

#endif




