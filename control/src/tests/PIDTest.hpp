/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef __PIDTest__
#define __PIDTest__

#include <cppunit/extensions/HelperMacros.h>

#include <FirstOrderLinearDS.hpp>
#include <LinearSensor.hpp>
#include <PID.hpp>

class PIDTest : public CppUnit::TestFixture {
 private:
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

  size_t _n;
  double _h;
  double _t0;
  double _T;
  double _tol;
  double _xFinal{0.};
  std::shared_ptr<siconos::modeling::FirstOrderLinearDS> _DS;
  std::shared_ptr<siconos::algebra::SiconosMatrix> _A;
  std::shared_ptr<siconos::algebra::SiconosVector> _b;
  std::shared_ptr<siconos::algebra::SiconosVector> _x0;
  std::shared_ptr<siconos::algebra::SiconosVector> _K;
  std::shared_ptr<siconos::control::LinearSensor> _sensor;
  std::shared_ptr<siconos::control::PID> _PIDcontroller;

 public:
  PIDTest() : _n(2), _h(0.05), _t0(0.0), _T(100.0), _tol(5e-12) {}
  void setUp();
  void tearDown();
};

#endif
