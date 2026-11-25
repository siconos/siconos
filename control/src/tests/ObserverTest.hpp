// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef __ObserverTest__
#define __ObserverTest__

#include <cppunit/extensions/HelperMacros.h>

#include <FirstOrderLinearDS.hpp>

#include "LinearSensor.hpp"
#include "PID.hpp"

class ObserverTest : public CppUnit::TestFixture {
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

  size_t _n{2};
  double _h{0.05};
  double _t0{0.};
  double _T{100};
  /** relative tolerance of the comparisons with the references (times max |data|): the ZOH
   *  matrices are computed with Lsodar at tolerances of a few machine epsilons, and the
   *  round-off differences between platforms (macOS arm64 / Linux x86_64) reach 8e-11 on data
   *  of size 400 (the former absolute 7e-11 only held on the platform of the references) */
  double _tol{1e-11};
  double _beta{0.2};

  std::shared_ptr<siconos::modeling::FirstOrderLinearDS> _DS{nullptr};
  std::shared_ptr<siconos::algebra::SiconosMatrix> _A{nullptr};
  std::shared_ptr<siconos::algebra::SiconosMatrix> _B{nullptr};
  std::shared_ptr<siconos::algebra::SiconosMatrix> _C{nullptr};
  std::shared_ptr<siconos::algebra::SiconosMatrix> _Csurface{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _b{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _x0{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _xHat0{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> _K{nullptr};
  std::shared_ptr<siconos::algebra::SiconosMatrix> _L{nullptr};
  std::shared_ptr<siconos::control::LinearSensor> _sensor{nullptr};
  std::shared_ptr<siconos::control::PID> _pid{nullptr};

 public:
  ObserverTest() = default;
  void setUp();
  void tearDown();
};

#endif
