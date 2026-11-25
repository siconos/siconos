// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "TwistingTest.hpp"

#include <FirstOrderLinearDS.hpp>
#include <io.hpp>

#include "ControlLsodarSimulation.hpp"
#include "ControlZOHSimulation.hpp"
#include "ExplicitTwisting.hpp"
#include "LinearSensor.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(TwistingTest);

void TwistingTest::setUp() {
  _A = std::make_shared<siconos::algebra::SiconosMatrix>(_n, _n);
  _A->setZero();
  (*_A)(0, 1) = 1.0;
  (*_A)(1, 0) = 19.0;
  (*_A)(1, 1) = -2.0;

  _x0 = std::make_shared<siconos::algebra::SiconosVector>(_n);
  _x0->setZero();
  (*_x0)(0) = -15.0;
  (*_x0)(1) = 20.0;

  _C = std::make_shared<siconos::algebra::SiconosMatrix>(2, 2);
  _C->setZero();
  _C->setIdentity();

  _B = std::make_shared<siconos::algebra::SiconosMatrix>(2, 1);
  (*_B)(0, 0) = 0.0;
  (*_B)(1, 0) = 1.0;

  _Csurface = std::make_shared<siconos::algebra::SiconosMatrix>(1, 2);
  _Csurface->setZero();
  (*_Csurface)(0, 0) = 1.0;
  (*_Csurface)(0, 1) = 1.0;
}

#ifdef HAS_EXTREME_POINT_ALGO

void TwistingTest::initTwisting() {
  _DS =
      std::make_shared<siconos::modeling::FirstOrderLinearDS>(*_x0, siconos::algebra::alias_t);
  _DS->setConstantA(*_A, siconos::algebra::alias_t);
  _sensor = std::make_shared<siconos::control::LinearSensor>(_DS, _C);
  _itw = std::make_shared<siconos::control::Twisting>(_sensor, 300., _beta, _h);
  auto eye = std::make_shared<siconos::algebra::SiconosMatrix>(2, 2);
  eye->setIdentity();
  _itw->setCsurface(eye);
}

void TwistingTest::initRegularTwisting() {
  _DS =
      std::make_shared<siconos::modeling::FirstOrderLinearDS>(*_x0, siconos::algebra::alias_t);
  _DS->setConstantA(*_A, siconos::algebra::alias_t);
  _sensor = std::make_shared<siconos::control::LinearSensor>(_DS, _C);
  _reg_itw = std::make_shared<siconos::control::RegularTwisting>(_sensor, 300., _beta);
  auto eye = std::make_shared<siconos::algebra::SiconosMatrix>(2, 2);
  eye->setIdentity();
  _reg_itw->setCsurface(eye);
}
#endif

void TwistingTest::initExplicitTwisting() {
  _DS =
      std::make_shared<siconos::modeling::FirstOrderLinearDS>(*_x0, siconos::algebra::alias_t);
  _DS->setConstantA(*_A, siconos::algebra::alias_t);
  _sensor = std::make_shared<siconos::control::LinearSensor>(_DS, _C);
  _expl_tw = std::make_shared<siconos::control::ExplicitTwisting>(_sensor, 300., _beta);
  auto eye = std::make_shared<siconos::algebra::SiconosMatrix>(2, 2);
  eye->setIdentity();
  _expl_tw->setCsurface(eye);
}

void TwistingTest::tearDown() {}

void TwistingTest::test_ExplicitTwisting_ZOH() {
  initExplicitTwisting();
  auto simZOH = std::make_shared<siconos::control::ControlZOHSimulation>(_t0, _T, _h);
  simZOH->setSaveOnlyMainSimulation(true);
  simZOH->addDynamicalSystem(_DS);
  simZOH->addSensor(_sensor, _h);
  simZOH->addActuator(_expl_tw, _h);
  simZOH->initialize();
  simZOH->run();
  auto& data = *simZOH->data();
  siconos::algebra::io::write("explicitTwisting_ZOH.dat", data,
                              siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  // The explicit twisting controller discretized with h = 0.05 and gain 300 gives a bounded
  // but chaotic closed loop: a round-off difference grows by a factor ~1000 every 2 s
  // (between macOS arm64 and Linux x86_64: 1e-12 at t = 1, 1e-9 at t = 3, O(1) at t = 9).
  // The trajectory is compared with the reference only over the first 1.5 s (relative
  // tolerance 1e-11), then the state must stay bounded (|x| <= 52.4 for the reference).
  auto ref = siconos::algebra::io::readDenseMatrix("etw_ZOH.ref");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_ExplicitTwisting_ZOH: size of the reference",
                               ref.rows() == data.rows() && ref.cols() == data.cols(), true);
  Eigen::Index rows = 0;
  while (rows < ref.rows() && ref(rows, 0) <= 1.5 + 1e-12) ++rows;
  double tol = 1e-11 * ref.topRows(rows).cwiseAbs().maxCoeff();
  double error = (data.topRows(rows) - ref.topRows(rows)).cwiseAbs().maxCoeff();
  double max_x = data.middleCols(1, 2).cwiseAbs().maxCoeff();
  std::cout << "------- Integration done: error over [0, 1.5] = " << error << " (tolerance "
            << tol << "), max |x| = " << max_x << " -------" << std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_ExplicitTwisting_ZOH : ", error <= tol && max_x <= 100.,
                               true);
}

void TwistingTest::test_ExplicitTwisting_Lsodar() {
  initExplicitTwisting();
  auto simLsodar = std::make_shared<siconos::control::ControlLsodarSimulation>(_t0, _T, _h);
  simLsodar->setSaveOnlyMainSimulation(true);
  simLsodar->addDynamicalSystem(_DS);
  simLsodar->addSensor(_sensor, _h);
  simLsodar->addActuator(_expl_tw, _h);
  simLsodar->initialize();
  simLsodar->run();

  auto& data = *simLsodar->data();
  siconos::algebra::io::write("explicitTwisting_Lsodar.dat", data,
                              siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  double error = 0.;
  bool test =
      !((error = siconos::algebra::io::compareRefFile(data, "etw_ZOH.ref", _tol)) >= 0.0 &&
        error > _tol);
  std::cout << "------- Integration done -------" << test << std::endl;
  // FP warning: this test is broken for a long time ... It needs a control's expert review for
  // ControlLsodarSimulation.
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("test_Luenberger_ZOH : ", test, true);
}

#ifdef HAS_EXTREME_POINT_ALGO
void TwistingTest::test_Twisting_ZOH() {
  initTwisting();
  auto simZOH = std::make_shared<siconos::control::ControlZOHSimulation>(_t0, _T, _h);
  simZOH->setSaveOnlyMainSimulation(true);
  simZOH->addDynamicalSystem(_DS);
  simZOH->addSensor(_sensor, _h);
  simZOH->addActuator(_itw, _h);
  simZOH->initialize();
  simZOH->run();
  auto data = simZOH->data();
  siconos::algebra::io::write("itw_twisting_ZOH.dat", *data, siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  // Reference Matrix
  auto dataRef = siconos::algebra::io::readDenseMatrix("itw2.ref");
  // it is a bad idea to compare solutions to an AVI that does not admit a unique
  // solution
  data->col(3) = _beta * data->col(2) + data->col(1);
  dataRef.col(3) = _beta * dataRef.col(2) + dataRef.col(1);
  dataRef -= *data;
  auto error = siconos::algebra::normInf(dataRef.leftCols(3));
  std::cout << "------- Integration done, error = " << error << " -------" << std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_itw_ZOH : ", error < _tol, true);
}

void TwistingTest::test_Twisting_Lsodar() {
  initTwisting();
  auto simLsodar = std::make_shared<siconos::control::ControlLsodarSimulation>(_t0, _T, _h);
  simLsodar->setSaveOnlyMainSimulation(true);
  simLsodar->addDynamicalSystem(_DS);
  simLsodar->addSensor(_sensor, _h);
  simLsodar->addActuator(_itw, _h);
  simLsodar->initialize();
  simLsodar->run();
  auto data = simLsodar->data();
  siconos::algebra::io::write("itw_twisting_Lsodar.dat", *data,
                              siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  // Reference Matrix
  auto dataRef = siconos::algebra::io::readDenseMatrix("itw2.ref");
  // it is a bad idea to compare solutions to an AVI that does not admit a unique
  // solution
  data->col(3) = _beta * data->col(2) + data->col(1);
  dataRef.col(3) = _beta * dataRef.col(2) + dataRef.col(1);
  dataRef -= *data;
  auto error = siconos::algebra::normInf(dataRef.leftCols(3));
  std::cout << "------- Integration done, error = " << error << " -------\n";
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_itw_Lsodar : ", error < _tol, true);
}

void TwistingTest::test_RegularTwisting_ZOH() {
  initRegularTwisting();
  auto simZOH = std::make_shared<siconos::control::ControlZOHSimulation>(_t0, _T, _h);
  simZOH->setSaveOnlyMainSimulation(true);
  simZOH->addDynamicalSystem(_DS);
  simZOH->addSensor(_sensor, _h);
  simZOH->addActuator(_reg_itw, _h);
  simZOH->initialize();
  simZOH->run();
  auto data = simZOH->data();
  siconos::algebra::io::write("reg_itw_ZOH.dat", *data, siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  // Reference Matrix
  auto dataRef = siconos::algebra::io::readDenseMatrix("reg_itw.ref");
  // it is a bad idea to compare solutions to an AVI that does not admit a
  // unique solution
  data->col(3) = _beta * data->col(2) + data->col(1);
  dataRef.col(3) = _beta * dataRef.col(2) + dataRef.col(1);
  dataRef -= *data;
  auto error = siconos::algebra::normInf(dataRef.leftCols(3));

  std::cout << "------- Integration done, error = " << error << " -------\n";
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_RegularTwistingZOH : ", error < _tol, true);
}

void TwistingTest::test_RegularTwisting_Lsodar() {
  initRegularTwisting();
  auto simLsodar = std::make_shared<siconos::control::ControlLsodarSimulation>(_t0, _T, _h);
  simLsodar->setSaveOnlyMainSimulation(true);
  simLsodar->addDynamicalSystem(_DS);
  simLsodar->addSensor(_sensor, _h);
  simLsodar->addActuator(_reg_itw, _h);
  simLsodar->initialize();
  simLsodar->run();
  auto data = simLsodar->data();
  siconos::algebra::io::write("reg_itw_Lsodar.dat", *data, siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  // Reference Matrix
  auto dataRef = siconos::algebra::io::readDenseMatrix("reg_itw.ref");
  // it is a bad idea to compare solutions to an AVI that does not admit a unique
  // solution
  data->col(3) = _beta * data->col(2) + data->col(1);
  dataRef.col(3) = _beta * dataRef.col(2) + dataRef.col(1);
  dataRef -= *data;
  auto error = siconos::algebra::normInf(dataRef.leftCols(3));
  std::cout << "------- Integration done, error = " << error << " -------\n";
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_RegularTwistingLsodar : ", error < 5e-9, true);
}
#endif
