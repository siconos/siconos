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
#include "TwistingTest.hpp"

#include <FirstOrderLinearTIDS.hpp>
#include <io.hpp>

#include "ControlLsodarSimulation.hpp"
#include "ControlZOHSimulation.hpp"
#include "ExplicitTwisting.hpp"
#include "LinearSensor.hpp"
#include "RegularTwisting.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "Twisting.hpp"

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
  (*_B)(1, 0) = 1.0;

  _Csurface = std::make_shared<siconos::algebra::SiconosMatrix>(1, 2);
  _Csurface->setZero();
  (*_Csurface)(0, 0) = 1.0;
  (*_Csurface)(0, 1) = 1.0;
}

#ifdef HAS_EXTREME_POINT_ALGO

void TwistingTest::initTwisting() {
  _DS = std::make_shared<siconos::modeling::FirstOrderLinearTIDS>(_x0, _A);
  _sensor = std::make_shared<siconos::control::LinearSensor>(_DS, _C);
  _itw = std::make_shared<siconos::control::Twisting>(_sensor, 300., _beta, _h);
  auto eye = std::make_shared<siconos::algebra::SiconosMatrix>(2, 2);
  eye->setIdentity();
  _itw->setCsurface(eye);
}

void TwistingTest::initRegularTwisting() {
  _DS = std::make_shared<siconos::modeling::FirstOrderLinearTIDS>(_x0, _A);
  _sensor = std::make_shared<siconos::control::LinearSensor>(_DS, _C);
  _reg_itw = std::make_shared<siconos::control::RegularTwisting>(_sensor, 300., _beta);
  auto eye = std::make_shared<siconos::algebra::SiconosMatrix>(2, 2);
  eye->setIdentity();
  _reg_itw->setCsurface(eye);
}
#endif

void TwistingTest::initExplicitTwisting() {
  _DS = std::make_shared<siconos::modeling::FirstOrderLinearTIDS>(_x0, _A);
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
  double error = 0.0;
  bool test =
      !((error = siconos::algebra::io::compareRefFile(data, "etw_ZOH.ref", _tol)) >= 0.0 &&
        error > _tol);
  std::cout << "------- Integration done -------" << test << std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_Luenberger_ZOH : ", test, true);
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
  auto error = siconos::algebra::io::compareRefFile(data, "etw_lsodar.ref", _tol);
  bool test = !(error > _tol);
  std::cout << "------- Integration done -------" << test << std::endl;
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
  siconos::algebra::io::write("itw_ZOH.dat", *data, siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  // Reference Matrix
  siconos::algebra::SiconosMatrix dataRef(data->rows(), data->cols());
  siconos::algebra::io::read("itw2.ref", dataRef);
  // it is a bad idea to compare solutions to an AVI that does not admit a unique
  // solution
  data->col(3) = _beta * data->col(2) + data->col(1);
  dataRef.col(3) = _beta * dataRef.col(2) + dataRef.col(1);
  dataRef -= *data;

  std::cout << "------- Integration done, error = " << dataRef.normInf() << " -------"
            << std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_itw_ZOH : ", dataRef.normInf() < _tol, true);
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
  siconos::algebra::io::write("itw_Lsodar.dat", *data, siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  // Reference Matrix
  siconos::algebra::SiconosMatrix dataRef(data->rows(), data->cols());
  siconos::algebra::io::read("itw2.ref", dataRef);
  // it is a bad idea to compare solutions to an AVI that does not admit a unique
  // solution
  data->col(3) = _beta * data->col(2) + data->col(1);
  dataRef.col(3) = _beta * dataRef.col(2) + dataRef.col(1);
  dataRef -= *data;
  std::cout << "------- Integration done, error = " << dataRef.normInf() << " -------\n";
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_itw_Lsodar : ", dataRef.normInf() < _tol, true);
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
  siconos::algebra::SiconosMatrix dataRef(data->rows(), data->cols());
  siconos::algebra::io::read("reg_itw.ref", dataRef);
  // it is a bad idea to compare solutions to an AVI that does not admit a
  // unique solution
  data->col(3) = _beta * data->col(2) + data->col(1);
  dataRef.col(3) = _beta * dataRef.col(2) + dataRef.col(1);
  dataRef -= *data;
  std::cout << "------- Integration done, error = " << dataRef.normInf() << " -------\n";
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_RegularTwistingZOH : ", dataRef.normInf() < _tol, true);
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
  siconos::algebra::io::write("reg_itw_Lsodar..dat", *data, siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  // Reference Matrix
  siconos::algebra::SiconosMatrix dataRef(data->rows(), data->cols());
  siconos::algebra::io::read("reg_itw.ref", dataRef);
  // it is a bad idea to compare solutions to an AVI that does not admit a unique
  // solution
  data->col(3) = _beta * data->col(2) + data->col(1);
  dataRef.col(3) = _beta * dataRef.col(2) + dataRef.col(1);
  dataRef -= *data;
  std::cout << "------- Integration done, error = " << dataRef.normInf() << " -------\n";
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_RegularTwistingLsodar : ", dataRef.normInf() < 5e-9,
                               true);
}
#endif
