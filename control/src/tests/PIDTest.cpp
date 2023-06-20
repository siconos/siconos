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
#include "PIDTest.hpp"

#include "ControlLsodarSimulation.hpp"
#include "ControlZOHSimulation.hpp"
#include "FirstOrderLinearTIDS.hpp"
#include "LinearSensor.hpp"
#include "PID.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "io.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(PIDTest);

void PIDTest::setUp() {
  _A = std::make_shared<siconos::algebra::SimpleMatrix>(_n, _n);
  (*_A)(0, 1) = 1.0;
  _x0 = std::make_shared<siconos::algebra::SiconosVector>(_n, 0);
  (*_x0)(0) = 10.0;
  (*_x0)(1) = 0.0;
  _xFinal = 0.0;

  _K = std::make_shared<siconos::algebra::SiconosVector>(3, 0);
  (*_K)(0) = .25;
  (*_K)(1) = .125;
  (*_K)(2) = 2.0;
}

void PIDTest::init() {
  _DS = std::make_shared<siconos::modeling::FirstOrderLinearTIDS>(_x0, _A);
  auto C = std::make_shared<siconos::algebra::SimpleMatrix>(1, 2, 0);
  (*C)(0, 0) = 1;
  _sensor = std::make_shared<siconos::control::LinearSensor>(_DS, C);
  auto B = std::make_shared<siconos::algebra::SimpleMatrix>(2, 1);
  (*B)(1, 0) = 1;
  _PIDcontroller = std::make_shared<siconos::control::PID>(_sensor, B);
  _PIDcontroller->setRef(_xFinal);
  _PIDcontroller->setK(_K);
  _PIDcontroller->setDeltaT(_h);
}

void PIDTest::tearDown() {}

void PIDTest::testPIDZOH() {
  init();
  auto simZOH = std::make_shared<siconos::control::ControlZOHSimulation>(_t0, _T, _h);
  simZOH->addDynamicalSystem(_DS);
  simZOH->addSensor(_sensor, _h);
  simZOH->addActuator(_PIDcontroller, _h);
  simZOH->initialize();
  simZOH->run();
  auto& data = *simZOH->data();
  siconos::algebra::io::write("PIDZOH.dat", data, siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  // Reference Matrix
  siconos::algebra::SimpleMatrix dataRef(data);
  dataRef.zero();
  siconos::algebra::io::read("PID.ref", dataRef);

  std::cout << "------- Integration done, error = " << (data - dataRef).normInf() << " -------"
            << std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testPIDZOH : ", (data - dataRef).normInf() < _tol, true);
}

void PIDTest::testPIDLsodar() {
  init();
  auto simLsodar = std::make_shared<siconos::control::ControlLsodarSimulation>(_t0, _T, _h);
  simLsodar->addDynamicalSystem(_DS);
  simLsodar->addSensor(_sensor, _h);
  simLsodar->addActuator(_PIDcontroller, _h);
  simLsodar->initialize();
  simLsodar->run();
  auto& data = *simLsodar->data();
  siconos::algebra::io::write("PIDLsodar.dat", data, siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  // Reference Matrix
  siconos::algebra::SimpleMatrix dataRef(data);
  dataRef.zero();
  siconos::algebra::io::read("PID.ref", dataRef);
  std::cout << "------- Integration done, error = " << (data - dataRef).normInf() << " -------"
            << std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testPIDLsodar : ", (data - dataRef).normInf() < _tol, true);
}
