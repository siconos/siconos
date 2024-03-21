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
#include "ObserverTest.hpp"

#include <FirstOrderLinearTIDS.hpp>

#include "ControlLsodarSimulation.hpp"
#include "ControlZOHSimulation.hpp"
#include "LinearSMC.hpp"
#include "LinearSensor.hpp"
#include "LuenbergerObserver.hpp"
#include "PID.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "SlidingReducedOrderObserver.hpp"
#include "io.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(ObserverTest);

void ObserverTest::setUp()
{
  _A = std::make_shared<siconos::algebra::SimpleMatrix>(_n, _n);
  _A->setZero();
  (*_A)(0, 1) = 1.0;
  (*_A)(1, 0) = -1.0;

  _x0 = std::make_shared<siconos::algebra::SiconosVector>(_n);
  _x0->setZero();
  (*_x0)(0) = 10.0;
  (*_x0)(1) = 0.0;

  _C = std::make_shared<siconos::algebra::SimpleMatrix>(1, 2);
  _C->setZero();
  (*_C)(0, 0) = 1.0;

  _B = std::make_shared<siconos::algebra::SimpleMatrix>(2, 1);
  (*_B)(1, 0) = 1.0;

  _Csurface = std::make_shared<siconos::algebra::SimpleMatrix>(1, 2);
  _Csurface->setZero();
  (*_Csurface)(0, 0) = 1.0;
  (*_Csurface)(0, 1) = 1.0;

  _L = std::make_shared<siconos::algebra::SimpleMatrix>(2, 1);
  (*_L)(0, 0) = -7.5125146;
  (*_L)(1, 0) = -50.04168751;

  _xHat0 = std::make_shared<siconos::algebra::SiconosVector>(2);
  (*_xHat0)(0) = (*_x0)(0);
  (*_xHat0)(1) = -5.0;

  _K = std::make_shared<siconos::algebra::SiconosVector>(3);
  (*_K)(0) = .25;
  (*_K)(1) = .125;
  (*_K)(2) = 2.0;
}

void ObserverTest::init()
{
  _DS = std::make_shared<siconos::modeling::FirstOrderLinearTIDS>(_x0, _A);
  _sensor = std::make_shared<siconos::control::LinearSensor>(_DS, _C);
  _pid = std::make_shared<siconos::control::PID>(_sensor);
  _pid->setRef(0.0);
  _pid->setK(_K);
  _pid->setB(_B);
}

void ObserverTest::tearDown() {}

void ObserverTest::test_SMO_ZOH()
{
  init();
  auto simZOH = std::make_shared<siconos::control::ControlZOHSimulation>(_t0, _T, _h);
  simZOH->addDynamicalSystem(_DS);
  simZOH->addSensor(_sensor, _h);
  auto smo = std::make_shared<siconos::control::SlidingReducedOrderObserver>(_sensor, *_xHat0,
                                                                             _C, _L);
  simZOH->addObserver(smo, _h);
  simZOH->addActuator(_pid, _h);
  simZOH->initialize();
  simZOH->run();
  auto& data = *simZOH->data();
  siconos::algebra::io::write("SMO_ZOH.dat", data, siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  auto error = siconos::algebra::io::compareRefFile(data, "SMO.ref", _tol);
  bool test = !(error  > _tol);
  std::cout << "------- Integration done -------" << test << std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_SMO_ZOH : ", test, true);
}

void ObserverTest::test_SMO_Lsodar()
{
  init();
  auto simLsodar = std::make_shared<siconos::control::ControlLsodarSimulation>(_t0, _T, _h);
  simLsodar->addDynamicalSystem(_DS);
  simLsodar->addSensor(_sensor, _h);
  auto smo = std::make_shared<siconos::control::SlidingReducedOrderObserver>(_sensor, *_xHat0,
                                                                             _C, _L);
  simLsodar->addObserver(smo, _h);
  simLsodar->addActuator(_pid, _h);
  simLsodar->initialize();
  simLsodar->run();
  auto& data = *simLsodar->data();
  siconos::algebra::io::write("SMO_Lsodar.dat", data, siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  double error = 0.0;
  bool test = !((error = siconos::algebra::io::compareRefFile(data, "SMO.ref", _tol)) > _tol);
  std::cout << "------- Integration done -------" << test << std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_SMO_Lsodar : ", test, true);
}

void ObserverTest::test_Luenberger_ZOH()
{
  init();
  auto simZOH = std::make_shared<siconos::control::ControlZOHSimulation>(_t0, _T, _h);
  simZOH->addDynamicalSystem(_DS);
  simZOH->addSensor(_sensor, _h);
  auto luenberger =
      std::make_shared<siconos::control::LuenbergerObserver>(_sensor, *_xHat0, _C, _L);
  simZOH->addObserver(luenberger, _h);
  simZOH->addActuator(_pid, _h);
  simZOH->initialize();
  simZOH->run();
  auto& data = *simZOH->data();
  siconos::algebra::io::write("Luenberger_ZOH..dat", data, siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  double error = 0.0;
  bool test =
      !((error = siconos::algebra::io::compareRefFile(data, "Luenberger.ref", _tol)) > _tol);
  std::cout << "------- Integration done -------" << test << std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_Luenberger_ZOH : ", test, true);
}

void ObserverTest::test_Luenberger_Lsodar()
{
  init();
  auto simLsodar = std::make_shared<siconos::control::ControlLsodarSimulation>(_t0, _T, _h);
  simLsodar->addDynamicalSystem(_DS);
  simLsodar->addSensor(_sensor, _h);
  auto luenberger =
      std::make_shared<siconos::control::LuenbergerObserver>(_sensor, *_xHat0, _C, _L);
  simLsodar->addObserver(luenberger, _h);
  simLsodar->addActuator(_pid, _h);
  simLsodar->initialize();
  simLsodar->run();
  auto& data = *simLsodar->data();
  siconos::algebra::io::write("Luenberger_Lsodar.dat", data, siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  double error = 0.0;
  bool test =
      !((error = siconos::algebra::io::compareRefFile(data, "Luenberger.ref", _tol)) > _tol);
  std::cout << "------- Integration done -------" << test << std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_Luenberger_Lsodar : ", test, true);
}
