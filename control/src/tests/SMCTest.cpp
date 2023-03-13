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
#include "SMCTest.hpp"

#include "ControlLsodarSimulation.hpp"
#include "ControlZOHSimulation.hpp"
#include "ExplicitLinearSMC.hpp"
#include "FirstOrderLinearTIDS.hpp"
#include "LinearSMC.hpp"
#include "LinearSensor.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "Twisting.hpp"
#include "io.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(SMCTest);

void SMCTest::setUp() {
  _A = std::make_shared<siconos::algebra::SimpleMatrix>(_n, _n, 0);
  (*_A)(0, 1) = 1.0;
  (*_A)(1, 0) = 19.0;
  (*_A)(1, 1) = -2.0;

  _x0 = std::make_shared<siconos::algebra::SiconosVector>(_n, 0);
  (*_x0)(0) = -15.0;
  (*_x0)(1) = 20.0;

  _C = std::make_shared<siconos::algebra::SimpleMatrix>(2, 2, 0);
  _C->eye();

  _B = std::make_shared<siconos::algebra::SimpleMatrix>(2, 1);
  (*_B)(1, 0) = 1.0;

  _Csurface = std::make_shared<siconos::algebra::SimpleMatrix>(1, 2, 0);
  (*_Csurface)(0, 0) = 1.0;
  (*_Csurface)(0, 1) = 1.0;
}

void SMCTest::init() {
  _DS = std::make_shared<siconos::modeling::FirstOrderLinearTIDS>(_x0, _A);
  _sensor = std::make_shared<siconos::control::LinearSensor>(_DS, _C);
  _iSMC = std::make_shared<siconos::control::LinearSMC>(_sensor, _B);
  _iSMC->setCsurface(_Csurface);

  _iSMC->display();
}

void SMCTest::init2() {
  _DS = std::make_shared<siconos::modeling::FirstOrderLinearTIDS>(_x0, _A);
  _sensor = std::make_shared<siconos::control::LinearSensor>(_DS, _C);
  _eSMC = std::make_shared<siconos::control::ExplicitLinearSMC>(_sensor, _B);
  _eSMC->setCsurface(_Csurface);
}

#ifdef HAS_EXTREME_POINT_ALGO
void SMCTest::initTwisting() {
  _DS = std::make_shared<siconos::modeling::FirstOrderLinearTIDS>(_x0, _A);
  _sensor = std::make_shared<siconos::control::LinearSensor>(_DS, _C);
  _itw = std::make_shared<siconos::control::Twisting>(_sensor, 300., _beta, _h);
  auto eye = std::make_shared<siconos::algebra::SimpleMatrix>(2, 2);
  eye->eye();
  _itw->setCsurface(eye);
}
#endif

void SMCTest::tearDown() {}
#include <EventsManager.hpp>
#include <Simulation.hpp>

void SMCTest::test_iSMC_ZOH() {
  init();
  auto simZOH = std::make_shared<siconos::control::ControlZOHSimulation>(_t0, _T, _h);
  simZOH->setSaveOnlyMainSimulation(true);
  simZOH->addDynamicalSystem(_DS);
  simZOH->addSensor(_sensor, _h);
  simZOH->addActuator(_iSMC, _h);

  simZOH->initialize();
  simZOH->run();
  auto& data = *simZOH->data();
  siconos::algebra::io::write("iSMC_ZOH.dat", data, siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);

  double error = 0.0;
  bool test =
      !((error = siconos::algebra::io::compareRefFile(data, "iSMC.ref", _tol)) >= 0.0 &&
        error > _tol);
  std::cout << "------- Integration done -------" << test << std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_Luenberger_ZOH : ", test, true);
}

void SMCTest::test_iSMC_Lsodar() {
  init();
  auto simLsodar = std::make_shared<siconos::control::ControlLsodarSimulation>(_t0, _T, _h);
  simLsodar->setSaveOnlyMainSimulation(true);
  simLsodar->addDynamicalSystem(_DS);
  simLsodar->addSensor(_sensor, _h);
  simLsodar->addActuator(_iSMC, _h);
  simLsodar->initialize();
  simLsodar->run();
  auto& data = *simLsodar->data();
  siconos::algebra::io::write("iSMC_Lsodar.dat", data, siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  double error = 0.0;
  bool test =
      !((error = siconos::algebra::io::compareRefFile(data, "iSMC.ref", _tol)) >= 0.0 &&
        error > _tol);
  std::cout << "------- Integration done -------" << test << std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_Luenberger_ZOH : ", test, true);
}

void SMCTest::test_eSMC_ZOH() {
  init2();
  auto simZOH = std::make_shared<siconos::control::ControlZOHSimulation>(_t0, _T, _h);
  simZOH->setSaveOnlyMainSimulation(true);
  simZOH->addDynamicalSystem(_DS);
  simZOH->addSensor(_sensor, _h);
  simZOH->addActuator(_eSMC, _h);
  simZOH->initialize();
  simZOH->run();
  auto& data = *simZOH->data();
  siconos::algebra::io::write("eSMC_ZOH.dat", data, siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  double error = 0.0;
  bool test =
      !((error = siconos::algebra::io::compareRefFile(data, "eSMC.ref", _tol)) >= 0.0 &&
        error > _tol);
  std::cout << "------- Integration done -------" << test << std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_Luenberger_ZOH : ", test, true);
}

void SMCTest::test_eSMC_Lsodar() {
  init2();
  auto simLsodar = std::make_shared<siconos::control::ControlLsodarSimulation>(_t0, _T, _h);
  simLsodar->setSaveOnlyMainSimulation(true);
  simLsodar->addDynamicalSystem(_DS);
  simLsodar->addSensor(_sensor, _h);
  simLsodar->addActuator(_eSMC, _h);
  simLsodar->initialize();
  simLsodar->run();

  auto& data = *simLsodar->data();
  siconos::algebra::io::write("eSMC_Lsodar.dat", data, siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  double error = 0.0;
  bool test =
      !((error = siconos::algebra::io::compareRefFile(data, "eSMC.ref", _tol)) >= 0.0 &&
        error > _tol);
  std::cout << "------- Integration done -------" << test << std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_Luenberger_ZOH : ", test, true);
}

#ifdef HAS_EXTREME_POINT_ALGO
void SMCTest::test_itw_ZOH() {
  initTwisting();
  auto simZOH = std::make_shared<siconos::control::ControlZOHSimulation>(_t0, _T, _h);
  simZOH->setSaveOnlyMainSimulation(true);
  simZOH->addDynamicalSystem(_DS);
  simZOH->addSensor(_sensor, _h);
  simZOH->addActuator(_itw, _h);
  simZOH->initialize();
  simZOH->run();
  auto& data = *simZOH->data();
  siconos::algebra::io::write("itw_ZOH.dat", data, siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  // Reference Matrix
  siconos::algebra::SimpleMatrix dataRef(data);
  dataRef.zero();
  siconos::algebra::io::read("itw.ref", dataRef);
  // it is a bad idea to compare solutions to an AVI that does not admit a unique solution
  siconos::algebra::SiconosVector lambda1 {data.size(0)};
  siconos::algebra::SiconosVector lambda2{data.size(0)};
  data.getCol(3, lambda1);
  data.getCol(4, lambda2);
  axpy(_beta, lambda2, lambda1);
  siconos::algebra::SiconosVector lambda1Ref{data.size(0)};
  siconos::algebra::SiconosVector lambda2Ref{data.size(0)};
  dataRef.getCol(3, lambda1Ref);
  dataRef.getCol(4, lambda2Ref);
  axpy(_beta, lambda2Ref, lambda1Ref);
  data.setCol(3, lambda1);
  dataRef.setCol(3, lambda1Ref);
  data.resize(data.size(0), 4);
  dataRef.resize(data.size(0), 4);
  std::cout << "------- Integration done, error = " << (data - dataRef).normInf()
	    << " -------" << std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_itw_ZOH : ", (data - dataRef).normInf() < _tol,
			       true);
}

void SMCTest::test_itw_Lsodar() {
  initTwisting();
  auto simLsodar = std::make_shared<siconos::control::ControlLsodarSimulation>(_t0, _T, _h);
  simLsodar->setSaveOnlyMainSimulation(true);
  simLsodar->addDynamicalSystem(_DS);
  simLsodar->addSensor(_sensor, _h);
  simLsodar->addActuator(_itw, _h);
  simLsodar->initialize();
  simLsodar->run();
  auto& data = *simLsodar->data();
  siconos::algebra::io::write("itw_Lsodar.dat", data, siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  // Reference Matrix
  siconos::algebra::SimpleMatrix dataRef(data);
  dataRef.zero();
  siconos::algebra::io::read("itw.ref", dataRef);
  // it is a bad idea to compare solutions to an AVI that does not admit a unique
  // solution
  siconos::algebra::SiconosVector lambda1{data.size(0)};
  siconos::algebra::SiconosVector lambda2{data.size(0)};
  data.getCol(3, lambda1);
  data.getCol(4, lambda2);
  axpy(_beta, lambda2, lambda1);
  siconos::algebra::SiconosVector lambda1Ref{data.size(0)};
  siconos::algebra::SiconosVector lambda2Ref{data.size(0)};
  dataRef.getCol(3, lambda1Ref);
  dataRef.getCol(4, lambda2Ref);
  axpy(_beta, lambda2Ref, lambda1Ref);
  data.setCol(3, lambda1);
  dataRef.setCol(3, lambda1Ref);
  data.resize(data.size(0), 4);
  dataRef.resize(data.size(0), 4);
  std::cout << "------- Integration done, error = " << (data - dataRef).normInf()
            << " -------" << std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_itw_Lsodar : ", (data - dataRef).normInf() < _tol,
                               true);
}
#endif
