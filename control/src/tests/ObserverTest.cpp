/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include "ControlZOHSimulation.hpp"
#include "ControlLsodarSimulation.hpp"
#include <FirstOrderLinearTIDS.hpp>
#include "LinearSensor.hpp"
#include "LinearSMC.hpp"
#include "PID.hpp"
#include "SlidingReducedOrderObserver.hpp"
#include "LuenbergerObserver.hpp"
#include "ioMatrix.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(ObserverTest);


void ObserverTest::setUp()
{
  _A.reset(new SimpleMatrix(_n, _n, 0));
  (*_A)(0, 1) = 1.0;
  (*_A)(1, 0) = -1.0;

  _x0.reset(new SiconosVector(_n, 0));
  (*_x0)(0) = 10.0;
  (*_x0)(1) = 0.0;

  _C.reset(new SimpleMatrix(1, 2, 0));
  (*_C)(0, 0) = 1.0;

  _B.reset(new SimpleMatrix(2, 1));
  (*_B)(1, 0) = 1.0;

  _Csurface.reset(new SimpleMatrix(1, 2, 0));
  (*_Csurface)(0, 0) = 1.0;
  (*_Csurface)(0, 1) = 1.0;

  _L.reset(new SimpleMatrix(2, 1));
  (*_L)(0, 0) = -7.5125146;
  (*_L)(1, 0) = -50.04168751;

  _xHat0.reset(new SiconosVector(2));
  (*_xHat0)(0) = (*_x0)(0);
  (*_xHat0)(1) = -5.0;

  _K.reset(new SiconosVector(3, 0));
  (*_K)(0) = .25;
  (*_K)(1) = .125;
  (*_K)(2) = 2.0;
}

void ObserverTest::init()
{
  _DS.reset(new FirstOrderLinearTIDS(_x0, _A));
  _sensor.reset(new LinearSensor(_DS, _C));
  _pid.reset(new PID(_sensor));
  _pid->setRef(0.0);
  _pid->setK(_K);
  _pid->setB(_B);
}

void ObserverTest::tearDown()
{}

void ObserverTest::test_SMO_ZOH()
{
  init();
  SP::ControlZOHSimulation simZOH(new ControlZOHSimulation(_t0, _T, _h));
  simZOH->addDynamicalSystem(_DS);
  simZOH->addSensor(_sensor, _h);
  simZOH->addActuator(_pid, _h);
  SP::Observer smo(new SlidingReducedOrderObserver(_sensor, *_xHat0, _C, _L));
  simZOH->addObserver(smo, _h);
  simZOH->initialize();
  simZOH->run();
  SimpleMatrix& data = *simZOH->data();
  ioMatrix::write("SMO_ZOH.dat", "ascii", data, "noDim");
  double error =0.0;
  bool test = !(ioMatrix::compareRefFile(data, "SMO.ref", _tol, error)&& error > _tol);
  std::cout << "------- Integration done -------" << test <<std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_SMO_ZOH : ", test, true);
}

void ObserverTest::test_SMO_Lsodar()
{
  init();
  SP::ControlLsodarSimulation simLsodar(new ControlLsodarSimulation(_t0, _T, _h));
  simLsodar->addDynamicalSystem(_DS);
  simLsodar->addSensor(_sensor, _h);
  simLsodar->addActuator(_pid, _h);
  SP::Observer smo(new SlidingReducedOrderObserver(_sensor, *_xHat0, _C, _L));
  simLsodar->addObserver(smo, _h);
  simLsodar->initialize();
  simLsodar->run();
  SimpleMatrix& data = *simLsodar->data();
  ioMatrix::write("SMO_Lsodar.dat", "ascii", data, "noDim");
  double error =0.0;
  bool test = !(ioMatrix::compareRefFile(data, "SMO.ref", _tol, error)&& error > _tol);
  std::cout << "------- Integration done -------" << test <<std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_SMO_Lsodar : ", test, true);
}

void ObserverTest::test_Luenberger_ZOH()
{
  init();
  SP::ControlZOHSimulation simZOH(new ControlZOHSimulation(_t0, _T, _h));
  simZOH->addDynamicalSystem(_DS);
  simZOH->addSensor(_sensor, _h);
  simZOH->addActuator(_pid, _h);
  SP::Observer luenberger(new LuenbergerObserver(_sensor, *_xHat0, _C, _L));
  simZOH->addObserver(luenberger, _h);
  simZOH->initialize();
  simZOH->run();
  SimpleMatrix& data = *simZOH->data();
  ioMatrix::write("Luenberger_ZOH.dat", "ascii", data, "noDim");
  double error =0.0;
  bool test = !(ioMatrix::compareRefFile(data, "Luenberger.ref", _tol, error)&& error > _tol);
  std::cout << "------- Integration done -------" << test <<std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_Luenberger_ZOH : ", test , true);
}

void ObserverTest::test_Luenberger_Lsodar()
{
  init();
  SP::ControlLsodarSimulation simLsodar(new ControlLsodarSimulation(_t0, _T, _h));
  simLsodar->addDynamicalSystem(_DS);
  simLsodar->addSensor(_sensor, _h);
  simLsodar->addActuator(_pid, _h);
  SP::Observer luenberger(new LuenbergerObserver(_sensor, *_xHat0, _C, _L));
  simLsodar->addObserver(luenberger, _h);
  simLsodar->initialize();
  simLsodar->run();
  SimpleMatrix& data = *simLsodar->data();
  ioMatrix::write("Luenberger_Lsodar.dat", "ascii", data, "noDim");
  double error =0.0;
  bool test = !(ioMatrix::compareRefFile(data, "Luenberger.ref", _tol, error)&& error > _tol);
  std::cout << "------- Integration done -------" << test <<std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_Luenberger_Lsodar : ", test , true);

}
