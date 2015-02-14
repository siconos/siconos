/* Siconos-Kernel, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
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
  // Reference Matrix
  SimpleMatrix dataRef(data);
  dataRef.zero();
  ioMatrix::read("SMO.ref", "ascii", dataRef);
  std::cout << "------- Integration done, error = " << (data - dataRef).normInf() << " -------" <<std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_SMO_ZOH : ", (data - dataRef).normInf() < _tol, true);
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
  // Reference Matrix
  SimpleMatrix dataRef(data);
  dataRef.zero();
  ioMatrix::read("SMO.ref", "ascii", dataRef);
  std::cout << "------- Integration done, error = " << (data - dataRef).normInf() << " -------" <<std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_SMO_Lsodar : ", (data - dataRef).normInf() < _tol, true);
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
  // Reference Matrix
  SimpleMatrix dataRef(data);
  dataRef.zero();
  ioMatrix::read("Luenberger.ref", "ascii", dataRef);
  std::cout << "------- Integration done, error = " << (data - dataRef).normInf() << " -------" <<std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_Luenberger_ZOH : ", (data - dataRef).normInf() < _tol, true);
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
  // Reference Matrix
  SimpleMatrix dataRef(data);
  dataRef.zero();
  ioMatrix::read("Luenberger.ref", "ascii", dataRef);
  std::cout << "------- Integration done, error = " << (data - dataRef).normInf() << " -------" <<std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test_Luenberger_Lsodar : ", (data - dataRef).normInf() < _tol, true);
}
