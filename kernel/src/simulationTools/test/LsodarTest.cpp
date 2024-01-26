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
#include "LsodarTest.hpp"
#include "EventsManager.hpp"
#include "SimpleMatrix.hpp"
#include "FirstOrderLinearTIDS.hpp"
#include "LsodarOSI.hpp"
#include "SiconosVector.hpp"


#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LsodarTest);

static void computef1(double, unsigned int, double*, double* f, unsigned int, double*)
{
  f[0] = -1.;
  f[1] = 0.;
}

static void computeA1(double, unsigned int, double*, double* A, unsigned int, double*)
{
  A[0] = 0.;
  A[1] = 0.;
  A[2] = 0.;
  A[3] = 0.;
}

void LsodarTest::setUp()
{
  _A = std::make_shared<siconos::algebra::SimpleMatrix>(_n, _n,0);// siconos::algebra::UBlasType::BLOCK);
  _b = std::make_shared<siconos::algebra::SiconosVector>(_n, 0);
  _x0 = std::make_shared<siconos::algebra::SiconosVector>(_n, 0);
}

void LsodarTest::init(bool initDS)
{
  if(initDS)
  {
    _DS = std::make_shared<siconos::modeling::FirstOrderLinearTIDS>(_x0, _A, _b);
  }

  _TD = std::make_shared<siconos::simulation::TimeDiscretisation>(_t0, _h);
  _model = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(_t0, _T);
  _sim = std::make_shared<siconos::simulation::EventDriven>(_model, _TD, 0);
  _Lsodar = std::make_shared<siconos::integrators::LsodarOSI>();
  _model->insertDynamicalSystem(_DS);
  _sim->associate(_Lsodar, _DS);
  _sim->initialize();
}

void LsodarTest::tearDown()
{}

void LsodarTest::testCstGradTIDS()
{
  std::cout << "===========================================" <<std::endl;
  std::cout << " ===== Lsodar tests start ... ===== " <<std::endl;
  std::cout << "===========================================" <<std::endl;
  std::cout << "------- Integrate a TIL system with constant gradients -------" <<std::endl;
  _b->setValue(0, -1.);
  _x0->setValue(0, 5.);
  _x0->setValue(1, 10);

  _DS = std::make_shared<siconos::modeling::FirstOrderLinearTIDS>(_x0, _A, _b);

  init(false);

  while(_sim->hasNextEvent())
  {
    _sim->advanceToEvent();
    _sim->processEvents();
  }

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testCstGradTIDS : ", fabs(_DS->x()->getValue(0) +  5.) < _tol, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testCstGradTIDS : ", fabs(_DS->x()->getValue(1) - 10.) < _tol, true);
  std::cout <<std::endl <<std::endl;
}

void LsodarTest::testCstGradDS()
{
  std::cout << "===========================================" <<std::endl;
  std::cout << " ===== Lsodar tests start ... ===== " <<std::endl;
  std::cout << "===========================================" <<std::endl;
  std::cout << "------- Integrate a L system with constant gradients -------" <<std::endl;
  _b->setValue(0, -1.);
  _x0->setValue(0, 5.);
  _x0->setValue(1, 10);

  _DS = std::make_shared<siconos::modeling::FirstOrderLinearDS>(_x0, _A, _b);

  init(false);

  while(_sim->hasNextEvent())
  {
    _sim->advanceToEvent();
    _sim->processEvents();
  }

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testCstGradDS : ", fabs(_DS->x()->getValue(0) +  5.) < _tol, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testCstGradDS : ", fabs(_DS->x()->getValue(1) - 10.) < _tol, true);
  std::cout <<std::endl <<std::endl;
}

void LsodarTest::testCstGradNLDS()
{
  std::cout << "===========================================" <<std::endl;
  std::cout << " ===== Lsodar tests start ... ===== " <<std::endl;
  std::cout << "===========================================" <<std::endl;
  std::cout << "------- Integrate a NL system with constant gradients -------" <<std::endl;
  _b->setValue(0, -1.);
  _x0->setValue(0, 5.);
  _x0->setValue(1, 10);

  _DS = std::make_shared<siconos::modeling::FirstOrderNonLinearDS>(_x0);
  auto& DSNL = static_cast<siconos::modeling::FirstOrderNonLinearDS&>(*_DS);
  DSNL.setComputeFFunction(&computef1);
  DSNL.setComputeJacobianfxFunction(&computeA1);

  init(false);

  while(_sim->hasNextEvent())
  {
    _sim->advanceToEvent();
    _sim->processEvents();
  }

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testCstGradNLDS : ", fabs(_DS->x()->getValue(0) +  5.) < _tol, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testCstGradNLDS : ", fabs(_DS->x()->getValue(1) - 10.) < _tol, true);
  std::cout <<std::endl <<std::endl;
}

