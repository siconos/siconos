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
#include "ZOHTest.hpp"

#include "EventsManager.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixVectorOp.hpp"
#include "SiconosVector.hpp"
#include "io.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(ZOHTest);

void ZOHTest::setUp() {
  _A = std::make_shared<siconos::algebra::SiconosMatrix>(_n, _n);
  _b = std::make_shared<siconos::algebra::SiconosVector>(_n);
  _x0 = std::make_shared<siconos::algebra::SiconosVector>(_n);
  _A->setZero();
  _b->setZero();
  _x0->setZero();
}

void ZOHTest::init() {
  _DS = std::make_shared<siconos::modeling::FirstOrderLinearDS>(*_x0, *_A, *_b);
  _TD = std::make_shared<siconos::simulation::TimeDiscretisation>(_t0, _h);
  _model = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(_t0, _T);
  _sim = std::make_shared<siconos::simulation::TimeStepping>(_model, _TD, 0);
  _ZOH = std::make_shared<siconos::integrators::ZeroOrderHoldOSI>();
  _model->insertDynamicalSystem(_DS);
  _sim->associate(_ZOH, _DS);
  _sim->initialize();
}

void ZOHTest::tearDown() {}

void ZOHTest::testMatrixExp0() {
  std::cout << "===========================================" << std::endl;
  std::cout << " ===== ZOH tests start ... ===== " << std::endl;
  std::cout << "===========================================" << std::endl;
  std::cout << "------- Compute matrix exponential of the identity matrix -------"
            << std::endl;
  _A->setIdentity();
  init();
  _sim->computeOneStep();
  _sim->nextStep();
  auto tmpM = std::make_shared<siconos::algebra::SiconosMatrix>(_n, _n);
  tmpM->setZero();
  tmpM->setIdentity();
  *tmpM = (*tmpM) * exp(_h);
  const auto& Phi = _ZOH->Ad(_DS);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixExp0 : ", Phi.rows() == _n, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixExp0 : ", Phi.cols() == _n, true);
  siconos::algebra::SiconosMatrix diffTmp = *tmpM - Phi;
  double diff = siconos::algebra::normInf(diffTmp);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixExp0 : ", diff < _tol, true);
  std::cout << "------- First computation ok, error = " << diff << " -------" << std::endl;
  std::cout << std::endl << std::endl;
}

void ZOHTest::testMatrixExp1() {
  std::cout << "===========================================" << std::endl;
  std::cout << " ===== ZOH tests start ... ===== " << std::endl;
  std::cout << "===========================================" << std::endl;
  std::cout << "------- Compute matrix exponential of a upper triangular matrix -------"
            << std::endl;
  _A->setZero();
  (*_A)(0, 1) = 1;
  init();
  _sim->computeOneStep();
  _sim->nextStep();
  auto tmpM = std::make_shared<siconos::algebra::SiconosMatrix>(_n, _n);
  tmpM->setZero();
  tmpM->setIdentity();
  (*tmpM)(0, 1) = _h;
  const auto& Phi = _ZOH->Ad(_DS);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixExp1 : ", Phi.rows() == _n, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixExp1 : ", Phi.cols() == _n, true);
  siconos::algebra::SiconosMatrix diff2Tmp = *tmpM - Phi;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixExp1 : ", siconos::algebra::normInf(diff2Tmp) < _tol, true);
  std::cout << "------- Second computation ok, error = " << siconos::algebra::normInf(diff2Tmp) << " -------"
            << std::endl;
  std::cout << std::endl << std::endl;
}

void ZOHTest::testMatrixIntegration1() {
  std::cout << "===========================================" << std::endl;
  std::cout << " ===== ZOH tests start ... ===== " << std::endl;
  std::cout << "===========================================" << std::endl;
  std::cout << "------- Integrate an oscillator -------" << std::endl;
  _A->setZero();
  (*_A)(0, 1) = 1;
  (*_A)(1, 0) = -1;
  _x0->setZero();
  (*_x0)(0) = 1;
  init();
  auto N = (unsigned)ceil((_T - _t0) / _h);
  siconos::algebra::SiconosMatrix dataPlot(N + 1, 3);
  auto& xProc = *_DS->x();
  unsigned int k = 0;
  dataPlot(0, 0) = _t0;
  dataPlot(0, 1) = (*_x0)(0);
  dataPlot(0, 2) = (*_x0)(1);
  while (_sim->hasNextEvent()) {
    _sim->computeOneStep();
    k++;
    dataPlot(k, 0) = _sim->nextTime();
    dataPlot(k, 1) = xProc(0);
    dataPlot(k, 2) = xProc(1);
    _sim->nextStep();
    //_sim->eventsManager()->display();
  }
  siconos::algebra::print(dataPlot);
  std::cout << std::endl << std::endl;
  siconos::algebra::io::write("testMatrixIntegration1.dat", dataPlot,
                              siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  // Reference Matrix
  siconos::algebra::SiconosMatrix dataPlotRef(dataPlot);
  dataPlotRef.setZero();
  // magic line to compute the following:
  //  python -c "import numpy as np; t = np.linspace(0, 9.9, 100);
  //  np.savetxt('testMatrixIntegration1.ref', np.transpose([t, np.cos(t), -np.sin(t)]))" &&
  //  sed -i "1i100 3" testMatrixIntegration1.ref
  siconos::algebra::io::read("testMatrixIntegration1.ref", dataPlotRef);
  siconos::algebra::SiconosMatrix diff = dataPlot - dataPlotRef;
  std::cout << "------- Integration Ok, error = " << siconos::algebra::normInf(diff) << " -------" << std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixExp1 : ", siconos::algebra::normInf(diff) < _tol, true);
}

void ZOHTest::testMatrixIntegration2() {
  std::cout << "===========================================" << std::endl;
  std::cout << " ===== ZOH tests start ... ===== " << std::endl;
  std::cout << "===========================================" << std::endl;
  std::cout << "------- Integrate x \\in -sgn(x)  -------" << std::endl;
  _A->setZero();
  _x0->setZero();
  (*_x0)(0) = 1;
  (*_x0)(1) = -1;
  auto B = std::make_shared<siconos::algebra::SiconosMatrix>(_n, _n);
  auto C = std::make_shared<siconos::algebra::SiconosMatrix>(_n, _n);
  B->setIdentity();
  C->setIdentity();
  auto rel = std::make_shared<siconos::modeling::FirstOrderLinearTIR>(*C, *B);
  auto D = std::make_shared<siconos::algebra::SiconosMatrix>(_n, _n);
  D->setZero();
  rel->setConstantD(*D);
  auto nslaw = std::make_shared<siconos::modeling::RelayNSL>(_n);
  _DS = std::make_shared<siconos::modeling::FirstOrderLinearDS>(*_x0, *_A, *_b);
  _TD = std::make_shared<siconos::simulation::TimeDiscretisation>(_t0, _h);
  _model = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(_t0, _T);
  auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, rel);
  _ZOH = std::make_shared<siconos::integrators::ZeroOrderHoldOSI>();
  _model->insertDynamicalSystem(_DS);
  _model->link(inter, _DS);
  _model->setControlProperty(inter, true);
  _sim = std::make_shared<siconos::simulation::TimeStepping>(_model, _TD, 1);
  _sim->associate(_ZOH, _DS);
  auto osnspb = std::make_shared<siconos::nonsmooth_formulations::Relay>();
  _sim->insertNonSmoothProblem(osnspb);
  _sim->initialize();
  auto N = (unsigned)ceil((_T - _t0) / _h);
  siconos::algebra::SiconosMatrix dataPlot(N + 1, 5);
  auto& xProc = *_DS->x();
  auto& lambda = *inter->lambda(0);
  unsigned int k = 0;
  dataPlot(0, 0) = _t0;
  dataPlot(0, 1) = (*_x0)(0);
  dataPlot(0, 2) = (*_x0)(1);
  dataPlot(0, 3) = 0;
  dataPlot(0, 4) = 0;

  while (_sim->hasNextEvent()) {
    _sim->computeOneStep();
    k++;
    dataPlot(k, 0) = _sim->nextTime();
    dataPlot(k, 1) = xProc(0);
    dataPlot(k, 2) = xProc(1);
    dataPlot(k, 3) = lambda(0);
    dataPlot(k, 4) = lambda(1);
    _sim->nextStep();
  }
  siconos::algebra::print(dataPlot);
  std::cout << std::endl << std::endl;
  siconos::algebra::io::write("testMatrixIntegration2.dat", dataPlot,
                              siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  // Reference Matrix
  siconos::algebra::SiconosMatrix dataPlotRef(dataPlot);
  siconos::algebra::io::read("testMatrixIntegration2.ref", dataPlotRef);
  siconos::algebra::SiconosMatrix diff = dataPlot - dataPlotRef;
  std::cout << "------- Integration Ok, error = " << siconos::algebra::normInf(diff) << " -------" << std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixExp2 : ", siconos::algebra::normInf(diff) < _tol, true);
}

void ZOHTest::testMatrixIntegration3() {
  std::cout << "===========================================" << std::endl;
  std::cout << " ===== ZOH tests start ... ===== " << std::endl;
  std::cout << "===========================================" << std::endl;
  std::cout << "------- Integrate Orlov's controller  -------" << std::endl;
  _h = .001;
  _A->setZero();
  (*_A)(0, 1) = 1;
  _x0->setZero();
  (*_x0)(0) = 1;
  (*_x0)(1) = 1;
  auto B = std::make_shared<siconos::algebra::SiconosMatrix>(_n, _n);
  B->setZero();
  auto C = std::make_shared<siconos::algebra::SiconosMatrix>(_n, _n);
  (*B)(1, 0) = 2;
  (*B)(1, 1) = 1;
  C->setIdentity();
  auto rel = std::make_shared<siconos::modeling::FirstOrderLinearTIR>(*C, *B);
  auto D = std::make_shared<siconos::algebra::SiconosMatrix>(_n, _n);
  D->setZero();
  rel->setConstantD(*D);
  auto nslaw = std::make_shared<siconos::modeling::RelayNSL>(_n);
  _DS = std::make_shared<siconos::modeling::FirstOrderLinearDS>(*_x0, *_A, *_b);
  _TD = std::make_shared<siconos::simulation::TimeDiscretisation>(_t0, _h);
  _model = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(_t0, _T);
  auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, rel);
  _ZOH = std::make_shared<siconos::integrators::ZeroOrderHoldOSI>();
  _model->insertDynamicalSystem(_DS);
  _model->link(inter, _DS);
  _model->setControlProperty(inter, true);
  _sim = std::make_shared<siconos::simulation::TimeStepping>(_model, _TD, 1);
  _sim->associate(_ZOH, _DS);
  auto osnspb = std::make_shared<siconos::nonsmooth_formulations::Relay>();
  _sim->insertNonSmoothProblem(osnspb);
  _sim->initialize();
  auto N = (unsigned)ceil((_T - _t0) / _h);
  siconos::algebra::SiconosMatrix dataPlot(N + 1, 7);
  auto& xProc = *_DS->x();
  auto& lambda = *inter->lambda(0);
  siconos::algebra::SiconosVector sampledControl(_n);
  unsigned int k = 0;
  dataPlot(0, 0) = _t0;
  dataPlot(0, 1) = (*_x0)(0);
  dataPlot(0, 2) = (*_x0)(1);
  dataPlot(0, 3) = 0;
  dataPlot(0, 4) = 0;
  dataPlot(0, 5) = 0;
  dataPlot(0, 6) = 0;
  while (_sim->hasNextEvent()) {
    _sim->computeOneStep();
    siconos::algebra::prod(*B, lambda, sampledControl, true);
    k++;
    dataPlot(k, 0) = _sim->nextTime();
    dataPlot(k, 1) = xProc(0);
    dataPlot(k, 2) = xProc(1);
    dataPlot(k, 3) = sampledControl(0);
    dataPlot(k, 4) = sampledControl(1);
    dataPlot(k, 5) = lambda(0);
    dataPlot(k, 6) = lambda(1);
    _sim->nextStep();
  }
  siconos::algebra::print(dataPlot);
  std::cout << std::endl << std::endl;
  siconos::algebra::io::write("testMatrixIntegration3.dat", dataPlot,
                              siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  // Reference Matrix
  siconos::algebra::SiconosMatrix dataPlotRef(dataPlot);
  siconos::algebra::io::read("testMatrixIntegration3.ref", dataPlotRef);
  siconos::algebra::SiconosMatrix diff = dataPlot - dataPlotRef;
  std::cout << "------- Integration Ok, error = " << siconos::algebra::normInf(diff) << " -------" << std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixExp3 : ", siconos::algebra::normInf(diff) < _tol, true);
}

void ZOHTest::testMatrixIntegration4() {
  std::cout << "===========================================" << std::endl;
  std::cout << " ===== ZOH tests start ... ===== " << std::endl;
  std::cout << "===========================================" << std::endl;
  std::cout << "------- Integrate Orlov's controller  -------" << std::endl;
  _h = .001;
  _A->setZero();
  (*_A)(0, 1) = 1;
  _x0->setZero();
  (*_x0)(0) = 1;
  (*_x0)(1) = 1;
  auto B = std::make_shared<siconos::algebra::SiconosMatrix>(_n, _n);
  B->setZero();
  auto C = std::make_shared<siconos::algebra::SiconosMatrix>(_n, _n);
  (*B)(1, 0) = 2;
  (*B)(1, 1) = 1;
  C->setIdentity();
  auto rel = std::make_shared<siconos::modeling::FirstOrderLinearTIR>(*C, *B);
  auto D = std::make_shared<siconos::algebra::SiconosMatrix>(_n, _n);
  D->setZero();
  rel->setConstantD(*D);
  auto nslaw = std::make_shared<siconos::modeling::RelayNSL>(_n);
  _DS = std::make_shared<siconos::modeling::FirstOrderLinearDS>(*_x0, *_A, *_b);
  _TD = std::make_shared<siconos::simulation::TimeDiscretisation>(_t0, _h);
  _model = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(_t0, _T);
  auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, rel);
  _ZOH = std::make_shared<siconos::integrators::ZeroOrderHoldOSI>();
  _model->insertDynamicalSystem(_DS);
  _model->link(inter, _DS);
  _model->setControlProperty(inter, true);
  _sim = std::make_shared<siconos::simulation::TimeStepping>(_model, _TD, 1);
  _sim->associate(_ZOH, _DS);
  auto osnspb = std::make_shared<siconos::nonsmooth_formulations::Relay>();
  _sim->insertNonSmoothProblem(osnspb);
  _sim->initialize();
  auto N = (unsigned)ceil((_T - _t0) / _h);
  siconos::algebra::SiconosMatrix dataPlot(N + 1, 7);
  auto& xProc = *_DS->x();
  auto& lambda = *inter->lambda(0);
  siconos::algebra::SiconosVector sampledControl(_n);
  unsigned int k = 0;
  dataPlot(0, 0) = _t0;
  dataPlot(0, 1) = (*_x0)(0);
  dataPlot(0, 2) = (*_x0)(1);
  dataPlot(0, 3) = 0;
  dataPlot(0, 4) = 0;
  dataPlot(0, 5) = 0;
  dataPlot(0, 6) = 0;
  while (_sim->hasNextEvent()) {
    _sim->computeOneStep();
    siconos::algebra::prod(*B, lambda, sampledControl, true);
    k++;
    dataPlot(k, 0) = _sim->nextTime();
    dataPlot(k, 1) = xProc(0);
    dataPlot(k, 2) = xProc(1);
    dataPlot(k, 3) = sampledControl(0);
    dataPlot(k, 4) = sampledControl(1);
    dataPlot(k, 5) = lambda(0);
    dataPlot(k, 6) = lambda(1);
    _sim->nextStep();
  }
  siconos::algebra::print(dataPlot);
  std::cout << std::endl << std::endl;
  siconos::algebra::io::write("testMatrixIntegration4.dat", dataPlot,
                              siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  // Reference Matrix
  siconos::algebra::SiconosMatrix dataPlotRef(dataPlot);
  dataPlotRef.setZero();
  siconos::algebra::io::read("testMatrixIntegration4.ref", dataPlotRef);
  siconos::algebra::SiconosMatrix diff = dataPlot - dataPlotRef;
  std::cout << "------- Integration Ok, error = " << siconos::algebra::normInf(diff) << " -------" << std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixExp4 : ", siconos::algebra::normInf(diff) < _tol, true);
}
