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
#include "testAVI.hpp"

#include "EventsManager.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "io.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(AVITest);

void AVITest::setUp() {}

void AVITest::tearDown() {}

void AVITest::testAVI() {
  std::cout << "===========================================" << std::endl;
  std::cout << " ===== AVI tests start ... ===== " << std::endl;
  std::cout << "===========================================" << std::endl;
  std::cout << "------- implicit Twisting relation  -------" << std::endl;
  _h = 1e-1;
  _T = 20.0;
  double G = 10.0;
  double beta = .3;
  _A = std::make_shared<siconos::algebra::SiconosMatrix>(_n, _n);
  _b = std::make_shared<siconos::algebra::SiconosVector>(_n);
  _x0 = std::make_shared<siconos::algebra::SiconosVector>(_n);
  _A->setZero();
  _b->setZero();
  (*_A)(0, 1) = 1.0;
  _x0->setZero();
  (*_x0)(0) = 10.0;
  (*_x0)(1) = 10.0;
  auto B = std::make_shared<siconos::algebra::SiconosMatrix>(_n, _n);
  B->setZero();
  auto C = std::make_shared<siconos::algebra::SiconosMatrix>(_n, _n);
  (*B)(1, 0) = G;
  (*B)(1, 1) = G * beta;
  C->setIdentity();
  auto rel = std::make_shared<siconos::modeling::FirstOrderLinearTIR>(*C, *B);
  // H-K representation: the feasible set is given by all the element λ such that Hλ ≥ K
  auto H = std::make_shared<siconos::algebra::SiconosMatrix>(4, 2);
  (*H)(0, 0) = 1.0;
  (*H)(1, 0) = -_h / 2.0;
  (*H)(2, 0) = -1.0;
  (*H)(3, 0) = _h / 2.0;
  (*H)(1, 1) = 1.0;
  (*H)(3, 1) = -1.0;
  auto K = std::make_shared<siconos::algebra::SiconosVector>(4);
  (*K)(0) = -1.0;
  (*K)(1) = -1.0;
  (*K)(2) = -1.0;
  (*K)(3) = -1.0;
  auto nslaw = std::make_shared<siconos::modeling::NormalConeNSL>(_n, H, K);
  _DS = std::make_shared<siconos::modeling::FirstOrderLinearDS>(*_x0, *_A, *_b);
  _TD = std::make_shared<siconos::simulation::TimeDiscretisation>(_t0, _h);
  _nsds = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(_t0, _T);
  auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, rel);
  _osi = std::make_shared<siconos::integrators::EulerMoreauOSI>(_theta);
  _nsds->insertDynamicalSystem(_DS);
  _nsds->link(inter, _DS);
  _sim = std::make_shared<siconos::simulation::TimeStepping>(_nsds, _TD);
  _sim->associate(_osi, _DS);
  auto osnspb = std::make_shared<siconos::nonsmooth_formulations::AVI>();
  _sim->insertNonSmoothProblem(osnspb);

  auto N = (unsigned)ceil((_T - _t0) / _h);
  siconos::algebra::SiconosMatrix dataPlot(N + 1, 5);
  auto& xProc = *_DS->x();
  auto& lambda = *inter->lambda(0);
  unsigned int k = 0;
  dataPlot(0, 0) = _t0;
  dataPlot(0, 1) = (*_x0)(0);
  dataPlot(0, 2) = (*_x0)(1);
  dataPlot(0, 3) = -1.0;
  dataPlot(0, 4) = -1.0;
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

  siconos::algebra::io::write("testAVI.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  // Reference Matrix
  siconos::algebra::SiconosMatrix dataPlotRef(dataPlot);
  dataPlotRef.setZero();
  siconos::algebra::io::read("testAVI.ref", dataPlotRef);
  siconos::algebra::SiconosVector err{dataPlot.cols()};
  siconos::algebra::normInfByColumn(dataPlot - dataPlotRef, err);
  siconos::algebra::print(err);

  double maxErr = err(0) > err(1) ? (err(0) > err(2) ? err(0) : err(2))
                                  : (err(1) > err(2) ? err(1) : err(2));

  std::cout << "------- Integration Ok, error = " << maxErr << " -------\n";
  if (maxErr > _tol) {
    siconos::algebra::print(dataPlot);
    siconos::algebra::SiconosMatrix diff = dataPlot - dataPlotRef;
    std::cout << "------- diff \n";
    siconos::algebra::print(diff);
  }

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAVI : ", maxErr < _tol, true);
}
