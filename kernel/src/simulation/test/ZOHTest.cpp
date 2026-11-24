// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "ZOHTest.hpp"

#include <cmath>
#include <iostream>
#include <string>

#include "EventsManager.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "StorageTools.hpp"
#include "io.hpp"

namespace {
/* Comparison of a smooth trajectory with its reference. The matrices of the ZOH scheme are
   computed by integrating an ODE with Lsodar at tolerances of a few machine epsilons: the
   round-off differences between platforms (macOS arm64 / Linux x86_64) accumulate to ~1e-12
   over the run, hence a tolerance relative to the size of the data. */
bool compareSmooth(const std::string& name, const siconos::algebra::SiconosMatrix& dataPlot,
                   const siconos::algebra::SiconosMatrix& dataPlotRef) {
  double error = siconos::algebra::normInf(dataPlot - dataPlotRef);
  double tol = 1e-10 * dataPlotRef.cwiseAbs().maxCoeff();
  std::cout << "------- " << name << ": error = " << error << " (tolerance " << tol
            << ") -------" << std::endl;
  return error < tol;
}

/* Orlov's relay controller: x0' = x1, x1' = -(2 sgn(x0) + sgn(x1)). Once the state reaches
   the switching surface x0 = 0 (t ~ 5.083), y0 = x0 is ~1e-18 and does not depend on the
   multiplier (the first row of B is zero): any lambda0 in [-1, 1] satisfies the relay law at
   this step, and the relay solver and the platform round-off select different (valid)
   values (Lemke and AVI Cao-Ferris differ by 1.7e-7, the enumerative solver more). The
   trajectory is therefore compared with the reference only until the switching surface is
   reached; afterwards, the properties of the controller are checked: the state stays in an
   O(h) neighbourhood of the sliding surface (|x0| <= 2.4e-5 with h = 1e-3) and converges to
   the origin (|x(T)| ~ 1e-8). */
bool compareRelay(const std::string& name, const siconos::algebra::SiconosMatrix& dataPlot,
                  const siconos::algebra::SiconosMatrix& dataPlotRef, double tol) {
  // first row of the reference on the switching surface
  Eigen::Index reach = 0;
  while (reach < dataPlotRef.rows() && std::abs(dataPlotRef(reach, 1)) > 1e-12) ++reach;
  double error =
      siconos::algebra::normInf(dataPlot.topRows(reach) - dataPlotRef.topRows(reach));
  double max_x0 = dataPlot.col(1).tail(dataPlot.rows() - reach).cwiseAbs().maxCoeff();
  // columns: t, x0, x1, control (2), lambda (2)
  double final_x = dataPlot.block(dataPlot.rows() - 1, 1, 1, 2).cwiseAbs().maxCoeff();
  std::cout << "------- " << name
            << ": switching surface reached at t = " << dataPlotRef(reach, 0)
            << ", error before = " << error << ", max |x0| after = " << max_x0
            << ", |x(T)| = " << final_x << " -------" << std::endl;
  return error < tol && max_x0 < 1e-4 && final_x < 1e-6;
}
}  // namespace

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
  _DS = std::make_shared<siconos::modeling::FirstOrderLinearDS>(*_x0, *_A, *_b,
                                                                siconos::algebra::alias_t);
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
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixExp1 : ", siconos::algebra::normInf(diff2Tmp) < _tol,
                               true);
  std::cout << "------- Second computation ok, error = " << siconos::algebra::normInf(diff2Tmp)
            << " -------" << std::endl;
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
  // siconos::algebra::print(dataPlot);
  std::cout << std::endl << std::endl;
  siconos::algebra::io::write("testMatrixIntegration1.dat", dataPlot,
                              siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  // Reference Matrix
  auto dataPlotRef = siconos::algebra::io::readDenseMatrix("testMatrixIntegration1.ref");
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testMatrixExp1 : ", compareSmooth("testMatrixIntegration1", dataPlot, dataPlotRef),
      true);
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
  _DS = std::make_shared<siconos::modeling::FirstOrderLinearDS>(*_x0, *_A, *_b,
                                                                siconos::algebra::alias_t);
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
  // ssiconos::algebra::print(dataPlot);
  std::cout << std::endl << std::endl;
  siconos::algebra::io::write("testMatrixIntegration2.dat", dataPlot,
                              siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  // Reference Matrix
  auto dataPlotRef = siconos::algebra::io::readDenseMatrix("testMatrixIntegration2.ref");
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testMatrixExp2 : ", compareSmooth("testMatrixIntegration2", dataPlot, dataPlotRef),
      true);
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
  _DS = std::make_shared<siconos::modeling::FirstOrderLinearDS>(*_x0, *_A, *_b,
                                                                siconos::algebra::alias_t);
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
    sampledControl.noalias() = *B * lambda;
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
  // siconos::algebra::print(dataPlot);
  std::cout << std::endl << std::endl;
  siconos::algebra::io::write("testMatrixIntegration3.dat", dataPlot,
                              siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  // Reference Matrix
  auto dataPlotRef = siconos::algebra::io::readDenseMatrix("testMatrixIntegration3.ref");
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testMatrixExp3 : ", compareRelay("testMatrixIntegration3", dataPlot, dataPlotRef, _tol),
      true);
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
  _DS = std::make_shared<siconos::modeling::FirstOrderLinearDS>(*_x0, *_A, *_b,
                                                                siconos::algebra::alias_t);
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
    sampledControl.noalias() = *B * lambda;
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
  // siconos::algebra::print(dataPlot);
  std::cout << std::endl << std::endl;
  siconos::algebra::io::write("testMatrixIntegration4.dat", dataPlot,
                              siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  // Reference Matrix
  auto dataPlotRef = siconos::algebra::io::readDenseMatrix("testMatrixIntegration4.ref");
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testMatrixExp4 : ", compareRelay("testMatrixIntegration4", dataPlot, dataPlotRef, _tol),
      true);
}
