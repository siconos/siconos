// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "EulerMoreauTest.hpp"

#include "EventsManager.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(EulerMoreauTest);

void EulerMoreauTest::setUp() {
  _A = std::make_shared<siconos::algebra::SiconosMatrix>(_n, _n);
  _b = std::make_shared<siconos::algebra::SiconosVector>(_n);
  _x0 = std::make_shared<siconos::algebra::SiconosVector>(_n);
  _A->setZero();
  _b->setZero();
  _x0->setZero();
}

void EulerMoreauTest::init(bool initDS) {
  if (initDS) {
    _DS = std::make_shared<siconos::modeling::FirstOrderLinearDS>(*_x0, *_A, *_b,
                                                                  siconos::algebra::alias_t);
  }

  _TD = std::make_shared<siconos::simulation::TimeDiscretisation>(_t0, _h);
  _model = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(_t0, _T);
  _sim = std::make_shared<siconos::simulation::TimeStepping>(_model, _TD, 0);
  _EulerMoreau = std::make_shared<siconos::integrators::EulerMoreauOSI>(.5);
  _model->insertDynamicalSystem(_DS);
  _sim->associate(_EulerMoreau, _DS);
  _sim->initialize();
}

void EulerMoreauTest::tearDown() {}

void EulerMoreauTest::testCstGradTIDS() {
  std::cout << "===========================================" << std::endl;
  std::cout << " ===== EulerMoreau tests start ... ===== " << std::endl;
  std::cout << "===========================================" << std::endl;
  std::cout << "------- Integrate a time-invariant coeff and linear system with constant "
               "gradients -------"
            << std::endl;
  (*_b)(0) = -1.;
  (*_x0)(0) = 5.;
  (*_x0)(1) = 10;

  _DS = std::make_shared<siconos::modeling::FirstOrderLinearDS>(*_x0, *_A, *_b,
                                                                siconos::algebra::alias_t);

  init(false);

  while (_sim->hasNextEvent()) {
    _sim->computeOneStep();
    _sim->nextStep();
  }
  siconos::algebra::SiconosVector xref{2};
  xref << 2., 10.;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testCstGradDS : ", (*_DS->x() - xref).norm() < 1e-6, true);
}

void EulerMoreauTest::testCstGradDS() {
  std::cout << "===========================================" << std::endl;
  std::cout << " ===== EulerMoreau tests start ... ===== " << std::endl;
  std::cout << "===========================================" << std::endl;
  std::cout << "------- Integrate a linear system with constant gradients -------"
            << std::endl;

  *_x0 << 5, 10;

  _DS =
      std::make_shared<siconos::modeling::FirstOrderLinearDS>(*_x0, siconos::algebra::alias_t);
  auto folds = std::dynamic_pointer_cast<siconos::modeling::FirstOrderLinearDS>(_DS);
  folds->setComputeAFunction(
      [](double time, Eigen::Ref<siconos::algebra::MapType> result) { result.setZero(); });

  folds->setComputebVectorFunction(
      [](double time, Eigen::Ref<siconos::algebra::MapVectorType> result) {
        result(0) = time;         //-1.;
        result(1) = time * time;  // 0.;
      });

  folds->setComputeMMatrixFunction(
      [](double time, Eigen::Ref<siconos::algebra::MapType> result) {
        result.setZero();
        result(0, 0) = 1;
        result(1, 1) = 1;
      });

  init(false);

  while (_sim->hasNextEvent()) {
    _sim->computeOneStep();
    _sim->nextStep();
  }

  siconos::algebra::SiconosVector xref{2};
  xref << 9.5, 19;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testCstGradDS : ", (*_DS->x() - xref).norm() < 1e-6, true);
  std::cout << std::endl << std::endl;
}

void EulerMoreauTest::testCstGradNLDS() {
  std::cout << "===========================================" << std::endl;
  std::cout << " ===== EulerMoreau tests start ... ===== " << std::endl;
  std::cout << "===========================================" << std::endl;
  std::cout << "------- Integrate a nonlinear system with constant gradients -------"
            << std::endl;
  (*_b)(0) = -1.;
  siconos::algebra::SiconosVector x0{2};
  x0 << 5, 10;
  _DS = std::make_shared<siconos::modeling::FirstOrderNonLinearDS>(x0,
                                                                   siconos::algebra::alias_t);

  auto& DSNL = static_cast<siconos::modeling::FirstOrderNonLinearDS&>(*_DS);

  DSNL.setComputefVectorFunction([](const Eigen::Ref<const siconos::algebra::SiconosVector>& x,
                                    double time,
                                    Eigen::Ref<siconos::algebra::MapVectorType> result) {
    result(0) = x(1);
    result(1) = -x(0);
  });

  DSNL.setComputeJacobianfOver_xFunction(
      [](const Eigen::Ref<const siconos::algebra::SiconosVector>& x, double time,
         Eigen::Ref<siconos::algebra::MapType> result) {
        result(0, 1) = 1.;
        result(1, 0) = -1;
      });

  DSNL.setComputeMMatrixFunction(
      [](double time, Eigen::Ref<siconos::algebra::MapType> result) {
        result.setZero();
        result(0, 0) = 1.;
        result(1, 1) = 1.;
      });
  init(false);
  while (_sim->hasNextEvent()) {
    _sim->computeOneStep();
    _sim->nextStep();
  }
  siconos::algebra::SiconosVector xref{2};
  xref(0) = x0(0) * std::cos(_T) + x0(1) * std::sin(_T);
  xref(1) = -x0(0) * std::sin(_T) + x0(1) * std::cos(_T);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testCstGradDS : ", (*_DS->x() - xref).norm() < 1e-5, true);
  std::cout << std::endl << std::endl;
}
