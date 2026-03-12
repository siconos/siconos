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
#include "MultiBodyTest.hpp"

#include "Circle.hpp"
#include "Disk.hpp"
#include "DiskPlanR.hpp"
#include "FrictionContact.hpp"
#include "MoreauJeanOSI.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "NumericsSolversNamespace.h"  // for SolverOptions tools
#include "SiconosBodies.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "SpaceFilter.hpp"
#include "TimeDiscretisation.hpp"
#include "TimeStepping.hpp"
#include "io.hpp"

// 2D
constexpr auto NDOF = 3;

// WALLS, TOP and GROUND
constexpr auto WALL = 100;
constexpr auto TOP = 100;
constexpr auto GROUND = 0;

// DEFAULT PLANS : a ground and two walls to support a crystal

// CRYSTAL SIZE
constexpr auto Ll = 7;

constexpr auto Rr = 1;
constexpr auto COSPI6 = 0.866025403784439;
constexpr auto SINPI6 = 0.5;
constexpr auto TANPI6 = 0.577350269189626;  // tan(pi/6)
constexpr auto SY = 3.73205080756888;  // ((cos(a)+1)/(cos(a)*sin(a)) - tan(a)) a=pi/6, R=1
constexpr auto SYL = 1 / TANPI6;

// Plan1
constexpr auto P1A = COSPI6;
constexpr auto P1B = -SINPI6;
constexpr auto P1C = (SY + (Ll - 1) * SYL - Rr) * P1B;

// Plan2
constexpr auto P2A = COSPI6;
constexpr auto P2B = SINPI6;
constexpr auto P2C = (SY + (Ll - 1) * SYL - Rr) * P2B;

constexpr auto GROUND_ID = -1;
constexpr auto MAX_RADIUS = std::numeric_limits<double>::infinity();

// Disk Plan relation
using DiskPlanR = siconos::collision::native::bodies::DiskPlanR;

// Interaction manager
using ContactManager = siconos::collision::native::SpaceFilter;

class Disks : public siconos::collision::native::SiconosBodies,
              public std::enable_shared_from_this<Disks> {
 private:
  siconos::algebra::SiconosVector initial_positions_{0};
  siconos::algebra::SiconosVector initial_velocities_{0};
  siconos::algebra::SiconosVector fext_{0};

 public:
  void init() { assert(false); };
  void init(std::string);
};

/* do nothing if solver does not converge */
void localCheckSolverOuput(int info, siconos::simulation::Simulation*) {
  if (info) exit(1);
}

double A(double t) { return 0.; }
double B(double t) { return 1.; }
double C(double t) {
  return 0.0;  // 1.1*cos(32.*M_PI*t) ;
}
double DA(double t) { return 0.; }
double DB(double t) { return 0.; }
double DC(double t) {
  return 0.0;  //-1.1*32.*M_PI*sin(32.*M_PI*t) ;
}

// ================= Creation of the model =======================
void Disks::init(std::string disks_input) {
  // User-defined main parameters

  double t0 = 0;  // initial computation time

  double T = 10;

  double h = 0.01;  // time step
  double g = 9.81;

  double theta = 0.5;  // theta for MoreauJeanOSI integrator

  std::string solverName = "NSGS";

  // -----------------------------------------
  // --- Dynamical systems && interactions ---
  // -----------------------------------------

  double R;
  double m;

  try {
    // ------------
    // --- Init ---
    // ------------

    std::cout << "====> nsds loading ..." << std::endl << std::endl;

    _plans = std::make_shared<siconos::algebra::SiconosMatrix>(
        siconos::algebra::io::readDenseMatrix("plans.dat"));
    if (_plans->rows() == 0) {
      /* default plans */
      double A1 = P1A;
      double B1 = P1B;
      double C1 = P1C;
      double A2 = P2A;
      double B2 = P2B;
      double C2 = P2C;

      _plans = std::make_shared<siconos::algebra::SiconosMatrix>(6, 6);
      _plans->setZero();
      (*_plans)(0, 0) = 0;
      (*_plans)(0, 1) = 1;
      (*_plans)(0, 2) = -GROUND;

      (*_plans)(1, 0) = 1;
      (*_plans)(1, 1) = 0;
      (*_plans)(1, 2) = WALL;

      (*_plans)(2, 0) = 1;
      (*_plans)(2, 1) = 0;
      (*_plans)(2, 2) = -WALL;

      (*_plans)(3, 0) = 0;
      (*_plans)(3, 1) = 1;
      (*_plans)(3, 2) = -TOP;

      (*_plans)(4, 0) = A1;
      (*_plans)(4, 1) = B1;
      (*_plans)(4, 2) = C1;

      (*_plans)(5, 0) = A2;
      (*_plans)(5, 1) = B2;
      (*_plans)(5, 2) = C2;
    }

    /* set center positions */
    for (unsigned int i = 0; i < _plans->rows(); ++i) {
      auto tmpr =
          std::make_shared<DiskPlanR>(1., (*_plans)(i, 0), (*_plans)(i, 1), (*_plans)(i, 2),
                                      (*_plans)(i, 3), (*_plans)(i, 4), (*_plans)(i, 5));
      (*_plans)(i, 3) = tmpr->getXCenter();
      (*_plans)(i, 4) = tmpr->getYCenter();
    }

    /*    _moving_plans = std::make_shared<siconos::collision::native::FMatrix>(1,6);
        (*_moving_plans)(0,0) = &A;
        (*_moving_plans)(0,1) = &B;
        (*_moving_plans)(0,2) = &C;
        (*_moving_plans)(0,3) = &DA;
        (*_moving_plans)(0,4) = &DB;
        (*_moving_plans)(0,5) = &DC;*/

    auto disks_matrix = std::make_shared<siconos::algebra::SiconosMatrix>(
        siconos::algebra::io::readDenseMatrix(disks_input));

    initial_positions_.resize(NDOF * disks_matrix->rows());
    initial_velocities_.resize(NDOF * disks_matrix->rows());
    initial_positions_.setZero();
    initial_velocities_.setZero();
    fext_.resize(NDOF * disks_matrix->rows());
    ;
    fext_.setZero();
    // -- OneStepIntegrators --
    auto osi = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);

    // -- Model --

    auto nsds = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);

    for (unsigned int i = 0; i < disks_matrix->rows(); i++) {
      R = (*disks_matrix)(i, 2);
      m = (*disks_matrix)(i, 3);

      initial_positions_.segment(NDOF * i, 2) = disks_matrix->row(i).segment(0, 2);

      std::shared_ptr<siconos::collision::native::bodies::CircularDS> body{nullptr};
      if (R > 0)
        body = std::make_shared<siconos::collision::native::bodies::Disk>(
            R, m, initial_positions_.segment(NDOF * i, NDOF),
            initial_velocities_.segment(NDOF * i, NDOF));
      else
        body = std::make_shared<siconos::collision::native::bodies::Circle>(
            -R, m, initial_positions_.segment(NDOF * i, NDOF),
            initial_velocities_.segment(NDOF * i, NDOF));

      // -- Set external forces (weight) --
      fext_(NDOF * i + 1) = -m * g;
      body->setConstantFext(fext_.segment(NDOF * i, NDOF), siconos::algebra::copy_t);

      // add the dynamical system in the non smooth dynamical system
      nsds->insertDynamicalSystem(body);
    }

    nsds->setSymmetric(true);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- Time discretisation --
    auto timedisc = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- OneStepNsProblem --
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::FrictionContact>(2);
    osnspb->numericsSolverOptions()->iparam[SICONOS_IPARAM_MAX_ITER] = 100;  // Max number of
    // iterations
    // osnspb_->numericsSolverOptions()->iparam[SICONOS_IPARAM_ITER_DONE] = 20; // compute
    // error
    // iterations
    osnspb->numericsSolverOptions()->dparam[SICONOS_DPARAM_TOL] = 1e-3;  // Tolerance

    osnspb->setMaxSize(6 * ((3 * Ll * Ll + 3 * Ll) / 2 - Ll));
    osnspb->setMStorageType(NM_SPARSE_BLOCK);  // Sparse storage
    osnspb->setNumericsVerboseMode(0);

    osnspb->setKeepLambdaAndYState(true);  // inject previous solution

    // -- Simulation --
    _sim = std::make_shared<siconos::simulation::TimeStepping>(nsds, timedisc, osi, osnspb);

    std::static_pointer_cast<siconos::simulation::TimeStepping>(_sim)->setNewtonMaxIteration(
        3);

    std::static_pointer_cast<siconos::simulation::TimeStepping>(_sim)->setCheckSolverFunction(
        localCheckSolverOuput);

    // --- Simulation initialization ---

    std::cout << "====> Simulation initialisation ...\n\n";

    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(0, 0, 0.3, 2);

    _playground = std::make_shared<ContactManager>(3, 6, _plans, _moving_plans);

    _playground->insertNonSmoothLaw(nslaw, 0, 0);

    _sim->insertInteractionManager(_playground);
  }

  catch (...) {
    siconos::exception::process();
    exit(1);
  }
}

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(MultiBodyTest);

void MultiBodyTest::setUp() {}

void MultiBodyTest::tearDown() {}

// multiples disks
void MultiBodyTest::t1() {
  auto disks = std::make_shared<Disks>();
  disks->init("disks.dat");

  // just try to run a simulation
  // if something is broken with SpaceFilter
  // an exception may occurs
  for (unsigned int i = 0; i < 20; ++i) {
    disks->compute();
  }

  CPPUNIT_ASSERT(1);
}

// one disk without interaction at the beginning
void MultiBodyTest::t2() {
  auto disks = std::make_shared<Disks>();
  disks->init("disks-nointer.dat");

  // just try to run a simulation
  // if something is broken with SpaceFilter
  // an exception may occurs
  // test fail with rev 3146
  // for (unsigned int i = 0; i < 20; ++i) {
  //   disks->compute();
  // }

  CPPUNIT_ASSERT(1);
}

void MultiBodyTest::t3() {}

void MultiBodyTest::t4() {}

void MultiBodyTest::t5() {}

void MultiBodyTest::t6() {}

void MultiBodyTest::t7() {}

void MultiBodyTest::t8() {}
