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

#include "NonSmoothDynamicalSystem.hpp"
#include "SiconosConfig.h"
#define WITH_SERIALIZATION

#ifdef HAVE_SICONOS_MECHANICS
#include "MechanicsIO.hpp"
#endif

#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/matrix_sparse.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/vector_sparse.hpp>

#include "KernelTest.hpp"
#include "SiconosKernel.hpp"

#define DEBUG_MESSAGES 1
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

#include "SiconosFull.hpp"

CPPUNIT_TEST_SUITE_REGISTRATION(KernelTest);

void KernelTest::setUp() { BBxml = "BouncingBall1.xml"; }

void KernelTest::tearDown(){};

void KernelTest::t0() {
  auto q = std::make_shared<siconos::algebra::SiconosVector>(3);
  auto q0 = std::make_shared<siconos::algebra::SiconosVector>(3);

  (*q)(0) = 1.0;
  (*q)(1) = 2.0;
  (*q)(2) = 2.0;

  std::ofstream ofs("Kernelt0.xml");
  {
    boost::archive::xml_oarchive oa(ofs);
    oa.register_type(static_cast<siconos::algebra::SiconosVector*>(nullptr));
    oa << NVP(q);
  }

  std::ifstream ifs("Kernelt0.xml");
  {
    boost::archive::xml_iarchive ia(ifs);
    ia.register_type(static_cast<siconos::algebra::SiconosVector*>(nullptr));
    ia >> NVP(q0);
  }

  CPPUNIT_ASSERT(*q0 == *q);
}

void KernelTest::t1() {
  auto m1 = std::make_shared<siconos::algebra::SimpleMatrix>(3, 3);
  auto m2 = std::make_shared<siconos::algebra::SimpleMatrix>(3, 3);

  m1->eye();
  (*m1)(1, 0) = 3.0;
  (*m1)(2, 1) = -7;

  std::ofstream ofs("Kernelt1.xml");
  {
    boost::archive::xml_oarchive oa(ofs);
    oa.register_type(static_cast<siconos::algebra::SimpleMatrix*>(nullptr));
    oa << NVP(m1);
  }

  std::ifstream ifs("Kernelt1.xml");
  {
    boost::archive::xml_iarchive ia(ifs);
    ia.register_type(static_cast<siconos::algebra::SimpleMatrix*>(nullptr));
    ia >> NVP(m2);
  }

  m1->display();
  m2->display();

  CPPUNIT_ASSERT(*m1 == *m2);
}

void KernelTest::t2() {
  auto m = std::make_shared<siconos::algebra::SimpleMatrix>(3, 3);
  auto v = std::make_shared<siconos::algebra::SiconosVector>(3);
  auto q = std::make_shared<siconos::algebra::SiconosVector>(3);

  m->eye();

  auto ds1 = std::make_shared<siconos::modeling::LagrangianDS>(q, v, m);
  auto ds2 = std::make_shared<siconos::modeling::LagrangianDS>(q, v, m);

  std::ofstream ofs("Kernelt2.xml");
  {
    boost::archive::xml_oarchive oa(ofs);
    oa.register_type(static_cast<siconos::algebra::SimpleMatrix*>(nullptr));
    oa.register_type(static_cast<siconos::algebra::SiconosVector*>(nullptr));
    oa.register_type(static_cast<siconos::modeling::LagrangianDS*>(nullptr));
    oa << NVP(ds1);
  }

  std::ifstream ifs("Kernelt2.xml");
  {
    boost::archive::xml_iarchive ia(ifs);
    ia.register_type(static_cast<siconos::algebra::SimpleMatrix*>(nullptr));
    ia.register_type(static_cast<siconos::algebra::SiconosVector*>(nullptr));
    ia.register_type(static_cast<siconos::modeling::LagrangianDS*>(nullptr));
    ia >> NVP(ds2);
  }

  CPPUNIT_ASSERT(*(std::static_pointer_cast<siconos::modeling::LagrangianDS>(ds1)->mass()) ==
                 *(std::static_pointer_cast<siconos::modeling::LagrangianDS>(ds2)->mass()));
  CPPUNIT_ASSERT(*(std::static_pointer_cast<siconos::modeling::LagrangianDS>(ds1)->q()) ==
                 *(std::static_pointer_cast<siconos::modeling::LagrangianDS>(ds2)->q()));
  CPPUNIT_ASSERT(
      *(std::static_pointer_cast<siconos::modeling::LagrangianDS>(ds1)->velocity()) ==
      *(std::static_pointer_cast<siconos::modeling::LagrangianDS>(ds2)->velocity()));
}

void KernelTest::t3() {
  std::shared_ptr<SolverOptions> so{solver_options_create(SICONOS_FRICTION_3D_NSN_AC),
                                    solver_options_delete};
  auto sor = std::make_shared<SolverOptions>();
  so->numberOfInternalSolvers = 1;
  so->internalSolvers = (SolverOptions**)calloc(1, sizeof(SolverOptions*));
  so->internalSolvers[0] = solver_options_create(SICONOS_FRICTION_3D_NSN_AC);

  std::ofstream ofs("SolverOptions.xml");
  {
    boost::archive::xml_oarchive oa(ofs);
    oa << NVP(so);
  }

  std::ifstream ifs("SolverOptions.xml");
  {
    boost::archive::xml_iarchive ia(ifs);
    ia >> NVP(sor);
  }

  CPPUNIT_ASSERT((so->iSize == sor->iSize));
}

// void KernelTest::t4()
// {
//   auto m = std::make_shared<siconos::algebra::SimpleMatrix>(3, 3));
//   auto v = std::make_shared<siconos::algebra::SiconosVector>(3));
//   auto q = std::make_shared<siconos::algebra::SiconosVector>(3));

//   m->eye();

//   auto ds1= std::make_shared<siconos::modeling::LagrangianDS>(q, v, m));
//   auto ds2= std::make_shared<siconos::modeling::LagrangianDS>(q, v, m));

//   SP::DynamicalSystemsSet dsset(new DynamicalSystemsSet());

//   dsset->insert(ds1);
//   dsset->insert(ds2);

//   std::ofstream ofs("t4.xml");
//   {
//     boost::archive::xml_oarchive oa(ofs);
//     oa.register_type(static_cast<SimpleMatrix*>(nullptr));
//     oa.register_type(static_cast<siconos::algebra::SiconosVector*>(nullptr));
//     oa.register_type(static_cast<siconos::modeling::LagrangianDS*>(nullptr));
//     oa << NVP(dsset);
//   }

//   SP::DynamicalSystemsSet dssetfromfile(new DynamicalSystemsSet());

//   std::ifstream ifs("t4.xml");
//   {
//     boost::archive::xml_iarchive ia(ifs);
//     ia.register_type(static_cast<siconos::algebra::SimpleMatrix*>(nullptr));
//     ia.register_type(static_cast<siconos::algebra::SiconosVector*>(nullptr));
//     ia.register_type(static_cast<siconos::modeling::LagrangianDS*>(nullptr));
//     ia >> NVP(dssetfromfile);
//   }

// }

#include "SiconosRestart.hpp"

using namespace std;
void KernelTest::t5() {
  // ================= Creation of the model =======================

  // User-defined main parameters
  unsigned int nDof = 3;       // degrees of freedom for the ball
  double t0 = 0;               // initial computation time
  double T = 10;               // final computation time
  double h = 0.005;            // time step
  double position_init = 1.0;  // initial position for lowest bead.
  double velocity_init = 0.0;  // initial velocity for lowest bead.
  double theta = 0.5;          // theta for MoreauJeanOSI integrator
  double R = 0.1;              // Ball radius
  double m = 1;                // Ball mass
  double g = 9.81;             // Gravity
  // -------------------------
  // --- Dynamical systems ---
  // -------------------------

  cout << "====> Model loading ..." << endl << endl;

  auto Mass = std::make_shared<siconos::algebra::SimpleMatrix>(nDof, nDof);
  (*Mass)(0, 0) = m;
  (*Mass)(1, 1) = m;
  (*Mass)(2, 2) = 3. / 5 * m * R * R;

  // -- Initial positions and velocities --
  auto q0 = std::make_shared<siconos::algebra::SiconosVector>(nDof);
  auto v0 = std::make_shared<siconos::algebra::SiconosVector>(nDof);
  (*q0)(0) = position_init;
  (*v0)(0) = velocity_init;

  // -- The dynamical system --
  auto ball = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(q0, v0, Mass);

  // -- Set external forces (weight) --
  auto weight = std::make_shared<siconos::algebra::SiconosVector>(nDof);
  (*weight)(0) = -m * g;
  ball->setFExtPtr(weight);

  // --------------------
  // --- Interactions ---
  // --------------------

  // -- nslaw --
  double e = 0.9;

  // Interaction ball-floor
  //
  auto H = std::make_shared<siconos::algebra::SimpleMatrix>(1, nDof);
  (*H)(0, 0) = 1.0;

  auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
  auto relation = std::make_shared<siconos::modeling::LagrangianLinearTIR>(H);

  auto inter = std::make_shared < Interaction(nslaw, relation);

  // -------------
  // --- Model ---
  // -------------
  auto bouncingBall = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);

  // add the dynamical system in the non smooth dynamical system
  bouncingBall->insertDynamicalSystem(ball);

  // link the interaction and the dynamical system
  bouncingBall->link(inter, ball);

  // ------------------
  // --- Simulation ---
  // ------------------

  // -- (1) OneStepIntegrators --
  auto OSI = std::make_shared<MoreauJeanOSI>(theta);

  // -- (2) Time discretisation --
  auto t = std::make_shared<TimeDiscretisation>(t0, h);

  // -- (3) one step non smooth problem
  auto osnspb = std::make_shared<LCP>();

  // -- (4) Simulation setup with (1) (2) (3)
  auto s = std::make_shared<TimeStepping>(bouncingBall, t, OSI, osnspb);
  s->associate(OSI, ball);

  // =========================== End of model definition ===========================

  // ================================= Computation =================================

  Siconos::save(s, BBxml);

  auto simFromFile = Siconos::load(BBxml);
  auto bouncingBallFromFile = simFromFile->nonSmoothDynamicalSystem();

  CPPUNIT_ASSERT((bouncingBallFromFile->t0() == bouncingBall->t0()));
  // in depth comparison?

  // now we should try to run the bouncing ball from file

  // BUT: non serialized members => must be initialized or serialized
}

void KernelTest::t6() {
  auto sim = Siconos::load(BBxml);

  try {
    auto bouncingBall = sim->nonSmoothDynamicalSystem();

    double T = bouncingBall->finalT();
    double t0 = bouncingBall->t0();
    double h = sim->timeStep();
    int N = (int)((T - t0) / h);  // Number of time steps

    auto dsg = bouncingBall->topology()->dSG(0);

    auto ball = std::static_pointer_cast<siconos::modeling::LagrangianDS>(
        dsg->bundle(*(dsg->begin())));

    auto s = std::static_pointer_cast<TimeStepping>(sim);
    std::shared_ptr<siconos::modeling::Interaction> inter;
    InteractionsGraph::VIterator ui, uiend;
    auto indexSet0 = bouncingBall->topology()->indexSet(0);
    for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
      inter = indexSet0->bundle(*ui);

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 5;
    siconos::algebra::SimpleMatrix dataPlot(N + 1, outputSize);

    auto q = ball->q();
    auto v = ball->velocity();
    auto p = ball->p(1);
    auto lambda = inter->lambda(1);

    dataPlot(0, 0) = bouncingBall->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*v)(0);
    dataPlot(0, 3) = (*p)(0);
    dataPlot(0, 4) = (*lambda)(0);
    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    while (s->hasNextEvent()) {
      s->computeOneStep();

      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0);
      dataPlot(k, 4) = (*lambda)(0);
      s->nextStep();
      k++;
    }
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write("result.dat", "ascii", dataPlot, "noDim");
    // Comparison with a reference file
    siconos::algebra::SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
    ioMatrix::read("result.ref", "ascii", dataPlotRef);

    if ((dataPlot - dataPlotRef).normInf() > 1e-12) {
      std::cout << "Warning. The results is rather different from the reference file :"
                << (dataPlot - dataPlotRef).normInf() << std::endl;
      CPPUNIT_ASSERT(false);
    }

  }

  catch (...) {
    Siconos::exception::process();
    CPPUNIT_ASSERT(false);
  }
}

#ifdef HAVE_SICONOS_MECHANICS
#include <Disk.hpp>

void KernelTest::t7() {
  std::shared_ptr < siconos::modeling::DynamicalSystem ds1, ds2;

  // Must be size=1, cannot deserialize a LagrangianDS with _ndof==0
  auto q = std::make_shared<siconos::algebra::SiconosVector>(1);
  auto v = std::make_shared<siconos::algebra::SiconosVector>(1);

  ds1.reset = std::make_shared<Disk(1, 1, q, v));

  ds2.reset = std::make_shared<Disk(2, 2, q, v));

  std::ofstream ofs("Kernelt7.xml");
  {
    boost::archive::xml_oarchive oa(ofs);
    siconos_io_register_Mechanics(oa);
    oa << NVP(ds1);
  }

  std::ifstream ifs("Kernelt7.xml");
  {
    boost::archive::xml_iarchive ia(ifs);
    siconos_io_register_Mechanics(ia);
    ia >> NVP(ds2);
  }

  CPPUNIT_ASSERT(std::static_pointer_cast<Disk>(ds1)->getRadius() ==
                 std::static_pointer_cast<Disk>(ds2)->getRadius());
}

void KernelTest::t8() {
  std::shared_ptr < siconos::modeling::DynamicalSystem ds1, ds2;

  auto q = std::make_shared<siconos::algebra::SiconosVector>(3);
  auto v = std::make_shared<siconos::algebra::SiconosVector>(3);

  (*q)(0) = 0.;
  (*q)(1) = 1.;
  (*q)(2) = 1.;

  (*v)(0) = 0;
  (*v)(1) = 0;
  (*v)(2) = 10.;

  auto nsds = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(0, 10);

  ds1 = std::make_shared<Disk(1, 1, q, v));
  ds2 = std::make_shared<Disk(2, 2, q, v));

  nsds->insertDynamicalSystem(ds1);
  nsds->insertDynamicalSystem(ds2);

  MechanicsIO IO;

  auto positions = IO.positions(*nsds);
  auto velocities = IO.velocities(*nsds);

  // ids
  CPPUNIT_ASSERT((*positions)(0, 0) == 0);
  CPPUNIT_ASSERT((*velocities)(0, 0) == 0);
  CPPUNIT_ASSERT((*positions)(0, 0) == ds1->number());
  CPPUNIT_ASSERT((*velocities)(0, 0) == ds1->number());

  CPPUNIT_ASSERT((*positions)(1, 0) == 1);
  CPPUNIT_ASSERT((*velocities)(1, 0) == 1);
  CPPUNIT_ASSERT((*positions)(1, 0) == ds2->number());
  CPPUNIT_ASSERT((*velocities)(1, 0) == ds2->number());

  CPPUNIT_ASSERT((*positions)(0, 1) == 0.);
  CPPUNIT_ASSERT((*velocities)(0, 1) == 0.);
  CPPUNIT_ASSERT((*positions)(0, 2) == 1.);
  CPPUNIT_ASSERT((*positions)(1, 2) == 1.);
  CPPUNIT_ASSERT((*velocities)(0, 3) == 10.);
  CPPUNIT_ASSERT((*velocities)(1, 3) == 10.);
}
#endif

void KernelTest::t9() {
  try {
    // Serialize and deserialize an NSDS with T=inf
    // (a possible failure case for xml archives)
    double t0 = 0.0;
    double T = std::numeric_limits<double>::infinity();
    auto  nsds1 = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem(t0, T));
    std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> nsds2;

    std::ofstream ofs("Kernelt9.xml");
    {
      boost::archive::xml_oarchive oa(ofs);
      siconos_io_register_Kernel(oa);
      oa << NVP(nsds1);
    }

    std::ifstream ifs("Kernelt9.xml");
    {
      boost::archive::xml_iarchive ia(ifs);
      siconos_io_register_Kernel(ia);
      ia >> NVP(nsds2);
    }

    CPPUNIT_ASSERT(nsds1->finalT() == nsds2->finalT());
  } catch (...) {
    Siconos::exception::process();
    CPPUNIT_ASSERT(false);
  }
}
