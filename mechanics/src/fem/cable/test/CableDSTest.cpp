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

#include "CableDSTest.hpp"

#include "CableTools.h"
#include "Rope.h"
#include "SiconosKernel.hpp"
#include "SiconosVector.hpp"
#include "SiconosVectorIterator.hpp"
#include "TransportCableProfil.h"
#include "io.hpp"

using namespace std;

#include <nlohmann/json.hpp>

// for convenience
using json = nlohmann::json;

#include <string>

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(CableDSTest);

void CableDSTest::setUp() {}

void CableDSTest::tearDown() {}

void CableDSTest::testReadModel() {
  std::cout << "Start testReadModel ...\n";
  std::string modelFile = "data/model.origin.json";
  json args;

  auto manager = std::make_shared<siconos::fem::cable::TransportCableManager>();
  auto res = manager->importModel(args, modelFile);
  CPPUNIT_ASSERT_NOT_EQUAL("testReadModel: cannot load the model", res == EXIT_FAILURE, true);
  std::cout << "End testReadModel ...\n";
}

void CableDSTest::testBuildInitialProfile() {
  std::cout << "Start testBuildInitialProfile ...\n";
  std::string outFile = "results/compute.origin.json";
  std::string modelFile = "data/model.origin.json";
  auto model = std::make_shared<siconos::fem::cable::TransportCableModel>(modelFile);

  auto results = std::make_shared<siconos::fem::cable::TransportCableResult>();

  auto profil = std::make_shared<siconos::fem::cable::TransportCableProfil>(*model, *results);

  int nb_nodes = 50;   // Catenary, number of nodes per rope span
  double tol = 1e-10;  // Tol. used in Newton-Raphson for catenary equation
  int nmax = 20;       // Newton-Raphson, max number of iterations
  profil->computeInitialProfil(nb_nodes, tol, nmax);

  // At this stage, positions are saved in each Rope, for all points.
  // To get them and for comparison, we use output to json and 'reload' into a
  // siconos::algebra::SiconosVector. Save ropeways variables into json file
  ojson out;
  results->to_json(out, "ropeway");
  auto q1 = siconos::algebra::io::readVectorFromJson(out["rope1"]["q"]);
  auto q2 = siconos::algebra::io::readVectorFromJson(out["rope2"]["q"]);
  // Read reference
  std::ifstream in("data/bouquetins_ref.json");
  json reader;
  in >> reader;
  auto qref1 = siconos::algebra::io::readVectorFromJson(reader["rope1"]["q"]);
  auto qref2 = siconos::algebra::io::readVectorFromJson(reader["rope2"]["q"]);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testBuildInitialProfile: check catenary", (qref1->size() == q1->size()),
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testBuildInitialProfile: check catenary", ((*qref1) == (*q1)),
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testBuildInitialProfile: check catenary", ((*qref2) == (*q2)),
                              true);

  nb_nodes = 1400;  // FEM number of nodes
  double eps = 0.1;
  tol = 1e-3;  // tolerance used to activate constraints
  profil->computeFEM(nb_nodes, eps, tol);

  // Now we dump positions into json ...
  results->to_json(out);
  // And build the dof vector
  auto positions = siconos::algebra::io::readVectorFromJson(out["q"]);

  // Read the reference and compare

  auto positions_ref = siconos::algebra::io::readVectorFromJson(reader["q"]);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testBuildInitialProfile:  check fem",
                               ((*positions_ref) == (*positions)), true);
  std::cout << "End testBuildInitialProfile ...\n";
  // // compare results to a reference
}

void CableDSTest::testComputeDS() {
  // std::cout << "Start testComputeDS ...\n";
  // json args;
  // std::string outFile = "results/compute.origin.json";
  // ojson out;
  // std::string modelFile = "data/model.origin.json";
  // auto manager = std::make_shared<TransportCableManager>();
  // auto res = manager->importModel(args, modelFile);

  // res = manager->computeFEM(args, outFile, out);
  // CPPUNIT_ASSERT_NOT_EQUAL(" testComputeDS: compute FAIL", res == EXIT_FAILURE, true);
  // std::cout << "End testComputeDS ...\n";
  // // // compare results to a reference

  // Load the model
  std::string modelFile = "data/model.origin.json";

  json args;
  auto manager = std::make_shared<siconos::fem::cable::TransportCableManager>();
  int res = manager->importModel(args, modelFile);
  CPPUNIT_ASSERT_NOT_EQUAL(" setUp: cannot load the model", res == EXIT_FAILURE, true);

  // Compute, simulation
  std::string outFile = "results/compute.origin.json";
  ojson out;
  res = manager->computeFEM(args, outFile, out);
  CPPUNIT_ASSERT_NOT_EQUAL(" testComputeDS: compute FAIL", res == EXIT_FAILURE, true);
}

void CableDSTest::testComputeBouncingBall() {
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

  cout << "====> Model loading ..." << endl;

  auto Mass = std::make_shared<siconos::algebra::SimpleMatrix>(nDof, nDof);
  (*Mass)(0, 0) = m;
  (*Mass)(1, 1) = m;
  (*Mass)(2, 2) = 2. / 5 * m * R * R;

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

  auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, relation);

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
  auto OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);

  // -- (2) Time discretisation --
  auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

  // -- (3) one step non smooth problem
  auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();

  // -- (4) Simulation setup with (1) (2) (3)
  auto s = std::make_shared<siconos::simulation::TimeStepping>(bouncingBall, t, OSI, osnspb);

  // =========================== End of model definition ===========================

  // ================================= Computation =================================

  int N = ceil((T - t0) / h);  // Number of time steps

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
  cout << "====> Start computation ... " << endl;
  // ==== Simulation loop - Writing without explicit event handling =====
  int k = 1;
  auto start = std::chrono::system_clock::now();
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
  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
  cout << "Computation time : " << elapsed << " ms" << endl;
}

void CableDSTest::testNoFext() {
  // auto cable = std::make_shared<siconos::mechanics::fem::CableDS>(q0, v0, mass, 14000, 1);

  // CPPUNIT_ASSERT_EQUAL_MESSAGE(" testNoFext: ", cable->fExt() == nullptr, true);
  //  ...
}

void CableDSTest::testConstantFext() {
  // auto externalForces = std::make_shared<siconos::algebra::SiconosVector>(ndof);
  //// fill external forces as you want ...
  // for (auto& val : *externalForces) {
  //   val = 12;
  // }

  // auto cable = std::make_shared<siconos::mechanics::fem::CableDS>(q0, v0, mass, 14000, 1);
  // cable->setFExtPtr(externalForces);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(" testCFext 1: ", cable->fExt() != nullptr, true);
  // auto fext = cable->fExt();
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(" testCFext 2: ", fext != nullptr, true);
  // for (auto val : *fext) {
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE(" testCFext 3: ", val == 12, true);
  // }
  //// ...
}

void CableDSTest::testVariableFext()

{
  //// A lambda function example, used to compute external forces
  // auto myforces = [](double time, std::shared_ptr<siconos::algebra::SiconosVector> result) {
  //   assert(result);
  //   for (auto& v : *result)
  //     v = cos(time);
  // };

  // auto cable = std::make_shared<siconos::mechanics::fem::CableDS>(q0, v0, mass, 14000, 1,
  // myforces); auto fext = cable->fExt(); CPPUNIT_ASSERT_EQUAL_MESSAGE(" testVFext 1: ", fext
  // != nullptr, true);

  // cable->computeFExt(3.);
  // for (auto val : *fext)
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE(" testVFext 2: ", val == cos(3.), true);

  // auto positions = cable->q();
  // auto velocities = cable->velocity();
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(" testVFext 4: ", positions->size() == cable->dimension(),
  //                              true);
  // cable->computeForces(5., positions, velocities);
  // for (auto val : *fext)
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE(" testVFext 5: ", val == cos(5.), true);

  // Example on how to proceed to compare matrices:
  // matrix is the new matrix to be checked. SomeReferenceFile.ref contains a former matrix
  // saved in a file and used as reference.
  // auto matrix = mass;
  // double eps = 1e-12;
  // auto error = ioMatrix::compareRefFile(*matrix, "SomeReferenceFile.ref", eps);
  // Check error ...
}
