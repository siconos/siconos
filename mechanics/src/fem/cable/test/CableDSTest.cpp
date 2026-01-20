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
#include "Interaction.hpp"
#include "LCP.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "LagrangianLinearTIR.hpp"
#include "MoreauJeanOSI.hpp"
#include "NewtonImpactNSL.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Rope.h"  // IWYU pragma: keep
#include "SiconosVector.hpp"
#include "StorageTools.hpp"
#include "TimeDiscretisation.hpp"
#include "TimeStepping.hpp"
#include "TransportCableManager.h"
#include "TransportCableProfile.h"
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
  std::string modelFile = "data/model.origin.json";
  siconos::fem::cable::TransportCableManager manager;
  manager.importModel(modelFile);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testReadModel: model loading ok", manager.isReady(), true);
  std::cout << "✅ test readModel passed.\n";
}

void CableDSTest::testBuildInitialProfile() {
  // ---- Build a profile from a TransportCableModel read in a json file ----
  //
  // Model from a file
  std::string modelFile = "data/model.origin.json";
  auto json_input = siconos::fem::cable::tools::load_json_file(modelFile);
  siconos::fem::cable::TransportCableModel model{json_input};
  siconos::fem::cable::TransportCableResult results;
  siconos::fem::cable::TransportCableProfile profile{model, results};

  // Computes a profile (applies catenary equation)
  int nb_nodes = 50;   // Catenary, number of nodes per piece of rope
  double tol = 1e-10;  // Tol. used in Newton-Raphson for catenary equation
  int nmax = 20;       // Newton-Raphson, max number of iterations
  profile.computeInitialProfile(nb_nodes, tol, nmax);

  // At this stage, positions are saved in each Rope, for all points.
  // To get them and for comparison, we use output to json and 'reload' into a
  // siconos::algebra::SiconosVector.
  // 1. Save ropeways variables into json file
  nlohmann::ordered_json out;
  results.to_json(out);
  // 2. Upload to q1 and q2
  auto q1 = siconos::algebra::io::readVectorFromJson(out["ropes_up"]["q"]);
  auto q2 = siconos::algebra::io::readVectorFromJson(out["ropes_down"]["q"]);
  std::ofstream file("tmp_out.json");
  file << out.dump(4);
  // 3. Read reference
  std::ifstream in("data/bouquetins_ref.json");
  json reader;
  in >> reader;
  auto qref1 = siconos::algebra::io::readVectorFromJson(reader["ropes_up"]["q"]);
  auto qref2 = siconos::algebra::io::readVectorFromJson(reader["ropes_down"]["q"]);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testBuildInitialProfile: check catenary",
                               (qref1.size() == q1.size()), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testBuildInitialProfile: check catenary", qref1.isApprox(q1),
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testBuildInitialProfile: check catenary", qref2.isApprox(q2),
                               true);
  std::cout << "✅ test build initial profile passed.\n";
}

void CableDSTest::testInitializeFEM() {
  // 1. Build a profile from a TransportCableModel read in a json file ----
  // 2. Initialize its profile
  // 3. Initialize the FE mesh and the variables at mesh points/elements

  std::string modelFile = "data/model.origin.json";
  auto json_input = siconos::fem::cable::tools::load_json_file(modelFile);
  siconos::fem::cable::TransportCableModel model{json_input};
  siconos::fem::cable::TransportCableResult results;
  siconos::fem::cable::TransportCableProfile profile{model, results};
  int nb_nodes = 50;   // Catenary, number of nodes per rope span
  double tol = 1e-10;  // Tol. used in Newton-Raphson for catenary equation
  int nmax = 20;       // Newton-Raphson, max number of iterations
  profile.computeInitialProfile(nb_nodes, tol, nmax);

  int numberOfElements = 1400;  // FEM number of nodes
  double eps = 0.1;
  tol = 1e-3;  // tolerance used to activate constraints
  profile.initializeFEM(numberOfElements, eps, tol);
  nlohmann::ordered_json out;
  // Now we dump positions into json ...
  results.to_json(out);
  std::ofstream file("tmp_out2.json");
  file << out.dump(4);

  // And read the resulting dof vector
  auto positions = siconos::algebra::io::readVectorFromJson(out["q"]);

  // Read the reference and compare
  std::ifstream in("data/bouquetins_ref.json");
  json reader;
  in >> reader;
  auto positions_ref = siconos::algebra::io::readVectorFromJson(reader["q"]);
  std::cout << "norm diff " << (positions - positions_ref).norm() << "\n";
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testBuildInitialProfile:  check fem initialisation",
                               positions.isApprox(positions_ref, 1e-10), true);
  std::cout << "✅ test initialize FEM passed.\n";

  // // compare results to a reference
}

void CableDSTest::testComputeDS() {
  // std::cout << "Start testComputeDS ...\n";
  // json args;
  // std::string outFile = "results/compute.origin.json";
  // nlohmann::ordered_json out;
  // std::string modelFile = "data/model.origin.json";
  // auto manager = std::make_shared<TransportCableManager>();
  // auto res = manager->importModel(args, modelFile);

  // res = manager->initializeFEM(args, outFile, out);
  // CPPUNIT_ASSERT_NOT_EQUAL(" testComputeDS: compute FAIL", res == EXIT_FAILURE, true);
  // std::cout << "End testComputeDS ...\n";
  // // // compare results to a reference

  // Load the model
  std::string modelFile = "data/model.origin.json";

  json args;
  siconos::fem::cable::TransportCableManager manager;
  manager.importModel(modelFile);

  // Compute, simulation
  std::string outFile = "results/compute.origin.json";
  nlohmann::ordered_json out;
  // manager.computeFEM(args, outFile, out);
  std::cout << "✅ test computeDS passed.\n";
}

void CableDSTest::testComputeBouncingBall() {
  // ================= Creation of the model =======================

  // User-defined main parameters
  int nDof = 3;                // degrees of freedom for the ball
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

  siconos::algebra::SiconosMatrix mass{nDof, nDof};
  mass(0, 0) = m;
  mass(1, 1) = m;
  mass(2, 2) = 2. / 5 * m * R * R;

  // -- Initial positions and velocities --
  siconos::algebra::SiconosVector q0{nDof};
  siconos::algebra::SiconosVector v0{nDof};
  q0(0) = position_init;
  v0(0) = velocity_init;

  // -- The dynamical system --
  auto ball = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(
      q0, v0, mass, siconos::algebra::alias_t);

  // -- Set external forces (weight) --
  siconos::algebra::SiconosVector weight{nDof};
  weight(0) = -m * g;
  ball->setConstantFext(weight, siconos::algebra::copy_t);

  // --------------------
  // --- Interactions ---
  // --------------------

  // -- nslaw --
  double e = 0.9;

  // Interaction ball-floor
  //
  auto H = std::make_shared<siconos::algebra::SiconosMatrix>(1, nDof);
  (*H)(0, 0) = 1.0;

  auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
  auto relation = std::make_shared<siconos::modeling::LagrangianLinearTIR>(*H);

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
  siconos::algebra::SiconosMatrix dataPlot(N + 1, outputSize);

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
  std::cout << "✅ test computeBouncingBall passed.\n";
}

void CableDSTest::testNoFext() {
  // auto cable = std::make_shared<siconos::mechanics::fem::CableDS>(q0, v0, mass, 14000, 1);

  // CPPUNIT_ASSERT_EQUAL_MESSAGE(" testNoFext: ", cable->fext() == nullptr, true);
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
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(" testCFext 1: ", cable->fext() != nullptr, true);
  // auto fext = cable->fext();
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
  // myforces); auto fext = cable->fext(); CPPUNIT_ASSERT_EQUAL_MESSAGE(" testVFext 1: ", fext
  // != nullptr, true);

  // cable->computeFext(3.);
  // for (auto val : *fext)
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE(" testVFext 2: ", val == cos(3.), true);

  // auto positions = cable->q();
  // auto velocities = cable->velocity();
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(" testVFext 4: ", positions->size() == cable->dimension(),
  //                              true);
  // cable->computeTotalForces(velocities, positions, 5.);
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
