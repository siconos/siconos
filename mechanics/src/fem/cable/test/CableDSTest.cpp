/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#include "SiconosVector.hpp"
#include "SiconosVectorIterator.hpp"
#include "TransportCableProfil.h"
#include "ioMatrix.hpp"
#include "ioVector.hpp"

#include <nlohmann/json.hpp>

// for convenience
using json = nlohmann::json;

#include <string>

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)                                       \
  if ((alpha) == (omega))                                                                     \
    CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(CableDSTest);

void CableDSTest::setUp() {}

void CableDSTest::tearDown() {}

void CableDSTest::testReadModel()
{
  std::cout << "Start testReadModel ...\n";
  std::string modelFile = "data/model.origin.json";
  json args;

  auto manager = std::make_shared<TransportCableManager>();
  auto res = manager->importModel(args, modelFile);
  CPPUNIT_ASSERT_NOT_EQUAL("testReadModel: cannot load the model", res == EXIT_FAILURE, true);
  std::cout << "End testReadModel ...\n";
}

void CableDSTest::testBuildInitialProfile()
{
  std::cout << "Start testBuildInitialProfile ...\n";
  std::string outFile = "results/compute.origin.json";
  std::string modelFile = "data/model.origin.json";
  auto model = std::make_shared<TransportCableModel>(modelFile);

  auto results = std::make_shared<TransportCableResult>();

  auto profil = std::make_shared<TransportCableProfil>(*model, *results);

  int nb_nodes = 50;  // Catenary, number of nodes per rope span
  double tol = 1e-10; // Tol. used in Newton-Raphson for catenary equation
  int nmax = 20;      // Newton-Raphson, max number of iterations
  profil->computeInitialProfil(nb_nodes, tol, nmax);

  // At this stage, positions are saved in each Rope, for all points.
  // To get them and for comparison, we use output to json and 'reload' into a SiconosVector.
  // Save ropeways variables into json file
  ojson out;
  results->to_json(out, "ropeway");
  auto q1 = ioVector::readVectorFromJson(out["rope1"]["q"]);
  auto q2 = ioVector::readVectorFromJson(out["rope2"]["q"]);
  // Read reference
  std::ifstream in("data/bouquetins_ref.json");
  json reader;
  in >> reader;
  auto qref1 = ioVector::readVectorFromJson(reader["rope1"]["q"]);
  auto qref2 = ioVector::readVectorFromJson(reader["rope2"]["q"]);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testBuildInitialProfile: check catenary", ((*qref1) == (*q1)) , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testBuildInitialProfile: check catenary", ((*qref2) == (*q2)) , true);
  
  nb_nodes = 1400; // FEM number of nodes
  double eps = 0.1;
  tol = 1e-3; // tolerance used to activate constraints
  profil->computeFEM(nb_nodes, eps, tol);

  // Now we dump positions into json ...
  results->to_json(out);
  // And build the dof vector
  auto positions = ioVector::readVectorFromJson(out["q"]);

  // Read the reference and compare

  auto positions_ref = ioVector::readVectorFromJson(reader["q"]);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testBuildInitialProfile:  check fem", ((*positions_ref) == (*positions)) , true);
  std::cout << "End testBuildInitialProfile ...\n";
  // // compare results to a reference
}

void CableDSTest::testComputeDS()
{
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
}

void CableDSTest::testNoFext()
{

  // auto cable = std::make_shared<siconos::mechanics::fem::CableDS>(q0, v0, mass, 14000, 1);

  // CPPUNIT_ASSERT_EQUAL_MESSAGE(" testNoFext: ", cable->fExt() == nullptr, true);
  //  ...
}

void CableDSTest::testConstantFext()
{
  // auto externalForces = std::make_shared<SiconosVector>(ndof);
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
  // auto myforces = [](double time, std::shared_ptr<SiconosVector> result) {
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
