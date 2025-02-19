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
#include "TransportCableManager.h"

#include <nlohmann/json.hpp>

#include "CableCollisionManager.h"
#include "CableDS.hpp"
#include "CableTools.h"
#include "FrictionContact.hpp"
#include "MoreauJeanOSI.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "OneStepNSProblem.hpp"
#include "Rope.h"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "Support.h"
#include "TimeDiscretisation.hpp"
#include "TimeStepping.hpp"
#include "TransportCableProfil.h"

using json = nlohmann::json;
using ojson = nlohmann::ordered_json;

siconos::fem::cable::TransportCableManager::TransportCableManager(
    const std::string &a_filename)
    : TransportCableManager{} {
  m_model.from_file(a_filename);
}

int siconos::fem::cable::TransportCableManager::importModel(const json &a_input,
                                                            const std::string &a_filename) {
  int res = EXIT_FAILURE;
  if (a_input.is_null()) {
    res = m_model.from_file(a_filename);
  } else {
    res = m_model.from_json(a_input);
  }
  return res;
}

void siconos::fem::cable::TransportCableManager::computeFEM(const json &a_args,
                                                            const std::string &a_outfile,
                                                            ojson &output) {
  // Ensures the model is valid
  assert(m_model.isLoaded());

  // Creates the profile
  TransportCableProfil profile(m_model, m_results);

  // Reads method from json ()
  std::string method = "all";  // default method if json is empty
  method = siconos::fem::cable::tools::getParam(a_args, "compute", method);

  // Initialize supports and ropes
  profile.computeInitialProfile(siconos::fem::cable::tools::getParam(a_args, "nb_node0", 50),
                                siconos::fem::cable::tools::getParam(a_args, "tol", 1e-20),
                                siconos::fem::cable::tools::getParam(a_args, "nmax", 20));

  //
  profile.computeFEM(siconos::fem::cable::tools::getParam(a_args, "nb_node", 1400),
                     siconos::fem::cable::tools::getParam(a_args, "eps", 0.1),
                     siconos::fem::cable::tools::getParam(a_args, "tol_contact", 1e-3));
#ifndef NSICONOS
  if (method == "dynamics") {
    computeDS();
  }
#endif
  exportTC(a_args, a_outfile, output);
}

int siconos::fem::cable::TransportCableManager::exportTC(const json &a_args,
                                                         const std::string &a_outfile,
                                                         ojson &output) {
  auto vOption = siconos::fem::cable::tools::getParam(a_args, "export", (std::string) "all");
  m_results.exportTC(a_outfile, output, vOption);

  return EXIT_SUCCESS;
}

int siconos::fem::cable::TransportCableManager::simulation(const json &a_model,
                                                           const json &a_args,
                                                           const std::string &a_filename,
                                                           const std::string &a_outfile,
                                                           ojson &output) {
  int vRet = importModel(a_model, a_filename);
  if (vRet == EXIT_SUCCESS) {
    computeFEM(a_args, a_outfile, output);
  }
  return vRet;
}

void siconos::fem::cable::TransportCableManager::computeDS(double a_tolContact, double a_mus,
                                                           double a_mup) {
  // model is loaded
  // q0 must be computed
  // q0 = q
  int ndof = m_results.q.size() * 3;
  if (not m_results.q0) m_results.q0 = std::make_shared<siconos::algebra::SiconosVector>(ndof);
  if (not m_results.v0) m_results.v0 = std::make_shared<siconos::algebra::SiconosVector>(ndof);

  siconos::fem::cable::tools::pointsToSiconosVector(m_results.q, m_results.q0);

  double rho = m_model.mechanicalProperties().get_rho();

  // compute mass -> results.mass
  compute_mass(m_results.elem_length, rho);

  // compute external forces -> results.fext
  compute_external_load(m_results.elem_length, rho);

  // create dynamics model
  auto cable = std::make_shared<CableDS>(
      Eigen::Ref<siconos::algebra::SiconosVector>(*m_results.q0),
      Eigen::Ref<siconos::algebra::SiconosVector>(*m_results.v0),
      Eigen::Ref<siconos::algebra::SiconosMatrix>(*m_results.mass),
      m_model.mechanicalProperties().get_EA(), m_results.elem_length);

  cable->setConstantFext(Eigen::Ref<siconos::algebra::SiconosVector>(*m_results.fext));

  // contact conditions
  auto collisions =
      std::make_shared<CableCollisionManager>(cable, m_results.supports, a_tolContact);

  // frictions
  int ns = m_results.supports.size();
  for (int i = 0; i < ns; i++) {
    if (i == m_results.puller12idx || i == m_results.puller21idx) {
      m_results.supports[i]->InitFriction(a_mup);
    } else {
      m_results.supports[i]->InitFriction(a_mus);
    }
  }

  // simulation
  double t0 = 0;       // initial computation time
  double T = 1e-02;    // final computation time
  double h = 1e-05;    // time step
  double theta = 1.0;  // theta for MoreauJeanOSI integrator

  // model
  auto dynamicModel = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);

  // add the dynamical system in the non smooth dynamical system
  dynamicModel->insertDynamicalSystem(cable);

  // -- (1) OneStepIntegrators --
  auto OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);
  OSI->setIsWSymmetricDefinitePositive(true);
  OSI->setGamma(0.0);

  // -- (2) Time discretisation --
  auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

  // -- (3) one step non smooth problem
  auto osnspb = std::make_shared<siconos::nonsmooth_formulations::FrictionContact>(2);
  // SolverOptions* options = osnspb->numericsSolverOptions().get();
  // options->dparam[SICONOS_DPARAM_TOL] = 1e-10;

  // -- (4) Simulation setup with (1) (2) (3)
  auto s = std::make_shared<siconos::simulation::TimeStepping>(dynamicModel, t, OSI, osnspb);
  s->insertInteractionManager(collisions);

  // --- Time loop ---
  // int k = 1;
  //  auto start = std::chrono::system_clock::now();

  while (s->hasNextEvent()) {
    s->computeOneStep();

    /*if (k % 1 == 0)
            writeDisplacementforPython(mesh, femodel, q, filename);*/

    s->nextStep();

    // k++;
  }
  // end = std::chrono::system_clock::now();
  // auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
}

void siconos::fem::cable::TransportCableManager::compute_mass(double a_length, double a_rho) {
  /*get_mass_damp(elem_length,elem_rho,damping)
        where
        elem_length is the initial length of each element
        elem_rho is the vector which e-th element is the linear density of element number e

  */
  int ndof = m_results.q0->size();
  if (not m_results.mass)
    m_results.mass =
        std::make_shared<siconos::algebra::SiconosMatrix>(ndof, ndof);  // FP: must be SPARSE
  double k = a_rho * a_length / 3.0;
  for (auto i = 0; i < ndof - 3; i++) {
    m_results.mass->setValue(i, i, 4 * k);
    m_results.mass->setValue(i, i + 3, k);
    m_results.mass->setValue(i + 3, i, k);
  }
  for (size_t i = 0; i < 3; i++) {
    m_results.mass->setValue(i + ndof - 3, i + ndof - 3, 4 * k);
    m_results.mass->setValue(i, i + ndof - 3, k);
    m_results.mass->setValue(i + ndof - 3, i, k);
  }
}

void siconos::fem::cable::TransportCableManager::compute_external_load(double a_length,
                                                                       double a_rho) {
  auto ndof = m_results.q0->size();
  if (not m_results.fext)
    m_results.fext = std::make_shared<siconos::algebra::SiconosVector>(ndof);

  double k = -9.81 * a_rho * a_length;
  for (auto i = 2; i < ndof; i += 3) {
    m_results.fext->setValue(i, k);
  }
}
