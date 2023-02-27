#include "TransportCableManager.h"

#include "TransportCableProfil.h"

#ifndef NSICONOS

#include "CableCollisionManager.h"
#include "CableDS.hpp"
#include "FrictionContact.hpp"
#include "MoreauJeanOSI.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "OneStepNSProblem.hpp"
#include "SimpleMatrix.hpp"
#include "TimeDiscretisation.hpp"
#include "TimeStepping.hpp"
#include "CableTools.h"
#endif
#include <nlohmann/json.hpp>

using json = nlohmann::json;
using ojson = nlohmann::ordered_json;

// #include "TCException.h"


TransportCableManager::TransportCableManager(const std::string &a_filename): TransportCableManager{}
{
  m_model.from_file(a_filename);
}

int TransportCableManager::importModel(const json &a_input, const string &a_filename)
{
  int res = EXIT_FAILURE;
  if (a_input.is_null()) {
    res = m_model.from_file(a_filename);
  }
  else {
    res = m_model.from_json(a_input);
  }
  return res;
}

int TransportCableManager::computeFEM(const json &a_args, const string &a_outfile,
                                      ojson &output)
{
  if (m_model.isLoaded()) {
    TransportCableProfil P(m_model, m_results);
    string method = "all";
    method = getParam(a_args, "compute", method);

    P.computeInitialProfil(getParam(a_args, "nb_node0", 50), getParam(a_args, "tol", 1e-20),
                           getParam(a_args, "nmax", 20));

    P.computeFEM(getParam(a_args, "nb_node", 1400), getParam(a_args, "eps", 0.1),
                 getParam(a_args, "tol_contact", 1e-3));
#ifndef NSICONOS
    if (method == "dynamics") {
      computeDS();
    }
#endif
    exportTC(a_args, a_outfile, output);

    return EXIT_SUCCESS;
  }
  else {
    // throw TCException("Load a model before compute");
  }
  return EXIT_FAILURE;
}

int TransportCableManager::exportTC(const json &a_args, const string &a_outfile, ojson &output)
{
  string vOption = getParam(a_args, "export", (string) "all");
  m_results.exportTC(a_outfile, output, vOption);

  return EXIT_SUCCESS;
}

int TransportCableManager::simulation(const json &a_model, const json &a_args,
                                      const string &a_filename, const string &a_outfile,
                                      ojson &output)
{
  int vRet = importModel(a_model, a_filename);
  if (vRet == EXIT_SUCCESS) {
    vRet = computeFEM(a_args, a_outfile, output);
  }
  return vRet;
}

#ifndef NSICONOS

void TransportCableManager::computeDS(double a_tolContact, double a_mus, double a_mup)
{
  // model is loaded
  // q0 must be computed
  // q0 = q
  int ndof = m_results.q.size() * 3;
  if(not m_results.q0)
    m_results.q0 = std::make_shared<SiconosVector>(ndof);
  if(not m_results.v0)
    m_results.v0 = std::make_shared<SiconosVector>(ndof);

  siconos::mechanics::cables::tools::pointsToSiconosVector(m_results.q, m_results.q0);

  double rho = m_model.get_cable().get_rho();

  // compute mass -> results.mass
  compute_mass(m_results.elem_length, rho);

  // compute external forces -> results.b
  compute_external_load(m_results.elem_length, rho);

  // create dynamics model
  auto cable = std::make_shared<siconos::mechanics::fem::CableDS>(
      m_results.q0, m_results.v0, m_results.mass, m_model.get_cable().get_EA(),
      m_results.elem_length);

  cable->setFExtPtr(m_results.b);

  // contact conditions
  auto collisions =
      std::make_shared<CableCollisionManager>(cable, m_results.supports, a_tolContact);

  // frictions
  int ns = m_results.supports.size();
  for (int i = 0; i < ns; i++) {
    if (i == m_results.puller12idx || i == m_results.puller21idx) {
      m_results.supports[i]->InitFriction(a_mup);
    }
    else {
      m_results.supports[i]->InitFriction(a_mus);
    }
  }

  // simulation
  double t0 = 0;      // initial computation time
  double T = 1e-02;   // final computation time
  double h = 1e-05;   // time step
  double theta = 1.0; // theta for MoreauJeanOSI integrator

  // model
  std::shared_ptr<NonSmoothDynamicalSystem> dynamicModel =
      std::make_shared<NonSmoothDynamicalSystem>(t0, T);

  // add the dynamical system in the non smooth dynamical system
  dynamicModel->insertDynamicalSystem(cable);

  // -- (1) OneStepIntegrators --
  std::shared_ptr<MoreauJeanOSI> OSI = std::make_shared<MoreauJeanOSI>(theta);
  OSI->setIsWSymmetricDefinitePositive(true);
  OSI->setGamma(0.0);

  // -- (2) Time discretisation --
  std::shared_ptr<TimeDiscretisation> t = std::make_shared<TimeDiscretisation>(t0, h);

  // -- (3) one step non smooth problem
  std::shared_ptr<OneStepNSProblem> osnspb = std::make_shared<FrictionContact>(2);
  // SolverOptions* options = osnspb->numericsSolverOptions().get();
  // options->dparam[SICONOS_DPARAM_TOL] = 1e-10;

  // -- (4) Simulation setup with (1) (2) (3)
  std::shared_ptr<TimeStepping> s =
      std::make_shared<TimeStepping>(dynamicModel, t, OSI, osnspb);
  s->insertInteractionManager(collisions);

  // --- Time loop ---
  int k = 1;
  //  auto start = std::chrono::system_clock::now();

  while (s->hasNextEvent()) {
    s->computeOneStep();

    /*if (k % 1 == 0)
            writeDisplacementforPython(mesh, femodel, q, filename);*/

    s->nextStep();

    k++;
  }
  // end = std::chrono::system_clock::now();
  // auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
}

void TransportCableManager::compute_mass(double a_length, double a_rho)
{
  /*get_mass_damp(elem_length,elem_rho,damping)
        where
        elem_length is the initial length of each element
        elem_rho is the vector which e-th element is the linear density of element number e

  */
  int ndof = m_results.q0->size();
  if (not m_results.mass)
    m_results.mass = std::make_shared<SimpleMatrix>(ndof, ndof, Siconos::UBLAS_TYPE::SPARSE);
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

void TransportCableManager::compute_external_load(double a_length, double a_rho)
{
  size_t ndof = m_results.q0->size();
  if (not m_results.b)
    m_results.b = std::make_shared<SiconosVector>(ndof);

  double k = -9.81 * a_rho * a_length;
  for (size_t i = 2; i < ndof; i += 3) {
    m_results.b->setValue(i, k);
  }
}

#endif
