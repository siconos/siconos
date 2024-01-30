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
#include "RollingFrictionContact.hpp"

#include "Interaction.hpp"
#include "NewtonImpactRollingFrictionNSL.hpp"
#include "NumericsSolversNamespace.h"  // solver_options stuff
#include "OSNSMatrix.hpp"
#include "SiconosVector.hpp"
#include "Simulation.hpp"
#include "Topology.hpp"

siconos::nonsmooth_formulations::RollingFrictionContact::RollingFrictionContact(int dimPb,
                                                                    int numericsSolverId)
    : RollingFrictionContact(dimPb,
                             std::shared_ptr<SolverOptions>(
                                 solver_options_create(numericsSolverId),
                                 solver_options_delete))
{
}

siconos::nonsmooth_formulations::RollingFrictionContact::RollingFrictionContact(
    int dimPb, std::shared_ptr<SolverOptions> options)
    : LinearOSNS(options), _contactProblemDim(dimPb)
{
  if (dimPb == 3 && options->solverId == SICONOS_ROLLING_FRICTION_3D_NSGS) {
    _numerics_solver_options.reset(solver_options_create(
                                       SICONOS_ROLLING_FRICTION_2D_NSGS),
                                   solver_options_delete);
  }

  if (dimPb == 5) {
    _rolling_frictionContact_driver = &rolling_fc3d_driver;
  }
  else if (dimPb == 3) {
    _rolling_frictionContact_driver = &rolling_fc2d_driver;
  }
  else
    THROW_EXCEPTION(
        "Wrong dimension value (only 5 (3D) or 3 (2D) are allowed for RollingFrictionContact "
        "constructor.");

  _mu = std::make_shared<std::vector<double>>();
  _muR = std::make_shared<std::vector<double>>();
}

void siconos::nonsmooth_formulations::RollingFrictionContact::initialize(
    std::shared_ptr<siconos::simulation::Simulation> sim)
{
  // - Checks memory allocation for main variables (M,q,w,z)
  // - Formalizes the problem if the topology is time-invariant

  // This function performs all steps that are time-invariant

  // General initialize for OneStepNSProblem
  LinearOSNS::initialize(sim);

  // Connect to the right function according to dim. of the problem

  // get topology
  auto topology = simulation()->nonSmoothDynamicalSystem()->topology();

  // Note that interactionBlocks is up to date since updateInteractionBlocks
  // has been called during OneStepNSProblem::initialize()

  // Fill vector of friction coefficients
  auto sizeMu = simulation()->nonSmoothDynamicalSystem()->topology()->indexSet(0)->size();
  _mu->reserve(sizeMu);
  _muR->reserve(sizeMu);

  // If the topology is TimeInvariant ie if M structure does not
  // change during simulation:

  if (topology->indexSet0()->size() > 0) {
    // Get index set from Simulation
    auto indexSet = simulation()->indexSet(indexSetLevel());
    siconos::graphs::InteractionsGraph::VIterator ui, uiend;
    for (std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui) {
      _mu->push_back(
          std::static_pointer_cast<siconos::modeling::NewtonImpactRollingFrictionNSL>(
              indexSet->bundle(*ui)->nonSmoothLaw())
              ->mu());
      _muR->push_back(
          std::static_pointer_cast<siconos::modeling::NewtonImpactRollingFrictionNSL>(
              indexSet->bundle(*ui)->nonSmoothLaw())
              ->muR());
    }
  }
}

void siconos::nonsmooth_formulations::RollingFrictionContact::updateMu()
{
  _mu->clear();
  _muR->clear();
  auto indexSet = simulation()->indexSet(indexSetLevel());
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  for (std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui) {
    _mu->push_back(std::static_pointer_cast<siconos::modeling::NewtonImpactRollingFrictionNSL>(
                       indexSet->bundle(*ui)->nonSmoothLaw())
                       ->mu());
    _muR->push_back(
        std::static_pointer_cast<siconos::modeling::NewtonImpactRollingFrictionNSL>(
            indexSet->bundle(*ui)->nonSmoothLaw())
            ->muR());
  }
}

std::shared_ptr<RollingFrictionContactProblem>
siconos::nonsmooth_formulations::RollingFrictionContact::frictionContactProblem()
{
  auto numerics_problem = std::make_shared<RollingFrictionContactProblem>();
  numerics_problem->dimension = _contactProblemDim;
  numerics_problem->numberOfContacts = _sizeOutput / _contactProblemDim;
  numerics_problem->M = &*_M->numericsMatrix();
  numerics_problem->q = &*_q->data();
  numerics_problem->mu = _mu->data();
  numerics_problem->mu_r = _muR->data();
  return numerics_problem;
}

// RollingFrictionContactProblem *
// siconos::nonsmooth_formulations::RollingFrictionContact::frictionContactProblemPtr()
// {
//   RollingFrictionContactProblem *numerics_problem = &_numerics_problem;
//   numerics_problem->dimension = _contactProblemDim;
//   numerics_problem->numberOfContacts = _sizeOutput / _contactProblemDim;
//   numerics_problem->M = &*_M->numericsMatrix();
//   numerics_problem->q = &*_q->getArray();
//   numerics_problem->mu = _mu->data();
//   numerics_problem->mu_r = _muR->data();
//   return numerics_problem;
// }

int siconos::nonsmooth_formulations::RollingFrictionContact::solve(
    std::shared_ptr<RollingFrictionContactProblem> problem)
{
  if (!problem) {
    problem = frictionContactProblem();
  }

  return (*_rolling_frictionContact_driver)(&*problem, &*_z->data(), &*_w->data(),
                                            &*_numerics_solver_options);
}

bool siconos::nonsmooth_formulations::RollingFrictionContact::checkCompatibleNSLaw(
    siconos::modeling::NonSmoothLaw &nslaw)
{
  float type_number =
      static_cast<float>(siconos::types::type_value(nslaw)) + 0.1 * nslaw.size();
  _nslawtype.insert(type_number);

  if (siconos::types::type_value(nslaw) !=
      siconos::modeling::Type::NewtonImpactRollingFrictionNSL) {
    THROW_EXCEPTION(
        "\nsiconos::nonsmooth_formulations::RollingFrictionContact::checkCompatibleNSLaw -  \n\
                      The chosen nonsmooth law is not compatible with FrictionalContact one step nonsmooth problem. \n\
                      Compatible siconos::modeling::NonSmoothLaw are: NewtonImpactRollingFrictionNSL (2D or 3D) \n");
    return false;
  }
  if (_nslawtype.size() > 1) {
    THROW_EXCEPTION(
        "\nsiconos::nonsmooth_formulations::RollingFrictionContact::checkCompatibleNSLaw -  \n\
                     Compatible siconos::modeling::NonSmoothLaw are: NewtonImpactRollingFrictionNSL (2D or 3D), but you cannot mix them \n");
    return false;
  }

  return true;
}

int siconos::nonsmooth_formulations::RollingFrictionContact::compute(double time)
{
  int info = 0;
  // --- Prepare data for RollingFrictionContact computing ---
  bool cont = preCompute(time);
  if (!cont) {
    return info;
  }
  // nothing to do
  if (indexSetLevel() == siconos::internal::LEVELMAX) {
    return info;
  }

  updateMu();

  // --- Call Numerics driver ---
  // Inputs:
  // - the problem (M,q ...)
  // - the unknowns (z,w)
  // - the options for the solver (name, max iteration number ...)
  // - the global options for Numerics (verbose mode ...)
  if (_sizeOutput != 0) {
    // Call Numerics Driver for RollingFrictionContact
    info = solve();
    postCompute();
  }

  return info;
}

void siconos::nonsmooth_formulations::RollingFrictionContact::display() const
{
  std::cout << "===== " << _contactProblemDim << "D Rolling Friction Contact Problem "
            << std::endl;
  std::cout << "of size " << _sizeOutput << "(ie " << _sizeOutput / _contactProblemDim
            << " contacts).\n";
  LinearOSNS::display();
}
