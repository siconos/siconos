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
#include "GlobalRollingFrictionContact.hpp"

#include "Interaction.hpp"
#include "LagrangianDS.hpp"
#include "LagrangianSparseDS.hpp"
#include "MoreauJeanGOSI.hpp"  // Numerics Header
#include "NewtonEulerDS.hpp"
#include "NewtonImpactRollingFrictionNSL.hpp"
#include "NumericsSolversNamespace.h"  // solver_options stuff
#include "OSNSMatrix.hpp"
// #include "SiconosVector.hpp"
#include "Simulation.hpp"
// #include "Topology.hpp"
//  #define DEBUG_NOCOLOR
//  #define DEBUG_STDOUT
//  #define DEBUG_MESSAGES
#include "siconos_debug.h"

// Constructor from solver id - Uses delegated constructor
siconos::nonsmooth_formulations::GlobalRollingFrictionContact::GlobalRollingFrictionContact(
    int dimPb, const int numericsSolverId)
    : GlobalRollingFrictionContact(
          dimPb, std::shared_ptr<SolverOptions>(solver_options_create(numericsSolverId),
                                                solver_options_delete)) {}

// Constructor based on a pre-defined solver options set.
siconos::nonsmooth_formulations::GlobalRollingFrictionContact::GlobalRollingFrictionContact(
    int dimPb, std::shared_ptr<SolverOptions> options)
    : GlobalFrictionContact(dimPb, options), _g_rolling_driver(&g_rolling_fc3d_driver) {
  // Only rolling fc3d for the moment.
  if (_contactProblemDim != 5)
    THROW_EXCEPTION("GlobalRollingFrictionContact No solver for 2 dimensional problems");

  // Reset default storage type for numerics matrices.
  _numericsMatrixStorageType = NM_SPARSE;
}

void siconos::nonsmooth_formulations::GlobalRollingFrictionContact::initialize(
    std::shared_ptr<siconos::simulation::Simulation> sim) {
  GlobalFrictionContact::initialize(sim);
  // get topology
  auto topology = simulation()->nonSmoothDynamicalSystem()->topology();

  // Fill vector of rolling friction coefficients
  auto I0 = topology->indexSet0();
  _mu_r = std::make_shared<std::vector<double>>();
  _mu_r->reserve(I0->size());
}

std::shared_ptr<GlobalRollingFrictionContactProblem> siconos::nonsmooth_formulations::
    GlobalRollingFrictionContact::globalRollingFrictionContactProblem() {
  auto numerics_problem = std::make_shared<GlobalRollingFrictionContactProblem>();
  numerics_problem->M = &*_W->numericsMatrix();
  numerics_problem->H = &*_H->numericsMatrix();
  numerics_problem->q = _q->data();
  numerics_problem->b = _b->data();
  numerics_problem->numberOfContacts = _sizeOutput / _contactProblemDim;
  numerics_problem->mu = _mu->data();
  numerics_problem->mu_r = _mu_r->data();
  numerics_problem->dimension = 5;
  return numerics_problem;
}

bool siconos::nonsmooth_formulations::GlobalRollingFrictionContact::checkCompatibleNSLaw(
    siconos::modeling::NonSmoothLaw& nslaw) {
  if (siconos::types::type_value(nslaw) !=
      siconos::modeling::Type::NewtonImpactRollingFrictionNSL) {
    THROW_EXCEPTION(
        "\nsiconos::nonsmooth_formulations::GlobalRollingFrictionContact::checkCompatibleNSLaw -  \n\
                      The chosen nonsmooth law is not compatible with GlobalRollingFrictionalContact one step nonsmooth problem. \n\
                      Compatible siconos::modeling::NonSmoothLaw is NewtonImpactRollingFrictionNSL (3D) \n");
    return false;
  }
  NSLTypeKey nsltype{siconos::types::type_value(nslaw), nslaw.size()};
  _nslawtype.insert(nsltype);

  if (_nslawtype.size() > 1) {
    THROW_EXCEPTION(
        "\nFrictionContact::checkCompatibleNSLaw -  \n\
                     Compatible siconos::modeling::NonSmoothLaw is : NewtonImpactRollingFrictionNSL (3D), but you cannot mix them \n");
    return false;
  }

  return true;
}

bool siconos::nonsmooth_formulations::GlobalRollingFrictionContact::preCompute(double time) {
  DEBUG_BEGIN(
      "siconos::nonsmooth_formulations::GlobalRollingFrictionContact::preCompute(double "
      "time)\n");
  // This function is used to prepare data for the GlobalRollingFrictionContact problem
  // - computation of M, H _tildeLocalVelocity and q
  // - set _sizeOutput, sizeLocalOutput

  // If the topology is time-invariant, only q needs to be computed at each time step.
  // M, _sizeOutput have been computed in initialize and are uptodate.

  // Get topology
  auto topology = simulation()->nonSmoothDynamicalSystem()->topology();
  DEBUG_PRINTF("indexSetLevel = %i\n", indexSetLevel());
  if (indexSetLevel() == siconos::internal::LEVELMAX) {
    DEBUG_END(
        "siconos::nonsmooth_formulations::GlobalRollingFrictionContact::preCompute(double "
        "time)\n");
    return false;
  }
  if (!_hasBeenUpdated) {
    auto& indexSet =
        *simulation()->nonSmoothDynamicalSystem()->topology()->indexSet(_indexSetLevel);
    auto& DSG0 = *simulation()->nonSmoothDynamicalSystem()->dynamicalSystems();

    // compute size and nnz of M and collect all matrices
    // compute nnz of H and collect H blocks
    // fill  mu

    // if (_sizeOutput == 0)
    // {
    //   DEBUG_END("siconos::nonsmooth_formulations::GlobalRollingFrictionContact::preCompute(double
    //   time)\n"); return false; }

    _mu->clear();
    _mu_r->clear();
    //    _mu.reserve(indexSet.size())

#if !defined(SICONOS_USE_MAP_FOR_HASH)
    using dsMatMap = std::unordered_map<std::shared_ptr<siconos::modeling::DynamicalSystem>,
                                        siconos::algebra::SiconosMatrix*>;
    using dsPosMap =
        std::unordered_map<std::shared_ptr<siconos::modeling::DynamicalSystem>, std::size_t>;

#else
    using dsMatMap = std::map<std::shared_ptr<siconos::modeling::DynamicalSystem>,
                              siconos::algebra::SiconosMatrix*>;
    using dsPosMap =
        std::map<std::shared_ptr<siconos::modeling::DynamicalSystem>, std::size_t>;

#endif
    dsMatMap dsMat;
    dsPosMap absPosDS;

    // fill _W
    _W->fillW(DSG0);
    _sizeGlobalOutput = _W->rows();
    DEBUG_PRINTF("sizeW = %lu \n", _sizeGlobalOutput);

    // fill _q
    if (_q->size() != _sizeGlobalOutput) _q->resize(_sizeGlobalOutput);

    siconos::algebra::Index offset = 0;
    siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
    for (std::tie(dsi, dsend) = DSG0.vertices(); dsi != dsend; ++dsi) {
      auto ds = DSG0.bundle(*dsi);
      auto dss = ds->dimension();
      DEBUG_PRINTF("offset = %lu \n", offset);

      auto Osi = DSG0.properties(DSG0.descriptor(ds)).osi;

      if (std::dynamic_pointer_cast<siconos::integrators::MoreauJeanGOSI>(Osi)) {
        auto& ds_work_vectors = *DSG0.properties(DSG0.descriptor(ds)).workVectors;

        if (std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds) ||
            std::dynamic_pointer_cast<siconos::modeling::LagrangianSparseDS>(ds)) {
          auto& vfree = *ds_work_vectors[tools::enum_to_index(
              siconos::integrators::MoreauJeanGOSI::wk_ds::vfree)];
          _q->segment(offset, dss) = vfree;
        } else if (std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
          auto& vfree = *ds_work_vectors[tools::enum_to_index(
              siconos::integrators::MoreauJeanGOSI::wk_ds::vfree)];
          _q->segment(offset, dss) = vfree;
        }
      } else {
        THROW_EXCEPTION(
            "siconos::nonsmooth_formulations::GlobalRollingFrictionContact::computeq. Only  "
            "implemented "
            "for MoreauJeanGOSI integrator.");
      }
      offset += dss;
    }
    DEBUG_EXPR(siconos::algebra::print(*_q););

    /************************************/

    // fill H
    _H->fillH(DSG0, indexSet);
    DEBUG_EXPR(NM_display(_H->numericsMatrix().get()););

    _sizeOutput = _H->cols();
    DEBUG_PRINTF("_sizeOutput = %i\n ", _sizeOutput);

    // fill _b
    if (_b->size() != _sizeOutput) _b->resize(_sizeOutput);

    siconos::graphs::InteractionsGraph::VIterator ui, uiend;
    for (std::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui) {
      auto inter = indexSet.bundle(*ui);

      assert(siconos::types::type_value(*(inter->nonSmoothLaw())) ==
             siconos::modeling::Type::NewtonImpactRollingFrictionNSL);
      _mu->push_back(
          std::static_pointer_cast<siconos::modeling::NewtonImpactRollingFrictionNSL>(
              inter->nonSmoothLaw())
              ->mu());
      _mu_r->push_back(
          std::static_pointer_cast<siconos::modeling::NewtonImpactRollingFrictionNSL>(
              inter->nonSmoothLaw())
              ->muR());

      auto ds1 = indexSet.properties(*ui).source;
      auto ds2 = indexSet.properties(*ui).target;
      auto mjgosi1 = std::dynamic_pointer_cast<siconos::integrators::MoreauJeanGOSI>(
          DSG0.properties(DSG0.descriptor(ds1)).osi);

      if (mjgosi1 and std::dynamic_pointer_cast<siconos::integrators::MoreauJeanGOSI>(
                          DSG0.properties(DSG0.descriptor(ds2)).osi)) {
        mjgosi1->NonSmoothLawContributionToOutput(inter, *this);
      } else {
        THROW_EXCEPTION(
            "siconos::nonsmooth_formulations::GlobalRollingFrictionContact::computeq. Only  "
            "implemented "
            "for MoreauJeanGOSI integrator.");
      }
      auto& osnsp_rhs = *(*indexSet.properties(*ui).workVectors)[tools::enum_to_index(
          siconos::integrators::MoreauJeanGOSI::wk_inter::osnsp_rhs)];
      auto pos = indexSet.properties(*ui).absolute_position;
      auto sizeY = inter->dimension();
      _b->segment(pos, sizeY) = osnsp_rhs;
    }
    DEBUG_EXPR(siconos::algebra::print(*_b););
    // Checks z and w sizes and reset if necessary
    if (_z->size() != _sizeOutput) {
      _z->resize(_sizeOutput, Eigen::NoChange);
      _z->setZero();
    }

    if (_w->size() != _sizeOutput) {
      _w->resize(_sizeOutput);
      _w->setZero();
    }

    if (_globalVelocities->size() != _sizeGlobalOutput) {
      _globalVelocities->resize(_sizeGlobalOutput);
      _globalVelocities->setZero();
    }
  }
  DEBUG_END(
      "siconos::nonsmooth_formulations::GlobalRollingFrictionContact::preCompute(double "
      "time)\n");
  return true;
}

int siconos::nonsmooth_formulations::GlobalRollingFrictionContact::compute(double time) {
  DEBUG_BEGIN(
      "siconos::nonsmooth_formulations::GlobalRollingFrictionContact::compute(double time)\n")
  int info = 0;
  // --- Prepare data for GlobalRollingFrictionContact computing ---
  bool cont = preCompute(time);
  if (!cont) return info;
  updateMu();
  updateMur();
  // --- Call Numerics solver ---
  info = solve();
  DEBUG_EXPR(display(););
  postCompute();
  DEBUG_END(
      "siconos::nonsmooth_formulations::GlobalRollingFrictionContact::compute(double time)\n")
  return info;
}

int siconos::nonsmooth_formulations::GlobalRollingFrictionContact::solve(
    std::shared_ptr<GlobalRollingFrictionContactProblem> problem) {
  if (!problem) {
    problem = globalRollingFrictionContactProblem();
  }
  return (*_g_rolling_driver)(&*problem, _z->data(), _w->data(), _globalVelocities->data(),
                              &*_numerics_solver_options);
}

void siconos::nonsmooth_formulations::GlobalRollingFrictionContact::updateMur() {
  _mu_r->clear();
  auto indexSet = simulation()->indexSet(indexSetLevel());
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  for (std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui) {
    _mu_r->push_back(
        std::static_pointer_cast<siconos::modeling::NewtonImpactRollingFrictionNSL>(
            indexSet->bundle(*ui)->nonSmoothLaw())
            ->muR());
  }
}
void siconos::nonsmooth_formulations::GlobalRollingFrictionContact::updateMu() {
  _mu_r->clear();
  auto indexSet = simulation()->indexSet(indexSetLevel());
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  for (std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui) {
    _mu_r->push_back(
        std::static_pointer_cast<siconos::modeling::NewtonImpactRollingFrictionNSL>(
            indexSet->bundle(*ui)->nonSmoothLaw())
            ->mu());
  }
}

void siconos::nonsmooth_formulations::GlobalRollingFrictionContact::display() const {
  std::cout << "===== " << _contactProblemDim << "D Global Rolling Friction Contact Problem "
            << std::endl;
  std::cout << "size (_sizeOutput) " << _sizeOutput << "(ie "
            << _sizeOutput / _contactProblemDim << " contacts).\n";
  std::cout << "and  size (_sizeGlobalOutput) " << _sizeGlobalOutput << std::endl;
  std::cout << "_numericsMatrixStorageType" << _numericsMatrixStorageType << std::endl;
  std::cout << " - Matrix M  : \n";
  // if (_W) siconos::algebra::print(*_W);
  // else std::cout << "-> nullptr" <<std::endl;
  auto* W_NM = _W->numericsMatrix().get();
  if (W_NM) {
    NM_display(W_NM);
  }
  std::cout << " - Matrix H : \n";
  // if (_H) siconos::algebra::print(*_H);
  // else std::cout << "-> nullptr" <<std::endl;
  auto* H_NM = _H->numericsMatrix().get();
  if (H_NM) {
    NM_display(H_NM);
  }

  std::cout << " - Vector q : \n";
  if (_q)
    siconos::algebra::print(*_q);
  else
    std::cout << "-> nullptr\n";
  std::cout << " - Vector b : \n";
  if (_b)
    siconos::algebra::print(*_b);
  else
    std::cout << "-> nullptr\n";

  std::cout << " - Vector z (reaction) : \n";
  if (_z)
    siconos::algebra::print(*_z);
  else
    std::cout << "-> nullptr\n";

  std::cout << " - Vector w (local velocities): \n";
  if (_w)
    siconos::algebra::print(*_w);
  else
    std::cout << "-> nullptr\n";

  std::cout << " - Vector globalVelocities : \n";
  if (_globalVelocities)
    siconos::algebra::print(*_globalVelocities);
  else
    std::cout << "-> nullptr\n";

  std::cout << " - Vector mu : \n";
  if (_mu) {
    for (auto coeff : *_mu) std::cout << coeff << " ";
    std::cout << "\n";
  } else
    std::cout << "-> nullptr\n";

  std::cout << " - Vector mu_r : \n";
  if (_mu_r) {
    for (auto coeff : *_mu_r) std::cout << coeff << " ";

    std::cout << std::endl;
  } else
    std::cout << "-> nullptr\n";
  std::cout << "============================================================\n";
}
