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
#include "GlobalFrictionContact.hpp"

#include <cstddef>

#include "Interaction.hpp"
#include "LagrangianDS.hpp"
#include "LagrangianSparseDS.hpp"
#include "MohrCoulombPlasticityNSL.hpp"
#include "MoreauJeanGOSI.hpp"  // Numerics Header
#include "NewtonEulerDS.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "NumericsSolversNamespace.h"  // solver_options stuff
#include "OSNSMatrix.hpp"
#include "OneStepIntegrator.hpp"
#include "SecondOrderDS.hpp"
#include "SiconosVector.hpp"
#include "Simulation.hpp"
// #include "Topology.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

// Constructor from solver id - Uses delegated constructor
siconos::nonsmooth_formulations::GlobalFrictionContact::GlobalFrictionContact(
    int dimPb, const int numericsSolverId)
    : GlobalFrictionContact(
          dimPb, std::shared_ptr<SolverOptions>(solver_options_create(numericsSolverId),
                                                solver_options_delete)) {}

// Constructor based on a pre-defined solver options set.
siconos::nonsmooth_formulations::GlobalFrictionContact::GlobalFrictionContact(
    int dimPb, std::shared_ptr<SolverOptions> options)
    : LinearOSNS(options, LinearOSNSAssemblyType::GLOBAL),
      _contactProblemDim(dimPb),
      _gfc_driver(&gfc3d_driver) {
  if (dimPb == 2) {
    _gfc_driver = &gfc2d_driver;
  } else if (dimPb == 3) {
    _gfc_driver = &gfc3d_driver;
  }

  // Reset default storage type for numerics matrices.
  _numericsMatrixStorageType = NM_SPARSE;
}

void siconos::nonsmooth_formulations::GlobalFrictionContact::initVectorsMemory() {
  // Memory allocation for reaction, and velocity
  LinearOSNS::initVectorsMemory();

  if (!_globalVelocities)
    _globalVelocities = std::make_shared<siconos::algebra::SiconosVector>(_maxSize);
  else {
    if (_globalVelocities->size() != _maxSize) _globalVelocities->resize(_maxSize);
  }

  if (!_b)
    _b = std::make_shared<siconos::algebra::SiconosVector>(_maxSize);
  else {
    if (_b->size() != _maxSize) _b->resize(_maxSize);
  }
}

void siconos::nonsmooth_formulations::GlobalFrictionContact::initialize(
    std::shared_ptr<siconos::simulation::Simulation> sim) {
  // - Checks memory allocation for main variables (M,q,w,z)
  // - Formalizes the problem if the topology is time-invariant

  // This function performs all steps that are time-invariant

  // General initialize for OneStepNSProblem
  OneStepNSProblem::initialize(sim);

  initVectorsMemory();

  // get topology
  auto topology = simulation()->nonSmoothDynamicalSystem()->topology();

  // Note that interactionBlocks is up to date since updateInteractionBlocks has been called
  // during OneStepNSProblem::initialize()

  // Fill vector of friction coefficients
  auto I0 = topology->indexSet0();
  _mu = std::make_shared<std::vector<double>>();
  _mu->reserve(I0->size());

  initOSNSMatrix();
}

std::shared_ptr<GlobalFrictionContactProblem>
siconos::nonsmooth_formulations::GlobalFrictionContact::globalFrictionContactProblem() {
  auto numerics_problem = std::make_shared<GlobalFrictionContactProblem>();
  numerics_problem->M = &*_W->numericsMatrix();
  numerics_problem->H = &*_H->numericsMatrix();
  numerics_problem->q = _q->data();
  numerics_problem->b = _b->data();
  numerics_problem->numberOfContacts = _sizeOutput / _contactProblemDim;
  numerics_problem->mu = _mu->data();
  numerics_problem->dimension = _contactProblemDim;
  if (_assemblyType == LinearOSNSAssemblyType::GLOBAL_REDUCED) {
    numerics_problem->M_inverse = &*_W_inverse->numericsMatrix();
  }
  return numerics_problem;
}

bool siconos::nonsmooth_formulations::GlobalFrictionContact::checkCompatibleNSLaw(
    siconos::modeling::NonSmoothLaw& nslaw) {
  float type_number =
      static_cast<float>(siconos::types::type_value(nslaw)) + 0.1 * nslaw.size();
  _nslawtype.insert(type_number);

  if ((siconos::types::type_value(nslaw) !=
       siconos::modeling::Type::NewtonImpactFrictionNSL) ||
      (siconos::types::type_value(nslaw) !=
       siconos::modeling::Type::MohrCoulombPlasticityNSL)) {
    THROW_EXCEPTION(
        "\nsiconos::nonsmooth_formulations::GlobalFrictionContact::checkCompatibleNSLaw -  \n\
                      The chosen nonsmooth law is not compatible with FrictionalContact one step nonsmooth problem. \n\
                      Compatible siconos::modeling::NonSmoothLaw is NewtonImpactFrictionNSL (2D or 3D) \n");
    return false;
  }
  if (_nslawtype.size() > 1) {
    THROW_EXCEPTION(
        "\nFrictionContact::checkCompatibleNSLaw -  \n\
                     Compatible siconos::modeling::NonSmoothLaw is : NewtonImpactFrictionNSL (2D or 3D), but you cannot mix them \n");
    return false;
  }

  return true;
}

// #define WITH_TIMER

void siconos::nonsmooth_formulations::GlobalFrictionContact::compute_q() {
  if (_q->size() != _sizeGlobalOutput) _q->resize(_sizeGlobalOutput);

  size_t offset = 0;
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  auto& DSG0 = *simulation()->nonSmoothDynamicalSystem()->dynamicalSystems();
  for (std::tie(dsi, dsend) = DSG0.vertices(); dsi != dsend; ++dsi) {
    auto ds = DSG0.bundle(*dsi);
    auto ds_size = ds->dimension();

    // OneStepIntegrator& Osi = *DSG0.properties(DSG0.descriptor(ds)).osi;
    // siconos::integrators::IntegratorType osiType = Osi.getType();
    auto mjgosi = std::dynamic_pointer_cast<siconos::integrators::MoreauJeanGOSI>(
        DSG0.properties(DSG0.descriptor(ds)).osi);
    if (mjgosi) {
      auto& ds_work_vectors = *DSG0.properties(DSG0.descriptor(ds)).workVectors;
      if ((std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) ||
          (std::dynamic_pointer_cast<siconos::modeling::LagrangianSparseDS>(ds)) ||
          (std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds))) {
        auto& vfree = *ds_work_vectors[siconos::integrators::MoreauJeanGOSI::FREE];
        _q->segment(offset, ds_size) = vfree.head(ds_size);
      }
    } else
      THROW_EXCEPTION(
          "siconos::mechanics::fem::nonsmooth_formulations::GlobalFrictionContact::"
          "computeq: only implemented for a MoreauJeanGOSI integrator.");

    offset += ds_size;
  }
}

bool siconos::nonsmooth_formulations::GlobalFrictionContact::preCompute(double time) {
  DEBUG_BEGIN(
      "siconos::nonsmooth_formulations::GlobalFrictionContact::preCompute(double time)\n");
  // This function is used to prepare data for the GlobalFrictionContact problem
  // - computation of M, H _tildeLocalVelocity and q
  // - set _sizeOutput, sizeLocalOutput

  // If the topology is time-invariant, only q needs to be computed at each time step.
  // M, _sizeOutput have been computed in initialize and are uptodate.

  // Get topology
  auto topology = simulation()->nonSmoothDynamicalSystem()->topology();
  DEBUG_PRINTF("indexSetLevel = %i\n", indexSetLevel());
  if (indexSetLevel() == siconos::internal::LEVELMAX) {
    DEBUG_END(
        "siconos::nonsmooth_formulations::GlobalFrictionContact::preCompute(double time)\n");
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
    //   DEBUG_END("siconos::nonsmooth_formulations::GlobalFrictionContact::preCompute(double
    //   time)\n"); return false; }

    _mu->clear();
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
#ifdef WITH_TIMER
    auto start = std::chrono::system_clock::now();
#endif
    // fill _W
    _W->fillW(DSG0);
    _sizeGlobalOutput = _W->rows();
    DEBUG_PRINTF("sizeW = %lu \n", _sizeGlobalOutput);
#ifdef WITH_TIMER
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "\nGlobalFrictionContact: fill W  " << elapsed << " ms" << std::endl;
#endif

    if (_assemblyType == LinearOSNSAssemblyType::GLOBAL_REDUCED) {
      // fill _W_inverse
      _W_inverse->fillWinverse(DSG0);
    }
#ifdef WITH_TIMER
    auto end_old = end;
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - end_old).count();
    std::cout << "GlobalFrictionContact: fillW inverse " << elapsed << " ms" << std::endl;
#endif

    // fill _q
    compute_q();
    DEBUG_EXPR(siconos::algebra::print(*_q););
#ifdef WITH_TIMER
    end_old = end;
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - end_old).count();
    std::cout << "GlobalFrictionContact: fill q " << elapsed << " ms" << std::endl;
#endif

    /************************************/

    // fill H
    _H->fillH(DSG0, indexSet);
    DEBUG_EXPR(NM_display(_H->numericsMatrix().get()););

    _sizeOutput = _H->cols();
    DEBUG_PRINTF("_sizeOutput = %i\n ", _sizeOutput);
#ifdef WITH_TIMER
    end_old = end;
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - end_old).count();
    std::cout << "GlobalFrictionContact: fill H " << elapsed << " ms" << std::endl;
#endif

    // fill _b
    if (_b->size() != _sizeOutput) _b->resize(_sizeOutput);

    siconos::graphs::InteractionsGraph::VIterator ui, uiend;
    for (std::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui) {
      auto inter = indexSet.bundle(*ui);

      // assert(siconos::types::type_value(*(inter->nonSmoothLaw())) ==
      //        siconos::modeling::Type::NewtonImpactFrictionNSL);
      if (siconos::types::type_value(*(inter->nonSmoothLaw())) ==
          siconos::modeling::Type::NewtonImpactFrictionNSL)
        _mu->push_back(std::static_pointer_cast<siconos::modeling::NewtonImpactFrictionNSL>(
                           inter->nonSmoothLaw())
                           ->mu());  // curious !!
      // if(siconos::types::type_value(*(inter->nonSmoothLaw())) ==
      //     siconos::modeling::Type::MohrCoulombPlasticityNSL)
      //     _mu->push_back(tan(std::static_pointer_cast<siconos::modeling::MohrCoulombPlasticityNSL>(
      //                      inter->nonSmoothLaw())
      //                      ->phi()));

      auto ds1 = indexSet.properties(*ui).source;
      auto ds2 = indexSet.properties(*ui).target;
      auto mjgosi1 = std::dynamic_pointer_cast<siconos::integrators::MoreauJeanGOSI>(
          DSG0.properties(DSG0.descriptor(ds1)).osi);

      if (mjgosi1 and std::dynamic_pointer_cast<siconos::integrators::MoreauJeanGOSI>(
                          DSG0.properties(DSG0.descriptor(ds2)).osi)) {
        // std::cout << "MoreauJeanGOSI case" << std::endl;
        mjgosi1->NonSmoothLawContributionToOutput(inter, *this);
      } else {
        THROW_EXCEPTION(
            "siconos::nonsmooth_formulations::GlobalFrictionContact::computeq. Not yet "
            "implemented for "
            "the given Integrator type ");
      }
      auto& osnsp_rhs = *(*indexSet.properties(*ui)
                               .workVectors)[siconos::integrators::MoreauJeanGOSI::OSNSP_RHS];
      auto pos = indexSet.properties(*ui).absolute_position;
      auto sizeY = inter->dimension();
      _b->segment(pos, sizeY) = osnsp_rhs.head(sizeY);
    }
    DEBUG_EXPR(siconos::algebra::print(*_b););
#ifdef WITH_TIMER
    end_old = end;
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - end_old).count();
    std::cout << "GlobalFrictionContact: fill b " << elapsed << " ms" << std::endl;
#endif

    // Checks z and w sizes and reset if necessary
    if (_z->size() != _sizeOutput) {
      _z->resize(_sizeOutput);
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
    // nothing to do (IsLinear and not changed)
#ifdef WITH_TIMER
    end_old = end;
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - end_old).count();
    std::cout << "GlobalFrictionContact: init w and z and V " << elapsed << " ms" << std::endl;
#endif
  }
  DEBUG_END(
      "siconos::nonsmooth_formulations::GlobalFrictionContact::preCompute(double time)\n");
  return true;
}

int siconos::nonsmooth_formulations::GlobalFrictionContact::compute(double time) {
  int info = 0;
  // --- Prepare data for GlobalFrictionContact computing ---
  bool cont = preCompute(time);
  if (!cont) return info;
  updateMu();

  // --- Call Numerics solver ---
  // if(_sizeGlobalOutput != 0)
  {
    info = solve();
    DEBUG_EXPR(display(););
    postCompute();
  }
  return info;
}

int siconos::nonsmooth_formulations::GlobalFrictionContact::solve() {
  //   std::shared_ptr<GlobalFrictionContactProblem> problem) {
  // if (!problem) {
  auto problem = globalFrictionContactProblem();

  return (*_gfc_driver)(&*problem, _z->data(), _w->data(), _globalVelocities->data(),
                        &*_numerics_solver_options);
}

void siconos::nonsmooth_formulations::GlobalFrictionContact::update_dynamicalsystems_state() {
  auto& DSG0 = *simulation()->nonSmoothDynamicalSystem()->dynamicalSystems();

  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = DSG0.vertices(); dsi != dsend; ++dsi) {
    auto ds = DSG0.bundle(*dsi);

    if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianDS>(ds)) {
      auto sizeDS = d->dimension();
      auto velocity = d->velocity();
      DEBUG_PRINTF("ds.number() : %i \n", ds.number());
      DEBUG_EXPR(siconos::algebra::print(*velocity););
      DEBUG_EXPR(siconos::algebra::print(*_globalVelocities););
      auto pos = DSG0.properties(*dsi).absolute_position;
      *velocity = _globalVelocities->segment(pos, sizeDS);
      DEBUG_EXPR(siconos::algebra::print(*velocity););
    } else if (auto d = std::dynamic_pointer_cast<siconos::modeling::LagrangianSparseDS>(ds)) {
      auto sizeDS = d->dimension();
      auto velocity = d->velocity();
      DEBUG_PRINTF("ds.number() : %i \n", ds.number());
      DEBUG_EXPR(siconos::algebra::print(*velocity););
      DEBUG_EXPR(siconos::algebra::print(*_globalVelocities););
      auto pos = DSG0.properties(*dsi).absolute_position;
      *velocity = _globalVelocities->segment(pos, sizeDS);
      DEBUG_EXPR(siconos::algebra::print(*velocity););
    } else if (auto neds = std::dynamic_pointer_cast<siconos::modeling::NewtonEulerDS>(ds)) {
      auto sizeDS = neds->dimension();
      auto twist = neds->twist();
      DEBUG_PRINTF("ds.number() : %i \n", ds.number());
      DEBUG_EXPR(siconos::algebra::print(*twist););
      DEBUG_EXPR(siconos::algebra::print(*_globalVelocities););
      auto pos = DSG0.properties(*dsi).absolute_position;
      *twist = _globalVelocities->segment(pos, sizeDS);
      DEBUG_EXPR(siconos::algebra::print(*twist););
    } else
      THROW_EXCEPTION(
          "siconos::nonsmooth_formulations::GlobalFrictionContact::postCompute() - Only "
          "implemented for "
          "Lagrangian or NewtonEuler systems.");
  }
}

void siconos::nonsmooth_formulations::GlobalFrictionContact::postCompute() {
  DEBUG_BEGIN(
      "siconos::nonsmooth_formulations::GlobalFrictionContact::postCompute(double time)\n");

  // This function is used to set y/lambda values using output from
  // primalfrictioncontact_driver Only Interactions (ie Interactions) of indexSet(leveMin) are
  // concerned.

  // === Get index set from Topology ===
  auto& indexSet =
      *simulation()->nonSmoothDynamicalSystem()->topology()->indexSet(_indexSetLevel);
  // y and lambda vectors
  //   // === Loop through "active" Interactions (ie present in indexSets[1]) ===

  size_t pos = 0;

  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  for (std::tie(ui, uiend) = indexSet.vertices(); ui != uiend;
       ++ui, pos += _contactProblemDim) {
    auto& inter = *indexSet.bundle(*ui);
    // Get Y and Lambda for the current Interaction
    auto y = inter.y(inputOutputLevel());
    auto lambda = inter.lambda(inputOutputLevel());
    // Copy _w/_z values, starting from index pos into y/lambda.

    lambda->segment(0, lambda->size()) = _z->segment(pos, lambda->size());
    DEBUG_EXPR(siconos::algebra::print(*lambda););
  }
  // global_velocities --> to dynamical systems
  update_dynamicalsystems_state();

  DEBUG_END(
      "siconos::nonsmooth_formulations::GlobalFrictionContact::postCompute(double time)\n");
}

void siconos::nonsmooth_formulations::GlobalFrictionContact::display() const {
  std::cout << "===== " << _contactProblemDim << "D Global Friction Contact Problem "
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

  std::cout << "============================================================\n";
}

void siconos::nonsmooth_formulations::GlobalFrictionContact::updateMu() {
  _mu->clear();
  auto indexSet = simulation()->indexSet(indexSetLevel());
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  for (std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui) {
    auto inter = indexSet->bundle(*ui);
    if (siconos::types::type_value(*(inter->nonSmoothLaw())) ==
        siconos::modeling::Type::NewtonImpactFrictionNSL)
      _mu->push_back(std::static_pointer_cast<siconos::modeling::NewtonImpactFrictionNSL>(
                         indexSet->bundle(*ui)->nonSmoothLaw())
                         ->mu());
    // if(siconos::types::type_value(*(inter->nonSmoothLaw())) ==
    //     siconos::modeling::Type::MohrCoulombPlasticityNSL)
    //   _mu->push_back(tan(std::static_pointer_cast<siconos::modeling::MohrCoulombPlasticityNSL>(
    //                          inter->nonSmoothLaw())
    //                          ->phi()));
  }
}
