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
#include "GenericMechanical.hpp"

#include "FrictionContactProblem.h"
#include "GenericMechanicalProblem.h"
#include "GenericMechanical_Solvers.h"
#include "Interaction.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "OSNSMatrix.hpp"
#include "RelayNSL.hpp"
#include "RelayProblem.h"
#include "Simulation.hpp"
#include "SolverOptions.h"
// #define DEBUG_BEGIN_END_ONLY
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
// #define DEBUG_WHERE_MESSAGES
#include "siconos_debug.h"

siconos::nonsmooth_formulations::GenericMechanical::GenericMechanical(int FC3D_Solver_Id)
    : GenericMechanical(std::shared_ptr<SolverOptions>(
          solver_options_create(SICONOS_GENERIC_MECHANICAL_NSGS), solver_options_delete)) {
  solver_options_update_internal(_numerics_solver_options.get(), 1, FC3D_Solver_Id);
}

siconos::nonsmooth_formulations::GenericMechanical::GenericMechanical(
    std::shared_ptr<SolverOptions> options)
    : LinearOSNS(options) {
  DEBUG_BEGIN(
      "siconos::nonsmooth_formulations::GenericMechanical::GenericMechanical(std::shared_ptr<"
      "siconos::"
      "numerics::SolverOptions> options)\n");
  // assert(options->solverId == SICONOS_GENERIC_MECHANICAL_NSGS); this will be checked in the
  // driver
  _numericsMatrixStorageType = NM_SPARSE_BLOCK;

  // auto deleter = [](GenericMechanicalProblem* pb) {
  //   genericMechanicalProblem_free(pb, GMP_FREE_GMP);
  // };

  //_pnumerics_GMP.reset(genericMechanicalProblem_new(), deleter);

  _pnumerics_GMP = genericMechanicalProblem_new();
  DEBUG_END(
      "siconos::nonsmooth_formulations::GenericMechanical::GenericMechanical(std::shared_ptr<"
      "siconos::"
      "numerics::SolverOptions> options)\n");
}

siconos::nonsmooth_formulations::GenericMechanical::~GenericMechanical() noexcept {
  if (_pnumerics_GMP) {
    _pnumerics_GMP->M = NULL;
    _pnumerics_GMP->q = NULL;
    // genericMechanicalProblem_free(_pnumerics_GMP, GMP_FREE_GMP);
  }
  _pnumerics_GMP = nullptr;
}

void siconos::nonsmooth_formulations::GenericMechanical::initialize(
    std::shared_ptr<siconos::simulation::Simulation> sim) {
  // - Checks memory allocation for main variables (M,q,w,z)
  // - Formalizes the problem if the topology is time-invariant

  // This function performs all steps that are time-invariant

  // General initialize for OneStepNSProblem
  LinearOSNS::initialize(sim);
}

bool siconos::nonsmooth_formulations::GenericMechanical::checkCompatibleNSLaw(
    siconos::modeling::NonSmoothLaw& nslaw) {
  // do nothoing since it is check in
  // siconos::nonsmooth_formulations::GenericMechanical::computeDiagonalInteractionBlock
  return true;
}

void siconos::nonsmooth_formulations::GenericMechanical::computeDiagonalInteractionBlock(
    const siconos::graphs::InteractionsGraph::VDescriptor& vd) {
  auto indexSet = simulation()->indexSet(indexSetLevel());
  // bool isTimeInvariant =
  // simulation()->nonSmoothDynamicalSystem()->topology()->isTimeInvariant();

  /*Build the corresponding numerics problems*/

  auto DS1 = indexSet->properties(vd).source;
  auto DS2 = indexSet->properties(vd).target;
  auto inter = indexSet->bundle(vd);

  DEBUG_PRINT(
      "siconos::nonsmooth_formulations::GenericMechanical::computeInteractionBlock: add "
      "problem of type ");

  if (!_hasBeenUpdated) {
    auto size = inter->nonSmoothLaw()->size();
    if (siconos::types::type_value(*(inter->nonSmoothLaw())) ==
        siconos::modeling::Type::EqualityConditionNSL) {
      gmp_add(_pnumerics_GMP, SICONOS_NUMERICS_PROBLEM_EQUALITY, size);
      DEBUG_PRINT("siconos::modeling::Type::EqualityConditionNSL\n");
      // pAux->size= inter->nonSmoothLaw()->size();
    } else if (siconos::types::type_value(*(inter->nonSmoothLaw())) ==
               siconos::modeling::Type::NewtonImpactNSL) {
      gmp_add(_pnumerics_GMP, SICONOS_NUMERICS_PROBLEM_LCP, size);
      DEBUG_PRINT(" siconos::modeling::Type::NewtonImpactNSL\n");
    } else if (siconos::types::type_value(*(inter->nonSmoothLaw())) ==
               siconos::modeling::Type::RelayNSL) {
      auto* pAux = static_cast<RelayProblem*>(
          gmp_add(_pnumerics_GMP, SICONOS_NUMERICS_PROBLEM_RELAY, size));
      auto nsLaw =
          std::static_pointer_cast<siconos::modeling::RelayNSL>(inter->nonSmoothLaw());
      for (decltype(size) i = 0; i < size; i++) {
        pAux->lb[i] = nsLaw->lb();
        pAux->ub[i] = nsLaw->ub();
      }
      DEBUG_PRINT(" siconos::modeling::Type::RelayNSL\n");
    } else if (siconos::types::type_value(*(inter->nonSmoothLaw())) ==
               siconos::modeling::Type::NewtonImpactFrictionNSL) {
      if (size == 3) {
        auto* pAux = static_cast<FrictionContactProblem*>(
            gmp_add(_pnumerics_GMP, SICONOS_NUMERICS_PROBLEM_FC3D, size));
        auto nsLaw = std::static_pointer_cast<siconos::modeling::NewtonImpactFrictionNSL>(
            inter->nonSmoothLaw());
        pAux->dimension = 3;
        pAux->numberOfContacts = 1;
        *(pAux->mu) = nsLaw->mu();

        DEBUG_PRINT(" siconos::modeling::Type::NewtonImpactFrictionNSL\n");
      } else if (size == 2) {
        auto* pAux = static_cast<FrictionContactProblem*>(
            gmp_add(_pnumerics_GMP, SICONOS_NUMERICS_PROBLEM_FC2D, size));
        auto nsLaw = std::static_pointer_cast<siconos::modeling::NewtonImpactFrictionNSL>(
            inter->nonSmoothLaw());
        pAux->dimension = 2;
        pAux->numberOfContacts = 1;
        *(pAux->mu) = nsLaw->mu();

        DEBUG_PRINT(" siconos::modeling::Type::NewtonImpactFrictionNSL\n");
      }
    } else {
      THROW_EXCEPTION(
          "siconos::nonsmooth_formulations::GenericMechanical::"
          "computeDiagonalInteractionBlock- not yet "
          "implemented for that NSLAW type");
    }
  }

  LinearOSNS::computeDiagonalInteractionBlock(vd);
}

void siconos::nonsmooth_formulations::GenericMechanical::computeInteractionBlock(
    const siconos::graphs::InteractionsGraph::EDescriptor& ed) {
  LinearOSNS::computeInteractionBlock(ed);
}

int siconos::nonsmooth_formulations::GenericMechanical::compute(double time) {
  DEBUG_BEGIN("siconos::nonsmooth_formulations::GenericMechanical::compute(double time)\n");
  int info = 0;
  // --- Prepare data for GenericMechanical computing ---
  bool cont = preCompute(time);

  if (!cont) return info;
  // MB: if _hasBeenUpdated is set true then :
  // LinearOSNS.cpp:602
  // position unitialized, pos get a wrong value then :
  // computeqBlock(inter, pos) -> SEGFAULT
  // so I comment this:
  // _hasBeenUpdated = true;
  /*
    La matrice _M est construite.  Ici, il faut construire les
    sous-problemes, c'est a dire completer les champs des
    NumericsProblem (_mu, _e, _en, les dimentions...).  Il faut aussi
    remplir la sous matrice M du sous-probleme.  Pour cela, on peut
    boucler sur les interactions et completer le membres
    _numerics_problem.problems[i] and
    _numerics_problem.problemsType[i].
   */

  //......

  // --- Call Numerics driver ---
  // Inputs:
  // - the problem (M,q ...)
  // - the unknowns (z,w)
  // - the options for the solver (name, max iteration number ...)
  // - the global options for Numerics (verbose mode ...)

  if (_sizeOutput != 0) {
    // The GenericMechanical Problem in Numerics format
    assert(_pnumerics_GMP);
    _pnumerics_GMP->M = &*_M->numericsMatrix();
    _pnumerics_GMP->q = &*_q->data();
    DEBUG_EXPR(display(););
    // Call Numerics Driver for GenericMechanical
    //    display();
    info = gmp_driver(_pnumerics_GMP, &*_z->data(), &*_w->data(), &*_numerics_solver_options);
    // printf("siconos::nonsmooth_formulations::GenericMechanical::compute : R:\n");
    // siconos::algebra::print(*_z);
    _pnumerics_GMP->M = NULL;
    _pnumerics_GMP->q = NULL;

    postCompute();
  } else {
    DEBUG_PRINT(
        "siconos::nonsmooth_formulations::GenericMechanical::compute : sizeoutput is null\n");
  }
  DEBUG_END("siconos::nonsmooth_formulations::GenericMechanical::compute(double time)\n");
  return info;
}

void siconos::nonsmooth_formulations::GenericMechanical::display() const {
  std::cout << "===== " << "Generic mechanical Problem " << std::endl;
  LinearOSNS::display();
}

void siconos::nonsmooth_formulations::GenericMechanical::updateInteractionBlocks() {
  if (!_hasBeenUpdated) {
    _pnumerics_GMP = genericMechanicalProblem_new();
  }
  LinearOSNS::updateInteractionBlocks();
}
