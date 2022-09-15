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
#include "GenericMechanical.hpp"

#include "Interaction.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "NumericsSolversNamespace.h"
#include "OSNSMatrix.hpp"
#include "RelayNSL.hpp"
#include "SiconosVector.hpp"
#include "Simulation.hpp"
// #define DEBUG_BEGIN_END_ONLY
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
// #define DEBUG_WHERE_MESSAGES
#include "siconos_debug.h"

siconos::simulation::GenericMechanical::GenericMechanical(int FC3D_Solver_Id)
    : GenericMechanical(std::shared_ptr<siconos::numerics::SolverOptions>(
          siconos::numerics::solver_options_create(
              siconos::numerics::SICONOS_GENERIC_MECHANICAL_NSGS),
          siconos::numerics::solver_options_delete))
{
  siconos::numerics::solver_options_update_internal(_numerics_solver_options.get(), 1,
                                                    FC3D_Solver_Id);
}

siconos::simulation::GenericMechanical::GenericMechanical(
    std::shared_ptr<siconos::numerics::SolverOptions> options)
    : LinearOSNS(options)
{
  DEBUG_BEGIN(
      "siconos::simulation::GenericMechanical::GenericMechanical(std::shared_ptr<siconos::"
      "numerics::SolverOptions> options)\n");
  // assert(options->solverId == SICONOS_GENERIC_MECHANICAL_NSGS); this will be checked in the
  // driver
  _numericsMatrixStorageType = siconos::numerics::NM_SPARSE_BLOCK;
  _pnumerics_GMP = siconos::numerics::genericMechanicalProblem_new();
  DEBUG_END(
      "siconos::simulation::GenericMechanical::GenericMechanical(std::shared_ptr<siconos::"
      "numerics::SolverOptions> options)\n");
}

siconos::simulation::GenericMechanical::~GenericMechanical()
{
  genericMechanicalProblem_free(_pnumerics_GMP, siconos::numerics::GMP_FREE_GMP);
  _pnumerics_GMP = nullptr;
}

void siconos::simulation::GenericMechanical::initialize(
    std::shared_ptr<siconos::simulation::Simulation> sim)
{
  // - Checks memory allocation for main variables (M,q,w,z)
  // - Formalizes the problem if the topology is time-invariant

  // This function performs all steps that are time-invariant

  // General initialize for OneStepNSProblem
  LinearOSNS::initialize(sim);
}

bool siconos::simulation::GenericMechanical::checkCompatibleNSLaw(
    siconos::modeling::NonSmoothLaw& nslaw)
{
  // do nothoing since it is check in
  // siconos::simulation::GenericMechanical::computeDiagonalInteractionBlock
  return true;
}

void siconos::simulation::GenericMechanical::computeDiagonalInteractionBlock(
    const siconos::graphs::InteractionsGraph::VDescriptor& vd)
{
  auto indexSet = simulation()->indexSet(indexSetLevel());
  // bool isTimeInvariant =
  // simulation()->nonSmoothDynamicalSystem()->topology()->isTimeInvariant();

  /*Build the corresponding numerics problems*/

  auto DS1 = indexSet->properties(vd).source;
  auto DS2 = indexSet->properties(vd).target;
  auto inter = indexSet->bundle(vd);

  DEBUG_PRINT(
      "siconos::simulation::GenericMechanical::computeInteractionBlock: add problem of type ");

  if (!_hasBeenUpdated) {
    auto size = inter->nonSmoothLaw()->size();
    if (siconos::types::type_value(*(inter->nonSmoothLaw())) ==
        siconos::modeling::Type::EqualityConditionNSL) {
      siconos::numerics::gmp_add(_pnumerics_GMP,
                                 siconos::numerics::SICONOS_NUMERICS_PROBLEM_EQUALITY, size);
      DEBUG_PRINT("siconos::modeling::Type::EqualityConditionNSL\n");
      // pAux->size= inter->nonSmoothLaw()->size();
    }
    else if (siconos::types::type_value(*(inter->nonSmoothLaw())) ==
             siconos::modeling::Type::NewtonImpactNSL) {
      siconos::numerics::gmp_add(_pnumerics_GMP,
                                 siconos::numerics::SICONOS_NUMERICS_PROBLEM_LCP, size);
      DEBUG_PRINT(" siconos::modeling::Type::NewtonImpactNSL\n");
    }
    else if (siconos::types::type_value(*(inter->nonSmoothLaw())) ==
             siconos::modeling::Type::RelayNSL) {
      auto* pAux = static_cast<siconos::numerics::RelayProblem*>(siconos::numerics::gmp_add(
          _pnumerics_GMP, siconos::numerics::SICONOS_NUMERICS_PROBLEM_RELAY, size));
      auto nsLaw =
          std::static_pointer_cast<siconos::modeling::RelayNSL>(inter->nonSmoothLaw());
      for (int i = 0; i < size; i++) {
        pAux->lb[i] = nsLaw->lb();
        pAux->ub[i] = nsLaw->ub();
      }
      DEBUG_PRINT(" siconos::modeling::Type::RelayNSL\n");
    }
    else if (siconos::types::type_value(*(inter->nonSmoothLaw())) ==
             siconos::modeling::Type::NewtonImpactFrictionNSL) {
      if (size == 3) {
        auto* pAux =
            static_cast<siconos::numerics::FrictionContactProblem*>(siconos::numerics::gmp_add(
                _pnumerics_GMP, siconos::numerics::SICONOS_NUMERICS_PROBLEM_FC3D, size));
        auto nsLaw = std::static_pointer_cast<siconos::modeling::NewtonImpactFrictionNSL>(
            inter->nonSmoothLaw());
        pAux->dimension = 3;
        pAux->numberOfContacts = 1;
        *(pAux->mu) = nsLaw->mu();

        DEBUG_PRINT(" siconos::modeling::Type::NewtonImpactFrictionNSL\n");
      }
      else if (size == 2) {
        auto* pAux =
            static_cast<siconos::numerics::FrictionContactProblem*>(siconos::numerics::gmp_add(
                _pnumerics_GMP, siconos::numerics::SICONOS_NUMERICS_PROBLEM_FC2D, size));
        auto nsLaw = std::static_pointer_cast<siconos::modeling::NewtonImpactFrictionNSL>(
            inter->nonSmoothLaw());
        pAux->dimension = 2;
        pAux->numberOfContacts = 1;
        *(pAux->mu) = nsLaw->mu();

        DEBUG_PRINT(" siconos::modeling::Type::NewtonImpactFrictionNSL\n");
      }
    }
    else {
      THROW_EXCEPTION(
          "siconos::simulation::GenericMechanical::computeDiagonalInteractionBlock- not yet "
          "implemented for that NSLAW type");
    }
  }

  LinearOSNS::computeDiagonalInteractionBlock(vd);
}

void siconos::simulation::GenericMechanical::computeInteractionBlock(
    const siconos::graphs::InteractionsGraph::EDescriptor& ed)
{
  LinearOSNS::computeInteractionBlock(ed);
}

int siconos::simulation::GenericMechanical::compute(double time)
{
  DEBUG_BEGIN("siconos::simulation::GenericMechanical::compute(double time)\n");
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

    _pnumerics_GMP->M = &*_M->numericsMatrix();
    _pnumerics_GMP->q = &*_q->getArray();
    DEBUG_EXPR(display(););
    // Call Numerics Driver for GenericMechanical
    //    display();
    info = siconos::numerics::gmp_driver(_pnumerics_GMP, &*_z->getArray(), &*_w->getArray(),
                                         &*_numerics_solver_options);
    // printf("siconos::simulation::GenericMechanical::compute : R:\n");
    //_z->display();
    postCompute();
  }
  else {
    DEBUG_PRINT("siconos::simulation::GenericMechanical::compute : sizeoutput is null\n");
  }
  DEBUG_END("siconos::simulation::GenericMechanical::compute(double time)\n");
  return info;
}

void siconos::simulation::GenericMechanical::display() const
{
  std::cout << "===== "
            << "Generic mechanical Problem " << std::endl;
  LinearOSNS::display();
}

void siconos::simulation::GenericMechanical::updateInteractionBlocks()
{
  if (!_hasBeenUpdated) {
    //    printf("siconos::simulation::GenericMechanical::updateInteractionBlocks : must be
    //    updated\n");
    genericMechanicalProblem_free(_pnumerics_GMP, siconos::numerics::GMP_FREE_GMP);
    _pnumerics_GMP = siconos::numerics::genericMechanicalProblem_new();
  }
  LinearOSNS::updateInteractionBlocks();
}
