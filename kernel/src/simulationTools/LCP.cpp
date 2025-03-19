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
#include "LCP.hpp"

#include "NonSmoothLaw.hpp"
#include "NumericsSolversNamespace.h"  // solver_options stuff
#include "OSNSMatrix.hpp"
#include "SiconosException.hpp"
#include "SiconosVector.hpp"
#include "TypeName.hpp"  // check nslaw type, should be replaced by dynamic_cast or variant ?
// #include "SolverOptions.h"
// #include "ComplementarityConditionNSL.hpp"
// // --- numerics headers ---
// #include "NonSmoothDrivers.h"
// #include "LCP_Solvers.h"

// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
// #define DEBUG_NOCOLOR
#include "siconos_debug.h"

siconos::nonsmooth_formulations::LCP::LCP(int numericsSolverId)
    : LCP(std::shared_ptr<SolverOptions>(solver_options_create(numericsSolverId),
                                         solver_options_delete)) {}

siconos::nonsmooth_formulations::LCP::LCP(std::shared_ptr<SolverOptions> options)
    : LinearOSNS(options),
      _numerics_problem(std::make_shared<LinearComplementarityProblem>()) {}

bool siconos::nonsmooth_formulations::LCP::checkCompatibleNSLaw(
    siconos::modeling::NonSmoothLaw& nslaw) {
  float type_number = (float)(siconos::types::type_value(nslaw));
  _nslawtype.insert(type_number);

  if (not(siconos::types::type_value(nslaw) ==
              siconos::modeling::Type::ComplementarityConditionNSL ||
          siconos::types::type_value(nslaw) == siconos::modeling::Type::NewtonImpactNSL ||
          siconos::types::type_value(nslaw) == siconos::modeling::Type::MultipleImpactNSL)) {
    THROW_EXCEPTION(
        "\nsiconos::nonsmooth_formulations::LCP::checkCompatibleNSLaw -  \n\
                      The chosen nonsmooth law is not compatible with LCP one step nonsmooth problem. \n \
                      Compatible NonSmoothLaw are: ComplementarityConditionNSL, MultipleImpactNSL or NewtonImpactNSL\n");
    return false;
  }

  return true;
}

int siconos::nonsmooth_formulations::LCP::solve() {
  // Note FP : wrap call to numerics solver inside this function
  // for python API (e.g. to allow profiling without C struct handling)

  // The LCP in Numerics format
  _numerics_problem->M = &*_M->numericsMatrix();
  _numerics_problem->q = _q->data();
  int info = 0;
  // const char * name = &*_numerics_solver_options->solverName;

  if (_numerics_solver_options->solverId == SICONOS_LCP_ENUM) {
    lcp_enum_init(&*_numerics_problem, &*_numerics_solver_options, 1);
  }

  // Call LCP Driver
  info = linearComplementarity_driver(&*_numerics_problem, _z->data(), _w->data(),
                                      &*_numerics_solver_options);

  if (_numerics_solver_options->solverId == SICONOS_LCP_ENUM) {
    lcp_enum_reset(&*_numerics_problem, &*_numerics_solver_options, 1);
  }
  return info;
}

int siconos::nonsmooth_formulations::LCP::compute(double time) {
  DEBUG_BEGIN("siconos::nonsmooth_formulations::LCP::compute(double time)\n");
  int info = 0;

  // --- Prepare data for LCP computing ---
  // And check if there is something to be done
  bool not_empty = preCompute(time);
  if (!not_empty) {
    DEBUG_PRINT("Nothing to compute\n");
    DEBUG_END("siconos::nonsmooth_formulations::LCP::compute(double time)\n");
    return info;
  }
  // --- Call Numerics driver ---
  // Inputs:
  // - the problem (M,q ...)
  // - the unknowns (z,w)
  // - the options for the solver (name, max iteration number ...)
  // - the global options for Numerics (verbose mode ...)
  DEBUG_PRINTF("LCP : sizeOutput=%d\n", _sizeOutput);
  DEBUG_PRINTF("_indexSetLevel = %i\n", _indexSetLevel);
  DEBUG_EXPR(display(););

  if (_sizeOutput != 0) {
    info = solve();
    // --- Recovering of the desired variables from LCP output ---
    postCompute();
    DEBUG_EXPR(display());
  }
  DEBUG_END("siconos::nonsmooth_formulations::LCP::compute(double time)\n");
  return info;
}
