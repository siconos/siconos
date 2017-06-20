/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include "OSNSMatrix.hpp"
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
// #define DEBUG_NOCOLOR
#include "debug.h"

// --- numerics headers ---
#include "NonSmoothDrivers.h"
#include "LCP_Solvers.h"

using namespace RELATION;

LCP::LCP(int numericsSolverId):
  LinearOSNS(numericsSolverId)
{
  _numerics_problem.reset(new LinearComplementarityProblem);
  linearComplementarity_setDefaultSolverOptions(NULL, &*_numerics_solver_options, _numerics_solver_id);
}

int LCP::compute(double time)
{

  int info = 0;

  // --- Prepare data for LCP computing ---
  // And check if there is something to be done
  bool cont = preCompute(time);
  if (!cont)
    return info;

  // --- Call Numerics driver ---
  // Inputs:
  // - the problem (M,q ...)
  // - the unknowns (z,w)
  // - the options for the solver (name, max iteration number ...)
  // - the global options for Numerics (verbose mode ...)
  DEBUG_PRINTF("LCP : sizeOutput=%d\n", _sizeOutput);
  DEBUG_PRINTF("_indexSetLevel = %i\n", _indexSetLevel);
  DEBUG_EXPR(display(););

  if (_sizeOutput != 0)
  {

    // The LCP in Numerics format
    _numerics_problem->M = &*_M->getNumericsMatrix();
    _numerics_problem->q = _q->getArray();
    _numerics_problem->size = _sizeOutput;

    //const char * name = &*_numerics_solver_options->solverName;
    if (_numerics_solver_options->solverId == SICONOS_LCP_ENUM)
    {
      lcp_enum_init(&*_numerics_problem, &*_numerics_solver_options, 1);


    }
    info = linearComplementarity_driver(&*_numerics_problem, _z->getArray() , _w->getArray() ,
                                        &*_numerics_solver_options);

    if (_numerics_solver_options->solverId == SICONOS_LCP_ENUM)
    {
      lcp_enum_reset(&*_numerics_problem, &*_numerics_solver_options, 1);


    }

    // --- Recovering of the desired variables from LCP output ---
    postCompute();

    DEBUG_EXPR(display());

  }
  return info;
}

LCP::~LCP()
{
  solver_options_delete(&*_numerics_solver_options);
}
