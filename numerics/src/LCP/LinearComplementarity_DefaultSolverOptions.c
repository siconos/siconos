/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include "LinearComplementarityProblem.h"
#include "LCP_Solvers.h"
#include "lcp_cst.h"
#include "SolverOptions.h"
#include "NumericsMatrix.h"


#include "numerics_verbose.h"

int linearComplementarity_setDefaultSolverOptions(LinearComplementarityProblem* problem, SolverOptions* options, int solverId)
{

  int info = -1;
  switch (solverId)
  {
  case SICONOS_LCP_NSGS_SBM:
  {
    info =    linearComplementarity_nsgs_SBM_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_LCP_QP:
  {
    info =    linearComplementarity_qp_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_LCP_NSQP:
  {
    info =    linearComplementarity_nsqp_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_LCP_CPG:
  {
    info =    linearComplementarity_cpg_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_LCP_PGS:
  {
    info =    linearComplementarity_pgs_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_LCP_RPGS:
  {
    info =    linearComplementarity_rpgs_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_LCP_PSOR:
  {
    info =    linearComplementarity_psor_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_LCP_LATIN:
  {
    info =    linearComplementarity_latin_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_LCP_LATIN_W:
  {
    info =    linearComplementarity_latin_w_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_LCP_LEMKE:
  {
    info =    linearComplementarity_lexicolemke_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_LCP_ENUM:
  {
    info =    linearComplementarity_enum_setDefaultSolverOptions(problem, options);
    break;
  }
  case SICONOS_LCP_NEWTONMIN:
  {
    info =    linearComplementarity_newton_min_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_LCP_CONVEXQP_PG:
  {
    info =    linearComplementarity_ConvexQP_ProjectedGradient_setDefaultSolverOptions(options);
    break;
  }
  
  case SICONOS_LCP_PATH:
  case SICONOS_LCP_AVI_CAOFERRIS:
  case SICONOS_LCP_BARD:
  case SICONOS_LCP_MURTY:
  case SICONOS_LCP_PATHSEARCH:
  case SICONOS_LCP_NEWTON_MINFBLSA:
  case SICONOS_LCP_NEWTON_FBLSA:
  case SICONOS_LCP_GAMS:
  case SICONOS_LCP_PIVOT:
  case SICONOS_LCP_PIVOT_LUMOD:
  {
    solver_options_set(options, solverId);
    info = 0;
    break;
  }
  default:
  {
    numerics_error("linearComplementarity_setDefaultSolverOptions", "Unknown Solver");

  }
  }


  return info;
}
