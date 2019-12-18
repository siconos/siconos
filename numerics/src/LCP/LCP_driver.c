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
#include <assert.h>                        // for assert
#include <stdio.h>                         // for printf, fprintf, stderr
#include <stdlib.h>                        // for exit, EXIT_FAILURE
#include "LCP_Solvers.h"                   // for lcp_pivot, lcp_compute_error
#include "LinearComplementarityProblem.h"  // for LinearComplementarityProblem
#include "NonSmoothDrivers.h"              // for linearComplementarity_driver
#include "NumericsFwd.h"                   // for SolverOptions, LinearCompl...
#include "NumericsMatrix.h"                // for NumericsMatrix
#include "SolverOptions.h"                 // for SolverOptions, solver_opti...
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"                         // for DEBUG_END, DEBUG_BEGIN
#include "lcp_cst.h"                       // for SICONOS_LCP_IPARAM_PIVOTIN...
#include "numerics_verbose.h"              // for numerics_error, verbose

const char* const   SICONOS_LCP_LEMKE_STR = "Lemke";
const char* const   SICONOS_LCP_NSGS_SBM_STR = "NSGS_SBM";
const char* const   SICONOS_LCP_PGS_STR = "PGS";
const char* const   SICONOS_LCP_CPG_STR = "CPG";
const char* const   SICONOS_LCP_LATIN_STR = "Latin";
const char* const   SICONOS_LCP_LATIN_W_STR = "Latin_w";
const char* const   SICONOS_LCP_QP_STR = "QP";
const char* const   SICONOS_LCP_NSQP_STR = "NSQP";
const char* const   SICONOS_LCP_NEWTONMIN_STR = "NewtonMin";
const char* const   SICONOS_LCP_NEWTON_FB_FBLSA_STR = "NewtonFB";
const char* const   SICONOS_LCP_NEWTON_MIN_FBLSA_STR = "NewtonMinFB";
const char* const   SICONOS_LCP_PSOR_STR = "PSOR";
const char* const   SICONOS_LCP_RPGS_STR = "RPGS";
const char* const   SICONOS_LCP_PATH_STR = "PATH";
const char* const   SICONOS_LCP_ENUM_STR = "ENUM";
const char* const   SICONOS_LCP_AVI_CAOFERRIS_STR = "AVI CaoFerris";
const char* const   SICONOS_LCP_PIVOT_STR = "Pivot based method";
const char* const   SICONOS_LCP_BARD_STR = "Bard-type pivoting method";
const char* const   SICONOS_LCP_MURTY_STR = "Murty's least index pivoting method";
const char* const   SICONOS_LCP_PATHSEARCH_STR = "For testing only: solver used in the Pathsearch algorithm";
const char* const   SICONOS_LCP_PIVOT_LUMOD_STR = "Pivot based method with BLU updates using LUMOD";
const char* const   SICONOS_LCP_GAMS_STR = "Using GAMS solvers";
const char* const   SICONOS_LCP_CONVEXQP_PG_STR = "Convex QP Projected Gradient";

static int lcp_driver_SparseBlockMatrix(LinearComplementarityProblem* problem, double *z, double *w, SolverOptions* options);

int lcp_driver_SparseBlockMatrix(LinearComplementarityProblem* problem, double *z, double *w, SolverOptions* options)
{
  DEBUG_BEGIN("lcp_driver_SparseBlockMatrix(...)\n");
  /* Checks storage type for the matrix M of the LCP */
  if(problem->M->storageType == 0)
    numerics_error("lcp_driver_SparseBlockMatrix", "forbidden type of storage for the matrix M of the LCP");

  /*
    The options for the global "block" solver are defined in options->\n
    options[i], for 0<i<numberOfSolvers-1 correspond to local solvers.
  */

  /* Output info. : 0: ok -  >0: problem (depends on solver) */
  int info = -1;

  /******************************************
   *  1 - Check for trivial solution
   ******************************************/
  int i = 0;
  int n = problem->size;
  double *q = problem->q;

  if(options->iparam[SICONOS_LCP_IPARAM_SKIP_TRIVIAL] == SICONOS_LCP_SKIP_TRIVIAL_NO)
  {
    while((i < (n - 1)) && (q[i] >= 0.)) i++;
    if((i == (n - 1)) && (q[n - 1] >= 0.))
    {
      /* TRIVIAL CASE : q >= 0
       * z = 0 and w = q is solution of LCP(q,M)
       */
      for(int j = 0 ; j < n; j++)
      {
        z[j] = 0.0;
        w[j] = q[j];
      }
      info = 0;
      options->iparam[SICONOS_IPARAM_ITER_DONE] = 0;   /* Number of iterations done */
      options->dparam[SICONOS_DPARAM_RESIDU] = 0.0; /* Error */
      if(verbose > 0)
        printf("LCP_driver_SparseBlockMatrix: found trivial solution for the LCP (positive vector q => z = 0 and w = q). \n");
      DEBUG_END("lcp_driver_SparseBlockMatrix(...)\n");
      return info;
    }
  }

  /*************************************************
   *  2 - Call specific solver (if no trivial sol.)
   *************************************************/

  /* Solver name */
  //  const char* const  name = options->solverName;
  if(verbose == 1)
    printf(" ========================== Call %s SparseBlockMatrix solver for Linear Complementarity problem ==========================\n", solver_options_id_to_name(options->solverId));

  /****** Gauss Seidel block solver ******/
  if((options->solverId) == SICONOS_LCP_NSGS_SBM)
    lcp_nsgs_SBM(problem, z, w, &info, options);
  else
  {
    fprintf(stderr, "LCP_driver_SparseBlockMatrix error: unknown solver named: %s\n", solver_options_id_to_name(options->solverId));
    exit(EXIT_FAILURE);
  }

  /*************************************************
   *  3 - Computes w = Mz + q and checks validity
   *************************************************/
  if(options->filterOn > 0)
    info = lcp_compute_error(problem, z, w, options->dparam[SICONOS_DPARAM_TOL], &(options->dparam[SICONOS_DPARAM_RESIDU]));
  DEBUG_END("lcp_driver_SparseBlockMatrix(...)\n");
  return info;

}

int lcp_driver_DenseMatrix(LinearComplementarityProblem* problem, double *z, double *w, SolverOptions* options)
{
  DEBUG_BEGIN("lcp_driver_DenseMatrix(...)\n")
  /* Note: inputs are not checked since it is supposed to be done in lcp_driver() function which calls the present one. */

  /* Checks storage type for the matrix M of the LCP */
  if(problem->M->storageType == 1)
    numerics_error("lcp_driver_DenseMatrix", "forbidden type of storage for the matrix M of the LCP");

  assert(options->isSet);

  if(verbose > 0)
    solver_options_print(options);

  /* Output info. : 0: ok -  >0: problem (depends on solver) */
  int info = -1;
  /******************************************
   *  1 - Check for trivial solution
   ******************************************/

  int i = 0;
  int n = problem->size;
  double *q = problem->q;
  if(options->iparam[SICONOS_LCP_IPARAM_SKIP_TRIVIAL] == SICONOS_LCP_SKIP_TRIVIAL_NO)
  {
    /*  if (!((options->solverId == SICONOS_LCP_ENUM) && (options->iparam[SICONOS_IPARAM_MAX_ITER] == 1 )))*/
    {
      while((i < (n - 1)) && (q[i] >= 0.)) i++;
      if((i == (n - 1)) && (q[n - 1] >= 0.))
      {
        /* TRIVIAL CASE : q >= 0
         * z = 0 and w = q is solution of LCP(q,M)
         */
        for(int j = 0 ; j < n; j++)
        {
          z[j] = 0.0;
          w[j] = q[j];
        }
        info = 0;
        options->dparam[SICONOS_DPARAM_RESIDU] = 0.0; /* Error */
        if(verbose > 0)
          printf("LCP_driver_DenseMatrix: found trivial solution for the LCP (positive vector q => z = 0 and w = q). \n");
        DEBUG_END("lcp_driver_DenseMatrix(...)\n")
        return info;
      }
    }

    // trivial solution : size-1 LCP
    if(n == 1)
    {
      double *M = problem->M->matrix0;
      w[0] = 0.;
      z[0] = -q[0] / M[0];
      info = 0;
      options->dparam[SICONOS_DPARAM_RESIDU] = 0.0; /* Error */
      if(verbose > 0)
        printf("LCP_driver_DenseMatrix: found trivial solution for the LCP (problem of size 1). \n");
      DEBUG_END("lcp_driver_DenseMatrix(...)\n")
      return info;
    }
  }

  /*************************************************
   *  2 - Call specific solver (if no trivial sol.)
   *************************************************/

  if(verbose == 1)
    printf(" ========================== Call %s solver for Linear Complementarity problem ==========================\n", solver_options_id_to_name(options->solverId));

  /****** Lemke algorithm ******/
  /* IN: itermax
     OUT: iter */
  int id = options->solverId;
  switch(id)
  {
  case SICONOS_LCP_LEMKE :
    lcp_lexicolemke(problem, z, w, &info, options);
    break;
  /****** PGS Solver ******/
  /* IN: itermax, tolerance
     OUT: iter, error */
  case SICONOS_LCP_PGS :
    lcp_pgs(problem, z, w, &info, options);
    break;
  /****** CPG Solver ******/
  /* IN: itermax, tolerance
     OUT: iter, error */
  case SICONOS_LCP_CPG:
    lcp_cpg(problem, z, w, &info, options);
    break;
  /****** Latin Solver ******/
  /* IN: itermax, tolerance, k_latin
     OUT: iter, error */
  case SICONOS_LCP_LATIN:
    lcp_latin(problem, z, w, &info, options);
    break;
  /****** Latin_w Solver ******/
  /* IN: itermax, tolerance, k_latin, relax
     OUT: iter, error */
  case SICONOS_LCP_LATIN_W:
    lcp_latin_w(problem, z, w, &info, options);
    break;
  /****** QP Solver ******/
  /* IN: tolerance
     OUT:
  */
  /* We assume that the LCP matrix M is symmetric*/
  case SICONOS_LCP_QP:
    lcp_qp(problem, z, w, &info, options);
    break;
  /****** NSQP Solver ******/
  /* IN: tolerance
     OUT:
  */
  case SICONOS_LCP_NSQP:
    lcp_nsqp(problem, z, w, &info, options);
    break;
  /****** Newton min ******/
  /* IN: itermax, tolerance
     OUT: iter, error
  */
  case SICONOS_LCP_NEWTONMIN:
    lcp_newton_min(problem, z, w, &info, options);
    break;
  /****** Newton Fischer-Burmeister ******/
  /* IN: itermax, tolerance
     OUT: iter, error
  */
  case SICONOS_LCP_NEWTON_FB_FBLSA:
    lcp_newton_FB(problem, z, w, &info, options);
    break;
  /****** Newton min + Fischer-Burmeister ******/
  /* IN: itermax, tolerance
     OUT: iter, error
  */
  case SICONOS_LCP_NEWTON_MIN_FBLSA:
    lcp_newton_minFB(problem, z, w, &info, options);
    break;
  /****** PSOR Solver ******/
  /* IN: itermax, tolerance, relax
     OUT: iter, error
  */
  case SICONOS_LCP_PSOR:
    lcp_psor(problem, z, w, &info, options);
    break;
  /****** RPGS (Regularized Projected Gauss-Seidel) Solver ******/
  /* IN: itermax, tolerance, rho
     OUT: iter, error
  */
  case SICONOS_LCP_RPGS:
    lcp_rpgs(problem, z, w, &info, options);
    break;
  /****** PATH (Ferris) Solver ******/
  /* IN: itermax, tolerance, rho
     OUT: iter, error
  */
  case SICONOS_LCP_PATH:
    lcp_path(problem, z, w, &info, options);
    break;
  /****** Enumeratif Solver ******/
  /* IN:  tolerance,
     OUT: key
  */
  case SICONOS_LCP_ENUM:
    lcp_enum(problem, z, w, &info, options);
    break;
  /****** Reformulate as AVI and use a solver by Cao and Ferris ******/
  /* IN:  tolerance, itermax
     OUT: iter, error
  */
  case SICONOS_LCP_AVI_CAOFERRIS:
    lcp_avi_caoferris(problem, z, w, &info, options);
    break;
  case SICONOS_LCP_PIVOT:
    lcp_pivot(problem, z, w, &info, options);
    break;
  case SICONOS_LCP_BARD:
    options->iparam[SICONOS_LCP_IPARAM_PIVOTING_METHOD_TYPE] = SICONOS_LCP_PIVOT_BARD;
    lcp_pivot(problem, z, w, &info, options);
    break;
  case SICONOS_LCP_MURTY:
    options->iparam[SICONOS_LCP_IPARAM_PIVOTING_METHOD_TYPE] = SICONOS_LCP_PIVOT_LEAST_INDEX;
    lcp_pivot(problem, z, w, &info, options);
    break;
  case SICONOS_LCP_PATHSEARCH:
    lcp_pathsearch(problem, z, w, &info, options);
    break;
  case SICONOS_LCP_PIVOT_LUMOD:
    lcp_pivot_lumod(problem, z, w, &info, options);
    break;
  /*error */
  case SICONOS_LCP_GAMS:
    lcp_gams(problem, z, w, &info, options);
    break;
  case SICONOS_LCP_CONVEXQP_PG:
    lcp_ConvexQP_ProjectedGradient(problem, z, w, &info, options);
    break;
  default:
  {
    fprintf(stderr, "lcp_driver_DenseMatrix error: unknown solver name: %s\n", solver_options_id_to_name(options->solverId));
    exit(EXIT_FAILURE);
  }
  }
  /*************************************************
   *  3 - Computes w = Mz + q and checks validity
   *************************************************/
  if(options->filterOn > 0)
  {
    int info_ = lcp_compute_error(problem, z, w, options->dparam[SICONOS_DPARAM_TOL], &(options->dparam[SICONOS_DPARAM_RESIDU]));
    if(info <= 0)  /* info was not set or the solver was happy */
      info = info_;
  }
  DEBUG_END("lcp_driver_DenseMatrix(...)\n")
  return info;

}

int linearComplementarity_driver(LinearComplementarityProblem* problem, double *z, double *w, SolverOptions* options)
{
  assert(options && "lcp_driver : null input for solver options");
  DEBUG_BEGIN("linearComplementarity_driver(...)\n");
  /* Checks inputs */
  assert(problem && z && w &&
         "lcp_driver : input for LinearComplementarityProblem and/or unknowns (z,w)");

  /* Output info. : 0: ok -  >0: problem (depends on solver) */
  int info = -1;
  /* Switch to DenseMatrix or SparseBlockMatrix solver according to the type of storage for M */
  /* Storage type for the matrix M of the LCP */

  int storageType = problem->M->storageType;
  DEBUG_PRINTF("storageType = %i\n", storageType);
  /* Sparse Block Storage */
  if(storageType == 1)
  {
    info = lcp_driver_SparseBlockMatrix(problem, z, w, options);
  }
  else
  {
    info = lcp_driver_DenseMatrix(problem, z, w, options);
  }
  DEBUG_END("linearComplementarity_driver(...)\n");
  return info;
}
