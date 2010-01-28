/* Siconos-Numerics, Copyright INRIA 2005-2010.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "Numerics_Options.h"
#include "NonSmoothDrivers.h"

#define LCP_LEMKE "Lemke"
#define LCP_NSGS_SBM "NSGS_SBM"
#define LCP_PGS "PGS"
#define LCP_CPG "CPG"
#define LCP_LATIN "Latin"
#define LCP_LATIN_W "Latin_w"
#define LCP_QP "QP"
#define LCP_NSQP "NSQP"
#define LCP_NEWTONMIN "NewtonMin"
#define LCP_NEWTONFB "NewtonFB"
#define LCP_PSOR "PSOR"
#define LCP_RPGS "RPGS"
#define LCP_PATH "PATH"
#define LCP_ENUM "ENUM"

int lcp_driver_SparseBlockMatrix(LinearComplementarity_Problem* problem, double *z , double *w, Solver_Options* options)
{
  /* Checks storage type for the matrix M of the LCP */
  if (problem->M->storageType == 0)
    numericsError("lcp_driver_SparseBlockMatrix", "forbidden type of storage for the matrix M of the LCP");

  /*
    The options for the global "block" solver are defined in options[0].\n
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
  while ((i < (n - 1)) && (q[i] >= 0.)) i++;
  if ((i == (n - 1)) && (q[n - 1] >= 0.))
  {
    /* TRIVIAL CASE : q >= 0
     * z = 0 and w = q is solution of LCP(q,M)
     */
    for (int j = 0 ; j < n; j++)
    {
      z[j] = 0.0;
      w[j] = q[j];
    }
    info = 0;
    options[0].iparam[1] = 0;   /* Number of iterations done */
    options[0].dparam[1] = 0.0; /* Error */
    if (verbose > 0)
      printf("LCP_driver_SparseBlockMatrix: found trivial solution for the LCP (positive vector q => z = 0 and w = q). \n");
    return info;
  }

  /*************************************************
  *  2 - Call specific solver (if no trivial sol.)
  *************************************************/

  /* Solver name */
  char * name = options[0].solverName;
  if (verbose == 1)
    printf(" ========================== Call %s SparseBlockMatrix solver for Linear Complementarity problem ==========================\n", name);

  /****** Gauss Seidel block solver ******/
  if (strcmp(name , LCP_NSGS_SBM) == 0)
    lcp_nsgs_SBM(problem, z , w , &info , options);
  else
  {
    fprintf(stderr, "LCP_driver_SparseBlockMatrix error: unknown solver named: %s\n", name);
    exit(EXIT_FAILURE);
  }

  /*************************************************
   *  3 - Computes w = Mz + q and checks validity
   *************************************************/
  if (options[0].filterOn > 0)
    info = lcp_compute_error(problem, z, w, options[0].dparam[0], &(options[0].dparam[1]));

  return info;

}

int lcp_driver_DenseMatrix(LinearComplementarity_Problem* problem, double *z , double *w, Solver_Options* options)
{
  /* Note: inputs are not checked since it is supposed to be done in lcp_driver() function which calls the present one. */

  /* Checks storage type for the matrix M of the LCP */
  if (problem->M->storageType == 1)
    numericsError("lcp_driver_DenseMatrix", "forbidden type of storage for the matrix M of the LCP");

  /* If the options for solver have not been set, read default values in .opt file */
  int NoDefaultOptions = options->isSet; /* true(1) if the Solver_Options structure has been filled in else false(0) */

  if (NoDefaultOptions == 0)
  {
    readSolverOptions(0, options);
    options->filterOn = 1;
  }

  if (verbose > 0)
    printSolverOptions(options);

  /* Output info. : 0: ok -  >0: problem (depends on solver) */
  int info = 1;
  /******************************************
   *  1 - Check for trivial solution
   ******************************************/

  int i = 0;
  int n = problem->size;
  double *q = problem->q;
  while ((i < (n - 1)) && (q[i] >= 0.)) i++;
  if ((i == (n - 1)) && (q[n - 1] >= 0.))
  {
    /* TRIVIAL CASE : q >= 0
     * z = 0 and w = q is solution of LCP(q,M)
     */
    for (int j = 0 ; j < n; j++)
    {
      z[j] = 0.0;
      w[j] = q[j];
    }
    info = 0;
    options->dparam[1] = 0.0; /* Error */
    if (verbose > 0)
      printf("LCP_driver_DenseMatrix: found trivial solution for the LCP (positive vector q => z = 0 and w = q). \n");
    return info;
  }

  /*************************************************
   *  2 - Call specific solver (if no trivial sol.)
   *************************************************/

  /* Solver name */
  char * name = options->solverName;

  if (verbose == 1)
    printf(" ========================== Call %s solver for Linear Complementarity problem ==========================\n", name);

  /****** Lemke algorithm ******/
  /* IN: itermax
     OUT: iter */
  if (strcmp(name , LCP_LEMKE) == 0 || strcmp(name , "LexicoLemke") == 0)
    lcp_lexicolemke(problem, z , w , &info , options);

  /****** PGS Solver ******/
  /* IN: itermax, tolerance
     OUT: iter, error */
  else   if (strcmp(name , LCP_PGS) == 0)
    lcp_pgs(problem, z , w , &info , options);

  /****** CPG Solver ******/
  /* IN: itermax, tolerance
     OUT: iter, error */
  else   if (strcmp(name , LCP_CPG) == 0)
    lcp_cpg(problem, z , w , &info , options);

  /****** Latin Solver ******/
  /* IN: itermax, tolerance, k_latin
     OUT: iter, error */
  else   if (strcmp(name , LCP_LATIN) == 0)
    lcp_latin(problem, z , w , &info , options);

  /****** Latin_w Solver ******/
  /* IN: itermax, tolerance, k_latin, relax
     OUT: iter, error */
  else   if (strcmp(name , LCP_LATIN_W) == 0)
    lcp_latin_w(problem, z , w , &info , options);

  /****** QP Solver ******/
  /* IN: tolerance
     OUT:
  */
  /* We assume that the LCP matrix M is symmetric*/
  else   if (strcmp(name , LCP_QP) == 0)
    lcp_qp(problem, z , w , &info , options);

  /****** NSQP Solver ******/
  /* IN: tolerance
     OUT:
  */
  else   if (strcmp(name , LCP_NSQP) == 0)
    lcp_nsqp(problem, z , w , &info , options);

  /****** Newton min ******/
  /* IN: itermax, tolerance
     OUT: iter, error
  */
  else   if (strcmp(name , LCP_NEWTONMIN) == 0)
    lcp_newton_min(problem, z , w , &info , options);

  /****** Newton Fischer-Burmeister ******/
  /* IN: itermax, tolerance
     OUT: iter, error
  */
  else   if (strcmp(name , LCP_NEWTONFB) == 0)
    lcp_newton_FB(problem, z , w , &info , options);

  /****** PSOR Solver ******/
  /* IN: itermax, tolerance, relax
     OUT: iter, error
  */
  else   if (strcmp(name , LCP_PSOR) == 0)
    lcp_psor(problem, z , w , &info , options);

  /****** RPGS (Regularized Projected Gauss-Seidel) Solver ******/
  /* IN: itermax, tolerance, rho
     OUT: iter, error
  */
  else   if (strcmp(name , LCP_RPGS) == 0)
    lcp_rpgs(problem, z , w , &info , options);

  /****** PATH (Ferris) Solver ******/
  /* IN: itermax, tolerance, rho
     OUT: iter, error
  */
  else   if (strcmp(name , LCP_PATH) == 0)
    lcp_path(problem, z , w , &info , options);

  /****** Enumeratif Solver ******/
  /* IN:  tolerance,
     OUT: key
  */
  else   if (strcmp(name , LCP_ENUM) == 0)
    lcp_enum(problem, z , w , &info , options);

  /*error */
  else
  {
    fprintf(stderr, "lcp_driver_DenseMatrix error: unknown solver name: %s\n", name);
    exit(EXIT_FAILURE);
  }

  /*************************************************
   *  3 - Computes w = Mz + q and checks validity
   *************************************************/
  if (options->filterOn > 0)
    lcp_compute_error(problem, z, w, options->dparam[0], &(options->dparam[1]));

  return info;

}

int linearComplementarity_driver(LinearComplementarity_Problem* problem, double *z , double *w, Solver_Options* options,  Numerics_Options* global_options)
{
  if (options == NULL || global_options == NULL)
    numericsError("lcp_driver", "null input for solver and/or global options");

  /* Set global options */
  setNumericsOptions(global_options);

  /* Checks inputs */
  if (problem == NULL || z == NULL || w == NULL)
    numericsError("lcp_driver", "null input for LinearComplementarity_Problem and/or unknowns (z,w)");

  /* Output info. : 0: ok -  >0: problem (depends on solver) */
  int info = -1;

  /* Switch to DenseMatrix or SparseBlockMatrix solver according to the type of storage for M */
  /* Storage type for the matrix M of the LCP */

  /*   LinearComplementarity_display(problem); */

  int storageType = problem->M->storageType;



  /* Sparse Block Storage */
  if (storageType == 1)
  {
    info = lcp_driver_SparseBlockMatrix(problem, z , w, options);
  }
  else
  {
    info = lcp_driver_DenseMatrix(problem, z , w, options);
  }
  return info;
}

