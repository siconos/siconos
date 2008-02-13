/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "Numerics_Options.h"
#include "NonSmoothDrivers.h"

int lcp_driver(LinearComplementarity_Problem* problem, double *z , double *w, Solver_Options* options, Numerics_Options* global_options)
{
  if (options == NULL || global_options == NULL)
    numericsError("lcp_driver", "null input for solver and/or global options");

  /* Set global options */
  setNumericsOptions(global_options);

  /* If the options for solver have not been set, read default values in .opt file */
  int NoDefaultOptions = options->isSet; /* true(1) if the Solver_Options structure has been filled in else false(0) */

  if (!NoDefaultOptions)
    readSolverOptions(0, options);

  if (verbose > 0)
    printSolverOptions(options);

  /* Output info. : 0: ok -  >0: problem (depends on solver) */
  int info = 1;

  /* Checks inputs */
  if (problem == NULL || z == NULL || w == NULL)
    numericsError("lcp_driver", "null input for LinearComplementarity_Problem and/or unknowns (z,w)");

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
    options->iparam[2] = 0;   /* Number of iterations done */
    options->dparam[2] = 0.0; /* Error */
    if (verbose > 0)
      printf("LCP_driver: found trivial solution for the LCP (positive vector q => z = 0 and w = q). \n");
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
  if (strcmp(name , "Lemke") == 0)
    lcp_lexicolemke(problem, z , w , &info , options);

  /****** PGS Solver ******/
  /* IN: itermax, tolerance
     OUT: iter, error */
  else   if (strcmp(name , "PGS") == 0)
    lcp_pgs(problem, z , w , &info , options);

  /****** CPG Solver ******/
  /* IN: itermax, tolerance
     OUT: iter, error */
  else   if (strcmp(name , "CPG") == 0)
    lcp_cpg(problem, z , w , &info , options);

  /****** Latin Solver ******/
  /* IN: itermax, tolerance, k_latin
     OUT: iter, error */
  else   if (strcmp(name , "Latin") == 0)
    lcp_latin(problem, z , w , &info , options);

  /****** Latin_w Solver ******/
  /* IN: itermax, tolerance, k_latin, relax
     OUT: iter, error */
  else   if (strcmp(name , "Latin_w") == 0)
    lcp_latin_w(problem, z , w , &info , options);

  /****** QP Solver ******/
  /* IN: tolerance
     OUT:
  */
  /* We assume that the LCP matrix M is symmetric*/
  else   if (strcmp(name , "QP") == 0)
    lcp_qp(problem, z , w , &info , options);

  /****** NSQP Solver ******/
  /* IN: tolerance
     OUT:
  */
  else   if (strcmp(name , "NSQP") == 0)
    lcp_nsqp(problem, z , w , &info , options);

  /****** Newton min ******/
  /* IN: itermax, tolerance
     OUT: iter, error
  */
  else   if (strcmp(name , "NewtonMin") == 0)
    lcp_newton_min(problem, z , w , &info , options);

  /****** Newton Fischer-Burmeister ******/
  /* IN: itermax, tolerance
     OUT: iter, error
  */
  else   if (strcmp(name , "NewtonFB") == 0)
    lcp_newton_FB(problem, z , w , &info , options);

  /****** PSOR Solver ******/
  /* IN: itermax, tolerance, relax
     OUT: iter, error
  */
  else   if (strcmp(name , "PSOR") == 0)
    lcp_psor(problem, z , w , &info , options);

  /****** RPGS (Regularized Projected Gauss-Seidel) Solver ******/
  /* IN: itermax, tolerance, rho
     OUT: iter, error
  */
  else   if (strcmp(name , "RPGS") == 0)
    lcp_rpgs(problem, z , w , &info , options);

  /****** PATH (Ferris) Solver ******/
  /* IN: itermax, tolerance, rho
     OUT: iter, error
  */
  else   if (strcmp(name , "Path") == 0)
    lcp_path(problem, z , w , &info , options);

  /*error */
  else
    printf("LCP_driver error: unknown solver named: %s\n", name);

  /*************************************************
   *  3 - Check solution validity
   *************************************************/

  /* Warning: it depends on the chosen solver */

  /* Not done for:  PGS, RPGS */
  if ((strcmp(name , "PGS") != 0) && (strcmp(name , "CPG") != 0))
  {
    double tolerance = options->dparam[0];
    info = filter_result_LCP(problem, z, w, tolerance);
  }
  return info;

}
