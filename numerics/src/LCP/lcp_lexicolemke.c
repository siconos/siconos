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

#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "pivot-utils.h"
#include "numerics_verbose.h"
#include "LinearComplementarityProblem.h"
#include "LCP_Solvers.h"
#include "lcp_cst.h"
#include "SolverOptions.h"
#include "NumericsMatrix.h"

#include "numerics_verbose.h"
//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES

//#define MAX_PIVOT
//#define INV_PIVOT

#include "debug.h"

void lcp_lexicolemke(LinearComplementarityProblem* problem, double *zlem , double *wlem , int *info , SolverOptions* options)
{
  /* matrix M of the lcp */
  double * M = problem->M->matrix0;
  assert(M);
  /* size of the LCP */
  int dim = problem->size;
  assert(dim>0);
  int dim2 = 2 * (dim + 1);

  int i, drive, block, Ifound;
  int ic, jc;
  int ITER;
  int nobasis;
  int itermax = options->iparam[0];

  i=0;
  int n = problem->size;
  double *q = problem->q;

  while ((i < (n - 1)) && (q[i] >= 0.)) 
    i++;

  if ((i == (n - 1)) && (q[n - 1] >= 0.))
  {

    /* TRIVIAL CASE : q >= 0
     * z = 0 and w = q is solution of LCP(q,M)
     */
    for (int j = 0 ; j < n; j++)
    {
      zlem[j] = 0.0;
      wlem[j] = q[j];
    }
    *info = 0;
    options->iparam[1] = 0;   /* Number of iterations done */
    options->dparam[1] = 0.0; /* Error */
    if (verbose > 0)
      printf("lcp_lexicolemke: found trivial solution for the LCP (positive vector q => z = 0 and w = q). \n");
    return ;
  }

  double z0, zb, delta_lexico;
  double pivot, tovip, ratio;
  double tmp;
  int *basis;
  double** A;

  /*output*/

  options->iparam[1] = 0;

  /* Allocation */

  unsigned* candidate_pivots_indx = (unsigned*)malloc(dim * sizeof(unsigned));
  basis = (int *)malloc(dim * sizeof(int));
  A = (double **)malloc(dim * sizeof(double*));

  for (ic = 0 ; ic < dim; ++ic)
    A[ic] = (double *)malloc(dim2 * sizeof(double));

  /* construction of A matrix such that
   * A = [ q | Id | -d | -M ] with d = (1,...1)
   */

  /* We need to init only the part corresponding to Id */
  for (ic = 0 ; ic < dim; ++ic)
    for (jc = 1 ; jc <= dim; ++jc)
      A[ic][jc] = 0.0;

  double max_elt_M = 0.;
  for (ic = 0 ; ic < dim; ++ic)
    for (jc = 0 ; jc < dim; ++jc)
    {
      A[ic][jc + dim + 2] = -M[dim * jc + ic];
      if (fabs(M[dim * jc + ic]) > max_elt_M) max_elt_M = fabs(M[dim * jc + ic]);
    }

  double lexico_tol_diff;
  double lexico_tol_elt;
  if (options->iparam[2] == 0)
  {
    lexico_tol_diff = max_elt_M*DBL_EPSILON;
    lexico_tol_elt = max_elt_M*DBL_EPSILON;
  }
  else
  {
    lexico_tol_diff = options->dparam[2];
    lexico_tol_elt = options->dparam[3];
  }
  assert(problem->q);

  for (ic = 0 ; ic < dim; ++ic) A[ic][0] = problem->q[ic];

  for (ic = 0 ; ic < dim; ++ic) A[ic][ic + 1 ] =  1.0;
  for (ic = 0 ; ic < dim; ++ic) A[ic][dim + 1] = -1.0;

  DEBUG_PRINT("total matrix\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("%1.2e ", A[i][j]) }
      DEBUG_PRINT("\n")});
  /* End of construction of A */

  Ifound = 0;


  for (ic = 0 ; ic < dim  ; ++ic) basis[ic] = ic + 1;

  drive = dim + 1;
  block = 0;
  z0 = A[block][0];
  ITER = 0;

  /* Start research of argmin lexico */
  /* With this first step the covering vector enter in the basis */

  for (ic = 1 ; ic < dim ; ++ic)
  {
    zb = A[ic][0];
    if (zb < z0)
    {
      z0    = zb;
      block = ic;
    }
    else if (zb == z0)
    {
      for (jc = 0 ; jc < dim ; ++jc)
      {
        delta_lexico = A[block][1 + jc] - A[ic][1 + jc];
        if (delta_lexico < 0.)
        {
          break;
        }
        else if (delta_lexico > 0)
        {
          block = ic;
          break;
        }
      }
    }
  }

  /* Stop research of argmin lexico */
  DEBUG_PRINTF("Pivoting %i and %i\n", block, drive);

  pivot = A[block][drive];
  tovip = 1.0 / pivot;

  /* Pivot < block , drive > */

#ifdef INV_PIVOT
  A[block][drive] = tovip;
#else
  A[block][drive] = 1.;
#endif
  for (ic = 0       ; ic < drive ; ++ic) A[block][ic] = A[block][ic] * tovip;
  for (ic = drive + 1 ; ic < dim2  ; ++ic) A[block][ic] = A[block][ic] * tovip;

  /* */

  for (ic = 0 ; ic < block ; ++ic)
  {
    tmp = A[ic][drive];
    for (jc = 0 ; jc < dim2 ; ++jc) A[ic][jc] -=  tmp * A[block][jc];
  }
  for (ic = block + 1 ; ic < dim ; ++ic)
  {
    tmp = A[ic][drive];
    for (jc = 0 ; jc < dim2 ; ++jc) A[ic][jc] -=  tmp * A[block][jc];
  }

   nobasis = basis[block];
  basis[block] = drive;

  DEBUG_EXPR_WE( DEBUG_PRINT("new basis: ")
      for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%i ", basis[i])}
      DEBUG_PRINT("\n"));
  DEBUG_PRINT("total matrix\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("%1.2e ", A[i][j]) }
      DEBUG_PRINT("\n")});
  DEBUG_PRINT("lexico_mat\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF(ANSI_COLOR_YELLOW "% 1.1e " ANSI_COLOR_RESET, A[i][drive]);
      for(unsigned int j = 1 ; j <= dim; ++j) {
      if (fabs(A[i][j]) > 2.2e-16) {DEBUG_PRINTF(ANSI_COLOR_YELLOW " % 2.f " ANSI_COLOR_RESET, A[i][j])}
      else if (A[i][j] == 0.) { DEBUG_PRINT(ANSI_COLOR_BLUE " . " ANSI_COLOR_RESET) }
      else { DEBUG_PRINT(ANSI_COLOR_RED " X " ANSI_COLOR_RESET) } }
      DEBUG_PRINT("\n")});

  while (ITER < itermax && !Ifound)
  {

    ++ITER;

    if (nobasis < dim + 1)      drive = nobasis + (dim + 1);
    else if (nobasis > dim + 1) drive = nobasis - (dim + 1);

    DEBUG_EXPR_WE( DEBUG_PRINT("basis= "); for (unsigned i = 0; i < dim; ++i) { DEBUG_PRINTF("%s%d ", basis_to_name(basis[i], dim), basis_to_number(basis[i], dim)); } DEBUG_PRINT("\n"));

    /* Start research of argmin lexico for minimum ratio test */
    ratio = 1e20;
    block = -1;

    unsigned nb_candidate = 0;
    DEBUG_EXPR_WE(unsigned max_pivot_helped = 0;)
    for (ic = 0 ; ic < dim ; ++ic)
    {
      zb = A[ic][drive];
      if (zb > 0.0)
      {
        z0 = A[ic][0] / zb;
        if (z0 > ratio) continue;
        if (z0 < ratio)
        {
          ratio = z0;
          block = ic;
          nb_candidate = 0;
          DEBUG_EXPR_WE(max_pivot_helped = 0;)
        }
        else
        {
          candidate_pivots_indx[nb_candidate++] = ic;
        }
      }
    }

    if (nb_candidate > 0)
    {
      DEBUG_PRINTF("pivot_selection_lemke :: lexicomin %d candidates, ratio = %e, drive = %d\n", nb_candidate, ratio, drive);
      for (unsigned k = 0; k < nb_candidate; ++k)
      {
        unsigned var = candidate_pivots_indx[k];
        double candidate_pivot = A[var][drive];
        for (jc = 1 ; jc < dim + 1 ; ++jc)
        {
          assert(block >=0 && "lcp_lexicolemke: block <0");
          //                dblock = (A[block][jc] / A[block][drive]) - (A[var][jc] / candidate_pivot);
          if (fabs(A[block][jc]) < lexico_tol_elt) // XXX TOL
          {
            if (A[var][jc] < -lexico_tol_elt)
            {
              /*  delta_lexico > 0 (since pivot are always >0., => new lexicomin  */
              block = var;
              DEBUG_EXPR_WE(max_pivot_helped = 0;)
              break;
            }
            else if (A[var][jc] > lexico_tol_elt)
              /* delta_lexico < 0 => lexicomin does not change */
            {
              DEBUG_EXPR_WE(max_pivot_helped = 0;)
              break;
            }
            else /* delta_lexico not conclusive => equality */
            {
              continue;
            }
          }
          else if (fabs(A[var][jc]) < lexico_tol_elt)
          {
            if (A[block][jc] > lexico_tol_elt)
            {
              /*  delta_lexico > 0 (since pivot are always >0., => new lexicomin  */
              block = var;
              DEBUG_EXPR_WE(max_pivot_helped = 0;)
              break;
            }
            else if (A[block][jc] < -lexico_tol_elt)
              /* delta_lexico < 0 => lexicomin does not change */
            {
              DEBUG_EXPR_WE(max_pivot_helped = 0;)
              break;
            }
            else /* delta_lexico not conclusive => equality */
            {
              continue;
            }
          }
          else /* really compute delta_lexico */
          {
            delta_lexico = (A[block][jc] * candidate_pivot) - (A[var][jc] * A[block][drive]);
            DEBUG_EXPR_WE(if ((delta_lexico != 0.) && (fabs(delta_lexico) < 1e-10) && (fabs(delta_lexico) > lexico_tol_diff)) { printf("pivot_selection_lemke :: very small difference in lexicomin: %2.2e\n", delta_lexico);
                unsigned block_number = basis_to_number(basis[block], dim); unsigned var_number =  basis_to_number(basis[var], dim);
                char* block_name = basis_to_name(basis[block], dim); char* var_name = basis_to_name(basis[var], dim);
                printf("lexicomin: A[%s%d][jc] / A[%s%d][drive] = %e / %e vs A[%s%d][jc] / A[%s%d][drive] = %e / %e\n", block_name, block_number, block_name, block_number, A[block][jc], A[block][drive], var_name, var_number, var_name, var_number, A[var][jc], A[var][drive]);});
            if (delta_lexico < -lexico_tol_diff) break;
            else if (delta_lexico > lexico_tol_diff)
            {
              DEBUG_PRINTF("pivot_selection_lemke :: lexicomin change var block changes %s%d from %s%d, delta_lexico = %2.2e, new pivot = %e\n", basis_to_name(basis[var], dim), basis_to_number(basis[var], dim), basis_to_name(basis[block], dim), basis_to_number(basis[block], dim), delta_lexico, A[var][drive]);
              block = var;
              DEBUG_EXPR_WE(max_pivot_helped = 0;)
              break;
            }
#ifdef MAX_PIVOT
            else if (delta_lexico != 0.)
            {
              if ((A[block][drive] < 1e-10) && (candidate_pivot > A[block][drive]))
              {
                DEBUG_PRINTF("pivot_selection_lemke :: lexicomin small difference %2.2e, taking largest pivot %e > %e (var %s%d vs %s%d)\n", delta_lexico, candidate_pivot, A[block][drive], basis_to_name(basis[var], dim), basis_to_number(basis[var], dim), basis_to_name(basis[block], dim), basis_to_number(basis[block], dim));
                block = var;
                DEBUG_EXPR_WE(max_pivot_helped = 1;)
                break;
              }

            }
#endif
          }
        }
      }
    }


    DEBUG_EXPR_WE(if (max_pivot_helped) {DEBUG_PRINT("pivot_selection_lemke :: lexicomin MAX PIVOT HELPED!\n");});
    if (block == -1)
    {
      Ifound = 1;
      DEBUG_PRINT("The pivot column is nonpositive !\n"
          "It either means that the algorithm failed or that the LCP is infeasible\n"
          "Check the class of the M matrix to find out the meaning of this\n");
      break;
    }

    DEBUG_PRINTF("leaving variable %s%d entering variable %s%d\n", basis_to_name(basis[block], dim), basis_to_number(basis[block], dim), basis_to_name(drive, dim), basis_to_number(drive, dim));
    if (basis[block] == dim + 1) Ifound = 1;

    /* Pivot < block , drive > */

    pivot = A[block][drive];
    tovip = 1.0 / pivot;
#ifdef INV_PIVOT
    A[block][drive] = tovip;
#else
    A[block][drive] = 1.;
#endif

    for (ic = 0       ; ic < drive ; ++ic) A[block][ic] = A[block][ic] * tovip;
    for (ic = drive + 1 ; ic < dim2  ; ++ic) A[block][ic] = A[block][ic] * tovip;

    /* */

    for (ic = 0 ; ic < block ; ++ic)
    {
      tmp = A[ic][drive];
      for (jc = 0 ; jc < dim2 ; ++jc) A[ic][jc] -=  tmp * A[block][jc];
    }
    for (ic = block + 1 ; ic < dim ; ++ic)
    {
      tmp = A[ic][drive];
      for (jc = 0 ; jc < dim2 ; ++jc) A[ic][jc] -=  tmp * A[block][jc];
    }

    nobasis = basis[block];
    basis[block] = drive;

    DEBUG_EXPR_WE( DEBUG_PRINT("new basis: ")
      for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%i ", basis[i])}
      DEBUG_PRINT("\n"));

    DEBUG_PRINT("total matrix\n");
    DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("%1.2e ", A[i][j]) }
      DEBUG_PRINT("\n")});
    DEBUG_PRINT("lexico_mat\n");
    DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
        { DEBUG_PRINTF(ANSI_COLOR_YELLOW "% 1.1e " ANSI_COLOR_RESET, A[i][drive]);
        for(unsigned int j = 1 ; j <= dim; ++j) {
        if (fabs(A[i][j]) > 2.2e-16) {DEBUG_PRINTF(ANSI_COLOR_YELLOW " % 2.f " ANSI_COLOR_RESET, A[i][j])}
        else if (A[i][j] == 0.) { DEBUG_PRINT(ANSI_COLOR_BLUE " . " ANSI_COLOR_RESET) }
        else { DEBUG_PRINT(ANSI_COLOR_RED " X " ANSI_COLOR_RESET) } }
        DEBUG_PRINT("\n")});


  } /* end while*/

  DEBUG_EXPR_WE( DEBUG_PRINT("new basis: ")
      for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%i ", basis[i])}
      DEBUG_PRINT("\n"));

  DEBUG_PRINT("total matrix\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("%1.2e ", A[i][j]) }
      DEBUG_PRINT("\n")});

  for (ic = 0 ; ic < dim; ++ic)
  {
    drive = basis[ic];
    if (drive < dim + 1)
    {
      zlem[drive - 1] = 0.0;
      wlem[drive - 1] = A[ic][0];
    }
    else if (drive > dim + 1)
    {
      zlem[drive - dim - 2] = A[ic][0];
      wlem[drive - dim - 2] = 0.0;
    }
  }

  options->iparam[SICONOS_IPARAM_ITER_DONE] = ITER;

  if (Ifound) *info = 0;
  else *info = 1;

  free(basis);

  for (i = 0 ; i < dim ; ++i) free(A[i]);
  free(A);
  free(candidate_pivots_indx);
}


int linearComplementarity_lexicolemke_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the Lemke Solver\n");
  }

  options->solverId = SICONOS_LCP_LEMKE;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  options->dWork = NULL;
  solver_options_nullify(options);
  options->dparam[0] = 1e-6;
  options->iparam[0] = 10000;
  return 0;
}
