#include <assert.h>  // for assert
#include <float.h>   // for DBL_EPSILON
#include <math.h>    // for fabs
#ifndef __cplusplus
#include <stdbool.h>  // for false
#endif
#include <stdio.h>   // for printf
#include <stdlib.h>  // for free, malloc

#include "LCP_Solvers.h"                   // for lcp_compute_error, lcp_pgs
#include "LinearComplementarityProblem.h"  // for LinearComplementarityProblem
#include "NumericsFwd.h"                   // for SolverOptions, LinearCompl...
#include "NumericsMatrix.h"                // for NM_get_value, NM_row_prod_...
#include "SiconosBlas.h"                   // for cblas_dcopy, cblas_dnrm2
#include "SolverOptions.h"                 // for SolverOptions, SICONOS_DPA...
#include "numerics_verbose.h"              // for verbose
#include "siconos_debug.h"  

#include "graph_tools.h"

void lcp_pgs_graph_permut_equitable_opti(LinearComplementarityProblem *problem, double *z, double *w, int *info,
                                         SolverOptions *options) {

  NumericsMatrix *M = problem->M;
  double *q = problem->q;

  assert(M);
  assert(q);
  assert(options->iSize > 1);
  assert(options->dSize > 1);

  int n = problem->size;

  /* Solver parameters */
  int itermax = options->iparam[SICONOS_IPARAM_MAX_ITER];
  double tol = options->dparam[SICONOS_DPARAM_TOL];

  options->iparam[SICONOS_IPARAM_ITER_DONE] = 0;
  options->dparam[SICONOS_DPARAM_RESIDU] = 0.0;

  /* Precompute all 1 / M_{i,i} */
  double *diag = (double *)malloc((size_t)n * sizeof(double));
  NM_get_invdiag(n, info, M, diag);

  /* Check if diagonal has a zero */
  if (*info == 2) {
    return;
  }

  /* Graph coloring + permutation */
  int err_not_symmetric;
  size_t n_colors = 0;
  size_t *sum_sizes = NULL;
  size_t *inv_permutation = (size_t *)malloc((size_t)n * sizeof(size_t));
  
  err_not_symmetric = color_graph_permut_equitable(n, M, &n_colors, &sum_sizes, inv_permutation);

  /* Matrix not symmetric */
  if (err_not_symmetric == 1) {
    // printf("ERROR: matrix is not symmetric.\n");
    free(inv_permutation);
    free(diag);
    return;
  }

  /* printf("Number of colors = %d\n", n_colors);
  printf("sum_sizes = [");
  for (int i = 0; i < n_colors + 1; i++) {
    printf(" %d ", sum_sizes[i]);
  }
  printf("]\n"); */
  size_t *permutation = (size_t *)malloc((size_t)n * sizeof(size_t));
  for (size_t i = 0; i < (size_t)n; i++) {
    permutation[inv_permutation[i]] = i;
  }

  /* Check coloring */
  /* int *colors = (int *)malloc(n * sizeof(int));
  int current_color = 0;
  for (int i = 0; i < n; i++) {
    if (i == sum_sizes[current_color + 1]) current_color += 1;
    colors[inv_permutation[i]] = current_color;
  } */

  /* printf("Colors = [");
  for (int i = 0; i < n; i++) {
    printf(" %d ", colors[i]);
  }
  printf("]\n"); */

  /* TO CHECK COLORING
  
  CSparseMatrix *sparse;
  if (M->matrix2->origin == NSM_CSR)
  {
    sparse = NM_csr(M);
  }
  else
  {
    sparse = NM_csc_trans(M);
  }

  CS_INT *Mp = sparse->p;
  CS_INT *Mi = sparse->i;
  double *Mx = sparse->x;

  size_t col;

  for (int row = 0; row < n; row++)
  {
    for (CS_INT p = Mp[row]; p < Mp[row + 1]; ++p)
    {
      if (Mi[p] == row) printf("Same vertex\n");
      else if (colors[row] == colors[Mi[p]]) printf("ERROR: %d %d\n", row, Mi[p]);
      
    }

  } */

  /* 2-norm of q, to normalize error */
  double norm_q = cblas_dnrm2(n, q, 1); 
  if (fabs(norm_q) <= DBL_EPSILON) norm_q = 1.;

  /* Solver variables */
  int iter = 0;
  double err = 1.;

  double zi;
  double right = 0.;
  size_t i, row;
  double *lefts = (double *)malloc((size_t)n * sizeof(double));
  int is_last = 0;

  // double time = omp_get_wtime();

  
  /* size_t n_threads = (size_t)omp_get_max_threads();
  size_t *start_lines = (size_t *)malloc(n_colors * (n_threads + 1) * sizeof(size_t));
  size_t n_lines;

  start_lines[0] = 0;

  for (size_t color = 0; color < n_colors; color++) {
    n_lines = sum_sizes[color + 1] - sum_sizes[color];
    printf(" %ld ", n_lines);

    // start_lines[color * (n_threads + 1)] = 0;
    for (size_t thread = 0; thread < n_threads; thread++) {

      if (thread < n_lines % n_threads) {
        start_lines[color * (n_threads + 1) + thread + 1] = start_lines[color * (n_threads + 1) + thread] +  n_lines / n_threads + 1;
      } else {
        start_lines[color * (n_threads + 1) + thread + 1] = start_lines[color * (n_threads + 1) + thread] +  n_lines / n_threads;
      }
    }
  }
  printf("\n"); */

  /* Compute which lines will be treated by each thread */
  size_t n_threads = (size_t)omp_get_max_threads();
  size_t *start_lines = (size_t *)malloc((n_colors * n_threads + 1) * sizeof(size_t));
  size_t n_lines;

  start_lines[0] = 0;

  for (size_t color = 0; color < n_colors; color++) {
    n_lines = sum_sizes[color + 1] - sum_sizes[color];
    // printf(" %ld ", n_lines);

    // start_lines[color * (n_threads + 1)] = 0;
    for (size_t thread = 1; thread <= n_threads; thread++) {

      if (thread - 1 < n_lines % n_threads) {
        start_lines[color * n_threads + thread] = start_lines[color * n_threads + thread - 1] +   n_lines / n_threads + 1;
      } else {
        start_lines[color * n_threads + thread] = start_lines[color * n_threads + thread - 1] +   n_lines / n_threads;
      }
    }
  }

  // for (size_t i = 0; i < (n_colors * n_threads + 1); i++) printf(" %ld ", start_lines[i]);

  // printf("\n");

  size_t thread;

  #pragma omp parallel default(none) private(zi, right, i, row, thread, is_last) shared(iter, itermax, q, z, start_lines, n_colors, norm_q, inv_permutation, permutation, diag, M, n, sum_sizes, n_threads, w, lefts, err, tol)
  {
    thread = (size_t)omp_get_thread_num();

    is_last = 0;
    iter = 0;
    err = 1.;

    // printf("[%ld]\n", thread);

    /* for (size_t color = 0; color < n_colors; color++) {

      printf("[%ld]", thread);
      for (size_t i = start_lines[color * n_threads + thread]; i < start_lines[color * n_threads + thread + 1]; i++) {
        printf(" %ld ", i); 
      }
      printf("\n");
    }

    printf("[%d] iter = %d, is_last = %d\n", thread, iter, is_last); */

    /* Start solving */
    while ((iter < itermax) && (is_last == 0)) {

      if (err < tol) {
        is_last = 1;
      }

      #pragma omp barrier

      #pragma omp single
      err = 0.;
    
      for (size_t color = 0; color < n_colors; color++) {

        // printf("[%d] Starting color %d\n", thread, color);

        for (size_t i = start_lines[color * n_threads + thread]; i < start_lines[color * n_threads + thread + 1]; i++) {

          row = inv_permutation[i];
          right = q[row];
          DEBUG_PRINTF("zi = %e\n", zi);

          w[row] = lefts[row];

          lefts[row] = 0.;

          // printf("BEFORE ::: [%d] Color %d, lefts[%d] = %e, right = %e, w[%d] = %e\n", thread, color, row, lefts[row], right, row, w[row]);

          NM_row_prod_graph((size_t)n, (int)row, row, sum_sizes[color], sum_sizes[color + 1], M, permutation, z, &lefts[row], &right, false);

          // printf("AFTER ::: [%d] Color %d, lefts[%d] = %e, right = %e, w[%d] = %e\n", thread, color, row, lefts[row], right, row, w[row]);


          w[row] += right + z[row] / diag[row];

          // printf("BEFORE UPDATE ERR ::: w[%d] = %e\n", row, w[row]);

          #pragma omp atomic update
          err += pow(z[row] - fmax(0, (z[row] - w[row])), 2);

          /* If err < tol, do not execute this (= do not update z)
          So that we have w = Mz + q when we exit the loop */
          if (is_last == 0) {
          zi = lefts[row] + right;      

          DEBUG_PRINTF("diag[i] = %e\t zi = %e\n", diag[i], zi);
          zi = -(zi)*diag[row];

          if (zi < 0)
              z[row] = 0.0;
          else
              z[row] = zi;
          }

        }
        #pragma omp barrier
      }

      #pragma omp single
      {
        err = sqrt(err) / norm_q;
        iter += 1;
      }
      
      // printf("[%d] Iter = %d, Err = %e\n", thread, iter, err);

      if (verbose == 2) {
      printf(" # i%d -- %g : ", iter, err);
      for (i = 0; i < (size_t)n; ++i) printf(" %g", z[i]);
      for (i = 0; i < (size_t)n; ++i) printf(" %g", w[i]);
      printf("\n");
      }

    }
  }

  // time = omp_get_wtime() - time;
  // printf("lcp_pgs_graph_permut: time in while loop = %e s\n", time);

  options->iparam[SICONOS_IPARAM_ITER_DONE] = iter;
  options->dparam[SICONOS_DPARAM_RESIDU] = err;

  if (err > tol) {
    if (verbose > 0) {
      printf("Siconos/Numerics: lcp_pgs: No convergence of PGS after %d iterations\n", iter);
      printf("The residue is : %g \n", err);
    }
    *info = 1;
  } else {
    if (verbose > 0) {
      printf("Siconos/Numerics: lcp_pgs: Convergence of PGS after %d iterations\n", iter);
      printf("The residue is : %g \n", err);
    }
    *info = 0;
  }

  free(sum_sizes);
  free(inv_permutation);
  free(permutation);
  free(diag);
  free(lefts);
  free(start_lines);

}