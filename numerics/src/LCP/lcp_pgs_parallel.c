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

#include <omp.h>

void lcp_pgs_parallel(LinearComplementarityProblem *problem, double *z, double *w, int *info,
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
  /* Initialize output */

  options->iparam[SICONOS_IPARAM_ITER_DONE] = 0;
  options->dparam[SICONOS_DPARAM_RESIDU] = 0.0;

  /* Preparation of the diagonal of the inverse matrix */
  double *diag = (double *)malloc(n * sizeof(double));
  NM_get_diag(n, info, M, diag);
  /* double diag_i = 0.0;
  for (int i = 0; i < n; ++i) {
    diag_i = NM_get_value(M, i, i);
    if (fabs(diag_i) < DBL_EPSILON) {
      if (verbose > 0) {
        printf("Numerics::lcp_pgs, error: vanishing diagonal term \n");
        printf(" The problem cannot be solved with this method \n");
      }

      *info = 2;
      free(diag);
      return;
    } else
      diag[i] = 1.0 / diag_i;
  } */

  /* Block of coordinates (i_g, j_g) will have size (block_sizes[i_g], block_sizes[j_g]) */
  int *block_sizes = NULL;
  /* Block of coordinates (i_g, j_g) starts at coordinates (start_indexes[i_g], start_indexes[j_g]) */
  int *start_indexes = NULL;
  /* Count the number of diagonal blocks computed (shared between threads)*/
  int counter = 0;
  /* Precise error (not light error) */
  double err = 0.;
  /* Number of OpenMP threads available */
  int g = omp_get_max_threads();

  /* Blocks will have size (n / g) or (n / g + 1) */
  block_sizes = (int *)malloc(g * sizeof(int));
  for (int i = 0; i < g; i++) {
    if (i < n % g) {
      block_sizes[i] = n / g + 1;
    } else {
      block_sizes[i] = n / g;
    }
  }

  start_indexes = (int *)malloc(g * sizeof(int));
  start_indexes[0] = 0;
  for (int i = 1; i < g; i++) {
    start_indexes[i] = start_indexes[i - 1] + block_sizes[i - 1];
  }

  /* 2-norm of q, to normalize error */
  double norm_q = cblas_dnrm2(n, problem->q, 1); 
  if (fabs(norm_q) <= DBL_EPSILON) norm_q = 1.;

  /* I have to call this once here, otherwise 
   * I get "double free or corruption" errors.
   * I think its because of parallel calls to
   * NM_block_prod which execute this in parallel.
   * This is not done in the graph version because
   * it is actually called once in color_graph_petsc
   * before the parallel section
   */
  CSparseMatrix* S = NULL;
  if (M->storageType == NM_SPARSE) {
    if (M->matrix2->origin == NSM_CSR) {
      S = NM_csr(M);
    } 
    else {
      S = NM_csc_trans(M);
    }
  }

  /* Flag to stop when convergence is achieved */
  int flag = 1;

  /* Start solving */
#pragma omp parallel default(none) shared(g, itermax, counter, err, flag, block_sizes, start_indexes, M, diag, tol, q, norm_q, w, z)
  {
    /* Thread number */
    int rank = omp_get_thread_num();
    /* Local block size (nb. of rows) */
    int size_i = block_sizes[rank];
    /* Line at which all blocks will start */
    int start_i = start_indexes[rank];

    /* We need two arrays to store the sums,
    in order to compute the precise error efficiently.    
    */

    /* Array to store sums, for blocks above the diagonal */
    double *t_right = NULL;
    t_right = (double *)malloc(size_i * sizeof(double));
    /* Array to store sums, for blocks below the diagonal */
    double *t_left = NULL;
    t_left = (double *)malloc(size_i * sizeof(double));

    /* Flag for last thread:
    when the last thread checks that err < tol,
    it needs to do one last iteration to synchronize
    with the other threads, so we need a specific flag for it
    */
    int flag_last_thread = 0; /* For all other threads, we don't need this flag */
    if (rank == g - 1) flag_last_thread = 1; /* Only last rhead needs it */

    /* Initialize left sums to 0 */
    for (int i = 0; i < size_i; i++) {
      t_left[i] = 0.;
    }

    /* To store z_i when diagonal block */
    double zsave = 0.;

    /* Number of iterations of current thread */
    int iter = 0;

    /* Main loop */
    while ((iter < itermax) && (flag || flag_last_thread)) {
      /* Initialize right sums with q */
      cblas_dcopy(size_i, q + start_i, 1, t_right, 1);

      /* Blocks after diagonal blocks */
      for (int j_g = rank + 1; j_g < g; j_g++) {
        /* Wait for diagonal block to be computed */
        while (counter <= (iter - 1) * g + j_g) {
        }

        NM_block_prod(start_i, start_indexes[j_g], size_i, block_sizes[j_g], M, z, t_right, 0);
      }

      /* Start computing w_{iter} */
      for (int i = start_i; i < start_i + size_i; i++) {
        w[i] = t_left[i - start_i] + z[i] / diag[i];
      }

      /* Initialize left sums to 0 */
      for (int i = 0; i < size_i; i++) {
        t_left[i] = 0.;
      }

      /* Blocks before diagonal blocks */
      for (int j_g = 0; j_g < rank; j_g++) {
        /* Wait for diagonal block to be computed */
        while (counter <= iter * g + j_g) {
        }

        NM_block_prod(start_i, start_indexes[j_g], size_i, block_sizes[j_g], M, z, t_left, 0);
      }

      /* Reset error when first group is back */
      if (rank == 0) {
        err = 0.;
      }

      /* Diagonal block */
      NM_block_prod_no_diag(start_i, size_i, M, z, t_right, &zsave, 0);

      if (flag == 0) {
        /* Last update of w and error, without 
        updating z so that w = Mz + q when we exit the while loop */
        for (int i = start_i; i < start_i + size_i; i++) {
          w[i] += t_right[i - start_i];
          err += pow(z[i] - fmax(0, (z[i] - w[i])), 2);
        }
      }
      else {
        /* Update z, w, and error */
        for (int i = start_i; i < start_i + size_i; i++) {

          /* Finish computing w_{iter} */
          w[i] += t_right[i - start_i];

          /* Update error between z_{iter} and w_{iter} */
          err += pow(z[i] - fmax(0, (z[i] - w[i])), 2);

          /* Compute z_{iter + 1} */
          z[i] = - (t_right[i - start_i] + t_left[i - start_i]) * diag[i];

          /* Local solver */
          if (z[i] < 0.) 
            z[i] = 0.;
        }
      }

      /* Last thread: z_{iter+1} is fully computed 
      check if error is good */
      if (rank == g - 1) {
        err = sqrt(err);
        err /= norm_q; /* Normalize error by norm of q */
        /* If flag is 0, the last thread already did its last iteration
        so we can update the second flag to stop this thread */
        if (flag == 0) flag_last_thread = 0;
        else if (err < tol) flag = 0; /* Set flag to 0 to stop all threads */
      }   

      iter += 1;

      #pragma omp atomic update
      counter += 1;
    }

    free(t_right);
    free(t_left);
  }

/* End of parallel section */

  options->iparam[SICONOS_IPARAM_ITER_DONE] = counter / g;
  options->dparam[SICONOS_DPARAM_RESIDU] = err;

  if (err > tol) {
    if (verbose > 0) {
      printf("Siconos/Numerics: lcp_pgs: No convergence of PGS after %d iterations\n", options->iparam[SICONOS_IPARAM_ITER_DONE]);
      printf("The residue is : %g \n", err);
    }
    *info = 1;
  } else {
    if (verbose > 0) {
      printf("Siconos/Numerics: lcp_pgs: Convergence of PGS after %d iterations\n", options->iparam[SICONOS_IPARAM_ITER_DONE]);
      printf("The residue is : %g \n", err);
    }
    *info = 0;
  }

  free(diag);
  free(block_sizes);
  free(start_indexes);

}