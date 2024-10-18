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
  double diag_i = 0.0;
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
  }

  /* Block of coordinates (i_g, j_g) will have size (block_sizes[i_g], block_sizes[j_g]) */
  int *block_sizes = NULL;
  /* Block of coordinates (i_g, j_g) starts at coordinates (start_indexes[i_g], start_indexes[j_g]) */
  int *start_indexes = NULL;
  /* Count the number of diagonal blocks computed (shared between threads)*/
  int counter = 0;
  /* Error */
  double err = 1.;
  /* Number of threads available */
  int g;
  /* Error factor */
  double light_error_factor = 1.;

  /* Initialize variables */
  g = omp_get_max_threads();
  // printf("Number of threads available: %d\n", g);

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

  // Flag to stop when convergence achieved
  int flag = 1;

  // Start solving
#pragma omp parallel shared(counter, err, flag, block_sizes, start_indexes, M, diag, tol, z, light_error_factor)
  {
    // Get thread num
    int i_g = omp_get_thread_num();

    int m_i = block_sizes[i_g];  // Local block size (nb. of rows)
    int m_j;  // Local block size (nb. of columns)

    int start_i = start_indexes[i_g]; // Line at which all blocks will start

    /* Local array to store Gauss-Seidel sums */
    double *t = NULL;
    t = (double *)malloc(m_i * sizeof(double));

    /* To store zi and compute error */
    double zsave = 0.;

    /* Number of diagonal blocks computed by current thread */
    int iter = 0;

    /* Main loop */
    while ((iter < itermax) && (flag)) {
      /* Initialize sums with q */
      cblas_dcopy(m_i, q + start_i, 1, t, 1);

      /* Blocks after diagonal blocks */
      for (int j_g = i_g + 1; j_g < g; j_g++) {
        /* Wait for diagonal block to be computed */
        while (counter <= (iter - 1) * g + j_g) {
        }

        m_j = block_sizes[j_g];
        NM_block_prod(start_i, start_indexes[j_g], m_i, m_j, M, z, t, 0);
      }

      /* Blocks before diagonal blocks */
      for (int j_g = 0; j_g < i_g; j_g++) {
        /* Wait for diagonal block to be computed */
        while (counter <= iter * g + j_g) {
        }

        m_j = block_sizes[j_g];
        NM_block_prod(start_i, start_indexes[j_g], m_i, m_j, M, z, t, 0);
      }

      /* Reset error when first group is back */
      if (i_g == 0)
        err = 0.;

      /* Diagonal block */
      NM_block_prod_no_diag(start_i, m_i, M, z, t, &zsave, 0);

      /* Update z */
      for (int i = start_i; i < start_i + m_i; i++) {

        zsave = z[i];

        z[i] = - t[i - start_i] * diag[i];

        /* Local solver */
        if (z[i] < 0.)
          z[i] = 0.;

        /* Update light error */
        err += fabs(z[i] - zsave);
      }

      /* Last group: z_{k+1} is fully computed 
      check if error is good */
      if ((i_g == g - 1) && (light_error_factor * err < tol)) {
        /* Compute precise error */
        lcp_compute_error(problem, z, w, tol, &err);
        if (err < tol)
          flag = 0;
        else
          light_error_factor *= 1e3;
      }  

      iter += 1;

#pragma omp atomic update
      counter += 1;
    }

    free(t);   
  }
/* End of parallel section */

  /* Compute precise error */
  lcp_compute_error(problem, z, w, tol, &err);

  /* Number of iterations = (counter + 1) / g because:
  Threads nb 0 to g - 2 increase counter "iter" times
  Thread nb g - 1 increases counter "iter - 1" times
  */
  options->iparam[SICONOS_IPARAM_ITER_DONE] = (counter + 1) / g;
  options->dparam[SICONOS_DPARAM_RESIDU] = err;

  if (err > tol) {
    if (verbose > 0) {
      printf("Siconos/Numerics: lcp_pgs: No convergence of PGS after %d iterations\n", (counter + 1) / g);
      printf("The residue is : %g \n", err);
    }
    *info = 1;
  } else {
    if (verbose > 0) {
      printf("Siconos/Numerics: lcp_pgs: Convergence of PGS after %d iterations\n", (counter + 1) / g);
      printf("The residue is : %g \n", err);
    }
    *info = 0;
  }

  free(diag);
  free(block_sizes);
  free(start_indexes);

}