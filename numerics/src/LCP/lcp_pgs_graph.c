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

void lcp_pgs_graph(LinearComplementarityProblem *problem, double *z, double *w, int *info,
                   SolverOptions *options) {  
  
  // double total_time, time;
  // total_time = omp_get_wtime();

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

  // time = omp_get_wtime();
  /* Preparation of the diagonal of the inverse matrix */
  double *diag = (double *)malloc(n * sizeof(double));
  NM_get_invdiag(n, info, M, diag);
  // time = omp_get_wtime() - time;
  // printf("Time to prepare diagonal: %fs\n", time);

  /* Check if diagonal has a zero */
  if (*info == 2) {
    return;
  }

  /* Graph coloring */
  long int n_colors = 0;
  size_t *partition_size = NULL;
  size_t **partitions = NULL;
  
  color_graph(n, M, &n_colors, &partition_size, &partitions);  

  /* printf("colors = [");
  for (int i = 0; i < n; i++) {
    printf(" %d ", colors[i]);
  }
  printf("]\n"); */
  
  /* Solver variables */
  int iter = 0;
  double err = 1.;

  double zi;
  // double left = 0., right = 0.;
  int i;
  // double time = 0.;

  /* Start solving */
  while ((iter < itermax) && (err > tol)) {

    // time = omp_get_wtime();
    for (int color = 0; color < n_colors; color++) {

      #pragma omp parallel for default(none) private(zi, i) shared(q, z, diag, M, n, partitions, partition_size, color)
      for (int v = 0; v < partition_size[color]; v++) {
        i = partitions[color][v];
        zi = q[i];
        DEBUG_PRINTF("zi = %e\n", zi);
        // NM_row_prod_leftright(n, i, i, M, z, &left, &right, false);
        NM_row_prod_no_diag1x1(n, i, i, M, z, &zi, false);
        DEBUG_PRINTF("diag[i] = %e\t zi = %e\n", diag[i], zi);
        zi = -(zi)*diag[i];

        if (zi < 0)
          z[i] = 0.0;
        else
          z[i] = zi;

        /**
         * Some threads may modify z while others are using it in 
         * NM_row_prod_no_diag1x1 but it does not matter as 
         * each line are independent inside a color
        */
      }
    }
    // time = omp_get_wtime() - time;
    // printf("Time in inner loop = %f\n", time);

    // time = omp_get_wtime();
    lcp_compute_error(problem, z, w, tol, &err);
    // time = omp_get_wtime() - time;
    // printf("Time in lcp_compute_error = %f\n", time);

    if (verbose == 2) {
      printf(" # i%d -- %g : ", iter, err);
      for (i = 0; i < n; ++i) printf(" %g", z[i]);
      for (i = 0; i < n; ++i) printf(" %g", w[i]);
      printf("\n");
    }

    iter += 1;

  }

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

  // time = omp_get_wtime() - time;
  // printf("Time to solve: %fs\n", time);

  free(partition_size);
  for (int i = 0; i < n_colors; i++) free(partitions[i]);
  free(partitions);
  free(diag);

  // total_time = omp_get_wtime() - total_time;
  // printf("Total time: %fs\n", total_time);

}
