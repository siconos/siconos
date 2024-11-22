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

#include "graph_tools_petsc.h"

void lcp_pgs_graph_petsc(LinearComplementarityProblem *problem, double *z, double *w, int *info,
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

  /* Time */
  double time;

  /* Graph coloring */
  long int n_colors = 0;
  size_t *partition_size = NULL;
  size_t **partitions = NULL;
  
  time = omp_get_wtime();
  color_graph_petsc(n, M, &n_colors, &partition_size, &partitions);  
  time = omp_get_wtime() - time;
  printf("Time to color graph: %fs\n", time);

  /* printf("colors = [");
  for (int i = 0; i < n; i++) {
    printf(" %d ", colors[i]);
  }
  printf("]\n"); */
  
  /* Solver variables */
  int iter = 0;
  double err = 1.;

  double zi;
  int i;

  time = omp_get_wtime();

  /* Start solving */
  while ((iter < itermax) && (err > tol)) {

    for (int color = 0; color < n_colors; color++) {

      #pragma omp parallel for default(none) private(zi, i) shared(q, z, diag, M, n, partitions, partition_size, color)
      for (int v = 0; v < partition_size[color]; v++) {
        i = partitions[color][v];
        zi = q[i];
        DEBUG_PRINTF("zi = %e\n", zi);
        NM_row_prod_no_diag1x1(n, i, i, M, z, &zi, false);
        DEBUG_PRINTF("diag[i] = %e\t zi = %e\n", diag[i], zi);
        zi = -(zi)*diag[i];

        if (zi < 0)
          z[i] = 0.0;
        else
          z[i] = zi;

      }
    }

    lcp_compute_error(problem, z, w, tol, &err);

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

  time = omp_get_wtime() - time;
  printf("Time to solve: %fs\n", time);

  free(partition_size);
  for (int i = 0; i < n_colors; i++) free(partitions[i]);
  free(partitions);
  free(diag);

}
