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
#include "graph_tools.h"
#include "lcp_cst.h"
#include "numerics_verbose.h"  // for verbose
#include "siconos_debug.h"

/* Solver registration system */
#include "numerics_errors.h"
#include "solver_registry.h"

void lcp_pgs_graph(LinearComplementarityProblem* problem, double* z, double* w, int* info,
                   SolverOptions* options) {
  NumericsMatrix* M = problem->M;
  double* q = problem->q;

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
  double* diag = (double*)malloc((size_t)n * sizeof(double));
  NM_get_invdiag(n, info, M, diag);

  /* Check if diagonal has a zero */
  if (*info == 2) {
    return;
  }

  /* Graph coloring */
  int err_not_symmetric;
  size_t n_colors = 0;
  size_t* partition_size = NULL;
  size_t** partitions = NULL;

  err_not_symmetric = color_graph(n, M, &n_colors, &partition_size, &partitions);

  /* Matrix not symmetric */
  if (err_not_symmetric == 1) {
    // printf("ERROR: matrix is not symmetric.\n");
    free(diag);
    return;
  }

  /* Solver variables */
  int iter = 0;
  double err = 1.;

  double zi;
  size_t i;

  /* Start solving */
  while ((iter < itermax) && (err > tol)) {
    for (size_t color = 0; color < n_colors; color++) {
#pragma omp parallel for default(none) private(zi, i) \
    shared(q, z, diag, M, n, partitions, partition_size, color)
      for (int v = 0; v < partition_size[color]; v++) {
        i = partitions[color][v];
        zi = q[i];
        DEBUG_PRINTF("zi = %e\n", zi);
        NM_row_prod_no_diag1x1((size_t)n, (int)i, i, M, z, &zi, false);
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
         * AS LONG AS M IS SYMMETRIC
         */
      }
    }

    lcp_compute_error(problem, z, w, tol, &err);

    if (verbose == 2) {
      printf(" # i%d -- %g : ", iter, err);
      for (int j = 0; j < n; ++j) printf(" %g", z[j]);
      for (int j = 0; j < n; ++j) printf(" %g", w[j]);
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

  free(partition_size);
  for (size_t i = 0; i < n_colors; i++) free(partitions[i]);
  free(partitions);
  free(diag);
}

static void lcp_pgs_set_default(SolverOptions* options) {
  /* No specific defaults needed */
  (void)options;
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 * This registers SICONOS_LCP_PGS in the global solver registry.
 */

static int lcp_pgs_init_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  lcp_pgs_set_default(options);
  return NUMERICS_OK;
}

static int lcp_pgs_solve_wrap(void* problem, double* z, double* w, SolverOptions* options) {
  int info = NUMERICS_OK;
  lcp_pgs_graph((LinearComplementarityProblem*)problem, z, w, &info, options);
  return info;
}

static void lcp_pgs_free_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
}

REGISTER_SOLVER(SICONOS_LCP_PGS_GRAPH, "SICONOS_LCP_PGS_GRAPH", "SICONOS_LCP_PGS_GRAPH",
                lcp_pgs_init_wrap, lcp_pgs_solve_wrap, lcp_pgs_free_wrap,
                NULL,                /* error function */
                lcp_pgs_set_default, /* set_default */
                1000,                /* default_max_iter */
                1e-6,                /* default_tol */
                0)                   /* is_local_solver */
