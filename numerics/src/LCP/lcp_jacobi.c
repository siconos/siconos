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

void lcp_jacobi(LinearComplementarityProblem *problem, double *z, double *w, int *info,
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
  double *diag = (double *)malloc((size_t)n * sizeof(double));
  NM_get_invdiag(n, info, M, diag);

  /* Check if diagonal has a zero */
  if (*info == 2) {
    return;
  }
  
  /* Solver variables */
  int iter = 0;
  double err = 1.;

  double *new_z = (double *)malloc((size_t)n * sizeof(double));

  /* Start solving */
  while ((iter < itermax) && (err > tol)) {

    // Compute new_z
    for (int i = 0; i < n; i++) {
      new_z[i] = q[i];
      DEBUG_PRINTF("zi = %e\n", new_z[i]);
      NM_row_prod_no_diag1x1((size_t)n, i, i, M, z, &new_z[i], false);
      DEBUG_PRINTF("diag[i] = %e\t zi = %e\n", diag[i], new_z[i]);

      new_z[i] = -(new_z[i])*diag[i];

      if (new_z[i] < 0)
        new_z[i] = 0.0;
    }

    cblas_dcopy(n, new_z, 1, z, 1);

    lcp_compute_error(problem, z, w, tol, &err);

    iter += 1;

    if ((iter % 100) == 0) printf(" %d %g \n", iter, err);

    if (verbose == 2) {
      printf(" # i%d -- %g : ", iter, err);
      for (size_t i = 0; i < (size_t)n; ++i) printf(" %g", z[i]);
      for (size_t i = 0; i < (size_t)n; ++i) printf(" %g", w[i]);
      printf("\n");
    }
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

  free(diag);
  free(new_z);

}