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

void lcp_jacobi_parallel(LinearComplementarityProblem *problem, double *z, double *w, int *info,
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
  /* 2-norm of q, to normalize error */
  double norm_q = cblas_dnrm2(n, q, 1); 
  if (fabs(norm_q) <= DBL_EPSILON) norm_q = 1.;

  double *new_z = (double *)malloc((size_t)n * sizeof(double));

  /* double *true_w = (double *)malloc(n * sizeof(double));
  double true_err;
 */
  /* 
  Use a znew to store results 
  */

  /* Start solving */
  while ((iter < itermax) && (err > tol)) {

    err = 0.;

    // Compute new_z
    #pragma omp parallel for reduction(+:err) default(none) shared(q, z, new_z, w, diag, M, n)
    for (int i = 0; i < n; i++) {
      w[i] = q[i];
      DEBUG_PRINTF("zi = %e\n", new_z[i]);
      NM_row_prod(n, 1, i, M, z, &w[i], false);
      DEBUG_PRINTF("diag[i] = %e\t zi = %e\n", diag[i], new_z[i]);
      err += pow(z[i] - fmax(0, (z[i] - w[i])), 2);

      new_z[i] = w[i] - z[i] / diag[i];
      new_z[i] = -(new_z[i])*diag[i];

      if (new_z[i] < 0)
        new_z[i] = 0.0;
    }

    err = sqrt(err) / norm_q;

    // Update z if not last
    if (err > tol) {
      // Not sure if paralle update is faster than  
      #pragma omp parallel for default(none) shared(z, new_z, n)
      for (int i = 0; i < n; i++) {
        z[i] = new_z[i];
      }
    }

    iter += 1;

    if (verbose == 2) {
      printf(" # i%d -- %g : ", iter, err);
      for (size_t i = 0; i < (size_t)n; ++i) printf(" %g", z[i]);
      for (size_t i = 0; i < (size_t)n; ++i) printf(" %g", w[i]);
      printf("\n");
    }

    /* printf(" Iter %d, Err %g\n", iter, err);
    printf("z = ");
    for (size_t i = 0; i < (size_t)n; ++i) printf(" %g", z[i]);
    printf("\n");
    printf("w = ");
    for (size_t i = 0; i < (size_t)n; ++i) printf(" %g", w[i]);
    printf("\n"); */

    // lcp_compute_error(problem, z, true_w, tol, &true_err);

    /* printf(" Iter %d, True Err %g\n", iter, true_err);
    printf("z = ");
    for (size_t i = 0; i < (size_t)n; ++i) printf(" %g", z[i]);
    printf("\n");
    printf("w = ");
    for (size_t i = 0; i < (size_t)n; ++i) printf(" %g", true_w[i]);
    printf("\n");

    printf("\n"); */

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
