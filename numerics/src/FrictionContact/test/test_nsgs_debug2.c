/* Debug test 2 - trace accept_local */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "FrictionContactProblem.h"
#include "FrictionContact_options.h"
#include "SolverOptions.h"
#include "fc3d_Solvers.h"

/* Custom accept wrapper with debug output */
static void debug_accept_local(void* local_problem, SolverOptions* options,
                               unsigned int block, int iter,
                               double* var_z_global, double* var_z_local) {
  (void)local_problem;
  
  double local_residual = options->dparam[SICONOS_DPARAM_RESIDU];
  printf("  Accept block %d, iter %d: local=[%.4e, %.4e, %.4e], residual=%.4e ",
         block, iter, var_z_local[0], var_z_local[1], var_z_local[2], local_residual);
  
  if (isnan(local_residual) || isinf(local_residual) || local_residual > 1.0) {
    printf("[DISCARDED]\n");
    return;
  }
  
  var_z_global[block * 3 + 0] = var_z_local[0];
  var_z_global[block * 3 + 1] = var_z_local[1];
  var_z_global[block * 3 + 2] = var_z_local[2];
  printf("[ACCEPTED]\n");
}

int main() {
  printf("NSGS Debug Test 2 - Accept tracing\n");
  printf("===================================\n\n");
  
  /* Load problem */
  FrictionContactProblem* problem = frictionContact_new_from_filename("./data/FC3D_Example1.dat");
  if (!problem) {
    printf("Failed to load problem\n");
    return 1;
  }
  
  int nc = problem->numberOfContacts;
  int n = nc * 3;
  
  /* Setup options for generic - use original fc3d_nsgs first to compare */
  SolverOptions* options = solver_options_create(SICONOS_FRICTION_3D_NSGS);
  options->dparam[SICONOS_DPARAM_TOL] = 1e-4;
  options->iparam[SICONOS_IPARAM_MAX_ITER] = 3;  /* Just 3 iterations for debug */
  solver_options_update_internal(options, 0, SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone);
  options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-2;
  
  /* Allocate and zero initial guess */
  double* r = (double*)calloc(n, sizeof(double));
  double* v = (double*)calloc(n, sizeof(double));
  
  printf("\n=== Running ORIGINAL fc3d_nsgs ===\n");
  int info = -1;
  fc3d_nsgs(problem, r, v, &info, options);
  printf("After %d iterations: reaction=[%.4e, %.4e, %.4e, ...]\n",
         options->iparam[SICONOS_IPARAM_ITER_DONE], r[0], r[1], r[2]);
  
  /* Reset for generic */
  memset(r, 0, n * sizeof(double));
  memset(v, 0, n * sizeof(double));
  
  printf("\n=== Running GENERIC fc3d_nsgs_generic ===\n");
  options->iparam[SICONOS_IPARAM_ITER_DONE] = 0;
  fc3d_nsgs_generic(problem, r, v, &info, options);
  printf("After %d iterations: reaction=[%.4e, %.4e, %.4e, ...]\n",
         options->iparam[SICONOS_IPARAM_ITER_DONE], r[0], r[1], r[2]);
  
  /* Cleanup */
  free(r); free(v);
  solver_options_delete(options);
  frictionContactProblem_free(problem);
  return 0;
}
