/* Simple test to compare fc3d_nsgs and generic behavior */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "FrictionContactProblem.h"
#include "FrictionContact_options.h"
#include "SolverOptions.h"
#include "fc3d_Solvers.h"

int main() {
  printf("Simple NSGS Comparison Test\n");
  printf("===========================\n\n");
  
  /* Load problem */
  FrictionContactProblem* problem = frictionContact_new_from_filename("./data/FC3D_Example1.dat");
  if (!problem) {
    printf("Failed to load problem\n");
    return 1;
  }
  
  int nc = problem->numberOfContacts;
  int n = nc * 3;
  printf("Problem: FC3D_Example1.dat (%d contacts, %d variables)\n\n", nc, n);
  
  /* Test with different tolerances */
  double tolerances[] = {1e-4, 1e-6, 1e-8};
  int n_tols = sizeof(tolerances) / sizeof(tolerances[0]);
  
  for (int t = 0; t < n_tols; t++) {
    double tol = tolerances[t];
    printf("Tolerance: %.2e\n", tol);
    printf("-----------------\n");
    
    /* Setup options for original */
    SolverOptions* opts_orig = solver_options_create(SICONOS_FRICTION_3D_NSGS);
    opts_orig->dparam[SICONOS_DPARAM_TOL] = tol;
    opts_orig->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
    solver_options_update_internal(opts_orig, 0, SICONOS_ONECONE_ProjectionOnCone);
    opts_orig->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = tol * 100;
    
    /* Copy for generic */
    SolverOptions* opts_gen = solver_options_copy(opts_orig);
    
    /* Allocate and zero initial guess */
    double* r_orig = (double*)calloc(n, sizeof(double));
    double* v_orig = (double*)calloc(n, sizeof(double));
    double* r_gen = (double*)calloc(n, sizeof(double));
    double* v_gen = (double*)calloc(n, sizeof(double));
    
    /* Run original */
    int info_orig = -1;
    fc3d_nsgs(problem, r_orig, v_orig, &info_orig, opts_orig);
    
    /* Run generic */
    int info_gen = -1;
    fc3d_nsgs_generic(problem, r_gen, v_gen, &info_gen, opts_gen);
    
    /* Print results */
    printf("  Original:  iter=%3d, info=%d, error=%.4e\n",
           opts_orig->iparam[SICONOS_IPARAM_ITER_DONE],
           info_orig,
           opts_orig->dparam[SICONOS_DPARAM_RESIDU]);
    printf("  Generic:   iter=%3d, info=%d, error=%.4e\n",
           opts_gen->iparam[SICONOS_IPARAM_ITER_DONE],
           info_gen,
           opts_gen->dparam[SICONOS_DPARAM_RESIDU]);
    
    /* Compute solution difference */
    double diff = 0.0;
    for (int i = 0; i < n; i++) {
      double d = r_orig[i] - r_gen[i];
      diff += d * d;
    }
    diff = sqrt(diff);
    printf("  Solution diff: %.4e\n\n", diff);
    
    /* Cleanup */
    free(r_orig); free(v_orig);
    free(r_gen); free(v_gen);
    solver_options_delete(opts_orig);
    solver_options_delete(opts_gen);
  }
  
  frictionContactProblem_free(problem);
  return 0;
}
