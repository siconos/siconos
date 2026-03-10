/* Debug test for NSGS generic */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "FrictionContactProblem.h"
#include "FrictionContact_options.h"
#include "SolverOptions.h"
#include "fc3d_Solvers.h"

extern void fc3d_nsgs_generic(FrictionContactProblem* problem, double* reaction, double* velocity,
                              int* info, SolverOptions* options);

int main() {
  printf("NSGS Debug Test\n");
  printf("===============\n\n");
  
  /* Load problem */
  FrictionContactProblem* problem = frictionContact_new_from_filename("./data/FC3D_Example1.dat");
  if (!problem) {
    printf("Failed to load problem\n");
    return 1;
  }
  
  int nc = problem->numberOfContacts;
  int n = nc * 3;
  printf("Problem: %d contacts, %d variables\n", nc, n);
  printf("q vector: [%.4e, %.4e, %.4e, ...]\n", problem->q[0], problem->q[1], problem->q[2]);
  printf("mu: [%.4e, %.4e, %.4e]\n", problem->mu[0], problem->mu[1], problem->mu[2]);
  
  /* Setup options for generic */
  SolverOptions* options = solver_options_create(SICONOS_FRICTION_3D_NSGS);
  options->dparam[SICONOS_DPARAM_TOL] = 1e-4;
  options->iparam[SICONOS_IPARAM_MAX_ITER] = 10;  /* Just 10 iterations for debug */
  solver_options_update_internal(options, 0, SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone);
  options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-2;
  
  /* Allocate and zero initial guess */
  double* r = (double*)calloc(n, sizeof(double));
  double* v = (double*)calloc(n, sizeof(double));
  
  printf("\nInitial reaction: [%.4e, %.4e, %.4e, ...]\n", r[0], r[1], r[2]);
  
  /* Run generic solver */
  int info = -1;
  fc3d_nsgs_generic(problem, r, v, &info, options);
  
  printf("\nAfter %d iterations:\n", options->iparam[SICONOS_IPARAM_ITER_DONE]);
  printf("Final reaction: [%.4e, %.4e, %.4e, ...]\n", r[0], r[1], r[2]);
  printf("Final velocity: [%.4e, %.4e, %.4e, ...]\n", v[0], v[1], v[2]);
  printf("Info: %d, Error: %.4e\n", info, options->dparam[SICONOS_DPARAM_RESIDU]);
  
  /* Cleanup */
  free(r); free(v);
  solver_options_delete(options);
  frictionContactProblem_free(problem);
  return 0;
}
