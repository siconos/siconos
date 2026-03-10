/* Debug test 4 - check local residual inside the loop */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "FrictionContactProblem.h"
#include "FrictionContact_options.h"
#include "SolverOptions.h"
#include "fc3d_Solvers.h"

int main() {
  printf("NSGS Debug Test 4 - Check local residual inside loop\n");
  printf("=====================================================\n\n");
  
  /* Load problem */
  FrictionContactProblem* problem = frictionContact_new_from_filename("./data/FC3D_Example1.dat");
  if (!problem) {
    printf("Failed to load problem\n");
    return 1;
  }
  
  int nc = problem->numberOfContacts;
  printf("Problem: %d contacts\n", nc);
  
  /* Setup options */
  SolverOptions* options = solver_options_create(SICONOS_FRICTION_3D_NSGS);
  options->dparam[SICONOS_DPARAM_TOL] = 1e-4;
  solver_options_update_internal(options, 0, SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone);
  options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-2;
  
  printf("Internal solver options: %p\n", (void*)options->internalSolvers[0]);
  printf("Initial residual: %.4e\n", options->internalSolvers[0]->dparam[SICONOS_DPARAM_RESIDU]);
  
  /* Allocate */
  double* r = (double*)calloc(nc * 3, sizeof(double));
  double* v = (double*)calloc(nc * 3, sizeof(double));
  
  /* Run 1 iteration of fc3d_nsgs */
  options->iparam[SICONOS_IPARAM_MAX_ITER] = 1;
  
  printf("\nRunning fc3d_nsgs for 1 iteration...\n");
  int info = -1;
  fc3d_nsgs(problem, r, v, &info, options);
  
  printf("After fc3d_nsgs:\n");
  printf("  Residual: %.4e\n", options->internalSolvers[0]->dparam[SICONOS_DPARAM_RESIDU]);
  printf("  Reaction: [%.4e, %.4e, %.4e, ...]\n", r[0], r[1], r[2]);
  
  /* Cleanup */
  free(r); free(v);
  solver_options_delete(options);
  frictionContactProblem_free(problem);
  return 0;
}
