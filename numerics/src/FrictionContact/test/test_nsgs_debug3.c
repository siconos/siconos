/* Debug test 3 - check local residual value */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "FrictionContactProblem.h"
#include "FrictionContact_options.h"
#include "SolverOptions.h"
#include "fc3d_Solvers.h"

extern void fc3d_nsgs_update(int contact, FrictionContactProblem* problem,
                             FrictionContactProblem* localproblem, double* reaction,
                             SolverOptions* options);

int main() {
  printf("NSGS Debug Test 3 - Check local residual\n");
  printf("=========================================\n\n");
  
  /* Load problem */
  FrictionContactProblem* problem = frictionContact_new_from_filename("./data/FC3D_Example1.dat");
  if (!problem) {
    printf("Failed to load problem\n");
    return 1;
  }
  
  int nc = problem->numberOfContacts;
  
  /* Setup options */
  SolverOptions* options = solver_options_create(SICONOS_FRICTION_3D_NSGS);
  options->dparam[SICONOS_DPARAM_TOL] = 1e-4;
  solver_options_update_internal(options, 0, SICONOS_ONECONE_ProjectionOnCone);
  options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-2;
  
  /* Set residual to 0 initially */
  options->internalSolvers[0]->dparam[SICONOS_DPARAM_RESIDU] = 0.0;
  
  /* Allocate */
  double* r = (double*)calloc(nc * 3, sizeof(double));
  double* v = (double*)calloc(nc * 3, sizeof(double));
  
  /* Run fc3d_nsgs for 1 iteration */
  options->iparam[SICONOS_IPARAM_MAX_ITER] = 1;
  int info = -1;
  fc3d_nsgs(problem, r, v, &info, options);
  
  printf("After 1 iteration of fc3d_nsgs:\n");
  printf("  Local solver residual: %.4e\n", options->internalSolvers[0]->dparam[SICONOS_DPARAM_RESIDU]);
  printf("  Reaction: [%.4e, %.4e, %.4e, ...]\n", r[0], r[1], r[2]);
  
  /* Now try with fresh options */
  SolverOptions* options2 = solver_options_create(SICONOS_FRICTION_3D_NSGS);
  options2->dparam[SICONOS_DPARAM_TOL] = 1e-4;
  solver_options_update_internal(options2, 0, SICONOS_ONECONE_ProjectionOnCone);
  options2->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-2;
  
  /* Set residual to NaN initially */
  options2->internalSolvers[0]->dparam[SICONOS_DPARAM_RESIDU] = NAN;
  
  double* r2 = (double*)calloc(nc * 3, sizeof(double));
  double* v2 = (double*)calloc(nc * 3, sizeof(double));
  
  options2->iparam[SICONOS_IPARAM_MAX_ITER] = 1;
  info = -1;
  fc3d_nsgs(problem, r2, v2, &info, options2);
  
  printf("\nAfter 1 iteration with NaN initial residual:\n");
  printf("  Local solver residual: %.4e\n", options2->internalSolvers[0]->dparam[SICONOS_DPARAM_RESIDU]);
  printf("  Reaction: [%.4e, %.4e, %.4e, ...]\n", r2[0], r2[1], r2[2]);
  
  /* Cleanup */
  free(r); free(v);
  free(r2); free(v2);
  solver_options_delete(options);
  solver_options_delete(options2);
  frictionContactProblem_free(problem);
  return 0;
}
