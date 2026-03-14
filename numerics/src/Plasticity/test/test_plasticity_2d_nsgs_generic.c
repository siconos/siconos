/* Test to compare legacy plasticity_2d_nsgs with plasticity_2d_nsgs_generic */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "PlasticityProblem.h"
#include "SolverOptions.h"
#include "plasticity_2d_solvers.h"
#include "Plasticity_options.h"
#include "numerics_verbose.h"
#include "SiconosBlas.h"


int main(int argc, char** argv) {
  (void)argc;
  (void)argv;
  
  printf("========================================\n");
  printf("Test: Legacy NSGS vs Generic NSGS\n");
  printf("========================================\n\n");
  
  /* Load test problem */
  const char* filename = "./data/plasticity_2d_example1.dat";
  PlasticityProblem* problem = plasticity2D_new_from_filename(filename);
  if (!problem) {
    printf("Failed to load problem from %s\n", filename);
    return 1;
  }
  
  printf("Problem loaded:\n");
  printf("  Number of cones: %d\n", problem->numberOfCones);
  printf("  Dimension: %d\n", problem->dimension);
  printf("  Model type: %s\n", 
         problem->model_type == PLASTICITY_MODEL_DRUCKER_PRAGER ? "Drucker-Prager" :
         problem->model_type == PLASTICITY_MODEL_VON_MISES ? "Von Mises" : "Unknown");
  
  int nc = problem->numberOfCones;
  int dim = problem->dimension;
  int n = nc * dim;
  
  /* Allocate solution vectors */
  double* stress_legacy = (double*)calloc(n, sizeof(double));
  double* strain_legacy = (double*)calloc(n, sizeof(double));
  double* stress_generic = (double*)calloc(n, sizeof(double));
  double* strain_generic = (double*)calloc(n, sizeof(double));
  
  /* ========================================
   * Test 1: Legacy NSGS
   * ======================================== */
  printf("\n--- Test 1: Legacy NSGS (PLASTICITY_2D_NSGS) ---\n");
  
  SolverOptions* options_legacy = solver_options_create(PLASTICITY_2D_NSGS);
  options_legacy->dparam[SICONOS_DPARAM_TOL] = 1e-6;
  options_legacy->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  
  /* Set internal solver to projection for consistency */
  solver_options_update_internal(options_legacy, 0, PLASTICITY_2D_ONECONE_ProjectionOnCone);
  
  int info_legacy;
  plasticity_2d_nsgs(problem, stress_legacy, strain_legacy, &info_legacy, options_legacy);
  
  printf("Legacy NSGS result:\n");
  printf("  Info: %d (%s)\n", info_legacy, info_legacy == 0 ? "converged" : "failed");
  printf("  Iterations: %d\n", options_legacy->iparam[SICONOS_IPARAM_ITER_DONE]);
  printf("  Residual: %e\n", options_legacy->dparam[SICONOS_DPARAM_RESIDU]);
  printf("  Stress norm: %e\n", cblas_dnrm2(n, stress_legacy, 1));
  
  /* ========================================
   * Test 2: Generic NSGS
   * ======================================== */
  printf("\n--- Test 2: Generic NSGS (PLASTICITY_2D_NSGS_GENERIC) ---\n");
  
  SolverOptions* options_generic = solver_options_create(PLASTICITY_2D_NSGS_GENERIC);
  options_generic->dparam[SICONOS_DPARAM_TOL] = 1e-6;
  options_generic->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  
  /* Set internal solver to projection for consistency */
  solver_options_update_internal(options_generic, 0, PLASTICITY_2D_ONECONE_ProjectionOnCone);
  
  int info_generic;
  plasticity_2d_nsgs_generic(problem, stress_generic, strain_generic, &info_generic, options_generic);
  
  printf("Generic NSGS result:\n");
  printf("  Info: %d (%s)\n", info_generic, info_generic == 0 ? "converged" : "failed");
  printf("  Iterations: %d\n", options_generic->iparam[SICONOS_IPARAM_ITER_DONE]);
  printf("  Residual: %e\n", options_generic->dparam[SICONOS_DPARAM_RESIDU]);
  printf("  Stress norm: %e\n", cblas_dnrm2(n, stress_generic, 1));
  
  /* ========================================
   * Compare results
   * ======================================== */
  printf("\n--- Comparison ---\n");
  
  double diff_norm = 0.0;
  for (int i = 0; i < n; i++) {
    double diff = stress_legacy[i] - stress_generic[i];
    diff_norm += diff * diff;
  }
  diff_norm = sqrt(diff_norm);
  
  printf("Difference between legacy and generic:\n");
  printf("  ||stress_legacy - stress_generic|| = %e\n", diff_norm);
  
  /* Relative difference */
  double legacy_norm = cblas_dnrm2(n, stress_legacy, 1);
  double rel_diff = (legacy_norm > 0) ? diff_norm / legacy_norm : diff_norm;
  printf("  Relative difference = %e\n", rel_diff);
  
  /* Print first few values */
  printf("\nFirst 9 stress values (legacy vs generic):\n");
  for (int i = 0; i < 9 && i < n; i++) {
    printf("  [%d] legacy=%12.8e  generic=%12.8e  diff=%12.8e\n", 
           i, stress_legacy[i], stress_generic[i], 
           fabs(stress_legacy[i] - stress_generic[i]));
  }
  
  /* ========================================
   * Cleanup
   * ======================================== */
  solver_options_delete(options_legacy);
  solver_options_delete(options_generic);
  plasticity2DProblem_free(problem);
  free(stress_legacy);
  free(strain_legacy);
  free(stress_generic);
  free(strain_generic);
  
  printf("\n========================================\n");
  printf("Test completed.\n");
  printf("========================================\n");
  
  return (info_legacy != 0 || info_generic != 0) ? 1 : 0;
}
