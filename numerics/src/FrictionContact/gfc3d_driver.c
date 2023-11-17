/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#include <assert.h>                        // for assert
#include <stdio.h>                         // for NULL, fprintf, stderr
#include <stdlib.h>                        // for exit, EXIT_FAILURE
#include "Friction_cst.h"                  // for SICONOS_GLOBAL_FRICTION_3D...
#include "GlobalFrictionContactProblem.h"  // for GlobalFrictionContactProblem
#include "NonSmoothDrivers.h"              // for gfc3d_driver
#include "NumericsFwd.h"                   // for SolverOptions, GlobalFrict...
#include "SolverOptions.h"                 // for SolverOptions, solver_opti...
#include "siconos_debug.h"                         // for DEBUG_EXPR
#include "gfc3d_Solvers.h"                 // for gfc3d_ACLMFixedPoint, gfc3...
#include "numerics_verbose.h"              // for numerics_printf_verbose
#include "gfc3d_balancing.h"
#include "gfc3d_compute_error.h"
#include "SiconosBlas.h"                         // for cblas_dcopy, cblas_dscal

#ifdef  DEBUG_MESSAGES
#include "NumericsVector.h"
#include "NumericsMatrix.h"
#endif

const char* const SICONOS_GLOBAL_FRICTION_3D_NSGS_WR_STR = "GFC3D_NSGS_WR";
const char* const SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR_STR = "GFC3D_NSN_AC_WR";
const char* const SICONOS_GLOBAL_FRICTION_3D_NSGSV_WR_STR = "GFC3D_NSGSV_WR";
const char* const SICONOS_GLOBAL_FRICTION_3D_PROX_WR_STR = "GFC3D_PROX_WR";
const char* const SICONOS_GLOBAL_FRICTION_3D_DSFP_WR_STR = "GFC3D_DSFP_WR";
const char* const SICONOS_GLOBAL_FRICTION_3D_TFP_WR_STR = "GFC3D_TFP_WR";
const char* const SICONOS_GLOBAL_FRICTION_3D_NSGS_STR = "GFC3D_NSGS";
const char* const SICONOS_GLOBAL_FRICTION_3D_NSN_AC_STR = "GFC3D_NSN_AC";
const char* const  SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH_STR = "GFC3D_GAMS_PATH";
const char* const  SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI_STR = "GFC3D_GAMS_PATHVI";
const char* const  SICONOS_GLOBAL_FRICTION_3D_VI_EG_STR = "GFC3D_VI_EG";
const char* const  SICONOS_GLOBAL_FRICTION_3D_ACLMFP_STR = "GFC3D_ACLMFP";
const char* const  SICONOS_GLOBAL_FRICTION_3D_VI_FPP_STR = "GFC3D_VI_FPP";
const char* const SICONOS_GLOBAL_FRICTION_3D_ADMM_WR_STR = "GFC3D_ADMM_WR";
const char* const SICONOS_GLOBAL_FRICTION_3D_IPM_WR_STR = "GFC3D_IPM_WR";
const char* const SICONOS_GLOBAL_FRICTION_3D_IPM_SNM_WR_STR = "GFC3D_IPM_SNM_WR";
const char* const SICONOS_GLOBAL_FRICTION_3D_IPM_SNM_SEP_STR = "GFC3D_IPM_SNM_SEP";

static int gfc3d_balancing_check_drift(GlobalFrictionContactProblem* balanced_problem,
                                       GlobalFrictionContactProblem* problem,
                                       double *reaction, double *velocity,
                                       double* globalVelocity,  SolverOptions* options)
{
  if(options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]>0)
  {
    size_t nc = problem->numberOfContacts;
    size_t n = problem->M->size0;
    size_t m = 3 * nc;

    double norm_b = cblas_dnrm2(m, balanced_problem->b, 1);
    double norm_q = cblas_dnrm2(n, balanced_problem->q, 1);
    double error_balancing = 1e24;
    double tolerance = options->dparam[SICONOS_DPARAM_TOL];
    gfc3d_compute_error(balanced_problem,  reaction, velocity, globalVelocity,
                        tolerance, options,
                        norm_q, norm_b,  &error_balancing);

    /* Come back to original variables */
    gfc3d_balancing_back_to_original_variables(balanced_problem, options,
                                               reaction, velocity, globalVelocity);

    norm_b = cblas_dnrm2(m, problem->b, 1);
    norm_q = cblas_dnrm2(n, problem->q, 1);
    double error =0.0;
    gfc3d_compute_error(problem,  reaction, velocity, globalVelocity,
                        tolerance, options,
                        norm_q, norm_b,  &error);

    numerics_printf_verbose(0,"error with balancing = %8.4e", error_balancing);
    numerics_printf_verbose(0,"error with original = %8.4e", error);
  }
  //else continue

  return 0;

}


int gfc3d_driver(GlobalFrictionContactProblem* problem, double *reaction, double *velocity,
                 double* globalVelocity,  SolverOptions* options, const char* problem_name)
{
  assert(options->isSet);
  DEBUG_EXPR(NV_display(globalVelocity,problem_ori->M->size0););
  if(verbose > 0)
    solver_options_print(options);

  /* Solver name */
  /*  const char* const  name = options->solverName;*/

  FILE *fileName = fopen("problem_name.res", "w");
  fprintf(fileName, "%s", problem_name);
  fclose(fileName);

  int info = -1 ;

  if(problem->dimension != 3)
    numerics_error("gfc3d_driver", "Dimension of the problem : problem-> dimension is not compatible or is not set");

  /* if there is no contact, we compute directly the global velocity as M^{-1}q */
  int m = problem->H->size1;
  if(m ==0)
  {
    numerics_printf_verbose(1,"---- GFC3D - DRIVER . No contact case. Direct computation of global velocity");
    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);
    return 0;
  }

  /* Non Smooth Gauss Seidel (NSGS) */
  switch(options->solverId)
  {
  case SICONOS_GLOBAL_FRICTION_3D_NSGS_WR:
  {

    numerics_printf_verbose(1," ========================== Call NSGS_WR solver with reformulation into Friction-Contact 3D problem ==========================\n");
    gfc3d_nsgs_wr(problem, reaction, velocity, globalVelocity, &info, options);
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_NSGSV_WR:
  {

    numerics_printf_verbose(1," ========================== Call NSGSV_WR solver with reformulation into Friction-Contact 3D problem ==========================\n");
    gfc3d_nsgs_velocity_wr(problem, reaction, velocity, globalVelocity, &info, options);
    break;
  }
  case SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR:
  {

    numerics_printf_verbose(1," ========================== Call NSN_AC_WR solver with reformulation into Friction-Contact 3D problem ==========================\n");
    gfc3d_nonsmooth_Newton_AlartCurnier_wr(problem, reaction, velocity, globalVelocity, &info, options);
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_PROX_WR:
  {

    numerics_printf_verbose(1," ========================== Call PROX_WR solver with reformulation into Friction-Contact 3D problem ==========================\n");
    gfc3d_proximal_wr(problem, reaction, velocity, globalVelocity, &info, options);
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_DSFP_WR:
  {

    numerics_printf_verbose(1," ========================== Call DSFP_WR solver with reformulation into Friction-Contact 3D problem ==========================\n");
    gfc3d_DeSaxceFixedPoint_wr(problem, reaction, velocity, globalVelocity, &info, options);
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_TFP_WR:
  {

    numerics_printf_verbose(1," ========================== Call TFP_WR solver with reformulation into Friction-Contact 3D problem ==========================\n");
    gfc3d_TrescaFixedPoint_wr(problem, reaction, velocity, globalVelocity, &info, options);
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_NSGS:
  {
    gfc3d_nsgs(problem, reaction, velocity, globalVelocity,
               &info, options);
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_NSN_AC:
  {
    /* Balancing */
    /* here, the balancing is done outside the solver */
    /* therfore the solver does not account for the possible drift of error measure between
       the balanced problem and the original one */

    GlobalFrictionContactProblem* balanced_problem = gfc3d_balancing_problem(problem,options);
    gfc3d_balancing_go_to_balanced_variables(balanced_problem, options,
                                             reaction, velocity, globalVelocity);
    /* Call the solver with balanced data */
    gfc3d_nonsmooth_Newton_AlartCurnier(balanced_problem, reaction, velocity,
                                        globalVelocity, &info, options);

    /* check if the drift is large */
    // int info_check_drift =
    gfc3d_balancing_check_drift(balanced_problem,problem, reaction, velocity, globalVelocity,
                                options);

    balanced_problem = gfc3d_balancing_free(balanced_problem, options);

    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH:
  {
    numerics_printf_verbose(1," ========================== Call PATH solver via GAMS for an AVI Friction-Contact 3D problem ==========================\n");
    gfc3d_AVI_gams_path(problem, reaction, velocity, &info, options);
    break;
  }
  case SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI:
  {
    numerics_printf_verbose(1," ========================== Call PATHVI solver via GAMS for an AVI Friction-Contact 3D problem ==========================\n");
    gfc3d_AVI_gams_pathvi(problem, reaction, globalVelocity, &info, options);
    break;
  }
  case SICONOS_GLOBAL_FRICTION_3D_VI_FPP:
  {
    gfc3d_VI_FixedPointProjection(problem, reaction, velocity,
                                  globalVelocity, &info, options);
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_VI_EG:
  {
    gfc3d_VI_ExtraGradient(problem, reaction, velocity,
                           globalVelocity, &info, options);
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_ACLMFP:
  {
    gfc3d_ACLMFixedPoint(problem, reaction, velocity,
                         globalVelocity, &info, options);
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_ADMM:
  {
    /* globalFrictionContact_rescaling(problem, 1.0/1.512808e-04, 1.0/1.407230e+01, 1.0); */
    gfc3d_ADMM(problem, reaction, velocity,
               globalVelocity, &info, options);
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_ADMM_WR:
  {

    numerics_printf_verbose(1," ========================== Call NSGS_WR solver with reformulation into Friction-Contact 3D problem ==========================\n");
    gfc3d_admm_wr(problem, reaction, velocity, globalVelocity, &info, options);
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_IPM:
  {

    GlobalFrictionContactProblem* balanced_problem = gfc3d_balancing_problem(problem,options);
    gfc3d_balancing_go_to_balanced_variables(balanced_problem, options,
                                             reaction, velocity, globalVelocity);

    gfc3d_IPM(balanced_problem, reaction, velocity,
              globalVelocity, &info, options, problem_name);


    gfc3d_balancing_check_drift(balanced_problem,problem, reaction, velocity, globalVelocity,
                                options);

    balanced_problem = gfc3d_balancing_free(balanced_problem, options);
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_IPM_WR:
  {

    gfc3d_ipm_wr(problem, reaction, velocity,
              globalVelocity, &info, options);
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_IPM_SNM:
  {
    gfc3d_IPM_SNM(problem, reaction, velocity,
              globalVelocity, &info, options, problem_name);
    break;
  }
  case SICONOS_GLOBAL_FRICTION_3D_IPM_SNM_WR:
  {
    gfc3d_ipm_snm_wr(problem, reaction, velocity,
              globalVelocity, &info, options);
    break;
  }
  case SICONOS_GLOBAL_FRICTION_3D_IPM_SNM_SEP:
  {
    #include <string.h>
    #include <time.h>
    #include "gfc3d_ipm.h"

    char *str = (char *) malloc(200);
    strcpy( str, problem_name );
    char * separators = "/";
    char *strToken_name = strtok( str, separators );
    for(int i=0; i<5; i++)
    {
      if(strToken_name != NULL) strToken_name = strtok ( NULL, separators );
    }
    strToken_name = strtok ( strToken_name, "." );

    char *line = NULL, *saveptr = NULL;
    int load_data = 0, len_prob_name = 0, n_blocks = 0.;
    size_t len = 0;
    FILE *varsep = fopen("varsep.res", "r");
    if (!varsep) printf("\n\n ERROR: varsep.res file is not available!!! \n\n");
    else
    {
      // Traverse the problem names in data file for a match
      for (int i=0; i<1091; i++)
      {
        if (getline(&line, &len, varsep))     // Read 1st line = problem name
        {
          len_prob_name = strlen(line);
          if (len_prob_name > 0 && line[len_prob_name - 1] == '\n')
          {
              line[len_prob_name - 1] = '\0';  // Replace the newline character with null terminator
          }
          if (strcmp(line, strToken_name) == 0) // Problem names are matched
          {
            load_data = 1;
            // printf("Matched: problem name = %s,", line);
            // break;
            // free(str);
          }
        }
        else
        {
          printf("ERROR: Error reading from varsep.res file.\n");
          break;
        }


        // Read No. of blocks of the current test
        if (getline(&line, &len, varsep))     // Read 2nd line = n_blocks
        {
          n_blocks = atoi(line);
          // printf("2nd line = %s\nn_blocks = %d\n", line, n_blocks);
        }

        // if (load_data) printf(", n_blocks = %d\n", n_blocks);

        if (load_data) break;

        // Go to the next problem name
        for (int i=0; i<n_blocks; i++)
          for (int j=0; j<4; j++) // each block has exaclty 4 lines
            getline(&line, &len, varsep);
      }
    }


    /*** Solve separate sub-problems ***/
    int blk_index = -1, count_contact_cp = 0, count_contact_full = 0, count_body_cp = 0, count_body_full = 0;
    int *permu_contact_cp = NULL, *permu_contact_full = NULL, *permu_body_cp = NULL, *permu_body_full = NULL;
    int Hc_m = -1, Hc_n = -1, Hc_rank = -1;
    ssize_t read;
    separators = " \t\n";

    // printf("BEFORE SOLVING:");
    // for(int i=0; i<problem->numberOfContacts; i++)
    // {
    //   printf("\n%2i-%2i: u = ", i*3, (i+1)*3-1);
    //   for (int j=0; j<3; j++) printf(" %.2e", velocity[i*3+j]);
    //   printf(",\tr = ");
    //   for (int j=0; j<3; j++) printf(" %.2e", reaction[i*3+j]);
    // }
    // printf("\n");

    // FILE *file = fopen("Hmat.m", "w");
    // fprintf(file,"H = [\n");
    // CSparseMatrix_print_in_Matlab_file(NM_triplet(NM_transpose(problem->H)), 0, file);
    // fprintf(file,"];\n");
    // fprintf(file,"H = sparse(int32(H(:,1)), int32(H(:,2)), H(:,3));\n");
    // fprintf(file,"M = [\n");
    // CSparseMatrix_print_in_Matlab_file(NM_triplet(problem->M), 0, file);
    // fprintf(file,"];\n");
    // fprintf(file,"M = sparse(int32(M(:,1)), int32(M(:,2)), M(:,3));\n");
    // fclose(file);

    // TESTING: change of friction coef.
    // for(int i = 0; i < problem->numberOfContacts ; i++) problem->mu[i]=0.3;

    // Timer
    long clk_tck = CLOCKS_PER_SEC;
    clock_t t1, t2;
    double time_sub_probs = 0, time_comple_prob = 0;

    // for getting pinfeas
    double mean_pinfeas = 0.;
    // options->solverData = calloc(1,sizeof(double));
    // double *mean_pinfeas = options->solverData;
    // *mean_pinfeas = 0.2;
    // *options->solverData = 0.2;

    int all_iterations = 0, n_succ = 0, n_fail = 0;
    for (int blk_i=0; blk_i<n_blocks; blk_i++)
    {
      // Block index
      getline(&line, &len, varsep);
      blk_index = atoi(line);

      // Permutation contact
      count_contact_cp = 0, count_contact_full = 0;
      read = getline(&line, &len, varsep);
      permu_contact_cp = malloc(read * sizeof(int)); // cp = compressed
      permu_contact_full = malloc(read * sizeof(int) * 3);
      char *strToken = strtok (line, separators);

      while (strToken != NULL)
      {
        // Convert the strToken to an integer and store it in permu_contact_cp array
        permu_contact_cp[count_contact_cp] = atoi(strToken)-1;
        for (int k=0; k<3; k++) permu_contact_full[count_contact_full++] = permu_contact_cp[count_contact_cp]*3+k;
        count_contact_cp++;

        // Get the next token
        strToken = strtok(NULL, separators);
      }

      // Permutation body
      count_body_cp = 0, count_body_full = 0;
      read = getline(&line, &len, varsep);
      permu_body_cp = malloc(read * sizeof(int));
      permu_body_full = malloc(read * sizeof(int) * 6);
      strToken = strtok (line, separators);
      while (strToken != NULL)
      {
        permu_body_cp[count_body_cp] = atoi(strToken)-1;
        for (int k=0; k<6; k++) permu_body_full[count_body_full++] = permu_body_cp[count_body_cp]*6+k;
        count_body_cp++;

        strToken = strtok(NULL, separators);
      }


      // Hc info
      getline(&line, &len, varsep);
      strToken = strtok (line, separators);
      strToken = strtok(NULL, separators);
      strToken = strtok(NULL, separators);
      Hc_rank = atoi(strToken);

      // // Print out for double-check
      // printf("Block %d:\npermutation contact compressed =", blk_index);
      // for (int i=0; i<count_contact_cp; i++) printf(" %d", permu_contact_cp[i]);
      // printf("\npermutation contact FULL =");
      // for (int i=0; i<count_contact_full; i++) printf(" %d", permu_contact_full[i]);
      // printf("\npermu_body compressed =");
      // for (int i=0; i<count_body_cp; i++) printf(" %d", permu_body_cp[i]);
      // printf("\npermu_body FULL =");
      // for (int i=0; i<count_body_full; i++) printf(" %d", permu_body_full[i]);
      // printf("\nH size nd x m = %d x %d\n\n", count_contact_cp*3, count_body_cp*6);

      // Create sub-problem
      GlobalFrictionContactProblem * sub_prob = globalFrictionContactProblem_new();
      sub_prob->dimension = problem->dimension;
      sub_prob->numberOfContacts = count_contact_cp;
      sub_prob->M = NM_extract(problem->M, count_body_full, permu_body_full, count_body_full, permu_body_full);
      sub_prob->env = NULL;

      // Attention: matrix H stored in hdf5 is transposed!
      sub_prob->H = NM_extract(problem->H, count_body_full, permu_body_full, count_contact_full, permu_contact_full);

      // vector f = problem->q
      sub_prob->q = (double*)malloc(count_body_full*sizeof(double));
      for (int i=0; i<count_body_full; i++) sub_prob->q[i] = problem->q[permu_body_full[i]];

      // vector w = problem->b
      sub_prob->b = (double*)malloc(count_contact_full*sizeof(double));
      for (int i=0; i<count_contact_full; i++) sub_prob->b[i] = problem->b[permu_contact_full[i]];

      // friction coef mu
      sub_prob->mu = (double*)malloc(count_contact_cp*sizeof(double));
      for (int i=0; i<count_contact_cp; i++) sub_prob->mu[i] = problem->mu[permu_contact_cp[i]];

      // Sub-solutions
      double *sub_reaction = (double*)calloc(count_contact_full, sizeof(double));
      double *sub_velocity = (double*)calloc(count_contact_full, sizeof(double));
      double *sub_globalVelocity = (double*)calloc(count_body_full, sizeof(double));

      printf("\n********** START sub-problem %d / %d, H (ndxm)/rank = %dx%d / %d *************************************************************\n", blk_i+1, n_blocks, count_contact_full, count_body_full, Hc_rank);
      t1 = clock();
      gfc3d_IPM_SNM(sub_prob, sub_reaction, sub_velocity,
              sub_globalVelocity, &info, options, problem_name);
      t2 = clock();
      time_sub_probs += (double)(t2-t1)/(double)clk_tck;
      all_iterations += options->iparam[SICONOS_IPARAM_ITER_DONE];
      double *pinfeas_ptr = options->solverData;
      mean_pinfeas += *pinfeas_ptr;
      free(options->solverData); options->solverData = NULL;
      if (info)
      {
        printf("test: failure\n");
        n_fail++;
      }
      else
      {
        printf("test: success\n");
        n_succ++;
      }
      printf("*************** END sub-problem %d / %d **************************************************************************************\n\n", blk_i+1, n_blocks);

      // Copy sub-solutions to main solutions
      for (int i=0; i<count_contact_full; i++)
      {
        reaction[permu_contact_full[i]] = sub_reaction[i];
        velocity[permu_contact_full[i]] = sub_velocity[i];
      }

      // for (int i=0; i<count_body_full; i++)
      // {
      //   globalVelocity[permu_body_full[i]] = sub_globalVelocity[i];
      // }


      free(sub_reaction); free(sub_velocity); free(sub_globalVelocity);
      globalFrictionContact_free(sub_prob);

      free(permu_contact_cp); free(permu_contact_full);
      free(permu_body_cp); free(permu_body_full);
    } // loop for n_blocks


    // printf("AFTER SOLVING:");
    // for(int i=0; i<problem->numberOfContacts; i++)
    // {
    //   printf("\n%2i-%2i: u = ", i*3, (i+1)*3-1);
    //   for (int j=0; j<3; j++) printf(" %.2e", velocity[i*3+j]);
    //   printf(",\tr = ");
    //   for (int j=0; j<3; j++) printf(" %.2e", reaction[i*3+j]);
    // }
    // printf("\n\n");


    // Double-check the result
    unsigned int nd = problem->H->size1;
    unsigned int d = problem->dimension;
    unsigned int n = problem->numberOfContacts;
    unsigned int m = problem->M->size0;

    NumericsMatrix *P_mu = NM_create(NM_SPARSE, nd, nd);
    NumericsMatrix *P_mu_inv = NM_create(NM_SPARSE, nd, nd);
    NM_triplet_alloc(P_mu, nd);
    NM_triplet_alloc(P_mu_inv, nd);
    P_mu->matrix2->origin = NSM_TRIPLET;
    P_mu_inv->matrix2->origin = NSM_TRIPLET;
    for(unsigned int i = 0; i < nd; ++i)
      if(i % d == 0)
      {
        NM_entry(P_mu, i, i, 1.);
        NM_entry(P_mu_inv, i, i, 1.);
      }
      else
      {
        NM_entry(P_mu, i, i, problem->mu[(int)(i/d)]);
        NM_entry(P_mu_inv, i, i, 1.0/problem->mu[(int)(i/d)]);
      }

    NumericsMatrix *M = problem->M;
    NumericsMatrix *H_tilde = NM_transpose(problem->H);
    double *w_tilde = problem->b;
    double *w = (double*)calloc(nd, sizeof(double));
    double *f = problem->q;

    double *velocity_tmp = (double*)calloc(nd, sizeof(double));
    double *reaction_tmp = (double*)calloc(nd, sizeof(double));
    double *primalConstraint = (double*)calloc(nd, sizeof(double));
    double *dualConstraint = (double*)calloc(m, sizeof(double));
    double *s = (double*)calloc(n, sizeof(double));

    double pinfeas, dinfeas, complem, udotr, nub, diff_u, diff_r;

    // Change of variable
    // Current velocity = [u0; mu*ub]
    NumericsMatrix *H = NM_multiply(P_mu, H_tilde);
    NM_gemv(1.0, P_mu, w_tilde, 0.0, w);


    globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);


    // Compute residuals
    // Current velocity = [u0; ub]
    // for (unsigned int i = 0; i<n; i++) s[i] = cblas_dnrm2(2, velocity+i*d+1, 1);
    // primalResidual_s(velocity, H, globalVelocity, w, s, primalConstraint, &pinfeas, options->dparam[SICONOS_DPARAM_TOL]);
    // dualResidual(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, options->dparam[SICONOS_DPARAM_TOL]);
    // complem = complemResidualNorm(velocity, reaction, nd, n);
    // udotr = cblas_ddot(nd, velocity, 1, reaction, 1);

    // Compute vars without friction coef.
    NM_gemv(1.0, P_mu, velocity, 0.0, velocity_tmp);
    NM_gemv(1.0, P_mu_inv, reaction, 0.0, reaction_tmp);
    // for (unsigned int i = 0; i<n; i++) s[i] = cblas_dnrm2(2, velocity_tmp+i*d+1, 1);
    // primalResidual_s(velocity_tmp, H, globalVelocity, w, s, primalConstraint, &pinfeas, options->dparam[SICONOS_DPARAM_TOL]);
    dualResidual(M, globalVelocity, H, reaction_tmp, f, dualConstraint, &dinfeas, options->dparam[SICONOS_DPARAM_TOL]);
    complem = complemResidualNorm(velocity_tmp, reaction_tmp, nd, n);
    udotr = cblas_ddot(nd, velocity_tmp, 1, reaction_tmp, 1);


    // Return to the original solutions to compute projection error
    // NM_gemv(1.0, P_mu_inv, velocity, 0.0, velocity_tmp);
    // NM_gemv(1.0, P_mu, reaction, 0.0, reaction_tmp);
    double norm_q = cblas_dnrm2(m, problem->q, 1);
    double norm_b = cblas_dnrm2(nd, problem->b, 1);
    double projerr, projerr_comple;
    // gfc3d_compute_error(problem, reaction_tmp, velocity_tmp, globalVelocity, options->dparam[SICONOS_DPARAM_TOL], options, norm_q, norm_b, &projerr);
    // gfc3d_compute_error(problem, reaction, velocity, globalVelocity, options->dparam[SICONOS_DPARAM_TOL], options, norm_q, norm_b, &projerr);
    projerr = projectionError(velocity_tmp, reaction_tmp, n, options->dparam[SICONOS_DPARAM_TOL]);


    // Store solutions obtained by solving sub-problems
    cblas_dcopy(nd, velocity, 1, velocity_tmp, 1);
    cblas_dcopy(nd, reaction, 1, reaction_tmp, 1);


    printf("\n\n============ Solve the complete problem by gfc3d_IPM_SNM ============\n");
    t1 = clock();
    gfc3d_IPM_SNM(problem, reaction, velocity, globalVelocity, &info, options, problem_name);
    t2 = clock();
    time_comple_prob = (double)(t2-t1)/(double)clk_tck;
    printf("Execution time: %.4f,\t\ttest: ", time_comple_prob);
    if (info) printf("failure\n"); else printf("success\n");

    // for comparison
    cblas_daxpy(nd, -1, velocity, 1, velocity_tmp, 1);
    cblas_daxpy(nd, -1, reaction, 1, reaction_tmp, 1);
    diff_u = cblas_dnrm2(nd, velocity_tmp, 1);
    diff_r = cblas_dnrm2(nd, reaction_tmp, 1);
    NM_gemv(1.0, P_mu, velocity, 0.0, velocity_tmp);
    NM_gemv(1.0, P_mu_inv, reaction, 0.0, reaction_tmp);
    projerr_comple = projectionError(velocity_tmp, reaction_tmp, n, options->dparam[SICONOS_DPARAM_TOL]);
    printf("=====================================================================\n\n");


    printf("====================== Result obtained by solving sub-problems ======================\n");
    numerics_printf_verbose(-1, "| its  | pinfeas | dinfeas | |u o r| |   u'r   | prj err |");
    numerics_printf_verbose(-1, "---------------------------------------------------------");
    numerics_printf_verbose(-1, "| %4i | %.1e | %.1e | %.1e | %.1e | %.1e |",
                          all_iterations, mean_pinfeas/n_blocks, dinfeas, complem, udotr, projerr);
    printf("Execution time: %.4f\n", time_sub_probs);
    printf("%d sub-problems, %4d success, %4d failure\n", n_blocks, n_succ, n_fail);
    printf("=====================================================================================\n\n");


    // Comparison

    printf("============ COMPARISON ============\n");
    printf("  | u_subs - u_comple | = %.2e\n", diff_u);
    printf("  | r_subs - r_comple | = %.2e\n", diff_r);
    printf("====================================\n\n\n");

    double min_val = mean_pinfeas/n_blocks;
    min_val = (min_val > dinfeas ? dinfeas : min_val) > complem ? complem : min_val;
    // min_val = min_val > dinfeas ? dinfeas : min_val;
    // min_val = min_val > complem ? dinfeas : min_val;
    printf("sumry total:   %s,\t  n nd x m = %d %d x %d\n",strToken_name, n, nd, m);
    printf("sumry total: %4d sub-prob, %4d succ, %4d fail,\t %.2e %.2e,\t %.4f (s),\t | u_subs - u_comple | = %.2e\n", n_blocks, n_succ, n_fail, min_val, projerr, time_sub_probs, diff_u);

    if (info)
      printf("sumry total:   comple prob,   failure           ,");
    else
      printf("sumry total:   comple prob,   success           ,");

    printf("\t %.2e %.2e,\t %.4f (s),\t | r_subs - r_comple | = %.2e\n",
            options->dparam[SICONOS_DPARAM_RESIDU], projerr_comple, time_comple_prob, diff_r);

    printf("sumry total: --------------------------------------------------------------------------------------------------------------------\n");


    if(H_tilde) {H_tilde = NM_free(H_tilde); H_tilde = NULL;}
    if(H) {H = NM_free(H); H = NULL;}
    if (w) {free(w); w = NULL;}
    if (s) {free(s); s = NULL;}
    if (velocity_tmp) {free(velocity_tmp); velocity_tmp = NULL;}
    if (reaction_tmp) {free(reaction_tmp); reaction_tmp = NULL;}
    if (primalConstraint) {free(primalConstraint); primalConstraint = NULL;}
    if (dualConstraint) {free(dualConstraint); dualConstraint = NULL;}
    free(str);
    fclose(varsep);
    break;
  }
  default:
  {
    fprintf(stderr, "Numerics, gfc3d_driver failed. Unknown solver %d.\n", options->solverId);
    exit(EXIT_FAILURE);

  }
  }

  return info;

}

int gfc3d_checkTrivialCaseGlobal(int n, double* q, double* velocity, double* reaction, double * globalVelocity, SolverOptions* options)
{
  /* norm of vector q */
  /*   double qs = cblas_dnrm2( n , q , 1 ); */
  /*   int i; */
  int info = -1;
  /*   if( qs <= DBL_EPSILON )  */
  /*     { */
  /*       // q norm equal to zero (less than DBL_EPSILON) */
  /*       // -> trivial solution: reaction = 0 and velocity = q */
  /*       for( i = 0 ; i < n ; ++i ) */
  /*  { */
  /*    velocity[i] = q[i]; */
  /*    reaction[i] = 0.; */
  /*  } */
  /*       iparam[2] = 0; */
  /*       iparam[4]= 0; */
  /*       dparam[1] = 0.0; */
  /*       dparam[3] = 0.0; */
  /*       info = 0; */
  /*       if(iparam[1]>0) */
  /*  printf("fc3d driver, norm(q) = 0, trivial solution reaction = 0, velocity = q.\n"); */
  /*     } */
  return info;
}
