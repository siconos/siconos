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
#include "fc3d_Solvers.h"                   // for fc3d_nsgs_set_default


#include <string.h>
#include <time.h>
#include "gfc3d_ipm.h"
#include <math.h>

#ifdef  DEBUG_MESSAGES
#include "NumericsVector.h"
#include "NumericsMatrix.h"
#endif

const char* const SICONOS_GLOBAL_FRICTION_3D_NSGS_WR_STR = "GFC3D_NSGS_WR";
const char* const SICONOS_GLOBAL_FRICTION_3D_NSGS_SEP_WR_STR = "GFC3D_NSGS_SEP_WR";
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
const char* const   SICONOS_GLOBAL_FRICTION_3D_IPM_SNM_PROX_STR = "GFC3D IPM SNM PROX";


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
                 double* globalVelocity,  SolverOptions* options)
{
  assert(options->isSet);
  DEBUG_EXPR(NV_display(globalVelocity,problem_ori->M->size0););
  if(verbose > 0)
    solver_options_print(options);

  /* Solver name */
  /*  const char* const  name = options->solverName;*/

  if (problem->name)
  {
    FILE *fileName = fopen("problem_name.res", "w");
    char * problem_name = problem->name;
    fprintf(fileName, "%s", problem_name);
    fclose(fileName);
  }
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
  case SICONOS_GLOBAL_FRICTION_3D_NSGS_SEP_WR:
  {
    numerics_printf_verbose(1," ========================== Call NSGS_SEP_WR solver with reformulation into Friction-Contact 3D problem ==========================\n");
    gfc3d_nsgs_wr(problem, reaction, velocity, globalVelocity, &info, options);

    // int nd = 33, m=30;
    // FILE *sol_NSGS = fopen("solNSGS.m", "w");
    // fprintf(sol_NSGS,"uNSGS = [");
    // for(int i = 0; i < nd; i++)
    // {
    //   fprintf(sol_NSGS, "%.20e, ", velocity[i]);
    // }
    // fprintf(sol_NSGS,"];\n");
    // fprintf(sol_NSGS,"rNSGS = [");
    // for(int i = 0; i < nd; i++)
    // {
    //   fprintf(sol_NSGS, "%.20e, ", reaction[i]);
    // }
    // fprintf(sol_NSGS,"];\n");
    // fprintf(sol_NSGS,"vNSGS = [");
    // for(int i = 0; i < m; i++)
    // {
    //   fprintf(sol_NSGS, "%.20e, ", globalVelocity[i]);
    // }
    // fprintf(sol_NSGS,"];\n");
    // fclose(sol_NSGS);

    // char *str = (char *) malloc(200);
    // if (problem->name)
    // {
    //   strcpy( str, problem->name );
    // }
    // else
    // {
    //   strcpy( str, "foo_" );
    // }
    // char * separators = "/";
    // char *strToken_name = strtok( str, separators );
    // for(int i=0; i<5; i++)
    // {
    //   if(strToken_name != NULL) strToken_name = strtok ( NULL, separators );
    // }
    // strToken_name = strtok ( strToken_name, "." );

    // char *line = NULL, *saveptr = NULL;
    // int load_data = 0, len_prob_name = 0, n_blocks = 0.;
    // size_t len = 0;
    // FILE *varsep = fopen("varsep.res", "r");
    // if (!varsep) printf("\n\n ERROR: varsep.res file is not available!!! \n\n");
    // else
    // {
    //   // Traverse the problem names in data file for a match
    //   for (int i=0; i<1081; i++)
    //   {
    //     if (getline(&line, &len, varsep))     // Read 1st line = problem name
    //     {
    //       len_prob_name = strlen(line);
    //       if (len_prob_name > 0 && line[len_prob_name - 1] == '\n')
    //       {
    //           line[len_prob_name - 1] = '\0';  // Replace the newline character with null terminator
    //       }
    //       if (strcmp(line, strToken_name) == 0) // Problem names are matched
    //       {
    //         load_data = 1;
    //       }
    //     }
    //     else
    //     {
    //       printf("ERROR: Error reading from varsep.res file.\n");
    //       break;
    //     }


    //     // Read No. of blocks of the current test
    //     if (getline(&line, &len, varsep))     // Read 2nd line = n_blocks
    //     {
    //       n_blocks = atoi(line);
    //     }


    //     if (load_data) break;

    //     // Go to the next problem name
    //     for (int i=0; i<n_blocks; i++)
    //       for (int j=0; j<4; j++) // each block has exaclty 4 lines
    //         getline(&line, &len, varsep);
    //   }
    // }



    // SolverOptions *options_ipm = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_IPM_SNM_SEP);
    // gfc3d_ipm_snm_set_default(options_ipm);



    // /*** Solve separate sub-problems ***/
    // int blk_index = -1, count_contact_cp = 0, count_contact_full = 0, count_body_cp = 0, count_body_full = 0;
    // int *permu_contact_cp = NULL, *permu_contact_full = NULL, *permu_body_cp = NULL, *permu_body_full = NULL;
    // int Hc_m = -1, Hc_n = -1, Hc_rank = -1;
    // ssize_t read;
    // separators = " \t\n";

    // for (int blk_i=0; blk_i<n_blocks; blk_i++)
    // {
    //   // Block index
    //   getline(&line, &len, varsep);
    //   blk_index = atoi(line);

    //   // Permutation contact
    //   count_contact_cp = 0, count_contact_full = 0;
    //   read = getline(&line, &len, varsep);
    //   permu_contact_cp = malloc(read * sizeof(int)); // cp = compressed
    //   permu_contact_full = malloc(read * sizeof(int) * 3);
    //   char *strToken = strtok (line, separators);

    //   while (strToken != NULL)
    //   {
    //     // Convert the strToken to an integer and store it in permu_contact_cp array
    //     permu_contact_cp[count_contact_cp] = atoi(strToken)-1;
    //     for (int k=0; k<3; k++) permu_contact_full[count_contact_full++] = permu_contact_cp[count_contact_cp]*3+k;
    //     count_contact_cp++;

    //     // Get the next token
    //     strToken = strtok(NULL, separators);
    //   }

    //   // Permutation body
    //   count_body_cp = 0, count_body_full = 0;
    //   read = getline(&line, &len, varsep);
    //   permu_body_cp = malloc(read * sizeof(int));
    //   permu_body_full = malloc(read * sizeof(int) * 6);
    //   strToken = strtok (line, separators);
    //   while (strToken != NULL)
    //   {
    //     permu_body_cp[count_body_cp] = atoi(strToken)-1;
    //     for (int k=0; k<6; k++) permu_body_full[count_body_full++] = permu_body_cp[count_body_cp]*6+k;
    //     count_body_cp++;

    //     strToken = strtok(NULL, separators);
    //   }


    //   // Hc info
    //   getline(&line, &len, varsep);
    //   strToken = strtok (line, separators);
    //   strToken = strtok(NULL, separators);
    //   strToken = strtok(NULL, separators);
    //   Hc_rank = atoi(strToken);


    //   // Create sub-problem
    //   GlobalFrictionContactProblem * sub_prob = globalFrictionContactProblem_new();
    //   sub_prob->dimension = problem->dimension;
    //   sub_prob->numberOfContacts = count_contact_cp;
    //   sub_prob->M = NM_extract(problem->M, count_body_full, permu_body_full, count_body_full, permu_body_full);
    //   sub_prob->env = NULL;

    //   /* symmetrization of the matrix M */
    //   if(!(NM_is_symmetric(sub_prob->M)))
    //   {
    //     printf("#################### SYMMETRIZATION ####################\n");
    //     NumericsMatrix *MT = NM_transpose(sub_prob->M);
    //     NumericsMatrix * MSym = NM_add(1 / 2., sub_prob->M, 1 / 2., MT);
    //     NM_free(sub_prob->M);
    //     sub_prob->M = MSym;
    //     NM_free(MT);
    //   }

    //   // Attention: matrix H stored in hdf5 is transposed!
    //   sub_prob->H = NM_extract(problem->H, count_body_full, permu_body_full, count_contact_full, permu_contact_full);

    //   // vector f = problem->q
    //   sub_prob->q = (double*)malloc(count_body_full*sizeof(double));
    //   for (int i=0; i<count_body_full; i++) sub_prob->q[i] = problem->q[permu_body_full[i]];

    //   // vector w = problem->b
    //   sub_prob->b = (double*)malloc(count_contact_full*sizeof(double));
    //   for (int i=0; i<count_contact_full; i++) sub_prob->b[i] = problem->b[permu_contact_full[i]];

    //   // friction coef mu
    //   sub_prob->mu = (double*)malloc(count_contact_cp*sizeof(double));
    //   for (int i=0; i<count_contact_cp; i++) sub_prob->mu[i] = problem->mu[permu_contact_cp[i]];

    //   // name
    //   sub_prob->name = problem->name;

    //   // Sub-solutions
    //   double *sub_reaction = (double*)calloc(count_contact_full, sizeof(double));
    //   double *sub_velocity = (double*)calloc(count_contact_full, sizeof(double));
    //   double *sub_globalVelocity = (double*)calloc(count_body_full, sizeof(double));
    //   double *sub_reaction_tmp = (double*)calloc(count_contact_full, sizeof(double));
    //   double *sub_velocity_tmp = (double*)calloc(count_contact_full, sizeof(double));


    //   // if (options) solver_options_delete(options); free(options); options = NULL;
    //   // options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_NSGS_SEP_WR);
    //   // fc3d_nsgs_set_default(options);
    //   // options->dparam[SICONOS_DPARAM_TOL] = 1e-3;
    //   // options->iparam[SICONOS_IPARAM_MAX_ITER] = 100;

    //   // Save the block number to append in the test name
    //   options->solverData = (char *)malloc(10*sizeof(char));
    //   char *blk_num_name = (char *)options->solverData;
    //   sprintf(blk_num_name, "-%d", blk_i);

    //   printf("\n********** START sub-problem %d / %d, H (ndxm)/rank = %dx%d / %d *************************************************************\n", blk_i+1, n_blocks, count_contact_full, count_body_full, Hc_rank);
    //   info = -1;
    //   gfc3d_nsgs_wr(sub_prob, sub_reaction, sub_velocity, sub_globalVelocity, &info, options);

    //   if (info)
    //   {
    //     printf("test: failure\n");
    //   }
    //   else
    //   {
    //     printf("test: success\n");
    //   }
    //   printf("*************** END sub-problem %d / %d **************************************************************************************\n\n", blk_i+1, n_blocks);






    //   // sub_prob->M = NM_extract(problem->M, count_body_full, permu_body_full, count_body_full, permu_body_full);
    //   // sub_prob->env = NULL;

    //   // /* symmetrization of the matrix M */
    //   // if(!(NM_is_symmetric(sub_prob->M)))
    //   // {
    //   //   printf("#################### SYMMETRIZATION ####################\n");
    //   //   NumericsMatrix *MT = NM_transpose(sub_prob->M);
    //   //   NumericsMatrix * MSym = NM_add(1 / 2., sub_prob->M, 1 / 2., MT);
    //   //   NM_free(sub_prob->M);
    //   //   sub_prob->M = MSym;
    //   //   NM_free(MT);
    //   // }

    //   // // Attention: matrix H stored in hdf5 is transposed!
    //   // sub_prob->H = NM_extract(problem->H, count_body_full, permu_body_full, count_contact_full, permu_contact_full);

    //   // // vector f = problem->q
    //   // sub_prob->q = (double*)malloc(count_body_full*sizeof(double));
    //   // for (int i=0; i<count_body_full; i++) sub_prob->q[i] = problem->q[permu_body_full[i]];

    //   // // vector w = problem->b
    //   // sub_prob->b = (double*)malloc(count_contact_full*sizeof(double));
    //   // for (int i=0; i<count_contact_full; i++) sub_prob->b[i] = problem->b[permu_contact_full[i]];

    //   // // friction coef mu
    //   // sub_prob->mu = (double*)malloc(count_contact_cp*sizeof(double));
    //   // for (int i=0; i<count_contact_cp; i++) sub_prob->mu[i] = problem->mu[permu_contact_cp[i]];

    //   // // name
    //   // sub_prob->name = problem->name;


    //   // if (options) solver_options_delete(options); free(options); options = NULL;
    //   // options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_IPM_SNM_SEP);
    //   // gfc3d_ipm_snm_set_default(options);
    //   // Save the block number to append in the test name
    //   options_ipm->solverData = (char *)malloc(10*sizeof(char));
    //   blk_num_name = (char *)options_ipm->solverData;
    //   sprintf(blk_num_name, "-%d", blk_i);
    //   gfc3d_IPM_SNM(sub_prob, sub_reaction, sub_velocity, sub_globalVelocity, &info, options_ipm);






    //   // Copy sub-solutions to main solutions
    //   for (int i=0; i<count_contact_full; i++)
    //   {
    //     reaction[permu_contact_full[i]] = sub_reaction[i];
    //     velocity[permu_contact_full[i]] = sub_velocity[i];
    //   }

    //   // for (int i=0; i<count_body_full; i++)
    //   // {
    //   //   globalVelocity[permu_body_full[i]] = sub_globalVelocity[i];
    //   // }


    //   free(sub_reaction); free(sub_velocity); free(sub_globalVelocity);
    //   free(sub_reaction_tmp); free(sub_velocity_tmp);
    //   globalFrictionContact_free(sub_prob);

    //   free(permu_contact_cp); free(permu_contact_full);
    //   free(permu_body_cp); free(permu_body_full);
    // } // loop for n_blocks

    // fclose(varsep);
    // free(str);

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
  case SICONOS_GLOBAL_FRICTION_3D_ADMM_SEP:
  {
    numerics_printf_verbose(1," ========================== Call ADMM_SEP solver with reformulation into Friction-Contact 3D problem ==========================\n");
    gfc3d_ADMM(problem, reaction, velocity,
               globalVelocity, &info, options);
    break;


    // char *str = (char *) malloc(200);
    // if (problem->name)
    // {
    //   strcpy( str, problem->name );
    // }
    // else
    // {
    //   strcpy( str, "foo_" );
    // }
    // char * separators = "/";
    // char *strToken_name = strtok( str, separators );
    // for(int i=0; i<5; i++)
    // {
    //   if(strToken_name != NULL) strToken_name = strtok ( NULL, separators );
    // }
    // strToken_name = strtok ( strToken_name, "." );

    // char *line = NULL, *saveptr = NULL;
    // int load_data = 0, len_prob_name = 0, n_blocks = 0.;
    // size_t len = 0;
    // FILE *varsep = fopen("varsep.res", "r");
    // if (!varsep) printf("\n\n ERROR: varsep.res file is not available!!! \n\n");
    // else
    // {
    //   // Traverse the problem names in data file for a match
    //   for (int i=0; i<1081; i++)
    //   {
    //     if (getline(&line, &len, varsep))     // Read 1st line = problem name
    //     {
    //       len_prob_name = strlen(line);
    //       if (len_prob_name > 0 && line[len_prob_name - 1] == '\n')
    //       {
    //           line[len_prob_name - 1] = '\0';  // Replace the newline character with null terminator
    //       }
    //       if (strcmp(line, strToken_name) == 0) // Problem names are matched
    //       {
    //         load_data = 1;
    //       }
    //     }
    //     else
    //     {
    //       printf("ERROR: Error reading from varsep.res file.\n");
    //       break;
    //     }


    //     // Read No. of blocks of the current test
    //     if (getline(&line, &len, varsep))     // Read 2nd line = n_blocks
    //     {
    //       n_blocks = atoi(line);
    //     }


    //     if (load_data) break;

    //     // Go to the next problem name
    //     for (int i=0; i<n_blocks; i++)
    //       for (int j=0; j<4; j++) // each block has exaclty 4 lines
    //         getline(&line, &len, varsep);
    //   }
    // }


    // /*** Solve separate sub-problems ***/
    // int blk_index = -1, count_contact_cp = 0, count_contact_full = 0, count_body_cp = 0, count_body_full = 0;
    // int *permu_contact_cp = NULL, *permu_contact_full = NULL, *permu_body_cp = NULL, *permu_body_full = NULL;
    // int Hc_m = -1, Hc_n = -1, Hc_rank = -1;
    // ssize_t read;
    // separators = " \t\n";

    // for (int blk_i=0; blk_i<n_blocks; blk_i++)
    // {
    //   // Block index
    //   getline(&line, &len, varsep);
    //   blk_index = atoi(line);

    //   // Permutation contact
    //   count_contact_cp = 0, count_contact_full = 0;
    //   read = getline(&line, &len, varsep);
    //   permu_contact_cp = malloc(read * sizeof(int)); // cp = compressed
    //   permu_contact_full = malloc(read * sizeof(int) * 3);
    //   char *strToken = strtok (line, separators);

    //   while (strToken != NULL)
    //   {
    //     // Convert the strToken to an integer and store it in permu_contact_cp array
    //     permu_contact_cp[count_contact_cp] = atoi(strToken)-1;
    //     for (int k=0; k<3; k++) permu_contact_full[count_contact_full++] = permu_contact_cp[count_contact_cp]*3+k;
    //     count_contact_cp++;

    //     // Get the next token
    //     strToken = strtok(NULL, separators);
    //   }

    //   // Permutation body
    //   count_body_cp = 0, count_body_full = 0;
    //   read = getline(&line, &len, varsep);
    //   permu_body_cp = malloc(read * sizeof(int));
    //   permu_body_full = malloc(read * sizeof(int) * 6);
    //   strToken = strtok (line, separators);
    //   while (strToken != NULL)
    //   {
    //     permu_body_cp[count_body_cp] = atoi(strToken)-1;
    //     for (int k=0; k<6; k++) permu_body_full[count_body_full++] = permu_body_cp[count_body_cp]*6+k;
    //     count_body_cp++;

    //     strToken = strtok(NULL, separators);
    //   }


    //   // Hc info
    //   getline(&line, &len, varsep);
    //   strToken = strtok (line, separators);
    //   strToken = strtok(NULL, separators);
    //   strToken = strtok(NULL, separators);
    //   Hc_rank = atoi(strToken);


    //   // Create sub-problem
    //   GlobalFrictionContactProblem * sub_prob = globalFrictionContactProblem_new();
    //   sub_prob->dimension = problem->dimension;
    //   sub_prob->numberOfContacts = count_contact_cp;
    //   sub_prob->M = NM_extract(problem->M, count_body_full, permu_body_full, count_body_full, permu_body_full);
    //   sub_prob->env = NULL;

    //   // /* symmetrization of the matrix M */
    //   // if(!(NM_is_symmetric(sub_prob->M)))
    //   // {
    //   //   printf("#################### SYMMETRIZATION ####################\n");
    //   //   NumericsMatrix *MT = NM_transpose(sub_prob->M);
    //   //   NumericsMatrix * MSym = NM_add(1 / 2., sub_prob->M, 1 / 2., MT);
    //   //   NM_free(sub_prob->M);
    //   //   sub_prob->M = MSym;
    //   //   NM_free(MT);
    //   // }

    //   // Attention: matrix H stored in hdf5 is transposed!
    //   sub_prob->H = NM_extract(problem->H, count_body_full, permu_body_full, count_contact_full, permu_contact_full);

    //   // vector f = problem->q
    //   sub_prob->q = (double*)malloc(count_body_full*sizeof(double));
    //   for (int i=0; i<count_body_full; i++) sub_prob->q[i] = problem->q[permu_body_full[i]];

    //   // vector w = problem->b
    //   sub_prob->b = (double*)malloc(count_contact_full*sizeof(double));
    //   for (int i=0; i<count_contact_full; i++) sub_prob->b[i] = problem->b[permu_contact_full[i]];

    //   // friction coef mu
    //   sub_prob->mu = (double*)malloc(count_contact_cp*sizeof(double));
    //   for (int i=0; i<count_contact_cp; i++) sub_prob->mu[i] = problem->mu[permu_contact_cp[i]];

    //   // name
    //   sub_prob->name = problem->name;

    //   // Sub-solutions
    //   double *sub_reaction = (double*)calloc(count_contact_full, sizeof(double));
    //   double *sub_velocity = (double*)calloc(count_contact_full, sizeof(double));
    //   double *sub_globalVelocity = (double*)calloc(count_body_full, sizeof(double));
    //   double *sub_reaction_tmp = (double*)calloc(count_contact_full, sizeof(double));
    //   double *sub_velocity_tmp = (double*)calloc(count_contact_full, sizeof(double));

    //   // Save the block number to append in the test name
    //   options->solverData = (char *)malloc(10*sizeof(char));
    //   char *blk_num_name = (char *)options->solverData;
    //   sprintf(blk_num_name, "-%d", blk_i);


    //   printf("\n********** START sub-problem %d / %d, H (ndxm)/rank = %dx%d / %d *************************************************************\n", blk_i+1, n_blocks, count_contact_full, count_body_full, Hc_rank);
    //   info = -1;
    //   gfc3d_ADMM(sub_prob, sub_reaction, sub_velocity, sub_globalVelocity, &info, options);
    //   if (info)
    //   {
    //     printf("test: failure\n");
    //   }
    //   else
    //   {
    //     printf("test: success\n");
    //   }
    //   printf("*************** END sub-problem %d / %d **************************************************************************************\n\n", blk_i+1, n_blocks);

    //   // Copy sub-solutions to main solutions
    //   for (int i=0; i<count_contact_full; i++)
    //   {
    //     reaction[permu_contact_full[i]] = sub_reaction[i];
    //     velocity[permu_contact_full[i]] = sub_velocity[i];
    //   }


    //   free(sub_reaction); free(sub_velocity); free(sub_globalVelocity);
    //   free(sub_reaction_tmp); free(sub_velocity_tmp);
    //   globalFrictionContact_free(sub_prob);

    //   free(permu_contact_cp); free(permu_contact_full);
    //   free(permu_body_cp); free(permu_body_full);
    // } // loop for n_blocks

    // fclose(varsep);
    // free(str);

    // break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_ADMM_WR:
  {

    numerics_printf_verbose(1," ========================== Call NSGS_WR solver with reformulation into Friction-Contact 3D problem ==========================\n");
    gfc3d_admm_wr(problem, reaction, velocity, globalVelocity, &info, options);
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_IPM:
  case SICONOS_GLOBAL_FRICTION_3D_IPM_SEP:
  {

    GlobalFrictionContactProblem* balanced_problem = gfc3d_balancing_problem(problem,options);
    gfc3d_balancing_go_to_balanced_variables(balanced_problem, options,
                                             reaction, velocity, globalVelocity);

    gfc3d_IPM(balanced_problem, reaction, velocity,
              globalVelocity, &info, options);


    gfc3d_balancing_check_drift(balanced_problem,problem, reaction, velocity, globalVelocity,
                                options);

    balanced_problem = gfc3d_balancing_free(balanced_problem, options);
    break;

  }
  // case SICONOS_GLOBAL_FRICTION_3D_IPM_SEP:
  // {
  //   // Get problem name (without its path)
  //   char *str = (char *) malloc(200);
  //   if (problem->name)
  //   {
  //     strcpy( str, problem->name );
  //   }
  //   else
  //   {
  //     strcpy( str, "foo_" );
  //   }
  //   char * separators = "/";
  //   char *strToken_name = strtok( str, separators );
  //   for(int i=0; i<5; i++)
  //   {
  //     if(strToken_name != NULL) strToken_name = strtok ( NULL, separators );
  //   }
  //   strToken_name = strtok ( strToken_name, "." );


  //   // Find problem name in separation data (file varsep.res), then get its separation info
  //   // Note: this is temporary, need a proper structure (reading hdf5) used in Siconos
  //   char *line = NULL, *saveptr = NULL;
  //   int load_data = 0, len_prob_name = 0, n_blocks = 0.;
  //   size_t len = 0;
  //   FILE *varsep = fopen("varsep.res", "r");
  //   if (!varsep) printf("\n\n ERROR: varsep.res file is not available!!! \n\n");
  //   else
  //   {
  //     // Traverse the problem names in data file for a match
  //     for (int i=0; i<1091; i++)
  //     {
  //       if (getline(&line, &len, varsep))     // Read 1st line = problem name
  //       {
  //         len_prob_name = strlen(line);
  //         if (len_prob_name > 0 && line[len_prob_name - 1] == '\n')
  //         {
  //             line[len_prob_name - 1] = '\0';  // Replace the newline character with null terminator
  //         }
  //         if (strcmp(line, strToken_name) == 0) // Problem names are matched
  //         {
  //           load_data = 1;
  //         }
  //       }
  //       else
  //       {
  //         printf("ERROR: Error reading from varsep.res file.\n");
  //         break;
  //       }


  //       // Read No. of blocks of the current test
  //       if (getline(&line, &len, varsep))     // Read 2nd line = n_blocks
  //       {
  //         n_blocks = atoi(line);
  //       }


  //       if (load_data) break;

  //       // Go to the next problem name
  //       for (int i=0; i<n_blocks; i++)
  //         for (int j=0; j<4; j++) // each block has exaclty 4 lines
  //           getline(&line, &len, varsep);
  //     }
  //   }


  //   /*** Solve separate sub-problems ***/
  //   int blk_index = -1, count_contact_cp = 0, count_contact_full = 0, count_body_cp = 0, count_body_full = 0;
  //   int *permu_contact_cp = NULL, *permu_contact_full = NULL, *permu_body_cp = NULL, *permu_body_full = NULL;
  //   int Hc_m = -1, Hc_n = -1, Hc_rank = -1;
  //   ssize_t read;
  //   separators = " \t\n";


  //   /* For statistiques */
  //   double somme[3], dur[3], cesp = 1e-8, ns, ndur;
  //   int nsc, nN, nB, nR, nT;
  //   int nNcomple_vec[problem->numberOfContacts], nBcomple_vec[problem->numberOfContacts], nRcomple_vec[problem->numberOfContacts], nTcomple_vec[problem->numberOfContacts];
  //   // FILE *stats = fopen("stats_ipm_sep.m", "a+");
  //   FILE *stats = fopen("stats_ipm_sep.m", "w");
  //   // ###########################


  //   // For performance profile
  //   FILE *pp_complete = fopen("pp_complete.pp", "a+");
  //   FILE *pp_subs = fopen("pp_subs.pp", "a+");
  //   int info_subs = 0;
  //   // ###########################


  //   // Timer
  //   long clk_tck = CLOCKS_PER_SEC;
  //   clock_t t1, t2;
  //   double time_sub_probs = 0, time_comple_prob = 0;

  //   // for getting NORM_INF of totalresidu and projection err for all sub-probs
  //   double max_residu_subprobs = -1., max_prjerr_subprobs = -1.;

  //   GlobalFrictionContactProblem* balanced_problem = NULL;
  //   int all_iterations = 0, n_succ = 0, n_fail = 0;
  //   for (int blk_i=0; blk_i<n_blocks; blk_i++)
  //   {
  //     // Block index
  //     getline(&line, &len, varsep);
  //     blk_index = atoi(line);

  //     // Permutation contact
  //     count_contact_cp = 0, count_contact_full = 0;
  //     read = getline(&line, &len, varsep);
  //     permu_contact_cp = malloc(read * sizeof(int)); // cp = compressed
  //     permu_contact_full = malloc(read * sizeof(int) * 3);
  //     char *strToken = strtok (line, separators);

  //     while (strToken != NULL)
  //     {
  //       // Convert the strToken to an integer and store it in permu_contact_cp array
  //       permu_contact_cp[count_contact_cp] = atoi(strToken)-1;
  //       for (int k=0; k<3; k++) permu_contact_full[count_contact_full++] = permu_contact_cp[count_contact_cp]*3+k;
  //       count_contact_cp++;

  //       // Get the next token
  //       strToken = strtok(NULL, separators);
  //     }

  //     // Permutation body
  //     count_body_cp = 0, count_body_full = 0;
  //     read = getline(&line, &len, varsep);
  //     permu_body_cp = malloc(read * sizeof(int));
  //     permu_body_full = malloc(read * sizeof(int) * 6);
  //     strToken = strtok (line, separators);
  //     while (strToken != NULL)
  //     {
  //       permu_body_cp[count_body_cp] = atoi(strToken)-1;
  //       for (int k=0; k<6; k++) permu_body_full[count_body_full++] = permu_body_cp[count_body_cp]*6+k;
  //       count_body_cp++;

  //       strToken = strtok(NULL, separators);
  //     }


  //     // Hc info
  //     getline(&line, &len, varsep);
  //     strToken = strtok (line, separators);
  //     strToken = strtok(NULL, separators);
  //     strToken = strtok(NULL, separators);
  //     Hc_rank = atoi(strToken);

  //     /* For statistiques */
  //     int nN_vec[count_contact_cp], nB_vec[count_contact_cp], nR_vec[count_contact_cp], nT_vec[count_contact_cp];
  //     // ###########################

  //     // Create sub-problem
  //     GlobalFrictionContactProblem * sub_prob = globalFrictionContactProblem_new();
  //     sub_prob->dimension = problem->dimension;
  //     sub_prob->numberOfContacts = count_contact_cp;
  //     sub_prob->M = NM_extract(problem->M, count_body_full, permu_body_full, count_body_full, permu_body_full);
  //     sub_prob->env = NULL;

  //     // Attention: matrix H stored in hdf5 is transposed!
  //     sub_prob->H = NM_extract(problem->H, count_body_full, permu_body_full, count_contact_full, permu_contact_full);

  //     // vector f = problem->q
  //     sub_prob->q = (double*)malloc(count_body_full*sizeof(double));
  //     for (int i=0; i<count_body_full; i++) sub_prob->q[i] = problem->q[permu_body_full[i]];

  //     // vector w = problem->b
  //     sub_prob->b = (double*)malloc(count_contact_full*sizeof(double));
  //     for (int i=0; i<count_contact_full; i++) sub_prob->b[i] = problem->b[permu_contact_full[i]];

  //     // friction coef mu
  //     sub_prob->mu = (double*)malloc(count_contact_cp*sizeof(double));
  //     for (int i=0; i<count_contact_cp; i++) sub_prob->mu[i] = problem->mu[permu_contact_cp[i]];

  //     // Sub-solutions
  //     double *sub_reaction = (double*)calloc(count_contact_full, sizeof(double));
  //     double *sub_velocity = (double*)calloc(count_contact_full, sizeof(double));
  //     double *sub_globalVelocity = (double*)calloc(count_body_full, sizeof(double));
  //     double *sub_reaction_tmp = (double*)calloc(count_contact_full, sizeof(double));
  //     double *sub_velocity_tmp = (double*)calloc(count_contact_full, sizeof(double));


  //     balanced_problem = gfc3d_balancing_problem(sub_prob,options);
  //     gfc3d_balancing_go_to_balanced_variables(balanced_problem, options,
  //                                            sub_reaction, sub_velocity, sub_globalVelocity);
  //     printf("\n********** START sub-problem %d / %d, H (ndxm)/rank = %dx%d / %d *************************************************************\n", blk_i+1, n_blocks, count_contact_full, count_body_full, Hc_rank);
  //     // options->dparam[SICONOS_DPARAM_TOL] = 1e-10*sub_prob->numberOfContacts/problem->numberOfContacts;
  //     t1 = clock();
  //     gfc3d_IPM(balanced_problem, sub_reaction, sub_velocity,
  //             sub_globalVelocity, &info, options);
  //     t2 = clock();
  //     time_sub_probs += (double)(t2-t1)/(double)clk_tck;
  //     all_iterations += options->iparam[SICONOS_IPARAM_ITER_DONE];
  //     max_residu_subprobs = fmax(options->dparam[SICONOS_DPARAM_RESIDU], max_residu_subprobs);
  //     double *ptr = options->solverData;
  //     max_prjerr_subprobs = fmax(*ptr, max_prjerr_subprobs);
  //     free(options->solverData); options->solverData = NULL;
  //     gfc3d_balancing_check_drift(balanced_problem,problem, reaction, velocity, globalVelocity,
  //                               options);
  //     balanced_problem = gfc3d_balancing_free(balanced_problem, options);


  //     /* For statistiques */
  //     NumericsMatrix *P_mu_sub = NM_create(NM_SPARSE, count_contact_full, count_contact_full);
  //     NumericsMatrix *P_mu_inv_sub = NM_create(NM_SPARSE, count_contact_full, count_contact_full);
  //     NM_triplet_alloc(P_mu_sub, count_contact_full);
  //     NM_triplet_alloc(P_mu_inv_sub, count_contact_full);
  //     P_mu_sub->matrix2->origin = NSM_TRIPLET;
  //     P_mu_inv_sub->matrix2->origin = NSM_TRIPLET;
  //     for(unsigned int i = 0; i < count_contact_full; ++i)
  //       if(i % 3 == 0)
  //       {
  //         NM_entry(P_mu_sub, i, i, 1.);
  //         NM_entry(P_mu_inv_sub, i, i, 1.);
  //       }
  //       else
  //       {
  //         NM_entry(P_mu_sub, i, i, sub_prob->mu[(int)(i/3)]);
  //         NM_entry(P_mu_inv_sub, i, i, 1.0/sub_prob->mu[(int)(i/3)]);
  //       }
  //     NM_gemv(1.0, P_mu_sub, sub_velocity, 0.0, sub_velocity_tmp);
  //     NM_gemv(1.0, P_mu_inv_sub, sub_reaction, 0.0, sub_reaction_tmp);
  //     double projerr_sub = projectionError(sub_velocity_tmp, sub_reaction_tmp, count_contact_cp, options->dparam[SICONOS_DPARAM_TOL]);
  //     nsc = 0, nN = 0, nB = 0, nR = 0, nT = 0;
  //     for (int i = 0; i < count_contact_cp; i++)
  //     {
  //       somme[0] = sub_velocity_tmp[3*i] + sub_reaction_tmp[3*i];
  //       somme[1] = sub_velocity_tmp[3*i+1] + sub_reaction_tmp[3*i+1];
  //       somme[2] = sub_velocity_tmp[3*i+2] + sub_reaction_tmp[3*i+2];
  //       dur[0] = sub_velocity_tmp[3*i] - sub_reaction_tmp[3*i];
  //       dur[1] = sub_velocity_tmp[3*i+1] - sub_reaction_tmp[3*i+1];
  //       dur[2] = sub_velocity_tmp[3*i+2] - sub_reaction_tmp[3*i+2];

  //       ns = somme[0] - cblas_dnrm2(2,somme+1,1);
  //       ndur = dur[0] - cblas_dnrm2(2,dur+1,1);
  //       if (ns > cesp*cblas_dnrm2(3, somme, 1))
  //       {
  //         nsc +=1;
  //         if (dur[0] >= cblas_dnrm2(2,dur+1,1))       { nN_vec[nN++] = i+1; }
  //         else if (-dur[0] >= cblas_dnrm2(2,dur+1,1)) { nB_vec[nB++] = i+1; }
  //         else                                        { nR_vec[nR++] = i+1; }
  //       }
  //       else
  //         nT_vec[nT++] = i+1;
  //     }

  //     fprintf(stats, "stats(end+1) = struct('type', 'sub', 'name', '%s', 'blk_id', %d, 'nBlk', %d, 'mu', %.2e, 'ite', %d, 'residu', %.2e, 'projerr', %.2e,",
  //                     strToken_name, blk_i+1, n_blocks, problem->mu[0], options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU], projerr_sub);

  //     // For performance profile
  //     // For all sub-problems, if there is only one failure,
  //     // solving these sub-problems is considered a failure
  //     info_subs += info;
  //     // ###########################

  //     if (info)
  //     {
  //       printf("test: failure\n");
  //       n_fail++;
  //       fprintf(stats,"'result', 'fail',");
  //     }
  //     else
  //     {
  //       printf("test: success\n");
  //       n_succ++;
  //       fprintf(stats,"'result', 'succ',");
  //     }

  //     fprintf(stats,"'nd', %d, 'm', %d, 'r_H', %d,", count_contact_full, count_body_full, Hc_rank);
  //     fprintf(stats,"'B', ["); for (int i=0; i<nB; i++) fprintf(stats, " %d", nB_vec[i]); fprintf(stats,"],");
  //     fprintf(stats,"'N', ["); for (int i=0; i<nN; i++) fprintf(stats, " %d", nN_vec[i]); fprintf(stats,"],");
  //     fprintf(stats,"'R', ["); for (int i=0; i<nR; i++) fprintf(stats, " %d", nR_vec[i]); fprintf(stats,"],");
  //     fprintf(stats,"'T', ["); for (int i=0; i<nT; i++) fprintf(stats, " %d", nT_vec[i]); fprintf(stats,"]);\n");
  //     // ###########################
  //     printf("*************** END sub-problem %d / %d **************************************************************************************\n\n", blk_i+1, n_blocks);

  //     // Copy sub-solutions to main solutions
  //     for (int i=0; i<count_contact_full; i++)
  //     {
  //       reaction[permu_contact_full[i]] = sub_reaction[i];
  //       velocity[permu_contact_full[i]] = sub_velocity[i];
  //     }

  //     // for (int i=0; i<count_body_full; i++)
  //     // {
  //     //   globalVelocity[permu_body_full[i]] = sub_globalVelocity[i];
  //     // }


  //     free(sub_reaction); free(sub_velocity); free(sub_globalVelocity);
  //     free(sub_reaction_tmp); free(sub_velocity_tmp);
  //     globalFrictionContact_free(sub_prob);

  //     free(permu_contact_cp); free(permu_contact_full);
  //     free(permu_body_cp); free(permu_body_full);

  //     if(P_mu_sub) {P_mu_sub = NM_free(P_mu_sub); P_mu_sub = NULL;}
  //     if(P_mu_inv_sub) {P_mu_inv_sub = NM_free(P_mu_inv_sub); P_mu_inv_sub = NULL;}
  //   } // loop for n_blocks



  //   // Double-check the result
  //   unsigned int nd = problem->H->size1;
  //   unsigned int d = problem->dimension;
  //   unsigned int n = problem->numberOfContacts;
  //   unsigned int m = problem->M->size0;

  //   NumericsMatrix *P_mu = NM_create(NM_SPARSE, nd, nd);
  //   NumericsMatrix *P_mu_inv = NM_create(NM_SPARSE, nd, nd);
  //   NM_triplet_alloc(P_mu, nd);
  //   NM_triplet_alloc(P_mu_inv, nd);
  //   P_mu->matrix2->origin = NSM_TRIPLET;
  //   P_mu_inv->matrix2->origin = NSM_TRIPLET;
  //   for(unsigned int i = 0; i < nd; ++i)
  //     if(i % d == 0)
  //     {
  //       NM_entry(P_mu, i, i, 1.);
  //       NM_entry(P_mu_inv, i, i, 1.);
  //     }
  //     else
  //     {
  //       NM_entry(P_mu, i, i, problem->mu[(int)(i/d)]);
  //       NM_entry(P_mu_inv, i, i, 1.0/problem->mu[(int)(i/d)]);
  //     }

  //   NumericsMatrix *M = problem->M;
  //   NumericsMatrix *H_tilde = NM_transpose(problem->H);
  //   double *w_tilde = problem->b;
  //   double *w = (double*)calloc(nd, sizeof(double));
  //   double *f = problem->q;

  //   double *velocity_tmp = (double*)calloc(nd, sizeof(double));
  //   double *reaction_tmp = (double*)calloc(nd, sizeof(double));
  //   double *primalConstraint = (double*)calloc(nd, sizeof(double));
  //   double *dualConstraint = (double*)calloc(m, sizeof(double));
  //   double *s = (double*)calloc(n, sizeof(double));

  //   double pinfeas, dinfeas, complem, udotr, nub, diff_u, diff_r;

  //   // Change of variable
  //   // Current velocity = [u0; mu*ub]
  //   NumericsMatrix *H = NM_multiply(P_mu, H_tilde);
  //   NM_gemv(1.0, P_mu, w_tilde, 0.0, w);


  //   globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);


  //   // Compute vars without friction coef.
  //   NM_gemv(1.0, P_mu, velocity, 0.0, velocity_tmp);
  //   NM_gemv(1.0, P_mu_inv, reaction, 0.0, reaction_tmp);
  //   // for (unsigned int i = 0; i<n; i++) s[i] = cblas_dnrm2(2, velocity_tmp+i*d+1, 1);
  //   // primalResidual_s(velocity_tmp, H, globalVelocity, w, s, primalConstraint, &pinfeas, options->dparam[SICONOS_DPARAM_TOL]);
  //   dualResidual(M, globalVelocity, H, reaction_tmp, f, dualConstraint, &dinfeas, options->dparam[SICONOS_DPARAM_TOL]);
  //   complem = complemResidualNorm(velocity_tmp, reaction_tmp, nd, n);
  //   udotr = cblas_ddot(nd, velocity_tmp, 1, reaction_tmp, 1);

  //   double norm_q = cblas_dnrm2(m, problem->q, 1);
  //   double norm_b = cblas_dnrm2(nd, problem->b, 1);
  //   double projerr, projerr_comple;
  //   projerr = projectionError(velocity_tmp, reaction_tmp, n, options->dparam[SICONOS_DPARAM_TOL]);


  //   // Store solutions obtained by solving sub-problems
  //   cblas_dcopy(nd, velocity, 1, velocity_tmp, 1);
  //   cblas_dcopy(nd, reaction, 1, reaction_tmp, 1);


  //   balanced_problem = gfc3d_balancing_problem(problem,options);
  //   gfc3d_balancing_go_to_balanced_variables(balanced_problem, options,
  //                                            reaction, velocity, globalVelocity);
  //   printf("\n\n============ Solve the complete problem by gfc3d_IPM ================\n");
  //   // options->dparam[SICONOS_DPARAM_TOL] = 1e-10;
  //   t1 = clock();
  //   gfc3d_IPM(balanced_problem, reaction, velocity, globalVelocity, &info, options);
  //   t2 = clock();
  //   time_comple_prob = (double)(t2-t1)/(double)clk_tck;
  //   printf("Execution time: %.4f,\t\ttest: ", time_comple_prob);
  //   // if (info) printf("failure\n"); else printf("success\n");
  //   gfc3d_balancing_check_drift(balanced_problem,problem, reaction, velocity, globalVelocity,
  //                               options);
  //   balanced_problem = gfc3d_balancing_free(balanced_problem, options);


  //   // for comparison
  //   cblas_daxpy(nd, -1, velocity, 1, velocity_tmp, 1);
  //   cblas_daxpy(nd, -1, reaction, 1, reaction_tmp, 1);
  //   diff_u = cblas_dnrm2(nd, velocity_tmp, 1);
  //   diff_r = cblas_dnrm2(nd, reaction_tmp, 1);
  //   NM_gemv(1.0, P_mu, velocity, 0.0, velocity_tmp);
  //   NM_gemv(1.0, P_mu_inv, reaction, 0.0, reaction_tmp);
  //   projerr_comple = projectionError(velocity_tmp, reaction_tmp, n, options->dparam[SICONOS_DPARAM_TOL]);



  //   /* For statistiques */
  //   nsc = 0, nN = 0, nB = 0, nR = 0, nT = 0;
  //   for (int i = 0; i < problem->numberOfContacts; i++)
  //   {
  //     somme[0] = velocity_tmp[3*i] + reaction_tmp[3*i];
  //     somme[1] = velocity_tmp[3*i+1] + reaction_tmp[3*i+1];
  //     somme[2] = velocity_tmp[3*i+2] + reaction[3*i+2];
  //     dur[0] = velocity_tmp[3*i] - reaction_tmp[3*i];
  //     dur[1] = velocity_tmp[3*i+1] - reaction_tmp[3*i+1];
  //     dur[2] = velocity_tmp[3*i+2] - reaction_tmp[3*i+2];

  //     ns = somme[0] - cblas_dnrm2(2,somme+1,1);
  //     ndur = dur[0] - cblas_dnrm2(2,dur+1,1);
  //     if (ns > cesp*cblas_dnrm2(3, somme, 1))
  //     {
  //       nsc +=1;
  //       if (dur[0] >= cblas_dnrm2(2,dur+1,1))       { nNcomple_vec[nN++] = i+1; }
  //       else if (-dur[0] >= cblas_dnrm2(2,dur+1,1)) { nBcomple_vec[nB++] = i+1; }
  //       else                                        { nRcomple_vec[nR++] = i+1; }
  //     }
  //     else
  //       nTcomple_vec[nT++] = i+1;
  //   }

  //   fprintf(stats, "stats(end+1) = struct('type', 'complete', 'name', '%s', 'blk_id', %d, 'nBlk', %d, 'mu', %.2e, 'ite', %d, 'residu', %.2e, 'projerr', %.2e,",
  //                   strToken_name, -1, n_blocks, problem->mu[0], options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU], projerr_comple);

  //   if (info)
  //   {
  //     printf("failure\n");
  //     // n_fail++;
  //     fprintf(stats,"'result', 'fail',");
  //   }
  //   else
  //   {
  //     printf("success\n");
  //     // n_succ++;
  //     fprintf(stats,"'result', 'succ',");
  //   }

  //   fprintf(stats,"'nd', %d, 'm', %d, 'r_H', %d,", nd, m, -1);
  //   fprintf(stats,"'B', ["); for (int i=0; i<nB; i++) fprintf(stats, " %d", nBcomple_vec[i]); fprintf(stats,"],");
  //   fprintf(stats,"'N', ["); for (int i=0; i<nN; i++) fprintf(stats, " %d", nNcomple_vec[i]); fprintf(stats,"],");
  //   fprintf(stats,"'R', ["); for (int i=0; i<nR; i++) fprintf(stats, " %d", nRcomple_vec[i]); fprintf(stats,"],");
  //   fprintf(stats,"'T', ["); for (int i=0; i<nT; i++) fprintf(stats, " %d", nTcomple_vec[i]); fprintf(stats,"]);\n");
  //   fclose(stats);
  //   // ###########################
  //   printf("=====================================================================\n\n");



  //   printf("sumry total:   %s,\t  n,  nd x m = %d,  %d x %d\n",strToken_name, n, nd, m);
  //   printf("sumry total: %4d sub-prob, %4d succ, %4d fail,\t %.2e %.2e,\t %.4f (s),\t | u_subs - u_comple | = %.2e\n", n_blocks, n_succ, n_fail, max_residu_subprobs, max_prjerr_subprobs, time_sub_probs, diff_u);

  //   if (info)
  //     printf("sumry total:   comple prob,   failure           ,");
  //   else
  //     printf("sumry total:   comple prob,   success           ,");

  //   printf("\t %.2e %.2e,\t %.4f (s),\t | r_subs - r_comple | = %.2e\n",
  //           options->dparam[SICONOS_DPARAM_RESIDU], projerr_comple, time_comple_prob, diff_r);

  //   printf("sumry total: --------------------------------------------------------------------------------------------------------------------\n\n");


  //   // For performance profile
  //   fprintf(pp_subs,"sumry: %d  %9.2e  %5i  %10.4f", info_subs, max_residu_subprobs, all_iterations/n_blocks, time_sub_probs);
  //   fprintf(pp_subs,"%3i %5i %5i     %s\n", d, n, m, strToken_name);
  //   fprintf(pp_complete,"sumry: %d  %9.2e  %5i  %10.4f", info, options->dparam[SICONOS_DPARAM_RESIDU], options->iparam[SICONOS_IPARAM_ITER_DONE], time_comple_prob);
  //   fprintf(pp_complete,"%3i %5i %5i     %s\n", d, n, m, strToken_name);
  //   fclose(pp_subs);
  //   fclose(pp_complete);
  //   // ###########################

  //   if(H_tilde) {H_tilde = NM_free(H_tilde); H_tilde = NULL;}
  //   if(H) {H = NM_free(H); H = NULL;}
  //   if(P_mu) {P_mu = NM_free(P_mu); P_mu = NULL;}
  //   if(P_mu_inv) {P_mu_inv = NM_free(P_mu_inv); P_mu_inv = NULL;}
  //   if (w) {free(w); w = NULL;}
  //   if (s) {free(s); s = NULL;}
  //   if (velocity_tmp) {free(velocity_tmp); velocity_tmp = NULL;}
  //   if (reaction_tmp) {free(reaction_tmp); reaction_tmp = NULL;}
  //   if (primalConstraint) {free(primalConstraint); primalConstraint = NULL;}
  //   if (dualConstraint) {free(dualConstraint); dualConstraint = NULL;}
  //   fclose(varsep);
  //   free(str);
  //   break;

  // }
  case SICONOS_GLOBAL_FRICTION_3D_IPM_WR:
  {

    /* gfc3d_ipm_wr(problem, reaction, velocity, */
    /*           globalVelocity, &info, options); */
    break;

  }
  case SICONOS_GLOBAL_FRICTION_3D_IPM_SNM:
  {
    gfc3d_IPM_SNM(problem, reaction, velocity,
              globalVelocity, &info, options);
    break;
  }
  case SICONOS_GLOBAL_FRICTION_3D_IPM_SNM_PROX:
  {
    verbose=1;
    gfc3d_IPM_SNM(problem, reaction, velocity,
              globalVelocity, &info, options);

    SolverOptions * nsn_options = options->internalSolvers[0];

    gfc3d_proximal_wr(problem, reaction, velocity, globalVelocity, &info, nsn_options);

    if (problem->name)
      {
    numerics_printf_verbose(1, "problem = %s --  NSN_PROX gfc3d_error = %e, prox iterations %i, cumulative newton iterations %i \n",problem->name, nsn_options->dparam[1], nsn_options->iparam[1], nsn_options->iparam[SICONOS_FRICTION_3D_PROXIMAL_IPARAM_CUMULATIVE_ITER_DONE]);
      }
    else
      {
	numerics_printf_verbose(1, " --  NSN_PROX gfc3d_error = %e, prox iterations %i, cumulative newton iterations %i \n", nsn_options->dparam[1], nsn_options->iparam[1], nsn_options->iparam[SICONOS_FRICTION_3D_PROXIMAL_IPARAM_CUMULATIVE_ITER_DONE]);
      }

    options->iparam[SICONOS_IPARAM_ITER_DONE]=options->iparam[0]+nsn_options->iparam[SICONOS_FRICTION_3D_PROXIMAL_IPARAM_CUMULATIVE_ITER_DONE];
    options->dparam[SICONOS_DPARAM_RESIDU] = nsn_options->dparam[SICONOS_DPARAM_RESIDU];

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
    gfc3d_IPM_SNM(problem, reaction, velocity,
              globalVelocity, &info, options);
    // char *str = (char *) malloc(200);
    // if (problem->name)
    // {
    //   strcpy( str, problem->name );
    // }
    // else
    // {
    //   strcpy( str, "foo_" );
    // }
    // char * separators = "/";
    // char *strToken_name = strtok( str, separators );
    // for(int i=0; i<5; i++)
    // {
    //   if(strToken_name != NULL) strToken_name = strtok ( NULL, separators );
    // }
    // strToken_name = strtok ( strToken_name, "." );


    // // // For performance profile
    // // FILE *pp_complete, *pp_subs;
    // // char pp_name_complete[20], pp_name_subs[20];
    // // // ###########################


    // // int info_subs = 0;
    // // ###########################

    // // // TESTING: change of friction coef.
    // // double mu_start = 0.;
    // // for (int index_mu=0; index_mu<10; index_mu++)
    // // {
    // //   mu_start += 0.1;
    // //   for(int i = 0; i < problem->numberOfContacts ; i++) problem->mu[i]= mu_start;
    // // for(int i = 0; i < problem->numberOfContacts ; i++) problem->mu[i]= 0.3;



    // char *line = NULL, *saveptr = NULL;
    // int load_data = 0, len_prob_name = 0, n_blocks = 0.;
    // size_t len = 0;
    // FILE *varsep = fopen("varsep.res", "r");
    // if (!varsep) printf("\n\n ERROR: varsep.res file is not available!!! \n\n");
    // else
    // {
    //   // Traverse the problem names in data file for a match
    //   for (int i=0; i<1081; i++)
    //   {
    //     if (getline(&line, &len, varsep))     // Read 1st line = problem name
    //     {
    //       len_prob_name = strlen(line);
    //       if (len_prob_name > 0 && line[len_prob_name - 1] == '\n')
    //       {
    //           line[len_prob_name - 1] = '\0';  // Replace the newline character with null terminator
    //       }
    //       if (strcmp(line, strToken_name) == 0) // Problem names are matched
    //       {
    //         load_data = 1;
    //       }
    //     }
    //     else
    //     {
    //       printf("ERROR: Error reading from varsep.res file.\n");
    //       break;
    //     }


    //     // Read No. of blocks of the current test
    //     if (getline(&line, &len, varsep))     // Read 2nd line = n_blocks
    //     {
    //       n_blocks = atoi(line);
    //     }


    //     if (load_data) break;

    //     // Go to the next problem name
    //     for (int i=0; i<n_blocks; i++)
    //       for (int j=0; j<4; j++) // each block has exaclty 4 lines
    //         getline(&line, &len, varsep);
    //   }
    // }


    // /*** Solve separate sub-problems ***/
    // int blk_index = -1, count_contact_cp = 0, count_contact_full = 0, count_body_cp = 0, count_body_full = 0;
    // int *permu_contact_cp = NULL, *permu_contact_full = NULL, *permu_body_cp = NULL, *permu_body_full = NULL;
    // int Hc_m = -1, Hc_n = -1, Hc_rank = -1;
    // ssize_t read;
    // separators = " \t\n";




    // // /* For statistiques */
    // // double somme[3], dur[3], cesp = 1e-8, ns, ndur;
    // // int nsc, nN, nB, nR, nT;
    // // int nNcomple_vec[problem->numberOfContacts], nBcomple_vec[problem->numberOfContacts], nRcomple_vec[problem->numberOfContacts], nTcomple_vec[problem->numberOfContacts];
    // // // FILE *stats = fopen("stats_Spheres1mm_4204.m", "w");
    // // // ###########################




    // // Timer
    // long clk_tck = CLOCKS_PER_SEC;
    // clock_t t1, t2;
    // double time_sub_probs = 0, time_comple_prob = 0;

    // // for getting NORM_INF of totalresidu and projection err for all sub-probs
    // double max_residu_subprobs = -1., max_prjerr_subprobs = -1.;


    // int all_iterations = 0, n_succ = 0, n_fail = 0;
    // for (int blk_i=0; blk_i<n_blocks; blk_i++)
    // {
    //   // Block index
    //   getline(&line, &len, varsep);
    //   blk_index = atoi(line);

    //   // Permutation contact
    //   count_contact_cp = 0, count_contact_full = 0;
    //   read = getline(&line, &len, varsep);
    //   permu_contact_cp = malloc(read * sizeof(int)); // cp = compressed
    //   permu_contact_full = malloc(read * sizeof(int) * 3);
    //   char *strToken = strtok (line, separators);

    //   while (strToken != NULL)
    //   {
    //     // Convert the strToken to an integer and store it in permu_contact_cp array
    //     permu_contact_cp[count_contact_cp] = atoi(strToken)-1;
    //     for (int k=0; k<3; k++) permu_contact_full[count_contact_full++] = permu_contact_cp[count_contact_cp]*3+k;
    //     count_contact_cp++;

    //     // Get the next token
    //     strToken = strtok(NULL, separators);
    //   }

    //   // Permutation body
    //   count_body_cp = 0, count_body_full = 0;
    //   read = getline(&line, &len, varsep);
    //   permu_body_cp = malloc(read * sizeof(int));
    //   permu_body_full = malloc(read * sizeof(int) * 6);
    //   strToken = strtok (line, separators);
    //   while (strToken != NULL)
    //   {
    //     permu_body_cp[count_body_cp] = atoi(strToken)-1;
    //     for (int k=0; k<6; k++) permu_body_full[count_body_full++] = permu_body_cp[count_body_cp]*6+k;
    //     count_body_cp++;

    //     strToken = strtok(NULL, separators);
    //   }


    //   // Hc info
    //   getline(&line, &len, varsep);
    //   strToken = strtok (line, separators);
    //   strToken = strtok(NULL, separators);
    //   strToken = strtok(NULL, separators);
    //   Hc_rank = atoi(strToken);


    //   // /* For statistiques */
    //   // int nN_vec[count_contact_cp], nB_vec[count_contact_cp], nR_vec[count_contact_cp], nT_vec[count_contact_cp];
    //   // // ###########################

    //   // Create sub-problem
    //   GlobalFrictionContactProblem * sub_prob = globalFrictionContactProblem_new();
    //   sub_prob->dimension = problem->dimension;
    //   sub_prob->numberOfContacts = count_contact_cp;
    //   sub_prob->M = NM_extract(problem->M, count_body_full, permu_body_full, count_body_full, permu_body_full);
    //   sub_prob->env = NULL;

    //   // Attention: matrix H stored in hdf5 is transposed!
    //   sub_prob->H = NM_extract(problem->H, count_body_full, permu_body_full, count_contact_full, permu_contact_full);

    //   // vector f = problem->q
    //   sub_prob->q = (double*)malloc(count_body_full*sizeof(double));
    //   for (int i=0; i<count_body_full; i++) sub_prob->q[i] = problem->q[permu_body_full[i]];

    //   // vector w = problem->b
    //   sub_prob->b = (double*)malloc(count_contact_full*sizeof(double));
    //   for (int i=0; i<count_contact_full; i++) sub_prob->b[i] = problem->b[permu_contact_full[i]];

    //   // friction coef mu
    //   sub_prob->mu = (double*)malloc(count_contact_cp*sizeof(double));
    //   for (int i=0; i<count_contact_cp; i++) sub_prob->mu[i] = problem->mu[permu_contact_cp[i]];

    //   // name
    //   sub_prob->name = problem->name;

    //   // Sub-solutions
    //   double *sub_reaction = (double*)calloc(count_contact_full, sizeof(double));
    //   double *sub_velocity = (double*)calloc(count_contact_full, sizeof(double));
    //   double *sub_globalVelocity = (double*)calloc(count_body_full, sizeof(double));
    //   double *sub_reaction_tmp = (double*)calloc(count_contact_full, sizeof(double));
    //   double *sub_velocity_tmp = (double*)calloc(count_contact_full, sizeof(double));

    //   // Save the block number to append in the test name
    //   options->solverData = (char *)malloc(10*sizeof(char));
    //   char *blk_num_name = (char *)options->solverData;
    //   sprintf(blk_num_name, "-%d", blk_i);

    //   printf("\n********** START sub-problem %d / %d, H (ndxm)/rank = %dx%d / %d *************************************************************\n", blk_i+1, n_blocks, count_contact_full, count_body_full, Hc_rank);
    //   t1 = clock();
    //   gfc3d_IPM_SNM(sub_prob, sub_reaction, sub_velocity,
    //           sub_globalVelocity, &info, options);
    //   t2 = clock();
    //   time_sub_probs += (double)(t2-t1)/(double)clk_tck;
    //   all_iterations += options->iparam[SICONOS_IPARAM_ITER_DONE];
    //   max_residu_subprobs = fmax(options->dparam[SICONOS_DPARAM_RESIDU], max_residu_subprobs);
    //   double *ptr = options->solverData;
    //   max_prjerr_subprobs = fmax(*ptr, max_prjerr_subprobs);
    //   free(options->solverData); options->solverData = NULL;


    //   // /* For statistiques */
    //   // NumericsMatrix *P_mu_sub = NM_create(NM_SPARSE, count_contact_full, count_contact_full);
    //   // NumericsMatrix *P_mu_inv_sub = NM_create(NM_SPARSE, count_contact_full, count_contact_full);
    //   // NM_triplet_alloc(P_mu_sub, count_contact_full);
    //   // NM_triplet_alloc(P_mu_inv_sub, count_contact_full);
    //   // P_mu_sub->matrix2->origin = NSM_TRIPLET;
    //   // P_mu_inv_sub->matrix2->origin = NSM_TRIPLET;
    //   // for(unsigned int i = 0; i < count_contact_full; ++i)
    //   //   if(i % 3 == 0)
    //   //   {
    //   //     NM_entry(P_mu_sub, i, i, 1.);
    //   //     NM_entry(P_mu_inv_sub, i, i, 1.);
    //   //   }
    //   //   else
    //   //   {
    //   //     NM_entry(P_mu_sub, i, i, sub_prob->mu[(int)(i/3)]);
    //   //     NM_entry(P_mu_inv_sub, i, i, 1.0/sub_prob->mu[(int)(i/3)]);
    //   //   }
    //   // NM_gemv(1.0, P_mu_sub, sub_velocity, 0.0, sub_velocity_tmp);
    //   // NM_gemv(1.0, P_mu_inv_sub, sub_reaction, 0.0, sub_reaction_tmp);
    //   // double projerr_sub = projectionError(sub_velocity_tmp, sub_reaction_tmp, count_contact_cp, options->dparam[SICONOS_DPARAM_TOL]);
    //   // nsc = 0, nN = 0, nB = 0, nR = 0, nT = 0;
    //   // for (int i = 0; i < count_contact_cp; i++)
    //   // {
    //   //   somme[0] = sub_velocity_tmp[3*i] + sub_reaction_tmp[3*i];
    //   //   somme[1] = sub_velocity_tmp[3*i+1] + sub_reaction_tmp[3*i+1];
    //   //   somme[2] = sub_velocity_tmp[3*i+2] + sub_reaction_tmp[3*i+2];
    //   //   dur[0] = sub_velocity_tmp[3*i] - sub_reaction_tmp[3*i];
    //   //   dur[1] = sub_velocity_tmp[3*i+1] - sub_reaction_tmp[3*i+1];
    //   //   dur[2] = sub_velocity_tmp[3*i+2] - sub_reaction_tmp[3*i+2];

    //   //   ns = somme[0] - cblas_dnrm2(2,somme+1,1);
    //   //   ndur = dur[0] - cblas_dnrm2(2,dur+1,1);
    //   //   if (ns > cesp*cblas_dnrm2(3, somme, 1))
    //   //   {
    //   //     nsc +=1;
    //   //     if (dur[0] >= cblas_dnrm2(2,dur+1,1))       { nN_vec[nN++] = i+1; }
    //   //     else if (-dur[0] >= cblas_dnrm2(2,dur+1,1)) { nB_vec[nB++] = i+1; }
    //   //     else                                        { nR_vec[nR++] = i+1; }
    //   //   }
    //   //   else
    //   //     nT_vec[nT++] = i+1;
    //   // }

    //   // fprintf(stats, "stats(end+1) = struct('type', 'sub', 'name', '%s', 'blk_id', %d, 'nBlk', %d, 'mu', %.2e, 'ite', %d, 'residu', %.2e, 'projerr', %.2e,",
    //   //                 strToken_name, blk_i+1, n_blocks, problem->mu[0], options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU], projerr_sub);


    //   // For performance profile
    //   // For all sub-problems, if there is only one failure,
    //   // solving these sub-problems is considered a failure
    //   // info_subs += info;
    //   // ###########################



    //   if (info)
    //   {
    //     printf("test: failure\n");
    //     n_fail++;
    //     // fprintf(stats,"'result', 'fail',");
    //   }
    //   else
    //   {
    //     printf("test: success\n");
    //     n_succ++;
    //     // fprintf(stats,"'result', 'succ',");
    //   }

    //   // fprintf(stats,"'nd', %d, 'm', %d, 'r_H', %d,", count_contact_full, count_body_full, Hc_rank);
    //   // fprintf(stats,"'B', ["); for (int i=0; i<nB; i++) fprintf(stats, " %d", nB_vec[i]); fprintf(stats,"],");
    //   // fprintf(stats,"'N', ["); for (int i=0; i<nN; i++) fprintf(stats, " %d", nN_vec[i]); fprintf(stats,"],");
    //   // fprintf(stats,"'R', ["); for (int i=0; i<nR; i++) fprintf(stats, " %d", nR_vec[i]); fprintf(stats,"],");
    //   // fprintf(stats,"'T', ["); for (int i=0; i<nT; i++) fprintf(stats, " %d", nT_vec[i]); fprintf(stats,"]);\n");
    //   // ###########################
    //   printf("*************** END sub-problem %d / %d **************************************************************************************\n\n", blk_i+1, n_blocks);

    //   // Copy sub-solutions to main solutions
    //   for (int i=0; i<count_contact_full; i++)
    //   {
    //     reaction[permu_contact_full[i]] = sub_reaction[i];
    //     velocity[permu_contact_full[i]] = sub_velocity[i];
    //   }

    //   // for (int i=0; i<count_body_full; i++)
    //   // {
    //   //   globalVelocity[permu_body_full[i]] = sub_globalVelocity[i];
    //   // }


    //   free(sub_reaction); free(sub_velocity); free(sub_globalVelocity);
    //   free(sub_reaction_tmp); free(sub_velocity_tmp);
    //   globalFrictionContact_free(sub_prob);

    //   free(permu_contact_cp); free(permu_contact_full);
    //   free(permu_body_cp); free(permu_body_full);

    //   // if(P_mu_sub) {P_mu_sub = NM_free(P_mu_sub); P_mu_sub = NULL;}
    //   // if(P_mu_inv_sub) {P_mu_inv_sub = NM_free(P_mu_inv_sub); P_mu_inv_sub = NULL;}
    // } // loop for n_blocks



    // // // Double-check the result
    // // unsigned int nd = problem->H->size1;
    // // unsigned int d = problem->dimension;
    // // unsigned int n = problem->numberOfContacts;
    // // unsigned int m = problem->M->size0;

    // // NumericsMatrix *P_mu = NM_create(NM_SPARSE, nd, nd);
    // // NumericsMatrix *P_mu_inv = NM_create(NM_SPARSE, nd, nd);
    // // NM_triplet_alloc(P_mu, nd);
    // // NM_triplet_alloc(P_mu_inv, nd);
    // // P_mu->matrix2->origin = NSM_TRIPLET;
    // // P_mu_inv->matrix2->origin = NSM_TRIPLET;
    // // for(unsigned int i = 0; i < nd; ++i)
    // //   if(i % d == 0)
    // //   {
    // //     NM_entry(P_mu, i, i, 1.);
    // //     NM_entry(P_mu_inv, i, i, 1.);
    // //   }
    // //   else
    // //   {
    // //     NM_entry(P_mu, i, i, problem->mu[(int)(i/d)]);
    // //     NM_entry(P_mu_inv, i, i, 1.0/problem->mu[(int)(i/d)]);
    // //   }

    // // NumericsMatrix *M = problem->M;
    // // NumericsMatrix *H_tilde = NM_transpose(problem->H);
    // // double *w_tilde = problem->b;
    // // double *w = (double*)calloc(nd, sizeof(double));
    // // double *f = problem->q;

    // // double *velocity_tmp = (double*)calloc(nd, sizeof(double));
    // // double *reaction_tmp = (double*)calloc(nd, sizeof(double));
    // // double *primalConstraint = (double*)calloc(nd, sizeof(double));
    // // double *dualConstraint = (double*)calloc(m, sizeof(double));
    // // double *s = (double*)calloc(n, sizeof(double));

    // // double pinfeas, dinfeas, complem, udotr, nub, diff_u, diff_r;

    // // // Change of variable
    // // // Current velocity = [u0; mu*ub]
    // // NumericsMatrix *H = NM_multiply(P_mu, H_tilde);
    // // NM_gemv(1.0, P_mu, w_tilde, 0.0, w);


    // // globalFrictionContact_computeGlobalVelocity(problem, reaction, globalVelocity);


    // // // Compute vars without friction coef.
    // // NM_gemv(1.0, P_mu, velocity, 0.0, velocity_tmp);
    // // NM_gemv(1.0, P_mu_inv, reaction, 0.0, reaction_tmp);
    // // // for (unsigned int i = 0; i<n; i++) s[i] = cblas_dnrm2(2, velocity_tmp+i*d+1, 1);
    // // // primalResidual_s(velocity_tmp, H, globalVelocity, w, s, primalConstraint, &pinfeas, options->dparam[SICONOS_DPARAM_TOL]);
    // // dualResidual(M, globalVelocity, H, reaction_tmp, f, dualConstraint, &dinfeas, options->dparam[SICONOS_DPARAM_TOL]);
    // // complem = complemResidualNorm(velocity_tmp, reaction_tmp, nd, n);
    // // udotr = cblas_ddot(nd, velocity_tmp, 1, reaction_tmp, 1);

    // // double norm_q = cblas_dnrm2(m, problem->q, 1);
    // // double norm_b = cblas_dnrm2(nd, problem->b, 1);
    // // double projerr, projerr_comple;
    // // projerr = projectionError(velocity_tmp, reaction_tmp, n, options->dparam[SICONOS_DPARAM_TOL]);


    // // // Store solutions obtained by solving sub-problems
    // // cblas_dcopy(nd, velocity, 1, velocity_tmp, 1);
    // // cblas_dcopy(nd, reaction, 1, reaction_tmp, 1);


    // // printf("\n\n============ Solve the complete problem by gfc3d_IPM_SNM ============\n");
    // // t1 = clock();
    // // gfc3d_IPM_SNM(problem, reaction, velocity, globalVelocity, &info, options);
    // // t2 = clock();
    // // time_comple_prob = (double)(t2-t1)/(double)clk_tck;
    // // printf("Execution time: %.4f,\t\ttest: ", time_comple_prob);
    // // // if (info) printf("failure\n"); else printf("success\n");


    // // // for comparison
    // // cblas_daxpy(nd, -1, velocity, 1, velocity_tmp, 1);
    // // cblas_daxpy(nd, -1, reaction, 1, reaction_tmp, 1);
    // // diff_u = cblas_dnrm2(nd, velocity_tmp, 1);
    // // diff_r = cblas_dnrm2(nd, reaction_tmp, 1);
    // // NM_gemv(1.0, P_mu, velocity, 0.0, velocity_tmp);
    // // NM_gemv(1.0, P_mu_inv, reaction, 0.0, reaction_tmp);
    // // projerr_comple = projectionError(velocity_tmp, reaction_tmp, n, options->dparam[SICONOS_DPARAM_TOL]);



    // // /* For statistiques */
    // // nsc = 0, nN = 0, nB = 0, nR = 0, nT = 0;
    // // for (int i = 0; i < problem->numberOfContacts; i++)
    // // {
    // //   somme[0] = velocity_tmp[3*i] + reaction_tmp[3*i];
    // //   somme[1] = velocity_tmp[3*i+1] + reaction_tmp[3*i+1];
    // //   somme[2] = velocity_tmp[3*i+2] + reaction[3*i+2];
    // //   dur[0] = velocity_tmp[3*i] - reaction_tmp[3*i];
    // //   dur[1] = velocity_tmp[3*i+1] - reaction_tmp[3*i+1];
    // //   dur[2] = velocity_tmp[3*i+2] - reaction_tmp[3*i+2];

    // //   ns = somme[0] - cblas_dnrm2(2,somme+1,1);
    // //   ndur = dur[0] - cblas_dnrm2(2,dur+1,1);
    // //   if (ns > cesp*cblas_dnrm2(3, somme, 1))
    // //   {
    // //     nsc +=1;
    // //     if (dur[0] >= cblas_dnrm2(2,dur+1,1))       { nNcomple_vec[nN++] = i+1; }
    // //     else if (-dur[0] >= cblas_dnrm2(2,dur+1,1)) { nBcomple_vec[nB++] = i+1; }
    // //     else                                        { nRcomple_vec[nR++] = i+1; }
    // //   }
    // //   else
    // //     nTcomple_vec[nT++] = i+1;
    // // }

    // // // fprintf(stats, "stats(end+1) = struct('type', 'complete', 'name', '%s', 'blk_id', %d, 'nBlk', %d, 'mu', %.2e, 'ite', %d, 'residu', %.2e, 'projerr', %.2e,",
    // // //                 strToken_name, -1, n_blocks, problem->mu[0], options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU], projerr_comple);

    // // if (info)
    // // {
    // //   printf("failure\n");
    // //   // fprintf(stats,"'result', 'fail',");
    // // }
    // // else
    // // {
    // //   printf("success\n");
    // //   // fprintf(stats,"'result', 'succ',");
    // // }

    // // // fprintf(stats,"'nd', %d, 'm', %d, 'r_H', %d,", nd, m, -1);
    // // // fprintf(stats,"'B', ["); for (int i=0; i<nB; i++) fprintf(stats, " %d", nBcomple_vec[i]); fprintf(stats,"],");
    // // // fprintf(stats,"'N', ["); for (int i=0; i<nN; i++) fprintf(stats, " %d", nNcomple_vec[i]); fprintf(stats,"],");
    // // // fprintf(stats,"'R', ["); for (int i=0; i<nR; i++) fprintf(stats, " %d", nRcomple_vec[i]); fprintf(stats,"],");
    // // // fprintf(stats,"'T', ["); for (int i=0; i<nT; i++) fprintf(stats, " %d", nTcomple_vec[i]); fprintf(stats,"]);\n");
    // // // fclose(stats);
    // // // ###########################
    // // printf("=====================================================================\n\n");




    // // printf("sumry total:   %s,\t  n nd x m = %d %d x %d\n",strToken_name, n, nd, m);
    // // printf("sumry total: %4d sub-prob, %4d succ, %4d fail,\t %.2e %.2e,\t %.4f (s),\t | u_subs - u_comple | = %.2e\n", n_blocks, n_succ, n_fail, max_residu_subprobs, max_prjerr_subprobs, time_sub_probs, diff_u);

    // // if (info)
    // //   printf("sumry total:   comple prob,   failure           ,");
    // // else
    // //   printf("sumry total:   comple prob,   success           ,");

    // // printf("\t %.2e %.2e,\t %.4f (s),\t | r_subs - r_comple | = %.2e\n",
    // //         options->dparam[SICONOS_DPARAM_RESIDU], projerr_comple, time_comple_prob, diff_r);

    // // printf("sumry total: --------------------------------------------------------------------------------------------------------------------\n\n");



    // // if(H_tilde) {H_tilde = NM_free(H_tilde); H_tilde = NULL;}
    // // if(H) {H = NM_free(H); H = NULL;}
    // // if(P_mu) {P_mu = NM_free(P_mu); P_mu = NULL;}
    // // if(P_mu_inv) {P_mu_inv = NM_free(P_mu_inv); P_mu_inv = NULL;}
    // // if (w) {free(w); w = NULL;}
    // // if (s) {free(s); s = NULL;}
    // // if (velocity_tmp) {free(velocity_tmp); velocity_tmp = NULL;}
    // // if (reaction_tmp) {free(reaction_tmp); reaction_tmp = NULL;}
    // // if (primalConstraint) {free(primalConstraint); primalConstraint = NULL;}
    // // if (dualConstraint) {free(dualConstraint); dualConstraint = NULL;}
    // fclose(varsep);



    // // }// End of TESTING: change of friction coef.
    // free(str);
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
