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

#define _XOPEN_SOURCE 700

#include <math.h>                          // for isfinite
#include <stdio.h>                         // for printf, fclose, fopen, FILE
#include <stdlib.h>                        // for calloc, free
#include "GlobalFrictionContactProblem.h"  // for GlobalFrictionContactProblem
#include "NonSmoothDrivers.h"              // for gfc3d_driver
#include "NumericsFwd.h"                   // for GlobalFrictionContactProblem
#include "NumericsMatrix.h"                // for NumericsMatrix
#include "SolverOptions.h"                 // for SolverOptions, SICONOS_DPA...
#include "frictionContact_test_utils.h"    // for globalFrictionContact_test...
#include "test_utils.h"                    // for TestCase
#include <time.h>
#include "SiconosConfig.h"                 // for HAVE_GAMS_C_API // IWYU pragma: keep
#include <string.h>
#include "SiconosLapack.h"

void print_problem_data_in_Matlab_file(GlobalFrictionContactProblem * problem, FILE * file);

void print_solution_in_Matlab_file(double * v, double * u, double * r, int d, int n, int m, FILE * file);

int globalFrictionContact_test_function(TestCase* current)
{

  int k, info = -1 ;
  GlobalFrictionContactProblem* problem = globalFrictionContact_new_from_filename(current->filename);

  /* globalFrictionContact_display(problem); */

  /* /\* print problem data in Matlab file *\/ */
  /* char *fname; */
  /* FILE * mfile; */
  /* char filename[60] = "cp iterates.m mdata/"; */
  /* char fsolname[60] = "mdata/"; */
  /* char probname[50] = "\0"; */
  /* fname = strrchr(current->filename, '/'); */
  /* fname++; */
  /* for(int i = 0; *fname != '.'; i++) */
  /* { */
  /*   probname[i] = (*fname == '-' ? '_' : *fname); */
  /*   fname++; */
  /* } */
  /* size_t long_probname = strlen(probname); */
  /* strcat(filename, probname); */
  /* strcat(filename, ".m"); */
  /* mfile = fopen(filename, "w"); */
  /* print_problem_data_in_Matlab_file(problem, mfile); */
  /* fclose(mfile); */

  /* strcat(fsolname, probname); */
  /* strcat(fsolname, "_sol.m"); */

  /* FILE * foutput  =  fopen("checkinput.dat", "w"); */
  /* info = globalFrictionContact_printInFile(problem, foutput); */


  int NC = problem->numberOfContacts;
  int dim = problem->dimension;
  int n = problem->M->size1;


  double *reaction = calloc(dim * NC, sizeof(double));
  double *velocity = calloc(dim * NC, sizeof(double));
  double *globalvelocity = calloc(n, sizeof(double));

  long clk_tck = CLOCKS_PER_SEC;


  // --- Extra setup for options when the solver belongs to GAMS family ---
#ifdef HAVE_GAMS_C_API
  // Do we really need this?
  frictionContact_test_gams_opts(current->options);
#endif

  clock_t t1 = clock();

  if(dim == 2)
  {

    info = gfc2d_driver(problem,
                        reaction, velocity, globalvelocity,
                        current->options);
  }
  else if(dim == 3)
  {
    info = gfc3d_driver(problem,
                        reaction, velocity, globalvelocity,
                        current->options);
    //    printf("info = %i\n", info);
  }

  clock_t t2 = clock();

  /* print final solution in Matlab file */
  /* mfile = fopen(fsolname, "w"); */
  /* print_solution_in_Matlab_file(globalvelocity, velocity, reaction, dim, NC, n, mfile); */
  /* fclose(mfile); */

  int print_size = 10;

  printf("Norm velocity:  %12.8e\n", cblas_dnrm2(NC*dim, velocity, 1));
  printf("Norm reaction:  %12.8e\n", cblas_dnrm2(NC*dim, reaction, 1));
  printf("Norm GlobalVe:  %12.8e\n", cblas_dnrm2(n, globalvelocity, 1));

  if(dim * NC >= print_size)
  {
    printf("First values (%i)\n", print_size);
    for(k = 0 ; k < print_size; k++)
    {
      printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k, reaction[k]);
    }
  }
  else
  {
    for(k = 0 ; k < dim * NC; k++)
    {
      printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k, reaction[k]);
    }
  }
  printf(" ..... \n");
  if(n >= print_size)
  {
    printf("First values (%i)\n", print_size);
    for(k = 0 ; k < print_size; k++)
    {
      printf("GlocalVelocity[%i] = %12.8e\n", k, globalvelocity[k]);
    }
  }
  else
  {
    for(k = 0 ; k < n; k++)
    {
      printf("GlocalVelocity[%i] = %12.8e\n", k, globalvelocity[k]);
    }
  }


  for(k = 0; k < dim * NC; ++k)
  {
    info = info == 0 ? !(isfinite(velocity[k]) && isfinite(reaction[k])) : info;
  }

  for(k = 0; k < n; ++k)
  {
    info = info == 0 ? !(isfinite(globalvelocity[k])) : info;
  }

  if(!info)
    printf("test: success\n");
  else
    printf("test: failure\n");

  printf("\nsumry: %d  %9.2e  %5i  %10.4f", info, current->options->dparam[SICONOS_DPARAM_RESIDU], current->options->iparam[SICONOS_IPARAM_ITER_DONE], (double)(t2-t1)/(double)clk_tck);
  printf("%3i %5i %5i     %s\n\n", dim, NC, n, current->filename);

  /* system(filename); */

  free(reaction);
  free(velocity);
  free(globalvelocity);
  /* fclose(foutput); */
  globalFrictionContact_free(problem);
  return info;

}

void print_problem_data_in_Matlab_file(GlobalFrictionContactProblem * problem, FILE * file)
{
  int d = problem->dimension;
  int n = problem->numberOfContacts;
  int m = problem->M->size0;

  fprintf(file,"d = %3i;\n", d);
  fprintf(file,"n = %6i;\n", n);
  fprintf(file,"m = %6i;\n", m);

  fprintf(file,"M = [\n");
  CSparseMatrix_print_in_Matlab_file(NM_triplet(problem->M), 0, file);
  fprintf(file,"    ];\n");
  fprintf(file,"M = sparse(int32(M(:,1)), int32(M(:,2)), M(:,3));\n");

  fprintf(file,"H = [\n");
  CSparseMatrix_print_in_Matlab_file(NM_triplet(problem->H), 0, file);
  fprintf(file,"    ];\n");
  fprintf(file,"H = sparse(int32(H(:,1)), int32(H(:,2)), H(:,3));\n");

  fprintf(file,"q = [");
  for(int i = 0; i < m; i++)
  {
    fprintf(file,"%22.14e; ",problem->q[i]);
  }
  fprintf(file,"];\n");

  fprintf(file,"b = [");
  for(int i = 0; i < n*d; i++)
  {
    fprintf(file,"%22.14e; ",problem->b[i]);
  }
  fprintf(file,"];\n");

  fprintf(file,"mu = [");
  for(int i = 0; i < n; i++)
  {
    fprintf(file,"%22.14e; ",problem->mu[i]);
  }
  fprintf(file,"];\n");

  return;
}


void print_solution_in_Matlab_file(double * v, double * u, double * r, int d, int n, int m, FILE * file)
{
  fprintf(file, "globalVelocity=[\n");
  for(int i = 0; i < m; i++)
  {
    fprintf(file, "%22.15e\n", v[i]);
  }
  fprintf(file, "               ];\n");
  fprintf(file, "velocity=[\n");
  for(int i = 0; i < n*d; i++)
  {
    fprintf(file, "%22.15e\n", u[i]);
  }
  fprintf(file, "         ];\n");
  fprintf(file, "reaction=[\n");
  for(int i = 0; i < n*d; i++)
  {
    fprintf(file, "%22.15e\n", r[i]);
  }
  fprintf(file, "         ];\n");
  return;
}
