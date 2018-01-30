/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include <string.h>

#if (__linux ||  __APPLE__)
#elif _MSC_VER
#define strdup _strdup
#else
static inline char* strdup(char* src)
{
  size_t len = strlen(src) + 1;
  char* dest = (char*)malloc(len * sizeof(char));
  strcpy(dest, src, len);
  return dest;
}
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "NonSmoothDrivers.h"
#include "lcp_test_function.h"
#include "GAMSlink.h"
#include "lcp_cst.h"
#include "LCP_Solvers.h"
#include "LinearComplementarityProblem.h"
#include "SolverOptions.h"
#include "SiconosCompat.h"

#ifdef __cplusplus
using namespace std;
#endif

int lcp_test_function(FILE * f, int solverId, char* filename)
{

  int i, info = 0 ;
  LinearComplementarityProblem* problem = (LinearComplementarityProblem *)malloc(sizeof(LinearComplementarityProblem));

  info = linearComplementarity_newFromFile(problem, f);

  FILE * foutput  =  fopen("./lcp_mmc.verif", "w");
  info = linearComplementarity_printInFile(problem, foutput);
  fclose(foutput);
  SolverOptions options;
  solver_options_set(&options, solverId);

#ifdef HAVE_GAMS_C_API
  if (solverId == SICONOS_LCP_GAMS)
  {
    SN_GAMSparams* GP = (SN_GAMSparams*)options.solverParameters;
    assert(GP);
    GP->model_dir = strdup(GAMS_MODELS_SOURCE_DIR);
    assert(filename);
    GP->filename = filename;
  }
#endif

  double * z = (double *)calloc(problem->size, sizeof(double));
  double * w = (double *)calloc(problem->size, sizeof(double));

  info = linearComplementarity_driver(problem, z , w, &options);

  for (i = 0 ; i < problem->size ; i++)
  {
    printf("z[%i] = %12.8e\t,w[%i] = %12.8e\n", i, z[i], i, w[i]);
    info = info == 0 ? !(isfinite(z[i]) && isfinite(w[i])): info;
  }


  if (!info)
  {
    printf("test succeeded err = %e \n", options.dparam[1]);
  }
  else
  {
    printf("test unsuccessful err =%e  \n", options.dparam[1]);
  }
  free(z);
  free(w);

  solver_options_delete(&options);

  if (solverId == SICONOS_LCP_GAMS)
  {
    free(options.solverParameters);
    options.solverParameters = NULL;
  }

  freeLinearComplementarityProblem(problem);
  printf("End of test.\n");


  return info;
}

int lcp_test_function_SBM(FILE * f, int solverId)
{

  int i, info = 0 ;
  LinearComplementarityProblem* problem = (LinearComplementarityProblem *)malloc(sizeof(LinearComplementarityProblem));

  info = linearComplementarity_newFromFile(problem, f);

  FILE * foutput  =  fopen("./lcp_mmc.verif", "w");
  info = linearComplementarity_printInFile(problem, foutput);
  SolverOptions * options = (SolverOptions *)malloc(sizeof(SolverOptions));



  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_NSGS_SBM);

  solver_options_set(options->internalSolvers, solverId);

#ifdef HAVE_GAMS_C_API
  if (solverId == SICONOS_LCP_GAMS)
  {
    // no testing for now
    solver_options_delete(options);
    free(options);
    freeLinearComplementarityProblem(problem);
    fclose(foutput);
    return 0;
  }
#endif



  double * z = (double *)calloc(problem->size, sizeof(double));
  double * w = (double *)calloc(problem->size, sizeof(double));

  info = linearComplementarity_driver(problem, z , w, options);

  for (i = 0 ; i < problem->size ; i++)
  {
    printf("z[%i] = %12.8e\t,w[%i] = %12.8e\n", i, z[i], i, w[i]);
  }

  if (!info)
  {
    printf("test succeeded err=%e \n", options->dparam[1]);
  }
  else
  {
    printf("test unsuccessful err =%e \n", options->dparam[1]);
  }
  free(z);
  free(w);
  // info = linearComplementarity_deleteDefaultSolverOptions(&options,solvername);

  solver_options_delete(options);
  free(options);

  freeLinearComplementarityProblem(problem);
  fclose(foutput);

  return info;


}


