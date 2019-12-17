/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include <stdio.h>                         // for printf
#include <stdlib.h>                        // for calloc, free, malloc
#include "LinearComplementarityProblem.h"  // for LinearComplementarityProblem
#include "NonSmoothDrivers.h"              // for linearComplementarity_driver
#include "NumericsFwd.h"                   // for LinearComplementarityProblem
#include "SolverOptions.h"                 // for SICONOS_DPARAM_RESIDU, Sol...
#include "lcp_test_utils.h"                // for lcp_test_function
#include "test_utils.h"                    // for TestCase
#include "SiconosConfig.h" // for HAVE_GAMS_C_API // IWYU pragma: keep

// --------- GAMS stuff ---------
#ifdef HAVE_GAMS_C_API
#include <string.h>
#if (__linux ||  __APPLE__)
#elif _MSC_VER
#define strdup _strdup
#else // to convert char to const char ... 
static inline char* strdup(char* src)
{
  size_t len = strlen(src) + 1;
  char* dest = (char*)malloc(len * sizeof(char));
  strcpy(dest, src, len);
  return dest;
}
#endif
#endif
// --------- End of GAMS stuff ---------

int lcp_test_function(TestCase * current)
{
  //numerics_set_verbose(2);
  int i, info = 0 ;
  LinearComplementarityProblem* problem = (LinearComplementarityProblem *)malloc(sizeof(LinearComplementarityProblem));
  info = linearComplementarity_newFromFilename(problem, current->filename);
  
#ifdef HAVE_GAMS_C_API
  if (current->options->solverId == SICONOS_LCP_GAMS)
  {
    SN_GAMSparams* GP = (SN_GAMSparams*)current->options->solverParameters;
    assert(GP);
    GP->model_dir = strdup(GAMS_MODELS_SOURCE_DIR);
    assert(current->filename);
    GP->filename = current->filename;
  }
#endif

  double * z = (double *)calloc(problem->size, sizeof(double));
  double * w = (double *)calloc(problem->size, sizeof(double));

  info = linearComplementarity_driver(problem, z , w, current->options);

  for (i = 0 ; i < problem->size ; i++)
  {
    printf("z[%i] = %12.8e\t,w[%i] = %12.8e\n", i, z[i], i, w[i]);
    info = info == 0 ? !(isfinite(z[i]) && isfinite(w[i])): info;
  }

  if (!info)
    printf("test succeeded err = %e \n", current->options->dparam[SICONOS_DPARAM_RESIDU]);
  else
    printf("test unsuccessful err =%e  \n", current->options->dparam[SICONOS_DPARAM_RESIDU]);
  free(z);
  free(w);
  freeLinearComplementarityProblem(problem);
  printf("End of test.\n");
  return info;
}


  
