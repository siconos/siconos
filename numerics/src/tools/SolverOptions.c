/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.

 * Copyright 2016 INRIA.

 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at

 * http://www.apache.org/licenses/LICENSE-2.0

 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#include "SolverOptions.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "mlcp_cst.h"
#include "MCP_cst.h"
#include "NCP_cst.h"
#include "lcp_cst.h"
#include "relay_cst.h"
#include "Friction_cst.h"
#include "AVI_cst.h"
#include "VI_cst.h"
#include "misc.h"

#include "Newton_Methods.h"
#include "PathSearch.h"
#include "VariationalInequality_Solvers.h"

#include "GAMSlink.h"

#include "SiconosNumerics_Solvers.h"

//#define DEBUG_MESSAGES 1
#include "debug.h"

#define MAX_ENV_SIZE 200

char * SICONOS_NUMERICS_PROBLEM_LCP_STR = "LCP";
char * SICONOS_NUMERICS_PROBLEM_MLCP_STR = "MLCP";
char * SICONOS_NUMERICS_PROBLEM_NCP_STR = "NCP";
char * SICONOS_NUMERICS_PROBLEM_MCP_STR = "MCP";
char * SICONOS_NUMERICS_PROBLEM_EQUALITY_STR = "EQUALITY";
char * SICONOS_NUMERICS_PROBLEM_FC2D_STR = "FC2D";
char * SICONOS_NUMERICS_PROBLEM_FC3D_STR = "FC3D";
char * SICONOS_NUMERICS_PROBLEM_VI_STR = "VI";
char * SICONOS_NUMERICS_PROBLEM_AVI_STR = "AVI";

static void recursive_printSolverOptions(SolverOptions* options, int level);

char * idProblemToChar(int id)
{
  switch (id)
  {
  case (SICONOS_NUMERICS_PROBLEM_LCP):
  {
    return SICONOS_NUMERICS_PROBLEM_LCP_STR;
  }
  case (SICONOS_NUMERICS_PROBLEM_MLCP):
  {
    return SICONOS_NUMERICS_PROBLEM_MLCP_STR;
  }
  case (SICONOS_NUMERICS_PROBLEM_NCP):
  {
    return SICONOS_NUMERICS_PROBLEM_NCP_STR;
  }
  case (SICONOS_NUMERICS_PROBLEM_MCP):
  {
    return SICONOS_NUMERICS_PROBLEM_MCP_STR;
  }
  case (SICONOS_NUMERICS_PROBLEM_EQUALITY):
  {
    return SICONOS_NUMERICS_PROBLEM_EQUALITY_STR;
  }
  case (SICONOS_NUMERICS_PROBLEM_FC2D):
  {
    return SICONOS_NUMERICS_PROBLEM_FC2D_STR;
  }
  case (SICONOS_NUMERICS_PROBLEM_FC3D):
  {
    return SICONOS_NUMERICS_PROBLEM_FC3D_STR;
  }
  case (SICONOS_NUMERICS_PROBLEM_VI):
  {
    return SICONOS_NUMERICS_PROBLEM_VI_STR;
  }
  case (SICONOS_NUMERICS_PROBLEM_AVI):
  {
    return SICONOS_NUMERICS_PROBLEM_AVI_STR;
  }
  default:
    printf("Numerics:idProblemToChar, id unknown : %d \n", id);
    return NULL;
  }

}

void readSolverOptions(int driverType, SolverOptions* options)
{
  /* To each problem, corresponds a XXX_parameters.opt file where default parameters can be read, XXX being the problem name (LCP, FrictionContact3D ...) */

  if (verbose > 0)
    printf("\n ========== Numerics Non Smooth Solver - Read default parameters for the solver.\n ==========");

  // Checks if NUMERICSSPATH is set.
  if (getenv("SICONOSPATH") == NULL)
  {
    fprintf(stderr, "Numerics, readSolverOptions error, SICONOSPATH environment variable not set. Can not find default solver options file.\n");
    return;
    //exit(EXIT_FAILURE);
  }

  FILE * ficin;
  /* Name of the default parameters file */
  char name[MAX_ENV_SIZE];
  char nameSuffix[] = "/include/siconos";

  strncpy(name, getenv("SICONOSPATH"), MAX_ENV_SIZE - sizeof(nameSuffix));
  strcat(name, nameSuffix);

  char buffer[64];
  char bufferName[64];
  /* Return value for reading */
  int nval;

  // set default size to 4 ...
  if (options->iparam == NULL)
    options->iparam = (int*)malloc(4 * sizeof(*options->iparam));
  if (options->dparam == NULL)
    options->dparam = (double*)malloc(4 * sizeof(*options->dparam));

  switch (driverType)
  {

  case 0:
    strcat(name, "LCP_parameters.opt");
    break;
  case 1:
    strcat(name, "dfc2D_parameters.opt");
    break;
  case 2:
    strcat(name, "FrictionContact2D_parameters.opt");
    break;
  case 3:
    strcat(name, "FrictionContact3D_parameters.opt");
    ficin = fopen(name, "r");
    if (verbose > 0)
      printf("The default-parameters file is: %s\n", name);
    if (!ficin)
    {
      printf("Numerics, readSolverOptions error. Can not open file %60s", name);
      exit(-1);
    }
    //nval = fscanf(ficin, "%c", &(options->solverName));
    CHECK_IO(fgets(buffer, 64, ficin));
    CHECK_IO(fgets(buffer, 64, ficin));
    CHECK_IO(fgets(buffer, 64, ficin));
    /* Solver name */
    CHECK_IO(fgets(bufferName , 64, ficin));
    options->solverId = nameToId(bufferName);
    CHECK_IO(fgets(buffer, 64, ficin));
    /* iparam */
    nval = fscanf(ficin, "%d%d", &(options->iparam[0]), &(options->iparam[1]));
    if (nval != 4)
    {
      fprintf(stderr, "Numerics, readSolverOptions error, wrong number of parameters for iparam.\n");
      exit(EXIT_FAILURE);

    }
    /* dparam */
    nval = fscanf(ficin, "%lf%lf%lf", &(options->dparam[0]), &(options->dparam[1]), &(options->dparam[2]));
    if (nval != 3)
    {
      fprintf(stderr, "Numerics, readSolverOptions error, wrong number of parameters for dparam.\n");
      exit(EXIT_FAILURE);

    }
    fclose(ficin);
    break;
  default:
    fprintf(stderr, "Numerics, readSolverOptions error, unknown problem type.\n");
    exit(EXIT_FAILURE);
  }
}

void recursive_printSolverOptions(SolverOptions* options, int level)
{
  char* marge;
  int i;
  marge = (char*) malloc((level + 1) * sizeof(char));
  for (i = 0; i < level; i++)
    marge[i] = ' ';
  marge[level] = '\0';

  printf("%s\n ========== Numerics Non Smooth Solver parameters: \n", marge);
  if (options->isSet == 0)
    printf("%sThe solver parameters have not been set. \t options->isSet = %i \n", marge, options->isSet);
  else
  {
    printf("%sThe solver parameters below have  been set \t options->isSet = %i\n", marge, options->isSet);
    printf("%sId of the solver\t\t\t\t options->solverId = %d \n", marge, options->solverId);
    printf("%sName of the solver\t\t\t\t  %s \n", marge, idToName(options->solverId));
    if (options->iparam != NULL)
    {
      printf("%sint parameters \t\t\t\t\t options->iparam\n", marge);
      printf("%ssize of the int parameters\t\t\t options->iSize = %i\n", marge, options->iSize);
      for (int i = 0; i < options->iSize; ++i)
        printf("%s\t\t\t\t\t\t options->iparam[%i] = %d\n", marge, i, options->iparam[i]);
    }
    if (options->dparam != NULL)
    {
      printf("%sdouble parameters \t\t\t\t options->dparam\n", marge);
      printf("%ssize of the double parameters\t\t\t options->dSize = %i\n", marge, options->dSize);
      for (int i = 0; i < options->dSize; ++i)
        printf("%s\t\t\t\t\t\t options->dparam[%i] = %.6le\n", marge, i, options->dparam[i]);
    }
  }
  if (options->iWork == NULL)
  {
    printf("%sinteger work array have not been allocated. \t options->iWork = NULL \n", marge);
  }
  else
  {
    printf("%sinteger work array have been allocated. \t options->iWork = %p \n", marge, options->iWork);
    printf("%sinteger work array size \t\t\t options->iSize = %i \n", marge, options->iSize);
  }
  if (options->dWork == NULL)
  {
    printf("%sdouble work array have not been allocated. \t options->dWork = NULL \n", marge);
  }
  else
  {
    printf("%sdouble work array have been allocated. \t options->dWork = %p \n", marge, options->dWork);
    printf("%sdouble work array size \t\t\t options->dSize = %i \n", marge, options->dSize);
  }




  printf("%sSee %s documentation for parameters definition)\n", marge, idToName(options->solverId));

  printf("\n");

  printf("%snumber of internal (or local) solvers \t\t options->numberOfInternalSolvers = %i\n", marge, options->numberOfInternalSolvers);
  for (i = 0; i < options->numberOfInternalSolvers; i++)
  {
    recursive_printSolverOptions(options->internalSolvers + i, level + 1);
  }
  free(marge);

}
void printSolverOptions(SolverOptions* options)
{
  recursive_printSolverOptions(options, 0);
}

void free_solver_specific_data(SolverOptions* options)
{
  int id = options->solverId;
  switch (id)
  {
    case SICONOS_NCP_PATHSEARCH:
      assert(options->solverData);
      free_solverData_PathSearch(options->solverData);
      free(options->solverData);
      options->solverData = NULL;
      break;
    case SICONOS_FRICTION_3D_GAMS_PATH:
    case SICONOS_FRICTION_3D_GAMS_PATHVI:
    case SICONOS_FRICTION_3D_GAMS_LCP_PATH:
    case SICONOS_FRICTION_3D_GAMS_LCP_PATHVI:
    case SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH:
    case SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI:
    {
      deleteGAMSparams((SN_GAMSparams *)options->solverParameters);
      options->solverParameters = NULL;
      break;
    }
    case SICONOS_VI_BOX_AVI_LSA:
    {
     vi_box_AVI_free_solverData(options);
     break;
    }
    default:
      {
       if (options->solverParameters)
       {
         free(options->solverParameters);
         options->solverParameters = NULL;
       }
      }
  }

  if (newton_LSA_check_solverId(id))
  {
    newton_LSA_free_solverOptions(options);
  }
}

void deleteSolverOptions(SolverOptions* op)
{
  if(op)
  {
    assert(op->numberOfInternalSolvers >= 0);
    for (int i = 0; i < op->numberOfInternalSolvers; i++)
      deleteSolverOptions(&(op->internalSolvers[i]));
    if (op->numberOfInternalSolvers && op->internalSolvers)
      free(op->internalSolvers);
    op->internalSolvers = NULL;
    if (op->iparam)
      free(op->iparam);
    op->iparam = NULL;
    if (op->dparam)
      free(op->dparam);
    op->dparam = NULL;
    if (op->iWork)
      free(op->iWork);
    op->iWork = NULL;
    if (op->dWork)
      free(op->dWork);
    op->dWork = NULL;
    if (op->callback)
    {
      // MB: Yoyo & BeadPlan ok now.
      free(op->callback);
      op->callback = NULL;
    }
    free_solver_specific_data(op);
  }
}

void null_SolverOptions(SolverOptions* options)
{
  options->dWork = NULL;
  options->iWork = NULL;
  options->callback = NULL;
  options->numericsOptions = NULL;
  options->internalSolvers = NULL;
  options->solverData = NULL;
  options->solverParameters = NULL;
}

void fill_SolverOptions(SolverOptions* options, int solverId, int iSize, int dSize, int iter_max, double tol)
{
  options->solverId = solverId;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = iSize;
  options->dSize = dSize;
  options->iparam = (int *)calloc(iSize, sizeof(int));
  options->dparam = (double *)calloc(dSize, sizeof(double));
  null_SolverOptions(options);
  /* we set those value, even if they don't make sense. If this is the case,
   * they should be +inf */
  options->iparam[0] = iter_max;
  options->dparam[0] = tol;

}

void set_SolverOptions(SolverOptions* options, int solverId)
{
 int iSize = 0;
 int dSize = 0;
 int iter_max = 0;
 double tol = 0.0;

 switch (solverId)
 {
  case SICONOS_LCP_PATH:
  {
    tol = 1e-12;
  }
  case SICONOS_LCP_PGS:
  case SICONOS_LCP_CPG:
  case SICONOS_LCP_NEWTONMIN:
  case SICONOS_LCP_QP:
  case SICONOS_LCP_NSQP:
  {
    iSize = 2;
    dSize = 2;
    iter_max = 1000;
    tol = tol == 0. ? 1e-8 : tol;
    fill_SolverOptions(options, solverId, iSize, dSize, iter_max, tol);
    break;
  }
  case SICONOS_LCP_RPGS:
  {
    iSize = 2;
    dSize = 3;
    iter_max = 1000;
    tol = 1e-8;
    fill_SolverOptions(options, solverId, iSize, dSize, iter_max, tol);
    options->dparam[2] = 1.0;
    break;
  }
  case SICONOS_LCP_LATIN:
  {
    iSize = 2;
    dSize = 3;
    iter_max = 1000;
    tol = 1e-8;
    fill_SolverOptions(options, solverId, iSize, dSize, iter_max, tol);
    options->dparam[2] = 0.3;
    break;
  }
  case SICONOS_LCP_LATIN_W:
  {
    iSize = 2;
    dSize = 4;
    iter_max = 1000;
    tol = 1e-8;
    fill_SolverOptions(options, solverId, iSize, dSize, iter_max, tol);
    options->dparam[2] = 0.3;
    options->dparam[3] = 1.0;
    break;
  }
  case SICONOS_LCP_ENUM:
  {
    iSize = 5;
    dSize = 2;
    iter_max = 0; /* this indicates the number of solutions ...  */
    tol = 1e-8;
    fill_SolverOptions(options, solverId, iSize, dSize, iter_max, tol);
    /*use dgels:*/
    options->iparam[4] = 0;
    break;
  }

  case SICONOS_NCP_NEWTON_FBLSA:
  case SICONOS_NCP_NEWTON_MINFBLSA:
  case SICONOS_MCP_NEWTON_FBLSA:
  case SICONOS_MCP_NEWTON_MINFBLSA:
  case SICONOS_LCP_NEWTON_FBLSA:
  case SICONOS_LCP_NEWTON_MINFBLSA:
  case SICONOS_VI_BOX_QI:
    iSize = 6;
    dSize = 3;
    iter_max = 1000;
    tol = 1e-12;
    fill_SolverOptions(options, solverId, iSize, dSize, iter_max, tol);
    newton_lsa_default_SolverOption(options);
    break;

  case SICONOS_NCP_PATHSEARCH:
    iSize = 9;
    dSize = 8;
    iter_max = 100;
    tol = 1e-12;
    fill_SolverOptions(options, solverId, iSize, dSize, iter_max, tol);
    pathsearch_default_SolverOption(options);
    break;

  case SICONOS_VI_BOX_AVI_LSA:
    iSize = 6;
    dSize = 3;
    iter_max = 100;
    tol = 1e-12;
    fill_SolverOptions(options, solverId, iSize, dSize, iter_max, tol);
    vi_box_AVI_extra_SolverOptions(options);
    break;

  case SICONOS_AVI_CAOFERRIS:
  case SICONOS_LCP_AVI_CAOFERRIS:
  case SICONOS_RELAY_AVI_CAOFERRIS:
  case SICONOS_RELAY_AVI_CAOFERRIS_TEST:
    iSize = 6;
    dSize = 3;
    iter_max = 10000;
    tol = 1e-12;
    fill_SolverOptions(options, solverId, iSize, dSize, iter_max, tol);
    break;

  case SICONOS_LCP_LEMKE:
  case SICONOS_LCP_BARD:
  case SICONOS_LCP_MURTY:
  case SICONOS_LCP_PIVOT:
  case SICONOS_LCP_PIVOT_LUMOD:
  case SICONOS_LCP_PATHSEARCH:
    iSize = 6;
    dSize = 3;
    iter_max = 10000;
    tol = 100*DBL_EPSILON;
    fill_SolverOptions(options, solverId, iSize, dSize, iter_max, tol);
    options->iparam[SICONOS_IPARAM_PIVOT_RULE] = SICONOS_LCP_PIVOT_LEMKE;
    break;

  case SICONOS_LCP_GAMS:
  {
    tol = 1e-12;
  }
  case SICONOS_FRICTION_3D_GAMS_PATH:
  case SICONOS_FRICTION_3D_GAMS_PATHVI:
  case SICONOS_FRICTION_3D_GAMS_LCP_PATH:
  case SICONOS_FRICTION_3D_GAMS_LCP_PATHVI:
  case SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH:
  case SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI:
  {
#ifdef HAVE_GAMS_C_API
    iSize = 5;
    dSize = 4; // stupid thing in checkTrivialCase in fc3d_driver.c
    iter_max = 10000;
    tol = tol == 0. ? 1e-9 : tol;
    fill_SolverOptions(options, solverId, iSize, dSize, iter_max, tol);
    if (!options->solverParameters)
    {
      options->solverParameters = createGAMSparams(GAMS_MODELS_SHARE_DIR, GAMS_DIR);
    }
    break;
#else
    printf("set_SolverOptions :: GAMS was not enabled, exiting!\n");
    exit(EXIT_FAILURE);
#endif
  }

  case SICONOS_NCP_PATH:
  case SICONOS_VI_BOX_PATH:
    iSize = 6;
    dSize = 3;
    iter_max = 10000;
    tol = 1e-12;
    fill_SolverOptions(options, solverId, iSize, dSize, iter_max, tol);
    break;
  default:
   printf("set_SolverOptions not supported for solver id %d named %s\n", solverId, idToName(solverId));
   exit(EXIT_FAILURE);
 }

}

char * idToName(int Id)
{
  switch (Id)
  {

#undef SICONOS_SOLVER_MACRO
#define SICONOS_SOLVER_MACRO(X) case X: return X ## _STR ;
SICONOS_REGISTER_SOLVERS()
  default:
    return SICONOS_NONAME_STR;
  }
}
int nameToId(char * pName)
{
#undef SICONOS_SOLVER_MACRO
#define SICONOS_SOLVER_MACRO(X) if (strcmp(X ## _STR, pName) == 0) return X;
SICONOS_REGISTER_SOLVERS()
  return 0;

}
