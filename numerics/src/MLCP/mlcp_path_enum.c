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

/*
|A C| |u| |a| |0|
|   |*| |+| |=| |
|D B| |v| |b| |w|
0<z*v>0
dim(u)=mm
dim(v)=nn

**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MLCP_Solvers.h"
#include "SiconosCompat.h"
#include <math.h>
#include "mlcp_path_enum.h"
#include "mlcp_enum.h"
#include "mlcp_tool.h"

static int sN;
static int sM;

static int * siWorkEnum = 0;
static int * siWorkPath = 0;
static double * sdWorkEnum = 0;
static double * sdWorkPath = 0;

int mixedLinearComplementarity_pathEnum_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pSolver)
{
  mixedLinearComplementarity_default_setDefaultSolverOptions(problem, pSolver);
  return 0;
}


int mlcp_path_enum_getNbIWork(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  return mlcp_enum_getNbIWork(problem, options);
}
int mlcp_path_enum_getNbDWork(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  return mlcp_enum_getNbDWork(problem, options);
}



/*
 *options->iparam[5] : n0 number of possible configuration.
 * dparam[5] : (in) a positive value, tolerane about the sign.
 *options->iWork : double work memory of  mlcp_direct_enum_getNbIWork() integers  (2(nn+mm))+((n + m)*(n0+1) + nO*m)
 *options->dWork : double work memory of mlcp_direct_enum_getNbDWork() doubles  ((nn+mm)*(nn+mm) + 3*(nn+mm))+(n + m + n0*(n+m)*(n+m))
 *
 *
 */

void mlcp_path_enum_init(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  sN = problem->n;
  sM = problem->m;
  int iOffset = 0;/* mlcp_path_getNbIWork(problem,options);*/
  int dOffset = 0;/*mlcp_path_getNbDWork(problem,options);*/
  siWorkEnum = options->iWork + iOffset;
  siWorkPath = options->iWork;
  sdWorkEnum = options->dWork + dOffset;
  sdWorkPath = options->dWork;
  /*  mlcp_path_init(problem, options);*/

}
void mlcp_path_enum_reset()
{
  /*mlcp_path_reset();*/
  siWorkEnum = 0;
  siWorkPath = 0;
  sdWorkEnum = 0;
  sdWorkPath = 0;
}

/*
 * The are no memory allocation in mlcp_direct, all necessary memory must be allocated by the user.
 *
 *options:
 * iparam[0] : (in) verbose.
 * dparam[0] : (in) a positive value, tolerane about the sign used by the enum algo.
 * iparam[5] : (in)  n0 number of possible configuration.
 * dparam[5] : (in) a positive value, tolerane about the sign.
 * dWork : working float zone size : n + m + n0*(n+m)*(n+m)  . MUST BE ALLOCATED BY THE USER.
 * iWork : working int zone size : (n + m)*(n0+1) + nO*m. MUST BE ALLOCATED BY THE USER.
 * double *z : size n+m
 * double *w : size n+m
 * info : output. info == 0 if success
 */
void mlcp_path_enum(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{
  if (!siWorkEnum)
  {
    *info = 1;
    printf("MLCP_PATH_ENUM error, call a non initialised method!!!!!!!!!!!!!!!!!!!!!\n");
    return;
  }
  /*First, try direct solver*/
  //  options->dWork = sdWorkDirect;
  //  options->iWork = siWorkDirect;
  if (options->iparam[8] == 761727)
  {
    mixedLinearComplementarity_display(problem);
  }
  mlcp_path(problem, z, w, info, options);
  if (*info)
  {
    printf("MLCP_PATH_ENUM: path failed, call enum\n");
    options->dWork = sdWorkEnum;
    options->iWork = siWorkEnum;
    /*solver direct failed, so run the enum solver.*/
    mlcp_enum(problem, z, w, info, options);
  }
}
