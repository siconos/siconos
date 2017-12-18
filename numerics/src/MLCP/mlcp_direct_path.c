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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MLCP_Solvers.h"
#include "SiconosCompat.h"
#include <math.h>
#include "mlcp_direct_path.h"
#include "mlcp_direct.h"
#include "mlcp_tool.h"

static int sN;
static int sM;

int mixedLinearComplementarity_directPath_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pSolver)
{
  mixedLinearComplementarity_default_setDefaultSolverOptions(problem, pSolver);
  return 0;
}


int mlcp_direct_path_getNbIWork(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  return mlcp_direct_getNbIWork(problem, options); //+mlcp_path_getNbIWork(problem,options);
}
int mlcp_direct_path_getNbDWork(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  return mlcp_direct_getNbDWork(problem, options); //+mlcp_path_getNbDWork(problem,options);
}



/*
 *options->iparam[5] : n0 number of possible configuration.
 * dparam[5] : (in) a positive value, tolerane about the sign.
 *options->iWork : double work memory of  mlcp_direct_path_getNbIWork() integers
 *options->dWork : double work memory of mlcp_direct_path_getNbDWork() doubles
 *
 *
 */

void mlcp_direct_path_init(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  sN = problem->n;
  sM = problem->m;
  mlcp_direct_init(problem, options);
  //mlcp_path_init(problem, options);

}
void mlcp_direct_path_reset()
{
  mlcp_direct_reset();
  //mlcp_path_reset();
}

/*
 * The are no memory allocation in mlcp_direct, all necessary memory must be allocated by the user.
 *
 *options:
 * iparam[0] : (in) verbose.
 * dparam[0] : (in) a positive value, tolerane about the sign used by the path algo.
 * iparam[5] : (in)  n0 number of possible configuration.
 * dparam[5] : (in) a positive value, tolerane about the sign.
 * dWork : working float zone size : n + m + n0*(n+m)*(n+m)  . MUST BE ALLOCATED BY THE USER.
 * iWork : working int zone size : (n + m)*(n0+1) + nO*m. MUST BE ALLOCATED BY THE USER.
 * double *z : size n+m
 * double *w : size n+m
 * info : output. info == 0 if success
 */
void mlcp_direct_path(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{
  /*First, try direct solver*/
  mlcp_direct(problem, z, w, info, options);
  if (*info)
  {
    /*solver direct failed, so run the path solver.*/
    mlcp_path(problem, z, w, info, options);
    if (!(*info))
    {
      /*       for (i=0;i<problem->n+problem->m;i++){ */
      /*  printf("w[%d]=%f z[%d]=%f\t",i,w[i],i,z[i]);  */
      /*       } */
      mlcp_direct_addConfigFromWSolution(problem, w + sN);
    }
  }
}
