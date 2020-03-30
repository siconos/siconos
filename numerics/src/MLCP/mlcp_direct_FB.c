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

/*
|A C| |u| |a| |0|
|   |*| |+| |=| |
|D B| |v| |b| |w|
0<z*v>0
dim(u)=mm
dim(v)=nn

**************************************************************************/

#include "mlcp_direct_FB.h"
#include "MLCP_Solvers.h"                       // for mixedLinearComplement...
#include "MixedLinearComplementarityProblem.h"  // for MixedLinearComplement...
#include "SolverOptions.h"                      // for SolverOptions
#include "mlcp_FB.h"                            // for mlcp_FB_getNbDWork
#include "mlcp_direct.h"                        // for mlcp_direct_getNbDWork

static int sN;
static int sM;

static int * siWorkFB = 0;
static int * siWorkDirect = 0;
static double * sdWorkFB = 0;
static double * sdWorkDirect = 0;

int mlcp_direct_FB_getNbIWork(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  int aux = mlcp_FB_getNbIWork(problem, options);
  return mlcp_direct_getNbIWork(problem, options) + aux;
}
int mlcp_direct_FB_getNbDWork(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  int aux = mlcp_FB_getNbDWork(problem, options);
  return mlcp_direct_getNbDWork(problem, options) + aux;
}


void mlcp_direct_FB_init(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  sN = problem->n;
  sM = problem->m;
  int iOffset = mlcp_direct_getNbIWork(problem, options);
  int dOffset = mlcp_direct_getNbDWork(problem, options);
  siWorkFB = options->iWork + iOffset;
  siWorkDirect = options->iWork;
  sdWorkFB = options->dWork + dOffset;
  sdWorkDirect = options->dWork;

  mlcp_direct_init(problem, options);
  options->dWork = sdWorkFB;
  options->iWork = siWorkFB;

  mlcp_FB_init(problem, options);
  options->dWork = sdWorkDirect;
  options->iWork = siWorkDirect;


}
void mlcp_direct_FB_reset()
{
  mlcp_direct_reset();
  mlcp_FB_reset();
}

void mlcp_direct_FB(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{
  /*First, try direct solver*/
  mlcp_direct(problem, z, w, info, options);
  if(*info)
  {
    /*solver direct failed, so run the path solver.*/
    mlcp_FB(problem, z, w, info, options);
    if(!(*info))
    {
      /*       for (i=0;i<problem->n+problem->m;i++){ */
      /*  printf("w[%d]=%f z[%d]=%f\t",i,w[i],i,z[i]);  */
      /*       } */
      mlcp_direct_addConfigFromWSolution(problem, w + sN);
    }
  }
}

