/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

*
* mlcp_direct_addConfigFromWSolution consists in precomputing the linear system from the current solution.
*
* mlcp_direct, the next LCP instance is solved using the current linear system, two cases can happen:
* 1) The complementarity constraints hold --> Success.
* 2) The complementarity constraints don't hold --> Failed.
*
**************************************************************************/

#include "mlcp_direct.h"
#ifndef __cplusplus
#include <stdbool.h>                       // for false
#endif
#include <stdio.h>                              // for printf
#include <stdlib.h>                             // for malloc, exit, free
#include "MLCP_Solvers.h"                       // for mlcp_direct, mlcp_dir...
#include "MixedLinearComplementarityProblem.h"  // for MixedLinearComplement...
#include "NumericsMatrix.h"                     // for NM_dense_display, Num...
#include "SiconosBlas.h"                        // for cblas_dgemv, CblasCol...
#include "SiconosLapack.h"                      // for lapack_int, DGETRF
#include "SolverOptions.h"                      // for SolverOptions
#include "mlcp_cst.h"                           // for SICONOS_IPARAM_MLCP_N...
#include "mlcp_tool.h"                          // for mlcp_buildM, mlcp_fil...
#include "numerics_verbose.h"                   // for verbose
#include "SiconosConfig.h"                      // for DIRECT_SOLVER_USE_DGETRI // IWYU pragma: keep

#define DEBUG_MESSAGES
#include "debug.h"



#define DIRECT_SOLVER_USE_DGETRI
double * sVBuf;

struct dataComplementarityConf
{
  int * zw; /*zw[i] == 0 means w null and z >=0*/
  double * M;
  struct dataComplementarityConf * next;
  struct dataComplementarityConf * prev;
  int Usable;
  int used;
  lapack_int* IPV;
};

static double * sp_curDouble = 0;
static int * sp_curInt = 0;
static double * sQ = 0;
static int s_numberOfCC = 0;
static int s_maxNumberOfCC = 0;
static struct dataComplementarityConf * spFirstCC = 0;
static struct dataComplementarityConf * sp_curCC = 0;
static double sTolneg = 0;
static double sTolpos = 0;
static int s_n;
static int s_m;
static int s_nbLines;
static int s_npM;
static int* spIntBuf;
static int sProblemChanged = 0;

static double * mydMalloc(int n);
static int * myiMalloc(int n);
static int internalPrecompute(MixedLinearComplementarityProblem* problem);
static int internalAddConfig(MixedLinearComplementarityProblem* problem, int * zw, int init);
static int solveWithCurConfig(MixedLinearComplementarityProblem* problem);
//static int nbConfig(struct dataComplementarityConf * pC);

double * mydMalloc(int n)
{
  double * aux = sp_curDouble;
  sp_curDouble = sp_curDouble + n;
  return aux;
}
int * myiMalloc(int n)
{
  int * aux = sp_curInt;
  sp_curInt = sp_curInt + n;
  return aux;
}
// XXX this is going to fail
static lapack_int * myiMalloc2(int n)
{
  lapack_int * aux = (lapack_int*)sp_curInt;
  sp_curInt = (int*)&aux[n];
  return aux;
}

int mlcp_direct_getNbIWork(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  return (problem->n + problem->m) * (options->iparam[SICONOS_IPARAM_MLCP_NUMBER_OF_CONFIGURATIONS] + 1)
    + options->iparam[SICONOS_IPARAM_MLCP_NUMBER_OF_CONFIGURATIONS] * problem->m;
}

int mlcp_direct_getNbDWork(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  return  problem->n + problem->m
    + (options->iparam[SICONOS_IPARAM_MLCP_NUMBER_OF_CONFIGURATIONS]) * ((problem->n + problem->m) * (problem->n + problem->m))
    + (problem->n + problem->m);
}

void mlcp_direct_init(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  sp_curDouble = options->dWork;
  sp_curInt = options->iWork;
  s_maxNumberOfCC = options->iparam[SICONOS_IPARAM_MLCP_NUMBER_OF_CONFIGURATIONS];
  sTolneg = options->dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_NEG];
  sTolpos = options->dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_POS];
  options->iparam[7] = 0;
  sProblemChanged = options->iparam[SICONOS_IPARAM_MLCP_UPDATE_REQUIRED];
  s_n = problem->n;
  s_m = problem->m;
  s_nbLines = problem->n + problem->m;
  if(problem->M->size0 != s_nbLines)
  {
    printf("mlcp_direct_init : M rectangular, not yet managed\n");
    exit(1);
  }

  if(verbose)
    printf("n= %d  m= %d /n sTolneg= %lf sTolpos= %lf \n", s_n, s_m, sTolneg, sTolpos);

  s_npM = s_n + s_m;

  // If the problem comes from the kernel (dynamical systems)
  // Then update is needed but no reset of the previous solutions
  // (This avoids some memory loss by the way)
  if(options->iparam[SICONOS_IPARAM_MLCP_UPDATE_REQUIRED]==0)
  {
    sp_curCC = 0;
    spFirstCC = 0;
    s_numberOfCC = 0;
  }

  sQ = mydMalloc(s_npM);
  sVBuf = mydMalloc(s_npM);
  spIntBuf = myiMalloc(s_npM);
}

void mlcp_direct_reset()
{
  struct dataComplementarityConf * aux;
  while(spFirstCC)
  {
    aux = spFirstCC;
    spFirstCC = spFirstCC->next;
    free(aux);
  }
}

int internalPrecompute(MixedLinearComplementarityProblem* problem)
{
  lapack_int INFO = 0;
  mlcp_buildM(spFirstCC->zw, spFirstCC->M, problem->M->matrix0, s_n, s_m, s_nbLines);
  if(verbose)
  {
    printf("mlcp_direct, precomputed M :\n");
    NM_dense_display(spFirstCC->M, s_npM, s_npM, 0);
  }
  if(!(spFirstCC->Usable))
    return 0;
  DGETRF(s_npM, s_npM, spFirstCC->M, s_npM, spFirstCC->IPV, &INFO);
  if(INFO)
  {
    spFirstCC->Usable = 0;
    printf("mlcp_direct, internalPrecompute  error, LU impossible\n");
    return 0;
  }
#ifdef DIRECT_SOLVER_USE_DGETRI
  DGETRI(s_npM, spFirstCC->M, s_npM, spFirstCC->IPV, &INFO);
  if(INFO)
  {
    spFirstCC->Usable = 1;
    printf("mlcp_direct error, internalPrecompute  DGETRI impossible\n");
    return 0;
  }
#endif
  return 1;
}

/*memory management about floatWorkingMem and intWorkingMem*/
int internalAddConfig(MixedLinearComplementarityProblem* problem, int * zw, int init)
{
  if(verbose)
  {
    printf("mlcp_direct internalAddConfig\n");
    printf("---------\n");
    for(int i = 0; i < problem->m; i++)
      printf("zw[%d]=%d\t", i, zw[i]);
    printf("\n");
  }
  if(init)
  {
    spFirstCC->zw = myiMalloc(s_m);
    spFirstCC->IPV = myiMalloc2(s_npM);
    spFirstCC->M = mydMalloc(s_npM * s_npM);
  }
  for(int i = 0; i < s_m; i++)
  {
    spFirstCC->zw[i] = zw[i];
  }
  return internalPrecompute(problem);
}
/*memory management about dataComplementarityConf*/
void mlcp_direct_addConfig(MixedLinearComplementarityProblem* problem, int * zw)
{
  if(s_numberOfCC < s_maxNumberOfCC)  /*Add a configuration*/
  {
    s_numberOfCC++;
    if(spFirstCC == 0)  /*first add*/
    {
      spFirstCC = (struct dataComplementarityConf *) malloc(sizeof(struct dataComplementarityConf));
      spFirstCC->Usable = 1;
      sp_curCC = spFirstCC;
      spFirstCC->prev = 0;
      spFirstCC->next = 0;
    }
    else
    {
      spFirstCC->prev = (struct dataComplementarityConf *) malloc(sizeof(struct dataComplementarityConf));
      spFirstCC->prev->Usable = 1;
      spFirstCC->prev->next = spFirstCC;
      spFirstCC = spFirstCC->prev;
      spFirstCC->prev = 0;
    }
    internalAddConfig(problem, zw, 1);
  }
  else /*Replace an old one*/
  {
    struct dataComplementarityConf * aux = spFirstCC;
    while(aux->next) aux = aux->next;
    if(aux->prev)
    {
      aux->prev->next = 0;
      spFirstCC->prev = aux;
      aux->prev = 0;
      aux->next = spFirstCC;
      spFirstCC = aux;
      spFirstCC->Usable = 1;
    }
    internalAddConfig(problem, zw, 0);
  }
}
void mlcp_direct_addConfigFromWSolution(MixedLinearComplementarityProblem* problem, double * wSol)
{
  int i;

  for(i = 0; i < s_m; i++)
  {
    if(wSol[i] > sTolpos)
      spIntBuf[i] = 1;
    else
      spIntBuf[i] = 0;
  }
  mlcp_direct_addConfig(problem, spIntBuf);
}



int solveWithCurConfig(MixedLinearComplementarityProblem* problem)
{
  int lin;
  lapack_int INFO = 0;
  double ALPHA = 1;
  double BETA = 0;
  int INCX = 1;
  int INCY = 1;
  double * solTest = 0;
  /*  printf("cur config ");
  for (int i=0;i<s_m;i++)
    if (sp_curCC->zw[i])
      printf("1");
    else
      printf("0");
      printf("\n");*/
  if(sProblemChanged)
    internalPrecompute(problem);
  if(!sp_curCC->Usable)
  {
    if(verbose)
      printf("solveWithCurConfig not usable\n");
    return 0;
  }
#ifdef DIRECT_SOLVER_USE_DGETRI
  cblas_dgemv(CblasColMajor,CblasNoTrans, s_npM, s_npM, ALPHA, sp_curCC->M, s_npM, sQ, INCX, BETA, sVBuf, INCY);
  solTest = sVBuf;
#else
  for(lin = 0; lin < s_npM; lin++)
    sQ[lin] =  - problem->q[lin];
  DGETRS(LA_NOTRANS, s_npM, one, sp_curCC->M, s_npM, sp_curCC->IPV, sQ, s_npM, &INFO);
  solTest = sQ;
#endif
  if(INFO)
  {
    if(verbose)
      printf("solveWithCurConfig DGETRS failed\n");
    return 0;
  }
  else
  {
    for(lin = 0 ; lin < s_m; lin++)
    {
      if(solTest[s_n + lin] < - sTolneg)
      {
        if(verbose)
          printf("solveWithCurConfig Sol not in the positive cone because %lf\n", solTest[s_n + lin]);
        return 0;
      }
    }
  }
  /*  if (verbose)
      printf("solveWithCurConfig Success\n");*/
  return 1;
}

void mlcp_direct(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{
  int find = 0;
  int lin = 0;
  if(!spFirstCC)
  {
    (*info) = 1;
  }
  else
  {
    sp_curCC = spFirstCC;
#ifdef DIRECT_SOLVER_USE_DGETRI
    for(lin = 0; lin < s_npM; lin++)
      sQ[lin] =  - problem->q[lin];
#endif
    do
    {
      find = solveWithCurConfig(problem);
      if(find)
      {
#ifdef DIRECT_SOLVER_USE_DGETRI
        mlcp_fillSolution(z, z + s_n, w, w + s_n, s_n, s_m, s_nbLines, sp_curCC->zw, sVBuf);
#else
        mlcp_fillSolution(z, z + s_n, w, w + s_n, s_n, s_m, s_nbLines, sp_curCC->zw, sQ);
#endif
        /*Current becomes first for the next step.*/
        if(sp_curCC != spFirstCC)
        {
          /*    nbConfig(spFirstCC);
          nbConfig(sp_curCC);
          printf("bidouille pour devenir 1\n");*/
          sp_curCC->prev->next = sp_curCC->next;
          if(sp_curCC->next)
            sp_curCC->next->prev = sp_curCC->prev;
          spFirstCC->prev = sp_curCC;
          sp_curCC->next = spFirstCC;
          spFirstCC = sp_curCC;
          spFirstCC->prev = 0;
          /*    nbConfig(spFirstCC);*/
        }
      }
      else
      {
        sp_curCC = sp_curCC->next;
      }
    }
    while(sp_curCC && !find);

    if(find)
    {
      *info = 0;
    }
    else
    {
      options->iparam[7]++;
      *info = 1;
    }
  }
}

void mlcp_direct_set_default(SolverOptions* options)
{
  options->dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_POS] = 1e-12;
  options->dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_NEG] = 1e-12;
  options->iparam[SICONOS_IPARAM_MLCP_NUMBER_OF_CONFIGURATIONS] = 3;
  options->iparam[SICONOS_IPARAM_MLCP_UPDATE_REQUIRED] = 0;
  options->filterOn = false;
}
