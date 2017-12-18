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

*
* mlcp_direct_addConfigFromWSolution consists in precomputing the linear system from the current solution.
*
* mlcp_direct, the next LCP instance is solved using the current linear system, two cases can happen:
* 1) The complementarity constraints hold --> Success.
* 2) The complementarity constraints don't hold --> Failed.
*
**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MLCP_Solvers.h"
#include "SiconosCompat.h"
#include <math.h>
#include "mlcp_direct.h"
#include "mlcp_tool.h"
#include "SiconosLapack.h"
#include "NumericsMatrix.h"
#include "numerics_verbose.h"

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

static double * spCurDouble = 0;
static int * spCurInt = 0;
static double * sQ = 0;
static int sNumberOfCC = 0;
static int sMaxNumberOfCC = 0;
static struct dataComplementarityConf * spFirstCC = 0;
static struct dataComplementarityConf * spCurCC = 0;
static double sTolneg = 0;
static double sTolpos = 0;
static int sN;
static int sM;
static int sNbLines;
static int sNpM;
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
  double * aux = spCurDouble;
  spCurDouble = spCurDouble + n;
  return aux;
}
int * myiMalloc(int n)
{
  int * aux = spCurInt;
  spCurInt = spCurInt + n;
  return aux;
}
// XXX this is going to fail
static lapack_int * myiMalloc2(int n)
{
  lapack_int * aux = (lapack_int*)spCurInt;
  spCurInt = (int*)&aux[n];
  return aux;
}

int mlcp_direct_getNbIWork(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  return (problem->n + problem->m) * (options->iparam[5] + 1) + options->iparam[5] * problem->m;
}
int mlcp_direct_getNbDWork(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  return  problem->n + problem->m + (options->iparam[5]) * ((problem->n + problem->m) * (problem->n + problem->m)) + (problem->n + problem->m);
}

/*
 *options->iparam[5] : n0 number of possible configuration
 *options->dparam[5] : tol
 *options->iparam[7] : number of failed
 *options->iparam[8] : mlcp problem hab been changed since the previous execution.
 *options->iWork : double work memory of size (n + m)*(n0+1) + nO*m
 *options->dWork : double work memory of size n + m + n0*(n+m)*(n+m)
 *
 *
 */

void mlcp_direct_init(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  spCurDouble = options->dWork;
  spCurInt = options->iWork;
  sMaxNumberOfCC = options->iparam[5];
  sTolneg = options->dparam[5];
  sTolpos = options->dparam[6];
  options->iparam[7] = 0;
  sProblemChanged = options->iparam[8];
  sN = problem->n;
  sM = problem->m;
  sNbLines = problem->n + problem->m;
  if (problem->M->size0 != sNbLines)
  {
    printf("mlcp_direct_init : M rectangular, not yet managed\n");
    exit(1);
  }

  if (verbose)
    printf("n= %d  m= %d /n sTolneg= %lf sTolpos= %lf \n", sN, sM, sTolneg, sTolpos);

  sNpM = sN + sM;
  spCurCC = 0;
  spFirstCC = 0;
  sNumberOfCC = 0;
  sQ = mydMalloc(sNpM);
  sVBuf = mydMalloc(sNpM);
  spIntBuf = myiMalloc(sNpM);
}
void mlcp_direct_reset()
{
  struct dataComplementarityConf * aux;
  while (spFirstCC)
  {
    aux = spFirstCC;
    spFirstCC = spFirstCC->next;
    free(aux);
  }
}
int internalPrecompute(MixedLinearComplementarityProblem* problem)
{
  lapack_int INFO = 0;
  mlcp_buildM(spFirstCC->zw, spFirstCC->M, problem->M->matrix0, sN, sM, sNbLines);
  if (verbose)
  {
    printf("mlcp_direct, precomputed M :\n");
    NM_dense_display(spFirstCC->M, sNpM, sNpM, 0);
  }
  if (!(spFirstCC->Usable))
    return 0;
  DGETRF(sNpM, sNpM, spFirstCC->M, sNpM, spFirstCC->IPV, &INFO);
  if (INFO)
  {
    spFirstCC->Usable = 0;
    printf("mlcp_direct, internalPrecompute  error, LU impossible\n");
    return 0;
  }
#ifdef DIRECT_SOLVER_USE_DGETRI
  DGETRI(sNpM, spFirstCC->M, sNpM, spFirstCC->IPV, &INFO);
  if (INFO)
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
  int i;
  if (verbose)
  {
    printf("mlcp_direct internalAddConfig\n");
    printf("---------\n");
    for (i = 0; i < problem->m; i++)
      printf("zw[%d]=%d\t", i, zw[i]);
    printf("\n");
  }
  if (init)
  {
    spFirstCC->zw = myiMalloc(sM);
    spFirstCC->IPV = myiMalloc2(sNpM);
    spFirstCC->M = mydMalloc(sNpM * sNpM);
  }
  for (i = 0; i < sM; i++)
  {
    spFirstCC->zw[i] = zw[i];
  }
  return internalPrecompute(problem);
}
/*memory management about dataComplementarityConf*/
void mlcp_direct_addConfig(MixedLinearComplementarityProblem* problem, int * zw)
{
  if (sNumberOfCC < sMaxNumberOfCC) /*Add a configuration*/
  {
    sNumberOfCC++;
    if (spFirstCC == 0) /*first add*/
    {
      spFirstCC = (struct dataComplementarityConf *) malloc(sizeof(struct dataComplementarityConf));
      spFirstCC->Usable = 1;
      spCurCC = spFirstCC;
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
    while (aux->next) aux = aux->next;
    if (aux->prev)
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

  for (i = 0; i < sM; i++)
  {
    if (wSol[i] > sTolpos)
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
  for (int i=0;i<sM;i++)
    if (spCurCC->zw[i])
      printf("1");
    else
      printf("0");
      printf("\n");*/
  if (sProblemChanged)
    internalPrecompute(problem);
  if (!spCurCC->Usable)
  {
    if (verbose)
      printf("solveWithCurConfig not usable\n");
    return 0;
  }
#ifdef DIRECT_SOLVER_USE_DGETRI
  cblas_dgemv(CblasColMajor,CblasNoTrans, sNpM, sNpM, ALPHA, spCurCC->M, sNpM, sQ, INCX, BETA, sVBuf, INCY);
  solTest = sVBuf;
#else
  for (lin = 0; lin < sNpM; lin++)
    sQ[lin] =  - problem->q[lin];
  DGETRS(LA_NOTRANS, sNpM, one, spCurCC->M, sNpM, spCurCC->IPV, sQ, sNpM, &INFO);
  solTest = sQ;
#endif
  if (INFO)
  {
    if (verbose)
      printf("solveWithCurConfig DGETRS failed\n");
    return 0;
  }
  else
  {
    for (lin = 0 ; lin < sM; lin++)
    {
      if (solTest[sN + lin] < - sTolneg)
      {
        if (verbose)
          printf("solveWithCurConfig Sol not in the positive cone because %lf\n", solTest[sN + lin]);
        return 0;
      }
    }
  }
  /*  if (verbose)
      printf("solveWithCurConfig Success\n");*/
  return 1;
}
/*
 * The are no memory allocation in mlcp_direct, all necessary memory must be allocated by the user.
 *
 *options:
 * iparam[5] : (in)  n0 number of possible configuration.
 * dparam[5] : (in) a positive value, tolerane about the sign.
 * dWork : working float zone size : n + m + n0*(n+m)*(n+m)  . MUST BE ALLOCATED BY THE USER.
 * iWork : working int zone size : (n + m)*(n0+1) + nO*m. MUST BE ALLOCATED BY THE USER.
 * double *z : size n+m
 * double *w : size n+m
 * info : output. info == 0 if success
 */
/*int nbConfig(struct dataComplementarityConf * pC)
{
  struct dataComplementarityConf * paux = pC;
  int nb = 0;
  while (paux)
  {
    nb++;
    paux = paux->next;
  }
  printf("number of conf :%d\n", nb);
  return nb;

}
*/

void mlcp_direct(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{
  int find = 0;
  int lin = 0;
  if (!spFirstCC)
  {
    (*info) = 1;
  }
  else
  {
    spCurCC = spFirstCC;
#ifdef DIRECT_SOLVER_USE_DGETRI
    for (lin = 0; lin < sNpM; lin++)
      sQ[lin] =  - problem->q[lin];
#endif
    do
    {
      find = solveWithCurConfig(problem);
      if (find)
      {
#ifdef DIRECT_SOLVER_USE_DGETRI
        mlcp_fillSolution(z, z + sN, w, w + sN, sN, sM, sNbLines, spCurCC->zw, sVBuf);
#else
        mlcp_fillSolution(z, z + sN, w, w + sN, sN, sM, sNbLines, spCurCC->zw, sQ);
#endif
        /*Current becomes first for the next step.*/
        if (spCurCC != spFirstCC)
        {
          /*    nbConfig(spFirstCC);
          nbConfig(spCurCC);
          printf("bidouille pour devenir 1\n");*/
          spCurCC->prev->next = spCurCC->next;
          if (spCurCC->next)
            spCurCC->next->prev = spCurCC->prev;
          spFirstCC->prev = spCurCC;
          spCurCC->next = spFirstCC;
          spFirstCC = spCurCC;
          spFirstCC->prev = 0;
          /*    nbConfig(spFirstCC);*/
        }
      }
      else
      {
        spCurCC = spCurCC->next;
      }
    }
    while (spCurCC && !find);

    if (find)
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
