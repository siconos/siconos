/* Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr

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
#include "LA.h"
#include "MLCP_Solvers.h"
#include <math.h>
#include "mlcp_direct.h"
#include "mlcp_tool.h"


struct dataComplementarityConf
{
  int * zw; /*zw[i] == 0 if w null and z >=0*/
  double * M;
  struct dataComplementarityConf * next;
  struct dataComplementarityConf * prev;
  int used;
  int* IPV;
};

static double * spCurDouble = 0;
static int * spCurInt = 0;
static double * sQ = 0;
static int sNumberOfCC = 0;
static int sMaxNumberOfCC = 0;
static struct dataComplementarityConf * spFirstCC = 0;
static struct dataComplementarityConf * spCurCC = 0;
static double sTol = 0;
static int sN;
static int sM;
static int sNpM;
static int* spIntBuf;
static int sVerbose = 0;

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
int mlcp_direct_getNbIWork(MixedLinearComplementarity_Problem* problem, Solver_Options* options)
{
  return (problem->n + problem->m) * (options->iparam[5] + 1) + options->iparam[5] * problem->m;
}
int mlcp_direct_getNbDWork(MixedLinearComplementarity_Problem* problem, Solver_Options* options)
{
  return  problem->n + problem->m + (options->iparam[5]) * (problem->n + problem->m) * (problem->n + problem->m);
}

/*
 *options->iparam[5] : n0 number of possible configuration
 *options->dparam[5] : tol
 *options->iWork : double work memory of size (n + m)*(n0+1) + nO*m
 *options->dWork : double work memory of size n + m + n0*(n+m)*(n+m)
 *
 *
 */

void mlcp_direct_init(MixedLinearComplementarity_Problem* problem, Solver_Options* options)
{
  spCurDouble = options->dWork;
  spCurInt = options->iWork;
  sMaxNumberOfCC = options->iparam[5];
  sVerbose = options->iparam[6];
  sTol = options->dparam[5];
  sN = problem->n;
  sM = problem->m;
  sNpM = sN + sM;
  spCurCC = 0;
  spFirstCC = 0;
  sNumberOfCC = 0;
  sQ = mydMalloc(sNpM);
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
/*memory management about floatWorkingMem and intWorkingMem*/
int internalAddConfig(MixedLinearComplementarity_Problem* problem, int * zw, int init)
{
  int i;
  int npm;
  int INFO;
  if (sVerbose)
  {
    printf("mlcp_direct internalAddConfig\n");
    printf("-----------------------------\n");
    for (i = 0; i < problem->m; i++)
      printf("zw[%d]=%d\t", i, zw[i]);
    printf("\n");
  }
  if (init)
  {
    spFirstCC->zw = myiMalloc(sM);
    spFirstCC->IPV = myiMalloc(sNpM);
    spFirstCC->M = mydMalloc(sNpM * sNpM);
  }
  for (i = 0; i < sM; i++)
  {
    spFirstCC->zw[i] = zw[i];
  }
  mlcp_buildM(zw, spFirstCC->M, problem->M->matrix0, sN, sM);
  DGETRF(sNpM, sNpM, spFirstCC->M, sNpM, spFirstCC->IPV, INFO);
  if (INFO)
  {
    printf("mlcp_direct error, LU impossible");
    return 0;
  }
  return 1;
}
/*memory management about dataComplementarityConf*/
void mlcp_direct_addConfig(MixedLinearComplementarity_Problem* problem, int * zw)
{
  if (sNumberOfCC < sMaxNumberOfCC) /*Add a configuration*/
  {
    sNumberOfCC++;
    if (spFirstCC == 0) /*first add*/
    {
      spFirstCC = (struct dataComplementarityConf *) malloc(sizeof(struct dataComplementarityConf));
      spCurCC = spFirstCC;
      spFirstCC->prev = 0;
      spFirstCC->next = 0;
    }
    else
    {
      spFirstCC->prev = (struct dataComplementarityConf *) malloc(sizeof(struct dataComplementarityConf));
      spFirstCC->prev->next = spFirstCC;
      spFirstCC = spFirstCC->prev;
      spFirstCC->prev = 0;
    }
    internalAddConfig(problem, zw, 1);
  }
  else  /*Replace an old one*/
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
    }
    internalAddConfig(problem, zw, 0);
  }
}
void mlcp_direct_addConfigFromWSolution(MixedLinearComplementarity_Problem* problem, double * wSol)
{
  int i;

  for (i = 0; i < sM; i++)
  {
    if (wSol[i] > 0)
      spIntBuf[i] = 1;
    else
      spIntBuf[i] = 0;
  }
  mlcp_direct_addConfig(problem, spIntBuf);
}



int solveWithCurConfig(MixedLinearComplementarity_Problem* problem)
{
  int one = 1;
  int lin;
  int INFO;
  for (lin = 0; lin < sNpM; lin++)
    sQ[lin] =  - problem->q[lin];
  DGETRS(LA_NOTRANS, sNpM, one, spCurCC->M, sNpM, spCurCC->IPV, sQ, sNpM, INFO);
  if (INFO)
  {
    if (sVerbose)
      printf("solveWithCurConfig DGETRS failed\n");
    return 0;
  }
  else
  {
    for (lin = 0 ; lin < sM; lin++)
    {
      if (sQ[sN + lin] < - sTol)
      {
        if (sVerbose)
          printf("solveWithCurConfig Sol not in the positive cone\n");
        return 0;
      }
    }
  }
  if (sVerbose)
    printf("solveWithCurConfig Success\n");
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
void mlcp_direct(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options)
{
  int find = 0;
  if (!spFirstCC)
  {
    (*info) = 1;
  }
  else
  {
    spCurCC = spFirstCC;
    do
    {
      find = solveWithCurConfig(problem);
      if (find)
      {
        mlcp_fillSolution(z, z + sN, w, w + sN, sN, sM, spCurCC->zw, sQ);
        /*Current becomes forst for the next step.*/
        if (spCurCC != spFirstCC)
        {
          spCurCC->prev->next = spCurCC->next;
          spFirstCC->prev = spCurCC;
          spCurCC->next = spFirstCC;
          spFirstCC = spCurCC;
          spFirstCC->prev = 0;
        }
      }
      else
      {
        spCurCC = spCurCC->next;
      }
    }
    while (spCurCC && !find);

    if (find)
      *info = 0;
    else
      *info = 1;
  }
}
