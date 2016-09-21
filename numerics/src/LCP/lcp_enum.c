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

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LinearComplementarityProblem.h"
#include "LCP_Solvers.h"
#include "lcp_cst.h"
#include "SolverOptions.h"
#include "NumericsMatrix.h"

#include "SiconosLapack.h"
#include "lcp_enum.h"
#include "numerics_verbose.h"

static unsigned long  int sCurrentEnum = 0;
static unsigned long  int sCmpEnum = 0;
static unsigned long  int sNbCase = 0;
static double sProgress = 0;
static int* sWZ = 0;
static double * sM = 0;
static double * sMref = 0;
static double * sQ = 0;
static double * sQref = 0;
static double * sColNul = 0;
static int sSize = 0;
static int LWORK = 0;


static void affectWZ();
static void lcp_buildM(int * zw, double * M, double * Mref, int size);
static void lcp_fillSolution(double*  z, double * w, int size, int* zw, double * Q);
static void lcp_initEnum();
static int lcp_nextEnum();
static void lcp_buildQ();

/*case defined with sCurrentEnum
 *if sWZ[i]==0
 *  w[i] null
 *else
 *  z[i] null
 */
void affectWZ()
{
  unsigned long  int aux = sCurrentEnum;
  for (int i = 0; i < sSize; i++)
  {
    sWZ[i] = aux & 1;
    aux = aux >> 1;
  }
}
void lcp_buildM(int * zw, double * M, double * Mref, int size)
{
  int col;
  double * Aux;
  double * AuxRef;
  Aux = M;
  AuxRef = Mref;
  for (col = 0; col < size; col++)
  {
    if (zw[col] == 0)
    {
      memcpy(Aux, AuxRef, size * sizeof(double));
    }
    else
    {
      /*for(i=0;i<size;i++) Aux[i]=0;*/
      memcpy(Aux, sColNul, sSize * sizeof(double));
      Aux[col] = -1;
      /*M[(n+col)*npm+col+n]=-1;*/
    }
    Aux = Aux + size;
    AuxRef = AuxRef + size;
  }
}
void   lcp_fillSolution(double*  z, double * w, int size, int* zw, double * Q)
{
  int lin;

  for (lin = 0; lin < size; lin++)
  {
    if (zw[lin] == 0)
    {
      w[lin] = 0;
      z[lin] = Q[lin];
    }
    else
    {
      z[lin] = 0;
      w[lin] = Q[lin];
    }
  }
}
void lcp_initEnum()
{
  int cmp;
  sCmpEnum = 0;
  sNbCase = 1;
  for (cmp = 0; cmp < sSize; cmp++)
    sNbCase = sNbCase << 1;
  sProgress = 0;
}
int lcp_nextEnum()
{
  if (sCmpEnum == sNbCase)
    return 0;
  if (sCurrentEnum >= sNbCase)
  {
    sCurrentEnum = 0;
  }
  if (verbose)
    printf("try enum :%d\n", (int)sCurrentEnum);
  affectWZ();
  sCurrentEnum++;
  sCmpEnum++;
  if (verbose && sCmpEnum > (unsigned long int)sProgress * sNbCase)
  {
    sProgress += 0.001;
    printf("lcp_enum progress %f %d \n", sProgress, (int) sCurrentEnum);
  }

  return 1;
}

void lcp_buildQ()
{
  memcpy(sQ, sQref, (sSize)*sizeof(double));
}
int lcp_enum_getNbIWork(LinearComplementarityProblem* problem, SolverOptions* options)
{
  return 2 * (problem->size);
}
int lcp_enum_getNbDWork(LinearComplementarityProblem* problem, SolverOptions* options)
{
  int aux = 3 * (problem->size) + (problem->size) * (problem->size);
  if (options->iparam[4])
  {
    LWORK = -1;
    //int info = 0;
    double dgelsSize = 0;
    //DGELS(problem->M->size0, problem->size , 1, 0, problem->M->size0, 0, problem->M->size0, &dgelsSize, LWORK, &info);
    aux += (int) dgelsSize;
    LWORK = (int) dgelsSize;
  }
  return aux;
}
void lcp_enum_init(LinearComplementarityProblem* problem, SolverOptions* options, int withMemAlloc)
{
  if (withMemAlloc)
  {
    options->dWork = (double *) malloc(lcp_enum_getNbDWork(problem, options) * sizeof(double));
    options->iWork = (int *) malloc(lcp_enum_getNbIWork(problem, options) * sizeof(int));
  }
}
void lcp_enum_reset(LinearComplementarityProblem* problem, SolverOptions* options, int withMemAlloc)
{
  if (withMemAlloc)
  {
    free(options->dWork);
    free(options->iWork);
  }
  options->dWork = NULL;
  solver_options_nullify(options);
}


void lcp_enum(LinearComplementarityProblem* problem, double *z, double *w, int *info , SolverOptions* options)
{
  *info = 1;
  double tol ;
  if (options->dWork == NULL)
  {
    lcp_enum_init(problem, options, 1);
  }
  double * workingFloat = options->dWork;
  int * workingInt = options->iWork;
  int lin;
  sSize = (problem->size);
  int NRHS = 1;
  int * ipiv;
  int check;
  int LAinfo = 0;
  int useDGELS = options->iparam[4];

  /*OUTPUT param*/
  sCurrentEnum = options->iparam[3];
  tol = options->dparam[0];
  int multipleSolutions = options->iparam[0];
  int numberofSolutions = 0;




  sMref = problem->M->matrix0;
  if (!sMref)
  {
    printf("lcp_enum failed, problem->M->matrix0 is null");

  }

  if (verbose)
    printf("lcp_enum begin, size %d tol %e\n", sSize, tol);

  sM = workingFloat;
  sQ = sM + sSize * sSize;
  sColNul = sQ + sSize;
  sQref = sColNul + sSize;
  for (lin = 0; lin < sSize; lin++)
  {
    sQref[lin] =  - problem->q[lin];
    sColNul[lin] = 0;
  }
  sWZ = workingInt;
  ipiv = sWZ + sSize;
  *info = 0;
  lcp_initEnum();
  while (lcp_nextEnum())
  {
    lcp_buildM(sWZ, sM, sMref, sSize);
    lcp_buildQ();
    /*     if (verbose) */
    /*       printCurrentSystem(); */
    if (useDGELS)
    {
      /* if (verbose) */
      /*   { */
      /*     printf("call dgels on ||AX-B||\n"); */
      /*     printf("A\n"); */
      /*     NM_dense_display(sM,sSize,sSize,0); */
      /*     printf("B\n"); */
      /*     NM_dense_display(sQ,sSize,1,0); */
      /*   } */

      DGELS(LA_NOTRANS,sSize, sSize, NRHS, sM, sSize, sQ, sSize,&LAinfo);
      if (verbose)
      {
        printf("Solution of dgels (info=%i)\n", LAinfo);
        NM_dense_display(sQ, sSize, 1, 0);
      }
    }
    else
    {
      DGESV(sSize, NRHS, sM, sSize, ipiv, sQ, sSize, &LAinfo);
    }
    if (!LAinfo)
    {
      if (useDGELS)
      {
        int cc = 0;
        int ii;
        printf("DGELS LAInfo=%i\n", LAinfo);
        for (ii = 0; ii < sSize; ii++)
        {
          if (isnan(sQ[ii]) || isinf(sQ[ii]))
          {
            printf("DGELS FAILED\n");
            cc = 1;
            break;
          }
        }
        if (cc)
          continue;
      }

      if (verbose)
      {
        printf("lcp_enum LU factorization succeeded:\n");
      }

      check = 1;
      for (lin = 0 ; lin < sSize; lin++)
      {
        if (sQ[lin] < - tol)
        {
          check = 0;
          break;/*out of the cone!*/
        }
      }
      if (!check)
        continue;
      else
      {
        numberofSolutions++;
        if (verbose || multipleSolutions)
        {
          printf("lcp_enum find %i solution with sCurrentEnum = %ld!\n", numberofSolutions, sCurrentEnum - 1);
        }


        lcp_fillSolution(z, w, sSize, sWZ, sQ);
        options->iparam[1] = (int) sCurrentEnum - 1;
        options->iparam[2] = numberofSolutions;
        if (!multipleSolutions)  return;
      }
    }
  }
  *info = 1;
  if (verbose)
    printf("lcp_enum has not found a solution!\n");
}

int linearComplementarity_enum_setDefaultSolverOptions(LinearComplementarityProblem* problem, SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the ENUM Solver\n");
  }


  solver_options_nullify(options);
  options->solverId = SICONOS_LCP_ENUM;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  for (i = 0; i < 5; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  if (problem)
  {
    options->dWork = (double*) malloc(lcp_enum_getNbDWork(problem, options) * sizeof(double));
    options->iWork = (int*) malloc(lcp_enum_getNbIWork(problem, options) * sizeof(int));
  }
  else
  {
    options->dWork = NULL;
    options->iWork = NULL;
  }

  options->dparam[0] = 1e-12;



  return 0;
}
