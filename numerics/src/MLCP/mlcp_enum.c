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

#include "mlcp_enum.h"
#include <assert.h>                             // for assert
#include <math.h>                               // for isinf, isnan
#ifndef __cplusplus
#include <stdbool.h>                       // for false
#endif
#include <stdio.h>                              // for printf
#include <string.h>                             // for memcpy
#include "MLCP_Solvers.h"                       // for mlcp_compute_error
#include "MixedLinearComplementarityProblem.h"  // for MixedLinearComplement...
#include "NumericsMatrix.h"                     // for NM_dense_display, Num...
#include "SiconosBlas.h"                        // for cblas_dnrm2
#include "SiconosLapack.h"                      // for DGELS, DGESV, lapack_int
#include "SolverOptions.h"                      // for SolverOptions, SICONO...
#include "mlcp_cst.h"                           // for SICONOS_IPARAM_MLCP_E...
#include "mlcp_enum_tool.h"                     // for initEnum, nextEnum
#include "mlcp_tool.h"                          // for mlcp_DisplaySolution
#include "numerics_verbose.h"                   // for verbose
#include "SiconosConfig.h"                      // for MLCP_DEBUG // IWYU pragma: keep

#ifdef MLCP_DEBUG
static int *sLastIWork;
static double *sLastDWork;
#endif

static double * sQ = 0;
static double * sColNul = 0;
static double * sM = 0;
static double * sMref = 0;
static double * sQref = 0;
/* double working memory for dgels*/
static double * sDgelsWork = 0;
static int LWORK = 0;

static int sNn = 0;
static int sMm = 0;
static int sMl = 0;
static int* sW2V = 0;

/*OUTPUT */
/*sW2 is a pointer on the output w*/
static double* sW2;
/*sW1 is a pointer on the output w*/
static double* sW1;
/*sV is a pointer on the output v*/
static double* sV;
/*sU is a pointer on the output u*/
static double* sU;

/** Local, static functions **/
static void buildQ()
{
  memcpy(sQ, sQref, sMl * sizeof(double));
}

static void printCurrentSystem()
{
  int npm = sNn + sMm;
  printf("printCurrentSystemM:\n");
  NM_dense_display(sM, sMl, npm, 0);
  printf("printCurrentSystemQ (ie -Q from mlcp because of linear system MZ=Q):\n");
  NM_dense_display(sQ, sMl, 1, 0);
}

static void printRefSystem()
{
  int npm = sNn + sMm;
  printf("ref M NbLines %d n %d  m %d :\n", sMl, sNn, sMm);
  NM_dense_display(sMref, sMl, npm, 0);
  printf("ref Q (ie -Q from mlcp because of linear system MZ=Q):\n");
  NM_dense_display(sQref, sMl, 1, 0);
}

/* An adaptation of the enum algorithm, to manage the case of MLCP-block formalization
 */
static void mlcp_enum_Block(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{
  double tol ;
  double * workingFloat = options->dWork;
  int * workingInt = options->iWork;
  int lin;
  int npm = (problem->n) + (problem->m);
  int NRHS = 1;
  lapack_int * ipiv;
  int * indexInBlock;
  int check;
  lapack_int LAinfo = 0;
  int useDGELS = options->iparam[SICONOS_IPARAM_MLCP_ENUM_USE_DGELS];
  *info = 0;
  assert(problem->M);
  assert(problem->M->matrix0);
  assert(problem->q);

  sMl = problem->M->size0;
  sNn = problem->n;
  sMm = problem->m;

  /*OUTPUT param*/
  sW1 = w;
  /*sW2=w+(sMl-problem->m); sW2 size :m */
  sU = z;
  tol = options->dparam[SICONOS_DPARAM_TOL];
  int itermax = options->iparam[SICONOS_IPARAM_MAX_ITER];

  sMref = problem->M->matrix0;
  /*  LWORK = 2*npm; LWORK >= max( 1, MN + max( MN, NRHS ) ) where MN = min(M,N)*/
  //  verbose=1;
  if (verbose)
    printf("mlcp_enum begin, n %d m %d tol %lf\n", sNn, sMm, tol);

  sM = workingFloat;
  /*  sQ = sM + npm*npm;*/
  sQ = sM + (sNn + sMm) * sMl;
  /*  sColNul = sQ + sMm +sNn;*/
  sColNul = sQ + sMl;
  /*  sQref = sColNul + sMm +sNn;*/
  sQref = sColNul + sMl;

  sDgelsWork = sQref + sMl;

  for (lin = 0; lin < sMl; lin++)
    sQref[lin] =  - problem->q[lin];
  for (lin = 0; lin < sMl; lin++)
    sColNul[lin] = 0;

  /*  printf("sColNul\n");
      NM_dense_display(sColNul,npm,1);*/
  if (verbose)
    printRefSystem();
  sW2V = workingInt;
  ipiv = sW2V + sMm;
  indexInBlock = ipiv + sMm + sNn;
  if (sMm == 0)
    indexInBlock = 0;
  *info = 0;
  mlcp_buildIndexInBlock(problem, indexInBlock);
  initEnum(problem->m);
  while (nextEnum(sW2V) && itermax-- > 0)
  {
    mlcp_buildM_Block(sW2V, sM, sMref, sNn, sMm, sMl, indexInBlock);
    buildQ();
    if (verbose)
      printCurrentSystem();
    if (useDGELS)
    {
      DGELS(LA_NOTRANS,sMl, npm, NRHS, sM, sMl, sQ, sMl,&LAinfo);
      if (verbose)
      {
        printf("Solution of dgels\n");
        NM_dense_display(sQ, sMl, 1, 0);
      }
    }
    else
    {
      DGESV(npm, NRHS, sM, npm, ipiv, sQ, npm, &LAinfo);
      if (verbose)
      {
        printf("Solution of dgesv\n");
        NM_dense_display(sQ, sMl, 1, 0);
      }
    }
    if (!LAinfo)
    {
      if (useDGELS)
      {
        int cc = 0;
        int ii;
        double rest = 0;
        for (ii = 0; ii < npm; ii++)
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

        if (sMl > npm)
        {
          rest = cblas_dnrm2(sMl - npm, sQ + npm, 1);

          if (rest > tol || isnan(rest) || isinf(rest))
          {
            if (verbose)
              printf("DGELS, optimal point doesn't satisfy AX=b, rest = %e\n", rest);
            continue;
          }
          if (verbose)
            printf("DGELS, optimal point rest = %e\n", rest);
        }
      }

      if (verbose)
      {
        printf("Solving linear system success, solution in cone?\n");
        NM_dense_display(sQ, sMl, 1, 0);
      }

      check = 1;
      for (lin = 0 ; lin < sMm; lin++)
      {
        if (sQ[indexInBlock[lin]] < - tol)
        {
          check = 0;
          break;/*out of the cone!*/
        }
      }
      if (!check)
        continue;
      else
      {
        double err;
        mlcp_fillSolution_Block(sU, sW1, sNn, sMm, sMl, sW2V, sQ, indexInBlock);
        mlcp_compute_error(problem, z, w, tol, &err);
        /*because it happens the LU leads to an wrong solution witout raise any error.*/
        if (err > 10 * tol)
        {
          if (verbose)
            printf("LU no-error, but mlcp_compute_error out of tol: %e!\n", err);
          continue;
        }
        if (verbose)
        {
          printf("mlcp_enum find a solution err = %e!\n", err);
          mlcp_DisplaySolution_Block(sU, sW1, sNn, sMm, sMl, indexInBlock);
        }
        // options->iparam[1]=sCurrentEnum-1;
        return;
      }
    }
    else
    {
      if (verbose)
      {
        printf("LU factorization failed:\n");
      }
    }
  }
  *info = 1;
  if (verbose)
    printf("mlcp_enum failed!\n");
}


/** End of static functions **/

int mlcp_enum_getNbIWork(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  if (!problem)
    return 0;
  return 2 * (problem->n + problem->m) + problem->m;
}

int mlcp_enum_getNbDWork(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  if (!problem)
    return 0;
  assert(problem->M);
  LWORK = 0;
  if (options->iparam[SICONOS_IPARAM_MLCP_ENUM_USE_DGELS])
  {
    LWORK = -1;
    //int info = 0;
    double dgelsSize = 0;
    //DGELS(problem->M->size0, problem->n + problem->m, 1, 0, problem->M->size0, 0, problem->M->size0, &dgelsSize, LWORK, &info);
    LWORK = (int) dgelsSize;
  }
  return LWORK + 3 * (problem->M->size0) + (problem->n + problem->m) * (problem->M->size0); 
}

void mlcp_enum(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{
  int nbSol = 0;
  if (problem->blocksRows)
  {
    mlcp_enum_Block(problem, z, w, info, options);
    return;
  }
  double tol ;
  double * workingFloat = options->dWork;
  int * workingInt = options->iWork;
  int lin;
  int npm = (problem->n) + (problem->m);
  int NRHS = 1;
  lapack_int * ipiv;
  int check;
  lapack_int LAinfo = 0;
  *info = 0;
  sMl = problem->M->size0;
  sNn = problem->n;
  sMm = problem->m;
  int useDGELS = options->iparam[SICONOS_IPARAM_MLCP_ENUM_USE_DGELS];
  /*OUTPUT param*/
  sW1 = w;
  sW2 = w + (sMl - problem->m); /*sW2 size :m */
  sU = z;
  sV = z + problem->n;
  tol = options->dparam[SICONOS_DPARAM_TOL];
  int itermax = options->iparam[SICONOS_IPARAM_MAX_ITER];

  sMref = problem->M->matrix0;
  /*  LWORK = 2*npm; LWORK >= max( 1, MN + max( MN, NRHS ) ) where MN = min(M,N)*/
  //  verbose=1;
  if (verbose)
    printf("mlcp_enum begin, n %d m %d tol %lf\n", sNn, sMm, tol);

  sM = workingFloat;
  /*  sQ = sM + npm*npm;*/
  sQ = sM + (sNn + sMm) * sMl;
  /*  sColNul = sQ + sMm +sNn;*/
  sColNul = sQ + sMl;
  /*  sQref = sColNul + sMm +sNn;*/
  sQref = sColNul + sMl;

  sDgelsWork = sQref + sMl;

  for (lin = 0; lin < sMl; lin++)
    sQref[lin] =  - problem->q[lin];
  for (lin = 0; lin < sMl; lin++)
    sColNul[lin] = 0;
  /*  printf("sColNul\n");
      NM_dense_display(sColNul,npm,1);*/
  if (verbose)
    printRefSystem();
  sW2V = workingInt;
  ipiv = sW2V + sMm;

  initEnum(problem->m);
  while (nextEnum(sW2V) && itermax-- > 0)
  {
    mlcp_buildM(sW2V, sM, sMref, sNn, sMm, sMl);
    buildQ();
    if (verbose)
      printCurrentSystem();
    if (useDGELS)
    {
      DGELS(LA_NOTRANS,sMl, npm, NRHS, sM, sMl, sQ, sMl, &LAinfo);
      if (verbose)
      {
        printf("Solution of dgels\n");
        NM_dense_display(sQ, sMl, 1, 0);
      }
    }
    else
    {
      DGESV(npm, NRHS, sM, npm, ipiv, sQ, npm, &LAinfo);
      if (verbose)
      {
        printf("Solution of dgesv\n");
        NM_dense_display(sQ, sMl, 1, 0);
      }
    }
    if (!LAinfo)
    {
      if (useDGELS)
      {
        int cc = 0;
        int ii;
        double rest = 0;
        for (ii = 0; ii < npm; ii++)
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

        if (sMl > npm)
        {
          rest = cblas_dnrm2(sMl - npm, sQ + npm, 1);

          if (rest > tol || isnan(rest) || isinf(rest))
          {
            if (verbose)
              printf("DGELS, optimal point doesn't satisfy AX=b, rest = %e\n", rest);
            continue;
          }
          if (verbose)
            printf("DGELS, optimal point rest = %e\n", rest);
        }
      }

      if (verbose)
      {
        printf("Solving linear system success, solution in cone?\n");
        NM_dense_display(sQ, sMl, 1, 0);
      }

      check = 1;
      for (lin = 0 ; lin < sMm; lin++)
      {
        if (sQ[sNn + lin] < - tol)
        {
          check = 0;
          break;/*out of the cone!*/
        }
      }
      if (!check)
        continue;
      else
      {
        double err;
        mlcp_fillSolution(sU, sV, sW1, sW2, sNn, sMm, sMl, sW2V, sQ);
        mlcp_compute_error(problem, z, w, tol, &err);
        /*because it happens the LU leads to an wrong solution witout raise any error.*/
        if (err > 10 * tol)
        {
          if (verbose)
            printf("LU no-error, but mlcp_compute_error out of tol: %e!\n", err);
          continue;
        }
        nbSol++;
        if (verbose)
        {
          printf("mlcp_enum find a solution, err=%e !\n", err);
          mlcp_DisplaySolution(sU, sV, sW1, sW2, sNn, sMm, sMl);
        }
        return;
      }
    }
    else
    {
      if (verbose)
      {
        printf("LU factorization failed:\n");
      }
    }
  }
  *info = 1;
  if (verbose)
    printf("mlcp_enum failed nbSol=%i!\n", nbSol);
}

void mlcp_enum_set_default(SolverOptions* options)
{
  options->dparam[SICONOS_IPARAM_MLCP_ENUM_USE_DGELS] = 0;
  options->filterOn = false;

}


