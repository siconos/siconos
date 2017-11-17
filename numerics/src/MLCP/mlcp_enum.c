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
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MLCP_Solvers.h"
#include "SiconosCompat.h"
#include "mlcp_tool.h"
#include "NumericsMatrix.h"
#include <math.h>
#include "mlcp_enum.h"
#include "mlcp_enum_tool.h"
#include "SiconosLapack.h"
#include "numerics_verbose.h"

//#ifdef HAVE_DGELS
//#define ENUM_USE_DGELS
//#endif

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

static void buildQ(void);
static void printCurrentSystem(void);
static void printRefSystem(void);

/*case defined with sCurrentEnum
 *if sW2V[i]==0
 *  v[i] not null w2[i] null
 *else
 *  v[i] null and w2[i] not null
 */
/* void affectW2V(){ */
/*   unsigned long  int aux = sCurrentEnum; */
/*   for (unsigned int i =0; i <sMm; i++){ */
/*     sW2V[i]=aux & 1; */
/*     aux = aux >> 1; */
/*   } */

/* } */
/* void initEnum(){ */
/*   int cmp; */
/*   sCmpEnum=0;   */
/*   sNbCase = 1; */
/*   for (cmp=0;cmp<sMm;cmp++) */
/*     sNbCase = sNbCase << 1; */
/*   sProgress=0; */
/* } */
/* int nextEnum(){ */
/*   if (sCmpEnum == sNbCase) */
/*     return 0; */
/*   if (sCurrentEnum >= sNbCase){ */
/*     sCurrentEnum=0; */
/*   } */
/*   if (verbose) */
/*     printf("try enum :%d\n",(int)sCurrentEnum); */
/*   affectW2V(); */
/*   sCurrentEnum++; */
/*   sCmpEnum++; */
/*   if (verbose && sCmpEnum > sProgress*sNbCase){ */
/*     sProgress+=0.001; */
/*     printf(" progress %f %d \n",sProgress,(int) sCurrentEnum); */
/*   } */

/*   return 1; */
/* } */


int mixedLinearComplementarity_enum_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pOPtionSolver)
{
  mixedLinearComplementarity_default_setDefaultSolverOptions(problem, pOPtionSolver);
  pOPtionSolver->dparam[0] = 1e-12;
  return 0;
}


void buildQ()
{
  memcpy(sQ, sQref, sMl * sizeof(double));
}
void printCurrentSystem()
{
  int npm = sNn + sMm;
  printf("printCurrentSystemM:\n");
  NM_dense_display(sM, sMl, npm, 0);
  printf("printCurrentSystemQ (ie -Q from mlcp because of linear system MZ=Q):\n");
  NM_dense_display(sQ, sMl, 1, 0);
}
void printRefSystem()
{
  int npm = sNn + sMm;
  printf("ref M NbLines %d n %d  m %d :\n", sMl, sNn, sMm);
  NM_dense_display(sMref, sMl, npm, 0);
  printf("ref Q (ie -Q from mlcp because of linear system MZ=Q):\n");
  NM_dense_display(sQref, sMl, 1, 0);
}
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
  if (options->iparam[4])
  {
    LWORK = -1;
    //int info = 0;
    double dgelsSize = 0;
    //DGELS(problem->M->size0, problem->n + problem->m, 1, 0, problem->M->size0, 0, problem->M->size0, &dgelsSize, LWORK, &info);
    LWORK = (int) dgelsSize;
  }
  return LWORK + 3 * (problem->M->size0) + (problem->n + problem->m) * (problem->M->size0); 
}

/*
 * The are no memory allocation in mlcp_enum, all necessary memory must be allocated by the user.
 *
 *options:
 * dparam[0] : (in) a positive value, tolerane about the sign.
 * dWork : working float zone size : (nn+mm)*(nn+mm) + 3*(nn+mm). MUST BE ALLOCATED BY THE USER.
 * iWork : working int zone size : 2(nn+mm). MUST BE ALLOCATED BY THE USER.
 * double *z : size n+m
 * double *w : size n+m
 */
void mlcp_enum_Block(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);
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
  int useDGELS = options->iparam[4];
  /*OUTPUT param*/
  sW1 = w;
  sW2 = w + (sMl - problem->m); /*sW2 size :m */
  sU = z;
  sV = z + problem->n;
  tol = options->dparam[0];
  int itermax = options->iparam[0];

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
    printf("mlcp_enum failed nbSol=%i!\n", nbSol);
}
/*
An adaptation of the previuos algorithm, to manage the case of MLCP-block formalization
 */
void mlcp_enum_Block(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
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
  int useDGELS = options->iparam[4];
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
  tol = options->dparam[0];
  int itermax = options->iparam[0];

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
int mlcp_enum_alloc_working_memory(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  if (options->iWork || options->dWork)
    return 0;
  options->iWork = (int *) malloc(mlcp_enum_getNbIWork(problem, options) * sizeof(int));
  options->dWork = (double *) malloc(mlcp_enum_getNbDWork(problem, options) * sizeof(double)); 
  return 1;
}
void mlcp_enum_free_working_memory(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  if (options->iWork)
    free(options->iWork);
  if (options->dWork)
    free(options->dWork);
  options->iWork = NULL;
  options->dWork = NULL;
}
