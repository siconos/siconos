/* Siconos-Numerics, Copyright INRIA 2005-2010.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "LA.h"
#include <math.h>
#include "LCP_Solvers.h"



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

/*case defined with sCurrentEnum
 *if sWZ[i]==0
 *  w[i] null
 *else
 *  z[i] null
 */
void affectWZ()
{
  unsigned long  int aux = sCurrentEnum;
  for (unsigned int i = 0; i < sSize; i++)
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
  if (verbose && sCmpEnum > sProgress * sNbCase)
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
int lcp_enum_getNbIWork(LinearComplementarity_Problem* problem, Solver_Options* options)
{
  return 2 * (problem->size);
}
int lcp_enum_getNbDWork(LinearComplementarity_Problem* problem, Solver_Options* options)
{
  return 3 * (problem->size) + (problem->size) * (problem->size);
}
void lcp_enum_init(LinearComplementarity_Problem* problem, Solver_Options* options, int withMemAlloc)
{
  if (withMemAlloc)
  {
    options->dWork = (double *) malloc(lcp_enum_getNbDWork(problem, options) * sizeof(double));
    options->iWork = (int *) malloc(lcp_enum_getNbIWork(problem, options) * sizeof(int));
  }
}
void lcp_enum_reset(LinearComplementarity_Problem* problem, Solver_Options* options, int withMemAlloc)
{
  if (withMemAlloc)
  {
    free(options->dWork);
    free(options->iWork);
  }
  options->dWork = NULL;
  options->iWork = NULL;
}


void lcp_enum(LinearComplementarity_Problem* problem, double *z, double *w, int *info , Solver_Options* options)
{
  *info = 1;
  double tol ;
  double * workingFloat = options->dWork;
  int * workingInt = options->iWork;
  int lin;
  sSize = (problem->size);
  int NRHS = 1;
  int * ipiv;
  int check;
  int DGESVinfo;

  /*OUTPUT param*/
  sCurrentEnum = 0;
  tol = options->dparam[0];

  sMref = problem->M->matrix0;
  if (!sMref)
  {
    printf("lcp_enum failed, problem->M->matrix0 is null");

  }

  if (verbose)
    printf("lcp_enum begin, size %d tol %lf\n", sSize, tol);

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
    DGESV(sSize, NRHS, sM, sSize, ipiv, sQ, sSize, DGESVinfo);

    if (!DGESVinfo)
    {
      if (verbose)
      {
        printf("lcp_enum LU foctorization success:\n");
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
        if (verbose)
          printf("lcp_enum find a solution!\n");
        lcp_fillSolution(z, w, sSize, sWZ, sQ);
        options->iparam[1] = sCurrentEnum - 1;
        return;
      }
    }
  }
  *info = 1;
  if (verbose)
    printf("lcp_enum failed!\n");
}

int linearComplementarity_enum_setDefaultSolverOptions(LinearComplementarity_Problem* problem, Solver_Options* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default Solver_Options for the ENUM Solver\n");
  }

  strcpy(options->solverName, "ENUM");

  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = (double*) malloc((3 * problem->size + problem->size * problem->size) * sizeof(double));
  options->iWork = (int*) malloc(2 * problem->size * sizeof(int));
  for (i = 0; i < 5; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->dparam[0] = 1e-12;




  return 0;
}
