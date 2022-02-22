/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#include "lcp_enum.h"
#include <math.h>                          // for isinf, isnan
#include <stdlib.h>                        // for malloc, free
#include <string.h>                        // for NULL, memcpy
#include "LCP_Solvers.h"                   // for lcp_enum, lcp_enum_init
#include "LinearComplementarityProblem.h"  // for LinearComplementarityProblem
#include "NumericsMatrix.h"                // for NM_dense_display, Numerics...
#include "SiconosLapack.h"                 // for DGELS, DGESV, lapack_int, LA_NOTRANS
#include "SolverOptions.h"                 // for SolverOptions, solver_opti...
#include "lcp_cst.h"                       // for SICONOS_LCP_IPARAM_ENUM_US...
#include "numerics_verbose.h"              // for numerics_printf, verbose
#include "enum_tool.h"


static void lcp_buildM(int * zw,
                       double * M,
                       double * Mref,
                       int size,
                       double * column_of_zero)
{
  int col;
  double * Aux;
  double * AuxRef;
  Aux = M;
  AuxRef = Mref;
  for(col = 0; col < size; col++)
  {
    if(zw[col] == 0)
    {
      memcpy(Aux, AuxRef, size * sizeof(double));
    }
    else
    {
      /*for(i=0;i<size;i++) Aux[i]=0;*/
      memcpy(Aux, column_of_zero, size * sizeof(double));
      Aux[col] = -1;
      /*M[(n+col)*npm+col+n]=-1;*/
    }
    Aux = Aux + size;
    AuxRef = AuxRef + size;
  }
}
static void   lcp_fillSolution(double*  z, double * w, int size, int* zw, double * Q)
{
  int lin;

  for(lin = 0; lin < size; lin++)
  {
    if(zw[lin] == 0)
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

int lcp_enum_getNbIWork(LinearComplementarityProblem* problem, SolverOptions* options)
{
  return 2 * (problem->size);
}
int lcp_enum_getNbDWork(LinearComplementarityProblem* problem, SolverOptions* options)
{
  int aux = 3 * (problem->size) + (problem->size) * (problem->size);
  if(options->iparam[SICONOS_LCP_IPARAM_ENUM_USE_DGELS])
  {
    // to be reviewed ...
    /* int LWORK = -1; */
    /* int info = 0; */
    /* double dgelsSize = 0; */
    /* DGELS(problem->M->size0, problem->size , 1, 0, problem->M->size0, 0, problem->M->size0, &dgelsSize, LWORK, &info); */
    /* aux += (int) dgelsSize; */
    /* LWORK = (int) dgelsSize; */
  }
  return aux;
}
void lcp_enum_init(LinearComplementarityProblem* problem, SolverOptions* options, int withMemAlloc)
{
  if(withMemAlloc)
  {
    options->dWork = (double *) malloc(lcp_enum_getNbDWork(problem, options) * sizeof(double));
    options->iWork = (int *) malloc(lcp_enum_getNbIWork(problem, options) * sizeof(int));
  }
}
void lcp_enum_reset(LinearComplementarityProblem* problem, SolverOptions* options, int withMemAlloc)
{
  if(withMemAlloc)
  {
    if(options->dWork)
      free(options->dWork);
    if(options->iWork)
      free(options->iWork);
  }
  options->dWork = NULL;
  options->iWork = NULL;
}


void lcp_enum(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{
  *info = 1;
  double tol ;
  if(options->dWork == NULL)
  {
    lcp_enum_init(problem, options, 1);
  }
  double * workingFloat = options->dWork;
  int * workingInt = options->iWork;

  int size = problem->size;
  int NRHS = 1;
  lapack_int * ipiv;
  int check;
  lapack_int LAinfo = 0;
  int useDGELS = options->iparam[SICONOS_LCP_IPARAM_ENUM_USE_DGELS];

  /*OUTPUT param*/

  tol = options->dparam[SICONOS_DPARAM_TOL];
  int multipleSolutions = options->iparam[SICONOS_LCP_IPARAM_ENUM_MULTIPLE_SOLUTIONS];
  int numberofSolutions = 0;


  if(! problem->M->matrix0)
  {
    numerics_printf("lcp_enum failed, problem->M->matrix0 is null");

  }

  if(verbose)
    numerics_printf("lcp_enum begin, size %d tol %e", size, tol);

  double * M_linear_system = workingFloat;
  double * q_linear_system =  M_linear_system + size * size;
  double * column_of_zero = q_linear_system + size;
  double * q_linear_systemref = column_of_zero + size;

  for(int row = 0; row < size; row++)
  {
    q_linear_systemref[row] =  - problem->q[row];
    column_of_zero[row] = 0;
  }
  //sWZ = workingInt;
  int * zw_indices  =  workingInt;
  ipiv = zw_indices + size;
  *info = 0;
  EnumerationStruct * enum_struct = enum_init(size);
  enum_struct->current = options->iparam[SICONOS_LCP_IPARAM_ENUM_SEED];
  while(enum_next(zw_indices, size, enum_struct))
  {
    lcp_buildM(zw_indices,  M_linear_system, problem->M->matrix0, size, column_of_zero);
    memcpy(q_linear_system, q_linear_systemref, (size)*sizeof(double));
    /*     if (verbose) */
    /*       printCurrentSystem(); */
    if(useDGELS)
    {
      /* if (verbose) */
      /*   { */
      /*     numerics_printf("call dgels on ||AX-B||\n"); */
      /*     numerics_printf("A\n"); */
      /*     NM_dense_display( M_linear_system,sSize,sSize,0); */
      /*     numerics_printf("B\n"); */
      /*     NM_dense_display(q_linear_system,sSize,1,0); */
      /*   } */

      DGELS(LA_NOTRANS, size, size, NRHS,  M_linear_system, size, q_linear_system, size, &LAinfo);
      if(verbose)
      {
        numerics_printf("Solution of dgels (info=%i)", LAinfo);
        NM_dense_display(q_linear_system, size, 1, 0);
      }
    }
    else
    {
      DGESV(size, NRHS,  M_linear_system, size, ipiv, q_linear_system, size, &LAinfo);
    }
    if(!LAinfo)
    {
      if(useDGELS)
      {
        int cc = 0;
        int ii;
        numerics_printf("DGELS LAInfo=%i", LAinfo);
        for(ii = 0; ii < size; ii++)
        {
          if(isnan(q_linear_system[ii]) || isinf(q_linear_system[ii]))
          {
            numerics_printf("DGELS FAILED");
            cc = 1;
            break;
          }
        }
        if(cc)
          continue;
      }

      if(verbose)
      {
        numerics_printf("lcp_enum LU factorization succeeded:");
      }

      check = 1;
      for(int row  = 0 ; row < size; row++)
      {
        if(q_linear_system[row] < - tol)
        {
          check = 0;
          break;/*out of the cone!*/
        }
      }
      if(!check)
        continue;
      else
      {
        numberofSolutions++;
        if(verbose || multipleSolutions)
        {
          numerics_printf("lcp_enum find %i solution with scurrent = %ld!", numberofSolutions, enum_struct->current - 1);
        }
        *info = 0;
        lcp_fillSolution(z, w, size, zw_indices, q_linear_system);
        options->iparam[SICONOS_LCP_IPARAM_ENUM_CURRENT_ENUM ] = (int) enum_struct->current - 1;
        options->iparam[SICONOS_LCP_IPARAM_ENUM_NUMBER_OF_SOLUTIONS] = numberofSolutions;
        if(!multipleSolutions)  return;
      }
    }
  }
  *info = 1;
  if(verbose)
    numerics_printf("lcp_enum has not found a solution!\n");
}

void lcp_enum_set_default(SolverOptions* options)
{
  options->iparam[SICONOS_LCP_IPARAM_ENUM_USE_DGELS] = 0;
  options->iparam[SICONOS_LCP_IPARAM_ENUM_SEED] = 0;
  options->iparam[SICONOS_LCP_IPARAM_ENUM_MULTIPLE_SOLUTIONS] = 0;
  // SICONOS_LCP_IPARAM_ENUM_CURRENT_ENUM (out)
  // SICONOS_LCP_IPARAM_ENUM_NUMBER_OF_SOLUTIONS (out)

}
