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
#include <time.h>
#include <math.h>
#include "LinearSystemProblem.h"
#include "SolverOptions.h"
#include "NonSmoothDrivers.h"
#include "NumericsMatrix.h"
#include "SiconosLapack.h"
#include "numerics_verbose.h"

static int LWORK = 0;
int SICONOS_LS_0 = 0;

static int solveLeastSquareProblem(LinearSystemProblem* problem, double *z ,  SolverOptions* options);
static int myLu(LinearSystemProblem* problem, double *z ,  SolverOptions* options);
#ifdef LINEARSYSTEM_DEBUG
static void displayLS(LinearSystemProblem* p);
#endif
//#define LINEARSYSTEM_DEBUG

/*
 *check M z + q =0
 *
 *
 */
double LinearSystem_computeError(LinearSystemProblem* problem, double *z)
{
  double * pM = problem->M->matrix0;
  double * pQ = problem->q;
  double error = 10;
  int n = problem->size;
  double * res = (double*)malloc(n * sizeof(double));
  memcpy(res, pQ, n * sizeof(double));
  cblas_dgemv(CblasColMajor,CblasNoTrans, n, n, 1.0, pM, n, z, 1, 1.0, res, 1);
  error = cblas_dnrm2(n, res, 1);
  free(res);
  return error;

}
/* return 0 if ok
* otherwise return !=0
*/
int solveLeastSquareProblem(LinearSystemProblem* problem, double *z ,  SolverOptions* options)
{


  /* Checks inputs */
  if (problem == NULL || z == NULL)
    numerics_error("EqualityProblem", "null input for EqualityProblem and/or unknowns (z)");

  /* Output info. : 0: ok -  >0: problem (depends on solver) */
  int info = -1;
  int n = problem->size;
  int n2 = n * n;

  double * Maux = 0;
  int * ipiv = 0;
  if (options && options->dWork)
    Maux = options->dWork;
  else
    Maux = (double *) malloc(LinearSystem_getNbDwork(problem, options) * sizeof(double));

  if (options && options->iWork)
    ipiv = options->iWork;
  else
    ipiv = (int *) malloc(LinearSystem_getNbIwork(problem, options) * sizeof(int));
  int LAinfo = 0;
  //displayLS(problem);

  assert(problem->M->matrix0);
  assert(problem->q);

  memcpy(Maux, problem->M->matrix0, n2 * sizeof(double));
  //  memcpy(z,problem->q,n*sizeof(double));
  for (int ii = 0; ii < n; ii++)
    z[ii] = -problem->q[ii];


  printf("LinearSystem : solveLeastSquareProblem LWORK is :%d\n", LWORK);


  DGELS(LA_NOTRANS,n, n, 1, Maux, n, z, n, &LAinfo);
  if (LAinfo)
  {
    printf("LinearSystem_driver: DGELS  failed:\n");
    goto __fin;
  }

  int ii;
  for (ii = 0; ii < n; ii++)
  {
    if (isnan(z[ii]) || isinf(z[ii]))
    {
      printf("DGELS FAILED\n");
      goto __fin;
    }
  }
  info = 0;
  printf("LinearSystem_driver: computeError of LinearSystem : %e\n", LinearSystem_computeError(problem, z));
__fin:
  if (!(options && options->dWork))
    free(Maux);

  if (!(options && options->iWork))
    free(ipiv);

  return info;


}
int LinearSystem_getNbDwork(LinearSystemProblem* problem, SolverOptions* options)
{
  int aux = problem->size * problem->size;
  if (options && options->iparam[4])
  {
    LWORK = -1;
    //int info = 0;
    double dgelsSize = 0;
    //DGELS(problem->size, problem->size , 1, 0, problem->size, 0, problem->size, &dgelsSize, LWORK, &info);
    aux += (int) dgelsSize;
    LWORK = (int) dgelsSize;
  }
  return aux;
}
int LinearSystem_getNbIwork(LinearSystemProblem* problem, SolverOptions* options)
{
  return problem->size;
}
int myLu(LinearSystemProblem* problem, double *z ,  SolverOptions* options)
{
  /* Checks inputs */
  if (problem == NULL || z == NULL)
    numerics_error("EqualityProblem", "null input for EqualityProblem and/or unknowns (z)");
  /* Output info. : 0: ok -  >0: problem (depends on solver) */
  int info = -1;
  int n = problem->size;
  int n2 = n * n;

  double * Maux = 0;
  int * ipiv = 0;
  if (options && options->dWork)
    Maux = options->dWork;
  else
    Maux = (double *) malloc(LinearSystem_getNbDwork(problem, options) * sizeof(double));

  if (options && options->iWork)
    ipiv = options->iWork;
  else
    ipiv = (int *) malloc(LinearSystem_getNbIwork(problem, options) * sizeof(int));
  int LAinfo = 0;
  //displayLS(problem);

  assert(problem->M->matrix0);
  assert(problem->q);

  memcpy(Maux, problem->M->matrix0, n2 * sizeof(double));
  //  memcpy(z,problem->q,n*sizeof(double));
  for (int ii = 0; ii < n; ii++)
    z[ii] = -problem->q[ii];

  DGESV(n, 1, Maux, n, ipiv, z, n, &LAinfo);
  if (!LAinfo)
  {
    info = 0;
  }
  else
  {
    printf("Equality_driver:: LU factorization failed:\n");
  }

  //printf("LinearSystem_driver: computeError of LinearSystem : %e\n",LinearSystemComputeError(problem,z));

  if (!(options && options->dWork))
    free(Maux);

  if (!(options && options->iWork))
    free(ipiv);

  return info;
}
/*
 * options->iWork containing the int work memory.
 * options->dWork containing the double work memory.
 * options->iparam[4] : use DGELS (1) or DGESV (0).
 */
int LinearSystem_driver(LinearSystemProblem* problem, double *z , double *w, SolverOptions* options)
{
  int i;
  assert(problem->M);
  if (problem->M->storageType == 1)
    numerics_error("LinearSystem_driver", "forbidden type of storage for the matrix M of the LCP");
#ifdef LINEARSYSTEM_DEBUG
  displayLS(problem);
#endif
  for (i = 0; i < problem->size; i++)
    w[i] = 0.0;

  int res = 0;
  if (options && options->iparam[4])
  {
    res = solveLeastSquareProblem(problem, z , options);
  }
  else
  {
    res = myLu(problem, z , options);
  }
  if (options && options->filterOn)
  {
    options->dparam[1] = LinearSystem_computeError(problem, z);
    printf("LinearSystem_driver solved with error = %e\n", LinearSystem_computeError(problem, z));
    printf("The solution is :\n");
    for (i = 0; i < problem->size; i++)
      printf(" %e", z[i]);
    printf("\n");
  }
  return res;
}
#ifdef LINEARSYSTEM_DEBUG
void displayLS(LinearSystemProblem* p)
{
  printf("Numerics LinearSystem DISPLAY:\n-------------\n");
  if (!p)
  {
    printf("p is null \n");
    return;
  }
  int size = p->size;

  printf("size :%d \n", size);
  if (p->M)
  {
    printf("M matrix:\n");
    NM_dense_display(p->M->matrix0, size, size, size);
  }
  else
    printf("No M matrix:\n");
  if (p->q)
  {
    printf("q matrix:\n");
    NM_dense_display(p->q, size, 1, 0);
  }
  else
    printf("No q:\n");

}
#endif
int LinearSystem_newFromFile(LinearSystemProblem* problem, FILE* file)
{
  int n = 0;
  int i;

  CHECK_IO(fscanf(file, "%d\n", &n));
  problem->size = n;
  problem->M = newNumericsMatrix();

  newFromFile(problem->M, file);

  problem->q = (double *) malloc(problem->M->size1 * sizeof(double));
  for (i = 0; i < problem->M->size1; i++)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->q[i])));
  }
  return 1;
}
int LinearSystem_setDefaultSolverOptions(LinearSystemProblem* problem, SolverOptions* options, int solverId)
{
  int info = 0;
  if (verbose > 0)
    printf("Set the Default SolverOptions for the LS Solver\n");


  options->solverId = SICONOS_LS_0;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  /*use dgels ?*/
  options->iparam[4] = 0;
  if (problem)
  {
    options->dWork = (double*) malloc(LinearSystem_getNbDwork(problem, options) * sizeof(double));
    options->iWork = (int*) malloc(LinearSystem_getNbIwork(problem, options) * sizeof(int));
  }
  else
  {
    options->dWork = NULL;
    options->iWork = NULL;
  }
  options->dparam[0] = 1e-12;
  return info;

}
void LinearSystem_freeProblem(LinearSystemProblem *problem)
{
  freeNumericsMatrix(problem->M);
  if (problem->M)
    free(problem->M);
  if (problem->q)
    free(problem->q);
  free(problem);
}
