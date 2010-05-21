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
#include <time.h>
#include <math.h>
#include "NumericsOptions.h"
#include "NonSmoothDrivers.h"
#include "LA.h"

/*
 *
 *
 *
 */
double LinearSystemComputeError(LinearSystemProblem* problem, double *z)
{
  double * pM = problem->M->matrix0;
  double * pQ = problem->q;
  double error = 10;
  int n = problem->size;
  double * res = (double*)malloc(n * sizeof(double));
  memcpy(res, pQ, n * sizeof(double));
  DGEMV(LA_NOTRANS, n, n, 1.0, pM, n, z, 1, 1.0, res, 1);
  error = DNRM2(n, res, 1);
  free(res);
  return error;

}
/* return 0 if ok
* otherwise return !=0
*/
int solveLeastSquareProblem(LinearSystemProblem* problem, double *z ,  SolverOptions* options)
{
  int info = -1;
  int n = problem->size;
  int NRHS = 1;

  int n2 = n * n;
  double * Maux = (double*)malloc(n2 * sizeof(double));
  int * ipiv = (int *) malloc(n * sizeof(int));
  int LAinfo;
  /* Checks inputs */
  if (problem == NULL || z == NULL)
    numericsError("EqualityProblem", "null input for EqualityProblem and/or unknowns (z)");
  //displayLS(problem);
  memcpy(Maux, problem->M->matrix0, n2 * sizeof(double));
  memcpy(z, problem->q, n * sizeof(double));


  int LWORK = -1;
  double dgelsSize = 0;
  DGELS(n, n, 1, 0, n, 0, n, &dgelsSize, LWORK, info);
  LWORK = (int) dgelsSize;
  printf("LWORK is :%d\n", LWORK);

  double * dgelsWork = (double *) malloc(LWORK * sizeof(double));

  DGELS(n, n, NRHS, Maux, n, z, n, dgelsWork, LWORK, LAinfo);
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
  printf("LinearSystem_driver: computeError of LinearSystem : %e\n", LinearSystemComputeError(problem, z));
__fin:
  free(Maux);
  free(ipiv);
  free(dgelsWork);
  return info;


}
int LinearSystem_driver_get_dwork(LinearSystemProblem* problem, SolverOptions* options)
{
  return problem->size * problem->size;
}
int LinearSystem_driver_get_iwork(LinearSystemProblem* problem, SolverOptions* options)
{
  return problem->size;
}
int myLu(LinearSystemProblem* problem, double *z ,  SolverOptions* options)
{
  /* Output info. : 0: ok -  >0: problem (depends on solver) */
  int info = -1;
  int n = problem->size;
  int n2 = n * n;

  double * Maux = 0;
  int * ipiv = 0;
  if (options && options->dWork)
    Maux = options->dWork;
  else
    Maux = (double *) malloc(n2 * sizeof(double));

  if (options && options->iWork)
    ipiv = options->iWork;
  else
    ipiv = (int *) malloc(n * sizeof(int));
  int LAinfo;
  /* Checks inputs */
  if (problem == NULL || z == NULL)
    numericsError("EqualityProblem", "null input for EqualityProblem and/or unknowns (z)");
  //displayLS(problem);
  memcpy(Maux, problem->M->matrix0, n2 * sizeof(double));
  memcpy(z, problem->q, n * sizeof(double));

  DGESV(n, 1, Maux, n, ipiv, z, n, LAinfo);
  if (!LAinfo)
  {
    info = 0;
  }
  else
  {
    printf("Equality_driver:: LU foctorization failed:\n");
  }
  for (int i = 0; i < problem->size; i++)
    z[i] = -z[i];

  //printf("LinearSystem_driver: computeError of LinearSystem : %e\n",LinearSystemComputeError(problem,z));

  if (!(options && options->dWork))
    free(Maux);

  if (!(options && options->iWork))
    free(ipiv);

  return info;
}

int LinearSystem_driver(LinearSystemProblem* problem, double *z , double *w, SolverOptions* options)
{
  int i;
  if (problem->M->storageType == 1)
    numericsError("lcp_driver_DenseMatrix", "forbidden type of storage for the matrix M of the LCP");
  for (i = 0; i < problem->size; i++)
    w[i] = 0.0;
  //
  return myLu(problem, z , options);
  //return solveLeastSquareProblem(problem, z , options);
}

void displayLS(LinearSystemProblem* p)
{
  printf("Numerics LinearSystem DISPLAY:\n-------------\n");
  if (!p)
    printf("p is null \n");
  int size = p->size;

  printf("size :%d \n", size);
  if (p->M)
  {
    printf("M matrix:\n");
    displayMat(p->M->matrix0, size, size, size);
  }
  else
    printf("No M matrix:\n");
  if (p->q)
  {
    printf("q matrix:\n");
    displayMat(p->q, size, 1, 0);
  }
  else
    printf("No q:\n");

}
