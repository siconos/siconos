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
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "Numerics_Options.h"
#include "NonSmoothDrivers.h"

/*
 *
 *
 *
 */

/* return 0 if ok
* otherwise return !=0
*/
int myLu(LinearSystem_Problem* problem, double *z ,  Solver_Options* options)
{
  /* Output info. : 0: ok -  >0: problem (depends on solver) */
  int info = -1;
  int n = problem->size;
  int n2 = n * n;
  double * Maux = (double*)malloc(n2 * sizeof(double));
  int * ipiv = (int *) malloc(n * sizeof(int));
  int LAinfo;
  /* Checks inputs */
  if (problem == NULL || z == NULL)
    numericsError("Equality_Problem", "null input for Equality_Problem and/or unknowns (z)");
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

  free(Maux);
  return info;
}

int LinearSystem_driver(LinearSystem_Problem* problem, double *z , double *w, Solver_Options* options)
{
  int i;
  if (problem->M->storageType == 1)
    numericsError("lcp_driver_DenseMatrix", "forbidden type of storage for the matrix M of the LCP");
  for (i = 0; i < problem->size; i++)
    w[i] = 0.0;
  return myLu(problem, z , options);
}

void displayLS(LinearSystem_Problem* p)
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
