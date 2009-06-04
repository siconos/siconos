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
#include <float.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "LA.h"

double cond(double * A, int n, int m)
{
  int LWORK = -1;
  double * WORK;
  WORK = malloc(sizeof(*WORK));
  int dimS = n;
  if (m < n) dimS = m;
  double * S =  malloc(dimS * sizeof(*S));

  char JOBU[1] = "N";
  int LDU = 1;
  double *U = NULL;
  char JOBVT[1] = "N";
  int LDVT = 1;
  double *VT = NULL;
  int size = n * m * sizeof(double);
  double *Atmp = (double *)malloc(size);
  memcpy(Atmp, A, size);

  int InfoDGSVD = -1;

  assert(WORK);
  DGESVD(JOBU, JOBVT, n, m, A, m, S, U, LDU, VT, LDVT, WORK, LWORK, InfoDGSVD);
  LWORK = (int)(WORK[0]);
  WORK = realloc(WORK, LWORK * sizeof * WORK);
  DGESVD(JOBU, JOBVT, n, m, A, m, S, U, LDU, VT, LDVT, WORK, LWORK, InfoDGSVD);

  printf("SVD of A :\n ");
  printf("[\t ");
  for (int i = 0; i < dimS ; i++)
  {
    printf("%14.7e\t", S[i]);
  }
  printf("]\n ");
  memcpy(A, Atmp, size);
  double conditioning =  S[0] / S[dimS - 1];

  free(Atmp);
  free(WORK);
  free(S);
  return conditioning;

}
