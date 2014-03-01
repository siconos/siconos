/* Siconos-Numerics, Copyright INRIA 2005-2011.
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
#include <float.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "SiconosLapack.h"
#include "pinv.h"
/**
 * n : row number
 * m : col number
 * Warnning n correspond to M in the LAPACK routine, and m to N.

 This routine computes the pseudo inverse of A and returns its conditionning.
 
 */
double pinv(double * A, int n, int m, double tolerance)
{
  int dimS = min(n,m);
  double * S =  (double*)malloc(dimS * sizeof(*S));
  int LDU = n;
  double *U = (double*)malloc(LDU * n * sizeof(double));
  int LDVT = m;
  double *VT = (double*)malloc(LDVT * m * sizeof(double));
  int InfoDGSVD = -1;
  double * superb = (double*)malloc((min(m, n) - 1)*sizeof(double));
  char JOBU = 'A', JOBVT = 'A'; 
  DGESVD(JOBU, JOBVT, n, m, A, n, S, U, LDU, VT, LDVT, superb, &InfoDGSVD);

  double conditioning =  S[0] / S[dimS - 1];
  int rank = 0;
  for (int i = 0; i < dimS ; i++)
  {
    if (S[i] > tolerance)
    {
      rank ++;
      S[i] = 1.0 / S[i];
    }
  }

  /*Compute the pseudo inverse */
  /* Costly version with full DGEMM*/
  double * Utranstmp = (double*)malloc(n * m * sizeof(double));
  for (int i = 0;  i < dimS; i++)
  {
    for (int j = 0;  j < n; j++)
    {
      Utranstmp[i + j * m] = S[i] * U[j + i * n];
    }
  }
  for (int i = dimS;  i < m; i++)
  {
    for (int j = 0;  j < n; j++)
    {
      Utranstmp[i + j * m] = 0.0;
    }
  }

  cblas_dgemm(CblasColMajor,CblasTrans, CblasNoTrans, m, n, m, 1.0, VT, m, Utranstmp, m, 0.0, A, m);
  /*     for (int i = 0;  i < n; i++){ */
  /*  for (int j = 0;  j < n; j++) */
  /*      { */
  /*   U[j+i*n] = S[i]*U[j+i*n]; */
  /*      } */
  /*     } */

  /*   for (int i = 0;  i < rank; i++){ */
  /*  for (int j = 0;  j < m; j++) */
  /*      { */

  /*   A[i+j*n] =0.0; */
  /*   for (int k = 0;  k < rank; k++){ */
  /*       A[i+j*n] += VT[k+i*n]*U[j+k*n]; */
  /*   } */
  /*      } */
  /*     } */

  free(U);
  free(VT);
  free(Utranstmp);
  free(S);
  free(superb);
  return conditioning;
}
