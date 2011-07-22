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
#include "LA.h"
#include "NumericsOptions.h"

#define NUMERICS_PINV_DEBUG
/**
 * n : row number
 * m : col number
 * Warnning n correspond to M in the LAPACK routine, and m to N.
 */
double pinv(double * A, int n, int m, double tolerance)
{
#ifdef COMPLETE_LAPACK_LIBRARIES

  int dimS = n;
  if (m < n) dimS = m;
  double * S =  malloc(dimS * sizeof(*S));
  char JOBU = 'A';
  int LDU = n;
  double *U = malloc(n * n * sizeof(double));
  char JOBVT = 'A'  ;
  int LDVT = m;
  double *VT = malloc(m * m * sizeof(double));
  /*    printf("Matrix A:\n "); */
  /*     for (int i = 0; i< n; i++){ */
  /*  for (int j = 0; j < m; j++){ */
  /*      printf("%8.6e\t", A[i+j*n]) ; */
  /*  } */
  /*  printf("\n"); */
  /*     }  */

  int InfoDGSVD = -1;
#ifdef USE_MKL
  double superb[min(m, n) - 1];
  DGESVD(JOBU, JOBVT, n, m, A, n, S, U, LDU, VT, LDVT, superb, 0, InfoDGSVD);
#else
  int LWORK = -1;
  double * WORK = malloc(sizeof(*WORK));
  assert(WORK);
  DGESVD(&JOBU, &JOBVT, n, m, A, n, S, U, LDU, VT, LDVT, WORK, LWORK, InfoDGSVD);
  LWORK = (int)(WORK[0]);
  WORK = realloc(WORK, LWORK * sizeof * WORK);
  DGESVD(&JOBU, &JOBVT, n, m, A, n, S, U, LDU, VT, LDVT, WORK, LWORK, InfoDGSVD);
  free(WORK);
#endif

  /*    printf("Matrix U:\n "); */
  /*     for (int i = 0; i< n; i++){ */
  /*  for (int j = 0; j < n; j++){ */
  /*      printf("%8.6e\t", U[i+j*n]) ; */
  /*  } */
  /*  printf("\n"); */
  /*     }  */
#ifdef NUMERICS_PINV_DEBUG
  printf("pinv: SVD :\n ");
  printf("[\t ");
  for (int i = 0; i < dimS ; i++)
  {
    printf("%14.7e\t", S[i]);
  }
  /*     printf("]\n "); */
  /*     printf("Matrix VT:\n "); */
  /*     for (int i = 0; i< m; i++){ */
  /*  for (int j = 0; j < m; j++){ */
  /*      printf("%8.6e\t", VT[i+j*m]) ; */
  /*  } */
  /*  printf("\n"); */
  /*     }  */

#endif





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
#ifdef NUMERICS_PINV_DEBUG
  printf("\n");
  printf("pinv: Rank of A :%i\n ", rank);
#endif
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
      Utranstmp[i + j * m] = 0.0 * U[j + i * n];
    }
  }
  /*     printf("Matrix Utranstmp:\n "); */
  /*     for (int i = 0; i< m; i++){ */
  /*  for (int j = 0; j < n; j++){ */
  /*      printf("%8.6e\t", Utranstmp[i+j*m]); */
  /*  } */
  /*  printf("\n"); */
  /*     }  */
  //DGEMM(LA_TRANS,LA_NOTRANS,m,n,m,1.0,VT,m,Utranstmp,m,0.0,A,m);
  DGEMM(LA_TRANS, LA_NOTRANS, n, m, n, 1.0, VT, n, Utranstmp, n, 0.0, A, n);

#ifdef NUMERICS_PINV_DEBUG
  printf("\nApinv=[ \n");
  for (int i = 0; i < n ; i++)
  {
    printf("[");
    for (int j = 0; j < m; j++)
      printf("%14.7e\t", A[i + j * n]);
    printf("];\n");
  }
  printf("];\n");
#endif
  /*    printf("Matrix Pseudo-Inverse of A:\n "); */
  /*     for (int i = 0; i< m; i++){ */
  /*  for (int j = 0; j < n; j++){ */
  /*      printf("%8.6e\t", A[i+j*m]); */
  /*  } */
  /*  printf("\n"); */
  /*     }  */

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
  return conditioning;
#else
#include <stdio.h>
  numericsError("pinv", "DGESVD not found");
  return 0.0;
#endif
}
