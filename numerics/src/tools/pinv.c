/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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
#include "pinv.h"
#include <stdlib.h>         // for free, malloc
#include "NSSTools.h"  // for min
#include "SiconosBlas.h"    // for cblas_dgemm, CblasColMajor, CblasNoTrans
#include "SiconosLapack.h"  // for DGESVD, lapack_int

/**
 * n : row number
 * m : col number
 * Warning n correspond to M in the LAPACK routine, and m to N.

 This routine computes the pseudo inverse of A and returns its conditionning.

 */
double pinv(double * A, int n, int m, double tolerance)
{
  int dimS = min(n,m);
  double * S = (double*)malloc(dimS * sizeof(double));
  int LDU = n;
  double *U = (double*)malloc(LDU * n * sizeof(double));
  int LDVT = m;
  double *VT = (double*)malloc(LDVT * m * sizeof(double));
  lapack_int InfoDGSVD = -1;
  double * superb = (double*)malloc((min(m, n) - 1)*sizeof(double));
  char JOBU = 'A', JOBVT = 'A';
  DGESVD(JOBU, JOBVT, n, m, A, n, S, U, LDU, VT, LDVT, superb, &InfoDGSVD);

  double conditioning =  S[0] / S[dimS - 1];
  int rank = 0;
  for(int i = 0; i < dimS ; i++)
  {
    if(S[i] > tolerance)
    {
      rank ++;
      S[i] = 1.0 / S[i];
    }
  }

  /*Compute the pseudo inverse */
  /* Costly version with full DGEMM*/
  double * Utranstmp = (double*)malloc(n * m * sizeof(double));
  for(int i = 0;  i < dimS; i++)
  {
    for(int j = 0;  j < n; j++)
    {
      Utranstmp[i + j * m] = S[i] * U[j + i * n];
    }
  }
  for(int i = dimS;  i < m; i++)
  {
    for(int j = 0;  j < n; j++)
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
