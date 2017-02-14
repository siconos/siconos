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
#include <float.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "SiconosLapack.h"

#include "cond.h"

double cond(double * A, int n, int m)
{
//#ifdef COMPLETE_LAPACK_LIBRARIES
  int dimS = m < n ? m : n;
  double * S =  (double *)malloc(dimS * sizeof(double));

  char JOBU = 'N';
  int LDU = 1;
  double *U = NULL;
  char JOBVT = 'N';
  int LDVT = 1;
  double *VT = NULL;
  size_t size = n * m * sizeof(double);
  double *Atmp = (double *)malloc(size);
  memcpy(Atmp, A, size);

  lapack_int InfoDGSVD = -1;

  double * superb =  (double *)malloc((min(m, n) - 1)* sizeof(double));
  DGESVD(JOBU, JOBVT, n, m, A, n, S, U, LDU, VT, LDVT, superb, &InfoDGSVD);

/* #else */
/*   int LWORK = -1; */
/*   double * WORK; */
/*   WORK = malloc(sizeof(*WORK)); */
/*   assert(WORK); */
/*   DGESVD(&JOBU, &JOBVT, n, m, A, m, S, U, LDU, VT, LDVT, WORK, LWORK, InfoDGSVD); */
/*   LWORK = (int)(WORK[0]); */
/*   WORK = realloc(WORK, LWORK * sizeof * WORK); */
/*   DGESVD(&JOBU, &JOBVT, n, m, A, m, S, U, LDU, VT, LDVT, WORK, LWORK, InfoDGSVD); */
/*   free(WORK); */
/* #endif */

  printf("SVD of A :\n ");
  printf("[\t ");
  for (int i = 0; i < dimS ; i++)
  {
    printf("%14.7e\t", S[i]);
  }
  printf("]\n ");
  memcpy(A, Atmp, size);

  double conditioning = S[0] / S[dimS - 1];

  free(superb);
  free(Atmp);
  free(S);

  return conditioning;
/* #else */
/*   fprintf(stderr, "Numerics. cond.c dgesvd not found\n"); */
/*   return 0.0; */
/* #endif */

}
