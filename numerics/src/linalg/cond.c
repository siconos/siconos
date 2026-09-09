/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#include "cond.h"

#include <stdio.h>   // for printf, NULL, size_t
#include <stdlib.h>  // for free, malloc
#include <string.h>  // for memcpy

#include "NSSTools.h"       // for min
#include "SiconosLapack.h"  // for DGESVD, lapack_int
#include "safe_casts.h"

double cond(double *A, size_t n, size_t m) {
  size_t dimS = m < n ? m : n;
  double *S = (double *)malloc(dimS * sizeof(double));

  char JOBU = 'N';
  int LDU = 1;
  char JOBVT = 'N';
  int LDVT = 1;
  size_t size = n * m * sizeof(double);
  double *Atmp = (double *)malloc(size);
  memcpy(Atmp, A, size);

  lapack_int InfoDGSVD = -1;

  double *superb = (double *)malloc((min(m, n) - 1) * sizeof(double));
  lapack_int n_la = to_blasint(n);
  lapack_int m_la = to_blasint(m);
  DGESVD(JOBU, JOBVT, n_la, m_la, Atmp, n_la, S, NULL, LDU, NULL, LDVT, superb, &InfoDGSVD);
  if (InfoDGSVD != 0) {
    fprintf(stderr, "cond: DGESVD failed (info=%d)\n", InfoDGSVD);
    free(S);
    free(Atmp);
    free(superb);
    return -1.0;
  }
  // printf("SVD of A :\n ");
  // printf("[\t ");
  // for (size_t i = 0; i < dimS; i++) {
  //   printf("%14.7e\t", S[i]);
  // }
  // printf("]\n ");
  // memcpy(A, Atmp, size);

  double conditioning = (S[dimS - 1] > 0.0) ? S[0] / S[dimS - 1] : -1.0;

  free(superb);
  free(Atmp);
  free(S);

  return conditioning;
}
