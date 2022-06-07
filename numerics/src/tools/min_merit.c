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

#include "min_merit.h"
#include "SiconosBlas.h"     // for cblas_dcopy
#include <stdio.h>           // for NULL
#include "NumericsMatrix.h"  // for NumericsMatrix
#include "assert.h"          // for assert

void F_min(int n1, int n2, double* restrict z, double* restrict F, double* restrict Fmin)
{
  assert(z != NULL);
  assert(F != NULL);
  assert(Fmin != NULL);

  for(int i = 0; i < n1; ++i)
  {
    Fmin[i] = F[i];
  }
  for(int i = n1 ; i < n1+n2 ; ++i)
  {
    Fmin[i] = z[i] <= F[i] ? z[i] : F[i];
  }
}

void Jac_F_min(int n1, int n2, double* restrict z, double* restrict F, NumericsMatrix* restrict nabla_F, NumericsMatrix* restrict H)
{

  double* nabla_F_dense = nabla_F->matrix0;
  double* H_dense = H->matrix0;

  int n = n1 + n2;
  if(n1 > 0)
  {
    //printf("Jac_F_min: mixed case has to validated -- xhub\n");
    cblas_dcopy(n*n, nabla_F_dense, 1, H_dense, 1);
  }
  // Facchinei--Pang p. 660 and 661
  // i \in alpha if F_i > z_i
  // i \in beta if F_i == z_i
  // i \in gamma if F_i < z_i
  // We made the choice here to have the trivial case z_i + d_i = 0 for alpha U
  // beta. See Facchinei--Pang for more details
  //
  // TODO implement sparse version where the trivial case can be fully exploited
  for(int i = n1; i < n; ++i)
  {
    if(z[i] <= F[i])  // i in beta U alpha
    {
      for(int j = 0; j < n; ++j)
      {
        H_dense[j * n + i] = 0.0;
      }
      H_dense[i * n + i] = 1.0;

    }
    else // i in gamma
    {
      for(int j = 0; j < n; ++j)
      {
        H_dense[j * n + i] = nabla_F_dense[j * n + i];
      }
    }
  }
}
