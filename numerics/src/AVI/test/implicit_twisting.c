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

/* Simple twisting example with the double integrator */

/* to have drand48 ...
 * hmhmhm ??? -- xhub */
#define _XOPEN_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "AVI_Solvers.h"
#include "AVI_cst.h"
#include "SiconosSets.h"
#include "NumericsMatrix.h"
#define TS 10e-3
#define NB_ITER 10000
#define TOL_NC 1e-12

int main(void)
{
  double x[2];
  unsigned short xsubi1[] = {0, 0, 0};
  unsigned short xsubi2[] = {0, 0, 0};
  x[0] = 50*erand48(xsubi1);
  x[1] = 50*erand48(xsubi2);
  /* column major */
  double Hdat[8] = {1.0, -TS/2.0, -1.0, TS/2.0, 0.0, 1.0, 0.0, -1.0};
  double K[4] = {-1.0, -1.0, -1.0, -1.0};

  NumericsMatrix H;
  NM_null(&H);
  NM_fill(&H, NM_DENSE, 4, 2, Hdat);

  double v1[] = {-1.0, -1.0 -TS/2.0};
  double v2[] = {-1.0, 1.0 -TS/2.0};
  double v3[] = {1.0, 1.0 + TS/2.0};
  double v4[] = {1.0, -1.0 + TS/2.0};

  polyhedron poly = { SICONOS_SET_POLYHEDRON, 4, 0, &H, K, NULL, NULL };

  /* twisting gain */
  double G = 10;
  double beta = .1;
  NumericsMatrix num_mat;
  double M[4] = { G*TS*TS/2.0, G*TS, beta*G*TS*TS/2.0, beta*G*TS };
  NM_null(&num_mat);
  NM_fill(&num_mat, NM_DENSE, 2, 2, M);

  double q[2] = { x[0] + TS*x[1], x[1] };

  AffineVariationalInequalities avi;
  avi.size = 2;
  avi.M = &num_mat;
  avi.q = q;
  avi.d = NULL;
  avi.poly.split = &poly;

  SolverOptions options;
  solver_options_set(&options, SICONOS_AVI_CAOFERRIS);

  _Bool c1, c2, c3, c4;
  unsigned N = 0;
  int info = 0;
  double lambda[2] = {0.};
  double sigma[2] = {0.};
  do
  {
    N++;
    lambda[0] = 0.0;
    lambda[1] = 0.0;
    sigma[0] = 0.0;
    sigma[1] = 0.0;

    info = avi_caoferris(&avi, lambda, sigma, &options);
    if (info) fprintf(stderr, "SOLVER FAILED!\tinfo=%d\n", info);

/* /    printf("x_k: %2.6e %2.6e\n", x[0], x[1]);
    printf("lambda: %2.6e %2.6e\n", lambda[0], lambda[1]);
    printf("u = %2.6e\n", lambda[0] + beta*lambda[1]);*/
    sigma[0] = x[0] + TS*x[1] + M[0]*lambda[0] + M[2]*lambda[1];
    sigma[1] = x[1] + M[1]*lambda[0] + M[3]*lambda[1];
/*    printf("x_{k+1}: %2.6e %2.6e\n", sigma[0], sigma[1]);*/

    /* let's check is it is correct */
    c1 = (sigma[0]*(v1[0] + lambda[0]) + sigma[1]*(v1[1] + lambda[1])) <= TOL_NC;
    c2 = (sigma[0]*(v2[0] + lambda[0]) + sigma[1]*(v2[1] + lambda[1])) <= TOL_NC;
    c3 = (sigma[0]*(v3[0] + lambda[0]) + sigma[1]*(v3[1] + lambda[1])) <= TOL_NC;
    c4 = (sigma[0]*(v4[0] + lambda[0]) + sigma[1]*(v4[1] + lambda[1])) <= TOL_NC;

    if (!c1)
    {
      fprintf(stderr, "ERROR in implicit_twisting.c, bad value of lambda. test c1 failed\n");
      fprintf(stderr, "<v1-lambda, sigma> = %2.10e\n", (sigma[0]*(v1[0] + lambda[0]) + sigma[1]*(v1[1] + lambda[1])));
      goto expose_failure;
    }
    if (!c2)
    {
      fprintf(stderr, "ERROR in implicit_twisting.c, bad value of lambda. test c2 failed\n");
      fprintf(stderr, "<v2-lambda, sigma> = %2.10e\n", (sigma[0]*(v2[0] + lambda[0]) + sigma[1]*(v2[1] + lambda[1])));
      goto expose_failure;
    }

    if (!c3)
    {
      fprintf(stderr, "ERROR in implicit_twisting.c, bad value of lambda. test c3 failed\n");
      fprintf(stderr, "<v3-lambda, sigma> = %2.10e\n", (sigma[0]*(v3[0] + lambda[0]) + sigma[1]*(v3[1] + lambda[1])));
      goto expose_failure;
    }

    if (!c4)
    {
      fprintf(stderr, "ERROR in implicit_twisting.c, bad value of lambda. test c4 failed\n");
      fprintf(stderr, "<v4-lambda, sigma> = %2.10e\n", (sigma[0]*(v4[0] + lambda[0]) + sigma[1]*(v4[1] + lambda[1])));
      goto expose_failure;
    }
    x[0] = sigma[0];
    x[1] = sigma[1];
    q[0] = x[0] + TS*x[1];
    q[1] = x[1];
  }
  while (!info && c1 && c2 && c3 && c4 && (N <= NB_ITER));

  printf("final x: %2.6e %2.6e\n", sigma[0], sigma[1]);

  solver_options_delete(&options);
  return info;

expose_failure:
  fprintf(stderr, "x = %2.10e %2.10e\n", x[0], x[1]);
  fprintf(stderr, "q = %2.10e %2.10e\n", q[0], q[1]);
  fprintf(stderr, "lambda = %2.10e %2.10e\n", lambda[0], lambda[1]);

  return 1;
}
