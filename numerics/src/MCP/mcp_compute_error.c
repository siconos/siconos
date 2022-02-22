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
#include <math.h>                         // for fmax, pow, sqrt
#include <stddef.h>                       // for NULL
#include "SiconosBlas.h"                        // for cblas_dnrm2
#include "MCP_Solvers.h"                  // for mcp_compute_error, mcp_old_...
#include "MixedComplementarityProblem.h"  // for MixedComplementarityProblem
#include "NumericsFwd.h"                  // for MixedComplementarityProblem
#include "numerics_verbose.h"             // for numerics_error

int mcp_compute_error(MixedComplementarityProblem* problem, double *z, double *w,  double * error)
{
  /* Checks inputs */
  if(problem == NULL || z == NULL || w == NULL)
    numerics_error("mcp_compute_error", "null input for problem and/or z and/or w");

  int size =  problem->n1 + problem->n2;
  /* Computes w = F(z) */
  problem->compute_Fmcp(problem->env, size, z, w);

  * error =0.0;
  /* compute the error of the inequality constraints */
  for(int i = problem->n1; i < size ; i++)
  {
    *error += pow(z[i] - fmax(0,(z[i] - w[i])),2);
  }
  *error =sqrt(*error);
  /* add the error of the equality constraints */
  * error += cblas_dnrm2(problem->n1, w, 1);

  /* int incx = 1, incy = 1; */
  /* unsigned int n = problem->size; */
  /* cblas_dcopy(n , problem->q , incx , w , incy);  // w <-q */
  /* NM_gemv(1.0, problem->M, z, 1.0, w); */
  /* double norm_q = cblas_dnrm2(n , problem->q , incx); */
  /* lcp_compute_error_only(n, z, w, error); */
  /* if (fabs(norm_q) > DBL_EPSILON) */
  /*   *error /= norm_q; */



  return 0;

}
int mcp_old_compute_error(MixedComplementarityProblem_old* problem, double *z, double *w,  double * error)
{
  /* Checks inputs */
  if(problem == NULL || z == NULL || w == NULL)
    numerics_error("mcp_old_compute_error", "null input for problem and/or z and/or w");

  int size =  problem->sizeEqualities + problem->sizeInequalities;
  /* Computes w = F(z) */
  problem->computeFmcp(size, z, w);

  * error =0.0;
  /* compute the error of the inequality constraints */
  for(int i = problem->sizeInequalities; i < size ; i++)
  {
    *error += pow(z[i] - fmax(0,(z[i] - w[i])),2);
  }
  *error =sqrt(*error);
  /* add the error of the equality constraints */
  * error += cblas_dnrm2(problem->sizeEqualities, w, 1);

  /* int incx = 1, incy = 1; */
  /* unsigned int n = problem->size; */
  /* cblas_dcopy(n , problem->q , incx , w , incy);  // w <-q */
  /* NM_gemv(1.0, problem->M, z, 1.0, w); */
  /* double norm_q = cblas_dnrm2(n , problem->q , incx); */
  /* lcp_compute_error_only(n, z, w, error); */
  /* if (fabs(norm_q) > DBL_EPSILON) */
  /*   *error /= norm_q; */



  return 0;

}
