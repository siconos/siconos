/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

#include <stdio.h>                              // for printf
#include <stdlib.h>                             // for malloc, exit, free
#include <string.h>
#include "MixedLinearComplementarityProblem.h"  // for MixedLinearComplement...
#include "LinearComplementarityProblem.h"       // for LinearComplement...
#include "NumericsMatrix.h"                     // for NM_...
#include "numerics_verbose.h"                   // for verbose */
#include "mlcp_to_lcp.h"                        // for mlcp_to_lcp

/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"

LinearComplementarityProblem*  mlcp_to_lcp(MixedLinearComplementarityProblem* problem)
{

  if(!problem->isStorageType2)
  {
    //mixedLinearComplementarity_display(problem);
    numerics_printf_verbose(0,"mlcp_pgs: Wrong Storage (!isStorageType2) for PGS solver\n");
    MixedLinearComplementarityProblem* mlcp_abcd =  mixedLinearComplementarity_fromMtoABCD(problem);
    //mixedLinearComplementarity_display(mlcp_abcd);
    problem = mlcp_abcd;
    //exit(EXIT_FAILURE);
  }

  int n =problem->n;
  int m =problem->m;

  /* solve A^{-1} C */


  double * A_copy = (double* )calloc(n*n,sizeof(double));
  memcpy(A_copy, problem->A, n*n*sizeof(double)); /* we preserve the original A */
  NumericsMatrix * A = NM_new();
  NM_fill(A, NM_DENSE, n, n, A_copy);
  double * AinvC = (double* )calloc(n*m,sizeof(double));
  memcpy(AinvC, problem->C, n*m*sizeof(double));
  int info =  NM_gesv_expert_multiple_rhs(A, AinvC, m, NM_KEEP_FACTORS);

  if (info)
  {
    numerics_printf_verbose(0,"");
    free(AinvC);
    return NULL;
  }

  /* solve A^{-1} a */
  double * Ainva = (double* )calloc(n,sizeof(double));
  memcpy(Ainva, problem->a, n*sizeof(double));
  NM_gesv_expert(A, Ainva, NM_KEEP_FACTORS);

  LinearComplementarityProblem* lcp= newLCP();

  lcp->size = m;
  lcp->M  = NM_new();
  double * Mdata = (double* )calloc(m*m,sizeof(double));
  NM_fill(lcp->M, NM_DENSE, m, m, Mdata);

  lcp->q = (double* )calloc(m,sizeof(double));

  /* compute B-DA^{-1}C */
  memcpy(lcp->M->matrix0, problem->B, m*m*sizeof(double));
  DEBUG_EXPR(NM_display(lcp->M));

  NumericsMatrix * D = NM_new();
  NM_fill(D, NM_DENSE, m, n, problem->D);
  DEBUG_EXPR(NM_display(D));

  NumericsMatrix * nm_AinvC = NM_new();
  NM_fill(nm_AinvC, NM_DENSE, n, m, AinvC);
  NM_gemm(-1.0, D, nm_AinvC, 1.0,  lcp->M );

  DEBUG_EXPR(NM_display(lcp->M));

  /* compute b - DA^{-1}a */
  memcpy(lcp->q, problem->b, m*sizeof(double));
  NM_gemv(-1.0, D, Ainva, 1.0, lcp->q);



  NM_clear(A);
  D->matrix0=NULL;
  NM_clear(D);

  NM_clear(nm_AinvC);
  return lcp;
}
