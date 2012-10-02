/* Siconos-Numerics, Copyright INRIA 2005-2012.
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "FrictionContact2D_Solvers.h"
#include "LCP_Solvers.h"
#include <assert.h>

int FrictionContact2D_tolcp(FrictionContactProblem* problem, LinearComplementarityProblem * lcp_problem)
{
  if (problem->dimension != 2)
  {
    numericsError("FrictionContact2D_tolcp", "Dimension of the problem : problem-> dimension is not compatible or is not set");
    return 1;
  }
  int nc = problem->numberOfContacts;
  lcp_problem->size = 3 * nc ;
  lcp_problem->M = (NumericsMatrix *)malloc(sizeof(NumericsMatrix));
  lcp_problem->M->size0 = 3 * nc ;
  lcp_problem->M->size1 = 3 * nc ;

  lcp_problem->M->storageType = 0;
  lcp_problem->M->matrix1 = NULL;
  lcp_problem->M->matrix0 = (double*)malloc(lcp_problem->size * lcp_problem->size * sizeof(double));;
  lcp_problem->q = (double*)malloc(lcp_problem->size * sizeof(double));
  int n = 2 * nc;
  int i, j;
  for (i = 0; i < nc; i++)
  {
    for (j = 0; j < nc; j++)
    {
      /* first Column */
      lcp_problem->M->matrix0[3 * i + 3 * j * lcp_problem->size] =
        problem->M->matrix0[2 * i + 2 * j * n] - problem->mu[i] * problem->M->matrix0[2 * i + (2 * j + 1) * n] ; // compute W_NN-mu W_NT

      lcp_problem->M->matrix0[3 * i + 1 + 3 * j * lcp_problem->size] =
        problem->M->matrix0[(2 * i + 1) + 2 * j * n] - problem->mu[i] * problem->M->matrix0[(2 * i + 1) + (2 * j + 1) * n] ; // compute W_TN-mu W_TT

      if (i == j)
      {
        lcp_problem->M->matrix0[3 * i + 2 + 3 * j * lcp_problem->size] = 2.0 * problem->mu[i];
      }
      else
      {
        lcp_problem->M->matrix0[3 * i + 2 + 3 * j * lcp_problem->size] = 0.0;
      }

      /* second Column */
      lcp_problem->M->matrix0[3 * i + (3 * j + 1)*lcp_problem->size] =
        problem->M->matrix0[2 * i + (2 * j + 1) * n] ; // set W_NT

      lcp_problem->M->matrix0[3 * i + 1 + (3 * j + 1)*lcp_problem->size] =
        problem->M->matrix0[(2 * i + 1) + (2 * j + 1) * n] ; // set WTT

      if (i == j)
      {
        lcp_problem->M->matrix0[3 * i + 2 + (3 * j + 1)*lcp_problem->size] =  -1.0;
      }
      else
      {
        lcp_problem->M->matrix0[3 * i + 2 + (3 * j + 1)*lcp_problem->size] =  0.0;
      }

      /* Third Column */
      lcp_problem->M->matrix0[3 * i + (3 * j + 2)*lcp_problem->size] =  0.0;


      if (i == j)
      {
        lcp_problem->M->matrix0[3 * i + 1 + (3 * j + 2)*lcp_problem->size] =  1.0;
      }
      else
      {
        lcp_problem->M->matrix0[3 * i + 1 + (3 * j + 2)*lcp_problem->size] =  0.0;
      }

      lcp_problem->M->matrix0[3 * i + 2 + (3 * j + 2)*lcp_problem->size] =  0.0;

    }
  }


  for (i = 0; i < nc; i++)
  {
    lcp_problem->q[3 * i  ] = problem->q[2 * i];
    lcp_problem->q[3 * i + 1] = problem->q[2 * i + 1];
    lcp_problem->q[3 * i + 2] = 0.0;;
  }

  return 0;

}
