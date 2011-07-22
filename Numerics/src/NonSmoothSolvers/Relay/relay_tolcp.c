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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "LA.h"
#include "Relay_Solvers.h"
#include "LCP_Solvers.h"
#include <assert.h>

void relay_tolcp(RelayProblem* problem, LinearComplementarityProblem * lcp_problem)
{
  lcp_problem->size = 2 * problem->size ;
  lcp_problem->M = (NumericsMatrix *)malloc(sizeof(NumericsMatrix));
  lcp_problem->M->size0 = 2 * problem->size ;
  lcp_problem->M->size1 = 2 * problem->size ;

  lcp_problem->M->storageType = 0;
  lcp_problem->M->matrix1 = NULL;
  lcp_problem->M->matrix0 = (double*)malloc(lcp_problem->size * lcp_problem->size * sizeof(double));;
  lcp_problem->q = (double*)malloc(lcp_problem->size * sizeof(double));

  int i, j;
  for (i = 0; i < problem->size; i++)
  {
    for (j = 0; j < problem->size; j++)
    {
      lcp_problem->M->matrix0[i + j * lcp_problem->size] =  problem->M->matrix0[i + j * problem->size];
    }
  }
  for (i = 0; i < problem->size; i++)
  {
    for (j = problem->size; j < 2 * problem->size; j++)
    {
      lcp_problem->M->matrix0[i + j * lcp_problem->size] =  0.0;
    }
    lcp_problem->M->matrix0[i + (i + problem->size)*lcp_problem->size] =  1.0;
  }
  for (i = problem->size; i < 2 * problem->size; i++)
  {
    for (j = 0; j < 2 * problem->size; j++)
    {
      lcp_problem->M->matrix0[i + j * lcp_problem->size] =  0.0;
    }
    lcp_problem->M->matrix0[i + (i - problem->size)*lcp_problem->size] =  -1.0;
  }

  for (i = 0; i < problem->size; i++)
  {
    lcp_problem->q[i] = problem->q[i];
    lcp_problem->q[i + problem->size] = problem->ub[i] - problem->lb[i];
    for (j = 0; j < problem->size; j++)
    {
      lcp_problem->q[i] -= problem->M->matrix0[i + j * (problem->size)] * problem->ub[i];
    }
  }



}
