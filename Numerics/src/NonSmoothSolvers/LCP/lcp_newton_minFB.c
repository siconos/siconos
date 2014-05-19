/* Siconos-Numerics, Copyright INRIA 2005-2014.
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
#include <math.h>
#include <float.h>

#include "LCP_Solvers.h"
#include "SiconosLapack.h"
#include "Newton_Methods.h"
#include "FischerBurmeister.h"
#include "min_merit.h"
#include "lcp_newton_FB.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

void lcp_min(void* data_opaque, double* z, double* F, double* Fmin)
{
  F_min(0, ((LinearComplementarityProblem *)data_opaque)->size, z, F, Fmin);
}

void min_compute_H_lcp(void* data_opaque, double* z, double* F, double* workV1, double* workV2, double* H)
{
  LinearComplementarityProblem* data = (LinearComplementarityProblem *)data_opaque;
  unsigned int n = data->size;
  assert(data->M);
  assert(data->M->matrix0);
  double* M = data->M->matrix0;

  Jac_F_min(0, n, z, F, M, H);
}

void lcp_newton_minFB(LinearComplementarityProblem* problem, double *z, double *w, int *info , SolverOptions* options)
{
  functions_FBLSA functions_minFBLSA_lcp;
  functions_minFBLSA_lcp.compute_F = &FB_compute_F_lcp;
  functions_minFBLSA_lcp.compute_F_merit = &lcp_FB;
  functions_minFBLSA_lcp.compute_H = &FB_compute_H_lcp;
  functions_minFBLSA_lcp.compute_error = &FB_compute_error_lcp;
  functions_minFBLSA_lcp.compute_F_desc = &lcp_min;
  functions_minFBLSA_lcp.compute_H_desc = &min_compute_H_lcp;

  newton_FBLSA(problem->size, z, w, info, (void *)problem, options, &functions_minFBLSA_lcp);
}
