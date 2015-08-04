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

#include "Newton_Methods.h"

#include "NCP_cst.h"
#include "NonlinearComplementarityProblem.h"
#include "FischerBurmeister.h"
#include "SiconosBlas.h"

#include "NCP_Solvers.h"

#include "ncp_newton_FBLSA.h"

void ncp_FB(void* data_opaque, double* z, double* F, double* F_FB)
{
  phi_FB(((NonlinearComplementarityProblem *)data_opaque)->n, z, F, F_FB);
}

void FB_compute_F_ncp(void* data_opaque, double* z, double* F)
{
  // Computation of the new value F(z)
  NonlinearComplementarityProblem* data = (NonlinearComplementarityProblem *)data_opaque;
  data->compute_F(data->env, data->n, z, F);
}

void FB_compute_H_ncp(void* data_opaque, double* z, double* F, double* workV1, double* workV2, double* H)
{
  NonlinearComplementarityProblem* data = (NonlinearComplementarityProblem *)data_opaque;

  data->compute_nabla_F(data->env, data->n, z, data->nabla_F->matrix0);

  Jac_F_FB(0, data->n, z, F, workV1, workV2, data->nabla_F->matrix0, H);
}

void FB_compute_error_ncp(void* data_opaque, double* z, double* w, double* Jac_F_merit, double tol, double* err)
{
  NonlinearComplementarityProblem* data = (NonlinearComplementarityProblem *)data_opaque;
  *err = cblas_dnrm2(data->n, Jac_F_merit, 1);
}

void ncp_newton_FBLSA(NonlinearComplementarityProblem* problem, double *z, double* F, int *info, SolverOptions* options)
{
  functions_LSA functions_FBLSA_ncp;
  init_lsa_functions(&functions_FBLSA_ncp, &FB_compute_F_ncp, &ncp_FB);
  functions_FBLSA_ncp.compute_H = &FB_compute_H_ncp;
  functions_FBLSA_ncp.compute_error = &FB_compute_error_ncp;

  newton_LSA(problem->n, z, F, info, (void *)problem, options, &functions_FBLSA_ncp);
}

/*
void ncp_newton_FBLSA_setDefaultSolverOptions(SolverOptions* options)
{
  fill_SolverOptions(options, SICONOS_NCP_NEWTON_FBLSA, 5, 5, 100, 1e-16);
}
*/
