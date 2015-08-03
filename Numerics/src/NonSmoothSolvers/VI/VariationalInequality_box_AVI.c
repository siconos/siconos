/* Siconos-Numerics, Copyright INRIA 2005-2014
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
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

#include "VariationalInequality_Solvers.h"
#include "VariationalInequality_computeError.h"
#include "Relay_Solvers.h"
#include "SiconosSets.h"
#include "SiconosBlas.h"
#include "Newton_Methods.h"
#include "VI_Newton.h"

void vi_compute_decent_dir_by_avi(void* problem, double* z, double* F, double* descent_dir, SolverOptions* options);
void vi_compute_decent_dir_by_avi(void* problem, double* z, double* F, double* descent_dir, SolverOptions* options)
{
  VariationalInequality* vi_pb = (VariationalInequality*) problem;
  int n = vi_pb->size;
  vi_pb->F(vi_pb->env, n, z, F);
  RelayProblem* relay_pb = (RelayProblem*) options->solverData;
  vi_pb->compute_nabla_F(vi_pb->env, n, z, relay_pb->M->matrix0);

  cblas_dcopy(n, F, 1, relay_pb->q, 1);
  prodNumericsMatrix(n, n, -1.0, relay_pb->M, z, 1.0, relay_pb->q);

  int local_info = 0;
  relay_avi_caoferris(relay_pb, descent_dir, F, &local_info, options->internalSolvers);

  cblas_daxpy(n, -1.0, z, 1, descent_dir, 1);

//  avi_caoferris_stage3(avi_options, x, s_vec, problem->size, A,
//  options->internalSolver);
//  for (unsigned int i = 0; i<n; ++i) x[i] = s_vec[A[i]-1] + problem->lb[i];
}

void * vi_get_set(void* problem);
void * vi_get_set(void* problem)
{
  return ((VariationalInequality*) problem)->set;
}

void vi_box_AVI_LSA(VariationalInequality* problem, double* z, double* F, int* info, SolverOptions* options)
{

  int n = problem->size;

  RelayProblem relay_pb;
  relay_pb.size = n;
  NumericsMatrix num_mat;
  relay_pb.M = &num_mat;
  fillNumericsMatrix(&num_mat, NM_DENSE, n, n, malloc(n * n * sizeof(double)));
  relay_pb.q = (double*) malloc(n * sizeof(double));

  box_constraints* box = (box_constraints*) problem->set;
  relay_pb.lb = box->lb;
  relay_pb.ub = box->ub;
  options->solverData = &relay_pb;

  functions_LSA functions_AVI_LSA;
  init_lsa_functions(&functions_AVI_LSA, &VI_compute_F, &VI_compute_F_box_Qi);
  functions_AVI_LSA.compute_H = &VI_compute_H_box_Qi;
  functions_AVI_LSA.compute_error = &VI_compute_error_box;
  functions_AVI_LSA.descent_direction = &vi_compute_decent_dir_by_avi;
  functions_AVI_LSA.get_set_from_problem_data = &vi_get_set;
  options->iparam[SICONOS_IPARAM_LSA_FORCE_ARCSEARCH] = 1;

  newton_LSA(problem->size, z, F, info, (void *)problem, options, &functions_AVI_LSA);


  freeNumericsMatrix(&num_mat);
  free(relay_pb.q);
}

void vi_box_AVI_extra_SolverOptions(SolverOptions* options)
{
  options->numberOfInternalSolvers = 1;
  options->internalSolvers = (SolverOptions *)malloc(sizeof(SolverOptions));
  relay_avi_caoferris_setDefaultSolverOptions(options->internalSolvers);
}
