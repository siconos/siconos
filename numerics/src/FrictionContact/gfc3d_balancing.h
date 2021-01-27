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
#ifndef GLOBALFRICTIONCONTACT3DBALANCING_H
#define GLOBALFRICTIONCONTACT3DBALANCING_H


struct GlobalFrictionContactProblem_balancing_data
{
  BalancingMatrices * B_for_M;
  BalancingMatrices * B_for_H;
  GlobalFrictionContactProblem * original_problem;
  double alpha;
  double beta;
  double gamma;
};


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  GlobalFrictionContactProblem*  gfc3d_balancing_problem(GlobalFrictionContactProblem* problem,
                                                       SolverOptions* options);

  void gfc3d_balancing_go_to_balanced_variables(GlobalFrictionContactProblem* balanced_problem,
                                  SolverOptions* options,
                                  double *r, double *u, double* v);
  
  void gfc3d_balancing_back_to_original_variables(GlobalFrictionContactProblem* balanced_problem,
                               SolverOptions* options,
                               double *r, double *u, double *v);
  
  GlobalFrictionContactProblem*  gfc3d_balancing_free(GlobalFrictionContactProblem* balanced_problem,
                            SolverOptions* options);

  GlobalFrictionContactProblem_balancing_data  * gfc3d_balancing_data_free(GlobalFrictionContactProblem_balancing_data * data);
  
  GlobalFrictionContactProblem_balancing_data  * gfc3d_balancing_data_new();
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
