// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef GLOBALFRICTIONCONTACT3DBALANCING_H
#define GLOBALFRICTIONCONTACT3DBALANCING_H

#include "GlobalFrictionContactProblem.h"
#include "NumericsMatrix.h"  // For balancing matrix

struct GlobalFrictionContactProblem_balancing_data {
  BalancingMatrices* B_for_M;
  BalancingMatrices* B_for_H;
  GlobalFrictionContactProblem* original_problem;
  double alpha;
  double beta;
  double gamma;
};

#if defined(__cplusplus)
extern "C" {
#endif
void gfc3d_rescaling(GlobalFrictionContactProblem* problem, double alpha, double beta,
                     double gamma);

void gfc3d_balancing_M(GlobalFrictionContactProblem* problem, BalancingMatrices* B_for_M);

void gfc3d_balancing_MHHT(GlobalFrictionContactProblem* problem, BalancingMatrices* B_for_M,
                          BalancingMatrices* B_for_H);

GlobalFrictionContactProblem* gfc3d_balancing_problem(GlobalFrictionContactProblem* problem,
                                                      SolverOptions* options);

void gfc3d_balancing_go_to_balanced_variables(GlobalFrictionContactProblem* balanced_problem,
                                              SolverOptions* options, double* r, double* u,
                                              double* v);

void gfc3d_balancing_back_to_original_variables(GlobalFrictionContactProblem* balanced_problem,
                                                SolverOptions* options, double* r, double* u,
                                                double* v);

GlobalFrictionContactProblem* gfc3d_balancing_free(
    GlobalFrictionContactProblem* balanced_problem, SolverOptions* options);

GlobalFrictionContactProblem_balancing_data* gfc3d_balancing_data_free(
    GlobalFrictionContactProblem_balancing_data* data);

GlobalFrictionContactProblem_balancing_data* gfc3d_balancing_data_new(void);

#if defined(__cplusplus)
}
#endif

#endif
