// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "SecondOrderConeLinearComplementarityProblem_as_VI.h"

#include "NumericsMatrix.h"                               // for NM_gemv
#include "SecondOrderConeLinearComplementarityProblem.h"  // for SecondOrder...
#include "VariationalInequality.h"                        // for Variational...
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "SiconosBlas.h"  // for cblas_dcopy
#include "numerics_errors.h"
#include "projectionOnCone.h"  // for projectionO...
#include "siconos_debug.h"     // for DEBUG_PRINT

void Function_VI_SOCLCP(void* self, int n_notused, double* x, double* F) {
  DEBUG_PRINT("Function_VI_FC3D(void * self, double *x, double *F)\n")
  VariationalInequality* vi = (VariationalInequality*)self;
  SecondOrderConeLinearComplementarityProblem_as_VI* pb =
      (SecondOrderConeLinearComplementarityProblem_as_VI*)vi->env;
  SecondOrderConeLinearComplementarityProblem* soclcp = pb->soclcp;
  // frictionContactProblem_display(fc3d);

  int n = soclcp->n;

  cblas_dcopy(n, soclcp->q, 1, F, 1);
  NM_gemv(1.0, soclcp->M, x, 1.0, F);
}

void Projection_VI_SOCLCP(void* viIn, double* x, double* PX) {
  DEBUG_PRINT("Projection_VI_SOCLCP(void *viIn, double *x, double *PX)\n")

  VariationalInequality* vi = (VariationalInequality*)viIn;
  SecondOrderConeLinearComplementarityProblem_as_VI* pb =
      (SecondOrderConeLinearComplementarityProblem_as_VI*)vi->env;
  SecondOrderConeLinearComplementarityProblem* soclcp = pb->soclcp;
  // SecondOrderConeLinearComplementarityProblem_display(soclcp);

  int cone = 0;
  int n = soclcp->n;
  cblas_dcopy(n, x, 1, PX, 1);
  int dim;
  for (cone = 0; cone < soclcp->nc; ++cone) {
    dim = soclcp->coneIndex[cone + 1] - soclcp->coneIndex[cone];
    projectionOnSecondOrderCone(&PX[soclcp->coneIndex[cone]], soclcp->tau[cone], dim);
  }
}
