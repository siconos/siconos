// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "fc3d_solvers_wr.h"

#include <stdlib.h>  // for malloc, calloc

#include "CSparseMatrix.h"                 // for CSparseMatrix
#include "FrictionContactProblem.h"        // for FrictionContactProblem
#include "GlobalFrictionContactProblem.h"  // for GlobalFrictionContactProblem
#include "NumericsMatrix.h"                // for NM_create, NM_display, NM_eye
#include "NumericsSparseMatrix.h"          // for NumericsSparseMatrix, NSM_...
#include "SiconosBlas.h"                   // for cblas_dcopy
#include "numerics_errors.h"

// #define TEST_COND

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

/* Global Variable for the reformulation of the problem */

GlobalFrictionContactProblem* fc3d_expand_to_global_problem(FrictionContactProblem* problem) {
  int dimension = problem->dimension;
  int nc = problem->numberOfContacts;
  int m = dimension * nc;

  GlobalFrictionContactProblem* globalproblem = globalFrictionContactProblem_new();

  globalproblem->numberOfContacts = problem->numberOfContacts;
  globalproblem->dimension = problem->dimension;

  globalproblem->mu = (double*)malloc(m * sizeof(double));
  cblas_dcopy(nc, problem->mu, 1, globalproblem->mu, 1);
  globalproblem->q = (double*)malloc(m * sizeof(double));
  cblas_dcopy(m, problem->q, 1, globalproblem->q, 1);
  globalproblem->b = (double*)calloc(m, sizeof(double));

  NumericsMatrix* M = problem->M;

  NumericsMatrix* Mglobal = NM_create(NM_SPARSE, m, m);
  CSparseMatrix* M_triplet = NM_triplet(M);
  Mglobal->sparse->triplet = M_triplet;
  Mglobal->sparse->origin = NSM_TRIPLET;
  globalproblem->M = Mglobal;

  NumericsMatrix* Hglobal = NM_eye(m);
  NM_display(Hglobal);
  globalproblem->H = M = Hglobal;

  return globalproblem;
}

/* Backward compatibility function name (deprecated) */
GlobalFrictionContactProblem* fc3d_reformulation_global_problem(
    FrictionContactProblem* problem) {
  return fc3d_expand_to_global_problem(problem);
}
