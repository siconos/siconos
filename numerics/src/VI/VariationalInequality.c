// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "VariationalInequality.h"

#include <assert.h>  // for assert
#include <stdlib.h>  // for free, malloc, exit, EXIT_FAILURE

#include "NumericsMatrix.h"  // for NM_clear
#include "numerics_errors.h"

void variationalInequality_display(VariationalInequality* vi) { assert(vi); }

int variationalInequality_printInFile(VariationalInequality* vi, FILE* file) {
  if (!vi) {
    CHECK_ARG(0, "Numerics, VariationalInequality printInFile failed, NULL input.\n");
  }

  return 0;
}

int variationalInequality_newFromFile(VariationalInequality* vi, FILE* file) { return 0; }

VariationalInequality* freeVariationalInequalityProblem(VariationalInequality* vi) {
  if (!vi) return NULL;
  if (vi->nabla_F) {
    vi->nabla_F = NM_free(vi->nabla_F);
  }
  vi->env = NULL;
  vi->F = NULL;
  vi->compute_nabla_F = NULL;
  vi->ProjectionOnX = NULL;
  vi->set = NULL;
  free(vi);
  return NULL;
}

void variationalInequality_clear(VariationalInequality* vi) {
  vi->size = 0;
  vi->env = NULL;
  vi->F = NULL;
  vi->compute_nabla_F = NULL;
  vi->ProjectionOnX = NULL;
  vi->normVI = 0.;
  vi->istheNormVIset = 0.;
  vi->set = NULL;
  vi->nabla_F = NULL;
}

VariationalInequality* variationalInequality_new(int size) {
  VariationalInequality* fvi = (VariationalInequality*)malloc(sizeof(VariationalInequality));
  variationalInequality_clear(fvi);
  fvi->size = size;

  return fvi;
}

VariationalInequality* newVI(void) {
  VariationalInequality* vi = (VariationalInequality*)malloc(sizeof(VariationalInequality));
  variationalInequality_clear(vi);

  return vi;
}
