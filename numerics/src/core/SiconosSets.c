// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "SiconosSets.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "NumericsMatrix.h"
#include "numerics_errors.h"

#ifdef __cplusplus
#undef restrict
#include <sys/cdefs.h>  // for __restrict
#define restrict __restrict
#endif

void project_on_set(int n, double* restrict x, void* restrict set) {
  if (!set) return;
  if (!x) return;

  assert(n > 0);

  int set_id = ((generic_set*)set)->id;
  switch (set_id) {
    case SICONOS_SET_POSITIVE_ORTHANT:
      for (int i = 0; i < n; ++i) {
        if (x[i] < 0.0) x[i] = 0.0;
      }
      break;
    case SICONOS_SET_BOX: {
      box_constraints* box = (box_constraints*)set;
      assert(box->lb);
      assert(box->ub);
      for (int i = 0; i < n; ++i)
        x[i] = x[i] >= box->lb[i] ? (x[i] <= box->ub[i] ? x[i] : box->ub[i]) : box->lb[i];
      break;
    }
    default:
      printf("project_on_set :: not implemented for set of type %d\n", set_id);
      exit(EXIT_FAILURE);
  }
}

void free_siconos_set(void* set) {
  if (!set) return;

  int set_id = ((generic_set*)set)->id;
  switch (set_id) {
    case SICONOS_SET_POSITIVE_ORTHANT:
      break;
    case SICONOS_SET_BOX:
      free_box((box_constraints*)set);
      break;
    case SICONOS_SET_POLYHEDRON:
      free_polyhedron((polyhedron*)set);
      break;
    case SICONOS_SET_POLYHEDRON_UNIFIED:
      free_polyhedron_unified((polyhedron_unified*)set);
      break;
    default:
      printf("free_siconos_set :: not implemented for set of type %d\n", set_id);
      exit(EXIT_FAILURE);
  }
}

void free_box(box_constraints* b) {
  if (!b) return;
  if (b->lb) {
    free(b->lb);
    b->lb = NULL;
  }
  if (b->ub) {
    free(b->ub);
    b->ub = NULL;
  }
}

void free_polyhedron(polyhedron* poly) {
  if (!poly) return;
  if (poly->H) {
    NM_clear(poly->H);
    free(poly->H);
    poly->H = NULL;
  }
  if (poly->K) {
    free(poly->K);
    poly->K = NULL;
  }
  if (poly->Heq) {
    NM_clear(poly->Heq);
    free(poly->Heq);
    poly->Heq = NULL;
  }
  if (poly->Keq) {
    free(poly->Keq);
    poly->Keq = NULL;
  }
}

void free_polyhedron_unified(polyhedron_unified* poly) {
  if (!poly) return;
  if (poly->A) {
    NM_clear(poly->A);
    free(poly->A);
    poly->A = NULL;
  }
  if (poly->b) {
    free(poly->b);
    poly->b = NULL;
  }
  if (poly->type) {
    free(poly->type);
    poly->type = NULL;
  }
}
