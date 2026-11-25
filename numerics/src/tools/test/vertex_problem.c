// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <stdio.h>

#include "NumericsMatrix.h"
#include "SiconosSets.h"
#include "vertex_extraction.h"

int main(void) {
  double Hdat[] = {1, 0, -1, 0, 0, 1, 0, -1};

  NumericsMatrix* H = NM_create_with_data(NM_DENSE, 4, 2, Hdat);

  double K[] = {-2, -3, -7, -8};

  polyhedron P = {SICONOS_SET_POLYHEDRON, 4, 0, H, K, NULL, NULL};
  int basis[11] = {0};

  siconos_find_vertex(&P, 2, basis);
  H->dense = NULL;
  NM_free(H);
  for (unsigned i = 0; i < 11; ++i) printf("%d ", basis[i]);
  printf("\n");
  return 0;
}
