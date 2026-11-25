// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <stdlib.h>  // for malloc

#include "test_utils.h"  // for data_collection

const char** data_collection() {
  int n_data_1 = 10;

  const char** data_collection_1 = (const char**)malloc(n_data_1 * sizeof(const char*));
  int n_data = 0;
  data_collection_1[n_data++] = "./data/cohesive_test_2x2.dat";
  data_collection_1[n_data++] = "./data/sphere_2x2_mu0.dat";
  data_collection_1[n_data++] = "./data/sphere_2x2.dat";
  /* 55 contacts + 55 cohesive points, hyperstatic (W has a 15-dim null space) */
  data_collection_1[n_data++] = "./data/diamond_sphere_sphere_pack.dat";

  data_collection_1[n_data++] = "---";

  return data_collection_1;
}
