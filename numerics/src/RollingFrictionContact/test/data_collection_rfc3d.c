// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <stdlib.h>  // for malloc

#include "test_utils.h"  // for data_collection

const char** data_collection() {
  int n_data_1 = 150;

  const char** data_collection_1 = (const char**)malloc(n_data_1 * sizeof(const char*));
  int n_data = 0;
  /* data_collection_1[n_data++] = "./data/RFC3D_sphere_1.dat"; */
  /* data_collection_1[n_data++] = "./data/RFC3D_sphere_2.dat"; */
  /* data_collection_1[n_data++] = "./data/RFC3D_cube_1.dat"; */
  /* data_collection_1[n_data++] = "./data/rfc3d_sphere_from_grfc3d.dat"; */
  /* data_collection_1[n_data++] = "./data/rfc3d_sphere_from_grfc3d_dense.dat"; */
  data_collection_1[n_data++] = "./data/rockfall_14_25792.dat";
  data_collection_1[n_data++] = "./data/rockfall_21_54029.dat";
  data_collection_1[n_data++] = "./data/spheres_pile_270_8290.dat";
  data_collection_1[n_data++] = "---";

  return data_collection_1;
}
