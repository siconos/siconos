/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#include <stdlib.h>  // for malloc

#include "test_utils.h"  // for data_collection

const char** data_collection() {
  int n_data_1 = 150;

  const char** data_collection_1 = (const char**)malloc(n_data_1 * sizeof(const char*));
  int n_data = 0;
  data_collection_1[n_data++] = "./data/GFC3D_Example0.dat";
  data_collection_1[n_data++] = "./data/GFC3D_Example00.dat";
  data_collection_1[n_data++] = "./data/GFC3D_Example00_badly_scaled.dat";
  data_collection_1[n_data++] = "./data/GFC3D_Example1.dat";
  data_collection_1[n_data++] = "./data/GFC3D_OneContact.dat";
  data_collection_1[n_data++] = "./data/GFC3D_TwoRods1.dat";
  data_collection_1[n_data++] = "---";

  return data_collection_1;
}
