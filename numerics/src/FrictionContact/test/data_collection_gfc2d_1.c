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
  data_collection_1[n_data++] = "./data/gfc2d_granularflowonwall_00000_00117.dat";
  data_collection_1[n_data++] = "./data/gfc2d_granularflowonwall_00069_00350.dat";
  data_collection_1[n_data++] = "./data/gfc2d_granularflowonwall_00230_00206.dat";
  data_collection_1[n_data++] = "./data/GFC3D_Example0_2d.dat";
  data_collection_1[n_data++] = "./data/GFC3D_Example1_2d.dat";
  data_collection_1[n_data++] = "./data/GFC3D_OneContact_2d.dat";
  data_collection_1[n_data++] = "---";
  return data_collection_1;
}
