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

  /* Small Drucker-Prager problems - only use small tests that complete quickly */
  data_collection_1[n_data++] = "./data/plasticity_2d_example1.dat";
  data_collection_1[n_data++] = "./data/plasticity_2d_example1_theta0.dat";
  data_collection_1[n_data++] = "./data/plasticity_2d_footing_1.dat";
  /* data_collection_1[n_data++] = "./data/plasticity_2d_footing_10.dat"; */
  /* data_collection_1[n_data++] = "./data/plasticity_2d_footing_100_theta0.05.dat";   */
  /* data_collection_1[n_data++] = "./data/plasticity_2d_footing_100.dat"; */
  /* data_collection_1[n_data++] = "./data/plasticity_2d_footing_100_theta0.dat"; */

  data_collection_1[n_data++] = "---";

  return data_collection_1;
}
