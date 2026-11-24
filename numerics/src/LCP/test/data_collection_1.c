// SPDX-FileCopyrightText: 2026 INRIA
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
  data_collection_1[n_data++] = "./data/lcp_mmc.dat";
  data_collection_1[n_data++] = "./data/lcp_deudeu.dat";
  data_collection_1[n_data++] = "./data/lcp_trivial.dat";
  data_collection_1[n_data++] = "./data/lcp_ortiz.dat";
  data_collection_1[n_data++] = "./data/lcp_exp_murty.dat";
  data_collection_1[n_data++] = "./data/lcp_exp_murty2.dat";

  data_collection_1[n_data++] = "---";

  return data_collection_1;
}
