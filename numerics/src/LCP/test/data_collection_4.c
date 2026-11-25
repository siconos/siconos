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
  data_collection_1[n_data++] = "./data/lcp_Pang_isolated_sol.dat";
  data_collection_1[n_data++] = "./data/lcp_Pang_isolated_sol_perturbed.dat";
  data_collection_1[n_data++] = "./data/lcp_inf_sol_perturbed.dat";

  data_collection_1[n_data++] = "---";

  return data_collection_1;
}
