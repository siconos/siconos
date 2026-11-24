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
  data_collection_1[n_data++] = "./data/FrictionContactProblem00374.dat";
  data_collection_1[n_data++] = "./data/FrictionContactProblem00394.dat";

  data_collection_1[n_data++] = "./data/FC2D_SliderCrankLagrangian00000.dat";
  data_collection_1[n_data++] = "./data/FC2D_SliderCrankLagrangian00001.dat";

  /* ## test 2D sparse block on 4 differents files */
  data_collection_1[n_data++] = "./data/FrictionContactProblem00727.dat";
  data_collection_1[n_data++] = "./data/FrictionContactProblem00031.dat";
  data_collection_1[n_data++] = "./data/FrictionContactProblem00071.dat";
  data_collection_1[n_data++] = "./data/FrictionContactProblem00237.dat"; /* 5 */

  /* ## test 2D sparse  on 4 differents files */
  data_collection_1[n_data++] = "./data/fc2d_granularflowonwall_00667_00000.dat";
  data_collection_1[n_data++] = "./data/fc2d_granularflowonwall_02020_00213.dat";
  data_collection_1[n_data++] = "./data/fc2d_granularflowonwall_07002_00801.dat";

  data_collection_1[n_data++] = "---";

  return data_collection_1;
}
