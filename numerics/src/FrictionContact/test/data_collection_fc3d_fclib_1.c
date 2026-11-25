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

  data_collection_1[n_data++] = "./data/Capsules-i125-1213.hdf5";
  data_collection_1[n_data++] = "./data/LMGC_100_PR_PerioBox-i00361-60-03000.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-library/Local/siconos/Chute_4000/Chute_4000-i0339-17651-52.hdf5";

  data_collection_1[n_data++] = "---";

  return data_collection_1;
}
