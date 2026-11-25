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

  data_collection_1[n_data++] = "./data/Rover1039.dat";
  data_collection_1[n_data++] = "./data/Rover1040.dat";
  data_collection_1[n_data++] = "./data/Rover1041.dat";
  data_collection_1[n_data++] = "./data/Rover3865.dat";
  data_collection_1[n_data++] = "./data/Rover4144.dat";
  data_collection_1[n_data++] = "./data/Rover4396.dat";
  data_collection_1[n_data++] = "./data/Rover4493.dat";
  data_collection_1[n_data++] = "./data/Rover4516.dat";
  data_collection_1[n_data++] = "./data/Rover4609.dat";
  data_collection_1[n_data++] = "./data/Rover4613.dat";
  data_collection_1[n_data++] = "./data/Rover4622.dat";
  data_collection_1[n_data++] = "./data/Rover9770.dat";
  data_collection_1[n_data++] = "./data/Rover11035.dat";
  data_collection_1[n_data++] = "./data/Rover11211.dat";

  data_collection_1[n_data++] = "---";

  return data_collection_1;
}
