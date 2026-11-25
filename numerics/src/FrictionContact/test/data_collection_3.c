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

  data_collection_1[n_data++] = "./data/FC3D_Example1_SBM.dat";
  data_collection_1[n_data++] = "./data/Capsules-i122-1617.dat";
  data_collection_1[n_data++] = "./data/Confeti-ex13-4contact-Fc3D-SBM.dat";
  data_collection_1[n_data++] = "./data/GFC3D_TwoRods1-condensed.dat";
  data_collection_1[n_data++] = "./data/FC3D_Example1.dat";
  data_collection_1[n_data++] = "./data/Confeti-ex03-Fc3D-SBM.dat"; /*5*/
  data_collection_1[n_data++] = "./data/BoxesStack1-i100000-32.hdf5.dat";
  data_collection_1[n_data++] = "./data/RockPile_tob1.dat";
  data_collection_1[n_data++] = "./data/KaplasTower-i1061-4.hdf5.dat";
  data_collection_1[n_data++] = "./data/OneObject-i100000-499.hdf5.dat";
  data_collection_1[n_data++] = "./data/Rover4396.dat"; /*10*/
  data_collection_1[n_data++] = "---";

  return data_collection_1;
}
