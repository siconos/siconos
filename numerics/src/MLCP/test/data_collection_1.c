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
  data_collection_1[n_data++] = "./data/deudeu_mlcp.dat";
  data_collection_1[n_data++] = "./data/RLCD_mlcp.dat";
  data_collection_1[n_data++] = "./data/RCD_mlcp.dat";

  data_collection_1[n_data++] = "./data/relay3_mlcp.dat";
  data_collection_1[n_data++] = "./data/relay2_mlcp.dat";

  data_collection_1[n_data++] = "./data/m2n1_mlcp.dat";
  data_collection_1[n_data++] = "./data/m3n2_mlcp.dat";

  data_collection_1[n_data++] = "./data/PDSym_mlcp.dat";
  data_collection_1[n_data++] = "./data/PD_mlcp.dat";

  /* data_collection_1[n_data++] = "./data/Rectangular_mlcp.dat"; */

  data_collection_1[n_data++] = "./data/diodeBridge20_mlcp.dat";
  /* data_collection_1[n_data++] = "./data/diodeBridge40_mlcp.dat"; */
  data_collection_1[n_data++] = "./data/diodeBridge_mlcp.dat";

  data_collection_1[n_data++] = "./data/Buck2_mlcp.dat";
  data_collection_1[n_data++] = "./data/BuckFirstStep_mlcp.dat";
  data_collection_1[n_data++] = "./data/BuckConverter_mlcp.dat";

  /* data_collection_1[n_data++] = "./data/deltasigma_mlcp.dat"; */
  /* data_collection_1[n_data++] = "./data/deltasigma2_mlcp.dat"; */
  data_collection_1[n_data++] = "---";

  return data_collection_1;
}
