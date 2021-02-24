/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <stdlib.h>      // for malloc
#include "test_utils.h"  // for data_collection

const char ** data_collection()
{

  int n_data_1=150;

  const char ** data_collection_1 = (const char **)malloc(n_data_1*sizeof(const char *));
  int n_data=0;
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
