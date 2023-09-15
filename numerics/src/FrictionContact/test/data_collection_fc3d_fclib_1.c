/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

  int n_data_1=100;

  const char ** data_collection_1 = (const char **)malloc(n_data_1*sizeof(const char *));
  int n_data=0;

  int listprob[2] = {1, 0};
  /* 0: LowWall_FEM      #50   problems */
  /* 1: ...              #.    problems */
  /*   TOTAL             #.    problems */

  if (listprob[0]==1)
  {
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00014-624-00001.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00016-675-00031.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00021-678-00014.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00022-688-00078.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00024-679-00032.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00024-688-00160.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00025-688-00122.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00025-688-00255.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00026-688-00054.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00026-688-00133.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00026-688-00277.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00027-688-00112.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00027-688-00310.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00028-688-00089.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00028-688-00172.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00028-688-00299.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00029-688-00161.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00029-688-00307.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00029-688-00472.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00030-688-00142.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00030-688-00210.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00030-688-00268.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00030-688-00381.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00031-688-00218.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00031-688-00325.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00031-688-00490.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00032-688-00191.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00032-688-00438.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00033-688-00189.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00033-688-00456.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00034-688-00226.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00034-688-00330.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00035-688-00173.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00035-688-00488.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00036-688-00308.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00037-688-00238.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00037-688-00419.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00038-688-00461.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00038-688-00498.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00039-688-00389.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00040-688-00350.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00040-688-00475.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00041-688-00474.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00042-688-00452.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00043-688-00429.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00044-688-00407.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00046-688-00344.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00047-688-00367.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00050-688-00222.hdf5";
    // data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00056-688-00272.hdf5";
  }

  else if (listprob[1]==1)
  {

  }

  data_collection_1[n_data++] = "---";
  return data_collection_1;
}
