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

  int listprob[2] = {0, 1};
  /* 0: LowWall_FEM      #50   problems */
  /* 1: Aqueduc          #10    problems */
  /*   TOTAL             #.    problems */

  if (listprob[0]==1)
  {
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00014-624-00001.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00016-675-00031.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00021-678-00014.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00022-688-00078.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00024-679-00032.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00024-688-00160.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00025-688-00122.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00025-688-00255.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00026-688-00054.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00026-688-00133.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00026-688-00277.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00027-688-00112.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00027-688-00310.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00028-688-00089.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00028-688-00172.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00028-688-00299.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00029-688-00161.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00029-688-00307.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00029-688-00472.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00030-688-00142.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00030-688-00210.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00030-688-00268.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00030-688-00381.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00031-688-00218.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00031-688-00325.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00031-688-00490.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00032-688-00191.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00032-688-00438.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00033-688-00189.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00033-688-00456.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00034-688-00226.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00034-688-00330.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00035-688-00173.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00035-688-00488.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00036-688-00308.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00037-688-00238.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00037-688-00419.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00038-688-00461.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00038-688-00498.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00039-688-00389.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00040-688-00350.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00040-688-00475.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00041-688-00474.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00042-688-00452.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00043-688-00429.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00044-688-00407.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00046-688-00344.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00047-688-00367.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00050-688-00222.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/LowWall_FEM/LMGC_LowWall_FEM-i00056-688-00272.hdf5";
  }

  else if (listprob[1]==1)
  {
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Aqueduct_PR/LMGC_AqueducPR-i00001-4387-00025.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Aqueduct_PR/LMGC_AqueducPR-i00001-4471-00013.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Aqueduct_PR/LMGC_AqueducPR-i00001-4476-00048.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Aqueduct_PR/LMGC_AqueducPR-i00001-4811-00001.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Aqueduct_PR/LMGC_AqueducPR-i04236-4455-00009.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Aqueduct_PR/LMGC_AqueducPR-i04320-4453-00010.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Aqueduct_PR/LMGC_AqueducPR-i04404-4467-00011.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Aqueduct_PR/LMGC_AqueducPR-i04489-4477-00012.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Aqueduct_PR/LMGC_AqueducPR-i04632-4441-00008.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Aqueduct_PR/LMGC_AqueducPR-i05000-4454-00007.hdf5";
  }

  else if (listprob[1]==2)
  {
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_2-i00001-5-00031.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_2-i00004-3-00024.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_2-i00005-3-00003.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_2-i00007-3-00043.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_2-i00008-3-00035.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_2-i00009-3-00080.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_2-i00011-2-00034.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_2-i00012-3-00020.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_2-i00012-4-00067.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_2-i00013-2-00079.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_2-i00013-3-00112.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_2-i00013-5-00065.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_2-i00014-3-00123.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_2-i00015-3-00122.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_2-i00019-5-00010.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i01249-361-00143.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i01688-363-00029.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i01851-361-00147.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i01910-361-00088.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02115-361-00056.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02228-363-00039.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02247-366-00140.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02344-361-00188.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02377-364-00007.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02411-388-00002.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02458-361-00123.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02482-361-00092.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02563-365-00057.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02613-361-00040.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02655-361-00189.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02696-363-00041.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02716-373-00087.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02726-365-00102.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02760-361-00176.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02767-363-00043.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02799-366-00103.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02836-365-00124.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02842-363-00003.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02884-361-00196.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02890-361-00084.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02908-361-00090.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02913-368-00095.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02933-369-00079.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02939-361-00018.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i02967-365-00116.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i03011-366-00138.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i03036-361-00153.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i03068-361-00149.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i03170-361-00014.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i03318-361-00074.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i03428-364-00182.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i03549-361-00183.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i03713-368-00162.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i03802-361-00131.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i03849-361-00194.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i03873-361-00177.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i03898-366-00150.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i03944-361-00165.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i04004-365-00104.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i04068-361-00115.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i04116-361-00167.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i04189-361-00022.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Cubes_H8/LMGC_Cubes_H8_20-i04284-361-00129.hdf5";
  }

  else if (listprob[1]==3)
  {
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i00157-76-00640.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i00258-71-00810.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i00265-72-00770.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i00290-71-00850.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i00340-70-00790.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i00344-70-00830.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i00668-72-00740.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i01548-108-00020.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i01573-108-00040.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i01635-108-00060.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i02219-90-00560.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i02864-106-00260.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i02940-106-00240.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i02960-106-00220.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i02972-106-00200.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i02982-106-00170.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i02987-106-00180.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i03005-106-00150.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i03121-106-00290.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i03276-106-00120.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i04009-72-00750.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-100-00470.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-100-00490.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-100-00510.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-100-00530.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-102-00430.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-102-00450.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-106-00100.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-106-00310.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-106-00330.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-106-00350.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-106-00370.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-106-00390.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-106-00410.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-108-00080.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-72-00700.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-72-00720.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-72-00880.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-74-00590.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-74-00610.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-74-00630.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-74-00900.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-74-00920.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-74-01000.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-76-00660.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-76-00680.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-76-00970.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-76-00990.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-78-00570.hdf5";
    data_collection_1[n_data++] = "./fclib-library/Local/LMGC90/Bridge_PR/LMGC_Bridge_PR-i05000-78-00930.hdf5";
  }





  data_collection_1[n_data++] = "---";
  return data_collection_1;
}
