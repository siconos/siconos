/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include <stdlib.h>  // for malloc

#include "test_utils.h"  // for data_collection

const char **data_collection() {
  int n_data_1 = 400;

  const char **data_collection_1 = (const char **)malloc(n_data_1 * sizeof(const char *));
  int n_data = 0;

  // List of rolling problems
  // int listprob[3] = {1, 0, 0};
  int listprob[3] = {1, 1, 1};
  /* 0: Chute          # 39 problems */  // 161 problems does not satisfy Slater's conditions
  /* 1: PrimitiveSoup  #158 problems */  //  22 problems does not satisfy Slater's conditions
  /* 2: SpheresPile    #188 problems */
  /*   TOTAL           #385 problems */

  if (listprob[0] == 1) {
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-768-nc-4-3.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-768-nc-18-28.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-768-nc-20-38.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-768-nc-27-50.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-768-nc-29-54.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-768-nc-31-61.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-768-nc-32-64.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-768-nc-33-67.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1152-nc-55-126.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1152-nc-58-118.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1152-nc-63-138.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1152-nc-69-150.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1152-nc-79-180.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1206-nc-88-195.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1236-nc-93-199.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1260-nc-90-203.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1362-nc-98-213.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1536-nc-155-283.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1536-nc-161-292.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1536-nc-164-298.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1536-nc-172-296.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1590-nc-178-305.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1722-nc-196-327.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1746-nc-200-331.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1758-nc-204-333.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1920-nc-203-355.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1920-nc-206-364.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1920-nc-214-366.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1920-nc-223-376.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1920-nc-227-375.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1920-nc-239-391.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1920-nc-267-428.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1920-nc-275-431.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-1974-nc-270-438.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-2064-nc-273-451.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-2154-nc-282-461.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-2166-nc-285-462.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-2232-nc-287-470.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-2304-nc-321-501.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-2304-nc-330-517.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-2304-nc-334-515.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-2304-nc-339-516.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-2358-nc-365-558.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-2388-nc-354-561.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-2394-nc-351-563.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-2478-nc-358-568.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-2622-nc-383-592.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-2622-nc-385-591.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-2688-nc-388-602.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-2688-nc-420-620.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-2688-nc-424-621.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-2688-nc-429-636.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-2688-nc-430-630.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-2688-nc-445-644.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-2718-nc-457-678.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-2886-nc-468-693.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-2928-nc-476-701.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3006-nc-488-708.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3036-nc-487-711.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3072-nc-506-729.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3072-nc-518-730.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3072-nc-519-757.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3072-nc-522-763.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3072-nc-526-749.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3072-nc-531-765.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3084-nc-566-798.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3234-nc-576-815.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3282-nc-579-819.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3360-nc-612-824.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3444-nc-599-829.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3450-nc-609-831.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3450-nc-610-832.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3456-nc-605-844.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3456-nc-627-871.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3456-nc-629-867.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3456-nc-641-864.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3588-nc-651-909.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3672-nc-677-918.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3690-nc-658-920.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3840-nc-709-959.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3840-nc-712-953.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3840-nc-727-974.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3840-nc-734-996.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3840-nc-737-980.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3840-nc-750-997.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-3990-nc-767-1045.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4026-nc-780-1053.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4212-nc-827-1073.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4224-nc-800-1094.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4224-nc-803-1082.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4224-nc-806-1080.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4224-nc-809-1104.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4224-nc-811-1089.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4224-nc-814-1100.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4224-nc-815-1093.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4224-nc-845-1143.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4224-nc-869-1159.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4296-nc-834-1172.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4464-nc-853-1208.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4608-nc-872-1233.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4608-nc-876-1230.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4608-nc-902-1255.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4608-nc-908-1269.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4608-nc-911-1250.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4608-nc-915-1271.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4608-nc-919-1267.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4608-nc-932-1277.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4608-nc-934-1286.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4638-nc-938-1311.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4782-nc-916-1324.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4920-nc-950-1341.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4992-nc-951-1424.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4992-nc-960-1374.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4992-nc-963-1375.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-4992-nc-989-1400.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5136-nc-969-1447.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5316-nc-993-1472.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5376-nc-994-1482.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5376-nc-998-1501.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5376-nc-1008-1495.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5376-nc-1012-1503.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5376-nc-1032-1525.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5376-nc-1033-1561.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5376-nc-1036-1534.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5376-nc-1073-1536.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5424-nc-1072-1576.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5622-nc-1101-1601.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5724-nc-1119-1614.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5760-nc-1108-1620.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5760-nc-1110-1625.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5760-nc-1132-1634.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5760-nc-1160-1658.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5760-nc-1173-1664.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5760-nc-1179-1654.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5760-nc-1181-1647.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5760-nc-1191-1666.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5760-nc-1192-1668.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5760-nc-1214-1681.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5760-nc-1233-1675.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5904-nc-1259-1707.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5940-nc-1249-1718.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-5946-nc-1234-1719.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-6120-nc-1289-1741.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-6144-nc-1281-1756.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-6144-nc-1302-1763.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-6144-nc-1316-1797.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-6144-nc-1322-1803.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-6240-nc-1326-1829.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-6246-nc-1329-1830.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-6288-nc-1298-1834.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-6438-nc-1341-1844.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-6528-nc-1353-1865.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-6528-nc-1360-1868.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-6528-nc-1361-1873.hdf5";  // pinfeas_2
    data_collection_1[n_data++] =
        "./fclib-library/GlobalRolling/siconos/Chute/Chute-ndof-6528-nc-1372-1870.hdf5";
  }
  if (listprob[0] == 1) {
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-768-nc-13-4.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-768-nc-47-29.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-768-nc-57-43.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-768-nc-58-40.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-984-nc-81-79.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-1152-nc-89-101.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-1152-nc-91-108.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-1152-nc-98-119.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-1152-nc-112-142.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-1152-nc-128-153.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-1182-nc-141-179.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-1536-nc-156-226.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-1536-nc-167-243.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-1536-nc-168-248.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-1536-nc-173-255.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-1536-nc-175-256.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-1536-nc-202-286.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-1536-nc-207-285.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-1746-nc-239-321.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-1770-nc-232-326.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-1842-nc-253-333.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-1920-nc-263-348.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-1920-nc-300-388.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-1920-nc-327-421.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-1920-nc-335-415.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-2052-nc-328-437.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-2148-nc-330-455.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-2214-nc-329-462.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-2274-nc-343-468.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-2304-nc-387-508.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-2334-nc-417-546.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-2334-nc-418-547.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-2340-nc-415-549.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-2394-nc-403-560.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-2688-nc-437-602.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-2688-nc-445-600.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-2688-nc-451-617.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-2688-nc-453-619.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-2688-nc-462-621.hdf5";

    // /* ========================================== Slater's conditions are not satisfied
    // ==========================================*/ data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-2718-nc-485-679.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-2910-nc-504-702.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3072-nc-538-733.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3072-nc-561-754.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3072-nc-565-762.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3072-nc-567-768.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3072-nc-568-758.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3072-nc-569-771.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3072-nc-574-767.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3072-nc-580-795.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3072-nc-585-787.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3072-nc-590-786.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3120-nc-612-812.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3120-nc-616-811.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3348-nc-645-834.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3456-nc-638-851.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3456-nc-660-861.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3456-nc-662-867.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3456-nc-664-885.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3456-nc-665-868.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3456-nc-706-912.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3840-nc-756-967.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3840-nc-763-964.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3840-nc-774-976.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3840-nc-777-975.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3840-nc-801-1035.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3864-nc-797-1045.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3900-nc-811-1052.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-3900-nc-815-1053.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4068-nc-851-1080.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4224-nc-840-1104.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4224-nc-859-1129.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4224-nc-865-1131.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4224-nc-899-1167.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4224-nc-905-1171.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4224-nc-911-1186.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4224-nc-912-1185.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4230-nc-916-1194.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4398-nc-936-1213.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4416-nc-937-1223.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4416-nc-953-1220.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4608-nc-975-1315.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4608-nc-976-1306.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4608-nc-978-1294.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4608-nc-984-1260.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4812-nc-1021-1354.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4854-nc-1025-1360.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4884-nc-1019-1366.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4920-nc-1042-1368.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4992-nc-1055-1385.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4992-nc-1059-1387.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4992-nc-1091-1413.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4992-nc-1092-1415.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4992-nc-1106-1410.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4992-nc-1114-1419.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4992-nc-1116-1416.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4992-nc-1136-1426.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4992-nc-1141-1433.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-4992-nc-1154-1435.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5166-nc-1165-1483.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5262-nc-1208-1501.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5292-nc-1193-1507.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5376-nc-1226-1515.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5376-nc-1229-1530.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5376-nc-1233-1539.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5376-nc-1285-1572.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5376-nc-1304-1550.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5376-nc-1312-1589.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5376-nc-1313-1590.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5376-nc-1327-1585.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5502-nc-1343-1612.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5508-nc-1353-1615.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5526-nc-1341-1623.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5538-nc-1340-1626.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5556-nc-1348-1627.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5760-nc-1363-1657.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5760-nc-1386-1666.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5760-nc-1391-1659.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5760-nc-1398-1672.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5760-nc-1401-1687.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5760-nc-1414-1700.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5760-nc-1415-1705.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5760-nc-1428-1723.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5760-nc-1436-1712.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5958-nc-1470-1744.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-5964-nc-1460-1745.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6120-nc-1533-1768.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6126-nc-1522-1770.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6144-nc-1504-1777.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6144-nc-1510-1773.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6144-nc-1511-1776.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6144-nc-1513-1785.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6144-nc-1517-1780.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6144-nc-1537-1821.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6144-nc-1553-1808.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6186-nc-1543-1840.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6312-nc-1600-1857.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6312-nc-1610-1858.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6324-nc-1605-1860.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6480-nc-1637-1876.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6528-nc-1683-1896.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6528-nc-1694-1905.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6528-nc-1699-1919.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6528-nc-1703-1920.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6528-nc-1709-1915.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6528-nc-1710-1931.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6528-nc-1718-1938.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6528-nc-1740-1954.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6528-nc-1766-1944.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6624-nc-1763-1986.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6624-nc-1771-1987.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6636-nc-1783-1988.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6702-nc-1788-1996.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6786-nc-1805-2002.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6852-nc-1853-2011.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6912-nc-1843-2044.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6912-nc-1855-2051.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6912-nc-1857-2052.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6912-nc-1858-2057.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6912-nc-1871-2059.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6912-nc-1877-2054.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6912-nc-1882-2070.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6912-nc-1890-2077.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6912-nc-1921-2082.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-6912-nc-1925-2089.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7224-nc-1922-2129.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7236-nc-1942-2133.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7284-nc-1933-2139.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7296-nc-1930-2153.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7296-nc-1951-2143.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7296-nc-1953-2187.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7296-nc-1975-2190.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7296-nc-2008-2194.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7326-nc-2049-2225.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7392-nc-2048-2232.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7464-nc-2058-2245.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7620-nc-2070-2262.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7668-nc-2079-2269.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7680-nc-2065-2271.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7680-nc-2094-2294.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7680-nc-2118-2309.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7680-nc-2122-2322.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7680-nc-2135-2329.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7680-nc-2140-2334.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7680-nc-2146-2318.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7680-nc-2166-2342.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7722-nc-2187-2348.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7830-nc-2200-2369.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7902-nc-2222-2375.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7914-nc-2189-2376.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-7950-nc-2215-2380.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-8040-nc-2252-2390.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-8064-nc-2237-2397.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-8064-nc-2244-2398.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-8064-nc-2254-2400.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-8064-nc-2262-2409.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-8064-nc-2307-2457.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-8064-nc-2317-2452.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-8064-nc-2324-2462.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-8064-nc-2335-2463.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/Chute/Chute-ndof-8064-nc-2343-2447.hdf5";
  }

  if (listprob[1] == 1) {
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-96-0.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-141-1.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-193-3.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-198-4.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-204-6.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-205-7.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-211-8.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-214-2.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-215-12.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-216-10.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-217-11.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-242-13.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-263-16.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-268-18.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-273-20.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-309-21.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-350-22.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-359-24.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-381-25.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-384-26.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-400-29.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-410-32.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-414-31.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-429-33.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-439-35.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-440-34.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-450-36.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-466-39.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-468-38.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-493-41.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-499-42.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-504-43.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-510-45.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-513-44.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-518-46.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-519-47.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-547-48.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-563-49.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-575-50.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-595-51.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-609-52.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-636-54.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-637-53.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-694-55.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-699-56.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-727-57.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-739-58.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-792-59.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-805-60.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-826-61.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-840-62.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-852-63.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-854-64.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-884-65.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-899-66.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-903-67.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-942-68.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-990-70.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-995-69.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1022-71.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1031-73.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1034-74.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1143-78.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1156-186.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1199-81.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1207-188.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1222-197.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1226-194.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1228-195.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1229-191.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1231-193.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1250-183.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1251-201.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1252-199.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1253-184.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1256-164.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1258-165.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1261-174.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1297-205.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1298-231.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1304-227.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1327-154.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1328-257.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1346-267.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1354-86.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1358-260.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1363-138.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1364-250.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1373-254.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1396-132.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1433-292.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1450-284.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1463-282.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1472-302.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1475-283.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1491-117.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1499-120.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1505-112.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1518-306.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1536-311.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1538-107.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1548-319.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1551-106.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1591-334.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1607-341.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1608-337.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1622-343.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1642-345.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1647-352.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1653-346.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1654-348.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1656-344.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1668-351.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1670-350.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1701-356.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1745-365.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1768-376.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1769-374.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1774-370.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1781-371.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1782-378.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1793-379.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1840-389.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1843-390.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1844-396.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1857-382.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1858-397.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1864-395.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1875-403.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1881-399.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1882-405.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1912-409.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1930-410.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1963-411.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-1970-412.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-2005-418.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-2017-421.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-2018-419.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-2064-438.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-2079-442.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-2093-444.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-2119-449.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-2122-445.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-2124-461.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-2128-446.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-2150-463.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-2163-470.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-2165-468.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-2166-473.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-2182-474.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-2192-476.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-2196-480.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-2199-477.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-2217-487.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-2237-490.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-2245-504.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-2252-501.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/"
        "PrimitiveSoup-ndof-6000-nc-2269-511.hdf5";

    // /* ========================================== Slater's conditions are not satisfied
    // ==========================================*/ data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1237-182.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2329-546.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2372-575.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2390-581.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2394-580.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2398-599.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2409-598.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2418-595.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2419-606.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2421-639.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2430-613.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2432-618.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2437-612.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2448-646.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2452-642.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2463-655.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2464-654.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2475-660.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2501-667.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2517-703.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2521-706.hdf5";
    // data_collection_1[n_data++] =
    // "./fclib-library/RollingGlobal/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2530-710.hdf5";
  }

  if (listprob[2] == 1) {
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-138-nc-25-0.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-168-nc-30-1.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-174-nc-34-2.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-198-nc-41-3.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-198-nc-43-4.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-204-nc-42-5.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-228-nc-48-6.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-246-nc-61-7.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-270-nc-59-15.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-270-nc-63-14.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-270-nc-66-13.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-276-nc-69-17.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-306-nc-78-26.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-306-nc-79-25.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-312-nc-80-29.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-312-nc-82-28.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-318-nc-83-32.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-318-nc-86-33.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-324-nc-88-34.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-330-nc-91-38.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-330-nc-93-40.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-342-nc-95-41.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-342-nc-100-42.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-360-nc-102-55.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-378-nc-96-62.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-390-nc-104-66.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-402-nc-111-71.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-402-nc-114-73.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-408-nc-120-75.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-444-nc-126-89.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-450-nc-130-96.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-462-nc-134-106.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-462-nc-136-104.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-480-nc-146-114.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-486-nc-141-116.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-492-nc-143-119.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-498-nc-144-121.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-504-nc-150-123.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-516-nc-155-128.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-522-nc-166-130.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-534-nc-145-134.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-546-nc-148-137.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-552-nc-158-139.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-558-nc-153-146.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-558-nc-157-147.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-564-nc-164-149.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-570-nc-167-154.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-576-nc-172-158.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-576-nc-175-160.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-582-nc-178-161.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-594-nc-174-168.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-600-nc-182-169.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-606-nc-186-171.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-618-nc-193-174.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-618-nc-200-177.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-624-nc-195-180.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-630-nc-190-185.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-630-nc-199-182.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-636-nc-197-190.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-636-nc-204-189.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-642-nc-202-191.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-648-nc-211-196.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-654-nc-209-201.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-660-nc-210-203.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-666-nc-213-209.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-666-nc-216-206.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-696-nc-218-219.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-720-nc-231-229.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-732-nc-237-235.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-768-nc-251-259.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-780-nc-250-264.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-780-nc-267-265.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-786-nc-252-267.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-786-nc-257-268.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-792-nc-261-273.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-798-nc-259-274.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-798-nc-270-275.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-810-nc-262-282.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-816-nc-265-285.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-828-nc-266-295.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-834-nc-272-297.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-840-nc-274-299.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-846-nc-279-301.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-864-nc-283-308.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-870-nc-277-311.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-876-nc-271-314.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-876-nc-275-316.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-882-nc-278-318.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-882-nc-282-322.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-888-nc-286-323.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-888-nc-288-324.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-894-nc-298-326.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-906-nc-306-330.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-912-nc-309-336.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-912-nc-310-333.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-912-nc-312-334.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-918-nc-304-340.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-918-nc-305-341.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-918-nc-314-342.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-930-nc-315-344.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-942-nc-307-350.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-948-nc-303-352.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-948-nc-311-354.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-972-nc-333-361.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-978-nc-334-366.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-978-nc-336-363.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-984-nc-328-370.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-990-nc-335-372.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1014-nc-348-387.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1020-nc-358-392.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1038-nc-353-402.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1056-nc-345-406.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1068-nc-359-411.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1080-nc-357-418.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1086-nc-363-424.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1086-nc-367-423.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1086-nc-371-420.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1110-nc-375-442.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1116-nc-366-443.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1134-nc-384-450.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1146-nc-378-458.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1158-nc-372-465.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1164-nc-383-469.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1164-nc-390-471.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1164-nc-397-468.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1170-nc-392-473.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1176-nc-387-476.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1176-nc-388-479.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1182-nc-380-481.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1206-nc-399-488.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1212-nc-396-495.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1212-nc-398-491.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1218-nc-406-499.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1218-nc-409-500.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1224-nc-394-505.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1230-nc-411-510.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1236-nc-405-512.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1248-nc-402-518.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1248-nc-408-517.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1254-nc-403-521.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1260-nc-416-525.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1284-nc-433-535.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1302-nc-441-542.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1314-nc-437-547.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1320-nc-444-549.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1344-nc-447-559.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1356-nc-445-562.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1362-nc-450-565.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1362-nc-459-563.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1362-nc-461-566.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1380-nc-449-574.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1380-nc-460-570.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1386-nc-466-579.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1386-nc-467-575.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1386-nc-469-577.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1392-nc-458-584.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1392-nc-472-583.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1398-nc-482-588.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1404-nc-478-592.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1416-nc-486-597.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1416-nc-490-596.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1422-nc-491-599.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1440-nc-493-607.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1440-nc-501-608.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1446-nc-495-611.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1446-nc-496-609.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1452-nc-511-616.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1452-nc-515-614.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1458-nc-510-618.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1470-nc-513-623.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1482-nc-492-626.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1488-nc-512-630.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1494-nc-521-632.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1500-nc-519-634.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1500-nc-533-640.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1500-nc-534-654.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1500-nc-536-645.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1500-nc-541-673.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1500-nc-542-661.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1500-nc-543-719.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1500-nc-545-675.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1500-nc-547-706.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1500-nc-548-734.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1500-nc-553-808.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1500-nc-555-761.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1500-nc-556-786.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1500-nc-559-799.hdf5";
    data_collection_1[n_data++] =
        "./fclib-library/RollingGlobal/siconos/SpheresPile/"
        "SpheresPile-ndof-1500-nc-560-816.hdf5";
  }

  data_collection_1[n_data++] = "---";
  return data_collection_1;
}
