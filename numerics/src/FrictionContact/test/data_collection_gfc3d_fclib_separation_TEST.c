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

#include <stdlib.h>  // for malloc

#include "test_utils.h"  // for data_collection

const char **data_collection() {
  int n_data_1 = 30000;

  const char **data_collection_1 = (const char **)malloc(n_data_1 * sizeof(const char *));
  int n_data = 0;

  data_collection_1
      [n_data++] =
          "./fclib-sub/Global/siconos/PrimitiveSoup/"
          "PrimitiveSoup-ndof-6000-nc-1087-183-68.hdf5";  // n=11
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-13824-nc-3722-4252-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-155-18120-39.hdf5";

  // data_collection_1[n_data++] =
  // "./fclib-library/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-103-729.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-103-729-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-library/Global/siconos/Chute/Chute-ndof-6816-nc-859-1963.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-library/Global/siconos/Chute/Chute-ndof-5376-nc-434-1550.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-6816-nc-859-1963-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3033-8373-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-244-20980-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-110-14908-26.hdf5"; // n=5, BNR
  // = 2/0/3 data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-2730-6-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4383-18-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-240-201-53.hdf5";

  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-18138-nc-6060-5654-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4627-1845-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4625-1864-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3194-2081-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4517-49149-1.hdf5";

  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-5-2299-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-2-2505-1.hdf5";

  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-42-10505-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-6528-nc-815-1897-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-393-710-119.hdf5";

  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-14.hdf5";    // n=10
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1230-223-62.hdf5"; //
  // n=2 data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-10752-nc-2385-3278-1.hdf5";        // n =
  // 2109 data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-5376-nc-434-1550-3.hdf5";

  // // RESTART
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-110-14908-26.hdf5"; // n=5, BNR
  // = 2/0/3 data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-130-16647-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-155-18120-39.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-156-17670-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-172-18659-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-179-19746-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-197-19889-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-243-21083-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-280-23762-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-283-23857-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-284-23843-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-325-32051-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-355-36003-34.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-103-729-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-466-2477-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-482-2537-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-489-2484-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-506-1976-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-544-2626-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-633-3770-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-651-3596-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-656-2866-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-671-2887-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-712-4000-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-715-11607-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-766-8284-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-597-48-94.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-897-132-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1003-161-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1087-183-68.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1122-204-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1264-227-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1475-295-79.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1504-298-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1556-306-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1611-330-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1649-336-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1729-343-98.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1768-366-93.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1781-357-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1787-367-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1798-364-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1825-369-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1857-389-105.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1860-386-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1909-382-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1934-391-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1962-398-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1995-412-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2045-431-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2045-431-32.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2077-427-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2082-445-64.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2088-437-64.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2090-440-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2174-457-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2214-488-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2222-461-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2239-485-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2250-473-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2256-472-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2273-499-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2279-479-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2283-502-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2332-518-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2369-523-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2371-526-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2388-531-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2438-544-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2460-578-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2470-580-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2488-575-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2489-568-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2491-564-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2498-586-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2517-590-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2521-594-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2537-611-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2537-611-120.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2545-622-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2597-636-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2612-641-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2616-648-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2633-643-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2639-651-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2708-701-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2720-703-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2722-717-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2724-713-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2734-712-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2893-938-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2902-990-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2920-1003-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2923-1005-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2925-1012-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2926-1006-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2999-1148-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3004-1135-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3047-1333-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3059-1383-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3066-1595-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3141-1907-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3144-1937-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3147-1925-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3154-1968-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3179-2042-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3190-2112-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3194-2081-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3199-2147-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3221-2232-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3235-2535-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3246-2536-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3249-2486-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3256-2494-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3271-2771-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3277-2770-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3282-2721-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-789-1998-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-796-2019-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-1124-3390-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-2141-6186-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-2390-6741-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-2434-6864-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-2486-6930-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-2658-7448-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-2773-7698-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-2840-7935-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-2963-8251-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-2967-8115-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3002-8309-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3013-8445-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3033-8373-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3038-8510-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3040-8470-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3045-8393-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3094-8679-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3105-8750-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3134-8838-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3143-8862-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3162-8944-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3172-8890-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3209-9146-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3234-9251-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3372-9824-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3381-9741-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3419-9867-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3432-9890-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3458-10018-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3468-10034-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3521-10193-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3525-10094-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3527-10170-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3564-10244-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3602-10372-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3605-10519-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3622-10418-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3700-10899-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3712-10989-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3717-10927-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3728-11062-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3737-11065-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3744-11009-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3812-11179-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3895-11737-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3918-11500-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3920-11602-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3928-11551-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-3960-11806-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4023-12160-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4070-12579-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4077-12622-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4085-13149-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4157-13366-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4159-13705-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4164-13707-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4168-13718-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4169-13598-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4199-13842-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4237-15402-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4250-15000-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4255-15460-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4271-15612-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4280-15670-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4334-16513-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4346-18262-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4371-17786-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4375-18220-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4435-26334-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4440-33607-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4454-32008-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4461-32433-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4464-27544-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4466-31772-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4482-38823-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4517-49149-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-2407-3-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-2444-5-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-2730-6-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-2941-7-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-3120-8-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-3284-9-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-3389-10-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-3459-11-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-3940-12-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4025-13-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4144-14-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4369-17-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4383-18-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4430-20-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4432-21-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4451-85-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4464-23-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4481-105-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4484-70-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4487-76-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4488-135-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4489-114-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4490-51-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4500-148-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4506-172-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4515-182-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4520-251-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4536-377-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4547-441-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4557-720-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4559-682-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4564-791-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4573-845-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4575-901-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4581-866-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4603-1534-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4610-1948-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4613-1698-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4615-1898-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4620-1879-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4625-1864-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4627-1845-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-5376-nc-434-1550-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-6144-nc-654-1784-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-6816-nc-859-1963-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-6912-nc-965-2050-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-6942-nc-979-2057-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-6942-nc-979-2057-31.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-7296-nc-1068-2169-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-7296-nc-1121-2138-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-7614-nc-1192-2224-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-7680-nc-1226-2312-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-7680-nc-1226-2312-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-7692-nc-1247-2315-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-8412-nc-1369-2481-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-8448-nc-1408-2507-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-8448-nc-1412-2509-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-8448-nc-1436-2501-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-8448-nc-1436-2501-32.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-8574-nc-1557-2575-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-8604-nc-1556-2580-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-9216-nc-1885-2782-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-9216-nc-1900-2778-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-9216-nc-1957-2807-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-9600-nc-2023-2861-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-9600-nc-2037-2879-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-9600-nc-2065-2894-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-10752-nc-2385-3278-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-10752-nc-2459-3293-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-10752-nc-2514-3305-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-10872-nc-2558-3334-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-11082-nc-2572-3351-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-11100-nc-2574-3352-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-11136-nc-2604-3371-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-11136-nc-2657-3373-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-11136-nc-2716-3407-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-11136-nc-2750-3414-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-11136-nc-2758-3424-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-11184-nc-2797-3439-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-11304-nc-2865-3458-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-11520-nc-2871-3523-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-11520-nc-2873-3515-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-11520-nc-2888-3508-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-11520-nc-2944-3551-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-11520-nc-2965-3556-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-11520-nc-3015-3505-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-11610-nc-3066-3579-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-11904-nc-3025-3652-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-11904-nc-3073-3675-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-12012-nc-3018-3692-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-12012-nc-3018-3692-39.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-12090-nc-3167-3702-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-12192-nc-3091-3717-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-12270-nc-3101-3728-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-12288-nc-3130-3742-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-12354-nc-3188-3832-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-12534-nc-3157-3853-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-12648-nc-3178-3869-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-12672-nc-3248-3892-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-12672-nc-3271-3893-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-12690-nc-3210-3939-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-13056-nc-3278-4017-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-13056-nc-3509-4052-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-13194-nc-3513-4089-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-13194-nc-3513-4089-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-13200-nc-3482-4093-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-13440-nc-3561-4133-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-13440-nc-3572-4190-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-13614-nc-3613-4210-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-13614-nc-3613-4210-36.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-13746-nc-3686-4220-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-13776-nc-3694-4222-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-13824-nc-3622-4289-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-13824-nc-3722-4252-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-13824-nc-3809-4269-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-13974-nc-3848-4325-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-13980-nc-3852-4329-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-14208-nc-3906-4373-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-14208-nc-3918-4388-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-14208-nc-3935-4385-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-14208-nc-3954-4382-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-14208-nc-4013-4415-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-14358-nc-4144-4452-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-14460-nc-4227-4461-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-14532-nc-4164-4470-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-14592-nc-4053-4513-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-14592-nc-4111-4501-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-14592-nc-4130-4492-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-14592-nc-4143-4519-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-14592-nc-4145-4528-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-14592-nc-4287-4557-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-14694-nc-4259-4564-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-14934-nc-4214-4587-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-14952-nc-4281-4592-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-14970-nc-4314-4597-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-14976-nc-4332-4626-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-14976-nc-4343-4629-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-14976-nc-4362-4637-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-14976-nc-4386-4610-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-14976-nc-4421-4638-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-14976-nc-4512-4650-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-15360-nc-4616-4709-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-15360-nc-4682-4771-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-15462-nc-4726-4801-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-15732-nc-4824-4841-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-15744-nc-4702-4862-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-15744-nc-4728-4853-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-15744-nc-4755-4857-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-15744-nc-4768-4907-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-15744-nc-4846-4894-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-16080-nc-4855-4950-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-16128-nc-4881-5010-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-16128-nc-4903-4978-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-16128-nc-4983-4973-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-16128-nc-4993-4991-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-16512-nc-5057-5067-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-16512-nc-5059-5068-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-16512-nc-5075-5085-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-16512-nc-5081-5111-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-16512-nc-5155-5092-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-16512-nc-5185-5089-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-16668-nc-5148-5159-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-16740-nc-5287-5167-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-16824-nc-5239-5183-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-16896-nc-5251-5270-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-16896-nc-5259-5228-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-16896-nc-5312-5213-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-17028-nc-5387-5293-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-17082-nc-5343-5299-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-17124-nc-5293-5305-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-17280-nc-5472-5349-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-17280-nc-5602-5411-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-17280-nc-5616-5391-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-17280-nc-5648-5395-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-17280-nc-5654-5402-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-17400-nc-5727-5423-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-17478-nc-5630-5428-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-17664-nc-5710-5497-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-17664-nc-5733-5514-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-17664-nc-5737-5485-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-17664-nc-5743-5509-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-17664-nc-5759-5488-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-17664-nc-5763-5461-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-17748-nc-5811-5528-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-17856-nc-5813-5541-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-18000-nc-5766-5553-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-18048-nc-5922-5600-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-18048-nc-5940-5601-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-18054-nc-5959-5635-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-18090-nc-6040-5640-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-18138-nc-6060-5654-1.hdf5";
  // // ======= End RESTART

  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-27.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-28.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-29.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-30.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-31-0-31.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-27.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-28.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-29.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-30.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-31.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-27.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-28.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-29.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-30.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-109-3-31.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-27.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-28.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-29.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-30.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-112-2-31.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-404-4-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-404-4-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-404-4-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-464-7-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-464-7-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-464-7-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-464-7-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-464-7-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-464-7-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-552-8-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-552-8-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-552-8-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-580-6-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-580-6-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-580-6-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-638-10-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-638-10-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-638-10-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-639-13-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-639-13-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-639-13-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-640-14-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-640-14-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-640-14-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-641-19-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-641-19-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-641-19-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-642-72-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-642-72-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-642-72-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-644-12-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-644-12-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-644-12-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-645-11-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-645-11-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-645-11-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-660-9-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-660-9-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-660-9-3.hdf5";

  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-2-2505-1.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-3-536-1.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-4-724-1.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-5-2299-1.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-5-2299-2.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-5-2299-3.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-11-4036-1.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-11-4036-2.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-11-4036-3.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-11-4036-4.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-11-4036-5.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-13-5130-1.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-13-5130-2.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-13-5130-3.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-13-5130-4.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-13-5130-5.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-14-5611-1.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-14-5611-2.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-14-5611-3.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-14-5611-4.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-14-5611-5.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-14-5611-6.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-14-5611-7.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-16-4589-1.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-16-4589-2.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-16-4589-3.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-16-4589-4.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-17-4723-1.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-17-4723-2.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-17-4723-3.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-17-4723-4.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-20-6978-1.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-20-6978-2.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-20-6978-3.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-20-6978-4.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-20-6978-5.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-20-6978-6.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-20-6978-7.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-20-6978-8.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-20-6978-9.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-20-6978-10.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-20-6978-11.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-22-7123-1.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-22-7123-2.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-22-7123-3.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-22-7123-4.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-22-7123-5.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-22-7123-6.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-22-7123-7.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-22-7123-8.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-22-7123-9.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-22-7123-10.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-22-7123-11.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-22-7123-12.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-22-7123-13.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-23-7519-1.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-23-7519-2.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-23-7519-3.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-23-7519-4.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-23-7519-5.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-23-7519-6.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-23-7519-7.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-23-7519-8.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-23-7519-9.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-23-7519-10.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-23-7519-11.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-24-7206-1.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-24-7206-2.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-24-7206-3.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-24-7206-4.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-24-7206-5.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-24-7206-6.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-24-7206-7.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-24-7206-8.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-24-7206-9.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-24-7206-10.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-24-7206-11.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-24-7206-12.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-24-7206-13.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-24-7206-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-26-7743-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-26-7743-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-26-7743-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-26-7743-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-26-7743-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-26-7743-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-26-7743-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-26-7743-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-26-7743-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-26-7743-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-26-7743-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-26-7743-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-26-7743-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-29-8238-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-29-8238-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-29-8238-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-29-8238-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-29-8238-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-29-8238-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-29-8238-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-29-8238-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-29-8238-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-29-8238-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-29-8238-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-29-8238-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-29-8238-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-29-8238-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-29-8238-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-29-8238-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-29-8238-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-30-9952-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-30-9952-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-30-9952-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-30-9952-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-30-9952-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-30-9952-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-30-9952-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-30-9952-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-30-9952-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-30-9952-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-30-9952-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-30-9952-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-30-9952-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-30-9952-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-30-9952-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-30-9952-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-31-8401-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-31-8401-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-31-8401-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-31-8401-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-31-8401-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-31-8401-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-31-8401-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-31-8401-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-31-8401-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-31-8401-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-31-8401-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-31-8401-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-31-8401-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-31-8401-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-31-8401-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-31-8401-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-31-8401-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-33-8433-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-33-8433-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-33-8433-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-33-8433-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-33-8433-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-33-8433-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-33-8433-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-33-8433-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-33-8433-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-33-8433-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-33-8433-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-33-8433-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-33-8433-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-33-8433-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-33-8433-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-33-8433-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-33-8433-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-34-9800-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-34-9800-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-34-9800-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-34-9800-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-34-9800-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-34-9800-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-34-9800-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-34-9800-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-34-9800-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-34-9800-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-34-9800-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-34-9800-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-34-9800-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-34-9800-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-34-9800-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-34-9800-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-34-9800-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-34-9800-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-34-9800-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-34-9800-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-35-10048-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-35-10048-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-35-10048-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-35-10048-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-35-10048-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-35-10048-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-35-10048-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-35-10048-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-35-10048-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-35-10048-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-35-10048-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-35-10048-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-35-10048-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-35-10048-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-35-10048-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-35-10048-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-35-10048-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-35-10048-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-37-10046-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-37-10046-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-37-10046-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-37-10046-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-37-10046-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-37-10046-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-37-10046-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-37-10046-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-37-10046-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-37-10046-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-37-10046-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-37-10046-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-37-10046-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-37-10046-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-37-10046-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-37-10046-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-37-10046-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-37-10046-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-37-10046-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-38-9816-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-38-9816-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-38-9816-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-38-9816-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-38-9816-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-38-9816-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-38-9816-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-38-9816-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-38-9816-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-38-9816-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-38-9816-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-38-9816-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-38-9816-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-38-9816-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-38-9816-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-38-9816-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-38-9816-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-38-9816-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-38-9816-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-38-9816-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-42-10505-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-42-10505-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-42-10505-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-42-10505-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-42-10505-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-42-10505-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-42-10505-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-42-10505-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-42-10505-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-42-10505-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-42-10505-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-42-10505-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-42-10505-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-42-10505-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-42-10505-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-42-10505-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-42-10505-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-42-10505-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-42-10505-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-42-10505-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-42-10505-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-42-10505-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-43-10319-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-43-10319-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-43-10319-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-43-10319-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-43-10319-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-43-10319-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-43-10319-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-43-10319-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-43-10319-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-43-10319-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-43-10319-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-43-10319-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-43-10319-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-43-10319-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-43-10319-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-43-10319-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-43-10319-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-43-10319-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-45-10749-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-45-10749-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-45-10749-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-45-10749-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-45-10749-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-45-10749-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-45-10749-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-45-10749-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-45-10749-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-45-10749-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-45-10749-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-45-10749-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-45-10749-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-45-10749-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-45-10749-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-45-10749-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-45-10749-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-45-10749-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-45-10749-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-45-10749-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-45-10749-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-46-10805-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-46-10805-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-46-10805-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-46-10805-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-46-10805-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-46-10805-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-46-10805-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-46-10805-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-46-10805-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-46-10805-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-46-10805-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-46-10805-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-46-10805-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-46-10805-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-46-10805-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-46-10805-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-46-10805-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-46-10805-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-46-10805-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-46-10805-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-46-10805-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-47-10807-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-47-10807-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-47-10807-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-47-10807-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-47-10807-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-47-10807-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-47-10807-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-47-10807-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-47-10807-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-47-10807-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-47-10807-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-47-10807-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-47-10807-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-47-10807-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-47-10807-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-47-10807-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-47-10807-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-47-10807-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-47-10807-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-47-10807-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-47-10807-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-47-10807-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-47-10807-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-47-10807-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-50-10987-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-50-10987-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-50-10987-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-50-10987-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-50-10987-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-50-10987-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-50-10987-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-50-10987-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-50-10987-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-50-10987-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-50-10987-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-50-10987-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-50-10987-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-50-10987-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-50-10987-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-50-10987-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-50-10987-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-50-10987-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-50-10987-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-50-10987-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-50-10987-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-51-10992-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-51-10992-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-51-10992-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-51-10992-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-51-10992-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-51-10992-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-51-10992-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-51-10992-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-51-10992-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-51-10992-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-51-10992-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-51-10992-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-51-10992-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-51-10992-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-51-10992-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-51-10992-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-51-10992-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-51-10992-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-51-10992-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-51-10992-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-51-10992-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-53-11042-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-53-11042-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-53-11042-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-53-11042-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-53-11042-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-53-11042-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-53-11042-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-53-11042-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-53-11042-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-53-11042-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-53-11042-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-53-11042-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-53-11042-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-53-11042-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-53-11042-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-53-11042-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-53-11042-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-53-11042-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-53-11042-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-53-11042-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-53-11042-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-53-11042-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-53-11042-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-54-11011-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-54-11011-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-54-11011-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-54-11011-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-54-11011-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-54-11011-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-54-11011-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-54-11011-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-54-11011-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-54-11011-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-54-11011-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-54-11011-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-54-11011-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-54-11011-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-54-11011-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-54-11011-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-54-11011-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-54-11011-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-54-11011-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-54-11011-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-54-11011-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-54-11011-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-54-11011-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-58-11165-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-58-11165-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-58-11165-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-58-11165-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-58-11165-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-58-11165-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-58-11165-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-58-11165-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-58-11165-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-58-11165-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-58-11165-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-58-11165-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-58-11165-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-58-11165-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-58-11165-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-58-11165-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-58-11165-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-58-11165-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-58-11165-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-58-11165-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-58-11165-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-58-11165-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-58-11165-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-58-11165-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-59-11576-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-59-11576-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-59-11576-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-59-11576-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-59-11576-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-59-11576-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-59-11576-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-59-11576-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-59-11576-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-59-11576-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-59-11576-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-59-11576-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-59-11576-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-59-11576-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-59-11576-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-59-11576-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-59-11576-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-59-11576-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-59-11576-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-59-11576-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-59-11576-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-59-11576-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-59-11576-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-59-11576-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-59-11576-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-60-11224-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-60-11224-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-60-11224-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-60-11224-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-60-11224-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-60-11224-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-60-11224-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-60-11224-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-60-11224-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-60-11224-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-60-11224-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-60-11224-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-60-11224-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-60-11224-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-60-11224-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-60-11224-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-60-11224-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-60-11224-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-60-11224-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-60-11224-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-60-11224-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-61-12252-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-27.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-28.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-62-11959-29.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-65-11333-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-65-11333-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-65-11333-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-65-11333-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-65-11333-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-65-11333-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-65-11333-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-65-11333-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-65-11333-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-65-11333-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-65-11333-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-65-11333-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-65-11333-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-65-11333-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-65-11333-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-65-11333-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-65-11333-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-65-11333-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-65-11333-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-65-11333-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-65-11333-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-65-11333-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-65-11333-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-65-11333-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-67-12960-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-67-12960-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-67-12960-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-67-12960-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-67-12960-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-67-12960-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-67-12960-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-67-12960-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-67-12960-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-67-12960-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-67-12960-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-67-12960-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-67-12960-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-67-12960-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-67-12960-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-67-12960-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-67-12960-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-67-12960-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-67-12960-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-67-12960-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-67-12960-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-67-12960-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-67-12960-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-67-12960-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-67-12960-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-69-12431-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-69-12431-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-69-12431-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-69-12431-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-69-12431-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-69-12431-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-69-12431-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-69-12431-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-69-12431-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-69-12431-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-69-12431-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-69-12431-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-69-12431-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-69-12431-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-69-12431-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-69-12431-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-69-12431-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-69-12431-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-69-12431-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-69-12431-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-69-12431-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-69-12431-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-69-12431-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-69-12431-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-69-12431-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-27.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-28.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-29.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-30.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-71-12680-31.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-27.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-28.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-29.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-73-13715-30.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-27.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-28.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-74-13707-29.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-27.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-75-13586-28.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-78-13483-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-78-13483-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-78-13483-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-78-13483-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-78-13483-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-78-13483-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-78-13483-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-78-13483-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-78-13483-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-78-13483-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-78-13483-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-78-13483-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-78-13483-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-78-13483-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-78-13483-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-78-13483-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-78-13483-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-78-13483-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-78-13483-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-78-13483-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-78-13483-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-78-13483-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-27.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-80-13581-28.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-27.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-28.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-29.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-30.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-82-14045-31.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-27.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-28.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-29.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-30.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-31.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-32.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-33.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-34.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-35.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-36.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-83-14290-37.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-27.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-28.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-29.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-30.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-84-13999-31.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-27.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-28.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-29.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-30.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-31.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-32.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-33.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-34.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-35.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-36.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-37.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-38.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-39.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-86-14336-40.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-27.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-28.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-29.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-30.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-31.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-32.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-33.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-34.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-35.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-36.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-37.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-88-14521-38.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-27.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-28.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-29.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-89-13818-30.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-27.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-28.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-29.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-30.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-31.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-32.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-33.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-34.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-91-14199-35.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-27.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-28.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-29.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-30.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-31.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-32.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-33.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-92-14159-34.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-27.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-28.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-29.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-30.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-31.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-32.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-33.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-93-14530-34.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-27.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-28.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-29.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-30.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-31.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-100-15174-32.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-27.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-28.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-29.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-30.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-31.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-32.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-33.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-34.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-35.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-36.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-37.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-38.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-39.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-101-14399-40.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-27.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-28.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-29.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-30.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-31.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-32.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-33.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-34.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-35.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-36.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-102-15592-37.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-27.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-28.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-29.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-104-14772-30.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-4.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-7.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-8.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-9.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-11.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-12.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-13.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-14.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-16.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-17.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-18.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-19.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-20.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-23.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-24.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-25.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-26.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-27.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-28.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-29.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-30.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-31.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-32.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-33.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-34.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-35.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-36.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-37.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-106-15594-38.hdf5";

  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-2-13-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-5-25-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-5-25-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-6-30-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-6-30-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-6-30-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-35-147-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-35-147-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-35-147-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-35-147-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-35-147-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-35-147-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-35-147-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-35-147-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-35-147-9.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-35-147-10.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-35-147-11.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-35-147-12.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-36-93-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-36-93-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-36-93-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-36-93-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-36-93-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-36-93-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-36-93-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-36-93-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-36-93-9.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-36-93-10.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-36-93-11.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-36-93-12.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-45-308-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-45-308-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-45-308-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-45-308-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-45-308-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-45-308-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-45-308-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-45-308-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-45-308-9.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-45-308-10.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-47-316-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-47-316-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-47-316-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-47-316-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-47-316-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-47-316-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-47-316-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-47-316-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-47-316-9.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-48-314-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-48-314-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-48-314-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-48-314-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-48-314-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-48-314-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-48-314-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-48-314-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-48-314-9.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-59-566-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-59-566-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-59-566-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-59-566-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-59-566-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-59-566-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-64-321-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-64-321-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-64-321-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-64-321-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-64-321-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-64-321-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-68-584-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-68-584-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-68-584-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-68-584-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-68-584-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-70-357-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-70-357-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-71-774-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-71-774-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-71-774-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-71-774-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-71-774-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-71-774-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-74-547-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-74-547-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-74-547-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-74-547-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-75-516-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-81-700-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-81-700-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-81-700-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-81-700-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-83-665-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-83-665-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-83-665-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-83-665-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-85-858-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-85-858-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-85-858-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-85-858-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-85-858-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-92-876-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-92-876-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-92-876-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-93-724-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-93-724-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-93-724-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-96-920-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-96-920-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-96-920-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-96-920-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-96-920-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-96-920-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-96-920-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-99-980-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-99-980-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-99-980-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-99-980-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-99-980-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-99-980-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-99-980-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-99-980-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-100-739-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-100-739-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-100-739-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-100-739-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-100-739-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-100-739-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-100-739-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-102-996-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-102-996-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-102-996-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-102-996-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-103-729-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-103-729-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-103-729-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-103-729-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-108-995-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-108-995-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-108-995-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-118-1007-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-118-1007-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-118-1007-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-119-1030-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-119-1030-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-119-1030-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-133-1202-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-133-1202-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-134-1087-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-134-1087-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-134-1087-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-136-1085-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-136-1085-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-136-1085-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-137-1201-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-137-1201-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-138-1199-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-138-1199-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-139-1176-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-139-1176-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-139-1176-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-139-1176-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-139-1176-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-146-1207-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-146-1207-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-146-1207-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-150-1170-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-150-1170-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-150-1170-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-152-1220-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-152-1220-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-152-1220-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-155-1259-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-155-1259-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-155-1259-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-155-1259-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-158-1134-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-160-1256-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-160-1256-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-165-1135-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-166-1131-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-166-1131-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-166-1131-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-177-1106-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-177-1106-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-177-1106-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-179-1262-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-179-1262-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-179-1262-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-179-1262-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-179-1262-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-179-1262-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-179-1262-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-190-1340-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-190-1340-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-190-1340-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-190-1340-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-190-1340-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-198-1476-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-198-1476-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-198-1476-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-198-1476-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-198-1476-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-200-1318-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-200-1318-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-204-1544-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-204-1544-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-204-1544-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-204-1544-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-205-1319-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-205-1319-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-208-1508-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-208-1508-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-208-1508-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-208-1508-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-216-1364-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-216-1364-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-217-1366-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-217-1366-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-217-1366-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-224-1363-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-224-1363-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-227-1290-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-227-1290-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-239-1447-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-239-1447-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-250-1285-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-250-1285-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-250-1285-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-253-1649-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-253-1649-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-254-1568-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-254-1568-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-254-1568-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-254-1568-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-257-1642-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-262-1439-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-262-1439-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-267-1421-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-267-1421-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-267-1421-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-285-1606-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-287-1604-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-293-1917-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-293-1917-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-293-1917-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-298-1578-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-298-1578-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-298-1578-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-302-1920-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-302-1920-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-302-1920-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-304-1955-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-304-1955-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-304-1955-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-307-1947-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-307-1947-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-307-1947-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-309-1909-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-309-1909-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-309-1909-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-310-1701-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-310-1701-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-310-1701-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-318-1772-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-318-1772-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-321-1898-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-321-1898-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-321-1898-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-321-1898-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-322-1884-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-322-1884-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-322-1884-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-322-1884-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-324-1768-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-329-1764-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-329-1764-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-334-1804-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-334-1804-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-334-1804-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-338-1806-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-338-1806-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-338-1806-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-340-1709-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-340-1709-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-340-1709-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-346-1816-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-346-1816-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-346-1816-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-348-2272-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-348-2272-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-348-2272-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-348-2272-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-353-1779-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-353-1779-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-355-1878-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-355-1878-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-355-1878-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-356-1835-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-356-1835-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-356-1835-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-360-2294-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-360-2294-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-360-2294-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-360-2294-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-361-1745-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-361-1745-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-361-1745-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-363-2283-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-363-2283-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-363-2283-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-363-2283-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-363-2283-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-366-1965-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-366-1965-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-366-1965-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-366-1965-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-366-1965-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-366-1965-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-385-1841-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-385-1841-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-385-1841-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-387-1966-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-387-1966-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-387-1966-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-387-1966-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-387-1966-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-387-1966-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-388-2308-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-388-2308-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-388-2308-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-388-2308-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-388-2308-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-389-1851-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-389-1851-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-389-1851-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-391-2237-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-391-2237-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-392-2328-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-392-2328-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-392-2328-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-392-2328-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-396-1734-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-401-2358-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-401-2358-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-401-2358-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-406-2197-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-406-2197-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-407-2214-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-407-2214-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-413-2135-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-417-2130-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-417-2130-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-423-2377-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-451-2449-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-451-2449-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-452-2070-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-452-2070-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-453-2075-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-453-2075-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-455-2063-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-466-2477-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-466-2477-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-468-2029-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-468-2029-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-469-2043-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-472-2053-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-472-2053-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-472-2053-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-477-2469-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-477-2469-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-482-2537-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-483-2561-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-487-2476-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-487-2476-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-488-2485-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-488-2485-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-489-2484-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-489-2484-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-490-2544-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-492-2575-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-494-2533-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-497-2525-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-497-2525-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-505-2611-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-505-2611-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-505-2611-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-506-1976-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-506-1976-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-506-1976-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-511-2596-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-512-2590-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-514-2597-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-516-2609-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-516-2609-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-530-2631-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-532-2616-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-544-2626-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-551-2628-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-561-2644-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-562-2625-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-563-1993-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-577-1991-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-579-2661-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-580-2692-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-581-3735-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-582-2698-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-584-2665-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-585-2694-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-589-2754-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-589-2754-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-593-2756-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-593-2756-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-600-2703-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-604-11438-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-604-11438-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-604-11438-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-604-11438-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-604-11438-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-604-11438-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-604-11438-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-604-11438-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-604-11438-9.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-610-3876-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-611-11429-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-611-11429-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-611-11429-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-611-11429-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-611-11429-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-611-11429-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-611-11429-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-611-11429-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-611-11429-9.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-612-3736-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-617-3796-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-620-3797-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-632-11394-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-632-11394-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-632-11394-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-632-11394-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-632-11394-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-632-11394-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-632-11394-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-633-3770-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-634-3807-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-641-11340-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-641-11340-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-641-11340-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-641-11340-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-641-11340-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-641-11340-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-641-11340-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-641-11340-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-644-9914-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-644-9914-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-644-9914-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-644-9914-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-644-9914-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-644-9914-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-645-2858-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-646-9824-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-646-9824-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-646-9824-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-646-9824-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-646-9824-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-646-9824-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-651-3596-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-652-3789-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-656-2866-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-660-3841-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-667-11516-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-667-11516-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-667-11516-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-667-11516-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-667-11516-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-667-11516-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-667-11516-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-671-2887-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-674-9764-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-674-9764-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-674-9764-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-674-9764-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-674-9764-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-674-9764-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-674-9764-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-674-9764-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-678-3088-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-693-3964-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-698-2920-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-702-2932-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-712-4000-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-715-11607-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-715-11607-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-715-11607-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-715-11607-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-715-11607-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-715-11607-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-715-11607-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-721-3188-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-724-3264-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-733-3431-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-736-4524-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-744-5584-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-753-4678-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-754-4930-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-755-5712-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-758-9340-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-758-9340-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-758-9340-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-758-9340-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-758-9340-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-758-9340-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-758-9340-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-758-9340-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-766-8284-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-766-8284-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-766-8284-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-768-6530-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-769-9320-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-769-9320-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-769-9320-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-769-9320-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-769-9320-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-769-9320-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-769-9320-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-771-6558-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-800-10583-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-800-10583-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-800-10583-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-800-10583-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-800-10583-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-800-10583-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-800-10583-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-800-10583-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-803-16221-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-803-16221-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-803-16221-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-803-16221-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-803-16221-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-803-16221-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-803-16221-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-803-16221-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-803-16221-9.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-804-17375-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-804-17375-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-804-17375-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-804-17375-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-804-17375-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-804-17375-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-804-17375-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-804-17375-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-804-17375-9.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-805-15948-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-805-15948-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-805-15948-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-805-15948-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-805-15948-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-805-15948-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-805-15948-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-805-15948-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-805-15948-9.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-807-15857-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-807-15857-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-807-15857-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-807-15857-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-807-15857-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-807-15857-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-807-15857-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-807-15857-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-807-15857-9.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-809-18144-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-809-18144-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-809-18144-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-809-18144-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-809-18144-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-809-18144-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-809-18144-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-809-18144-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-809-18144-9.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-810-16130-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-810-16130-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-810-16130-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-810-16130-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-810-16130-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-810-16130-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-810-16130-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-810-16130-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-810-16130-9.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-811-10878-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-811-10878-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-811-10878-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-811-10878-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-811-10878-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-811-10878-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-811-10878-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-811-10878-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-812-16086-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-812-16086-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-812-16086-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-812-16086-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-812-16086-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-812-16086-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-812-16086-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-812-16086-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-812-16086-9.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-818-17116-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-818-17116-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-818-17116-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-818-17116-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-818-17116-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-818-17116-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-818-17116-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-818-17116-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-818-17116-9.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-821-11143-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-821-11143-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-821-11143-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-821-11143-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-821-11143-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-821-11143-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-821-11143-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-821-11143-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-823-17059-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-823-17059-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-823-17059-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-823-17059-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-823-17059-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-823-17059-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-823-17059-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-823-17059-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-823-17059-9.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-824-17047-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-824-17047-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-824-17047-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-824-17047-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-824-17047-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-824-17047-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-824-17047-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-824-17047-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-824-17047-9.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-826-11242-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-826-11242-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-826-11242-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-826-11242-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-826-11242-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-826-11242-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-826-11242-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-826-11242-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-829-11323-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-829-11323-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-829-11323-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-829-11323-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-829-11323-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-829-11323-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-829-11323-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-829-11323-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-837-11278-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-837-11278-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-837-11278-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-837-11278-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-837-11278-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-837-11278-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-837-11278-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-837-11278-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-839-11280-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-839-11280-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-839-11280-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-839-11280-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-839-11280-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-839-11280-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-839-11280-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-839-11280-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-841-11259-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-841-11259-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-841-11259-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-841-11259-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-841-11259-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-841-11259-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-841-11259-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-841-11259-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-845-11306-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-845-11306-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-845-11306-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-845-11306-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-845-11306-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-845-11306-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-845-11306-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-845-11306-8.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-851-11250-1.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-851-11250-2.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-851-11250-3.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-851-11250-4.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-851-11250-5.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-851-11250-6.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-851-11250-7.hdf5";
  //     data_collection_1[n_data++] =
  //     "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-851-11250-8.hdf5";

  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-73-1-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/BoxStacks/BoxStacks-ndof-450-nc-464-7-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-2-2505-1.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-3-536-1.hdf5";
  //   data_collection_1[n_data++] =
  //   "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-4-724-1.hdf5";

  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-259-21436-15.hdf5"; // CP:
  // 0.5,0.8,33 // CP+PRJ: 0.5, 0.8, 60

  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-262-1439-2.hdf5"; // CP:
  // 0.5, 0.8, 99  // CP+PRJ: 0.5, 0.8, 133 data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-505-2611-1.hdf5"; // CP:
  // 0.5, 0.8, 69 data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-563-1993-1.hdf5"; // CP:
  // 0.5, 0.7, 802   // CP+PRJ: 0.5, 0.9, 1417 or 0.5, 0.7, 857 or 0.49, 0.9, 265
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-633-3770-1.hdf5";

  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1353-243-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2009-420-1.hdf5"; //
  // CP: 0.5, 0.8, 398 data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2612-641-1.hdf5"; //
  // CP: 0.5, 0.8, 52 data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3190-2112-1.hdf5"; //
  // CP: 0.5, 0.8, 41 data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3179-2042-1.hdf5"; //
  // CP: 0.5, 0.8, 47

  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4250-15000-1.hdf5"; // CP: 0.5,
  // 0.7, 53 data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4371-17786-1.hdf5"; // CP: 0.5,
  // 0.6, 56

  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-2730-6-1.hdf5";
  // // CP: 0.5, 0.7, 72 data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-2941-7-2.hdf5";
  // // CP: 0.5, 0.8, 56 data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4432-21-1.hdf5";
  // // CP: 0.5, 0.6, 96 data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4489-114-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4520-251-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4536-377-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4557-720-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4573-845-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4575-901-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4603-1534-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4610-1948-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4613-1698-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4615-1898-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4625-1864-1.hdf5";

  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-12288-nc-3130-3742-1.hdf5";

  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-155-18120-39.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-355-36003-34.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-172-18659-21.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-224-20891-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-130-16647-10.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-259-21436-15.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-156-17670-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-325-32051-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-284-23843-3.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-325-32051-22.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-179-19746-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-244-20980-6.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-197-19889-5.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-280-23762-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Capsules/Capsules-ndof-600-nc-243-21083-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-12288-nc-3130-3742-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-16080-nc-4855-4950-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Chute/Chute-ndof-17856-nc-5813-5541-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-262-1439-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-392-2328-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-482-2537-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-489-2484-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-505-2611-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-532-2616-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-563-1993-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-579-2661-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-584-2665-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-633-3770-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-651-3596-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-674-9764-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-671-2887-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-712-4000-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/KaplasTower/KaplasTower-ndof-864-nc-766-8284-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1353-243-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1611-330-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1768-366-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2009-420-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2250-473-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2256-472-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2282-500-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2438-544-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2470-580-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2517-590-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2612-641-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2633-643-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2708-701-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2724-713-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2925-1012-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2920-1003-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2923-1005-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-2926-1006-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3047-1333-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3141-1907-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3190-2112-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3194-2081-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-3179-2042-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4237-15402-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4250-15000-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres/Spheres-ndof-12000-nc-4371-17786-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-2407-3-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-2730-6-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-2941-7-2.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-3389-10-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-3940-12-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4025-13-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4144-14-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4369-17-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4383-18-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4432-21-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4451-85-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4464-23-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4481-105-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4484-70-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4489-114-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4488-135-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4500-148-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4515-182-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4520-251-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4536-377-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4547-441-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4559-682-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4557-720-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4564-791-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4573-845-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4575-901-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4581-866-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4603-1534-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4610-1948-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4613-1698-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4615-1898-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4625-1864-1.hdf5";
  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/Spheres1mmScaled/Spheres1mmScaled-ndof-12000-nc-4627-1845-1.hdf5";

  // data_collection_1[n_data++] =
  // "./fclib-sub/Global/siconos/PrimitiveSoup/PrimitiveSoup-ndof-6000-nc-1087-183-68.hdf5";

  data_collection_1[n_data++] = "---";
  return data_collection_1;
}
