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

  // List of rolling problems
  int listprob[1] = {1}; // There is only 1 type now
  /* 0: RollingSpheres      #28  problems */
  /*   TOTAL          #1169 problems */

  if (listprob[0]==1)
  {
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-18-nc-2-0.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-1008-nc-337-7427.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-102-nc-17-205.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-1020-nc-341-7540.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-1026-nc-340-7579.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-1026-nc-343-7561.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-1026-nc-345-7595.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-1032-nc-342-7654.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-1038-nc-344-7690.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-1038-nc-346-7694.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-1044-nc-350-7753.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-1050-nc-348-7830.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-1050-nc-356-7786.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-1050-nc-357-7778.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-1050-nc-360-7775.hdf5";
    data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-1050-nc-361-7776.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-1056-nc-351-7839.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-108-nc-21-214.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-114-nc-22-219.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-114-nc-23-223.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-126-nc-25-247.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-138-nc-29-319.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-138-nc-33-309.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-144-nc-35-348.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-150-nc-34-365.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-150-nc-36-376.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-156-nc-37-398.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-162-nc-40-460.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-174-nc-39-517.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-186-nc-43-582.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-192-nc-45-585.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-192-nc-46-611.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-198-nc-48-634.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-198-nc-51-640.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-198-nc-52-649.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-210-nc-47-697.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-222-nc-55-731.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-222-nc-57-774.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-228-nc-59-802.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-228-nc-60-799.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-228-nc-61-804.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-228-nc-62-821.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-240-nc-63-827.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-240-nc-64-826.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-246-nc-65-846.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-252-nc-68-876.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-252-nc-70-869.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-282-nc-75-1057.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-288-nc-76-1076.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-294-nc-81-1127.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-294-nc-83-1139.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-306-nc-85-1220.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-312-nc-87-1230.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-312-nc-90-1232.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-318-nc-91-1309.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-324-nc-86-1336.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-330-nc-89-1412.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-330-nc-94-1421.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-336-nc-93-1436.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-36-nc-4-28.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-360-nc-104-1608.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-360-nc-99-1631.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-366-nc-101-1655.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-372-nc-100-1686.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-372-nc-97-1705.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-384-nc-106-1805.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-390-nc-108-1846.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-390-nc-109-1856.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-390-nc-111-1870.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-396-nc-113-1914.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-396-nc-115-1958.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-408-nc-118-2027.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-408-nc-120-2051.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-408-nc-124-2058.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-420-nc-117-2118.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-420-nc-119-2121.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-420-nc-123-2128.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-426-nc-127-2217.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-432-nc-132-2267.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-444-nc-129-2386.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-444-nc-130-2332.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-450-nc-135-2433.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-456-nc-139-2467.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-462-nc-136-2536.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-462-nc-140-2526.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-462-nc-141-2500.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-468-nc-138-2553.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-468-nc-146-2583.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-474-nc-142-2615.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-474-nc-147-2657.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-48-nc-5-75.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-48-nc-6-96.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-486-nc-149-2766.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-492-nc-150-2784.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-492-nc-152-2781.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-492-nc-156-2822.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-492-nc-158-2816.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-498-nc-151-2835.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-528-nc-159-3149.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-534-nc-163-3207.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-540-nc-166-3274.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-540-nc-167-3252.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-546-nc-171-3300.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-552-nc-170-3371.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-552-nc-173-3349.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-552-nc-174-3362.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-558-nc-172-3426.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-558-nc-178-3437.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-558-nc-179-3446.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-558-nc-181-3438.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-576-nc-183-3597.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-588-nc-184-3675.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-588-nc-189-3697.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-594-nc-190-3745.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-594-nc-193-3751.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-60-nc-7-117.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-600-nc-195-3856.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-600-nc-199-3846.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-606-nc-194-3879.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-606-nc-197-3882.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-612-nc-198-3920.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-618-nc-206-4016.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-624-nc-212-4100.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-630-nc-210-4128.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-636-nc-201-4175.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-636-nc-207-4178.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-642-nc-211-4248.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-642-nc-213-4256.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-654-nc-215-4357.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-654-nc-217-4345.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-654-nc-218-4361.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-66-nc-9-121.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-660-nc-226-4393.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-672-nc-220-4466.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-678-nc-232-4555.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-684-nc-230-4604.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-690-nc-231-4643.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-696-nc-216-4718.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-702-nc-221-4789.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-702-nc-225-4763.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-726-nc-224-4963.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-732-nc-233-5014.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-738-nc-237-5076.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-738-nc-238-5072.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-744-nc-240-5177.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-750-nc-244-5216.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-756-nc-245-5274.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-756-nc-246-5266.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-762-nc-242-5326.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-762-nc-249-5297.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-762-nc-250-5331.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-762-nc-253-5333.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-774-nc-257-5426.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-78-nc-10-135.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-78-nc-12-147.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-780-nc-256-5491.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-780-nc-259-5509.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-780-nc-260-5470.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-786-nc-265-5544.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-798-nc-268-5597.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-798-nc-274-5582.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-816-nc-262-5775.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-822-nc-258-5782.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-822-nc-264-5787.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-822-nc-272-5832.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-828-nc-270-5851.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-828-nc-277-5856.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-834-nc-280-5930.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-834-nc-281-5957.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-840-nc-283-5962.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-852-nc-282-6072.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-858-nc-275-6106.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-864-nc-287-6161.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-870-nc-288-6260.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-876-nc-294-6270.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-876-nc-295-6299.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-882-nc-297-6328.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-888-nc-303-6428.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-894-nc-292-6440.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-894-nc-293-6470.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-894-nc-296-6445.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-90-nc-14-182.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-906-nc-302-6577.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-912-nc-299-6633.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-912-nc-300-6637.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-924-nc-310-6734.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-930-nc-311-6786.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-936-nc-315-6821.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-936-nc-316-6842.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-936-nc-317-6844.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-942-nc-314-6857.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-948-nc-319-6922.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-96-nc-18-186.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-960-nc-321-7029.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-960-nc-322-7030.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-960-nc-330-7010.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-966-nc-323-7052.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-972-nc-328-7120.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-990-nc-327-7308.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-990-nc-329-7310.hdf5";
    // data_collection_1[n_data++] = "/Users/minhnguyen/fclib-library/Global/siconos/Rolling/RollingSpheres/RollingSpheres-ndof-996-nc-333-7321.hdf5";
  }

  data_collection_1[n_data++] = "---";
  return data_collection_1;

}

