/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

#include <stdio.h>                       // for NULL
#include <stdlib.h>                      // for malloc
#include "Friction_cst.h"                // for SICONOS_FRICTION_3D_NSN_HYBR...
#include "SolverOptions.h"               // for SICONOS_DPARAM_TOL, SICONOS_...
#include "frictionContact_test_utils.h"  // for build_friction_test, build_t...
#include "test_utils.h"                  // for TestCase

TestCase * build_test_collection(int n_data, const char ** data_collection, int* number_of_tests)
{
  *number_of_tests = 9;//n_data * n_solvers;
  TestCase * collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));


  // "External" solver parameters
  // -> same values for all tests.
  // The differences between tests are only for internal solvers and input data.
  int topsolver = SICONOS_FRICTION_3D_NSGS;
  int dpos[] = {1, SICONOS_DPARAM_TOL};  // ipos = [number of values in parameters list, indices]
  double dparam[] = {1e-5};
  int ipos[] = {1, SICONOS_IPARAM_MAX_ITER};  // ipos = [number of values in parameters list, indices]
  int iparam[] = {2000};

  int current = 0;

  

  {
    int d = 6; // BoxesStack1-i100000-32.hdf5.dat
    // nonsmooth newton 'damped', Moreau-Jean formulation. Default for other parameters
    int intern_ipos[] = {1, SICONOS_FRICTION_3D_NSN_FORMULATION};
    int internal_iparam[] = {SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD};
    dparam[SICONOS_DPARAM_TOL] = 1e-5;
    iparam[SICONOS_IPARAM_MAX_ITER] = 1500;
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_NSN_GP, NULL,  NULL, intern_ipos, internal_iparam,
               &collection[current++]);
    collection[current-1].will_fail=1;
  }
  {
    int d = 6; // BoxesStack1-i100000-32.hdf5.dat
    // nonsmooth newton 'damped', change hybrid strategy. Default for other parameters
    int intern_ipos[] = {1, SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY};
    int internal_iparam[] = {SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP};
    dparam[SICONOS_DPARAM_TOL] = 1e-5;
    iparam[SICONOS_IPARAM_MAX_ITER] = 1500;
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_NSN_GP, NULL,  NULL, intern_ipos, internal_iparam,
               &collection[current++]);
    // expected to fail
    collection[current-1].will_fail=1;

  }
  {
    int d = 8; // KaplasTower-i1061-4.hdf5.dat";
    
    // Nonsmooth Newton "damped", default values.
    dparam[SICONOS_DPARAM_TOL] = 1e-5;
    iparam[SICONOS_IPARAM_MAX_ITER] = 2000;
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_NSN_GP, NULL, NULL, NULL, NULL,
               &collection[current++]);
  }
  {
    int d = 8; // KaplasTower-i1061-4.hdf5.dat";
   
    // nonsmooth newton 'damped', Moreau-Jean formulation. Default for other parameters
    int intern_ipos[] = {1, SICONOS_FRICTION_3D_NSN_FORMULATION};
    int internal_iparam[] = {SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD};
    dparam[SICONOS_DPARAM_TOL] = 1e-5;
    iparam[SICONOS_IPARAM_MAX_ITER] = 1500;
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_NSN_GP, NULL,  NULL, intern_ipos, internal_iparam,
               &collection[current++]);

  }
  {
    int d = 8; // KaplasTower-i1061-4.hdf5.dat";
    // nonsmooth newton 'damped', change hybrid strategy. Default for other parameters
    int intern_ipos[] = {1, SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY};
    int internal_iparam[] = {SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP};
    dparam[SICONOS_DPARAM_TOL] = 1e-5;
    iparam[SICONOS_IPARAM_MAX_ITER] = 1500;
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_NSN_GP, NULL,  NULL, intern_ipos, internal_iparam,
               &collection[current++]);
  }


  {
    int d = 9; // OneObject-i100000-499.hdf5.dat
    // Nonsmooth Newton "damped", default values.
    dparam[SICONOS_DPARAM_TOL] = 1e-5;
    iparam[SICONOS_IPARAM_MAX_ITER] = 100000;
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_NSN_GP, NULL, NULL, NULL, NULL,
               &collection[current++]);
    // expected to fail
    collection[current-1].will_fail=1;
  }
  {
    int d = 9; // OneObject-i100000-499.hdf5.dat
    // nonsmooth newton 'damped', Moreau-Jean formulation. Default for other parameters
    int intern_ipos[] = {1, SICONOS_FRICTION_3D_NSN_FORMULATION};
    int internal_iparam[] = {SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD};
    dparam[SICONOS_DPARAM_TOL] = 1e-5;
    iparam[SICONOS_IPARAM_MAX_ITER] = 1500;
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_NSN_GP, NULL,  NULL, intern_ipos, internal_iparam,
               &collection[current++]);
    // expected to fail
    collection[current-1].will_fail=1;

  }  
  {
    int d = 9; // OneObject-i100000-499.hdf5.dat
    // nonsmooth newton 'damped', change hybrid strategy. Default for other parameters
    int intern_ipos[] = {1, SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY};
    int internal_iparam[] = {SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP};
    dparam[SICONOS_DPARAM_TOL] = 1e-5;
    iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_NSN_GP, NULL,  NULL, intern_ipos, internal_iparam,
               &collection[current++]);
    // expected to fail
    collection[current-1].will_fail=1;
  }

  {
    int d = 9; // OneObject-i100000-499.hdf5.dat
    // Projection on cone with local iteration, set tol and max iter.
    dparam[SICONOS_DPARAM_TOL] = 1e-5;
    iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    double internal_dparam[] = {1e-12};
    int internal_iparam[] = {10};
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration,
               dpos, internal_dparam, ipos, internal_iparam,
               &collection[current++]);
    // expected to fail
    collection[current-1].will_fail=1;
    
  }
  
  *number_of_tests = current;
  return collection;

}
