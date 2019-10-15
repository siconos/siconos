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
  *number_of_tests = 25;//n_data * n_solvers;
  TestCase * tests_list = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));


  // "External" solver parameters
  // -> same values for all tests.
  // The differences between tests are only for internal solvers and input data.
  int topsolver = SICONOS_FRICTION_3D_NSGS;
  int dpos[] = {1, SICONOS_DPARAM_TOL};  // ipos = [number of values in parameters list, indices]
  double dparam[] = {1e-16};
  int ipos[] = {1, SICONOS_IPARAM_MAX_ITER};  // ipos = [number of values in parameters list, indices]
  int iparam[] = {10000};

  int current = 0;

  
  {
    int d = 0; // FC3D_Example1_SBM.dat
    
    // Projection on cone, default values.
    dparam[SICONOS_DPARAM_TOL] = 1e-16;
    iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone, NULL, NULL, NULL, NULL,
               &tests_list[current++]);
  }

  {
    int d = 0; // FC3D_Example1_SBM.dat
    // Projection on cone with diagonalization, default value.
    dparam[SICONOS_DPARAM_TOL] = 1e-16;
    iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithDiagonalization, NULL, NULL, NULL, NULL,
               &tests_list[current++]);
  }
  
  {
    int d = 0; // FC3D_Example1_SBM.dat
    // Projection on cone with local iteration, set tol and max iter.
    dparam[SICONOS_DPARAM_TOL] = 1e-16;
    iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    double internal_dparam[] = {1e-3};
    int internal_iparam[] = {10};
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration,
               dpos, internal_dparam, ipos, internal_iparam,
               &tests_list[current++]);
  }

  {
    int d= 0; // FC3D_Example1_SBM.dat
    dparam[SICONOS_DPARAM_TOL] = 1e-16;
    iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    // Projection on cone with regularization, set rho
    int intern_dpos[] = {1, SICONOS_FRICTION_3D_NSN_RHO};
    double internal_dparam[] = {0.1}; // rho value
    
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization,
               intern_dpos, internal_dparam,  NULL, NULL,
               &tests_list[current++]);
  }

  {
    int d=1; // "./data/Capsules-i122-1617.dat"

    // Projection on cone with local iteration, set tol, itermax, d[9], i[8]
    dparam[SICONOS_DPARAM_TOL] = 1e-16;
    iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    int intern_dpos[] = {2, SICONOS_DPARAM_TOL, 9};
    double internal_dparam[] = {1e-16, 1.};
    int intern_ipos[] = {2, SICONOS_IPARAM_MAX_ITER, 8};
    int internal_iparam[] = {20, 1};
    
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration,
               intern_dpos, internal_dparam, intern_ipos, internal_iparam,
               &tests_list[current++]);
    // Expected to fail ...
    tests_list[current - 1].will_fail = 1;
  }
  
  {
    int d = 2; // Confeti-ex13-4contact-Fc3D-SBM.dat
    // Projection on cone set d[9], i[8]
    dparam[SICONOS_DPARAM_TOL] = 1e-5;
    iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    int intern_dpos[] = {1, 9};
    double internal_dparam[] = { 1.};
    int intern_ipos[] = {1, 8};
    int internal_iparam[] = { 1};
    
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone,
               intern_dpos, internal_dparam, intern_ipos, internal_iparam,
               &tests_list[current++]);
    // Expected to fail ...
    tests_list[current - 1].will_fail = 1;
  }

  {
    int d = 2; // Confeti-ex13-4contact-Fc3D-SBM.dat
    // nonsmooth newton. Set tol, max iter and i[1]. Default for other parameters
    
    dparam[SICONOS_DPARAM_TOL] = 1e-12;
    iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    double internal_dparam[] = {1e-18};
    int intern_ipos[] = {2, SICONOS_IPARAM_MAX_ITER, 1};  // current iteration number ???
    int internal_iparam[] = {10, 1};
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_NSN,
               dpos, internal_dparam, intern_ipos, internal_iparam,
               &tests_list[current++]);
  }

  
  {
    int d = 2; // Confeti-ex13-4contact-Fc3D-SBM.dat
    // Projection on cone with local iteration, set tol and maxiter
    dparam[SICONOS_DPARAM_TOL] = 1e-12;
    iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    double internal_dparam[] = { 1e-6};
    int internal_iparam[] = { 100};
    
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration,
               dpos, internal_dparam, ipos, internal_iparam,
               &tests_list[current++]);
  }

  
  {
    int d = 2; // Confeti-ex13-4contact-Fc3D-SBM.dat
    // Projection on cone with regularization, default values.
    dparam[SICONOS_DPARAM_TOL] = 1e-12;
    iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization,
               NULL, NULL, NULL, NULL,
               &tests_list[current++]);
  }
  
  {
    int d = 2; // Confeti-ex13-4contact-Fc3D-SBM.dat
    // Projection on cone, default values.
    dparam[SICONOS_DPARAM_TOL] = 1e-2;
    iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone,
               NULL, NULL, NULL, NULL,
               &tests_list[current++]);
  }

  {
    int d = 2; // Confeti-ex13-4contact-Fc3D-SBM.dat
    // nonsmooth newton. Set tol and i[1]. Default for other parameters
    dparam[SICONOS_DPARAM_TOL] = 1e-5;
    iparam[SICONOS_IPARAM_MAX_ITER]= 1000;
    double internal_dparam[] = {1e-16};
    int internal_iparam[] = {10};
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_NSN,
               dpos, internal_dparam, ipos, internal_iparam,
               &tests_list[current++]);
    }

  {
    int d = 2; // Confeti-ex13-4contact-Fc3D-SBM.dat
    // Projection on cone with local iteration, set tol and maxiter
    dparam[SICONOS_DPARAM_TOL] = 1e-12;
    iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    double internal_dparam[] = { 1e-6};
    int internal_iparam[] = { 100};
    
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration,
               dpos, internal_dparam, ipos, internal_iparam,
               &tests_list[current++]);

    internal_dparam[SICONOS_DPARAM_TOL] = 1e-16;
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration,
               dpos, internal_dparam, ipos, internal_iparam,
               &tests_list[current++]);
  }

  {
    int d = 3; // GFC3D_TwoRods1-condensed.dat
    // nonsmooth newton. Set tol and i[1]. Default for other parameters
    dparam[SICONOS_DPARAM_TOL] = 1e-12;
    iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    double internal_dparam[] = {1e-18};
    int intern_ipos[] = {2, SICONOS_IPARAM_MAX_ITER, 1};  // current iteration number ???
    int internal_iparam[] = {10, 1};
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_NSN,
               dpos, internal_dparam, intern_ipos, internal_iparam,
               &tests_list[current++]);
  }
  

  {
    int d  = 4; // FC3D_Example1.dat
    
    // nonsmooth newton. Set tol and max iter. Default for other parameters
    dparam[SICONOS_DPARAM_TOL] = 1e-12;
    iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    double internal_dparam[] = {1e-18};
    int internal_iparam[] = {10};
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_NSN,
               dpos, internal_dparam, ipos, internal_iparam,
               &tests_list[current++]);
  }

  {
    int d = 5; // Confeti-ex03-Fc3D-SBM.dat
    // Projection on cone, default values.
    dparam[SICONOS_DPARAM_TOL] = 1e-5;
    iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone,
               NULL, NULL, NULL, NULL,
               &tests_list[current++]);
    // Expected to fail ...
    tests_list[current - 1].will_fail = 1;
  }

  {
    int d = 5;
    // nonsmooth newton. Set tol and max iter. Default for other parameters
    dparam[SICONOS_DPARAM_TOL] = 1e-5;
    iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    double internal_dparam[] = {1e-16};
    int internal_iparam[] = {10};
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_NSN,
               dpos, internal_dparam, ipos, internal_iparam,
               &tests_list[current++]);
    // Expected to fail ...
    tests_list[current - 1].will_fail = 1;
  }

  {
    int d = 5;
    // Projection on cone with local iteration, set tol and maxiter
    dparam[SICONOS_DPARAM_TOL] = 1e-5;
    iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    double internal_dparam[] = { 1e-12};
    int internal_iparam[] = { 10};
    
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration,
               dpos, internal_dparam, ipos, internal_iparam,
               &tests_list[current++]);
  }
 

  {
    int d = 5;
    // Projection on cone with regularization, set tol and maxiter
    dparam[SICONOS_DPARAM_TOL] = 1e-5;
    iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
    double internal_dparam[] = { 1e-8};
    int internal_iparam[] = { 10};
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization,
               dpos, internal_dparam, ipos, internal_iparam,
               &tests_list[current++]);
    // Expected to fail ...
    tests_list[current - 1].will_fail = 1;
  }


  
  {
    int d = 7;
    // nonsmooth newton 'damped'. Set tol and max iter. Default for other parameters
    dparam[SICONOS_DPARAM_TOL] = 1e-3;
    iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
    double internal_dparam[] = {1e-16};
    int internal_iparam[] = {100};
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_NSN_GP,
               dpos, internal_dparam, ipos, internal_iparam,
               &tests_list[current++]);
    // Expected to fail ...
    tests_list[current - 1].will_fail = 1;

    internal_iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_NSN_GP,
               dpos, internal_dparam, ipos, internal_iparam,
               &tests_list[current++]);
    // Expected to fail ...
    tests_list[current - 1].will_fail = 1;
  }

  {
    int d = 7;
    // nonsmooth newton. Set tol and max iter. Default for other parameters
    dparam[SICONOS_DPARAM_TOL] = 1e-3;
    iparam[SICONOS_IPARAM_MAX_ITER] = 2000;
    double internal_dparam[] = {1e-16};
    int internal_iparam[] = {100};
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_NSN,
               dpos, internal_dparam, ipos, internal_iparam,
               &tests_list[current++]);
    // Expected to fail ...
    tests_list[current - 1].will_fail = 1;

    internal_iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_NSN,
               dpos, internal_dparam, ipos, internal_iparam,
               &tests_list[current++]);
    // Expected to fail ...
    tests_list[current - 1].will_fail = 1;
  }

  {
    int d = 7;
    // Projection on cone with local iteration, set tol and maxiter
    dparam[SICONOS_DPARAM_TOL] = 1e-3;
    iparam[SICONOS_IPARAM_MAX_ITER] = 2000;
    double internal_dparam[] = { 1e-6};
    int internal_iparam[] = { 100};
    
    build_friction_test(data_collection[d],
               topsolver, dpos, dparam, ipos, iparam,
               SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration,
               dpos, internal_dparam, ipos, internal_iparam,
               &tests_list[current++]);
    // Expected to fail ...
    tests_list[current - 1].will_fail = 1;
  }

  *number_of_tests = current;
  return tests_list;

}
