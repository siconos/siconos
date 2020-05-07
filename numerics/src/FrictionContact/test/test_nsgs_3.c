/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

#include <stdlib.h>                      // for malloc
#include "Friction_cst.h"                // for SICONOS_FRICTION_3D_ONECONTA...
#include "NumericsFwd.h"                 // for SolverOptions
#include "SolverOptions.h"               // for SolverOptions, solver_option...
#include "frictionContact_test_utils.h"  // for build_test_collection
#include "test_utils.h"                  // for TestCase

TestCase * build_test_collection(int n_data, const char ** data_collection, int* number_of_tests)
{
  *number_of_tests = 22;//n_data * n_solvers;
  TestCase * collection = (TestCase*)malloc((*number_of_tests) * sizeof(TestCase));


  // "External" solver parameters
  // -> same values for all tests.
  // The differences between tests are only for internal solvers and input data.
  int topsolver = SICONOS_FRICTION_3D_NSGS;
  int current = 0;


  {
    // Projection on cone, default values.
    int d = 0; // FC3D_Example1_SBM.dat
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-16;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;

    solver_options_update_internal(collection[current].options, 0, SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone);

    current++;
  }

  {
    int d = 0; // FC3D_Example1_SBM.dat
    // Projection on cone with diagonalization, default value.
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-16;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;

    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithDiagonalization);

    current++;
  }

  {
    int d = 0; // FC3D_Example1_SBM.dat
    // Projection on cone with local iteration, set tol and max iter.
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-16;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;

    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);
    collection[current].options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-3;
    collection[current].options->internalSolvers[0]->iparam[SICONOS_IPARAM_MAX_ITER] = 10;
    current++;
  }

  {
    int d= 0; // FC3D_Example1_SBM.dat
    // Projection on cone with regularization, set rho
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-16;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;

    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization);
    collection[current].options->internalSolvers[0]->dparam[SICONOS_FRICTION_3D_NSN_RHO] = 0.1;
    current++;
  }


  {
    int d=1; // "./data/Capsules-i122-1617.dat"

    // Projection on cone with local iteration, set tol, itermax, d[9], i[8]
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-16;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;

    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);

    collection[current].options->internalSolvers[0]->dparam[SICONOS_IPARAM_MAX_ITER] = 20;
    collection[current].options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-16;
    collection[current].options->internalSolvers[0]->dparam[9] = 1.; // ???
    collection[current].options->internalSolvers[0]->iparam[8] = 1;  // ???
    // Expected to fail ...
    collection[current].will_fail = 1;
    current++;
  }

  {
    int d = 2; // Confeti-ex13-4contact-Fc3D-SBM.dat
    // Projection on cone set d[9], i[8]
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;

    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone);

    collection[current].options->internalSolvers[0]->dparam[9] = 1.; // ???
    collection[current].options->internalSolvers[0]->iparam[8] = 1;  // ???
    // Expected to fail ...
    collection[current].will_fail = 1;
    current++;
  }

  {
    int d = 2; // Confeti-ex13-4contact-Fc3D-SBM.dat
    // nonsmooth newton. Set tol, max iter and i[1]. Default for other parameters
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-12;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;

    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_NSN);

    collection[current].options->internalSolvers[0]->dparam[SICONOS_IPARAM_MAX_ITER] = 10;
    collection[current].options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-18;
    collection[current].options->internalSolvers[0]->iparam[1] = 1;  // ???
    current++;
  }


  {
    int d = 2; // Confeti-ex13-4contact-Fc3D-SBM.dat
    // Projection on cone with local iteration, set tol and maxiter
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-12;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;

    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);

    collection[current].options->internalSolvers[0]->dparam[SICONOS_IPARAM_MAX_ITER] = 100;
    collection[current].options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-6;
    current++;
  }


  {
    int d = 2; // Confeti-ex13-4contact-Fc3D-SBM.dat
    // Projection on cone with regularization, default values.
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-12;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;

    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization);

    current++;
  }

  {
    int d = 2; // Confeti-ex13-4contact-Fc3D-SBM.dat
    // Projection on cone, default values.
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-2;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;

    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone);

    current++;
  }

  {
    int d = 2; // Confeti-ex13-4contact-Fc3D-SBM.dat
    // nonsmooth newton. Set tol and i[1]. Default for other parameters
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;

    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_NSN);
    collection[current].options->internalSolvers[0]->dparam[SICONOS_IPARAM_MAX_ITER] = 10;
    collection[current].options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-16;
    current++;
  }

  {
    int d = 2; // Confeti-ex13-4contact-Fc3D-SBM.dat
    // Projection on cone with local iteration, set tol and maxiter
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-12;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;

    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);
    collection[current].options->internalSolvers[0]->dparam[SICONOS_IPARAM_MAX_ITER] = 100;
    collection[current].options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-6;
    current++;


    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-12;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;

    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);
    collection[current].options->internalSolvers[0]->dparam[SICONOS_IPARAM_MAX_ITER] = 100;
    collection[current].options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-16;
    current++;
  }

  {
    int d = 3; // GFC3D_TwoRods1-condensed.dat
    // nonsmooth newton. Set tol and i[1]. Default for other parameters
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-12;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;

    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_NSN);
    collection[current].options->internalSolvers[0]->dparam[SICONOS_IPARAM_MAX_ITER] = 10;
    collection[current].options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-18;
    current++;
  }


  {
    int d  = 4; // FC3D_Example1.dat

    // nonsmooth newton. Set tol and max iter. Default for other parameters
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-12;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;

    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_NSN);
    collection[current].options->internalSolvers[0]->dparam[SICONOS_IPARAM_MAX_ITER] = 10;
    collection[current].options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-18;
    current++;
  }

  {
    int d = 5; // Confeti-ex03-Fc3D-SBM.dat
    // Projection on cone, default values.
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;

    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone);
    // Expected to fail ...
    collection[current].will_fail = 1;
    current++;
  }

  {
    int d = 5;
    // nonsmooth newton. Set tol and max iter. Default for other parameters
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;

    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_NSN);
    collection[current].options->internalSolvers[0]->dparam[SICONOS_IPARAM_MAX_ITER] = 10;
    collection[current].options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-16;
    // Expected to fail ...
    collection[current].will_fail = 1;
    current++;
  }

  {
    int d = 5;
    // Projection on cone with local iteration, set tol and maxiter
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;

    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);
    collection[current].options->internalSolvers[0]->dparam[SICONOS_IPARAM_MAX_ITER] = 10;
    collection[current].options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-12;
    current++;
  }


  {
    int d = 5;
    // Projection on cone with regularization, set tol and maxiter
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-5;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;

    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization);
    collection[current].options->internalSolvers[0]->dparam[SICONOS_IPARAM_MAX_ITER] = 10;
    collection[current].options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-8;
    // Expected to fail ...
    collection[current].will_fail = 1;
    current++;
  }



  {
    int d = 7;
    // nonsmooth newton 'damped'. Set tol and max iter. Default for other parameters
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-3;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;

    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_NSN_GP);
    collection[current].options->internalSolvers[0]->dparam[SICONOS_IPARAM_MAX_ITER] = 1000;
    collection[current].options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-16;
    // Expected to fail ...
    collection[current].will_fail = 1;
    current++;
  }

  {
    int d = 7;
    // nonsmooth newton. Set tol and max iter. Default for other parameters
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-3;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 2000;

    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_NSN);
    collection[current].options->internalSolvers[0]->dparam[SICONOS_IPARAM_MAX_ITER] = 1000;
    collection[current].options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-16;
    // Expected to fail ...
    collection[current].will_fail = 1;
    current++;
  }

  {
    int d = 7;
    // Projection on cone with local iteration, set tol and maxiter
    collection[current].filename = data_collection[d];
    collection[current].options = solver_options_create(topsolver);
    collection[current].options->dparam[SICONOS_DPARAM_TOL] = 1e-3;
    collection[current].options->iparam[SICONOS_IPARAM_MAX_ITER] = 2000;

    solver_options_update_internal(collection[current].options, 0,
                                   SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);
    collection[current].options->internalSolvers[0]->dparam[SICONOS_IPARAM_MAX_ITER] = 100;
    collection[current].options->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-6;
    // Expected to fail ...
    collection[current].will_fail = 1;
    current++;
  }

  *number_of_tests = current;
  return collection;

}
