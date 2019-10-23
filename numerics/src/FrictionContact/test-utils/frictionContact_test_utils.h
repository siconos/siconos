/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#ifndef FRICTIONCONTACT_TEST_UTILS_H
#define FRICTIONCONTACT_TEST_UTILS_H

#include "SiconosConfig.h" // for BUILD_AS_CPP
#include "GAMSlink.h"
#include "Friction_cst.h"
#include "test_utils.h" // for TestCase

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Solve lcp using parameters and reference from a pre-defined TestCase
      return 1 if the test has succeeded.
  */
  int frictionContact_test_function(TestCase*);


  
  void frictionContact_test_gams_opts(SN_GAMSparams* GP, int solverId);

#if defined(WITH_FCLIB)
  int frictionContact_test_function_hdf5(const char * path, SolverOptions * options);
  int gfc3d_test_function_hdf5(const char* path, SolverOptions* options);
#endif
 
  TestCase * build_test_collection(int n_data, const char ** data_collection, int*);

  void build_friction_test(const char *,
                  int solver_id, int* d_ind, double* dparam, int * i_ind, int* iparam,
                  int internal_solver_id, int * i_d_ind, double * internal_dparam, int * i_i_ind, int * internal_iparam,
                  TestCase* testname);

  void build_gfc3d_test(const char * filename,
                        int solver_id, int* d_ind, double* dparam, int * i_ind, int* iparam,
                        int internal_solver_id, int * i_d_ind, double * internal_dparam, int * i_i_ind, int * internal_iparam,
                        TestCase* testname);

  void build_rfc3d_test(const char * filename,
                        int solver_id, int* d_ind, double* dparam, int * i_ind, int* iparam,
                        int internal_solver_id, int * i_d_ind, double * internal_dparam, int * i_i_ind, int * internal_iparam,
                        TestCase* testname);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP) 
}
#endif

#endif


