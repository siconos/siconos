/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.

 * Copyright 2016 INRIA.

 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at

 * http://www.apache.org/licenses/LICENSE-2.0

 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/


/** structure for use in PATH */
typedef struct
{
  int n; /**< Size of the problem */
  int nnz; /**< number of non-zero elements */

  double *z; /**< current iterate */
  double *f; /**< current value of the function */

  double *lb; /**< lower bounds */
  double *ub; /**< upper bounds */

  void* problem; /**< problem data */
} SN_generic_path_env;

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Interface to the PATH solver
   * \param mcp_interface data structure for PATH
   * \param z solution
   * \param F value of the function at the solution
   * \param info return code from PATH
   * \param dparam set of double parameters
   * \param iparam set of integer parameters
   */
 void SN_path_interface(MCP_Interface* mcp_interface, double* z, double* F, int* info , double* dparam, int* iparam);


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

