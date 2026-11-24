// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/** structure for use in PATH */
typedef struct {
  int n;   /**< Size of the problem */
  int nnz; /**< number of non-zero elements */

  double* z; /**< current iterate */
  double* f; /**< current value of the function */

  double* lb; /**< lower bounds */
  double* ub; /**< upper bounds */

  void* problem; /**< problem data */
} SN_generic_path_env;

#if defined(__cplusplus)
extern "C" {
#endif

/** Interface to the PATH solver
 * \param mcp_interface data structure for PATH
 * \param z solution
 * \param F value of the function at the solution
 * \param info return code from PATH
 * \param dparam set of double parameters
 * \param iparam set of integer parameters
 */
void SN_path_interface(MCP_Interface* mcp_interface, double* z, double* F, int* info,
                       double* dparam, int* iparam);

#if defined(__cplusplus)
}
#endif
