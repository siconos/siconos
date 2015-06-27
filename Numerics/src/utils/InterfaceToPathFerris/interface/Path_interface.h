/* Siconos-Numerics, Copyright INRIA 2005-2015
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
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

