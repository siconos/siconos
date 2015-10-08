/* Siconos-Numerics, Copyright INRIA 2005-2014
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
/*!\file PathSearch.h
 * \brief Path search related functions and data
 */

#ifndef PATHSEARCH_H
#define PATHSEARCH_H

#include "NMS.h"
#include "Newton_Methods.h"

/** struct ncp_pathsearch_data NCP_PathSearch.h
 * solver specific data
 */
typedef struct
{
  NMS_data* data_NMS; /**< struct for the NMS scheme */
  functions_LSA* lsa_functions; /**< functions for the search */
} pathsearch_data;

#include "SiconosConfig.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** set some default value for the solver option when the path search
   * algorithm is used
   * \param options the structure to be modified
   */
  void pathsearch_default_SolverOption(SolverOptions* options);

  /** free solverData specific for the path search
   * \param solverData the struct to free
   */
  void free_solverData_PathSearch(void* solverData);


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
