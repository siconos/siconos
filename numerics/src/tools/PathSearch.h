// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*!\file PathSearch.h
 * \brief Path search related functions and data
 */

#ifndef PATHSEARCH_H
#define PATHSEARCH_H

#include "NMS.h"             // for NMS_data
#include "Newton_methods.h"  // for functions_LSA
#include "NumericsFwd.h"     // for SolverOptions

/** struct ncp_pathsearch_data NCP_PathSearch.h
 * solver specific data
 */
typedef struct {
  NMS_data* data_NMS;           /**< struct for the NMS scheme */
  functions_LSA* lsa_functions; /**< functions for the search */
} pathsearch_data;

#if defined(__cplusplus)
extern "C" {
#endif

/** set some default value for the solver option when the path search
 * algorithm is used
 * \param options the structure to be modified
 */
void pathsearch_set_default(SolverOptions* options);

#if defined(__cplusplus)
}
#endif

#endif
