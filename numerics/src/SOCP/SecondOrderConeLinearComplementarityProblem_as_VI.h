/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef SOCLCP_AS_VI_H
#define SOCLCP_AS_VI_H

/*!\file SecondOrderConeLinearComplementarityProblem_as_VI.h
  \brief Definition of a structure to handle with SOCLCP problems.
*/

#include "NumericsFwd.h"  // for SecondOrderConeLinearComplementarityProblem

/** \struct SecondOrderConeLinearComplementarityProblem_as_VI
 * SecondOrderConeLinearComplementarityProblem_as_VI.h
 *
 */
struct SecondOrderConeLinearComplementarityProblem_as_VI {
  /* the VI associated with the FC3D probelem */
  VariationalInequality* vi;
  /* the SOCLCP associated with the VI  */
  SecondOrderConeLinearComplementarityProblem* soclcp;
};

#if defined(__cplusplus)
extern "C" {
#endif

void Function_VI_SOCLCP(void* self, int n, double* x, double* F);

void Projection_VI_SOCLCP(void* viIn, double* x, double* PX);

#if defined(__cplusplus)
}
#endif

#endif
