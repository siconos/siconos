// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*!\file gfc3d_Gams.h
  \brief Various structure to work with GAMS
*/

#ifndef gfc3d_Gams_h
#define gfc3d_Gams_h

#include "NumericsMatrix.h"

/** Structure needed to solve the global friction contact problem */
typedef struct {
  NumericsMatrix* Mlu; /**< M in LU factored form*/
} GFC3D_Gams;

#endif
