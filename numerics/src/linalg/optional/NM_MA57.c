/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#ifdef WITH_MA57
#include "NM_MA57.h"

#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"

/*#define DEBUG_MESSAGES*/
#include "numerics_errors.h"
#include "siconos_debug.h"

void NM_MA57_free(void* p) {
  NSM_linear_solver_params* params = (NSM_linear_solver_params*)p;
  if (params->linear_solver_data) {
    LBL_Data* lbl = (LBL_Data*)params->linear_solver_data;
    FILE* logfile = lbl->ma57->logfile;
    LBL_Finalize(lbl);
    fclose(logfile);
  }
  // free(params->linear_solver_data);
  params->linear_solver_data = NULL;
}

#endif
