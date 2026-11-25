// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "NCP_Path.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Standalone_Path.h"
#include "numerics_errors.h"

static int fill_structure; /* Do we need to fill in the structure of    */
/* the Jacobian?                             */

int NCP_Path(int n, double* z, FuncEvalPtr F, JacEvalPtr jacobianF, int* iparam,
             double* dparam) {
  /* Lower bounds on the variables = 0 for NCP */
  double* lb = (double*)malloc(sizeof(double) * n);
  /* Upper bounds on the variables = +inf for NCP */
  double* ub = (double*)malloc(sizeof(double) * n);
  /* Function evaluation */
  double* f = (double*)malloc(sizeof(double) * n);
  ;
  /* Number of nonzeros */
  int nnz = n * n;
  /* Termination status from PATH */
  int status;
  int j;

  /* Connect F and its jacobian to input functions */
  setFuncEval(F);
  setJacEval(jacobianF);

  /**************************************/
  /* Fill in the lower and upper bounds */
  /**************************************/

  for (j = 0; j < n; j++) {
    lb[j] = 0;
    ub[j] = 1e20;
  }

  /************************************************************************/
  /* Set fill_structure to true (we ALWAYS need to fill in the structure  */
  /* of the jacobian at least once per call to PATH).                     */
  /************************************************************************/

  fill_structure = 1;

  /************************************************************************/
  /* Call PATH.                                                           */
  /************************************************************************/

  pathMain(n, nnz, &status, z, f, lb, ub);

  /************************************************************************/
  /* Deallocate memory.                                                   */
  /************************************************************************/

  free(lb);
  free(ub);
  free(f);
  return 0;
}
