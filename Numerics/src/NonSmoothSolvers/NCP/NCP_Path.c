/* Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Standalone_Path.h"

static int fill_structure;  /* Do we need to fill in the structure of    */
/* the Jacobian?                             */

int NCP_Path(int n, double* z, FuncEvalPtr* F, JacEvalPtr* jacobianF, int* iparam, double* dparam)
{
  /* Lower bounds on the variables = 0 for NCP */
  double *lb  = (double *)malloc(sizeof(double) * n);
  /* Upper bounds on the variables = +inf for NCP */
  double *ub = (double *)malloc(sizeof(double) * n);
  /* Function evaluation */
  double *f = (double *)malloc(sizeof(double) * n);;
  /* Number of nonzeros */
  int nnz = n * n;
  /* Termination status from PATH */
  int status;
  int i, j;

  /* Connect F and its jacobian to input functions */
  setFuncEval(*F);
  setJacEval(*jacobianF);

  /**************************************/
  /* Fill in the lower and upper bounds */
  /**************************************/

  for (j = 0; j < n; j++)
  {
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

