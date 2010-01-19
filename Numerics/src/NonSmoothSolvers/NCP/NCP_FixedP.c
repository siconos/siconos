/* Siconos-Numerics, Copyright INRIA 2005-2010.
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
#include "NonSmoothNewton.h"
#include "NCP_Path.h"
#include "FrictionContact3D_NCPGlockerFixedPoint.h"
#include "Numerics_Options.h"
#include "LA.h"
#include "FrictionContact3D2NCP_Glocker.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"



/*============================ Fixed point Solver ==================================*/

int Fixe(int n, double* z, int* iparam, double* dparam)
{

  int itermax = iparam[0]; // maximum number of iterations allowed
  int niter = 0; // current iteration number
  double tolerance = dparam[0];
  if (verbose > 0)
  {
    printf(" ============= Starting of Newton process =============\n");
    printf(" - tolerance: %f\n - maximum number of iterations: %i\n", tolerance, itermax);
  }

  int i;

  /* Connect F and its jacobian to input functions */
  //  setFuncEval(F);

  /* Memory allocation for phi and its jacobian */
  double * www = (double*)malloc(sizeof(double) * n);

  double terminationCriterion = 1;

  /** Iterations ... */
  while ((niter < itermax) && (terminationCriterion > tolerance))
  {
    ++niter;

    //printf(" ============= Fixed Point Iteration ============= %i\n",niter);
    for (i = 0; i < n ; ++i)
      compute_Z_GlockerFixedP(i, www);

    terminationCriterion = DNRM2(n, www, 1);
    //printf(" error = %14.7e\n", terminationCriterion);
    if (verbose > 0)
    {
      printf("Non Smooth Newton, iteration number %i, error equal to %14.7e .\n", niter, terminationCriterion);
      printf(" -----------------------------------------------------------------------");
    }
  }

  /* Total number of iterations */
  iparam[1] = niter;
  /* Final error */
  dparam[1] = terminationCriterion;

  /** Free memory*/
  free(www);

  if (verbose > 0)
  {
    if (dparam[1] > tolerance)
      printf("Non Smooth Newton warning: no convergence after %i iterations\n" , niter);

    else
      printf("Non Smooth Newton: convergence after %i iterations\n" , niter);
    printf(" The residue is : %e \n", dparam[1]);
  }

  if (dparam[1] > tolerance)
    return 1;
  else return 0;
}
