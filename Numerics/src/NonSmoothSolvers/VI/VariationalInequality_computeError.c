/* Siconos-Numerics, Copyright INRIA 2005-2012.
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


#include "NumericsOptions.h" // for global options
#include "VariationalInequality.h"
#include "SolverOptions.h"
#include "VariationalInequality_computeError.h"
#include "SiconosLapack.h"

#include <math.h>
#include <assert.h>

/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"

int variationalInequality_computeError(
  VariationalInequality* problem,
  double *z , double *w, double tolerance,
  SolverOptions * options, double * error)
{

  assert(problem);
  assert(z);
  assert(w);
  assert(error);
  
  int incx = 1;

  int n = problem->size;

  *error = 0.;
  double *ztmp = (double*)malloc(n* sizeof(double));

  for (int i=0; i<n ; i++)
  {
    ztmp[i]=0.0;
  }
  problem->F(problem,ztmp,w);
  double normq = cblas_dnrm2(n , w , incx);
  DEBUG_PRINTF("normq = %12.8e\n", normq);
  cblas_dcopy(n , z , 1 , ztmp, 1);

  problem->F(problem,ztmp,w);
 
  
  cblas_daxpy(n, -1.0, w , 1, ztmp , 1) ;
  
  problem->ProjectionOnX(problem,ztmp,w);
  
  cblas_dcopy(n , z , 1 , ztmp, 1);

  cblas_daxpy(n, -1.0, w , 1, ztmp , 1) ;

  *error = cblas_dnrm2(n , ztmp , incx);
  free(ztmp);
  
  problem->F(problem,z,w);


  /* Computes error */
  *error = *error / (normq + 1.0);
  if (*error > tolerance)
  {
    if (verbose > 1)
      printf(" Numerics - variotionalInequality_compute_error: error = %g > tolerance = %g.\n",
             *error, tolerance);
    return 1;
  }
  else
    return 0;
}




int variationalInequality_computeError_wait(
  VariationalInequality* problem,
  double *z , double *w, double tolerance,
  SolverOptions * options, double * error)
{

  assert(problem);
  assert(z);
  assert(w);
  assert(error);
  
  int incx = 1;

  int n = problem->size;

  *error = 0.;
  double *ztmp = (double*)malloc(n* sizeof(double));
  cblas_dcopy(n , z , 1 , ztmp, 1);

  problem->F(problem,ztmp,w);
  double normq = cblas_dnrm2(n , w , incx);
  
  
  cblas_daxpy(n, -1.0, w , 1, ztmp , 1) ;
  
  problem->ProjectionOnX(problem,ztmp,w);
  
  cblas_dcopy(n , z , 1 , ztmp, 1);

  cblas_daxpy(n, -1.0, w , 1, ztmp , 1) ;

  *error = cblas_dnrm2(n , ztmp , incx);
  free(ztmp);
  
  problem->F(problem,z,w);


  
  /* Computes error */
  *error = *error / (normq + 1.0);
  if (*error > tolerance)
  {
    if (verbose > 1)
      printf(" Numerics - variotionalInequality_compute_error: error = %g > tolerance = %g.\n",
             *error, tolerance);
    return 1;
  }
  else
    return 0;
}
