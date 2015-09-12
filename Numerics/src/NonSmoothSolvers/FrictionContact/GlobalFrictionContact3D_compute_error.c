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

#include "GlobalFrictionContact3D_compute_error.h"

#include "NumericsOptions.h" // for global options
#include "GlobalFrictionContactProblem.h"
#include "GlobalFrictionContact3D_Solvers.h"
#include "projectionOnCone.h"
#include "SiconosLapack.h"
#include <math.h>
#include <assert.h>
#include <float.h>
extern int *Global_ipiv;
extern int  Global_MisInverse;
extern int  Global_MisLU;

int GlobalFrictionContact3D_compute_error(GlobalFrictionContactProblem* problem, double *reaction , double *velocity, double* globalVelocity, double tolerance, double * error)
{
  /* Checks inputs */
  if (problem == NULL || reaction == NULL || velocity == NULL || globalVelocity == NULL)
    numericsError("GlobalFrictionContact3D_compute_error", "null input");

  /* Computes error = dnorm2( GlobalVelocity -M^-1( q + H reaction)*/
  int incx = 1;
  int nc = problem->numberOfContacts;
  int m = nc * 3;
  int n = problem->M->size0;
  double *mu = problem->mu;
  double *q = problem->q;
  NumericsMatrix *H = problem->H;
  NumericsMatrix *M = problem->M;

  double* qtmp = (double*)malloc(n * sizeof(double));
  double* globalVelocitytmp = (double*)malloc(n * sizeof(double));
  cblas_dcopy(n, q, 1, qtmp, 1);

  double alpha = 1.0;
  double beta = 1.0;
  prodNumericsMatrix(m, n, alpha, H, reaction , beta, qtmp);

  /* dense */
  if (M->storageType == 0)
  {
    int infoDGETRS = -1;
    cblas_dcopy(n, qtmp, 1, globalVelocitytmp, 1);
    assert(Global_MisLU);
    DGETRS(LA_NOTRANS, n, 1,  M->matrix0, n, Global_ipiv, globalVelocitytmp , n, &infoDGETRS);
    assert(!infoDGETRS);
  }
  /* SBM */
  else if (M->storageType == 1)
  {
    beta = 0.0;
    if (!Global_MisInverse)
    {
      inverseDiagSBM(M->matrix1);
      Global_MisInverse = 1;
    }
    prodNumericsMatrix(n, n, alpha, M, qtmp , beta, globalVelocitytmp);

  }
  /* coordinate */
  else if (M->storageType == 2)
  {
    cblas_dcopy(n, qtmp, 1, globalVelocitytmp, 1);
    cs_lusol(1, NM_triplet(M), globalVelocitytmp, DBL_EPSILON);
  }


  cblas_daxpy(n , -1.0 , globalVelocity , 1 , globalVelocitytmp, 1);
  *error =   cblas_dnrm2(n , globalVelocitytmp , 1);
  free(qtmp);
  free(globalVelocitytmp);

  cblas_dcopy(m, problem->b, 1, velocity, 1);
  if (H->storageType == 1)
  {
    beta = 1.0;
    SparseBlockStructuredMatrix *Htrans = (SparseBlockStructuredMatrix*)malloc(sizeof(SparseBlockStructuredMatrix));
    transposeSBM(H->matrix1, Htrans);
    prodSBM(n, m, alpha, Htrans, globalVelocity , beta, velocity);
    freeSBM(Htrans);
    free(Htrans);
  }
  else if (H->storageType == 0)
  {
    cblas_dgemv(CblasColMajor,CblasTrans, n, m, 1.0, H->matrix0 , n, globalVelocity , 1, 1.0, velocity, 1);
  }


  double worktmp[3];
  double normUT;
  double rho = 1.0;
  for (int ic = 0 ; ic < nc ; ic++)
  {
    /* Compute the modified local velocity */
    normUT = sqrt(velocity[ic * 3 + 1] * velocity[ic * 3 + 1] + velocity[ic * 3 + 2] * velocity[ic * 3 + 2]);
    worktmp[0] = reaction[ic * 3] - rho * (velocity[ic * 3] + mu[ic] * normUT);
    worktmp[1] = reaction[ic * 3 + 1] - rho * velocity[ic * 3 + 1] ;
    worktmp[2] = reaction[ic * 3 + 2] - rho * velocity[ic * 3 + 2] ;
    projectionOnCone(worktmp, mu[ic]);
    worktmp[0] = reaction[ic * 3] -  worktmp[0];
    worktmp[1] = reaction[ic * 3 + 1] -  worktmp[1];
    worktmp[2] = reaction[ic * 3 + 2] -  worktmp[2];
    *error +=  worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2];
  }
  /*   *error = sqrt(*error); */

  /* Computes error */
  double normq = cblas_dnrm2(n , problem->q , incx);
  *error = *error / (normq + 1.0);

  if (*error > tolerance)
  {
    /*       if (verbose > 0) printf(" Numerics - GlobalFrictionContact3D_compute_error failed: error = %g > tolerance = %g.\n",*error, tolerance); */
    return 1;
  }
  else
    return 0;
}
