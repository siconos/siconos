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
#include <stdlib.h>
#include <assert.h>
#include "FrictionContactProblem_as_VI.h"

#include "projectionOnCone.h"
#include "misc.h"
#include "SiconosBlas.h"

/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"

void Function_VI_FC3D(void * self, double *x, double *F)
{
  DEBUG_PRINT("Function_VI_FC3D(void * self, double *x, double *F)\n")
  VariationalInequality * vi = (VariationalInequality *) self;
  FrictionContactProblem_as_VI* pb = (FrictionContactProblem_as_VI*)vi->env;
  FrictionContactProblem * fc3d = pb->fc3d;
  //frictionContact_display(fc3d);

  int nLocal =  fc3d->dimension;
  int n = fc3d->numberOfContacts *  fc3d->dimension;

  cblas_dcopy(n , fc3d->q , 1 , F, 1);
  prodNumericsMatrix(n, n, 1.0, fc3d->M, x, 1.0, F);
  int contact =0;

  for (contact = 0 ; contact <  fc3d->numberOfContacts ; ++contact)
  {
    double  normUT = sqrt(F[contact * nLocal + 1] * F[contact * nLocal + 1]
                          + F[contact * nLocal + 2] * F[contact * nLocal + 2]);
    F[contact * nLocal] +=  (fc3d->mu[contact] * normUT);
  }
}


void Projection_VI_FC3D(void *viIn, double *x, double *PX)
{
  DEBUG_PRINT("Projection_VI_FC3D(void *viIn, double *x, double *PX)\n")

  VariationalInequality * vi = (VariationalInequality *) viIn;
  FrictionContactProblem_as_VI* pb = (FrictionContactProblem_as_VI*)vi->env;
  FrictionContactProblem * fc3d = pb->fc3d;
  //frictionContact_display(fc3d);

  int contact =0;
  int nLocal =  fc3d->dimension;
  int n = fc3d->numberOfContacts* nLocal;
  cblas_dcopy(n , x , 1 , PX, 1);
  for (contact = 0 ; contact < fc3d->numberOfContacts  ; ++contact)
  {
    projectionOnCone(&PX[ contact * nLocal ], fc3d->mu[contact]);
  }
}
