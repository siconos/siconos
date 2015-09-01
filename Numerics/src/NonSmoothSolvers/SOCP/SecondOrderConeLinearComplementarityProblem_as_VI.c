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
#include <math.h>
#include "SecondOrderConeLinearComplementarityProblem_as_VI.h"

#include "projectionOnCone.h"
#include "misc.h"
#include "SiconosBlas.h"

/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"

void Function_VI_SOCLCP(void * self, int n_notused, double *x, double *F)
{
  DEBUG_PRINT("Function_VI_FC3D(void * self, double *x, double *F)\n")
  VariationalInequality * vi = (VariationalInequality *) self;
  SecondOrderConeLinearComplementarityProblem_as_VI* pb = (SecondOrderConeLinearComplementarityProblem_as_VI*)vi->env;
  SecondOrderConeLinearComplementarityProblem * soclcp = pb->soclcp;
  //frictionContact_display(fc3d);

  int n =   soclcp->n;

  cblas_dcopy(n , soclcp->q , 1 , F, 1);
  prodNumericsMatrix(n, n, 1.0, soclcp->M, x, 1.0, F);
}


void Projection_VI_SOCLCP(void *viIn, double *x, double *PX)
{
  DEBUG_PRINT("Projection_VI_SOCLCP(void *viIn, double *x, double *PX)\n")

  VariationalInequality * vi = (VariationalInequality *) viIn;
  SecondOrderConeLinearComplementarityProblem_as_VI* pb = (SecondOrderConeLinearComplementarityProblem_as_VI*)vi->env;
  SecondOrderConeLinearComplementarityProblem * soclcp = pb->soclcp;
  //SecondOrderConeLinearComplementarityProblem_display(soclcp);

  int cone =0;
  int n = soclcp->n;
  cblas_dcopy(n , x , 1 , PX, 1);
  int dim;
  for (cone = 0 ; cone < soclcp->nc  ; ++cone)
  {
    dim=soclcp->coneIndex[cone+1]-soclcp->coneIndex[cone];
    projectionOnSecondOrderCone(&PX[soclcp->coneIndex[cone]], soclcp->mu[cone], dim);
  }
}
