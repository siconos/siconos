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

|A C| |u| |a| |0|
|   |*| |+| |=| |
|D B| |v| |b| |w|
0<z*v>0
dim(u)=mm
dim(v)=nn

**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "LA.h"
#include "MLCP_Solvers.h"
#include <math.h>

#ifdef HAVE_MLCPSIMPLEX
/*import external implementation*/
#include "external_mlcp_simplex.h"

static int sIsInitialize = 0;
#endif

void mlcp_simplex_initialize(MixedLinearComplementarity_Problem* problem)
{
#ifdef HAVE_MLCPSIMPLEX
  int nn = problem->n;
  int mm = problem->m;
  extern_mlcp_simplex_init_with_M(&nn , &mm, problem->M->matrix0);
  sIsInitialize = 1;
#endif
}
void mlcp_simplex_reset()
{
#ifdef HAVE_MLCPSIMPLEX
  extern_mlcp_simplex_stop();
  sIsInitialize = 0;
#endif
}
/*  tolVar =options->dparam[0];      tolerance to consider that a var is null
 *  tolComp = options->dparam[1];      tolerance to consider that complementarity holds
 *  tolNegVar = options->dparam[2];     tolerance to consider a value is negative
 *  nIterMax = options->iparam[0];     max number of nodes to consider in tree search
*/
void mlcp_simplex(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options)
{
#ifdef HAVE_MLCPSIMPLEX
  double tol ;
  double * workingFloat = options->floatWorkingMem;
  int * workingInt = options->intWorkingMem;
  int lin;
  int npm = (problem->n) + (problem->m);
  int npm2 = npm * npm;
  int NRHS = 1;
  int one = 1;
  int * ipiv;
  int check;
  int DGESVinfo = 1;
  int nn = problem->n;
  int mm = problem->m;

  if (!sIsInitialize)
    extern_mlcp_simplex_init_with_M(&nn , &mm, problem->M->matrix0);

  extern_mlcp_simplex(problem->q, problem->q + nn, z, z + nn, w + nn , info ,  options->iparam , options->dparam);
  for (lin = 0; lin < nn; lin++)
    w[lin] = 0;

  if (!sIsInitialize)
    extern_mlcp_simplex_stop();
#endif
}
