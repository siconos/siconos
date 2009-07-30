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
#include "LA.h"
#include "MLCP_Solvers.h"
#include <math.h>
#include "mlcp_path_enum.h"
#include "mlcp_enum.h"
#include "mlcp_tool.h"

static int sN;
static int sM;

static int * siWorkEnum = 0;
static int * siWorkPath = 0;
static double * sdWorkEnum = 0;
static double * sdWorkPath = 0;



int mlcp_path_enum_getNbIWork(MixedLinearComplementarity_Problem* problem, Solver_Options* options)
{
  return mlcp_enum_getNbIWork(problem, options);
}
int mlcp_path_enum_getNbDWork(MixedLinearComplementarity_Problem* problem, Solver_Options* options)
{
  return mlcp_enum_getNbDWork(problem, options);
}



/*
 *options->iparam[5] : n0 number of possible configuration.
 * dparam[5] : (in) a positive value, tolerane about the sign.
 *options->iWork : double work memory of  mlcp_direct_enum_getNbIWork() integers  (2(nn+mm))+((n + m)*(n0+1) + nO*m)
 *options->dWork : double work memory of mlcp_direct_enum_getNbDWork() doubles  ((nn+mm)*(nn+mm) + 3*(nn+mm))+(n + m + n0*(n+m)*(n+m))
 *
 *
 */

void mlcp_path_enum_init(MixedLinearComplementarity_Problem* problem, Solver_Options* options)
{
  sN = problem->n;
  sM = problem->m;
  int iOffset = 0;/* mlcp_path_getNbIWork(problem,options);*/
  int dOffset = 0;/*mlcp_path_getNbDWork(problem,options);*/
  siWorkEnum = options->iWork + iOffset;
  siWorkPath = options->iWork;
  sdWorkEnum = options->dWork + dOffset;
  sdWorkPath = options->dWork;
  /*  mlcp_path_init(problem, options);*/

}
void mlcp_path_enum_reset()
{
  /*mlcp_path_reset();*/
  siWorkEnum = 0;
  siWorkPath = 0;
  sdWorkEnum = 0;
  sdWorkPath = 0;
}

/*
 * The are no memory allocation in mlcp_direct, all necessary memory must be allocated by the user.
 *
 *options:
 * iparam[0] : (in) verbose.
 * dparam[0] : (in) a positive value, tolerane about the sign used by the enum algo.
 * iparam[5] : (in)  n0 number of possible configuration.
 * dparam[5] : (in) a positive value, tolerane about the sign.
 * dWork : working float zone size : n + m + n0*(n+m)*(n+m)  . MUST BE ALLOCATED BY THE USER.
 * iWork : working int zone size : (n + m)*(n0+1) + nO*m. MUST BE ALLOCATED BY THE USER.
 * double *z : size n+m
 * double *w : size n+m
 * info : output. info == 0 if success
 */
void mlcp_path_enum(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options)
{
  if (!siWorkEnum)
  {
    *info = 1;
    printf("MLCP_PATH_ENUM error, call a non initialised method!!!!!!!!!!!!!!!!!!!!!\n");
    return;
  }
  /*First, try direct solver*/
  //  options->dWork = sdWorkDirect;
  //  options->iWork = siWorkDirect;
  mlcp_path(problem, z, w, info, options);
  if (*info)
  {
    printf("MLCP_PATH_ENUM: path failed, call enum\n");
    options->dWork = sdWorkEnum;
    options->iWork = siWorkEnum;
    /*solver direct failed, so run the enum solver.*/
    mlcp_enum(problem, z, w, info, options);
  }
}
