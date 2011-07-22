/* Siconos-Numerics, Copyright INRIA 2005-2011.
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include "LA.h"
#include "NumericsOptions.h"
#include "LCP_Solvers.h"
#include "NonSmoothDrivers.h"

int linearComplementarity_setDefaultSolverOptions(LinearComplementarityProblem* problem, SolverOptions* options, int solverId)
{

  int info = -1;
  switch (solverId)
  {
  case SICONOS_LCP_NSGS_SBM:
  {
    info =    linearComplementarity_nsgs_SBM_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_LCP_QP:
  {
    info =    linearComplementarity_qp_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_LCP_NSQP:
  {
    info =    linearComplementarity_nsqp_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_LCP_CPG:
  {
    info =    linearComplementarity_cpg_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_LCP_PGS:
  {
    info =    linearComplementarity_pgs_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_LCP_RPGS:
  {
    info =    linearComplementarity_rpgs_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_LCP_PSOR:
  {
    info =    linearComplementarity_psor_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_LCP_LATIN:
  {
    info =    linearComplementarity_latin_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_LCP_LATIN_W:
  {
    info =    linearComplementarity_latin_w_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_LCP_LEMKE:
  {
    info =    linearComplementarity_lexicolemke_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_LCP_PATH:
  {
    info =    linearComplementarity_path_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_LCP_ENUM:
  {
    info =    linearComplementarity_enum_setDefaultSolverOptions(problem, options);
    break;
  }
  case SICONOS_LCP_NEWTONMIN:
  {
    info =    linearComplementarity_newton_min_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_LCP_NEWTONFB:
  {
    info =    linearComplementarity_newton_FB_setDefaultSolverOptions(options);
    break;
  }
  default:
  {
    numericsError("linearComplementarity_setDefaultSolverOptions", "Unknown Solver");

  }
  }


  return info;
}
