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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include "LA.h"
#include "Numerics_Options.h"
#include "LCP_Solvers.h"
#include "NonSmoothDrivers.h"

int linearComplementarity_setDefaultSolverOptions(LinearComplementarity_Problem* problem, Solver_Options** arrayOfSolver_Options, char *solvername)
{
  int info = -1;
  if (strcmp(solvername , "GaussSeidel_SBM") == 0)
  {
    info =    linearComplementarity_qp_setDefaultSolverOptions(arrayOfSolver_Options);
  }
  else if (strcmp(solvername, "QP") == 0)
  {
    info =    linearComplementarity_qp_setDefaultSolverOptions(arrayOfSolver_Options);
  }
  else if (strcmp(solvername, "NSQP") == 0)
  {
    info =    linearComplementarity_nsqp_setDefaultSolverOptions(arrayOfSolver_Options);
  }
  else if (strcmp(solvername, "CPG") == 0)
  {
    info =    linearComplementarity_cpg_setDefaultSolverOptions(arrayOfSolver_Options);
  }
  else if (strcmp(solvername, "PGS") == 0)
  {
    info =    linearComplementarity_pgs_setDefaultSolverOptions(arrayOfSolver_Options);
  }
  else if (strcmp(solvername, "RPGS") == 0)
  {
    info =    linearComplementarity_rpgs_setDefaultSolverOptions(arrayOfSolver_Options);
  }
  else if (strcmp(solvername, "PSOR") == 0)
  {
    info =    linearComplementarity_psor_setDefaultSolverOptions(arrayOfSolver_Options);
  }
  else if (strcmp(solvername, "Latin") == 0)
  {
    info =    linearComplementarity_latin_setDefaultSolverOptions(arrayOfSolver_Options);
  }
  else if (strcmp(solvername, "Latin_w") == 0)
  {
    info =    linearComplementarity_latin_w_setDefaultSolverOptions(arrayOfSolver_Options);
  }
  else if (strcmp(solvername , "Lemke") == 0 || strcmp(solvername , "LexicoLemke") == 0)
  {
    info =    linearComplementarity_lexicolemke_setDefaultSolverOptions(arrayOfSolver_Options);
  }
  else if (strcmp(solvername, "Path") == 0)
  {
    info =    linearComplementarity_path_setDefaultSolverOptions(arrayOfSolver_Options);
  }
  else if (strcmp(solvername, "ENUM") == 0)
  {
    info =    linearComplementarity_enum_setDefaultSolverOptions(problem, arrayOfSolver_Options);
  }
  else if (strcmp(solvername, "NewtonMin") == 0)
  {
    info =    linearComplementarity_newton_min_setDefaultSolverOptions(arrayOfSolver_Options);
  }
  else if (strcmp(solvername, "NewtonFB") == 0)
  {
    info =    linearComplementarity_newton_FB_setDefaultSolverOptions(arrayOfSolver_Options);
  }
  else
  {
    numericsError("linearComplementarity_setDefaultSolverOptions", "Unknow Solver");

  }


  return info;
}

int linearComplementarity_deleteDefaultSolverOptions(Solver_Options** arrayOfSolver_Options, char *solvername)
{
  int info = -1;
  if (strcmp(solvername , "GaussSeidel_SBM") == 0)
  {
    info =    linearComplementarity_qp_deleteDefaultSolverOptions(arrayOfSolver_Options);
  }
  else if (strcmp(solvername, "QP") == 0)
  {
    info =    linearComplementarity_qp_deleteDefaultSolverOptions(arrayOfSolver_Options);
  }
  else if (strcmp(solvername, "NSQP") == 0)
  {
    info =    linearComplementarity_nsqp_deleteDefaultSolverOptions(arrayOfSolver_Options);
  }
  else if (strcmp(solvername, "CPG") == 0)
  {
    info =    linearComplementarity_cpg_deleteDefaultSolverOptions(arrayOfSolver_Options);
  }
  else if (strcmp(solvername, "PGS") == 0)
  {
    info =    linearComplementarity_pgs_deleteDefaultSolverOptions(arrayOfSolver_Options);
  }
  else if (strcmp(solvername, "RPGS") == 0)
  {
    info =    linearComplementarity_rpgs_deleteDefaultSolverOptions(arrayOfSolver_Options);
  }
  else if (strcmp(solvername, "PSOR") == 0)
  {
    info =    linearComplementarity_psor_deleteDefaultSolverOptions(arrayOfSolver_Options);
  }
  else if (strcmp(solvername, "Latin") == 0)
  {
    info =    linearComplementarity_latin_deleteDefaultSolverOptions(arrayOfSolver_Options);
  }
  else if (strcmp(solvername, "Latin_w") == 0)
  {
    info =    linearComplementarity_latin_w_deleteDefaultSolverOptions(arrayOfSolver_Options);
  }
  else if (strcmp(solvername , "Lemke") == 0 || strcmp(solvername , "LexicoLemke") == 0)
  {
    info =    linearComplementarity_lexicolemke_deleteDefaultSolverOptions(arrayOfSolver_Options);
  }
  else if (strcmp(solvername, "Path") == 0)
  {
    info =    linearComplementarity_path_deleteDefaultSolverOptions(arrayOfSolver_Options);
  }
  else if (strcmp(solvername, "ENUM") == 0)
  {
    info =    linearComplementarity_enum_deleteDefaultSolverOptions(arrayOfSolver_Options);
  }
  else if (strcmp(solvername, "NewtonMin") == 0)
  {
    info =    linearComplementarity_newton_min_deleteDefaultSolverOptions(arrayOfSolver_Options);
  }
  else if (strcmp(solvername, "NewtonFB") == 0)
  {
    info =    linearComplementarity_newton_FB_deleteDefaultSolverOptions(arrayOfSolver_Options);
  }
  else
  {
    numericsError("linearComplementarity_deleteDefaultSolverOptions", "Unknow Solver");

  }


  return info;
}
