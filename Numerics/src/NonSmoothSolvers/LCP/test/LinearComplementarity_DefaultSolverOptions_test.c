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
#include "NonSmoothDrivers.h"

int main(void)
{
  printf("\n Start of test on Default Solver_Options\n");
  int info = 0 ;
  Solver_Options * options = (Solver_Options *)malloc(sizeof(Solver_Options));

  FILE * finput  =  fopen("./data/lcp_mmc.dat", "r");
  LinearComplementarity_Problem* problem = (LinearComplementarity_Problem *)malloc(sizeof(LinearComplementarity_Problem));

  info = linearComplementarity_newFromFile(problem, finput);



  fclose(finput);


  info = linearComplementarity_setDefaultSolverOptions(problem, options, "PGS_SBM");
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, "PGS");
  printSolverOptions(options);
  deleteSolverOptions(options);


  info = linearComplementarity_setDefaultSolverOptions(problem, options, "RPGS");
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, "QP");
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, "NSQP");
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, "CPG");
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, "PSOR");
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, "Latin");
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, "Latin_w");
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, "Lemke");
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, "Path");
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, "ENUM");
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, "NewtonMin");
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, "NewtonFB");
  printSolverOptions(options);
  deleteSolverOptions(options);




  freeLinearComplementarity_problem(problem);
  free(options);


  printf("\n End of test on Default Solver_Options\n");
  return info;
}
