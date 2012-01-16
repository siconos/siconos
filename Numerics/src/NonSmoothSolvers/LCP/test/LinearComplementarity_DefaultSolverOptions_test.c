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
#include "NonSmoothDrivers.h"

int main(void)
{
  printf("\n Start of test on Default SolverOptions\n");
  int info = 0 ;
  SolverOptions * options = (SolverOptions *)malloc(sizeof(SolverOptions));

  FILE * finput  =  fopen("./data/lcp_mmc.dat", "r");
  LinearComplementarityProblem* problem = (LinearComplementarityProblem *)malloc(sizeof(LinearComplementarityProblem));

  info = linearComplementarity_newFromFile(problem, finput);

  fclose(finput);


  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_NSGS_SBM);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_PGS);
  printSolverOptions(options);
  deleteSolverOptions(options);


  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_RPGS);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_QP);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_NSQP);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_CPG);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_PSOR);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_LATIN);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_LATIN_W);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_LEMKE);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_PATH);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_ENUM);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_NEWTONMIN);
  printSolverOptions(options);
  deleteSolverOptions(options);

  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_NEWTONFB);
  printSolverOptions(options);
  deleteSolverOptions(options);




  freeLinearComplementarityProblem(problem);
  free(options);


  printf("\n End of test on Default SolverOptions\n");
  return info;
}
