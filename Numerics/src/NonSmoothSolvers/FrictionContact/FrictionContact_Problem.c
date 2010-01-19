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
#include <stdlib.h>
#include "FrictionContact_Problem.h"


int frictionContact_printInFile(FrictionContact_Problem*  problem, FILE* file)
{
  if (! problem)
  {
    fprintf(stderr, "Numerics, FrictionContact_Problem printInFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  int i;

  int d  = problem->dimension;
  fprintf(file, "%d\n", d);
  int nc = problem->numberOfContacts;
  fprintf(file, "%d\n", nc);
  printInFile(problem->M, file);
  for (i = 0; i < problem->M->size1; i++)
  {
    fprintf(file, "%32.24e ", problem->q[i]);
  }
  fprintf(file, "\n");
  for (i = 0; i < nc; i++)
  {
    fprintf(file, "%32.24e ", problem->mu[i]);
  }
  fprintf(file, "\n");
  fprintf(file, "%d\n", problem->isComplete);
  return 0;
}

int frictionContact_newFromFile(FrictionContact_Problem* problem, FILE* file)
{
  int nc = 0, d = 0;
  int i;
  fscanf(file, "%d\n", &d);
  problem->dimension = d;
  fscanf(file, "%d\n", &nc);
  problem->numberOfContacts = nc;
  problem->M = (NumericsMatrix *)malloc(sizeof(NumericsMatrix));

  newFromFile(problem->M, file);

  problem->q = (double *) malloc(problem->M->size1 * sizeof(double));
  for (i = 0; i < problem->M->size1; i++)
  {
    fscanf(file, "%lf ", &(problem->q[i]));
  }

  fscanf(file, "\n");
  problem->mu = (double *) malloc(nc * sizeof(double));
  for (i = 0; i < nc; i++)
  {
    fscanf(file, "%lf ", &(problem->mu[i]));
  }
  fscanf(file, "\n");
  fscanf(file, "%d\n", &(problem->isComplete));
  return 0;
}

void freeFrictionContact_problem(FrictionContact_Problem* problem)
{

  freeNumericsMatrix(problem->M);
  free(problem->M);
  free(problem->mu);
  free(problem->q);
  free(problem);
  problem = NULL;

}
