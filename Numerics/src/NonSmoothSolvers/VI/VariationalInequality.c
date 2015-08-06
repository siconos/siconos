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
#include "VariationalInequality.h"
#include "misc.h"


void variationalInequality_display(VariationalInequality* problem)
{

  assert(problem);

}

int variationalInequality_printInFile(VariationalInequality*  problem, FILE* file)
{
  if (! problem)
  {
    fprintf(stderr, "Numerics, VariationalInequality printInFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }

  return 0;
}

int variationalInequality_newFromFile(VariationalInequality* problem, FILE* file)
{

  return 0;
}

void freeVariationalInequalityProblem(VariationalInequality* problem)
{
  assert(problem);
}

/* VariationalInequality* variationalInequalityProblem_new( int dim, void (* F)(void *vi, double *,double *)  ) */
/* { */
/*   VariationalInequality* fvi = (VariationalInequality*) malloc(sizeof(VariationalInequality)); */

/*   fvi->size = dim; */
/*   fvi->F = Callback; */
/*   return fvi; */
/* } */
