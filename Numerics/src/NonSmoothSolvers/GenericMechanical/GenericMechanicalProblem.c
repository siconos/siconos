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



#include "GenericMechanicalProblem.h"
#include "NonSmoothDrivers.h"
/* void * solverFC3D; */
/* void * solverEquality; */
/* void * solverLCP; */
/* void * solverMLCP; */

/** GenericMechanicalProblem
    \param numberOfBlockLine The number of line of blocks.
    \param M A sparse blocks matrix.
    \param q A dense vector.
    \param problems An array of pointer on Numerics problem, either FrictionContactProblem, LinearSystemProblem, MixedLinearComplementarityProblem or LinearComplementarityProblem.
    Each of them contain only the parameters value that can not be decuced from the corresponding block lines (n and m for the MLCP, mu and e for the Friction problem). The size is numberOfBlockLine.
    \param problemsType
 */

/** Remark:

    The M and q contains de matrices of the problem. The sub problems (problems) has also a M and q member usfull for the computation of the local error.


 */
GenericMechanicalProblem * buildEmptyGenericMechanicalProblem()
{
  GenericMechanicalProblem * paux = malloc(sizeof(GenericMechanicalProblem));
  paux->firstListElem = 0;
  paux->lastListElem = 0;
  paux->size = 0;
  return paux;
}

void freeGenericMechanicalProblem(GenericMechanicalProblem * pGMP)
{
  if (!pGMP)
    return;
  while (pGMP->lastListElem)
  {
    listNumericsProblem * pElem = pGMP->lastListElem;
    switch (pElem->type)
    {
    case SICONOS_NUMERICS_PROBLEM_EQUALITY:
    {
      free(((LinearSystemProblem *)(pElem->problem))->M);
      free(((LinearSystemProblem *)(pElem->problem))->q);
      break;
    }
    case SICONOS_NUMERICS_PROBLEM_FC3D:
    {
      free(((FrictionContactProblem*)(pElem->problem))->M);
      free(((FrictionContactProblem*)(pElem->problem))->mu);
      printf("Numerics : freeGenericMechanicalProblem missing free\n");
      break;
    }
    default:
      printf("Numerics : freeGenericMechanicalProblem case %d not managed.\n", pElem->type);
    }

    free(pElem->problem);
    pGMP->lastListElem = pElem->prevProblem;
    free(pElem);
  }
  free(pGMP);
}
void * addProblem(GenericMechanicalProblem * pGMP, int problemType, int size)
{
  listNumericsProblem * newProblem = (listNumericsProblem*) malloc(sizeof(listNumericsProblem));
  newProblem->nextProblem = 0;
  newProblem->type = problemType;
  pGMP->size += size;
  if (!pGMP->lastListElem)
  {
    pGMP->firstListElem = newProblem;
    pGMP->lastListElem = newProblem;
    newProblem->prevProblem = 0;
  }
  else
  {
    pGMP->lastListElem->nextProblem = newProblem;
    newProblem->prevProblem =  pGMP->lastListElem;
    pGMP->lastListElem = newProblem;
  }
  switch (problemType)
  {
  case (SICONOS_NUMERICS_PROBLEM_LCP):
  {
    newProblem->problem = (void *) malloc(sizeof(LinearComplementarityProblem));
    printf("NUMERICS::::Warning1, missing code to allocate q.\n");
    break;
  }
  case (SICONOS_NUMERICS_PROBLEM_EQUALITY):
  {
    newProblem->problem = (void *) malloc(sizeof(LinearSystemProblem));
    ((LinearSystemProblem*)(newProblem->problem))->size = size;
    ((LinearSystemProblem*)(newProblem->problem))->M = (NumericsMatrix*) malloc(sizeof(NumericsMatrix));
    ((LinearSystemProblem*)(newProblem->problem))->q = (double*) malloc(size * sizeof(double));
    break;
  }
  case (SICONOS_NUMERICS_PROBLEM_FC3D):
  {
    newProblem->problem = (void *) malloc(sizeof(FrictionContactProblem));
    ((FrictionContactProblem*)(newProblem->problem))->mu = (double*) malloc(sizeof(double));
    ((FrictionContactProblem*)(newProblem->problem))->M = (NumericsMatrix*) malloc(sizeof(NumericsMatrix));
    printf("NUMERICS::::Warning2, missing code to allocate q.\n");
    break;
  }
  default:
    printf("GenericMechanicalProblem.h addProblem : problemType unknown: %d . \n", problemType);
  }
  return 0;
}


void displayGMP(GenericMechanicalProblem * pGMP)
{
  listNumericsProblem * pElem = pGMP->firstListElem;
  int ii;
  printf("\nBEGIN Display a GenericMechanicalProblem:\n");

  while (pElem)
  {
    printf("-->An sub-problem %s.\n", idProblemToChar(pElem->type));
    pElem = pElem->nextProblem;
  }
  printf("The sparce block matrice is :\n");
  display(pGMP->M);
  printf("The q vector is :\n");
  for (ii = 0; ii < pGMP->size; ii++)
    printf("%e ", pGMP->q[ii]);

  //printSBM(pGMP->M->matrix1);
  printf("\nEND Display a GenericMechanicalProblem:\n");
}
