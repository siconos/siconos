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

#define GMP_DEBUG

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
    free(pElem->q);
    switch (pElem->type)
    {
    case SICONOS_NUMERICS_PROBLEM_EQUALITY:
    {
      free(((LinearSystemProblem *)(pElem->problem))->M);
      //  free(((LinearSystemProblem *)(pElem->problem))->q);
      break;
    }
    case SICONOS_NUMERICS_PROBLEM_LCP:
    {
      free(((LinearComplementarityProblem *)(pElem->problem))->M);
      //  free(((LinearComplementarityProblem *)(pElem->problem))->q);
      break;
    }
    case SICONOS_NUMERICS_PROBLEM_FC3D:
    {
      free(((FrictionContactProblem*)(pElem->problem))->M);
      free(((FrictionContactProblem*)(pElem->problem))->mu);
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
  if (size > GMP_MAX_SIZE_LOCAL)
  {
    printf("GMP: addProblem local size to big, set GMP_MAX_SIZE_LOCAL\n");
    exit(1);
  }
  listNumericsProblem * newProblem = (listNumericsProblem*) malloc(sizeof(listNumericsProblem));
  newProblem->nextProblem = 0;
  newProblem->type = problemType;
  newProblem->size = size;
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
    LinearComplementarityProblem * pLCP = (LinearComplementarityProblem*)newProblem->problem;
    pLCP->M = (NumericsMatrix*) malloc(sizeof(NumericsMatrix));
    pLCP->q = (double*) malloc(size * sizeof(double));
    newProblem->q = pLCP->q;
    pLCP->M->storageType = 0; /*local prb is dense*/
    pLCP->M->size0 = size;
    pLCP->M->size1 = size;
    pLCP->size = size;

    break;
  }
  case (SICONOS_NUMERICS_PROBLEM_EQUALITY):
  {
    newProblem->problem = (void *) malloc(sizeof(LinearSystemProblem));
    LinearSystemProblem* pLS = (LinearSystemProblem*) newProblem->problem;
    pLS->size = size;
    pLS->M = (NumericsMatrix*) malloc(sizeof(NumericsMatrix));
    pLS->q = (double*) malloc(size * sizeof(double));
    newProblem->q = pLS->q;
    pLS->M->storageType = 0; /*local prb is dense*/
    pLS->M->size0 = size;
    pLS->M->size1 = size;
    break;
  }
  case (SICONOS_NUMERICS_PROBLEM_FC3D):
  {
    newProblem->problem = (void *) malloc(sizeof(FrictionContactProblem));
    FrictionContactProblem* pFC3D = (FrictionContactProblem*) newProblem->problem;
    pFC3D->mu = (double*) malloc(sizeof(double));
    pFC3D->M = (NumericsMatrix*) malloc(sizeof(NumericsMatrix));
    pFC3D->M->storageType = 0; /*Local prb is dense*/
    pFC3D->M->size0 = size;
    pFC3D->M->size1 = size;
    pFC3D->numberOfContacts = 1;
    pFC3D->q = (double*) malloc(size * sizeof(double));
    pFC3D->dimension = 3;
    newProblem->q = pFC3D->q;
    break;
  }
  default:
    printf("GenericMechanicalProblem.h addProblem : problemType unknown: %d . \n", problemType);
  }
  return  newProblem->problem;
}


void displayGMP(GenericMechanicalProblem * pGMP)
{
  listNumericsProblem * pElem = pGMP->firstListElem;
  int ii;
  printf("\nBEGIN Display a GenericMechanicalProblem(Numerics):\n");

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

void genericMechnical_printInFile(GenericMechanicalProblem*  pGMP, FILE* file)
{
  listNumericsProblem * curProblem = pGMP->firstListElem;
  /*Print M*/
  printInFile(pGMP->M, file);
  fprintf(file, "\n");
  /*Print Q*/
  for (int ii = 0; ii < pGMP->size; ii++)
    fprintf(file, "%e\n", pGMP->q[ii]);
  fprintf(file, "\n");
  /*Print lthe type and options (mu)*/
  while (curProblem)
  {
    fprintf(file, "%d\n", curProblem->type);
    if (curProblem->type == SICONOS_NUMERICS_PROBLEM_FC3D)
      fprintf(file, "%e\n", ((FrictionContactProblem*)curProblem->problem)->mu[0]);
    curProblem = curProblem->nextProblem;
  }
}

GenericMechanicalProblem * genericMechnical_newFromFile(FILE* file)
{
  int nsubProb = 0;
  int prbType = 0;
  int i, posInX, localSize;
  void * prb;

  GenericMechanicalProblem*  pGMP = buildEmptyGenericMechanicalProblem();

  //fscanf(file,"%d\n",&nsubProb);

  pGMP->M = (NumericsMatrix *)malloc(sizeof(NumericsMatrix));
  newFromFile(pGMP->M, file);
  SparseBlockStructuredMatrix* m = pGMP->M->matrix1;

  pGMP->q = (double *) malloc(pGMP->M->size1 * sizeof(double));
  for (i = 0; i < pGMP->M->size1; i++)
  {
    fscanf(file, "%lf ", pGMP->q + i);
  }
  nsubProb = m->filled1 - 1;
  posInX = 0;
  for (int ii = 0; ii < nsubProb; ii++)
  {
    if (ii)
      posInX = m->blocksize0[ii - 1];
    localSize = m->blocksize0[ii] - posInX;
    fscanf(file, "%d\n", &prbType);
    prb = addProblem(pGMP, prbType, localSize);
    if (prbType == SICONOS_NUMERICS_PROBLEM_FC3D)
    {
      fscanf(file, "%lf ", ((FrictionContactProblem*)prb)->mu);
    }
  }



#ifdef GMP_DEBUG
  displayGMP(pGMP);
#endif

  return pGMP;
}
