/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/



#include "GenericMechanicalProblem.h"
#include "NonSmoothDrivers.h"
#include "numerics_verbose.h"
#include "FrictionContactProblem.h"
#include "LinearComplementarityProblem.h"
#include "RelayProblem.h"
#include "GenericMechanical_Solvers.h"
#include "NumericsMatrix.h"
#include "SparseBlockMatrix.h"
//#define GMP_DEBUG

/* void * solverFC3D; */
/* void * solverEquality; */
/* void * solverLCP; */
/* void * solverMLCP; */

GenericMechanicalProblem * buildEmptyGenericMechanicalProblem()
{
  GenericMechanicalProblem * paux = (GenericMechanicalProblem *)malloc(sizeof(GenericMechanicalProblem));
  paux->firstListElem = 0;
  paux->lastListElem = 0;
  paux->size = 0;
  paux->maxLocalSize = 0;
  return paux;
}

void freeGenericMechanicalProblem(GenericMechanicalProblem * pGMP, unsigned int level)
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
      break;
    }
    case SICONOS_NUMERICS_PROBLEM_LCP:
    {
      free(((LinearComplementarityProblem *)(pElem->problem))->M);
      //  free(((LinearComplementarityProblem *)(pElem->problem))->q);
      break;
    }
    case SICONOS_NUMERICS_PROBLEM_RELAY:
    {
      free(((RelayProblem *)(pElem->problem))->M);
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
  if (level & NUMERICS_GMP_FREE_MATRIX)
  {
    assert(pGMP->M);
    int storageType = pGMP->M->storageType;
    if (storageType == 0)
      free(pGMP->M->matrix0);
    else
      SBMfree(pGMP->M->matrix1, NUMERICS_SBM_FREE_BLOCK | NUMERICS_SBM_FREE_SBM);
    free(pGMP->q);
    free(pGMP->M);
  }
  if (level & NUMERICS_GMP_FREE_GMP)
    free(pGMP);

}
void * addProblem(GenericMechanicalProblem * pGMP, int problemType, int size)
{
  listNumericsProblem * newProblem = (listNumericsProblem*) malloc(sizeof(listNumericsProblem));
  newProblem->nextProblem = 0;
  newProblem->type = problemType;
  newProblem->size = size;
  newProblem->error = 0;
  pGMP->size += size;
  if (size > pGMP->maxLocalSize)
    pGMP->maxLocalSize = size;
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
    pLCP->M = NM_new();
    pLCP->q = (double*) malloc(size * sizeof(double));
    newProblem->q = pLCP->q;
    pLCP->M->storageType = 0; /*local prb is dense*/
    pLCP->M->size0 = size;
    pLCP->M->size1 = size;
    pLCP->size = size;

    break;
  }
  case (SICONOS_NUMERICS_PROBLEM_RELAY):
  {
    newProblem->problem = (void *) malloc(sizeof(RelayProblem));
    RelayProblem * pRelay = (RelayProblem*)newProblem->problem;
    pRelay->M = NM_new();
    pRelay->q = (double*) malloc(size * sizeof(double));
    newProblem->q = pRelay->q;
    pRelay->M->storageType = 0; /*local prb is dense*/
    pRelay->M->size0 = size;
    pRelay->M->size1 = size;
    pRelay->size = size;
    pRelay->lb = (double*)malloc(size * sizeof(double));
    pRelay->ub = (double*)malloc(size * sizeof(double));

    break;
  }
  case (SICONOS_NUMERICS_PROBLEM_EQUALITY):
  {
    newProblem->problem = NULL;
    newProblem->q = (double*) malloc(size * sizeof(double));;
    break;
  }
  case (SICONOS_NUMERICS_PROBLEM_FC3D):
  {
    newProblem->problem = (void *) malloc(sizeof(FrictionContactProblem));
    FrictionContactProblem* pFC3D = (FrictionContactProblem*) newProblem->problem;
    pFC3D->mu = (double*) malloc(sizeof(double));
    pFC3D->M = NM_new();
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
    exit(EXIT_FAILURE);
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
    printf("-->An sub-problem %s.\n", ns_problem_id_to_name(pElem->type));
    pElem = pElem->nextProblem;
  }
  printf("The sparce block matrice is :\n");
  NM_display(pGMP->M);
  printf("The q vector is :\n");
  for (ii = 0; ii < pGMP->size; ii++)
    printf("%e ", pGMP->q[ii]);

  //SBM_print(pGMP->M->matrix1);
  printf("\nEND Display a GenericMechanicalProblem:\n");
}

void genericMechanical_printInFile(GenericMechanicalProblem*  pGMP, FILE* file)
{
  listNumericsProblem * curProblem = pGMP->firstListElem;
  /*Print M*/
  NM_write_in_file(pGMP->M, file);
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

GenericMechanicalProblem * genericMechanical_newFromFile(FILE* file)
{
  size_t nsubProb = 0;
  int prbType = 0;
  int i, posInX, localSize;
  void * prb;

  GenericMechanicalProblem*  pGMP = buildEmptyGenericMechanicalProblem();

  //fscanf(file,"%d\n",&nsubProb);

  pGMP->M = NM_new();
  NM_new_from_file(pGMP->M, file);
  SparseBlockStructuredMatrix* m = pGMP->M->matrix1;

  pGMP->q = (double *) malloc(pGMP->M->size1 * sizeof(double));
  for (i = 0; i < pGMP->M->size1; i++)
  {
    CHECK_IO(fscanf(file, "%lf ", pGMP->q + i));
  }
  nsubProb = m->filled1 - 1;
  posInX = 0;
  for (size_t ii = 0; ii < nsubProb; ii++)
  {
    if (ii)
      posInX = m->blocksize0[ii - 1];
    localSize = m->blocksize0[ii] - posInX;
    CHECK_IO(fscanf(file, "%d\n", &prbType));
    prb = addProblem(pGMP, prbType, localSize);
    if (prbType == SICONOS_NUMERICS_PROBLEM_FC3D)
    {
      CHECK_IO(fscanf(file, "%lf ", ((FrictionContactProblem*)prb)->mu));
    }
  }



#ifdef GMP_DEBUG
  displayGMP(pGMP);
#endif

  return pGMP;
}
