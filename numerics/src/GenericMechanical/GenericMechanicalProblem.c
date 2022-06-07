/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#include <assert.h>                        // for assert
#include <stdlib.h>                        // for malloc, free, exit, EXIT_F...
#include "FrictionContactProblem.h"        // for FrictionContactProblem
#include "GenericMechanical_Solvers.h"     // for NUMERICS_GMP_FREE_GMP, NUM...
#include "LinearComplementarityProblem.h"  // for LinearComplementarityProblem
#include "NumericsMatrix.h"                // for NumericsMatrix, NM_new
#include "RelayProblem.h"                  // for RelayProblem
#include "SparseBlockMatrix.h"             // for SBMfree, SparseBlockStruct...
#include "numerics_verbose.h"              // for CHECK_IO

GenericMechanicalProblem * genericMechanicalProblem_new()
{
  GenericMechanicalProblem * paux = (GenericMechanicalProblem *)malloc(sizeof(GenericMechanicalProblem));
  paux->firstListElem = 0;
  paux->lastListElem = 0;
  paux->size = 0;
  paux->maxLocalSize = 0;
  return paux;
}

void genericMechanicalProblem_free(GenericMechanicalProblem * pGMP, unsigned int level)
{
  if(!pGMP)
    return;
  while(pGMP->lastListElem)
  {
    listNumericsProblem * pElem = pGMP->lastListElem;
    free(pElem->q);
    switch(pElem->type)
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
    case SICONOS_NUMERICS_PROBLEM_FC2D:
    case SICONOS_NUMERICS_PROBLEM_FC3D:
    {
      free(((FrictionContactProblem*)(pElem->problem))->M);
      free(((FrictionContactProblem*)(pElem->problem))->mu);
      break;
    }
    default:
      printf("Numerics : genericMechanicalProblem_free case %d not managed.\n", pElem->type);
    }

    free(pElem->problem);
    pGMP->lastListElem = pElem->prevProblem;
    free(pElem);
  }
  if(level & NUMERICS_GMP_FREE_MATRIX)
  {
    assert(pGMP->M);
    NM_types storageType = pGMP->M->storageType;
    if(storageType == NM_DENSE)
      free(pGMP->M->matrix0);
    else
      SBMfree(pGMP->M->matrix1, NUMERICS_SBM_FREE_BLOCK | NUMERICS_SBM_FREE_SBM);
    free(pGMP->q);
    free(pGMP->M);
  }

  if(level & NUMERICS_GMP_FREE_GMP)
    free(pGMP);
}
void * gmp_add(GenericMechanicalProblem * pGMP, int problemType, int size)
{
  listNumericsProblem * newProblem = (listNumericsProblem*) malloc(sizeof(listNumericsProblem));
  newProblem->nextProblem = 0;
  newProblem->type = problemType;
  newProblem->size = size;
  newProblem->error = 0;
  pGMP->size += size;
  if(size > pGMP->maxLocalSize)
    pGMP->maxLocalSize = size;
  if(!pGMP->lastListElem)
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
  switch(problemType)
  {
  case(SICONOS_NUMERICS_PROBLEM_LCP):
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
  case(SICONOS_NUMERICS_PROBLEM_RELAY):
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
  case(SICONOS_NUMERICS_PROBLEM_EQUALITY):
  {
    newProblem->problem = NULL;
    newProblem->q = (double*) malloc(size * sizeof(double));;
    break;
  }
  case(SICONOS_NUMERICS_PROBLEM_FC3D):
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
  case(SICONOS_NUMERICS_PROBLEM_FC2D):
  {
    newProblem->problem = (void *) malloc(sizeof(FrictionContactProblem));
    FrictionContactProblem* pFC2D = (FrictionContactProblem*) newProblem->problem;
    pFC2D->mu = (double*) malloc(sizeof(double));
    pFC2D->M = NM_new();
    pFC2D->M->storageType = 0; /*Local prb is dense*/
    pFC2D->M->size0 = size;
    pFC2D->M->size1 = size;
    pFC2D->numberOfContacts = 1;
    pFC2D->q = (double*) malloc(size * sizeof(double));
    pFC2D->dimension = 3;
    newProblem->q = pFC2D->q;
    break;
  }
  default:
    printf("GenericMechanicalProblem.h gmp_add : problemType unknown: %d . \n", problemType);
    exit(EXIT_FAILURE);
  }
  return  newProblem->problem;
}


void genericMechanicalProblem_display(GenericMechanicalProblem * pGMP)
{
  listNumericsProblem * pElem = pGMP->firstListElem;
  int ii;
  printf("\nBEGIN Display a GenericMechanicalProblem(Numerics):\n");

  while(pElem)
  {
    printf("-->An sub-problem %s.\n", ns_problem_id_to_name(pElem->type));
    pElem = pElem->nextProblem;
  }
  printf("The sparce block matrice is :\n");
  NM_display(pGMP->M);
  printf("The q vector is :\n");
  for(ii = 0; ii < pGMP->size; ii++)
    printf("%e ", pGMP->q[ii]);

  //SBM_print(pGMP->M->matrix1);
  printf("\nEND Display a GenericMechanicalProblem:\n");
}

void genericMechanicalProblem_printInFile(GenericMechanicalProblem*  pGMP, FILE* file)
{
  listNumericsProblem * curProblem = pGMP->firstListElem;
  /*Print M*/
  NM_write_in_file(pGMP->M, file);
  fprintf(file, "\n");
  /*Print Q*/
  for(int ii = 0; ii < pGMP->size; ii++)
    fprintf(file, "%e\n", pGMP->q[ii]);
  fprintf(file, "\n");
  /*Print the type and options (mu)*/
  while(curProblem)
  {
    fprintf(file, "%d\n", curProblem->type);
    if(curProblem->type == SICONOS_NUMERICS_PROBLEM_FC3D)
      fprintf(file, "%e\n", ((FrictionContactProblem*)curProblem->problem)->mu[0]);
    curProblem = curProblem->nextProblem;
  }
}

GenericMechanicalProblem * genericMechanical_newFromFile(FILE* file)
{
  GenericMechanicalProblem* problem = genericMechanicalProblem_new();
  size_t nsubProb = 0;
  int prbType = 0;
  int i, posInX, localSize;
  void * prb;

  problem->M = NM_new_from_file(file);
  SparseBlockStructuredMatrix* m = problem->M->matrix1;

  problem->q = (double *) malloc(problem->M->size1 * sizeof(double));
  for(i = 0; i < problem->M->size1; i++)
  {
    CHECK_IO(fscanf(file, "%lf ", problem->q + i));
  }
  nsubProb = m->filled1 - 1;
  posInX = 0;
  for(size_t ii = 0; ii < nsubProb; ii++)
  {
    if(ii)
      posInX = m->blocksize0[ii - 1];
    localSize = m->blocksize0[ii] - posInX;
    CHECK_IO(fscanf(file, "%d\n", &prbType));
    prb = gmp_add(problem, prbType, localSize);
    if(prbType == SICONOS_NUMERICS_PROBLEM_FC3D)
    {
      CHECK_IO(fscanf(file, "%lf ", ((FrictionContactProblem*)prb)->mu));
    }
  }



#ifdef GMP_DEBUG
  genericMechanicalProblem_display(pGMP);
#endif
  return problem;
}

GenericMechanicalProblem * genericMechanical_new_from_filename(const char* filename)
{
  GenericMechanicalProblem * problem = NULL;
  FILE * file = fopen(filename, "r");
  if(file == NULL)
  {
    printf("Error! Could not open filename %s\n", filename);
    exit(EXIT_FAILURE);
  }

  problem = genericMechanical_newFromFile(file);
  fclose(file);
  return problem;
}

/** return nonsmooth problem formulation name, from its id number. */
const char * ns_problem_id_to_name(enum SICONOS_NUMERICS_PROBLEM_TYPE id)
{
  switch(id)
  {
  case(SICONOS_NUMERICS_PROBLEM_LCP):
  {
    return "LCP";
  }
  case(SICONOS_NUMERICS_PROBLEM_MLCP):
  {
    return "MLCP";
  }
  case(SICONOS_NUMERICS_PROBLEM_NCP):
  {
    return "NCP";
  }
  case(SICONOS_NUMERICS_PROBLEM_MCP):
  {
    return "MCP";
  }
  case(SICONOS_NUMERICS_PROBLEM_EQUALITY):
  {
    return "EQUALITY";
  }
  case(SICONOS_NUMERICS_PROBLEM_FC2D):
  {
    return "FC2D";
  }
  case(SICONOS_NUMERICS_PROBLEM_FC3D):
  {
    return "FC3D";
  }
  case(SICONOS_NUMERICS_PROBLEM_VI):
  {
    return "VI";
  }
  case(SICONOS_NUMERICS_PROBLEM_AVI):
  {
    return "AVI";
  }
  case(SICONOS_NUMERICS_PROBLEM_RELAY):
  {
    return "RELAY";
  }
  default:
    printf("Numerics:ns_problem_id_to_name, id unknown : %d \n", id);
    return NULL;
  }
}


