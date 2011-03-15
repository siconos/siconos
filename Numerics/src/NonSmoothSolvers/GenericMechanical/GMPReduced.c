#include "NonSmoothDrivers.h"
#include "pinv.h"
#include "LA.h"
//#define GMP_DEBUG_REDUCED
//#define GMP_DEBUG_GMPREDUCED_SOLVE

void _GMPReducedEquality(GenericMechanicalProblem* pInProblem, double * reducedProb, double * Qreduced, int * Me_size, int* Mi_size);
void _GMPReducedGetSizes(GenericMechanicalProblem* pInProblem, int * Me_size, int* Mi_size);


void printDenseMatrice(char* name, FILE * titi, double * m, int N, int M)
{
  if (titi)
  {
    fprintf(titi, "%s=[ \n", name);
    for (int i = 0; i < N; i++)
    {
      fprintf(titi, "[");
      for (int j = 0; j < M; j++)
      {
        fprintf(titi, "%e\t  ", m[i + j * N]);
      }
      fprintf(titi, "];\n");
    }
    fprintf(titi, "];\n");
  }
  else
  {
    printf("%s=[ \n", name);
    for (int i = 0; i < N; i++)
    {
      printf("[");
      for (int j = 0; j < M; j++)
      {
        printf("%e\t  ", m[i + j * N]);
      }
      printf("];\n");
    }
    printf("];\n");
  }
}

void GMPReducedSolToSol(GenericMechanicalProblem* pInProblem, double * reaction, double * velocity,
                        double * Re, double * Rreduced, double * Vreduced)
{
  listNumericsProblem * curProblem = 0;
  //int curRowE=0;
  //int curRowI=0;
  curProblem =  pInProblem->firstListElem;
  int curSize = 0;
  while (curProblem)
  {
    curSize = curProblem->size;
    switch (curProblem->type)
    {
    case SICONOS_NUMERICS_PROBLEM_EQUALITY:
    {
      memcpy(reaction, Re, curSize * sizeof(double));
      for (int i = 0; i < curSize; i++)
        velocity[i] = 0.0;
      Re += curSize;
      break;
    }
    case SICONOS_NUMERICS_PROBLEM_LCP:
    case SICONOS_NUMERICS_PROBLEM_FC3D:
    {
      memcpy(reaction, Rreduced, curSize * sizeof(double));
      memcpy(velocity, Vreduced, curSize * sizeof(double));
      Rreduced += curSize;
      Vreduced += curSize;
      break;
    }
    default:
      printf("GMPReducedSolToSol Numerics : genericMechanicalProblem_GS unknown problem type %d.\n", curProblem->type);

    }
    reaction += curSize;
    velocity += curSize;
    curProblem = curProblem->nextProblem;
  }
}
void _GMPReducedGetSizes(GenericMechanicalProblem* pInProblem, int * Me_size, int* Mi_size)
{
  listNumericsProblem * curProblem = 0;
  (* Me_size) = 0;
  (* Mi_size) = 0;
  curProblem =  pInProblem->firstListElem;
  while (curProblem)
  {
    if (curProblem->type == SICONOS_NUMERICS_PROBLEM_EQUALITY)
    {
      (*Me_size) += curProblem->size;;
    }
    else
    {
      (*Mi_size) += curProblem->size;
    }
    curProblem = curProblem->nextProblem;
  }
}

/*mem loc done */
void buildReducedGMP(GenericMechanicalProblem* pInProblem, double * Me, double * Mi, double * Qe, double * Qi, int * Me_Size, int* Mi_Size)
{
  assert(pInProblem->M->storageType);
  //#ifdef TYTYFCRR
  SparseBlockStructuredMatrix* m = pInProblem->M->matrix1;
#ifdef GMP_DEBUG_REDUCED
  FILE * titi  = fopen("buildReducedGMP_input.txt", "w");
  printInFileSBMForScilab(m, titi);
  fclose(titi);
#endif
  int currentRowNumber = 0;
  int  posInX = 0;
  int curSize = 0;
  //  int *newIndexOfBlockI;
  // int NbCol=pInProblem->size;
  int nbBlockCol = m->blocknumber1;
  int * newIndexOfCol = (int*) malloc(nbBlockCol * sizeof(int));

  /*Me building*/
  int MeRow = 0;
  int MiRow = 0;

  /**size of Me */
  listNumericsProblem * curProblem = 0;
  int nbBlockRowE = 0;
  int nbBlockRowI = 0;
  int numBlockRow = 0;
  curProblem =  pInProblem->firstListElem;
  while (curProblem)
  {
    if (numBlockRow)
      curSize = m->blocksize0[numBlockRow] - m->blocksize0[numBlockRow - 1];
    else
      curSize = m->blocksize0[numBlockRow];

    if (curProblem->type == SICONOS_NUMERICS_PROBLEM_EQUALITY)
    {
      nbBlockRowE++;
      MeRow += curSize;
    }
    else
    {
      nbBlockRowI++;
      MiRow += curSize;
    }
    curProblem = curProblem->nextProblem;
    numBlockRow++;
  }
  numBlockRow = 0;
  int numRowE = 0;
  int numRowI = 0;
  int numRow = 0;
  curProblem =  pInProblem->firstListElem;
  while (curProblem)
  {
    if (curProblem->type == SICONOS_NUMERICS_PROBLEM_EQUALITY)
    {
      newIndexOfCol[numRow] = numRowE;
      numRowE++;
    }
    else
    {
      newIndexOfCol[numRow] = numRowI + nbBlockRowE;
      numRowI++;
    }
    numRow++;
    curProblem = curProblem->nextProblem;
  }
#ifdef GMP_DEBUG_REDUCED
  printf("buildReducedGMP nb of block of eq=%i. nb of iq=%i\n", numRowE, numRowI);
#endif
  /*building of the permutation matrices*/
  SparseBlockStructuredMatrix Maux;
  ColPermutationSBM(newIndexOfCol, m, &Maux);
  SparseBlockStructuredMatrix Morder;
  RowPermutationSBM(newIndexOfCol, &Maux, &Morder);
  SBMfree(&Maux, 0);
  /*
    get the permutation indices of col (and row).

   */
  curProblem =  pInProblem->firstListElem;




  /**mem alloc for Me and Mi*/
  //int nbCol=MeRow+MiRow;
  *Me_Size = MeRow;
  *Mi_Size = MiRow;
  /*  Me=(double *) malloc(MeRow*nbCol*sizeof(double));
  Qe=(double *) malloc(MeRow*sizeof(double));
  Mi=(double *) malloc(MiRow*nbCol*sizeof(double));
  Qi=(double *) malloc(MiRow*sizeof(double));*/
  /** copi in Me*/
  int curPos = 0;
  for (int numBlockRow = 0; numBlockRow < nbBlockRowE; numBlockRow++)
  {
    SBMRowToDense(&Morder, numBlockRow, Me, curPos, MeRow);
    curPos = Morder.blocksize1[numBlockRow];
  }
  curPos = 0;
  int firtMiLine = 0;
  if (numBlockRow < nbBlockRowE + nbBlockRowI)
    firtMiLine = Morder.blocksize1[nbBlockRowE];

  for (int numBlockRow = nbBlockRowE; numBlockRow < nbBlockRowE + nbBlockRowI; numBlockRow++)
  {
    curPos = Morder.blocksize1[numBlockRow] - firtMiLine;
    SBMRowToDense(&Morder, numBlockRow, Mi, curPos, MiRow);
  }
  SBMfree(&Morder, 0);

  curProblem =  pInProblem->firstListElem;
  currentRowNumber = 0;
  posInX = 0;
  int curBlock = 0;
  int curPosIq = 0;
  int curPosEq = 0;
  double *curQ = pInProblem->q;
  double *curQe = Qe;
  double *curQi = Qi;
  curBlock = 0;
  while (curProblem)
  {
    if (curBlock)
    {
      curSize = m->blocksize0[curBlock] - m->blocksize0[curBlock - 1];
    }
    else
    {
      curSize = m->blocksize0[curBlock];
    }

    switch (curProblem->type)
    {
    case SICONOS_NUMERICS_PROBLEM_EQUALITY:
    {
      /** copy the current line block in Me*/
      memcpy(curQe, curQ, curSize * sizeof(double));
      curPosEq += curSize;
      curQe += curSize;
      break;
    }
    case SICONOS_NUMERICS_PROBLEM_LCP:
    case SICONOS_NUMERICS_PROBLEM_FC3D:
    {
      memcpy(curQi, curQ, curSize * sizeof(double));
      curPosIq += curSize;
      curQi += curSize;
      break;
    }
    default:
      printf("GMPReduced  buildReducedGMP: problemType unknown: %d . \n", curProblem->type);
    }
    curProblem = curProblem->nextProblem;
    curQ += curSize;
    curBlock++;
  }
#ifdef GMP_DEBUG_REDUCED
  printf("\\The Me matrix is:\n");
  printf("Me=[ \n");
  for (int i = 0; i < MeRow; i++)
  {
    printf("[");
    for (int j = 0; j < nbCol; j++)
    {
      printf("%e\t  ", Me[i + j * MeRow]);
    }
    printf("];\n");
  }
  printf("];\n");
  printf("Qe= [ \n");
  for (int i = 0; i < MeRow; i++)
    printf("%e\t  ", Qe[i]);
  printf("];\n");
  printf("\\The Mi matrix is:\n");
  printf("Mi=[ \n");
  for (int i = 0; i < MiRow; i++)
  {
    printf("[");
    for (int j = 0; j < nbCol; j++)
    {
      printf("%e\t  ", Mi[i + j * MiRow]);
    }
    printf("];\n");
  }
  printf("];\n");
  printf("Qi= [ \n");
  for (int i = 0; i < MiRow; i++)
    printf("%e\t  ", Qi[i]);
  printf("];\n");
#endif
  free(newIndexOfCol);
  //#endif
}
/*
 * The equalities are assamblate in one block.
 *
 *0=(Me_1 Me_2)(Re Ri)' + Qe
 *Vi=(Mi_1 Mi_2)(Re Ri)' + Qi
 *
 *and GS.
 */
void GMPReducedEqualitySolve(GenericMechanicalProblem* pInProblem, double *reaction , double *velocity, int * info, SolverOptions* options)
{

  SparseBlockStructuredMatrix* m = pInProblem->M->matrix1;
  int nbRow = m->blocksize0[m->blocknumber0 - 1];
  int nbCol = m->blocksize1[m->blocknumber1 - 1];

  int Me_size;
  int Mi_size;
  double * reducedProb = (double *)malloc(nbRow * nbCol * sizeof(double));
  double * Qreduced = (double *)malloc(nbRow * sizeof(double));
  _GMPReducedEquality(pInProblem, reducedProb, Qreduced, &Me_size, &Mi_size);

  if (Me_size == 0)
  {
    genericMechanicalProblem_GS(pInProblem, reaction, velocity, info, options);
    free(reducedProb);
    free(Qreduced);
    return;
  }
  listNumericsProblem * curProblem = 0;
  GenericMechanicalProblem * _pnumerics_GMP = buildEmptyGenericMechanicalProblem();
  curProblem =  pInProblem->firstListElem;
  if (Me_size)
    addProblem(_pnumerics_GMP, SICONOS_NUMERICS_PROBLEM_EQUALITY, Me_size);
  while (curProblem)
  {
    switch (curProblem->type)
    {
    case SICONOS_NUMERICS_PROBLEM_EQUALITY:
    {
      break;
    }
    case SICONOS_NUMERICS_PROBLEM_LCP:
    {
      addProblem(_pnumerics_GMP, curProblem->type, curProblem->size);
      break;
    }
    case SICONOS_NUMERICS_PROBLEM_FC3D:
    {
      FrictionContactProblem* pFC3D = (FrictionContactProblem*)addProblem(_pnumerics_GMP, curProblem->type, curProblem->size);
      *(pFC3D->mu) = *(((FrictionContactProblem*)curProblem->problem)->mu);
      break;
    }
    default:
      printf("GMPReduced  buildReducedGMP: problemType unknown: %d . \n", curProblem->type);
    }
    curProblem = curProblem->nextProblem;
  }

#ifdef GMP_DEBUG_GMPREDUCED_SOLVE
  printDenseMatrice("newPrb", titi, reducedProb, nbRow, nbCol);
  printDenseMatrice("newQ", titi, Qreduced, nbRow, 1);
#endif
  NumericsMatrix numM;
  numM.storageType = 0;
  numM.matrix0 = reducedProb;
  numM.matrix1 = 0;
  numM.size0 = nbRow;
  numM.size1 = nbCol;
  _pnumerics_GMP->M = &numM;
  _pnumerics_GMP->q = Qreduced;
  double *Rreduced = (double *) calloc(nbCol, sizeof(double));
  double *Vreduced = (double *) calloc(nbRow, sizeof(double));
  genericMechanicalProblem_GS(_pnumerics_GMP, Rreduced, Vreduced, info, options);
#ifdef GMP_DEBUG_GMPREDUCED_SOLVE
  if (*info)
  {
    printf("\nGMPREduced2 failed!\n");
  }
  else
  {
    printf("\nGMPREduced2 succed!\n");
    printDenseMatrice("R", titi, Rreduced, nbRow, 1);
    printDenseMatrice("V", titi, Vreduced, nbRow, 1);
  }
#endif
  if (*info)
    goto END_GMP2;
  GMPReducedSolToSol(pInProblem, reaction, velocity, Rreduced, Rreduced + Me_size, Vreduced + Me_size);
#ifdef GMP_DEBUG_GMPREDUCED_SOLVE
  printDenseMatrice("R2", titi, reaction, nbRow, 1);
  printDenseMatrice("V2", titi, velocity, nbRow, 1);
#endif
  double err;
  int tolViolate = GenericMechanical_compute_error(pInProblem, reaction, velocity, options->dparam[0], options, &err);
  if (tolViolate)
  {
    printf("GMPReduced2, warnning, reduced problem solved, but error of intial probleme violated tol = %e, err= %e\n", options->dparam[0], err);
  }


END_GMP2:
  ;
#ifdef GMP_DEBUG_GMPREDUCED_SOLVE
  fclose(titi);
#endif
  free(Rreduced);
  free(Vreduced);
  freeGenericMechanicalProblem(_pnumerics_GMP, NUMERICS_GMP_FREE_GMP);
  free(Qreduced);
  free(reducedProb);
}

/*
 * The equalities are eliminated.
 *
 *0=(Me_1 Me_2)(Re Ri)' + Qe
 *Vi=(Mi_1 Mi_2)(Re Ri)' + Qi
 *
 *Re=-Me_1^{-1}(Me_2Ri+Qe)
 *
 *Vi=(Mi_2-Mi_1 Me_1^{-1} Me_2)Ri+Qi-Mi1 Me_1^{-1} Qe
 *
 */
void GMPReducedSolve(GenericMechanicalProblem* pInProblem, double *reaction , double *velocity, int * info, SolverOptions* options)
{

  SparseBlockStructuredMatrix* m = pInProblem->M->matrix1;
  int nbRow = m->blocksize0[m->blocknumber0 - 1];
  int nbCol = m->blocksize1[m->blocknumber1 - 1];
  double *Me = (double *) malloc(nbRow * nbCol * sizeof(double));
  double *Qe = (double *) malloc(nbRow * sizeof(double));
  double *Mi = (double *) malloc(nbRow * nbCol * sizeof(double));
  double *Qi = (double *) malloc(nbRow * sizeof(double));
  int Me_size;
  int Mi_size;
  buildReducedGMP(pInProblem, Me, Mi, Qe, Qi, &Me_size, &Mi_size);
  if ((Me_size == 0 || Mi_size == 0))
  {
    genericMechanicalProblem_GS(pInProblem, reaction, velocity, info, options);
    free(Me);
    free(Qe);
    free(Mi);
    free(Qi);
    return;
  }
  double * pseduInvMe1 = (double *)malloc(Me_size * Me_size * sizeof(double));
  if (Mi_size == 0)
  {
    DGEMV(LA_NOTRANS, Me_size, Me_size, -1.0, pseduInvMe1, Me_size, pInProblem->q , 1, 0.0, reaction, 1);
    for (int i = 0; i < Me_size; i++)
      velocity[i] = 0.0;

    free(Me);
    free(Qe);
    free(Mi);
    free(Qi);
    free(pseduInvMe1);
    return;

  }
  memcpy(pseduInvMe1, Me, Me_size * Me_size * sizeof(double));
  pinv(pseduInvMe1, Me_size, Me_size, 1e-16);
  double *Mi2 = Mi + Mi_size * Me_size;
  double *Mi1 = Mi;
  double *Me2 = Me + Me_size * Me_size;
#ifdef GMP_DEBUG_GMPREDUCED_SOLVE
  double *Me1 = Me;
  FILE * titi  = fopen("buildReducedGMP_output.txt", "w");
  printf("GMPReducedsolve\n");
  printDenseMatrice("Me1", titi, Me1, Me_size, Me_size);
  printDenseMatrice("Me2", titi, Me2, Me_size, Mi_size);
  printDenseMatrice("Mi1", titi, Mi1, Mi_size, Me_size);
  printDenseMatrice("Mi2", titi, Mi2, Mi_size, Mi_size);
  printDenseMatrice("Qe", titi, Qe, Me_size, 1);
  printDenseMatrice("Qi", titi, Qi, Mi_size, 1);
  printDenseMatrice("Me1inv", titi, pseduInvMe1, Me_size, Me_size);
#endif


  double * reducedProb = (double *)malloc(Mi_size * Mi_size * sizeof(double));
  memcpy(reducedProb, Mi2, Mi_size * Mi_size * sizeof(double));

  double * Mi1pseduInvMe1 = (double *)malloc(Mi_size * Me_size * sizeof(double));
  DGEMM(LA_NOTRANS, LA_NOTRANS, Mi_size, Me_size, Me_size, -1.0, Mi1, Mi_size, pseduInvMe1, Me_size, 0.0, Mi1pseduInvMe1, Mi_size);
#ifdef GMP_DEBUG_GMPREDUCED_SOLVE
  printDenseMatrice("minusMi1pseduInvMe1", titi, Mi1pseduInvMe1, Mi_size, Me_size);
  fprintf(titi, "_minusMi1pseduInvMe1=-Mi1*Me1inv;\n");
#endif
  DGEMV(LA_NOTRANS, Mi_size, Me_size, 1.0, Mi1pseduInvMe1, Mi_size, Qe, 1, 1.0, Qi, 1);
#ifdef GMP_DEBUG_GMPREDUCED_SOLVE
  printDenseMatrice("newQi", titi, Qi, Mi_size, 1);
  fprintf(titi, "_newQi=Qi+_minusMi1pseduInvMe1*Qe;\n");
#endif
  DGEMM(LA_NOTRANS, LA_NOTRANS, Mi_size, Mi_size, Me_size, 1.0, Mi1pseduInvMe1, Mi_size, Me2, Me_size, 1.0, reducedProb, Mi_size);
#ifdef GMP_DEBUG_GMPREDUCED_SOLVE
  printDenseMatrice("W", titi, reducedProb, Mi_size, Mi_size);
  fprintf(titi, "_W=Mi2+_minusMi1pseduInvMe1*Me2;\n");

#endif
  listNumericsProblem * curProblem = 0;
  GenericMechanicalProblem * _pnumerics_GMP = buildEmptyGenericMechanicalProblem();
  curProblem =  pInProblem->firstListElem;
  while (curProblem)
  {
    switch (curProblem->type)
    {
    case SICONOS_NUMERICS_PROBLEM_EQUALITY:
    {
      break;
    }
    case SICONOS_NUMERICS_PROBLEM_LCP:
    {
      addProblem(_pnumerics_GMP, curProblem->type, curProblem->size);
      break;
    }
    case SICONOS_NUMERICS_PROBLEM_FC3D:
    {
      FrictionContactProblem* pFC3D = (FrictionContactProblem*)addProblem(_pnumerics_GMP, curProblem->type, curProblem->size);
      *(pFC3D->mu) = *(((FrictionContactProblem*)curProblem->problem)->mu);
      break;
    }
    default:
      printf("GMPReduced  buildReducedGMP: problemType unknown: %d . \n", curProblem->type);
    }
    curProblem = curProblem->nextProblem;
  }
  NumericsMatrix numM;
  numM.storageType = 0;
  numM.matrix0 = reducedProb;
  numM.matrix1 = 0;
  numM.size0 = Mi_size;
  numM.size1 = Mi_size;
  _pnumerics_GMP->M = &numM;
  _pnumerics_GMP->q = Qi;
  double *Rreduced = (double *) malloc(Mi_size * sizeof(double));
  double *Vreduced = (double *) malloc(Mi_size * sizeof(double));
  genericMechanicalProblem_GS(_pnumerics_GMP, Rreduced, Vreduced, info, options);
#ifdef GMP_DEBUG_GMPREDUCED_SOLVE
  if (*info)
  {
    printf("\nGMPREduced failed!\n");
  }
  else
  {
    printf("\nGMPREduced succed!\n");
    printDenseMatrice("Ri", titi, Rreduced, Mi_size, 1);
    printDenseMatrice("Vi", titi, Vreduced, Mi_size, 1);
  }
#endif
  if (*info)
    goto END_GMP1;
  /*Re computation*/
  double * Re = (double*)malloc(Me_size * sizeof(double));
  double * Rbuf = (double*)malloc(Me_size * sizeof(double));
  memcpy(Rbuf, Qe, Me_size * sizeof(double));
  DGEMV(LA_NOTRANS, Me_size, Mi_size, 1.0, Me2, Me_size, Rreduced, 1, 1.0, Rbuf, 1);
  DGEMV(LA_NOTRANS, Me_size, Me_size, -1.0, pseduInvMe1, Me_size, Rbuf, 1, 0.0, Re, 1);
#ifdef GMP_DEBUG_GMPREDUCED_SOLVE
  fprintf(titi, "_Re=-Me1inv*(Me2*Ri+Qe);\n");
  printDenseMatrice("Re", titi, Re, Me_size, 1);
#endif
  GMPReducedSolToSol(pInProblem, reaction, velocity, Re, Rreduced, Vreduced);
  double err;
  int tolViolate = GenericMechanical_compute_error(pInProblem, reaction, velocity, options->dparam[0], options, &err);
  if (tolViolate)
  {
    printf("GMPReduced, warnning, reduced problem solved, but error of intial probleme violated tol = %e, err= %e\n", options->dparam[0], err);
  }
  free(Re);
  free(Rbuf);
END_GMP1:
  ;
#ifdef GMP_DEBUG_GMPREDUCED_SOLVE
  fclose(titi);
#endif
  free(Rreduced);
  free(Vreduced);
  freeGenericMechanicalProblem(_pnumerics_GMP, NUMERICS_GMP_FREE_GMP);
  free(Me);
  free(Mi);
  free(Qe);
  free(Qi);
  free(pseduInvMe1);
  free(reducedProb);
  free(Mi1pseduInvMe1);
  //  GenericMechanicalProblem GMPOutProblem;
  //  SparseBlockStructuredMatrix mOut;

}

void _GMPReducedEquality(GenericMechanicalProblem* pInProblem, double * reducedProb, double * Qreduced, int * Me_size, int* Mi_size)
{

  SparseBlockStructuredMatrix* m = pInProblem->M->matrix1;
  int nbRow = m->blocksize0[m->blocknumber0 - 1];
  int nbCol = m->blocksize1[m->blocknumber1 - 1];

  _GMPReducedGetSizes(pInProblem, Me_size, Mi_size);

  double *Me = (*Me_size) ? (double *) malloc((*Me_size) * nbCol * sizeof(double)) : 0;
  double *Mi = (*Mi_size) ? (double *) malloc((*Mi_size) * nbCol * sizeof(double)) : 0;
  double *Qi = (double *) malloc(nbRow * sizeof(double));
  buildReducedGMP(pInProblem, Me, Mi, Qreduced, Qi, Me_size, Mi_size);

#ifdef GMP_DEBUG_GMPREDUCED_SOLVE
  double *Me1 = Me;
  double *Me2 = Me + (*Me_size) * (*Me_size);
  double *Mi1 = Mi;
  double *Mi2 = Mi + (*Mi_size) * (*Me_size);
  FILE * titi  = fopen("buildReduced2GMP_output.txt", "w");
  printf("GMP2Reducedsolve\n");
  printDenseMatrice("Me1", titi, Me1, Me_size, Me_size);
  printDenseMatrice("Me2", titi, Me2, Me_size, Mi_size);
  printDenseMatrice("Mi1", titi, Mi1, Mi_size, Me_size);
  printDenseMatrice("Mi2", titi, Mi2, Mi_size, Mi_size);
  printDenseMatrice("Qe", titi, Qreduced, Me_size, 1);
  printDenseMatrice("Qi", titi, Qi, Mi_size, 1);
#endif
  for (int numCol = 0; numCol < nbCol; numCol++)
  {
    if (*Me_size)
      memcpy(reducedProb + numCol * nbRow, Me + numCol * (*Me_size), (*Me_size)*sizeof(double));
    if (*Mi_size)
      memcpy(reducedProb + numCol * nbRow + (*Me_size), Mi + numCol * (*Mi_size), (*Mi_size)*sizeof(double));
  }
  if (*Mi_size)
    memcpy(Qreduced + (*Me_size), Qi, (*Mi_size)*sizeof(double));
  free(Me);
  free(Mi);
  free(Qi);
}

void GMPasMLCP(GenericMechanicalProblem* pInProblem, double *reaction , double *velocity, int* info, SolverOptions* options)
{

  /*First, we don't manage FC3D.*/
  listNumericsProblem * curProblem = 0;
  curProblem =  pInProblem->firstListElem;
  while (curProblem)
  {
    switch (curProblem->type)
    {
    case SICONOS_NUMERICS_PROBLEM_EQUALITY:
    case SICONOS_NUMERICS_PROBLEM_LCP:
      break;
    case SICONOS_NUMERICS_PROBLEM_FC3D:
    {
      printf("GMPasMLCP Numerics ERROR: GMPasMLCP doesn't deal with FC3D.\n");
      *info = 1;
      return;
      break;
    }
    default:
      printf("GMPasMLCP Numerics : genericMechanicalProblem_GS unknown problem type %d.\n", curProblem->type);
    }
    curProblem = curProblem->nextProblem;
  }
  int Me_size;
  int Mi_size;

  SparseBlockStructuredMatrix* m = pInProblem->M->matrix1;
  int nbRow = m->blocksize0[m->blocknumber0 - 1];
  int nbCol = m->blocksize1[m->blocknumber1 - 1];

  double * reducedProb = (double *)malloc(nbRow * nbCol * sizeof(double));
  double * Qreduced = (double *)malloc(nbRow * sizeof(double));
  _GMPReducedEquality(pInProblem, reducedProb, Qreduced, &Me_size, &Mi_size);

  if (!Me_size)
  {
    /*it is a lcp.*/
    LinearComplementarityProblem aLCP;
    SolverOptions aLcpOptions;
    NumericsMatrix M;
    M.storageType = 0;
    M.size0 = Mi_size;
    M.size1 = Mi_size;
    M.matrix0 = reducedProb;
    M.matrix1 = 0;
    aLCP.size = Mi_size;
    aLCP.q = Qreduced;
    aLCP.M = &M;
    linearComplementarity_setDefaultSolverOptions(&aLCP, &aLcpOptions, SICONOS_LCP_ENUM);
    lcp_enum_init(&aLCP, &aLcpOptions, 1);
    *info = linearComplementarity_driver(&aLCP, reaction, velocity, &aLcpOptions, 0);
    lcp_enum_reset(&aLCP, &aLcpOptions, 1);
    goto END_GMP3;
  }
  if (!Mi_size)
  {
    /*it is a linear system.*/
    LinearSystemProblem aLS;
    NumericsMatrix M;
    M.storageType = 0;
    M.size0 = Mi_size;
    M.size1 = Mi_size;
    M.matrix0 = reducedProb;
    M.matrix1 = 0;
    aLS.size = Me_size;
    aLS.q = Qreduced;
    aLS.M = &M;
    *info = LinearSystem_driver(&aLS, reaction, velocity, 0);
    goto END_GMP3;
  }
  /*it is a MLCP*/
  MixedLinearComplementarityProblem aMLCP;
  SolverOptions aMlcpOptions;
  aMLCP.n = Me_size;
  aMLCP.m = Mi_size;
  aMLCP.blocksLine = 0;
  aMLCP.blocksIsComp = 0;
  aMLCP.problemType = 0;
  aMLCP.A = 0;
  aMLCP.B = 0;
  aMLCP.C = 0;
  aMLCP.D = 0;
  aMLCP.a = 0;
  aMLCP.b = 0;
  aMLCP.q = Qreduced;
  NumericsMatrix M;
  M.storageType = 0;
  M.size0 = Mi_size + Me_size;
  M.size1 = Mi_size + Me_size;
  M.matrix0 = reducedProb;
  M.matrix1 = 0;
  aMLCP.M = &M;
  aMlcpOptions.solverId = SICONOS_MLCP_ENUM;
  mixedLinearComplementarity_setDefaultSolverOptions(&aMLCP, &aMlcpOptions);
  mlcp_driver_init(&aMLCP, &aMlcpOptions);
  aMlcpOptions.dparam[0] = options->dparam[0];
  *info = mlcp_driver(&aMLCP, reaction, velocity, &aMlcpOptions, 0);

  mlcp_driver_reset(&aMLCP, &aMlcpOptions);
  mixedLinearComplementarity_deleteDefaultSolverOptions(&aMLCP, &aMlcpOptions);
END_GMP3:
  ;

  free(reducedProb);
  free(Qreduced);

}
