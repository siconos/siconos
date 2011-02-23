#include "GenericMechanical_Solvers.h"
#include "pinv.h"
#include "LA.h"
#define GMP_DEBUG_REDUCED


void printDenseMatrice(char* name, FILE * titi, double * m, int N, int M)
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

/*mem loc done */
void buildReducedGMP(GenericMechanicalProblem* pInProblem, double * Me, double * Mi, double * Qe, double * Qi, int * Me_Size, int* Mi_Size)
{

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
  int nbCol = MeRow + MiRow;
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
  for (int numBlockRow = nbBlockRowE; numBlockRow < nbBlockRowE + nbBlockRowI; numBlockRow++)
  {
    SBMRowToDense(&Morder, numBlockRow, Mi, curPos, MiRow);
    curPos = Morder.blocksize1[numBlockRow];
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

  //#endif
}

#define GMP_DEBUG_GMPREDUCED_SOLVE
void GMPReducedsolve(GenericMechanicalProblem* pInProblem, double *reaction , double *velocity, SolverOptions* options)
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
  double * pseduInvMe1 = (double *)malloc(Me_size * Me_size * sizeof(double));
  memcpy(pseduInvMe1, Me, Me_size * Me_size * sizeof(double));
  pinv(pseduInvMe1, Me_size, Me_size, 1e-16);
  double *Me1 = Me;
  double *Me2 = Me + Me_size * Me_size;
  double *Mi1 = Mi;
  double *Mi2 = Mi + Mi_size * Me_size;
#ifdef GMP_DEBUG_GMPREDUCED_SOLVE
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
  printf("-Mi1*pseduInvMe1");
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
  fclose(titi);
#endif
  //  GenericMechanicalProblem GMPOutProblem;
  //  SparseBlockStructuredMatrix mOut;

}

void convertReducedSolToSol(double *pVin, double *pRin, double *pVout, double* pRout, double * matL, int* corresp)
{
}

void deleteReducedGMP(GenericMechanicalProblem* pGMP)
{
}
