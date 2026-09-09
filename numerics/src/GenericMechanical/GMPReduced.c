/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include "GMPReduced.h"

#include <assert.h>  // for assert

#include "safe_casts.h"
#ifndef __cplusplus
#include <stdbool.h>  // for true
#endif
#include <stdio.h>   // for printf, size_t, NULL
#include <stdlib.h>  // for free, malloc, calloc
#include <string.h>  // for memcpy

#include "FrictionContactProblem.h"        // IWYU pragma: keep
#include "GenericMechanicalProblem.h"      // for listNumericsProblem
#include "GenericMechanical_Solvers.h"     // for gmp_gauss_seidel, gmp...
#include "LCP_Solvers.h"                   // for lcp_enum_init, lcp_en...
#include "LinearComplementarityProblem.h"  // IWYU pragma: keep
#include "MLCP_Solvers.h"                  // for mixedLinearComplement...
#include "NonSmoothDrivers.h"              // for linearComplementarity...
#include "NumericsMatrix.h"                // for NumericsMatrix, NM_fill
#include "SiconosBlas.h"                   // for cblas_dgemv, CblasNoT...
#include "SolverOptions.h"                 // for SICONOS_NUMERICS_PROB...
#include "SparseBlockMatrix.h"             // for SparseBlockStructured...
#include "lcp_cst.h"                       // for SICONOS_LCP_ENUM
#include "mlcp_cst.h"                      // for SICONOS_MLCP_ENUM
#include "pinv.h"                          // for pinv

#ifdef GMP_DEBUG_REDUCED
static void printDenseMatrice(char *name, FILE *file, double *m, int N, int M) {
  if (file) {
    fprintf(file, "%s=[ \n", name);
    for (int i = 0; i < N; i++) {
      fprintf(file, "[");
      for (int j = 0; j < M; j++) {
        fprintf(file, "%e\t  ", m[i + j * N]);
      }
      fprintf(file, "];\n");
    }
    fprintf(file, "];\n");
  } else {
    printf("%s=[ \n", name);
    for (int i = 0; i < N; i++) {
      printf("[");
      for (int j = 0; j < M; j++) {
        printf("%e\t  ", m[i + j * N]);
      }
      printf("];\n");
    }
    printf("];\n");
  }
}
#endif

void gmp_reduced_convert_solution(GenericMechanicalProblem *pInProblem, double *reaction,
                                  double *velocity, double *Re, double *Rreduced,
                                  double *Vreduced) {
  listNumericsProblem *curProblem = 0;
  // int curRowE=0;
  // int curRowI=0;
  curProblem = pInProblem->firstListElem;
  size_t curSize = 0;
  while (curProblem) {
    curSize = to_size_t(curProblem->size);
    switch (curProblem->type) {
      case SICONOS_NUMERICS_PROBLEM_EQUALITY: {
        memcpy(reaction, Re, curSize * sizeof(double));
        for (size_t i = 0; i < curSize; i++) velocity[i] = 0.0;
        Re += curSize;
        break;
      }
      case SICONOS_NUMERICS_PROBLEM_LCP:
      case SICONOS_NUMERICS_PROBLEM_FC3D: {
        memcpy(reaction, Rreduced, curSize * sizeof(double));
        memcpy(velocity, Vreduced, curSize * sizeof(double));
        Rreduced += curSize;
        Vreduced += curSize;
        break;
      }
      default:
        printf(
            "gmp_reduced_convert_solution Numerics : gmp_gauss_seidel unknown problem type "
            "%d.\n",
            curProblem->type);
    }
    reaction += curSize;
    velocity += curSize;
    curProblem = curProblem->nextProblem;
  }
}
static void _GMPReducedGetSizes(GenericMechanicalProblem *pInProblem, size_t *Me_size,
                                size_t *Mi_size) {
  listNumericsProblem *curProblem = 0;
  (*Me_size) = 0;
  (*Mi_size) = 0;
  curProblem = pInProblem->firstListElem;
  while (curProblem) {
    if (curProblem->type == SICONOS_NUMERICS_PROBLEM_EQUALITY) {
      (*Me_size) += to_size_t(curProblem->size);
      ;
    } else {
      (*Mi_size) += to_size_t(curProblem->size);
    }
    curProblem = curProblem->nextProblem;
  }
}

/*mem loc done */
static void buildReducedGMP(GenericMechanicalProblem *pInProblem, double *Me, double *Mi,
                            double *Qe, double *Qi, size_t *Me_Size, size_t *Mi_Size) {
  assert(pInProblem->M->storageType);
  // #ifdef TYTYFCRR
  SparseBlockStructuredMatrix *m = pInProblem->M->matrix1;
#ifdef GMP_DEBUG_REDUCED
  FILE *file = fopen("buildReducedGMP_input.txt", "w");
  SBM_write_in_fileForScilab(m, file);
  fclose(file);
#endif
  size_t curSize = 0;
  //  int *newIndexOfBlockI;
  // int NbCol=pInProblem->size;
  size_t nbBlockCol = m->blocknumber1;
  size_t *newIndexOfCol = (size_t *)malloc(nbBlockCol * sizeof(size_t));

  /*Me building*/
  size_t MeRow = 0;
  size_t MiRow = 0;

  /**size of Me */
  listNumericsProblem *curProblem = 0;
  size_t nbBlockRowE = 0;
  size_t nbBlockRowI = 0;
  size_t numBlockRow = 0;
  curProblem = pInProblem->firstListElem;
  while (curProblem) {
    if (numBlockRow)
      curSize = m->blocksize0[numBlockRow] - m->blocksize0[numBlockRow - 1];
    else
      curSize = m->blocksize0[numBlockRow];

    if (curProblem->type == SICONOS_NUMERICS_PROBLEM_EQUALITY) {
      nbBlockRowE++;
      MeRow += curSize;
    } else {
      nbBlockRowI++;
      MiRow += curSize;
    }
    curProblem = curProblem->nextProblem;
    numBlockRow++;
  }
  numBlockRow = 0;
  size_t numRowE = 0;
  size_t numRowI = 0;
  size_t numRow = 0;
  curProblem = pInProblem->firstListElem;
  while (curProblem) {
    if (curProblem->type == SICONOS_NUMERICS_PROBLEM_EQUALITY) {
      newIndexOfCol[numRow] = numRowE;
      numRowE++;
    } else {
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
  SparseBlockStructuredMatrix *Maux = SBM_new();
  SBM_column_permutation(newIndexOfCol, m, Maux);
  SparseBlockStructuredMatrix *Morder = SBM_new();
  SBM_row_permutation(newIndexOfCol, Maux, Morder);
  // free Maux but not the blocks, since they are shared with m and Morder
  Maux = SBM_free(Maux, SBM_FREE_KEEP_BLOCKS);

  /*
    get the permutation indices of col (and row).

   */
  curProblem = pInProblem->firstListElem;

  /**mem alloc for Me and Mi*/
  // int nbCol=MeRow+MiRow;
  *Me_Size = MeRow;
  *Mi_Size = MiRow;
  /*  Me=(double *) malloc(MeRow*nbCol*sizeof(double));
  Qe=(double *) malloc(MeRow*sizeof(double));
  Mi=(double *) malloc(MiRow*nbCol*sizeof(double));
  Qi=(double *) malloc(MiRow*sizeof(double));*/
  /** copi in Me*/
  size_t curPos = 0;
  for (size_t numBlockRow = 0; numBlockRow < nbBlockRowE; numBlockRow++) {
    SBM_row_to_dense(Morder, numBlockRow, Me, curPos, MeRow);
    curPos = Morder->blocksize1[numBlockRow];
  }
  curPos = 0;
  size_t firtMiLine = 0;
  if (nbBlockRowI > 0) firtMiLine = Morder->blocksize1[nbBlockRowE];

  for (size_t numBlockRow = nbBlockRowE; numBlockRow < nbBlockRowE + nbBlockRowI;
       numBlockRow++) {
    curPos = Morder->blocksize1[numBlockRow] - firtMiLine;
    SBM_row_to_dense(Morder, numBlockRow, Mi, curPos, MiRow);
  }
  Morder = SBM_free(Morder, SBM_FREE_KEEP_BLOCKS);

  curProblem = pInProblem->firstListElem;
  int curBlock = 0;
  // int curPosIq = 0;
  // int curPosEq = 0;
  double *curQ = pInProblem->q;
  double *curQe = Qe;
  double *curQi = Qi;
  curBlock = 0;
  while (curProblem) {
    if (curBlock) {
      curSize = m->blocksize0[curBlock] - m->blocksize0[curBlock - 1];
    } else {
      curSize = m->blocksize0[curBlock];
    }

    switch (curProblem->type) {
      case SICONOS_NUMERICS_PROBLEM_EQUALITY: {
        /** copy the current line block in Me*/
        memcpy(curQe, curQ, curSize * sizeof(double));
        // curPosEq += curSize;
        curQe += curSize;
        break;
      }
      case SICONOS_NUMERICS_PROBLEM_LCP:
      case SICONOS_NUMERICS_PROBLEM_FC3D: {
        memcpy(curQi, curQ, curSize * sizeof(double));
        //        curPosIq += curSize;
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
  int nbCol = MeRow + MiRow;
  printf("\\The Me matrix is:\n");
  printf("Me=[ \n");
  for (int i = 0; i < MeRow; i++) {
    printf("[");
    for (int j = 0; j < nbCol; j++) {
      printf("%e\t  ", Me[i + j * MeRow]);
    }
    printf("];\n");
  }
  printf("];\n");
  printf("Qe= [ \n");
  for (int i = 0; i < MeRow; i++) printf("%e\t  ", Qe[i]);
  printf("];\n");
  printf("\\The Mi matrix is:\n");
  printf("Mi=[ \n");
  for (int i = 0; i < MiRow; i++) {
    printf("[");
    for (int j = 0; j < nbCol; j++) {
      printf("%e\t  ", Mi[i + j * MiRow]);
    }
    printf("];\n");
  }
  printf("];\n");
  printf("Qi= [ \n");
  for (int i = 0; i < MiRow; i++) printf("%e\t  ", Qi[i]);
  printf("];\n");
#endif
  free(newIndexOfCol);
  // #endif
}

static void _GMPReducedEquality(GenericMechanicalProblem *pInProblem, double *reducedProb,
                                double *Qreduced, size_t *Me_size, size_t *Mi_size) {
  SparseBlockStructuredMatrix *m = pInProblem->M->matrix1;
  size_t nbRow = m->blocksize0[m->blocknumber0 - 1];
  size_t nbCol = m->blocksize1[m->blocknumber1 - 1];

  _GMPReducedGetSizes(pInProblem, Me_size, Mi_size);
  if (*Me_size == 0) {
    memcpy(Qreduced, pInProblem->q, (*Mi_size) * sizeof(double));
    SBM_to_dense(m, reducedProb);
    return;
  }

  double *Me = (*Me_size) ? (double *)malloc((*Me_size) * nbCol * sizeof(double)) : 0;
  double *Mi = (*Mi_size) ? (double *)malloc((*Mi_size) * nbCol * sizeof(double)) : 0;
  double *Qi = (double *)malloc(nbRow * sizeof(double));
  buildReducedGMP(pInProblem, Me, Mi, Qreduced, Qi, Me_size, Mi_size);

#ifdef GMP_DEBUG_GMPREDUCED_SOLVE
  double *Me1 = Me;
  double *Me2 = Me + (*Me_size) * (*Me_size);
  double *Mi1 = Mi;
  double *Mi2 = Mi + (*Mi_size) * (*Me_size);
  FILE *file = fopen("buildReduced2GMP_output.txt", "w");
  printf("GMP2Reducedsolve\n");
  printDenseMatrice("Me1", file, Me1, *Me_size, *Me_size);
  printDenseMatrice("Me2", file, Me2, *Me_size, *Mi_size);
  printDenseMatrice("Mi1", file, Mi1, *Mi_size, *Me_size);
  printDenseMatrice("Mi2", file, Mi2, *Mi_size, *Mi_size);
  printDenseMatrice("Qe", file, Qreduced, *Me_size, 1);
  printDenseMatrice("Qi", file, Qi, *Mi_size, 1);
  fclose(file);
#endif
  for (size_t numCol = 0; numCol < nbCol; numCol++) {
    if (*Me_size)
      memcpy(reducedProb + numCol * nbRow, Me + numCol * (*Me_size),
             (*Me_size) * sizeof(double));
    if (*Mi_size)
      memcpy(reducedProb + numCol * nbRow + (*Me_size), Mi + numCol * (*Mi_size),
             (*Mi_size) * sizeof(double));
  }
  if (*Mi_size) memcpy(Qreduced + (*Me_size), Qi, (*Mi_size) * sizeof(double));
  free(Me);
  free(Mi);
  free(Qi);
}

/*
 * The equalities are assamblate in one block.
 *
 *0=(Me_1 Me_2)(Re Ri)' + Qe
 *Vi=(Mi_1 Mi_2)(Re Ri)' + Qi
 *
 *and GS.
 */
void gmp_reduced_equality_solve(GenericMechanicalProblem *pInProblem, double *reaction,
                                double *velocity, int *info, SolverOptions *options) {
  SparseBlockStructuredMatrix *m = pInProblem->M->matrix1;
  size_t nbRow = m->blocksize0[m->blocknumber0 - 1];
  size_t nbCol = m->blocksize1[m->blocknumber1 - 1];

  size_t Me_size;
  size_t Mi_size;
  double *reducedProb = (double *)malloc(nbRow * nbCol * sizeof(double));
  double *Qreduced = (double *)malloc(nbRow * sizeof(double));
  double *Rreduced = (double *)calloc(nbCol, sizeof(double));
  double *Vreduced = (double *)calloc(nbRow, sizeof(double));

  _GMPReducedEquality(pInProblem, reducedProb, Qreduced, &Me_size, &Mi_size);

  if (Me_size == 0) {
    gmp_gauss_seidel(pInProblem, reaction, velocity, info, options);
    free(reducedProb);
    free(Qreduced);
    free(Rreduced);
    free(Vreduced);
    return;
  }
  listNumericsProblem *curProblem = 0;
  GenericMechanicalProblem *_pnumerics_GMP = genericMechanicalProblem_new();
  curProblem = pInProblem->firstListElem;
  if (Me_size) gmp_add(_pnumerics_GMP, SICONOS_NUMERICS_PROBLEM_EQUALITY, Me_size);
  size_t curPos = 0;
  size_t curPosEq = 0;
  size_t curPosInq = Me_size;
  while (curProblem) {
    size_t curSize = to_size_t(curProblem->size);
    switch (curProblem->type) {
      case SICONOS_NUMERICS_PROBLEM_EQUALITY: {
        memcpy(Vreduced + curPosEq, velocity + curPos, curSize * sizeof(double));
        memcpy(Rreduced + curPosEq, reaction + curPos, curSize * sizeof(double));
        curPosEq += curSize;
        curPos += curSize;
        break;
      }
      case SICONOS_NUMERICS_PROBLEM_LCP: {
        memcpy(Vreduced + curPosInq, velocity + curPos, curSize * sizeof(double));
        memcpy(Rreduced + curPosInq, reaction + curPos, curSize * sizeof(double));
        curPosInq += curSize;
        curPos += curSize;
        gmp_add(_pnumerics_GMP, curProblem->type, to_size_t(curProblem->size));
        break;
      }
      case SICONOS_NUMERICS_PROBLEM_FC3D: {
        memcpy(Vreduced + curPosInq, velocity + curPos, curSize * sizeof(double));
        memcpy(Rreduced + curPosInq, reaction + curPos, curSize * sizeof(double));
        curPosInq += curSize;
        curPos += curSize;
        FrictionContactProblem *pFC3D = (FrictionContactProblem *)gmp_add(
            _pnumerics_GMP, curProblem->type, to_size_t(curProblem->size));
        *(pFC3D->mu) = *(((FrictionContactProblem *)curProblem->problem)->mu);
        break;
      }
      default:
        printf("GMPReduced  buildReducedGMP: problemType unknown: %d . \n", curProblem->type);
    }
    curProblem = curProblem->nextProblem;
  }

#ifdef GMP_DEBUG_GMPREDUCED_SOLVE
  //  printDenseMatrice("newPrb",file,reducedProb,nbRow,nbCol);
  //  printDenseMatrice("newQ",file,Qreduced,nbRow,1);
#endif
  NumericsMatrix numM;
  NM_null(&numM);
  numM.storageType = 0;
  numM.matrix0 = reducedProb;
  numM.matrix1 = 0;
  numM.size0 = to_int(nbRow);
  numM.size1 = to_int(nbCol);
  _pnumerics_GMP->M = &numM;
  _pnumerics_GMP->q = Qreduced;
  gmp_gauss_seidel(_pnumerics_GMP, Rreduced, Vreduced, info, options);
#ifdef GMP_DEBUG_GMPREDUCED_SOLVE
  if (*info) {
    printf("\nGMPREduced2 failed!\n");
  } else {
    printf("\nGMPREduced2 succed!\n");
    //    printDenseMatrice("R",file,Rreduced,nbRow,1);
    //    printDenseMatrice("V",file,Vreduced,nbRow,1);
  }
#endif
  if (!*info) {
    gmp_reduced_convert_solution(pInProblem, reaction, velocity, Rreduced, Rreduced + Me_size,
                                 Vreduced + Me_size);
#ifdef GMP_DEBUG_GMPREDUCED_SOLVE
    //    printDenseMatrice("R2",file,reaction,nbRow,1);
    //    printDenseMatrice("V2",file,velocity,nbRow,1);
#endif
    double err;
    int tolViolate = gmp_compute_error(pInProblem, reaction, velocity,
                                       options->dparam[SICONOS_DPARAM_TOL], options, &err);
    if (tolViolate) {
      printf(
          "GMPReduced2, warnning, reduced problem solved, but error of initial probleme "
          "violated tol = %e, err= %e\n",
          options->dparam[SICONOS_DPARAM_TOL], err);
    }
  }

#ifdef GMP_DEBUG_GMPREDUCED_SOLVE
  //    fclose(file);
#endif
  free(Rreduced);
  free(Vreduced);
  genericMechanicalProblem_free(_pnumerics_GMP, GMP_FREE_GMP);
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
void gmp_reduced_solve(GenericMechanicalProblem *pInProblem, double *reaction,
                       double *velocity, int *info, SolverOptions *options) {
  SparseBlockStructuredMatrix *m = pInProblem->M->matrix1;
  size_t nbRow = m->blocksize0[m->blocknumber0 - 1];
  size_t nbCol = m->blocksize1[m->blocknumber1 - 1];
  double *Me = (double *)malloc(nbRow * nbCol * sizeof(double));
  double *Qe = (double *)malloc(nbRow * sizeof(double));
  double *Mi = (double *)malloc(nbRow * nbCol * sizeof(double));
  double *Qi = (double *)malloc(nbRow * sizeof(double));
  size_t Mesize;
  size_t Misize;
  buildReducedGMP(pInProblem, Me, Mi, Qe, Qi, &Mesize, &Misize);

  if ((Mesize == 0 || Misize == 0)) {
    gmp_gauss_seidel(pInProblem, reaction, velocity, info, options);
    free(Me);
    free(Qe);
    free(Mi);
    free(Qi);
    return;
  }

  const size_t Me_size = Mesize;
  const size_t Mi_size = Misize;
  double *pseduInvMe1 = (double *)malloc(Me_size * Me_size * sizeof(double));
  memcpy(pseduInvMe1, Me, Me_size * Me_size * sizeof(double));
  pinv(pseduInvMe1, to_int(Me_size), to_int(Me_size), 1e-16);
  double *Mi2 = Mi + Mi_size * Me_size;
  double *Mi1 = Mi;
  double *Me2 = Me + Me_size * Me_size;
#ifdef GMP_DEBUG_GMPREDUCED_SOLVE
  double *Me1 = Me;
  FILE *file = fopen("buildReducedGMP_output.txt", "w");
  printf("GMPReducedsolve\n");
  printDenseMatrice("Me1", file, Me1, Me_size, Me_size);
  printDenseMatrice("Me2", file, Me2, Me_size, Mi_size);
  printDenseMatrice("Mi1", file, Mi1, Mi_size, Me_size);
  printDenseMatrice("Mi2", file, Mi2, Mi_size, Mi_size);
  printDenseMatrice("Qe", file, Qe, Me_size, 1);
  printDenseMatrice("Qi", file, Qi, Mi_size, 1);
  printDenseMatrice("Me1inv", file, pseduInvMe1, Me_size, Me_size);
#endif

  double *reducedProb = (double *)malloc(Mi_size * Mi_size * sizeof(double));
  memcpy(reducedProb, Mi2, Mi_size * Mi_size * sizeof(double));

  double *Mi1pseduInvMe1 = (double *)malloc(Mi_size * Me_size * sizeof(double));
  blasint mesize = to_blasint(Me_size);
  blasint misize = to_blasint(Mi_size);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, misize, mesize, mesize, -1.0, Mi1,
              misize, pseduInvMe1, mesize, 0.0, Mi1pseduInvMe1, misize);
#ifdef GMP_DEBUG_GMPREDUCED_SOLVE
  printDenseMatrice("minusMi1pseduInvMe1", file, Mi1pseduInvMe1, Mi_size, Me_size);
  fprintf(file, "_minusMi1pseduInvMe1=-Mi1*Me1inv;\n");
#endif
  cblas_dgemv(CblasColMajor, CblasNoTrans, misize, mesize, 1.0, Mi1pseduInvMe1, misize, Qe, 1,
              1.0, Qi, 1);
#ifdef GMP_DEBUG_GMPREDUCED_SOLVE
  printDenseMatrice("newQi", file, Qi, Mi_size, 1);
  fprintf(file, "_newQi=Qi+_minusMi1pseduInvMe1*Qe;\n");
#endif
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, misize, misize, mesize, 1.0,
              Mi1pseduInvMe1, misize, Me2, mesize, 1.0, reducedProb, misize);
#ifdef GMP_DEBUG_GMPREDUCED_SOLVE
  printDenseMatrice("W", file, reducedProb, Mi_size, Mi_size);
  fprintf(file, "_W=Mi2+_minusMi1pseduInvMe1*Me2;\n");

#endif
  listNumericsProblem *curProblem = 0;
  GenericMechanicalProblem *_pnumerics_GMP = genericMechanicalProblem_new();
  curProblem = pInProblem->firstListElem;
  while (curProblem) {
    switch (curProblem->type) {
      case SICONOS_NUMERICS_PROBLEM_EQUALITY: {
        break;
      }
      case SICONOS_NUMERICS_PROBLEM_LCP: {
        gmp_add(_pnumerics_GMP, curProblem->type, to_size_t(curProblem->size));
        break;
      }
      case SICONOS_NUMERICS_PROBLEM_FC3D: {
        FrictionContactProblem *pFC3D = (FrictionContactProblem *)gmp_add(
            _pnumerics_GMP, curProblem->type, to_size_t(curProblem->size));
        *(pFC3D->mu) = *(((FrictionContactProblem *)curProblem->problem)->mu);
        break;
      }
      default:
        printf("GMPReduced  buildReducedGMP: problemType unknown: %d . \n", curProblem->type);
    }
    curProblem = curProblem->nextProblem;
  }
  NumericsMatrix numM;
  NM_null(&numM);
  numM.storageType = 0;
  numM.matrix0 = reducedProb;
  numM.matrix1 = 0;
  numM.size0 = to_int(Mi_size);
  numM.size1 = numM.size0;
  _pnumerics_GMP->M = &numM;
  _pnumerics_GMP->q = Qi;
  double *Rreduced = (double *)malloc(Mi_size * sizeof(double));
  double *Vreduced = (double *)malloc(Mi_size * sizeof(double));
  gmp_gauss_seidel(_pnumerics_GMP, Rreduced, Vreduced, info, options);
#ifdef GMP_DEBUG_GMPREDUCED_SOLVE
  if (*info) {
    printf("\nGMPREduced failed!\n");
  } else {
    printf("\nGMPREduced succed!\n");
    printDenseMatrice("Ri", file, Rreduced, Mi_size, 1);
    printDenseMatrice("Vi", file, Vreduced, Mi_size, 1);
  }
#endif
  if (!*info) {
    /*Re computation*/
    double *Re = (double *)malloc(Me_size * sizeof(double));
    double *Rbuf = (double *)malloc(Me_size * sizeof(double));
    memcpy(Rbuf, Qe, Me_size * sizeof(double));
    cblas_dgemv(CblasColMajor, CblasNoTrans, mesize, misize, 1.0, Me2, mesize, Rreduced, 1,
                1.0, Rbuf, 1);
    cblas_dgemv(CblasColMajor, CblasNoTrans, mesize, mesize, -1.0, pseduInvMe1, mesize, Rbuf,
                1, 0.0, Re, 1);
#ifdef GMP_DEBUG_GMPREDUCED_SOLVE
    fprintf(file, "_Re=-Me1inv*(Me2*Ri+Qe);\n");
    printDenseMatrice("Re", file, Re, Me_size, 1);
#endif
    gmp_reduced_convert_solution(pInProblem, reaction, velocity, Re, Rreduced, Vreduced);
    double err;
    int tolViolate = gmp_compute_error(pInProblem, reaction, velocity,
                                       options->dparam[SICONOS_DPARAM_TOL], options, &err);
    if (tolViolate) {
      printf(
          "GMPReduced, warnning, reduced problem solved, but error of initial probleme "
          "violated tol = %e, err= %e\n",
          options->dparam[SICONOS_DPARAM_TOL], err);
    }
    free(Re);
    free(Rbuf);
  }

#ifdef GMP_DEBUG_GMPREDUCED_SOLVE
  fclose(file);
#endif
  free(Rreduced);
  free(Vreduced);
  genericMechanicalProblem_free(_pnumerics_GMP, GMP_FREE_GMP);
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

void gmp_as_mlcp(GenericMechanicalProblem *pInProblem, double *reaction, double *velocity,
                 int *info, SolverOptions *options) {
  /*First, we don't manage FC3D.*/
  listNumericsProblem *curProblem = 0;
  curProblem = pInProblem->firstListElem;
  while (curProblem) {
    switch (curProblem->type) {
      case SICONOS_NUMERICS_PROBLEM_EQUALITY:
      case SICONOS_NUMERICS_PROBLEM_LCP:
        break;
      case SICONOS_NUMERICS_PROBLEM_FC3D: {
        printf("gmp_as_mlcp Numerics ERROR: gmp_as_mlcp doesn't deal with FC3D.\n");
        *info = 1;
        return;
      }
      default:
        printf("gmp_as_mlcp Numerics : gmp_gauss_seidel unknown problem type %d.\n",
               curProblem->type);
    }
    curProblem = curProblem->nextProblem;
  }
  size_t Me_size;
  size_t Mi_size;

  SparseBlockStructuredMatrix *m = pInProblem->M->matrix1;
  size_t nbRow = m->blocksize0[m->blocknumber0 - 1];
  size_t nbCol = m->blocksize1[m->blocknumber1 - 1];

  double *reducedProb = (double *)malloc(nbRow * nbCol * sizeof(double));
  double *Qreduced = (double *)malloc(nbRow * sizeof(double));
  _GMPReducedEquality(pInProblem, reducedProb, Qreduced, &Me_size, &Mi_size);

  if (!Me_size) {
    /*it is a lcp.*/
    LinearComplementarityProblem aLCP;
    SolverOptions *aLcpOptions = solver_options_create(SICONOS_LCP_ENUM);
    NumericsMatrix M;
    NM_null(&M);
    M.storageType = 0;
    M.size0 = to_int(Mi_size);
    M.size1 = to_int(Mi_size);
    M.matrix0 = reducedProb;
    M.matrix1 = 0;
    aLCP.size = to_int(Mi_size);
    aLCP.q = Qreduced;
    aLCP.M = &M;
    lcp_enum_init(&aLCP, aLcpOptions, 1);
    *info = linearComplementarity_driver(&aLCP, reaction, velocity, aLcpOptions);
    lcp_enum_reset(&aLCP, aLcpOptions, 1);
    solver_options_delete(aLcpOptions);
    aLcpOptions = NULL;

    goto END_GMP3;
  }
  if (!Mi_size) {
    /*it is a linear system.*/
    for (size_t i = 0; i < (size_t)Me_size; ++i) reaction[i] = -Qreduced[i];
    NumericsMatrix M;
    NM_null(&M);
    NM_fill(&M, NM_DENSE, to_int(Me_size), to_int(Me_size), reducedProb);
    // *info = NM_gesv(&M, reaction, true);
    *info = NM_LU_solve(&M, reaction, 1);
    M.matrix0 = NULL;
    NM_clear(&M);
    goto END_GMP3;
  }
  /*it is a MLCP*/
  MixedLinearComplementarityProblem aMLCP;
  SolverOptions *aMlcpOptions = solver_options_create(SICONOS_MLCP_ENUM);
  aMLCP.n = to_int(Me_size);
  aMLCP.m = to_int(Mi_size);
  aMLCP.blocksRows = 0;
  aMLCP.blocksIsComp = 0;
  aMLCP.isStorageType1 = 1;
  aMLCP.isStorageType2 = 0;

  aMLCP.A = 0;
  aMLCP.B = 0;
  aMLCP.C = 0;
  aMLCP.D = 0;
  aMLCP.a = 0;
  aMLCP.b = 0;
  aMLCP.q = Qreduced;
  NumericsMatrix M;
  NM_null(&M);
  M.storageType = 0;
  M.size0 = to_int(Mi_size + Me_size);
  M.size1 = to_int(Mi_size + Me_size);
  M.matrix0 = reducedProb;
  M.matrix1 = 0;
  aMLCP.M = &M;
  mlcp_driver_init(&aMLCP, aMlcpOptions);
  aMlcpOptions->dparam[SICONOS_DPARAM_TOL] = options->dparam[SICONOS_DPARAM_TOL];
  *info = mlcp_driver(&aMLCP, reaction, velocity, aMlcpOptions);

  mlcp_driver_reset(&aMLCP, aMlcpOptions);
  solver_options_delete(aMlcpOptions);
  aMlcpOptions = NULL;
END_GMP3:;

  free(reducedProb);
  free(Qreduced);
}
