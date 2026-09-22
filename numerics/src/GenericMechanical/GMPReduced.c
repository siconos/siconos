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

/*! \file GMPReduced.c
 *  \brief Reduced-strategy solvers for GenericMechanicalProblem.
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
#include "GenericMechanicalProblem.h"      // for GMP_LocalProblem
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
#include "numerics_verbose.h"              // for numerics_error, numerics_warning
#include "pinv.h"                          // for pinv
#include "siconos_debug.h"                 // for DEBUG_PRINTF, DEBUG_EXPR_WE, DEBUG_PRINT_MAT_STR

/** Print a dense column-major matrix in Scilab-compatible format.
 *  If file is NULL the matrix is printed to stdout.
 *  \param[in] name variable name to print
 *  \param[in] file output file, or NULL for stdout
 *  \param[in] m column-major matrix data
 *  \param[in] N number of rows
 *  \param[in] M number of columns
 */
static void printDenseMatrix(char *name, FILE *file, double *m, int N, int M) {
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

/** Convert a reduced solution back to the original GenericMechanicalProblem layout.
 *  Equality blocks take their value from Re, while LCP/FC3D blocks take (Rreduced, Vreduced).
 *  \param[in] problem original problem description
 *  \param[out] reaction global reaction vector
 *  \param[out] velocity global velocity vector
 *  \param[in] Re equality-block reactions
 *  \param[in] Rreduced inequality-block reactions
 *  \param[in] Vreduced inequality-block velocities
 */
void gmp_reduced_convert_solution(GenericMechanicalProblem *problem, double *reaction,
                                  double *velocity, double *Re, double *Rreduced,
                                  double *Vreduced) {
  GMP_LocalProblem *local = 0;
  local = problem->firstLocal;
  size_t local_size = 0;
  while (local) {
    local_size = to_size_t(local->size);
    switch (local->type) {
      case SICONOS_NUMERICS_PROBLEM_EQUALITY: {
        memcpy(reaction, Re, local_size * sizeof(double));
        for (size_t i = 0; i < local_size; i++) velocity[i] = 0.0;
        Re += local_size;
        break;
      }
      case SICONOS_NUMERICS_PROBLEM_LCP:
      case SICONOS_NUMERICS_PROBLEM_FC3D: {
        memcpy(reaction, Rreduced, local_size * sizeof(double));
        memcpy(velocity, Vreduced, local_size * sizeof(double));
        Rreduced += local_size;
        Vreduced += local_size;
        break;
      }
      default:
        (void)numerics_warning("gmp_reduced_convert_solution", "unknown problem type");
    }
    reaction += local_size;
    velocity += local_size;
    local = local->next;
  }
}
/** Compute the total size of equality blocks (Me) and inequality blocks (Mi).
 *  \param[in] problem the GenericMechanicalProblem
 *  \param[out] Me_size total size of equality blocks
 *  \param[out] Mi_size total size of inequality blocks
 */
/** Build a temporary GenericMechanicalProblem for the reduced Gauss-Seidel solve.
 *
 *  The temporary problem contains an optional leading equality block followed by all
 *  inequality blocks (LCP/FC3D) of the original problem. FC3D friction coefficients are
 *  copied from the original problem.
 *
 *  \param[in] problem the original GenericMechanicalProblem
 *  \param[in] equality_size size of the leading equality block (0 to omit it)
 *  \return a new GenericMechanicalProblem, or NULL on allocation failure
 */
static GenericMechanicalProblem *gmp_build_reduced_gs_problem(GenericMechanicalProblem *problem,
                                                              size_t equality_size) {
  GenericMechanicalProblem *reduced = genericMechanicalProblem_new();
  if (!reduced) return NULL;

  if (equality_size) {
    if (!gmp_add(reduced, SICONOS_NUMERICS_PROBLEM_EQUALITY, equality_size)) {
      genericMechanicalProblem_free(reduced, GMP_FREE_GMP);
      return NULL;
    }
  }

  GMP_LocalProblem *local = problem->firstLocal;
  while (local) {
    if (local->type != SICONOS_NUMERICS_PROBLEM_EQUALITY) {
      void *prb = gmp_add(reduced, local->type, to_size_t(local->size));
      if (!prb) {
        genericMechanicalProblem_free(reduced, GMP_FREE_GMP);
        return NULL;
      }
      if (local->type == SICONOS_NUMERICS_PROBLEM_FC3D) {
        FrictionContactProblem *pFC3D = (FrictionContactProblem *)prb;
        *pFC3D->mu = *(((FrictionContactProblem *)local->problem)->mu);
      }
    }
    local = local->next;
  }
  return reduced;
}

static void _GMPReducedGetSizes(GenericMechanicalProblem *problem, size_t *Me_size,
                                size_t *Mi_size) {
  GMP_LocalProblem *local = 0;
  (*Me_size) = 0;
  (*Mi_size) = 0;
  local = problem->firstLocal;
  while (local) {
    if (local->type == SICONOS_NUMERICS_PROBLEM_EQUALITY) {
      (*Me_size) += to_size_t(local->size);
      ;
    } else {
      (*Mi_size) += to_size_t(local->size);
    }
    local = local->next;
  }
}

/** Build the reduced matrices Me, Mi and reduced right-hand sides Qe, Qi.
 *
 *  The rows of the global sparse-block matrix are permuted so that equality rows come first,
 *  followed by inequality rows. The output is stored in column-major dense format.
 *
 *  \param[in] problem the GenericMechanicalProblem (must use sparse-block storage)
 *  \param[out] Me equality rows of the permuted matrix, size Me_Size x (Me_Size + Mi_Size)
 *  \param[out] Mi inequality rows of the permuted matrix, size Mi_Size x (Me_Size + Mi_Size)
 *  \param[out] Qe equality part of the permuted right-hand side
 *  \param[out] Qi inequality part of the permuted right-hand side
 *  \param[out] Me_Size number of equality rows
 *  \param[out] Mi_Size number of inequality rows
 */
static int buildReducedGMP(GenericMechanicalProblem *problem, double *Me, double *Mi,
                           double *Qe, double *Qi, size_t *Me_Size, size_t *Mi_Size) {
  assert(problem->M->storageType == NM_SPARSE_BLOCK);
  SparseBlockStructuredMatrix *m = problem->M->matrix1;
  DEBUG_EXPR_WE({
    FILE *file = fopen("buildReducedGMP_input.txt", "w");
    if (file) {
      SBM_write_in_fileForScilab(m, file);
      fclose(file);
    }
  });
  size_t local_size = 0;
  size_t numberOfBlockColumns = m->blocknumber1;
  size_t *newIndexOfCol = (size_t *)malloc(numberOfBlockColumns * sizeof(size_t));
  if (!newIndexOfCol) {
    (void)numerics_error("buildReducedGMP", "memory allocation failed");
    return 1;
  }

  /*Me building*/
  size_t MeRow = 0;
  size_t MiRow = 0;

  /**size of Me */
  GMP_LocalProblem *local = 0;
  size_t numberOfEqualityBlockRows = 0;
  size_t numberOfInequalityBlockRows = 0;
  size_t blockRowIndex = 0;
  local = problem->firstLocal;
  while (local) {
    if (blockRowIndex)
      local_size = m->blocksize0[blockRowIndex] - m->blocksize0[blockRowIndex - 1];
    else
      local_size = m->blocksize0[blockRowIndex];

    if (local->type == SICONOS_NUMERICS_PROBLEM_EQUALITY) {
      numberOfEqualityBlockRows++;
      MeRow += local_size;
    } else {
      numberOfInequalityBlockRows++;
      MiRow += local_size;
    }
    local = local->next;
    blockRowIndex++;
  }
  blockRowIndex = 0;
  size_t equalityRowIndex = 0;
  size_t inequalityRowIndex = 0;
  size_t rowIndex = 0;
  local = problem->firstLocal;
  while (local) {
    if (local->type == SICONOS_NUMERICS_PROBLEM_EQUALITY) {
      newIndexOfCol[rowIndex] = equalityRowIndex;
      equalityRowIndex++;
    } else {
      newIndexOfCol[rowIndex] = inequalityRowIndex + numberOfEqualityBlockRows;
      inequalityRowIndex++;
    }
    rowIndex++;
    local = local->next;
  }
  DEBUG_PRINTF("buildReducedGMP equality row index=%i. inequality row index=%i\n", (int)equalityRowIndex, (int)inequalityRowIndex);

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
  local = problem->firstLocal;

  /**mem alloc for Me and Mi*/
  // int numberOfColumns=MeRow+MiRow;
  *Me_Size = MeRow;
  *Mi_Size = MiRow;

  /** copy rows into Me and Mi */
  size_t currentPosition = 0;
  for (size_t blockRowIndex = 0; blockRowIndex < numberOfEqualityBlockRows; blockRowIndex++) {
    SBM_row_to_dense(Morder, blockRowIndex, Me, currentPosition, MeRow);
    currentPosition = Morder->blocksize1[blockRowIndex];
  }
  currentPosition = 0;
  size_t firtMiLine = 0;
  if (numberOfInequalityBlockRows > 0) firtMiLine = Morder->blocksize1[numberOfEqualityBlockRows];

  for (size_t blockRowIndex = numberOfEqualityBlockRows; blockRowIndex < numberOfEqualityBlockRows + numberOfInequalityBlockRows;
       blockRowIndex++) {
    currentPosition = Morder->blocksize1[blockRowIndex] - firtMiLine;
    SBM_row_to_dense(Morder, blockRowIndex, Mi, currentPosition, MiRow);
  }
  Morder = SBM_free(Morder, SBM_FREE_KEEP_BLOCKS);

  local = problem->firstLocal;
  int currentBlockIndex = 0;
  double *curQ = problem->q;
  double *curQe = Qe;
  double *curQi = Qi;
  currentBlockIndex = 0;
  while (local) {
    if (currentBlockIndex) {
      local_size = m->blocksize0[currentBlockIndex] - m->blocksize0[currentBlockIndex - 1];
    } else {
      local_size = m->blocksize0[currentBlockIndex];
    }

    switch (local->type) {
      case SICONOS_NUMERICS_PROBLEM_EQUALITY: {
        /** copy the current line block in Me*/
        memcpy(curQe, curQ, local_size * sizeof(double));
        curQe += local_size;
        break;
      }
      case SICONOS_NUMERICS_PROBLEM_LCP:
      case SICONOS_NUMERICS_PROBLEM_FC3D: {
        memcpy(curQi, curQ, local_size * sizeof(double));
        curQi += local_size;
        break;
      }
      default:
        (void)numerics_warning("buildReducedGMP", "unknown problem type");
    }
    local = local->next;
    curQ += local_size;
    currentBlockIndex++;
  }
  DEBUG_EXPR_WE({
    size_t numberOfColumns = MeRow + MiRow;
    DEBUG_PRINT("The Me matrix is:\n");
    DEBUG_PRINT_MAT_STR("Me", Me, (unsigned)MeRow, (unsigned)numberOfColumns);
    DEBUG_PRINT_MAT_STR("Qe", Qe, (unsigned)MeRow, 1);
    DEBUG_PRINT("The Mi matrix is:\n");
    DEBUG_PRINT_MAT_STR("Mi", Mi, (unsigned)MiRow, (unsigned)numberOfColumns);
    DEBUG_PRINT_MAT_STR("Qi", Qi, (unsigned)MiRow, 1);
  });
  free(newIndexOfCol);
  return 0;
}

/** Assemble the reduced problem with equalities grouped in one block.
 *
 *  Output matrices are stored in column-major dense format:
 *    - reducedProb is numberOfRows x numberOfColumns with equality rows on top,
 *    - Qreduced is the permuted right-hand side.
 *
 *  \param[in] problem the original GenericMechanicalProblem
 *  \param[out] reducedProb dense reduced matrix (must be pre-allocated to numberOfRows * numberOfColumns)
 *  \param[out] Qreduced dense reduced right-hand side (must be pre-allocated to numberOfRows)
 *  \param[out] Me_size number of equality rows
 *  \param[out] Mi_size number of inequality rows
 */
static int _GMPReducedEquality(GenericMechanicalProblem *problem, double *reducedProb,
                               double *Qreduced, size_t *Me_size, size_t *Mi_size) {
  SparseBlockStructuredMatrix *m = problem->M->matrix1;
  size_t numberOfRows = m->blocksize0[m->blocknumber0 - 1];
  size_t numberOfColumns = m->blocksize1[m->blocknumber1 - 1];

  _GMPReducedGetSizes(problem, Me_size, Mi_size);
  if (*Me_size == 0) {
    memcpy(Qreduced, problem->q, (*Mi_size) * sizeof(double));
    SBM_to_dense(m, reducedProb);
    return 0;
  }

  double *Me = (*Me_size) ? (double *)malloc((*Me_size) * numberOfColumns * sizeof(double)) : NULL;
  double *Mi = (*Mi_size) ? (double *)malloc((*Mi_size) * numberOfColumns * sizeof(double)) : NULL;
  double *Qi = (double *)malloc(numberOfRows * sizeof(double));
  if ((*Me_size && !Me) || (*Mi_size && !Mi) || !Qi) {
    (void)numerics_error("_GMPReducedEquality", "memory allocation failed");
    free(Me);
    free(Mi);
    free(Qi);
    return 1;
  }
  if (buildReducedGMP(problem, Me, Mi, Qreduced, Qi, Me_size, Mi_size)) {
    free(Me);
    free(Mi);
    free(Qi);
    return 1;
  }

  DEBUG_EXPR_WE({
    double *Me1 = Me;
    double *Me2 = Me + (*Me_size) * (*Me_size);
    double *Mi1 = Mi;
    double *Mi2 = Mi + (*Mi_size) * (*Me_size);
    FILE *file = fopen("buildReduced2GMP_output.txt", "w");
    DEBUG_PRINT("GMP2Reducedsolve\n");
    if (file) {
      printDenseMatrix("Me1", file, Me1, (int)(*Me_size), (int)(*Me_size));
      printDenseMatrix("Me2", file, Me2, (int)(*Me_size), (int)(*Mi_size));
      printDenseMatrix("Mi1", file, Mi1, (int)(*Mi_size), (int)(*Me_size));
      printDenseMatrix("Mi2", file, Mi2, (int)(*Mi_size), (int)(*Mi_size));
      printDenseMatrix("Qe", file, Qreduced, (int)(*Me_size), 1);
      printDenseMatrix("Qi", file, Qi, (int)(*Mi_size), 1);
      fclose(file);
    }
  });
  for (size_t columnIndex = 0; columnIndex < numberOfColumns; columnIndex++) {
    if (*Me_size)
      memcpy(reducedProb + columnIndex * numberOfRows, Me + columnIndex * (*Me_size),
             (*Me_size) * sizeof(double));
    if (*Mi_size)
      memcpy(reducedProb + columnIndex * numberOfRows + (*Me_size), Mi + columnIndex * (*Mi_size),
             (*Mi_size) * sizeof(double));
  }
  if (*Mi_size) memcpy(Qreduced + (*Me_size), Qi, (*Mi_size) * sizeof(double));
  free(Me);
  free(Mi);
  free(Qi);
  return 0;
}

/** Solve a GenericMechanicalProblem by assembling all equalities into a single block.
 *
 *  The reduced problem is solved with a Non-Smooth Gauss-Seidel sweep. Equality variables
 *  are then recovered and the solution is converted back to the original layout.
 *
 *  \param[in] problem the GenericMechanicalProblem to solve
 *  \param[out] reaction global reaction vector
 *  \param[out] velocity global velocity vector
 *  \param[out] info solver status (0 on success)
 *  \param[in,out] options solver options
 */
void gmp_reduced_equality_solve(GenericMechanicalProblem *problem, double *reaction,
                                double *velocity, int *info, SolverOptions *options) {
  SparseBlockStructuredMatrix *m = problem->M->matrix1;
  size_t numberOfRows = m->blocksize0[m->blocknumber0 - 1];
  size_t numberOfColumns = m->blocksize1[m->blocknumber1 - 1];

  size_t Me_size = 0;
  size_t Mi_size = 0;
  double *reducedProb = NULL;
  double *Qreduced = NULL;
  double *Rreduced = NULL;
  double *Vreduced = NULL;
  GenericMechanicalProblem *_pnumerics_GMP = NULL;

  reducedProb = (double *)malloc(numberOfRows * numberOfColumns * sizeof(double));
  Qreduced = (double *)malloc(numberOfRows * sizeof(double));
  Rreduced = (double *)calloc(numberOfColumns, sizeof(double));
  Vreduced = (double *)calloc(numberOfRows, sizeof(double));
  if (!reducedProb || !Qreduced || !Rreduced || !Vreduced) {
    *info = 1;
    (void)numerics_error("gmp_reduced_equality_solve", "memory allocation failed");
    goto cleanup;
  }

  if (_GMPReducedEquality(problem, reducedProb, Qreduced, &Me_size, &Mi_size)) {
    *info = 1;
    goto cleanup;
  }

  if (Me_size == 0) {
    gmp_gauss_seidel(problem, reaction, velocity, info, options);
    goto cleanup;
  }

  _pnumerics_GMP = gmp_build_reduced_gs_problem(problem, Me_size);
  if (!_pnumerics_GMP) {
    *info = 1;
    (void)numerics_error("gmp_reduced_equality_solve", "failed to build reduced problem");
    goto cleanup;
  }

  /* Copy initial guesses for the reduced problem. */
  GMP_LocalProblem *local = problem->firstLocal;
  size_t currentPosition = 0;
  size_t equalityPosition = 0;
  size_t inequalityPosition = Me_size;
  while (local) {
    size_t local_size = to_size_t(local->size);
    switch (local->type) {
      case SICONOS_NUMERICS_PROBLEM_EQUALITY: {
        memcpy(Vreduced + equalityPosition, velocity + currentPosition, local_size * sizeof(double));
        memcpy(Rreduced + equalityPosition, reaction + currentPosition, local_size * sizeof(double));
        equalityPosition += local_size;
        currentPosition += local_size;
        break;
      }
      case SICONOS_NUMERICS_PROBLEM_LCP:
      case SICONOS_NUMERICS_PROBLEM_FC3D: {
        memcpy(Vreduced + inequalityPosition, velocity + currentPosition, local_size * sizeof(double));
        memcpy(Rreduced + inequalityPosition, reaction + currentPosition, local_size * sizeof(double));
        inequalityPosition += local_size;
        currentPosition += local_size;
        break;
      }
      default:
        (void)numerics_warning("gmp_reduced_equality_solve", "unknown problem type");
    }
    local = local->next;
  }

  NumericsMatrix numM;
  NM_null(&numM);
  numM.storageType = NM_DENSE;
  numM.matrix0 = reducedProb;
  numM.matrix1 = NULL;
  numM.size0 = to_int(numberOfRows);
  numM.size1 = to_int(numberOfColumns);
  _pnumerics_GMP->M = &numM;
  _pnumerics_GMP->q = Qreduced;
  gmp_gauss_seidel(_pnumerics_GMP, Rreduced, Vreduced, info, options);
  DEBUG_PRINTF("GMPREduced2 %s\n", *info ? "failed" : "succed");
  if (!*info) {
    gmp_reduced_convert_solution(problem, reaction, velocity, Rreduced, Rreduced + Me_size,
                                 Vreduced + Me_size);
    double error;
    int toleranceViolation = gmp_compute_error(problem, reaction, velocity,
                                       options->dparam[SICONOS_DPARAM_TOL], options, &error);
    if (toleranceViolation) {
      numerics_printf(
          "GMPReduced_equality_solve: reduced problem solved, but original error violated "
          "tolerance = %e, error = %e\n",
          options->dparam[SICONOS_DPARAM_TOL], error);
    }
  }

cleanup:
  free(Rreduced);
  free(Vreduced);
  if (_pnumerics_GMP) genericMechanicalProblem_free(_pnumerics_GMP, GMP_FREE_GMP);
  free(Qreduced);
  free(reducedProb);
}

/** Solve a GenericMechanicalProblem by eliminating equality blocks.
 *
 *  Equalities are eliminated using the pseudo-inverse of Me_1:
 *    Re = -Me_1^{-1} (Me_2 Ri + Qe)
 *    Vi = (Mi_2 - Mi_1 Me_1^{-1} Me_2) Ri + Qi - Mi_1 Me_1^{-1} Qe
 *  The remaining inequality-only problem is solved with NSGS, then Re is recovered.
 *
 *  \param[in] problem the GenericMechanicalProblem to solve
 *  \param[out] reaction global reaction vector
 *  \param[out] velocity global velocity vector
 *  \param[out] info solver status (0 on success)
 *  \param[in,out] options solver options
 */
void gmp_reduced_solve(GenericMechanicalProblem *problem, double *reaction,
                       double *velocity, int *info, SolverOptions *options) {
  SparseBlockStructuredMatrix *m = problem->M->matrix1;
  size_t numberOfRows = m->blocksize0[m->blocknumber0 - 1];
  size_t numberOfColumns = m->blocksize1[m->blocknumber1 - 1];

  size_t Mesize = 0;
  size_t Misize = 0;
  double *Me = NULL;
  double *Qe = NULL;
  double *Mi = NULL;
  double *Qi = NULL;
  double *pseudoInvMe1 = NULL;
  double *reducedProb = NULL;
  double *Mi1pseudoInvMe1 = NULL;
  double *Rreduced = NULL;
  double *Vreduced = NULL;
  double *Re = NULL;
  double *Rbuf = NULL;
  GenericMechanicalProblem *_pnumerics_GMP = NULL;

  Me = (double *)malloc(numberOfRows * numberOfColumns * sizeof(double));
  Qe = (double *)malloc(numberOfRows * sizeof(double));
  Mi = (double *)malloc(numberOfRows * numberOfColumns * sizeof(double));
  Qi = (double *)malloc(numberOfRows * sizeof(double));
  if (!Me || !Qe || !Mi || !Qi) {
    *info = 1;
    (void)numerics_error("gmp_reduced_solve", "memory allocation failed");
    goto cleanup;
  }

  if (buildReducedGMP(problem, Me, Mi, Qe, Qi, &Mesize, &Misize)) {
    *info = 1;
    goto cleanup;
  }

  if (Mesize == 0 || Misize == 0) {
    gmp_gauss_seidel(problem, reaction, velocity, info, options);
    goto cleanup;
  }

  const size_t Me_size = Mesize;
  const size_t Mi_size = Misize;

  pseudoInvMe1 = (double *)malloc(Me_size * Me_size * sizeof(double));
  reducedProb = (double *)malloc(Mi_size * Mi_size * sizeof(double));
  Mi1pseudoInvMe1 = (double *)malloc(Mi_size * Me_size * sizeof(double));
  Rreduced = (double *)malloc(Mi_size * sizeof(double));
  Vreduced = (double *)malloc(Mi_size * sizeof(double));
  if (!pseudoInvMe1 || !reducedProb || !Mi1pseudoInvMe1 || !Rreduced || !Vreduced) {
    *info = 1;
    (void)numerics_error("gmp_reduced_solve", "memory allocation failed");
    goto cleanup;
  }

  double *Mi2 = Mi + Mi_size * Me_size;
  double *Mi1 = Mi;
  double *Me2 = Me + Me_size * Me_size;

  memcpy(pseudoInvMe1, Me, Me_size * Me_size * sizeof(double));
  pinv(pseudoInvMe1, to_int(Me_size), to_int(Me_size), 1e-16);
  memcpy(reducedProb, Mi2, Mi_size * Mi_size * sizeof(double));

  DEBUG_EXPR_WE({
    FILE *file = fopen("buildReducedGMP_output.txt", "w");
    DEBUG_PRINT("GMPReducedsolve\n");
    if (file) {
      printDenseMatrix("Me1", file, Me, (int)Me_size, (int)Me_size);
      printDenseMatrix("Me2", file, Me2, (int)Me_size, (int)Mi_size);
      printDenseMatrix("Mi1", file, Mi1, (int)Mi_size, (int)Me_size);
      printDenseMatrix("Mi2", file, Mi2, (int)Mi_size, (int)Mi_size);
      printDenseMatrix("Qe", file, Qe, (int)Me_size, 1);
      printDenseMatrix("Qi", file, Qi, (int)Mi_size, 1);
      printDenseMatrix("Me1inv", file, pseudoInvMe1, (int)Me_size, (int)Me_size);
      fclose(file);
    }
  });

  blasint mesize = to_blasint(Me_size);
  blasint misize = to_blasint(Mi_size);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, misize, mesize, mesize, -1.0, Mi1,
              misize, pseudoInvMe1, mesize, 0.0, Mi1pseudoInvMe1, misize);
  DEBUG_EXPR_WE({
    FILE *file = fopen("buildReducedGMP_output.txt", "a");
    if (file) {
      printDenseMatrix("minusMi1pseudoInvMe1", file, Mi1pseudoInvMe1, (int)Mi_size, (int)Me_size);
      fprintf(file, "_minusMi1pseudoInvMe1=-Mi1*Me1inv;\n");
      fclose(file);
    }
  });
  cblas_dgemv(CblasColMajor, CblasNoTrans, misize, mesize, 1.0, Mi1pseudoInvMe1, misize, Qe, 1,
              1.0, Qi, 1);
  DEBUG_EXPR_WE({
    FILE *file = fopen("buildReducedGMP_output.txt", "a");
    if (file) {
      printDenseMatrix("newQi", file, Qi, (int)Mi_size, 1);
      fprintf(file, "_newQi=Qi+_minusMi1pseudoInvMe1*Qe;\n");
      fclose(file);
    }
  });
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, misize, misize, mesize, 1.0,
              Mi1pseudoInvMe1, misize, Me2, mesize, 1.0, reducedProb, misize);
  DEBUG_EXPR_WE({
    FILE *file = fopen("buildReducedGMP_output.txt", "a");
    if (file) {
      printDenseMatrix("W", file, reducedProb, (int)Mi_size, (int)Mi_size);
      fprintf(file, "_W=Mi2+_minusMi1pseudoInvMe1*Me2;\n");
      fclose(file);
    }
  });

  _pnumerics_GMP = gmp_build_reduced_gs_problem(problem, 0);
  if (!_pnumerics_GMP) {
    *info = 1;
    (void)numerics_error("gmp_reduced_solve", "failed to build reduced problem");
    goto cleanup;
  }

  NumericsMatrix numM;
  NM_null(&numM);
  numM.storageType = NM_DENSE;
  numM.matrix0 = reducedProb;
  numM.matrix1 = NULL;
  numM.size0 = to_int(Mi_size);
  numM.size1 = numM.size0;
  _pnumerics_GMP->M = &numM;
  _pnumerics_GMP->q = Qi;
  gmp_gauss_seidel(_pnumerics_GMP, Rreduced, Vreduced, info, options);
  DEBUG_EXPR_WE({
    FILE *file = fopen("buildReducedGMP_output.txt", "a");
    if (file) {
      if (*info) {
        fprintf(file, "\nGMPREduced failed!\n");
      } else {
        fprintf(file, "\nGMPREduced succed!\n");
        printDenseMatrix("Ri", file, Rreduced, (int)Mi_size, 1);
        printDenseMatrix("Vi", file, Vreduced, (int)Mi_size, 1);
      }
      fclose(file);
    }
  });
  if (!*info) {
    Re = (double *)malloc(Me_size * sizeof(double));
    Rbuf = (double *)malloc(Me_size * sizeof(double));
    if (!Re || !Rbuf) {
      *info = 1;
      (void)numerics_error("gmp_reduced_solve", "memory allocation failed");
      goto cleanup;
    }
    memcpy(Rbuf, Qe, Me_size * sizeof(double));
    cblas_dgemv(CblasColMajor, CblasNoTrans, mesize, misize, 1.0, Me2, mesize, Rreduced, 1,
                1.0, Rbuf, 1);
    cblas_dgemv(CblasColMajor, CblasNoTrans, mesize, mesize, -1.0, pseudoInvMe1, mesize, Rbuf,
                1, 0.0, Re, 1);
    DEBUG_EXPR_WE({
      FILE *file = fopen("buildReducedGMP_output.txt", "a");
      if (file) {
        fprintf(file, "_Re=-Me1inv*(Me2*Ri+Qe);\n");
        printDenseMatrix("Re", file, Re, (int)Me_size, 1);
        fclose(file);
      }
    });
    gmp_reduced_convert_solution(problem, reaction, velocity, Re, Rreduced, Vreduced);
    double error;
    int toleranceViolation = gmp_compute_error(problem, reaction, velocity,
                                       options->dparam[SICONOS_DPARAM_TOL], options, &error);
    if (toleranceViolation) {
      numerics_printf(
          "GMPReduced_solve: reduced problem solved, but original error violated "
          "tolerance = %e, error = %e\n",
          options->dparam[SICONOS_DPARAM_TOL], error);
    }
  }

cleanup:
  free(Re);
  free(Rbuf);
  free(Rreduced);
  free(Vreduced);
  if (_pnumerics_GMP) genericMechanicalProblem_free(_pnumerics_GMP, GMP_FREE_GMP);
  free(Me);
  free(Mi);
  free(Qe);
  free(Qi);
  free(pseudoInvMe1);
  free(reducedProb);
  free(Mi1pseudoInvMe1);
}

/** Solve a GenericMechanicalProblem by reformulating it as an MLCP.
 *
 *  Equalities are assembled into one block. If no equality exists the problem is an LCP;
 *  if no inequality exists it is a linear system. Otherwise it is solved as a mixed LCP.
 *  FC3D blocks are not supported by this solver.
 *
 *  \param[in] problem the GenericMechanicalProblem to solve
 *  \param[out] reaction global reaction vector
 *  \param[out] velocity global velocity vector
 *  \param[out] info solver status (0 on success)
 *  \param[in,out] options solver options
 */
static void gmp_as_mlcp_lcp(double *reducedProb, double *Qreduced, size_t Mi_size,
                            double *reaction, double *velocity, int *info) {
  LinearComplementarityProblem aLCP;
  SolverOptions *aLcpOptions = solver_options_create(SICONOS_LCP_ENUM);
  NumericsMatrix M;
  NM_null(&M);
  NM_fill(&M, NM_DENSE, to_int(Mi_size), to_int(Mi_size), reducedProb);
  aLCP.size = to_int(Mi_size);
  aLCP.q = Qreduced;
  aLCP.M = &M;
  lcp_enum_init(&aLCP, aLcpOptions, 1);
  *info = linearComplementarity_driver(&aLCP, reaction, velocity, aLcpOptions);
  lcp_enum_reset(&aLCP, aLcpOptions, 1);
  solver_options_delete(aLcpOptions);
}

static void gmp_as_mlcp_linear_system(double *reducedProb, double *Qreduced, size_t Me_size,
                                      double *reaction, int *info) {
  for (size_t i = 0; i < (size_t)Me_size; ++i) reaction[i] = -Qreduced[i];
  NumericsMatrix M;
  NM_null(&M);
  NM_fill(&M, NM_DENSE, to_int(Me_size), to_int(Me_size), reducedProb);
  *info = NM_LU_solve(&M, reaction, 1);
  M.matrix0 = NULL;  /* reducedProb is owned by the caller, do not free it. */
  NM_clear(&M);
}

static void gmp_as_mlcp_mlcp(double *reducedProb, double *Qreduced, size_t Me_size,
                             size_t Mi_size, double *reaction, double *velocity, int *info,
                             SolverOptions *options) {
  MixedLinearComplementarityProblem aMLCP;
  SolverOptions *aMlcpOptions = solver_options_create(SICONOS_MLCP_ENUM);
  aMLCP.n = to_int(Me_size);
  aMLCP.m = to_int(Mi_size);
  aMLCP.blocksRows = NULL;
  aMLCP.blocksIsComp = NULL;
  aMLCP.isStorageType1 = 1;
  aMLCP.isStorageType2 = 0;
  aMLCP.A = NULL;
  aMLCP.B = NULL;
  aMLCP.C = NULL;
  aMLCP.D = NULL;
  aMLCP.a = NULL;
  aMLCP.b = NULL;
  aMLCP.q = Qreduced;

  NumericsMatrix M;
  NM_null(&M);
  NM_fill(&M, NM_DENSE, to_int(Mi_size + Me_size), to_int(Mi_size + Me_size), reducedProb);
  aMLCP.M = &M;

  mlcp_driver_init(&aMLCP, aMlcpOptions);
  aMlcpOptions->dparam[SICONOS_DPARAM_TOL] = options->dparam[SICONOS_DPARAM_TOL];
  *info = mlcp_driver(&aMLCP, reaction, velocity, aMlcpOptions);
  mlcp_driver_reset(&aMLCP, aMlcpOptions);
  solver_options_delete(aMlcpOptions);
}

void gmp_as_mlcp(GenericMechanicalProblem *problem, double *reaction, double *velocity,
                 int *info, SolverOptions *options) {
  /* FC3D blocks are not supported by the MLCP reformulation. */
  GMP_LocalProblem *local = problem->firstLocal;
  while (local) {
    switch (local->type) {
      case SICONOS_NUMERICS_PROBLEM_EQUALITY:
      case SICONOS_NUMERICS_PROBLEM_LCP:
        break;
      case SICONOS_NUMERICS_PROBLEM_FC3D: {
        (void)numerics_error("gmp_as_mlcp", "gmp_as_mlcp doesn't deal with FC3D");
        *info = 1;
        return;
      }
      default:
        (void)numerics_warning("gmp_as_mlcp", "unknown problem type");
    }
    local = local->next;
  }
  size_t Me_size;
  size_t Mi_size;

  SparseBlockStructuredMatrix *m = problem->M->matrix1;
  size_t numberOfRows = m->blocksize0[m->blocknumber0 - 1];
  size_t numberOfColumns = m->blocksize1[m->blocknumber1 - 1];

  double *reducedProb = (double *)malloc(numberOfRows * numberOfColumns * sizeof(double));
  double *Qreduced = (double *)malloc(numberOfRows * sizeof(double));
  if (!reducedProb || !Qreduced) {
    *info = 1;
    (void)numerics_error("gmp_as_mlcp", "memory allocation failed");
    free(reducedProb);
    free(Qreduced);
    return;
  }
  if (_GMPReducedEquality(problem, reducedProb, Qreduced, &Me_size, &Mi_size)) {
    *info = 1;
    free(reducedProb);
    free(Qreduced);
    return;
  }

  if (!Me_size) {
    gmp_as_mlcp_lcp(reducedProb, Qreduced, Mi_size, reaction, velocity, info);
  } else if (!Mi_size) {
    gmp_as_mlcp_linear_system(reducedProb, Qreduced, Me_size, reaction, info);
  } else {
    gmp_as_mlcp_mlcp(reducedProb, Qreduced, Me_size, Mi_size, reaction, velocity, info, options);
  }

  free(reducedProb);
  free(Qreduced);
}
