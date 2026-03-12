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

/*!\file test_op5x5.c
 * \brief Test file for 5x5 matrix operations in op5x5.h
 *
 * Matrices are stored in COLUMN-MAJOR order (like BLAS/LAPACK):
 * For a 5x5 matrix A where A[row][col]:
 *   A[0][0] = a[0],  A[1][0] = a[1],  A[2][0] = a[2],  A[3][0] = a[3],  A[4][0] = a[4] (column
 * 0) A[0][1] = a[5],  A[1][1] = a[6],  A[2][1] = a[7],  A[3][1] = a[8],  A[4][1] = a[9]
 * (column 1) A[0][2] = a[10], A[1][2] = a[11], A[2][2] = a[12], A[3][2] = a[13], A[4][2] =
 * a[14] (column 2) A[0][3] = a[15], A[1][3] = a[16], A[2][3] = a[17], A[3][3] = a[18], A[4][3]
 * = a[19] (column 3) A[0][4] = a[20], A[1][4] = a[21], A[2][4] = a[22], A[3][4] = a[23],
 * A[4][4] = a[24] (column 4)
 *
 * In memory: a[0]=A[0][0], a[1]=A[1][0], a[2]=A[2][0], a[3]=A[3][0], a[4]=A[4][0],
 *            a[5]=A[0][1], a[6]=A[1][1], ...
 */

#include "CSparseMatrix.h"
#define _XOPEN_SOURCE 700

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "SiconosBlas.h"
#include "SiconosLapack.h"
#include "op5x5.h"

/* Helper: check if two 5x5 matrices are equal (column-major) */
static int equal5x5(double* a, double* b) {
  for (int i = 0; i < 25; i++) {
    if (fabs(a[i] - b[i]) > 1e-10) {
      return 0;
    }
  }
  return 1;
}

/* Helper: check if two 5D vectors are equal */
static int equal5(double* a, double* b) {
  for (int i = 0; i < 5; i++) {
    if (fabs(a[i] - b[i]) > 1e-10) {
      return 0;
    }
  }
  return 1;
}

/* Helper: print 5x5 matrix in human-readable row-major format */
static void print5x5_rm(double* a) {
  for (int row = 0; row < 5; row++) {
    printf("[");
    for (int col = 0; col < 5; col++) {
      printf(" %6.2f", a[col * 5 + row]);  // Convert column-major to row-major for display
    }
    printf(" ]\n");
  }
}

int main() {
  printf("Testing op5x5.h - 5x5 matrix operations (column-major)\n");
  printf("==========================================================\n\n");
  double FOO[25] = {
      1,  2,  3,  4,  5,   // column 0: A[0][0]=1, A[1][0]=2, A[2][0]=3, A[3][0]=4, A[4][0]=5
      6,  7,  8,  9,  10,  // column 1: A[0][1]=6, A[1][1]=7, A[2][1]=8, A[3][1]=9, A[4][1]=10
      11, 12, 13, 14, 15,  // column 2
      16, 17, 18, 19, 20,  // column 3
      21, 22, 23, 24, 25   // column 4
  };
  /* Test 1: Identity matrix */
  printf("Test 1: eye5x5\n");
  double _I[25];
  eye5x5_for_loop(_I);
  // Column-major: diagonal at positions 0, 6, 12, 18, 24
  double I_expected[25] = {
      1, 0, 0, 0, 0,  // column 0: A[0][0]=1, A[1][0]=0, A[2][0]=0, A[3][0]=0, A[4][0]=0
      0, 1, 0, 0, 0,  // column 1: A[0][1]=0, A[1][1]=1, A[2][1]=0, A[3][1]=0, A[4][1]=0
      0, 0, 1, 0, 0,  // column 2
      0, 0, 0, 1, 0,  // column 3
      0, 0, 0, 0, 1   // column 4
  };
  cpy5x5_for_loop(FOO, _I);
  eye5x5(_I);
  assert(equal5x5(_I, I_expected));
  printf("  PASSED\n\n");

  /* Test 2: Zero matrix */

  double Z[25];
  cpy5x5_for_loop(FOO, Z);
  zero5x5_for_loop(Z);
  for (int i = 0; i < 25; i++) assert(Z[i] == 0.0);
  zero5x5(Z);
  for (int i = 0; i < 25; i++) assert(Z[i] == 0.0);
  printf("Test 2: zero5x5\n");
  printf("  PASSED\n\n");

  /* Test 3: Copy matrix (cpy5x5) */
  printf("Test 3: cpy5x5\n");
  // Column-major: matrix with elements 1,2,3,...25 in column-major order
  double A[25] = {
      1,  2,  3,  4,  5,   // column 0: A[0][0]=1, A[1][0]=2, A[2][0]=3, A[3][0]=4, A[4][0]=5
      6,  7,  8,  9,  10,  // column 1: A[0][1]=6, A[1][1]=7, A[2][1]=8, A[3][1]=9, A[4][1]=10
      11, 12, 13, 14, 15,  // column 2
      16, 17, 18, 19, 20,  // column 3
      21, 22, 23, 24, 25   // column 4
  };
  double B[25];
  cpy5x5_for_loop(A, B);
  assert(equal5x5(A, B));
  cpy5x5_for_loop(FOO, B);
  cpy5x5(A, B);
  assert(equal5x5(A, B));
  printf("  PASSED\n\n");

  /* Test 4: Matrix addition (add5x5) */
  printf("Test 4: add5x5\n");
  double C[25];
  cpy5x5_for_loop(A, C);
  add5x5_for_loop(B, C);  // C = A + B = A + A = 2*A
  for (int i = 0; i < 25; i++) assert(fabs(C[i] - 2 * A[i]) < 1e-10);

  cpy5x5(FOO, C);
  cpy5x5(A, C);
  add5x5(B, C);  // C = A + B = A + A = 2*A
  for (int i = 0; i < 25; i++) assert(fabs(C[i] - 2 * A[i]) < 1e-10);

  printf("  PASSED\n\n");

  /* Test 5: Matrix subtraction (sub5x5) */
  printf("Test 5: sub5x5\n");
  sub5x5_for_loop(C, C);  // C = C - C = 0
  assert(equal5x5(C, Z));
  cpy5x5(FOO, C);
  sub5x5(C, C);  // C = C - C = 0
  assert(equal5x5(C, Z));
  printf("  PASSED\n\n");

  /* Test 6: Scalar multiplication (scal5x5) */
  printf("Test 6: scal5x5\n");
  double D[25];
  cpy5x5_for_loop(A, D);
  scal5x5_for_loop(2.0, D);
  for (int i = 0; i < 25; i++) assert(fabs(D[i] - 2 * A[i]) < 1e-10);

  cpy5x5(FOO, D);
  cpy5x5(A, D);
  scal5x5(2.0, D);
  for (int i = 0; i < 25; i++) assert(fabs(D[i] - 2 * A[i]) < 1e-10);
  printf("  PASSED\n\n");

  /* Test 7: Matrix-vector product (mvp5x5) */
  printf("Test 7: mvp5x5\n");
  // w = A * v + w with A = I
  double v[5] = {1, 2, 3, 4, 5};
  double w[5] = {1, 2, 3, 4, 5};
  double test[5] = {2, 4, 6, 8, 10};

  eye5x5_for_loop(_I);
  mvp5x5_for_loop(_I, v, w);  // _I * v + w = w
  assert(equal5(w, test));

  double v_[5] = {1, 2, 3, 4, 5};
  double w_[5] = {1, 2, 3, 4, 5};
  mvp5x5(_I, v_, w_);  // _I * v + w = w
  assert(equal5(test, w_));

  printf("  PASSED\n\n");

  /* Test 8: Matrix-matrix multiplication (mm5x5) */
  printf("Test 8: mm5x5\n");
  double E[25];
  eye5x5_for_loop(_I);
  mm5x5_for_loop(_I, A, E);  // E = _I * A = A
  assert(equal5x5(A, E));

  mm5x5_for_loop(A, _I, E);  // E = A * _I = A
  assert(equal5x5(A, E));

  cpy5x5(FOO, E);

  mm5x5(_I, A, E);  // E = _I * A = A
  assert(equal5x5(A, E));

  mm5x5(A, _I, E);  // E = A * _I = A
  assert(equal5x5(A, E));

  printf("  PASSED\n\n");

  /* Test 9: Transpose (transpose5x5) */
  printf("Test 9: transpose5x5\n");
  // Create symmetric matrix in column-major
  double Sym[25] = {
      1, 2, 3,  4,  5,   // column 0
      2, 6, 7,  8,  9,   // column 1: symmetric with column 0
      3, 7, 10, 11, 12,  // column 2
      4, 8, 11, 13, 14,  // column 3
      5, 9, 12, 14, 15   // column 4
  };
  double SymT[25];
  cpy5x5_for_loop(Sym, SymT);
  transpose5x5_for_loop(SymT);  // Should be unchanged (symmetric)
  assert(equal5x5(Sym, SymT));
  printf("  PASSED\n\n");

  /* Test 10: Vector dot product (dot5) */
  printf("Test 10: dot5\n");
  double u[5] = {1, 2, 3, 4, 5};
  double v_dot[5] = {1, 2, 3, 4, 5};
  double dot = dot5_for_loop(u, v_dot);  // 1+4+9+16+25 = 55
  assert(fabs(dot - 55.0) < 1e-10);
  printf("  PASSED\n\n");

  /* Test 11: Vector norm (norm5) */
  printf("Test 11: norm5\n");
  double n = norm5_for_loop(u);  // sqrt(55)
  assert(fabs(n - sqrt(55.0)) < 1e-10);
  printf("  PASSED\n\n");

  /* Test 12: Normalize vector (normalize5) */
  printf("Test 12: normalize5\n");
  double u2[5];
  cpy5_for_loop(u, u2);
  normalize5_for_loop(u2);
  n = norm5_for_loop(u2);
  assert(fabs(n - 1.0) < 1e-10);
  printf("  PASSED\n\n");

  /* Test 13: Solve linear system (solve_5x5_gepp) */
  printf("Test 13: solve_5x5_gepp (identity matrix)\n");
  double A_solve[25];
  eye5x5_for_loop(A_solve);
  double b[5] = {1, 2, 3, 4, 5};
  double x[5];
  cpy5_for_loop(b, x);
  int info = solve_5x5_gepp_for_loop(A_solve, x);
  assert(info == 0);
  assert(equal5(x, b));
  printf("  PASSED\n\n");

  /* Test 14: Solve with diagonal matrix */
  printf("Test 14: solve_5x5_gepp (diagonal matrix)\n");
  // A = 2*I, solve A*x = b -> x = b/2
  for (int i = 0; i < 5; i++) A_solve[i * 5 + i] = 2.0;
  double b2[5] = {2, 4, 6, 8, 10};
  cpy5_for_loop(b2, x);
  info = solve_5x5_gepp_for_loop(A_solve, x);
  double x_expected[5] = {1, 2, 3, 4, 5};
  assert(info == 0);
  assert(equal5(x, x_expected));
  printf("  PASSED\n\n");

  /* Test 15: Matrix inverse (inv_5x5_gepp) */
  printf("Test 15: inv_5x5_gepp\n");
  double A_inv[25];
  eye5x5_for_loop(A_solve);  // A = _I
  cpy5x5_for_loop(A_solve, A_inv);
  info = inv_5x5_gepp_for_loop(A_inv);
  assert(info == 0);
  assert(equal5x5(A_inv, I_expected));
  printf("  PASSED\n\n");

  /* Test 16: Trace (trace5x5) */
  printf("Test 16: trace5x5\n");
  // A has diagonal: 1, 7, 13, 19, 25 -> sum = 65
  double tr = trace5x5_for_loop(A);
  assert(fabs(tr - 65.0) < 1e-10);
  printf("  PASSED\n\n");

  /* Test 17: AXPY (axpy5) */
  printf("Test 17: axpy5\n");
  double y[5] = {1, 1, 1, 1, 1};
  double y_expected[5] = {3, 5, 7, 9, 11};  // y = 2*u + y
  axpy5_for_loop(2.0, u, y);
  assert(equal5(y, y_expected));
  printf("  PASSED\n\n");

  /* Test 18: Display functions */
  printf("Test 18: display5 and display5x5\n");
  printf("  Vector u: ");
  display5_for_loop(u);
  printf("  Matrix A (in row-major format):\n");
  print5x5_rm(A);
  printf("  PASSED\n\n");

  /* Test 19: Eigenvalues of diagonal matrix */
  printf("Test 19: eigvals_sym5x5_for_loop (diagonal matrix)\n");
  {
    double A_diag[25] = {0};
    A_diag[0] = 5.0;
    A_diag[6] = 3.0;
    A_diag[12] = 8.0;
    A_diag[18] = 1.0;
    A_diag[24] = 6.0;
    double eigvals[5];
    int info = eigvals_sym5x5_for_loop(A_diag, eigvals);
    double expected[5] = {1.0, 3.0, 5.0, 6.0, 8.0};
    assert(info == 0);
    assert(equal5(eigvals, expected));
    printf("  PASSED\n\n");
  }

  /* Test 20: Eigenvalues of 2x2 block matrix */
  printf("Test 20: eigvals_sym5x5_for_loop (2x2 block)\n");
  {
    double A_block[25] = {0};
    A_block[0] = 4.0;
    A_block[1] = 1.0;
    A_block[5] = 1.0;
    A_block[6] = 4.0;
    A_block[12] = 2.0;
    A_block[18] = 3.0;
    A_block[24] = 5.0;
    double eigvals[5];
    int info = eigvals_sym5x5_for_loop(A_block, eigvals);
    double expected[5] = {2.0, 3.0, 3.0, 5.0, 5.0};
    assert(info == 0);
    for (int i = 0; i < 5; i++) {
      assert(fabs(eigvals[i] - expected[i]) < 1e-8);
    }
    printf("  PASSED\n\n");
  }

  /* Test 21: Eigenvalues and eigenvectors */
  printf("Test 21: eigsym5x5_for_loop (eigenvectors)\n");
  {
    double A_full[25] = {2, 1, 0, 0, 0, 1, 2, 1, 0, 0, 0, 1, 2,
                         1, 0, 0, 0, 1, 2, 1, 0, 0, 0, 1, 2};
    double eigvals[5];
    double eigvecs[25];
    int info = eigsym5x5_for_loop(A_full, eigvals, eigvecs);
    assert(info == 0);

    /* Verify A*v = lambda*v for each eigenpair */
    for (int k = 0; k < 5; ++k) {
      double Av[5] = {0};
      double v[5];
      for (int i = 0; i < 5; ++i) v[i] = eigvecs[k * 5 + i];

      /* Compute A*v */
      for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
          Av[i] += A_full[j * 5 + i] * v[j];
        }
      }

      /* Check A*v = lambda*v */
      double error = 0;
      for (int i = 0; i < 5; ++i) {
        error += fabs(Av[i] - eigvals[k] * v[i]);
      }
      assert(error < 1e-8);
    }
    printf("  PASSED\n\n");
  }

  /* Test 22: Trace equals sum of eigenvalues */
  printf("Test 22: Trace equals sum of eigenvalues\n");
  {
    double A_test[25] = {5, 2, 1, 0, 0, 2, 4, 1, 1, 0, 1, 1, 6,
                         2, 1, 0, 1, 2, 3, 1, 0, 0, 1, 1, 4};
    double trace = A_test[0] + A_test[6] + A_test[12] + A_test[18] + A_test[24];
    double eigvals[5];
    int info = eigvals_sym5x5_for_loop(A_test, eigvals);
    double sum_eigs = eigvals[0] + eigvals[1] + eigvals[2] + eigvals[3] + eigvals[4];
    assert(info == 0);
    assert(fabs(trace - sum_eigs) < 1e-8);
    printf("  PASSED\n\n");
  }

  printf("==========================================================\n");
  printf("All tests PASSED!\n");

  return 0;
}
