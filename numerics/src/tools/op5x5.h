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

/*!\file op5x5.h
 * \brief linear algebra operations in 5D (5x5 matrices and 5D vectors)*/

#ifndef _op5x5_h_
#define _op5x5_h_

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288   /* pi             */
#define M_PI_2 1.57079632679489661923132169163975144 /* pi/2           */
#endif

#ifndef MAYBE_UNUSED
#ifdef __GNUC__
#define MAYBE_UNUSED __attribute__((unused))
#else
#define MAYBE_UNUSED
#endif
#endif

#ifndef WARN_RESULT_IGNORED
#ifdef __GNUC__
#define WARN_RESULT_IGNORED __attribute__((warn_unused_result))
#else
#define WARN_RESULT_IGNORED
#endif
#endif

#ifdef __cplusplus
#undef restrict
#include <sys/cdefs.h>  // for __restrict
#define restrict __restrict
#endif

/** OP5X5(EXPR) do EXPR 25 times
 * \param EXPR a C expression that should contains self incrementing
 *        pointers on arrays[25] */
#define OP5X5(EXPR) \
  do {              \
    EXPR; EXPR; EXPR; EXPR;  EXPR; \
    EXPR; EXPR; EXPR; EXPR;  EXPR; \
    EXPR; EXPR; EXPR; EXPR;  EXPR; \
    EXPR; EXPR; EXPR; EXPR;  EXPR; \
    EXPR; EXPR; EXPR; EXPR;  EXPR; \
  } while (0)

/** OP5(EXPR) do EXPR 5 times
 * \param EXPR a C expression that should contains self incrementing
 *        pointers on arrays[5] */
#define OP5(EXPR) \
  do {            \
    EXPR; EXPR; EXPR; EXPR;  EXPR; \
  } while (0)

/** SET5X5(V) set a 5x5 matrix to V column-major order
 * \param V a double[25] */
#define SET5X5(V)                   \
  double* V##00 MAYBE_UNUSED = V++; \
  double* V##10 MAYBE_UNUSED = V++; \
  double* V##20 MAYBE_UNUSED = V++; \
  double* V##30 MAYBE_UNUSED = V++; \
  double* V##40 MAYBE_UNUSED = V++; \
  double* V##01 MAYBE_UNUSED = V++; \
  double* V##11 MAYBE_UNUSED = V++; \
  double* V##21 MAYBE_UNUSED = V++; \
  double* V##31 MAYBE_UNUSED = V++; \
  double* V##41 MAYBE_UNUSED = V++; \
  double* V##02 MAYBE_UNUSED = V++; \
  double* V##12 MAYBE_UNUSED = V++; \
  double* V##22 MAYBE_UNUSED = V++; \
  double* V##32 MAYBE_UNUSED = V++; \
  double* V##42 MAYBE_UNUSED = V++; \
  double* V##03 MAYBE_UNUSED = V++; \
  double* V##13 MAYBE_UNUSED = V++; \
  double* V##23 MAYBE_UNUSED = V++; \
  double* V##33 MAYBE_UNUSED = V++; \
  double* V##43 MAYBE_UNUSED = V++; \
  double* V##04 MAYBE_UNUSED = V++; \
  double* V##14 MAYBE_UNUSED = V++; \
  double* V##24 MAYBE_UNUSED = V++; \
  double* V##34 MAYBE_UNUSED = V++; \
  double* V##44 MAYBE_UNUSED = V++;

/** SET5(V) set a 5D vector
 * \param V a double[5] */
#define SET5(V)                    \
  double* V##0 MAYBE_UNUSED = V++; \
  double* V##1 MAYBE_UNUSED = V++; \
  double* V##2 MAYBE_UNUSED = V++; \
  double* V##3 MAYBE_UNUSED = V++; \
  double* V##4 MAYBE_UNUSED = V++;

/** SET5MAYBE : set pointers on a vector5 v (*v0 *v1 *v2 *v3 *v4) only if v is
 * non null.
 * Warning: the pointer v is modified and is ready for a next SET5
 * use *v0 if you need *v
 */
#define SET5MAYBE(V)             \
  double* V##0 MAYBE_UNUSED = 0; \
  double* V##1 MAYBE_UNUSED = 0; \
  double* V##2 MAYBE_UNUSED = 0; \
  double* V##3 MAYBE_UNUSED = 0; \
  double* V##4 MAYBE_UNUSED = 0; \
  if (V) {                       \
    V##0 = V++;                  \
    V##1 = V++;                  \
    V##2 = V++;                  \
    V##3 = V++;                  \
    V##4 = V++;                  \
  }

/*============= 5x5 matrices and 5D vectors =============*/

/** copy a 5x5 matrix or a vector[25] */
static inline void cpy5x5(double* restrict a, double* restrict b) { OP5X5(*b++ = *a++); }

/** add a 5x5 matrix or a vector[9] */
static inline void add5x5(double a[25], double b[25]) { OP5X5(*b++ += *a++); }

/** sub a 5x5 matrix or a vector[25] */
static inline void sub5x5(double a[25], double b[25]) { OP5X5(*b++ -= *a++); }

/** copy a 5D vector */
static inline void cpy5(double a[5], double b[5]) { OP5(*b++ = *a++); }

/** add two 5D vectors */
static inline void add5(double a[5], double b[5]) { OP5(*b++ += *a++); }

/** scalar multiplication of a matrix5x5 */
static inline void scal5x5(double scal, double m[25]) { OP5X5(*m++ *= scal); }

/** set to zero a matrix5x5 */
static inline void zero5x5(double m[25]) { OP5X5(*m++ = 0.); }

/** set to identity a matrix5x5
 * \param[in,out] m  m[25]
 */
static inline void eye5x5(double m[25]) {
  SET5X5(m);
  *m00 = 1.;
  *m01 = 0.;
  *m02 = 0.;
  *m03 = 0.;
  *m04 = 0.;
  *m10 = 0.;
  *m11 = 1.;
  *m12 = 0.;
  *m13 = 0.;
  *m14 = 0.;
  *m20 = 0.;
  *m21 = 0.;
  *m22 = 1.;
  *m23 = 0.;
  *m24 = 0.;
  *m30 = 0.;
  *m31 = 0.;
  *m32 = 0.;
  *m33 = 1.;
  *m34 = 0.;
  *m40 = 0.;
  *m41 = 0.;
  *m42 = 0.;
  *m43 = 0.;
  *m44 = 1.;
}

/** scalar multiplication of a 5D vector */
static inline void scal5(double scal, double* v) { OP5(*v++ *= scal); }

static inline void display5(double* restrict a) {
  SET5(a);
  printf("[ %6.4e\t  %6.4e\t %6.4e\t  %6.4e\t %6.4e ]\n", *a0, *a1, *a2, *a3, *a4);
}
static inline void display5x5(double* restrict a) {
  SET5X5(a);
  printf("[ %6.4e\t  %6.4e\t %6.4e\t  %6.4e\t %6.4e ]\n", *a00, *a01, *a02, *a03, *a04);
  printf("[ %6.4e\t  %6.4e\t %6.4e\t  %6.4e\t %6.4e ]\n", *a10, *a11, *a12, *a13, *a14);
  printf("[ %6.4e\t  %6.4e\t %6.4e\t  %6.4e\t %6.4e ]\n", *a20, *a21, *a22, *a23, *a24);
  printf("[ %6.4e\t  %6.4e\t %6.4e\t  %6.4e\t %6.4e ]\n", *a30, *a31, *a32, *a33, *a34);
  printf("[ %6.4e\t  %6.4e\t %6.4e\t  %6.4e\t %6.4e ]\n", *a40, *a41, *a42, *a43, *a44);
}



static inline void mvp5x5(const double* restrict a, const double* restrict v,
                          double* restrict r) {
  double* pr;

  pr = r;

  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v++;

  pr = r;

  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v++;

  pr = r;

  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v++;

  pr = r;

  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v++;

  pr = r;

  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v++;
}



/** matrix matrix multiplication : c = a * b
 * \param[in] a  a[25]
 * \param[in] b  b[25]
 * \param[out] c c[25]
 */
static inline void mm5x5(double* restrict a, double* restrict b, double* restrict c) {
  SET5X5(a);
  SET5X5(b);
  SET5X5(c);

  *c00 = *a00 * *b00 + *a01 * *b10 + *a02 * *b20 + *a03 * *b30 + *a04 * *b40;
  *c01 = *a00 * *b01 + *a01 * *b11 + *a02 * *b21 + *a03 * *b31 + *a04 * *b41;
  *c02 = *a00 * *b02 + *a01 * *b12 + *a02 * *b22 + *a03 * *b32 + *a04 * *b42;
  *c03 = *a00 * *b03 + *a01 * *b13 + *a02 * *b23 + *a03 * *b33 + *a04 * *b43;
  *c04 = *a00 * *b04 + *a01 * *b14 + *a02 * *b24 + *a03 * *b34 + *a04 * *b44;

  *c10 = *a10 * *b00 + *a11 * *b10 + *a12 * *b20 + *a13 * *b30 + *a14 * *b40;
  *c11 = *a10 * *b01 + *a11 * *b11 + *a12 * *b21 + *a13 * *b31 + *a14 * *b41;
  *c12 = *a10 * *b02 + *a11 * *b12 + *a12 * *b22 + *a13 * *b32 + *a14 * *b42;
  *c13 = *a10 * *b03 + *a11 * *b13 + *a12 * *b23 + *a13 * *b33 + *a14 * *b43;
  *c14 = *a10 * *b04 + *a11 * *b14 + *a12 * *b24 + *a13 * *b34 + *a14 * *b44;

  *c20 = *a20 * *b00 + *a21 * *b10 + *a22 * *b20 + *a23 * *b30 + *a24 * *b40;
  *c21 = *a20 * *b01 + *a21 * *b11 + *a22 * *b21 + *a23 * *b31 + *a24 * *b41;
  *c22 = *a20 * *b02 + *a21 * *b12 + *a22 * *b22 + *a23 * *b32 + *a24 * *b42;
  *c23 = *a20 * *b03 + *a21 * *b13 + *a22 * *b23 + *a23 * *b33 + *a24 * *b43;
  *c24 = *a20 * *b04 + *a21 * *b14 + *a22 * *b24 + *a23 * *b34 + *a24 * *b44;

  *c30 = *a30 * *b00 + *a31 * *b10 + *a32 * *b20 + *a33 * *b30 + *a34 * *b40;
  *c31 = *a30 * *b01 + *a31 * *b11 + *a32 * *b21 + *a33 * *b31 + *a34 * *b41;
  *c32 = *a30 * *b02 + *a31 * *b12 + *a32 * *b22 + *a33 * *b32 + *a34 * *b42;
  *c33 = *a30 * *b03 + *a31 * *b13 + *a32 * *b23 + *a33 * *b33 + *a34 * *b43;
  *c34 = *a30 * *b04 + *a31 * *b14 + *a32 * *b24 + *a33 * *b34 + *a34 * *b44;

  *c40 = *a40 * *b00 + *a41 * *b10 + *a42 * *b20 + *a43 * *b30 + *a44 * *b40;
  *c41 = *a40 * *b01 + *a41 * *b11 + *a42 * *b21 + *a43 * *b31 + *a44 * *b41;
  *c42 = *a40 * *b02 + *a41 * *b12 + *a42 * *b22 + *a43 * *b32 + *a44 * *b42;
  *c43 = *a40 * *b03 + *a41 * *b13 + *a42 * *b23 + *a43 * *b33 + *a44 * *b43;
  *c44 = *a40 * *b04 + *a41 * *b14 + *a42 * *b24 + *a43 * *b34 + *a44 * *b44;
}



static inline double hypot5(double* a) {
  double r;

  r = *a * *a;
  a++;
  r += *a * *a;
  a++;
  r += *a * *a;
  a++;
  r += *a * *a;
  a++;
  r += *a * *a;

  return sqrt(r);
}
static inline void zero5(double* v) { OP5(*v++ = 0.); }


/*============= 5x5 matrices and 5D vectors =============*/

/** copy a 5x5 matrix or a vector[25] */
static inline void cpy5x5_for_loop(double* restrict a, double* restrict b) {
  for (int i = 0; i < 25; ++i) b[i] = a[i];
}

/** add a 5x5 matrix or a vector[25] */
static inline void add5x5_for_loop(double a[25], double b[25]) {
  for (int i = 0; i < 25; ++i) b[i] += a[i];
}

/** sub a 5x5 matrix or a vector[25] */
static inline void sub5x5_for_loop(double a[25], double b[25]) {
  for (int i = 0; i < 25; ++i) b[i] -= a[i];
}

/** copy a 5D vector */
static inline void cpy5_for_loop(double a[5], double b[5]) {
  for (int i = 0; i < 5; ++i) b[i] = a[i];
}

/** add two 5D vectors */
static inline void add5_for_loop(double a[5], double b[5]) {
  for (int i = 0; i < 5; ++i) b[i] += a[i];
}

/** scalar multiplication of a 5D vector */
static inline void scal5_for_loop(double scal, double* v) {
  for (int i = 0; i < 5; ++i) v[i] *= scal;
}

/** scalar multiplication of a matrix5x5 */
static inline void scal5x5_for_loop(double scal, double m[25]) {
  for (int i = 0; i < 25; ++i) m[i] *= scal;
}

/** set to zero a matrix5x5 */
static inline void zero5x5_for_loop(double m[25]) {
  for (int i = 0; i < 25; ++i) m[i] = 0.0;
}

/** set to identity a matrix5x5
 * \param[in,out] m  m[25]
 */
static inline void eye5x5_for_loop(double m[25]) {
  zero5x5_for_loop(m);
  m[0] = 1.0;
  m[6] = 1.0;
  m[12] = 1.0;
  m[18] = 1.0;
  m[24] = 1.0;
}

/** set to zero a 5D vector */
static inline void zero5_for_loop(double* v) {
  for (int i = 0; i < 5; ++i) v[i] = 0.0;
}

/** display a 5D vector */
static inline void display5_for_loop(double* restrict a) {
  printf("%12.5f  %12.5f  %12.5f  %12.5f  %12.5f \n", a[0], a[1], a[2], a[3], a[4]);
}

/** display a 5x5 matrix (column-major storage, displayed in row-major format) */
static inline void display5x5_for_loop(double* restrict a) {
  printf("\n");
  for (int i = 0; i < 5; ++i) {
    printf("%12.5f  %12.5f  %12.5f  %12.5f  %12.5f \n", a[0 * 5 + i], a[1 * 5 + i],
           a[2 * 5 + i], a[3 * 5 + i], a[4 * 5 + i]);
  }
  printf("\n");
}

/* matrix-vector product of a 5x5 matrix (column-major) and a 5D vector */
static inline void mvp5x5_for_loop(const double* restrict a, const double* restrict v,
                                   double* restrict w) {
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      w[i] += a[j * 5 + i] * v[j]; /* column-major: a[j][i] = a[j*5+i] */
    }
  }
}

/* matrix-vector product of a 5x5 matrix (transposed) and a 5D vector */
static inline void mvp5x5T_for_loop(const double* restrict a, const double* restrict v,
                                    double* restrict w) {
  for (int i = 0; i < 5; ++i) {
    w[i] = 0;
    for (int j = 0; j < 5; ++j) {
      w[i] += a[j * 5 + i] * v[j];
    }
  }
}

/** matrix matrix multiplication : c = a * b
 * \param[in] a  a[25]
 * \param[in] b  b[25]
 * \param[out] c c[25]
 */
static inline void mm5x5_for_loop(double* restrict a, double* restrict b, double* restrict c) {
  /* c = a * b */
  /* a, b, c are 5x5 column-major: a[row][col] = a[col*5+row] */
  for (int i = 0; i < 5; ++i) {   /* row i of result */
    for (int j = 0; j < 5; ++j) { /* column j of result */
      c[j * 5 + i] = 0;
      for (int k = 0; k < 5; ++k) { /* sum over k: a[i][k] * b[k][j] */
        c[j * 5 + i] += a[k * 5 + i] * b[j * 5 + k];
      }
    }
  }
}

/* multiply a and b (column-major), both 5x5 matrices, and add to c */
static inline void mm5x5_acc_for_loop(double* restrict a, double* restrict b,
                                      double* restrict c) {
  /* c += a * b */
  /* a, b, c are 5x5 column-major: a[row][col] = a[col*5+row] */
  for (int i = 0; i < 5; ++i) {     /* row i of result */
    for (int j = 0; j < 5; ++j) {   /* column j of result */
      for (int k = 0; k < 5; ++k) { /* sum over k: a[i][k] * b[k][j] */
        c[j * 5 + i] += a[k * 5 + i] * b[j * 5 + k];
      }
    }
  }
}

/* multiply a (5x5) and b (5x5) column-major, subtract result from c */
static inline void mm5x5_sub_for_loop(double* restrict a, double* restrict b,
                                      double* restrict c) {
  /* c -= a * b */
  /* a, b, c are 5x5 column-major: a[row][col] = a[col*5+row] */
  for (int i = 0; i < 5; ++i) {     /* row i of result */
    for (int j = 0; j < 5; ++j) {   /* column j of result */
      for (int k = 0; k < 5; ++k) { /* sum over k: a[i][k] * b[k][j] */
        c[j * 5 + i] -= a[k * 5 + i] * b[j * 5 + k];
      }
    }
  }
}

/* c = a^T * b (column-major) */
static inline void mTm5x5_for_loop(double* restrict a, double* restrict b,
                                   double* restrict c) {
  /* a, b, c are 5x5 column-major */
  /* (a^T)[i][k] = a[k][i], so c[i][j] = sum_k a[k][i] * b[k][j] */
  for (int i = 0; i < 5; ++i) {   /* row i of result (col i of a) */
    for (int j = 0; j < 5; ++j) { /* column j of result */
      c[j * 5 + i] = 0;
      for (int k = 0; k < 5; ++k) { /* sum over k */
        c[j * 5 + i] += a[i * 5 + k] * b[j * 5 + k];
      }
    }
  }
}

/* a <- a^T, swap in place the non diagonal entries (column-major) */
static inline void transpose5x5_for_loop(double* a) {
  /* In column-major: a[i][j] = a[j*5+i], a[j][i] = a[i*5+j] */
  double tmp;
  for (int i = 0; i < 5; ++i) {
    for (int j = i + 1; j < 5; ++j) {
      tmp = a[j * 5 + i];          /* a[i][j] */
      a[j * 5 + i] = a[i * 5 + j]; /* swap with a[j][i] */
      a[i * 5 + j] = tmp;
    }
  }
}

/* c = a * b^T (column-major) */
static inline void mmt5x5_for_loop(double* restrict a, double* restrict b,
                                   double* restrict c) {
  /* a, b, c are 5x5 column-major */
  /* c[i][j] = sum_k a[i][k] * b[j][k] (since (b^T)[k][j] = b[j][k]) */
  for (int i = 0; i < 5; ++i) {   /* row i of result */
    for (int j = 0; j < 5; ++j) { /* column j of result */
      c[j * 5 + i] = 0;
      for (int k = 0; k < 5; ++k) { /* sum over k */
        c[j * 5 + i] += a[k * 5 + i] * b[k * 5 + j];
      }
    }
  }
}

/* 5x5 matrix - 5x5 matrix hadamard (element-wise) product */
static inline void hadamard5x5_for_loop(double* restrict a, double* restrict b,
                                        double* restrict c) {
  for (int i = 0; i < 25; ++i) {
    c[i] = a[i] * b[i];
  }
}

/* a <- a + b, with b a 5x5 matrix */
static inline void add5x5_to_for_loop(double* a, double* b) {
  for (int i = 0; i < 25; ++i) {
    a[i] += b[i];
  }
}

static inline double trace5x5_for_loop(const double* a) {
  return a[0] + a[6] + a[12] + a[18] + a[24];
}

/* 5D vector dot product */
static inline double dot5_for_loop(double* a, double* b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3] + a[4] * b[4];
}

/* compute the Euclidean norm of a 5D vector */
static inline double norm5_for_loop(double* a) {
  return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2] + a[3] * a[3] + a[4] * a[4]);
}

/* compute the squared Euclidean norm of a 5D vector */
static inline double norm2_5_for_loop(double* a) {
  return a[0] * a[0] + a[1] * a[1] + a[2] * a[2] + a[3] * a[3] + a[4] * a[4];
}

/* compute the Euclidean norm of a 5D vector, in place and fast */
static inline void normalize5_for_loop(double* a) {
  double norm = norm5_for_loop(a);
  a[0] /= norm;
  a[1] /= norm;
  a[2] /= norm;
  a[3] /= norm;
  a[4] /= norm;
}

static inline void axpy5_for_loop(double a, double* x, double* y) {
  /* y <- a*x + y */
  for (int i = 0; i < 5; ++i) {
    y[i] = a * x[i] + y[i];
  }
}

static inline void axpby5_for_loop(double a, double b, double* x, double* y) {
  /* y <- a*x + b*y */
  for (int i = 0; i < 5; ++i) {
    y[i] = a * x[i] + b * y[i];
  }
}

/* symmetric part of a 5x5 matrix (column-major): (a + a^T)/2 */
static inline void sym5x5_for_loop(double* a, double* sym) {
  /* a[i][j] = a[j*5+i], a[j][i] = a[i*5+j] */
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      sym[j * 5 + i] = 0.5 * (a[j * 5 + i] + a[i * 5 + j]);
    }
  }
}

/* skew-symmetric part of a 5x5 matrix (column-major): (a - a^T)/2 */
static inline void skew5x5_for_loop(double* a, double* skew) {
  /* a[i][j] = a[j*5+i], a[j][i] = a[i*5+j] */
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      skew[j * 5 + i] = 0.5 * (a[j * 5 + i] - a[i * 5 + j]);
    }
  }
}

/* 5x5 matrix-vector product (column-major), result added to existing vector: y += a * x */
static inline void mvp5x5_add_for_loop(double* a, double* x, double* y) {
  /* a[row][col] = a[col*5+row] */
  for (int i = 0; i < 5; ++i) {   /* row i */
    for (int j = 0; j < 5; ++j) { /* column j */
      y[i] += a[j * 5 + i] * x[j];
    }
  }
}

/* Infinity norm (max row sum) for column-major 5x5 matrix */
static inline double infnorm5x5_for_loop(double* a) {
  /* a[row][col] = a[col*5+row] */
  double max = 0;
  for (int i = 0; i < 5; ++i) { /* for each row i */
    double row_sum = 0;
    for (int j = 0; j < 5; ++j) { /* sum across columns */
      row_sum += fabs(a[j * 5 + i]);
    }
    if (row_sum > max) max = row_sum;
  }
  return max;
}

static inline double infnorm5_for_loop(double* v) {
  double max = fabs(v[0]);
  for (int i = 1; i < 5; ++i) {
    if (fabs(v[i]) > max) max = fabs(v[i]);
  }
  return max;
}

/* Column sum norm for column-major 5x5 matrix */
static inline void column_sum_norm5x5_for_loop(double* a, double* norm) {
  /* a[row][col] = a[col*5+row], so column j is contiguous */
  for (int j = 0; j < 5; ++j) {
    norm[j] = 0;
    for (int i = 0; i < 5; ++i) {
      norm[j] += fabs(a[j * 5 + i]);
    }
  }
}

static inline double hypot5_for_loop(double* a) {
  return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2] + a[3] * a[3] + a[4] * a[4]);
}

/* 5x5 matrix solve - Gaussian elimination with partial pivoting (column-major) */
static inline int solve_5x5_gepp_for_loop(const double* restrict a, double* restrict b) {
  double A[25];
  double x[5];
  for (int i = 0; i < 25; ++i) A[i] = a[i];
  for (int i = 0; i < 5; ++i) x[i] = b[i];

  int ipiv[5];
  for (int i = 0; i < 5; ++i) ipiv[i] = i;

  for (int i = 0; i < 5; ++i) {
    /* Find pivot - column i has elements at A[i*5 + row] in column-major */
    int max_row = i;
    double max_val = fabs(A[i * 5 + i]);
    for (int k = i + 1; k < 5; ++k) {
      if (fabs(A[i * 5 + k]) > max_val) {
        max_val = fabs(A[i * 5 + k]);
        max_row = k;
      }
    }

    if (max_val < 1e-16) {
      return i + 1; /* Singular matrix */
    }

    /* Swap rows in column-major: swap elements in each column */
    if (max_row != i) {
      for (int j = 0; j < 5; ++j) {
        double tmp = A[j * 5 + i];
        A[j * 5 + i] = A[j * 5 + max_row];
        A[j * 5 + max_row] = tmp;
      }
      double tmp = x[i];
      x[i] = x[max_row];
      x[max_row] = tmp;
    }

    /* Eliminate column */
    for (int k = i + 1; k < 5; ++k) {
      double factor = A[i * 5 + k] / A[i * 5 + i];
      for (int j = i; j < 5; ++j) {
        A[j * 5 + k] -= factor * A[j * 5 + i];
      }
      x[k] -= factor * x[i];
    }
  }

  /* Back substitution */
  for (int i = 4; i >= 0; --i) {
    double sum = 0;
    for (int j = i + 1; j < 5; ++j) {
      sum += A[j * 5 + i] * x[j];
    }
    x[i] = (x[i] - sum) / A[i * 5 + i];
  }

  for (int i = 0; i < 5; ++i) b[i] = x[i];
  return 0;
}

/* 5x5 matrix inverse (column-major) */
static inline int inv_5x5_gepp_for_loop(double* a) {
  double A[25];
  double B[25];

  /* Create augmented matrix [A|I] in column-major */
  for (int i = 0; i < 25; ++i) {
    A[i] = a[i];
    B[i] = 0.0;
  }
  /* Set identity in B (diagonal elements) */
  for (int i = 0; i < 5; ++i) {
    B[i * 5 + i] = 1.0; /* B[i][i] in column-major */
  }

  for (int i = 0; i < 5; ++i) {
    /* Find pivot - look down column i */
    int max_row = i;
    double max_val = fabs(A[i * 5 + i]);
    for (int k = i + 1; k < 5; ++k) {
      if (fabs(A[i * 5 + k]) > max_val) {
        max_val = fabs(A[i * 5 + k]);
        max_row = k;
      }
    }

    if (max_val < 1e-16) {
      return 1; /* Singular matrix */
    }

    /* Swap rows in column-major: swap elements in each column */
    if (max_row != i) {
      for (int j = 0; j < 5; ++j) {
        double tmp = A[j * 5 + i];
        A[j * 5 + i] = A[j * 5 + max_row];
        A[j * 5 + max_row] = tmp;
        tmp = B[j * 5 + i];
        B[j * 5 + i] = B[j * 5 + max_row];
        B[j * 5 + max_row] = tmp;
      }
    }

    /* Scale pivot row */
    double pivot = A[i * 5 + i];
    for (int j = 0; j < 5; ++j) {
      A[j * 5 + i] /= pivot;
      B[j * 5 + i] /= pivot;
    }

    /* Eliminate column */
    for (int k = 0; k < 5; ++k) {
      if (k != i) {
        double factor = A[i * 5 + k];
        for (int j = 0; j < 5; ++j) {
          A[j * 5 + k] -= factor * A[j * 5 + i];
          B[j * 5 + k] -= factor * B[j * 5 + i];
        }
      }
    }
  }

  /* Copy inverse back */
  for (int i = 0; i < 25; ++i) a[i] = B[i];
  return 0;
}

#endif /* _op5x5_h_ */

/* NOTE: This file (op5x5.h) contains 5x5 matrix and 5D vector operations.
 * It does NOT include op3x3.h - include both headers if you need both:
 *
 * Example usage:
 *   #include "op3x3.h"       // for 3x3 matrices and 3D vectors
 *   #include "op5x5.h"       // for 5x5 matrices and 5D vectors
 */
