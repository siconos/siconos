/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

/*!\file op3x3.h
 * \brief linear algebra operations in 3D*/

#ifndef _op3x3_h_
#define _op3x3_h_

#include <math.h>
#include <float.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288   /* pi             */
#define M_PI_2      1.57079632679489661923132169163975144   /* pi/2           */
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
#define WARN_RESULT_IGNORED __attribute__ ((warn_unused_result))
#else
#define WARN_RESULT_IGNORED
#endif
#endif

#ifdef __cplusplus
#undef restrict
#include <sys/cdefs.h>  // for __restrict
#define restrict __restrict
#endif

#define SOLVE_3X3_GEPP(MAT, X) \
{ \
  int res = solve_3x3_gepp(mat, X); \
  if (res) { printf("%s :: Call solve_3x3_gepp(" #MAT ", " #X " failed!. "\
                    "No pivot could be selected for column %d\n", __func__, res); } }


/** OP3X3(EXPR) do EXPR 9 times
 * \param EXPR a C expression that should contains self incrementing
 *        pointers on arrays[9] */
#define OP3X3(EXPR)                             \
  do                                            \
  {                                             \
    EXPR;                                       \
    EXPR;                                       \
    EXPR;                                       \
    EXPR;                                       \
    EXPR;                                       \
    EXPR;                                       \
    EXPR;                                       \
    EXPR;                                       \
    EXPR;                                       \
  } while(0)                                    \

/** OP3(EXPR) do EXPR 3 times
 * \param EXPR a C expression that should contains self incrementing
 *        pointers on arrays[9] */
#define OP3(EXPR)                               \
  do                                            \
  {                                             \
    EXPR;                                       \
    EXPR;                                       \
    EXPR;                                       \
  } while(0)                                    \

/** SET3X3 : set pointers on a 3x3 matrix a (*a00 *a01 *a10 etc.)
 * warning the pointer a is modified (use a00 instead) and is ready
 * for a next SET3X3
 */
#define SET3X3(V)                                                       \
  double* V##00 MAYBE_UNUSED = V++;                                     \
  double* V##10 MAYBE_UNUSED = V++;                                     \
  double* V##20 MAYBE_UNUSED = V++;                                     \
  double* V##01 MAYBE_UNUSED = V++;                                     \
  double* V##11 MAYBE_UNUSED = V++;                                     \
  double* V##21 MAYBE_UNUSED = V++;                                     \
  double* V##02 MAYBE_UNUSED = V++;                                     \
  double* V##12 MAYBE_UNUSED = V++;                                     \
  double* V##22 MAYBE_UNUSED = V++;
#define SET3X3MAYBE(V)                                       \
  double* V##00 MAYBE_UNUSED = 0;                            \
  double* V##10 MAYBE_UNUSED = 0;                            \
  double* V##20 MAYBE_UNUSED = 0;                            \
  double* V##01 MAYBE_UNUSED = 0;                            \
  double* V##11 MAYBE_UNUSED = 0;                            \
  double* V##21 MAYBE_UNUSED = 0;                            \
  double* V##02 MAYBE_UNUSED = 0;                            \
  double* V##12 MAYBE_UNUSED = 0;                            \
  double* V##22 MAYBE_UNUSED = 0;                            \
  if (V)                                        \
  {                                             \
    V##00 = V++;                                \
    V##10 = V++;                                \
    V##20 = V++;                                \
    V##01 = V++;                                \
    V##11 = V++;                                \
    V##21 = V++;                                \
    V##02 = V++;                                \
    V##12 = V++;                                \
    V##22 = V++;                                \
  }

/** SET3 : set pointers on a vector3 v (*v0 *v1 *v2)
 * Warning: the pointer v is modified and is ready for a next SET3
 * use *v0 if you need *v
 */
#define SET3(V)                                 \
  double* V##0 MAYBE_UNUSED = V++;                           \
  double* V##1 MAYBE_UNUSED = V++;                           \
  double* V##2 MAYBE_UNUSED = V++;

/** SET3MAYBE : set pointers on a vector3 v (*v0 *v1 *v2) only if v is
 * non null.
 * Warning: the pointer v is modified and is ready for a next SET3
 * use *v0 if you need *v
 */
#define SET3MAYBE(V)                                 \
  double* V##0 MAYBE_UNUSED = 0;                                  \
  double* V##1 MAYBE_UNUSED = 0;                                  \
  double* V##2 MAYBE_UNUSED = 0;                                  \
  if (V)                                             \
  {                                                  \
    V##0 = V++;                                      \
    V##1 = V++;                                      \
    V##2 = V++;                                      \
  }

/** copy a 3x3 matrix or a vector[9]
 * \param[in] a  a[9]
 * \param[out] b  b[9]*/
static inline void cpy3x3(double* restrict a, double* restrict b)
{
  OP3X3(*b++ = *a++);
}

/** add a 3x3 matrix or a vector[9]
 *\param[in] a a[9]
 *\param[in,out]  b b[9]*/
static inline void add3x3(double a[9], double b[9])
{
  OP3X3(*b++ += *a++);
}

/** sub a 3x3 matrix or a vector[9]
 *\param[in] a a[9]
 *\param[in,out] b b[9]*/
static inline void sub3x3(double a[9], double b[9])
{
  OP3X3(*b++ -= *a++);
}

/** copy a vector[3]
 * \param[in] a a[3]
 * \param[out] b b[3]
 */
static inline void cpy3(double a[3], double b[3])
{
  OP3(*b++ = *a++);
}

/** add a vector[3]
 * \param[in] a a[3]
 * \param[in,out] b b[3]*/
static inline void add3(double a[3], double b[3])
{
  OP3(*b++ += *a++);
}

/** sub a vector[3]
 * \param[in] a a[3]
 * \param[in,out] b  b[3]
 */
static inline void sub3(double a[3], double b[3])
{
  OP3(*b++ -= *a++);
}

/** scalar multiplication of a matrix3x3
 * \param[in] scal double scalar
 * \param[in,out] m  m[9]
 */
static inline void scal3x3(double scal, double m[9])
{
  OP3X3(*m++ *= scal);
}

/** diagonal scaling of a vector
 * \param[in] scal_coeffs diagonal part of a matrix
 * \param[in,out] v a 3D vector
 */
static inline void diag_scal3(double* restrict scal_coeffs, double* restrict v)
{
  OP3(*v++ *= *scal_coeffs++);
}

/** scalar multiplication of a vector3
 * \param[in] scal double scalar
 * \param[in,out] v v[3]
 */
static inline void scal3(double scal, double* v)
{
  OP3(*v++ *= scal);
}


/** copy & transpose a matrix
 * \param[in] a *a
 * \param[out] b transpose(*a)
 */
static inline void cpytr3x3(double* restrict a, double* restrict b)
{
  SET3X3(a);
  SET3X3(b);
  *b00 = *a00;
  *b10 = *a01;
  *b20 = *a02;
  *b01 = *a10;
  *b11 = *a11;
  *b21 = *a12;
  *b02 = *a20;
  *b12 = *a21;
  *b22 = *a22;
}

/** matrix vector multiplication
 * \param[in] a 3 by 3 matrix in col-major
 * \param[in] v 3 dimensional vector
 * \param[out] r \f$r = a*v\f$
 */
static inline void mv3x3(double* restrict a, double* restrict v, double* restrict r)
{

  double* pr;

  pr = r;

  *pr++ = *a++ * *v;
  *pr++ = *a++ * *v;
  *pr++ = *a++ * *v++;

  pr = r;

  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v++;

  pr = r;

  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v++;
}

/** transpose(matrix) vector multiplication
 * \param[in] a 3 by 3 matrix in col-major
 * \param[in] v 3 dimensional vector
 * \param[out] r \f$r = a^T*v\f$
 */
static inline void mtv3x3(double* restrict a, double* restrict v, double* restrict r)
{

  double* pv = v;

  *r = *a++ * *pv++;
  *r += *a++ * *pv++;
  *r++ += *a++ * *pv;

  pv = v;

  *r = *a++ * *pv++;
  *r += *a++ * *pv++;
  *r++ += *a++ * *pv;

  pv = v;

  *r = *a++ * *pv++;
  *r += *a++ * *pv++;
  *r += *a * *pv;
}

/** add a matrix vector multiplication
 * \param[in] a a[4]
 * \param[in] v v[2]
 * \param[out] r r[2] the result of r += av
 */
static inline void mvp2x2(const double* restrict a, const double* restrict v, double* restrict r)
{

  double* pr;

  pr = r;

  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v++;

  pr = r;

  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v++;
}

/** add a matrix vector multiplication
 * \param[in] a a[9]
 * \param[in] v v[3]
 * \param[out] r r[3] the result of r += av
 */
static inline void mvp3x3(const double* restrict a, const double* restrict v, double* restrict r)
{

  double* pr;

  pr = r;

  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v++;

  pr = r;

  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v++;

  pr = r;

  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v;
  *pr++ += *a++ * *v++;
}

static inline void mvp5x5(const double* restrict a, const double* restrict v, double* restrict r)
{

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

/** add a matrix vector multiplication scaled by alpha 
 * \param[in] alpha scalar coeff
 * \param[in] a a[9]
 * \param[in] v v[3]
 * \param[out] r r[3] the result of r += av
 */
static inline void mvp_alpha3x3(double alpha, const double* restrict a, const double* restrict v, double* restrict r)
{

  double* pr;

  pr = r;

  *pr++ += alpha * *a++ * *v;
  *pr++ += alpha * *a++ * *v;
  *pr++ += alpha * *a++ * *v++;

  pr = r;

  *pr++ += alpha *  *a++ * *v;
  *pr++ += alpha *  *a++ * *v;
  *pr++ += alpha *  *a++ * *v++;

  pr = r;

  *pr++ += alpha *  *a++ * *v;
  *pr++ += alpha *  *a++ * *v;
  *pr++ += alpha *  *a++ * *v++;
}
/** minux the result a matrix vector multiplication
 * \param[in] a matrix
 * \param[in] v vector
 * \param[out] r the result of r -= av
 */
static inline void mvm3x3(double* restrict a, double* restrict v, double* restrict r)
{

  double* pr;

  pr = r;

  *pr++ -= *a++ * *v;
  *pr++ -= *a++ * *v;
  *pr++ -= *a++ * *v++;

  pr = r;

  *pr++ -= *a++ * *v;
  *pr++ -= *a++ * *v;
  *pr++ -= *a++ * *v++;

  pr = r;

  *pr++ -= *a++ * *v;
  *pr++ -= *a++ * *v;
  *pr++ -= *a++ * *v++;
}


/** transpose(matrix) vector multiplication
 * \param[in] a 3 by 3 matrix in col-major
 * \param[in] v 3 dimensional vector
 * \param[out] r \f$r = a^T*v\f$
 */
static inline void mtvm3x3(double* restrict a, double* restrict v, double* restrict r)
{

  double* pv = v;

  *r -= *a++ * *pv++;
  *r -= *a++ * *pv++;
  *r++ -= *a++ * *pv;

  pv = v;

  *r -= *a++ * *pv++;
  *r -= *a++ * *pv++;
  *r++ -= *a++ * *pv;

  pv = v;

  *r -= *a++ * *pv++;
  *r -= *a++ * *pv++;
  *r -= *a * *pv;
}
/** matrix matrix multiplication : c = a * b
 * \param[in] a  a[9]
 * \param[in] b b[9]
 * \param[out] c c[9]
 */
static inline void mm3x3(double* restrict a, double* restrict b, double* restrict c)
{

  SET3X3(a);
  SET3X3(b);
  SET3X3(c);

  *c00 = *a00 * *b00 + *a01 * *b10 + *a02 * *b20;
  *c01 = *a00 * *b01 + *a01 * *b11 + *a02 * *b21;
  *c02 = *a00 * *b02 + *a01 * *b12 + *a02 * *b22;

  *c10 = *a10 * *b00 + *a11 * *b10 + *a12 * *b20;
  *c11 = *a10 * *b01 + *a11 * *b11 + *a12 * *b21;
  *c12 = *a10 * *b02 + *a11 * *b12 + *a12 * *b22;

  *c20 = *a20 * *b00 + *a21 * *b10 + *a22 * *b20;
  *c21 = *a20 * *b01 + *a21 * *b11 + *a22 * *b21;
  *c22 = *a20 * *b02 + *a21 * *b12 + *a22 * *b22;

}

/** add a matrix matrix multiplication : c += a*b
 * \param[in] a a[9]
 * \param[in] b b[9]
 * \param[out] c c[9]
 */
static inline void mmp3x3(double* restrict a, double* restrict b, double* restrict c)
{

  SET3X3(a);
  SET3X3(b);
  SET3X3(c);

  *c00 += *a00 * *b00 + *a01 * *b10 + *a02 * *b20;
  *c01 += *a00 * *b01 + *a01 * *b11 + *a02 * *b21;
  *c02 += *a00 * *b02 + *a01 * *b12 + *a02 * *b22;

  *c10 += *a10 * *b00 + *a11 * *b10 + *a12 * *b20;
  *c11 += *a10 * *b01 + *a11 * *b11 + *a12 * *b21;
  *c12 += *a10 * *b02 + *a11 * *b12 + *a12 * *b22;

  *c20 += *a20 * *b00 + *a21 * *b10 + *a22 * *b20;
  *c21 += *a20 * *b01 + *a21 * *b11 + *a22 * *b21;
  *c22 += *a20 * *b02 + *a21 * *b12 + *a22 * *b22;

}

/** sub a matrix matrix multiplication : c -= a*b
 * \param[in] a a[9]
 * \param[in] b b[9]
 * \param[out] c c[9]
 */
static inline void mmm3x3(double* restrict a, double* restrict b, double* restrict c)
{

  SET3X3(a);
  SET3X3(b);
  SET3X3(c);

  *c00 -= *a00 * *b00 + *a01 * *b10 + *a02 * *b20;
  *c01 -= *a00 * *b01 + *a01 * *b11 + *a02 * *b21;
  *c02 -= *a00 * *b02 + *a01 * *b12 + *a02 * *b22;

  *c10 -= *a10 * *b00 + *a11 * *b10 + *a12 * *b20;
  *c11 -= *a10 * *b01 + *a11 * *b11 + *a12 * *b21;
  *c12 -= *a10 * *b02 + *a11 * *b12 + *a12 * *b22;

  *c20 -= *a20 * *b00 + *a21 * *b10 + *a22 * *b20;
  *c21 -= *a20 * *b01 + *a21 * *b11 + *a22 * *b21;
  *c22 -= *a20 * *b02 + *a21 * *b12 + *a22 * *b22;

}

/** determinant
 * \param[in] a double* a
 * \return the value of the determinant
 */
static inline double det3x3(double* a)
{
  SET3X3(a);

  return
    *a00 * *a11 * *a22 + *a01 * *a12 * *a20 + *a02 * *a10 * *a21 -
    *a00 * *a12 * *a21 - *a01 * *a10 * *a22 - *a02 * *a11 * *a20;
}


/** system resolution : x <- sol(Ax = b)
 * \param[in] a double a[9]
 * \param[out] x double x[3]
 * \param[in] b double b[3]
 * \return 0 if success, 1 if failed
 */
WARN_RESULT_IGNORED
static inline int solv3x3(double* restrict a, double* restrict x, double* restrict b)
{

  SET3X3(a);
  double* b0 = b++;
  double* b1 = b++;
  double* b2 = b;

  double det =
    *a00 * *a11 * *a22 + *a01 * *a12 * *a20 + *a02 * *a10 * *a21 -
    *a00 * *a12 * *a21 - *a01 * *a10 * *a22 - *a02 * *a11 * *a20;

  if (fabs(det) > DBL_EPSILON)
  {
    double idet = 1.0 / det;
    *x++ = idet * (*a01 * *a12 * *b2 +  *a02 * *a21 * *b1 +
                   *a11 * *a22 * *b0 -  *a01 * *a22 * *b1 -
                   *a02 * *a11 * *b2 -  *a12 * *a21 * *b0);
    *x++ = idet * (*a00 * *a22 * *b1 +  *a02 * *a10 * *b2 +
                   *a12 * *a20 * *b0 -  *a00 * *a12 * *b2 -
                   *a02 * *a20 * *b1 -  *a10 * *a22 * *b0);
    *x   = idet * (*a00 * *a11 * *b2 +  *a01 * *a20 * *b1 +
                   *a10 * *a21 * *b0 -  *a00 * *a21 * *b1 -
                   *a01 * *a10 * *b2 -  *a11 * *a20 * *b0);
  }
  else
  {
    *x++ = NAN;
    *x++ = NAN;
    *x   = NAN;
    return 1;
  }
  return 0;
}

/** check equality : a[9] == b[9]
 * \param[in] a double a[9]
 * \param[in] b double b[9]
 * \return 1 if true, 0 if false
 */
static inline int equal3x3(double* restrict a, double* restrict b)
{
  return *a++ == *b++ &&
         *a++ == *b++ &&
         *a++ == *b++ &&
         *a++ == *b++ &&
         *a++ == *b++ &&
         *a++ == *b++ &&
         *a++ == *b++ &&
         *a++ == *b++ &&
         *a == *b;
}

/** check equality : a[3] == b[3]
 * \param[in] a double a[3]
 * \param[in] b double b[3]
 * \return 1 if true, 0 if false
 */
static inline int equal3(double* restrict a, double* restrict b)
{
  return *a++ == *b++ &&
         *a++ == *b++ &&
         *a == *b;
}

/** scalar product : c <- a.b
 * \param[in] a double a[3]
 * \param[in] b double b[3]
 * \return the scalar product
 */
static inline double dot3(double* restrict a, double* restrict b)
{
  double r;
  r = *a++ * * b++;
  r += *a++ * * b++;
  r += *a * *b;
  return r;
}

/** cross product : c <- a x b
 * \param[in] a double a[3]
 * \param[in] b double b[3]
 * \param[out] c double c[3]
 */
static inline void cross3(double* restrict a, double* restrict b, double* restrict c)
{
  double* a0 = a++;
  double* a1 = a++;
  double* a2 = a;
  double* b0 = b++;
  double* b1 = b++;
  double* b2 = b;

  *c++ = *a1 * *b2 - *a2 * *b1;
  *c++ = *a2 * *b0 - *a0 * *b2;
  *c   = *a0 * *b1 - *a1 * *b0;
}

/** norm : || a ||
 *  may underflow & overflow
 * \param[in] a a[2]
 * \return the norm
 */
static inline double hypot2(double* a)
{
  double r;

  r = *a * *a;
  a++;
  r += *a * *a;
  return sqrt(r);
}


/** norm : || a ||
 *  may underflow & overflow
 * \param[in] a a[3]
 * \return the norm
 */
static inline double hypot3(double* a)
{
  double r;

  r = *a * *a;
  a++;
  r += *a * *a;
  a++;
  r += *a * *a;
  return sqrt(r);
}
static inline double hypot5(double* a)
{
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
static inline double hypot9(double* a)
{
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


/** extract3x3 : copy a sub 3x3 matrix of *a into *b */
/* \param[in] n row numbers of matrix a
 * \param[in] i0 row of first element
 * \param[in] j0 column of first element
 * \param[in] a a[n,n] matrix
 * \param[out] b b[9] a 3x3 matrix */
static inline void extract3x3(int n, int i0, int j0, double* restrict a, double* restrict b)
{

  int k0 = i0 + n * j0;

  int nm3 = n - 3;

  a += k0;

  *b++ = *a++;
  *b++ = *a++;
  *b++ = *a++;
  a += nm3;
  *b++ = *a++;
  *b++ = *a++;
  *b++ = *a++;
  a += nm3;
  *b++ = *a++;
  *b++ = *a++;
  *b   = *a;
}

/** insert3x3 : insert a 3x3 matrix *b into *a */
/* \param[in] n row numbers of matrix a
 * \param[in] i0 row of first element
 * \param[in] j0 column of first element
 * \param[in,out] a  a[n,n] matrix
 * \param[in] b b[9] a 3x3 matrix */
static inline void insert3x3(int n, int i0, int j0, double* restrict a, double* restrict b)
{

  int k0 = i0 + n * j0;

  int nm3 = n - 3;

  a += k0;

  *a++ = *b++;
  *a++ = *b++;
  *a++ = *b++;
  a += nm3;
  *a++ = *b++;
  *a++ = *b++;
  *a++ = *b++;
  a += nm3;
  *a++ = *b++;
  *a++ = *b++;
  *a = *b;
}

/** print a matrix
 * \param mat the matrix
 */
void print3x3(double* mat);


/** print a vector
 * \param[in] v the vector
 */
void print3(double* v);

/** orthoBaseFromVector : From a vector A, build a matrix (A,A1,A2) such that it is an
 * orthonormal.
 * \param[in,out] Ax first component of the vector A
 * \param[in,out] Ay second component of the vector A
 * \param[in,out] Az third component of the vector A
 * \param[out] A1x first component of the vector A
 * \param[out] A1y second component of the vector A
 * \param[out] A1z third component of the vector A
 * \param[out] A2x first component of the vector A
 * \param[out] A2y second component of the vector A
 * \param[out] A2z third component of the vector A
 *
 * \return 0 if success, 1 if there is a problem
*/
WARN_RESULT_IGNORED
static inline int orthoBaseFromVector(double *Ax, double *Ay, double *Az,
                                      double *A1x, double *A1y, double *A1z,
                                      double *A2x, double *A2y, double *A2z)
{

  double normA = sqrt((*Ax) * (*Ax) + (*Ay) * (*Ay) + (*Az) * (*Az));
  if (normA == 0.)
  {
    (*Ax) = NAN;
    (*Ay) = NAN;
    (*Az) = NAN;
    (*A1x) = NAN;
    (*A1y) = NAN;
    (*A1z) = NAN;
    (*A2x) = NAN;
    (*A2y) = NAN;
    (*A2z) = NAN;
    return 1;
  }
  (*Ax) /= normA;
  (*Ay) /= normA;
  (*Az) /= normA;
  /*build (*A1*/
  if (fabs(*Ax) > fabs(*Ay))
  {
    if (fabs((*Ax)) > fabs((*Az))) /*(*Ax is the bigest*/
    {
      (*A1x) = (*Ay);
      (*A1y) = -(*Ax);
      (*A1z) = 0;
    }
    else /*(*Az is the biggest*/
    {
      (*A1z) = (*Ax);
      (*A1y) = -(*Az);
      (*A1x) = 0;
    }
  }
  else if (fabs(*Ay) > fabs(*Az)) /*(*Ay is the bigest*/
  {
    (*A1y) = (*Ax);
    (*A1x) = -(*Ay);
    (*A1z) = 0;
  }
  else /*(*Az is the biggest*/
  {
    (*A1z) = (*Ax);
    (*A1y) = -(*Az);
    (*A1x) = 0;
  }
  double normA1 = sqrt((*A1x) * (*A1x) + (*A1y) * (*A1y) + (*A1z) * (*A1z));
  (*A1x) /= normA1;
  (*A1y) /= normA1;
  (*A1z) /= normA1;

  (*A2x) = *Ay * *A1z - *Az * *A1y;
  (*A2y) = *Az * *A1x - *Ax * *A1z;
  (*A2z) = *Ax * *A1y - *Ay * *A1x;
  return 0;
}


/** solve Ax = b by partial pivoting Gaussian elimination. This function is 10
 * to 20 times faster than calling LAPACK (tested with netlib).
 *
 * \param a column-major matrix (not modified)
 * \param[in,out] b on input, the right-hand side; on output the solution x
 *
 * \return 0 if ok, otherwise the column where no pivot could be selected
 */
WARN_RESULT_IGNORED
static inline int solve_3x3_gepp(const double* restrict a, double* restrict b)
{
  double lp0, lp1, lp2, lm1, lm2, ln1, ln2;
  double bl, bm, bn;
  double factor1, factor2;
  double sol0, sol1, sol2;
  double alp0;
  double a0 = a[0];
  double a1 = a[1];
  double a2 = a[2];
  double aa0 = fabs(a0);
  double aa1 = fabs(a1);
  double aa2 = fabs(a2);
  /* do partial search of pivot */
  int pivot2, pivot1 = aa0 >= aa1 ? aa0 >= aa2 ? 0 : 20 : aa1 >= aa2 ? 10 : 20;
  int info = 0;

  /* We are going to put the system in the form
   * | lp0 lp1 lp2 ; bl |
   * |  0  lm1 lm2 ; bm |
   * |  0  ln1 ln2 ; bn |
   */
  switch (pivot1)
  {
    case 0: /* first element is pivot, first line does not change */
      factor1 = a1/a0;
      factor2 = a2/a0;
      lp0 = a0;
      alp0 = fabs(a0);
      lp1 = a[3];
      lp2 = a[6];
      lm1 = a[4] - factor1*lp1;
      lm2 = a[7] - factor1*lp2;
      ln1 = a[5] - factor2*lp1;
      ln2 = a[8] - factor2*lp2;
      bl = b[0];
      bm = b[1] - factor1*bl;
      bn = b[2] - factor2*bl;
      break;
    case 10: /* first element is pivot, first line does not change */
      factor1 = a0/a1;
      factor2 = a2/a1;
      lp0 = a1;
      alp0 = fabs(a1);
      lp1 = a[4];
      lp2 = a[7];
      lm1 = a[3] - factor1*lp1;
      lm2 = a[6] - factor1*lp2;
      ln1 = a[5] - factor2*lp1;
      ln2 = a[8] - factor2*lp2;
      bl = b[1];
      bm = b[0] - factor1*bl;
      bn = b[2] - factor2*bl;
      break;
    case 20: /* first element is pivot, first line does not change */
      factor1 = a0/a2;
      factor2 = a1/a2;
      lp0 = a2;
      alp0 = fabs(a2);
      lp1 = a[5];
      lp2 = a[8];
      lm1 = a[3] - factor1*lp1;
      lm2 = a[6] - factor1*lp2;
      ln1 = a[4] - factor2*lp1;
      ln2 = a[7] - factor2*lp2;
      bl = b[2];
      bm = b[0] - factor1*bl;
      bn = b[1] - factor2*bl;
      break;
    default:
      exit(EXIT_FAILURE);
  }

  if (alp0 <= DBL_EPSILON)
  {
    info = 1;
    return info;
  }

  double alm1 = fabs(lm1);
  double aln1 = fabs(ln1);
  pivot2 = alm1 >= aln1 ? 0 : 1;

  /* We now solve the system
   * | lm1 lm2 ; bm |
   * | ln1 ln2 ; bn |
   */
  switch (pivot2)
  {
    case 0:
      if (alm1 < DBL_EPSILON)
      {
        info = 1;
        return info;
      }
      factor1 = ln1/lm1;
      sol2 = (bn - factor1*bm)/(ln2 - factor1*lm2);
      sol1 = (bm - lm2*sol2)/lm1;
      break;
    case 1:
      if (aln1 < DBL_EPSILON)
      {
        info = 1;
        return info;
      }
      factor1 = lm1/ln1;
      sol2 = (bm - factor1*bn)/(lm2 - factor1*ln2);
      sol1 = (bn - ln2*sol2)/ln1;
      break;
    default:
      exit(EXIT_FAILURE);
  }
  sol0 = (bl - sol1*lp1 - sol2*lp2)/lp0;

  b[0] = sol0;
  b[1] = sol1;
  b[2] = sol2;

  return info;
}


#define mat_elem(a, y, x, n) (a + ((y) * (n) + (x)))

static void swap_row(double *a, double *b, int r1, int r2, int n)
{
	double tmp, *p1, *p2;
	int i;

	if (r1 == r2) return;
	for (i = 0; i < n; i++) {
		p1 = mat_elem(a, r1, i, n);
		p2 = mat_elem(a, r2, i, n);
		tmp = *p1, *p1 = *p2, *p2 = tmp;
	}
	tmp = b[r1], b[r1] = b[r2], b[r2] = tmp;
}

static inline  void solve_nxn_gepp(int n, double *a, double *b, double *x)
{
#define A(y, x) (*mat_elem(a, y, x, n))
	int  j, col, row, max_row,dia;
	double max, tmp;

	for (dia = 0; dia < n; dia++) {
		max_row = dia, max = A(dia, dia);

		for (row = dia + 1; row < n; row++)
			if ((tmp = fabs(A(row, dia))) > max)
				max_row = row, max = tmp;

		swap_row(a, b, dia, max_row, n);

		for (row = dia + 1; row < n; row++) {
			tmp = A(row, dia) / A(dia, dia);
			for (col = dia+1; col < n; col++)
				A(row, col) -= tmp * A(dia, col);
			A(row, dia) = 0;
			b[row] -= tmp * b[dia];
		}
	}
	for (row = n - 1; row >= 0; row--) {
		tmp = b[row];
		for (j = n - 1; j > row; j--)
			tmp -= x[j] * A(row, j);
		x[row] = tmp / A(row, row);
	}
#undef A
}

/** Computation of the eigenvalues of a symmetric 3x3 real matrix
 *
 * \param a symmetric column-major matrix (not modified)
 * \param b a 3x3 work matrix
 * \param[in,out] eig eigenvalie in decreasing order
 *
 * \return 0 all the time
 */
WARN_RESULT_IGNORED
static inline int eig_3x3(double* restrict a, double* restrict b, double* restrict eig)
{
  SET3X3(a);
  SET3X3(b);
  double pi = M_PI;
  double p1 =  *a01* *a01 +  *a02* *a02 +  *a12* *a12;
  if (p1 == 0)
  {
    /* A is diagonal. */
    eig[0] =  *a00;
    eig[1] =  *a11;
    eig[2] =  *a22;
  }
  else
  {
    double q = ( *a00+ *a11+ *a22)/3.0;
    double  p2 = ( *a00 - q)*( *a00 - q) + ( *a11 - q)*( *a11 - q) + ( *a22 - q)*( *a22 - q) + 2 * p1;
    double p = sqrt(p2 / 6.0);
    *b00 = (1 / p) * ( *a00 - q);
    *b11 = (1 / p) * ( *a11 - q);
    *b22 = (1 / p) * ( *a22 - q);
    *b01 = (1 / p) * ( *a01);
    *b02 = (1 / p) * ( *a02);
    *b10 = (1 / p) * ( *a10);
    *b12 = (1 / p) * ( *a12);
    *b20 = (1 / p) * ( *a20);
    *b21 = (1 / p) * ( *a21);
    //B = (1 / p) * ( *A - q * I)       % I is the identity matrix
    double det =
      *b00 * *b11 * *b22 + *b01 * *b12 * *b20 + *b02 * *b10 * *b21 -
      *b00 * *b12 * *b21 - *b01 * *b10 * *b22 - *b02 * *b11 * *b20;
    double r = det/2.0;

    /* % In exact arithmetic for a symmetric matrix  -1 <= r <= 1 */
    /* % but computation error can leave it slightly outside this range. */
    double phi;
    if (r <= -1)
      phi = pi / 3.0;
    else if (r >= 1)
      phi = 0.0;
    else
      phi = acos(r) / 3;

   /* % the eigenvalues satisfy eig3 <= eig2 <= eig1 */
    eig[0] = q + 2 * p * cos(phi);
    eig[2] = q + 2 * p * cos(phi + (2*pi/3));
    eig[1] = 3 * q - eig[0] - eig[2];    /* % since trace(A) = eig1 + eig2 + eig3; */
  }

  return 0;
}



#endif
