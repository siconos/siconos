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

#ifndef _op3x3_h_
#define _op3x3_h_

#include <math.h>
#include <float.h>

/** OP3X3(EXPR) do EXPR 9 times
 * \param a C expression that should contains self incrementing
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
 * \param a C expression that should contains self incrementing
 *        pointers on arrays[9] */
#define OP3(EXPR)                               \
  do                                            \
  {                                             \
    EXPR;                                       \
    EXPR;                                       \
    EXPR;                                       \
  } while(0)                                    \
 

#if defined(OP3X3_C_STORAGE)
#define SET3X3(V)                             \
  double* V##00 = V++;                        \
  double* V##01 = V++;                        \
  double* V##02 = V++;                        \
  double* V##10 = V++;                        \
  double* V##11 = V++;                        \
  double* V##12 = V++;                        \
  double* V##20 = V++;                        \
  double* V##21 = V++;                        \
  double* V##22 = V++;
#define SET3X3MAYBE(V)                          \
  double* V##00 = 0;                            \
  double* V##01 = 0;                            \
  double* V##02 = 0;                            \
  double* V##10 = 0;                            \
  double* V##11 = 0;                            \
  double* V##12 = 0;                            \
  double* V##20 = 0;                            \
  double* V##21 = 0;                            \
  double* V##22 = 0;                            \
  if (V)                                        \
  {                                             \
    V##00 = V++;                                \
    V##01 = V++;                                \
    V##02 = V++;                                \
    V##10 = V++;                                \
    V##11 = V++;                                \
    V##12 = V++;                                \
    V##20 = V++;                                \
    V##21 = V++;                                \
    V##22 = V++;                                \
  }

#else // fortran storage

/** SET3X3 : set pointers on a 3x3 matrix a (*a00 *a01 *a10 etc.)
 * warning the pointer a is modified (use a00 instead) and is ready
 * for a next SET3X3
 */
#define SET3X3(V)                                                       \
  double* V##00 = V++;                                                  \
  double* V##10 = V++;                                                  \
  double* V##20 = V++;                                                  \
  double* V##01 = V++;                                                  \
  double* V##11 = V++;                                                  \
  double* V##21 = V++;                                                  \
  double* V##02 = V++;                                                  \
  double* V##12 = V++;                                                  \
  double* V##22 = V++;
#define SET3X3MAYBE(V)                          \
  double* V##00 = 0;                            \
  double* V##10 = 0;                            \
  double* V##20 = 0;                            \
  double* V##01 = 0;                            \
  double* V##11 = 0;                            \
  double* V##21 = 0;                            \
  double* V##02 = 0;                            \
  double* V##12 = 0;                            \
  double* V##22 = 0;                            \
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
#endif

/** SET3 : set pointers on a vector3 v (*v0 *v1 *v2)
 * Warning: the pointer v is modified and is ready for a next SET3
 * use *v0 if you need *v
 */
#define SET3(V)                                 \
  double* V##0 = V++;                           \
  double* V##1 = V++;                           \
  double* V##2 = V++;

/** SET3MAYBE : set pointers on a vector3 v (*v0 *v1 *v2) only if v is
 * non null.
 * Warning: the pointer v is modified and is ready for a next SET3
 * use *v0 if you need *v
 */
#define SET3MAYBE(V)                                 \
  double* V##0 = 0;                                  \
  double* V##1 = 0;                                  \
  double* V##2 = 0;                                  \
  if (V)                                             \
  {                                                  \
    V##0 = V++;                                      \
    V##1 = V++;                                      \
    V##2 = V++;                                      \
  }

/** copy a 3x3 matrix or a vector[9]
 *\param[in] a[9]
 *\param[out] b[9]*/
static inline void cpy3x3(double a[9], double b[9])
{
  OP3X3(*b++ = *a++);
}

/** add a 3x3 matrix or a vector[9]
 *\param[in] a[9]
 *\param[in,out] b[9]*/
static inline void add3x3(double a[9], double b[9])
{
  OP3X3(*b++ += *a++);
}

/** sub a 3x3 matrix or a vector[9]
 *\param[in] a[9]
 *\param[in,out] b[9]*/
static inline void sub3x3(double a[9], double b[9])
{
  OP3X3(*b++ -= *a++);
}

/** copy a vector[3]
 *\param[in] a[3]
 *\param[out] b[3]*/
static inline void cpy3(double a[3], double b[3])
{
  OP3(*b++ = *a++);
}

/** add a vector[3]
 *\param[in] a[3]
 *\param[in,out] b[3]*/
static inline void add3(double a[3], double b[3])
{
  OP3(*b++ += *a++);
}

/** sub a vector[3]
 *\param[in] a[3]
 *\param[in,out] b[3]*/
static inline void sub3(double a[3], double b[3])
{
  OP3(*b++ -= *a++);
}

/** scalar multiplication of a matrix3x3
 * \param[in] double scalar
 * \param[in,out] b[9]
 */
static inline void scal3x3(double scal, double m[9])
{
  OP3X3(*m++ *= scal);
}

/** scalar multiplication of a vector3
 * \param[in] double scalar
 * \param[in,out] v[3]
 */
static inline void scal3(double scal, double v[3])
{
  OP3(*v++ *= scal);
}


/** copy & transpose a matrix
 * \param[in] *a
 * \param[out] transpose(*a)
 */
static inline void cpytr3x3(double a[9], double b[9])
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
};

/** matrix vector multiplication
 * \param[in] a[9]
 * \param[in] v[3]
 * \param[out] r[3]
 */
static inline void mv3x3(double a[9], double v[9], double r[9])
{

#if defined(OP3X3_C_STORAGE)
  double* pv;

  pv = v;
  *r++ = *a++ * *pv++ + *a++ * *pv++ + *a++ * *pv++;

  pv = v;
  *r++ = *a++ * *pv++ + *a++ * *pv++ + *a++ * *pv++;

  pv = v;
  *r++ = *a++ * *pv++ + *a++ * *pv++ + *a++ * *pv++;
#else
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

#endif

}


/** add a matrix vector multiplication
 * \param[in] a[9]
 * \param[in] v[3]
 * \param[out] r[3]
 */
static inline void mvp3x3(double a[9], double v[9], double r[9])
{

#if defined(OP3X3_C_STORAGE)
  double* pv;

  pv = v;
  *r++ += *a++ * *pv++ + *a++ * *pv++ + *a++ * *pv++;

  pv = v;
  *r++ += *a++ * *pv++ + *a++ * *pv++ + *a++ * *pv++;

  pv = v;
  *r++ += *a++ * *pv++ + *a++ * *pv++ + *a++ * *pv++;
#else
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

#endif

}

/** matrix matrix multiplication : c = a * b
 * \param[in] a[9]
 * \param[in] b[9]
 * \param[out] c[9]
 */
static inline void mm3x3(double a[9], double b[9], double c[9])
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
 * \param[in] a[9]
 * \param[in] b[9]
 * \param[out] c[9]
 */
static inline void mmp3x3(double a[9], double b[9], double c[9])
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
 * \param[in] a[9]
 * \param[in] b[9]
 * \param[out] c[9]
 */
static inline void mmm3x3(double a[9], double b[9], double c[9])
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
 * \param[in] double* a
 */
static inline double det3x3(double a[9])
{
  SET3X3(a);

  return
    *a00 * *a11 * *a22 + *a01 * *a12 * *a20 + *a02 * *a10 * *a21 -
    *a00 * *a12 * *a21 - *a01 * *a10 * *a22 - *a02 * *a11 * *a20;
}


/** system resolution : x <- sol(Ax = b)
 * \param[in] double a[9]
 * \param[out] double x[3]
 * \param[in] double b[3]
 */
static inline void solv3x3(double a[9], double x[3], double b[3])
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
  }
}

/** check equality : a[9] == b[9]
 * \param[in] double a[9]
 * \param[in] double b[9]
 */
static inline int equal3x3(double a[9], double b[9])
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
 * \param[in] double a[3]
 * \param[in] double b[3]
 */
static inline int equal3(double a[3], double b[3])
{
  return *a++ == *b++ &&
         *a++ == *b++ &&
         *a == *b;
}

/** scalar product : c <- a.b
 * \param[in] double a[3]
 * \param[in] double b[3]
 * \return double
 */
static inline double dot3(double a[3], double b[3])
{
  double r;
  r = *a++ * * b++;
  r += *a++ * * b++;
  r += *a * *b;
  return r;
};

/** cross product : c <- a x b
 * \param[in] double a[3]
 * \param[in] double b[3]
 * \param[out] double c[3]
 */
static inline void cross3(double* a, double* b, double* c)
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
 * \param[in] a[3]
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


/* check nan of component
 * \param double* a
 */
static inline int isnan3(double* a)
{
  double* a0 = a++;
  double* a1 = a++;
  double* a2 = a;

  return isnan(*a0) || isnan(*a1) || isnan(*a2);
}

/** extract3x3 : copy a sub 3x3 matrix of *a into *b */
/* \param[in] n row numbers of matrix a
 * \param[in] i0 row of first element
 * \param[in] j0 column of first element
 * \param[in] a[n,n] matrix
 * \param[out] b[9] a 3x3 matrix */
static inline void extract3x3(int n, int i0, int j0,
                              double* a, double* b)
{
#if defined(OP3X3_C_STORAGE)
  int k0 = n * i0 + j0;
#else
  int k0 = i0 + n * j0;
#endif
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
 * \param[in,out] a[n,n] matrix
 * \param[in] b[9] a 3x3 matrix */
static inline void insert3x3(int n, int i0, int j0,
                             double* a, double* b)
{

#if defined(OP3X3_C_STORAGE)
  int k0 = n * i0 + j0;
#else
  int k0 = i0 + n * j0;
#endif
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
 * \param double* a
 */
void print3x3(double* mat);


/** print a vector
 * \param[in] double* v
 */
void print3(double* v);
/** orthoBaseFromVector : From A, build (A,A1,A2) such that it was an orthonormal base.
 * \param[in,out] Ax
 * \param[in,out] Ay
 * \param[in,out] Az
 * \param[out] A1x
 * \param[out] A1y
 * \param[out] A1z
 * \param[out] A2x
 * \param[out] A2y
 * \param[out] A2z
*/
static inline void orthoBaseFromVector(double *Ax, double *Ay, double *Az, double *A1x, double *A1y, double *A1z, double *A2x, double *A2y, double *A2z)
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
    return;
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
    else  /*(*Az is the biggest*/
    {
      (*A1z) = (*Ax);
      (*A1y) = -(*Az);
      (*A1x) = 0;
    }
  }
  else if (fabs(*Ay) > fabs(*Az))  /*(*Ay is the bigest*/
  {
    (*A1y) = (*Ax);
    (*A1x) = -(*Ay);
    (*A1z) = 0;
  }
  else  /*(*Az is the biggest*/
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
}


#endif
