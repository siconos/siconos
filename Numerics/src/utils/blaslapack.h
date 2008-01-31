/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */

#ifndef BLASLAPACK_H
#define BLASLAPACK_H

#if __GNUC__ >= 4
#include "f2c.h"
#else
#include <g2c.h>
#endif

//#define _F2C_INCLUDE_H // to avoid f2c.h include => conflict with g2c.h
#if  defined(RIOS) && !defined(CLAPACK)
#define F77NAME(x) x
#else
#define F77NAME(x) x##_
#endif



#ifndef _BLAS1_H_

#define _BLAS1_H_

extern  double F77NAME(dasum)(const integer *n, const double *dx, const integer *incx);


extern void F77NAME(daxpy)(const integer *n, const double *da, const double *dx,
                           const integer *incx, double *dy, const integer *incy);

extern void F77NAME(dcopy)(const integer *n, const double *dx, const integer *incx, double *dy,
                           const integer *incy);


extern double F77NAME(ddot)(const integer *n, const double *dx, const integer *incx,
                            const double *dy, const integer *incy);

extern double F77NAME(dnrm2)(const integer *n, const double *dx, const integer *incx);

extern void F77NAME(drot)(const integer *n, double *dx, const integer *incx, double *dy,
                          const integer *incy, const double *c, const double *s);

extern void F77NAME(drotg)(double *da, double *db, double *c, double *s);

extern void F77NAME(dscal)(const integer *n, double *da, double *dx, const integer *incx);

extern void F77NAME(dswap)(const integer *n, double *dx, const integer *incx, double *dy,
                           const integer *incy);

extern integer F77NAME(idamax)(const integer *n, const double *dx, const integer *incx);

#ifdef LA_COMPLEX_SUPPORT

extern double F77NAME(zdotc)(doublecomplex *c, const integer *n,
                             const doublecomplex *cx,
                             const integer *incx, const doublecomplex *cy, const integer *incy);

extern double F77NAME(zdotu)(doublecomplex *c, const integer *n,
                             const doublecomplex *cx, const integer *incx,
                             const doublecomplex *cy, const integer *incy);

extern void F77NAME(zaxpy)(const integer *n, const doublecomplex *da,
                           const doublecomplex *dx,
                           const integer *incx, doublecomplex *dy,
                           const integer *incy);

extern void F77NAME(zcopy)(const integer *n, const doublecomplex *dx, const integer *incx,
                           doublecomplex *dy, const integer *incy);

extern double  F77NAME(dzasum)(const integer *n, const doublecomplex *dx, const integer *incx);

extern double  F77NAME(dznrm2)(const integer *n, const doublecomplex *dx, const integer *incx);

extern void F77NAME(zdscal)(const integer *n, const double *da, doublecomplex *dx,
                            const integer *incx);

extern void F77NAME(zscal)(const integer *n, const doublecomplex *da, doublecomplex *dx,
                           const integer *incx);

extern integer F77NAME(izamax)(const integer *n, const doublecomplex *dx, const integer *incx);

extern void F77NAME(zswap)(const integer *n, doublecomplex *dx, const integer *incx,
                           doublecomplex *dy, integer *incy);

#endif // LA_COMPLEX_SUPPORT
#endif // _BLAS1_H_


#ifndef _BLAS2_H_
#define _BLAS2_H_





extern void F77NAME(dgemv)(char* trans, integer* M, integer* N, double* alpha,
                           const double* A, integer* lda, const double* dx,
                           integer* incx, double* beta, double* dy, integer* incy);

extern void F77NAME(dgbmv)(char* trans, integer* M, integer* N, integer* kl,
                           integer* ku, double* alpha, const double* A, integer* lda,
                           const double* dx, integer* incx, double* beta,
                           double* dy, integer* incy);

extern void F77NAME(dsymv)(char* uplo, integer* N, double* alpha, const double* A,
                           integer* lda, const double* dx, integer* incx, double* beta,
                           double* dy, integer* incy);

extern void F77NAME(dsbmv)(char* uplo, integer* N, integer* k, double* alpha,
                           const double* A, integer* lda, const double* dx,
                           integer* incx, double* beta, double* dy, integer* incy);

extern void F77NAME(dspmv)(char* uplo, integer* N, double* alpha, double* AP,
                           const double* dx, integer* incx, double* beta, double* dy,
                           integer* incy);

extern void F77NAME(dtrmv)(char* uplo, char* trans, char* diag, const integer* N,
                           const double* A, integer* lda, const double* dx,
                           integer* incx);

// currently not implemented.
//F77NAME(dtbmv) ( UPLO, TRANS, DIAG, N, K, A, LDA, dx, INCX )

extern  void F77NAME(dtrsv)(char* uplo, char* trans, char* diag, const integer* N,
                            double* A, integer* lda, double* dx, integer* incx);

// currently not implemented.
//F77NAME(dtbsv) ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )

// currently not implemented.
//F77NAME(dtpsv) ( UPLO, TRANS, DIAG, N, AP, X, INCX )

extern  void F77NAME(dger)(integer* M, integer* N, const double* alpha,
                           const double* dx, integer* incx,
                           const double* dy, integer* incy,
                           double* A, integer* lda);

extern void F77NAME(dsyr)(char* uplo, integer* N, const double* alpha,
                          const double* dx, integer* incx,
                          double* A, integer* lda);

extern  void F77NAME(dspr)(char* uplo, integer* N, const double* alpha,
                           const double* dx, integer* incx,
                           double* AP);

extern  void F77NAME(dsyr2)(char* uplo, integer* N, const double* alpha,
                            const double* dx, integer* incx,
                            const double* dy, integer* incy, double* A,
                            integer* lda);

extern  void F77NAME(dspr2)(char* uplo, integer* N, const double* alpha,
                            const double* dx, integer* incx,
                            const double* dy, integer* incy,
                            double* AP);


#ifdef LA_COMPLEX_SUPPORT

extern  void F77NAME(zgemv)(const char* trans,
                            const integer* M, const integer* N,
                            const doublecomplex* alpha,
                            const doublecomplex* A, const integer* lda,
                            const doublecomplex* dx, const integer* incx,
                            const doublecomplex* beta,
                            doublecomplex* dy, const integer* incy);

extern  void F77NAME(zgerc)(const integer* M, const integer* N,
                            const doublecomplex* alpha,
                            const doublecomplex* dx, const integer* incx,
                            const doublecomplex* dy, const integer* incy,
                            doublecomplex* A, const integer* lda);
#endif // LA_COMPLEX_SUPPORT

#endif // _BLAS2_H_

#ifndef _BLAS3_H
#define _BLAS3_H


extern void F77NAME(dgemm)(char *transa, char *transb, integer *m, integer *n, integer *k,
                           double *alpha, const double *a, integer *lda, const double *b,
                           integer *ldb, double *beta, double *c, integer *ldc);

extern void F77NAME(dtrsm)(char *side, char *uplo, char *transa, char *diag,
                           integer *m, integer *n, double *alpha, const double *A, integer *lda,
                           const double *B, integer *ldb);

extern void F77NAME(dtrmm)(char *side, char *uplo, char *transa, char *diag,
                           integer *m, integer *n, double *alpha, const double *A, integer *lda,
                           const double *B, integer *ldb);

extern void F77NAME(dsymm)(char *side, char *uplo, integer *m, integer *n,
                           double *alpha, const double *A, integer *lda, const double *B,
                           integer *ldb, double *beta, double *C, integer *ldc);

extern void F77NAME(dsyrk)(char *uplo, char *transa, integer *n, integer *k,
                           double *alpha, double *A, integer *lda, double *beta, double *C,
                           integer *ldc);

extern void F77NAME(dsyr2k)(char *uplo, char *transa, integer *n, integer *k,
                            double *alpha, double *A, integer *lda, double *B, integer *ldb,
                            double *beta, double *C, integer *ldc);

#ifdef LA_COMPLEX_SUPPORT
extern void F77NAME(zgemm)(char *transa, char *transb, integer *m, integer *n, integer *k,
                           doublecomplex *alpha, const doublecomplex *a, integer *lda, const doublecomplex *b,
                           integer *ldb, doublecomplex *beta, doublecomplex *c, integer *ldc);
#endif // LA_COMPLEX_SUPPORT



#endif // _BLAS3_H

#ifndef _LAPACKC_H_
#define _LAPACKC_H_



// *************************** Utility Routines **********************

extern float F77NAME(slamch)(char *t);

extern doublecomplex F77NAME(zlamch)(char *t);


//void F77NAME(zswap)(integer *n, doublecomplex *x, integer *incx, doublecomplex *y, integer *incy);

extern void F77NAME(zgesv)(integer *n, integer *k, doublecomplex *A, integer *lda, integer *ipiv,
                           doublecomplex *X, integer *ldx, integer *info);

extern void F77NAME(zposv)(char *uplo, integer *m, integer *k , doublecomplex *A, integer *lda,
                           doublecomplex *X, integer *ldx, integer *info);

extern void F77NAME(zgels)(char *trans, integer *m, integer *n, integer *nrhs, doublecomplex *A,
                           integer *lda, doublecomplex *B, integer *ldb, doublecomplex *work, integer *lwork, integer *info);

extern void F77NAME(ztimmg)(integer *iflag, integer *m, integer *n, doublecomplex *A, integer *lda,
                            integer *kl, integer *ku);

extern void F77NAME(zlaswp)(integer *n, doublecomplex *A, integer *lda, integer *k1, integer *k2,
                            integer *ipiv, integer *incx);

extern doublecomplex F77NAME(zopla)(char *subname, integer *m, integer *n, integer *kl, integer *ku,
                                    integer *nb);

// ******************* LU Factorization Routines **********************

extern void F77NAME(zgetrf)(integer *m, integer *n, doublecomplex *A, integer *lda, integer *ipiv,
                            integer *info);

extern void F77NAME(zgetf2)(integer *m, integer *n, doublecomplex *A, integer *lda, integer *ipiv,
                            integer *info);

extern void F77NAME(zgbtrf)(integer *m, integer *n, integer *KL, integer *KU, doublecomplex *BM,
                            integer *LDBM, integer *ipiv, integer *info);

extern void F77NAME(zgttrf)(integer *N, doublecomplex *DL, doublecomplex *D, doublecomplex *DU,
                            doublecomplex *DU2, integer *ipiv, integer *info);

extern void F77NAME(zpotrf)(char *UPLO, integer *N, doublecomplex *SM, integer *LDSM,
                            integer *info);

extern void F77NAME(zsytrf)(char *UPLO, integer *N, doublecomplex *SM, integer *LDSM,
                            integer *ipiv, doublecomplex *WORK, integer *LWORK, integer *info);

extern void F77NAME(zpbtrf)(char *UPLO, integer *N, integer *KD, doublecomplex *SBM,
                            integer *LDSM, integer *info);

extern void F77NAME(zpttrf)(integer *N, doublecomplex *D, doublecomplex *E, integer *info);

// ********************* LU Solve Routines ***************************

extern void F77NAME(zgetrs)(char *trans, integer *N, integer *nrhs, doublecomplex *A, integer *lda,
                            integer * ipiv, doublecomplex *b, integer *ldb, integer *info);

extern void F77NAME(zgbtrs)(char *trans, integer *N, integer *kl, integer *ku, integer *nrhs,
                            doublecomplex *AB, integer *ldab, integer *ipiv, doublecomplex *b, integer *ldb, integer *info);

extern void F77NAME(zsytrs)(char *uplo, integer *N, integer *nrhs, doublecomplex *A, integer *lda,
                            integer *ipiv, doublecomplex *b, integer *ldb, integer *info);

extern void F77NAME(zgttrs)(char *trans, integer *N, integer *nrhs, doublecomplex *DL,
                            doublecomplex *D, doublecomplex *DU, doublecomplex *DU2, integer *ipiv, doublecomplex *b,
                            integer *ldb, integer *info);

extern void F77NAME(zpotrs)(char *UPLO, integer *N, integer *nrhs, doublecomplex *A, integer *LDA,
                            doublecomplex *b, integer *ldb, integer *info);

extern void F77NAME(zpttrs)(integer *N, integer *nrhs, doublecomplex *D, doublecomplex *E,
                            doublecomplex *b, integer *ldb, integer *info);

extern void F77NAME(zpbtrs)(char *UPLO, integer *N, integer *KD, integer *nrhs, doublecomplex *AB,
                            integer *LDAB, doublecomplex *b, integer *ldb, integer *info);

// ******************** QR factorizations

extern void F77NAME(zgeqrf)(integer *m, integer *n, doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *work, integer *lwork, integer *info);
extern void F77NAME(zungqr)(integer *m, integer *n, integer *k, doublecomplex *a, integer *lda, const doublecomplex *tau, doublecomplex *work, integer *lwork, integer *info);
extern void F77NAME(zunmqr)(char *side, char *trans, integer *m, integer *n, integer *k, const doublecomplex *a, integer *lda, const doublecomplex *tau, doublecomplex *c, integer *ldc, doublecomplex *work, integer *lwork, integer *info);

// ********************* Eigenvalue/Singular Value Decomposition Drivers

extern void F77NAME(zheevd)(char *jobz, char *uplo, integer *n, doublecomplex *a, integer *lda, double *w, doublecomplex *work, integer *lwork, double *rwork, integer *lrwork, integer *info);
extern void F77NAME(zheevr)(char *jobz, char *range, char *uplo, integer *n, doublecomplex *a, integer *lda, double *vl, double *vu, integer *il, integer *iu, double *abstol, integer *m, double *w, doublecomplex *z, integer *ldz, integer *isuppz, integer *info);

extern void F77NAME(zgesvd)(char *jobu, char *jobvt, integer *m, integer *n, doublecomplex *a, integer *lda, double *sing, doublecomplex *u, integer *ldu, doublecomplex *vt, integer *ldvt, integer *info);
extern void F77NAME(zgesdd)(char *jobz, integer *m, integer *n, doublecomplex *a, integer *lda, double *s, doublecomplex *u, integer *ldu, doublecomplex *vt, integer *ldvt, doublecomplex *work, integer *lwork, double *rwork, integer *iwork, integer *info);


// *******************************


#endif

//      Double precision Lapack routines

#ifndef _DLAPACK_H_
#define _DLAPACK_H_



// *************************** Utility Routines **********************


extern double F77NAME(dlamch)(char *t);



//******************** Linear Equation Solvers *************************
extern void F77NAME(dgesv)(integer *n, integer *k, doublereal *A, integer *lda, integer *ipiv,
                           doublereal *X, integer *ldx, integer *info);

extern void F77NAME(dposv)(char *uplo, integer *m, integer *k , doublereal *A, integer *lda,
                           doublereal *X, integer *ldx, integer *info);

extern  void F77NAME(dsysv)(const char *uplo, integer *n, integer *nrhs,
                            doublereal *A, integer *lda,
                            integer *ipiv,
                            doublereal *B, integer *ldb,
                            doublereal *work, integer *lwork, integer *info);

extern void F77NAME(dtrtrs)(char *uplo, char *trans, char *diag, integer *n,
                            integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *info);

//******************** Lapack Utility Routines ************************

extern void F77NAME(dgels)(char *trans, integer *m, integer *n, integer *nrhs, doublereal *A,
                           integer *lda, doublereal *B, integer *ldb, doublereal *work, integer *lwork, integer *info);

extern void F77NAME(dtimmg)(integer *iflag, integer *m, integer *n, doublereal *A, integer *lda,
                            integer *kl, integer *ku);

extern void F77NAME(dlaswp)(integer *n, doublereal *A, integer *lda, integer *k1, integer *k2,
                            integer *ipiv, integer *incx);

extern  doublereal F77NAME(dopla)(char *subname, integer *m, integer *n, integer *kl, integer *ku,
                                  integer *nb);

// ******************* LU Factorization Routines **********************

extern void F77NAME(dgetrf)(integer *m, integer *n, doublereal *A, integer *lda, integer *ipiv,
                            integer *info);

extern void F77NAME(dgetri)(integer *n, doublereal* A, integer *lda, integer* ipiv, doublereal* work, integer *lwork,
                            integer *info);

extern void F77NAME(dgetf2)(integer *m, integer *n, doublereal *A, integer *lda, integer *ipiv,
                            integer *info);

extern void F77NAME(dgbtrf)(integer *m, integer *n, integer *KL, integer *KU, doublereal *BM,
                            integer *LDBM, integer *ipiv, integer *info);

extern void F77NAME(dgttrf)(integer *N, doublereal *DL, doublereal *D, doublereal *DU,
                            doublereal *DU2, integer *ipiv, integer *info);

extern void F77NAME(dpotrf)(char *UPLO, integer *N, doublereal *SM, integer *LDSM,
                            integer *info);

extern void F77NAME(dsytrf)(char *UPLO, integer *N, doublereal *SM, integer *LDSM,
                            integer *ipiv, doublereal *WORK, integer *LWORK, integer *info);

extern void F77NAME(dpbtrf)(char *UPLO, integer *N, integer *KD, doublereal *SBM,
                            integer *LDSM, integer *info);

extern  void F77NAME(dpttrf)(integer *N, doublereal *D, doublereal *E, integer *info);

extern void F77NAME(dgecon)(char *norm, integer *n, double *a,
                            integer *lda, double *anorm, double *rcond,
                            double *work, integer *iwork,
                            integer *info);
// ********************* LU Solve Routines ***************************

extern void F77NAME(dgetrs)(char *trans, integer *N, integer *nrhs, doublereal *A, integer *lda,
                            integer * ipiv, doublereal *b, integer *ldb, integer *info);

extern void F77NAME(dgbtrs)(char *trans, integer *N, integer *kl, integer *ku, integer *nrhs,
                            doublereal *AB, integer *ldab, integer *ipiv, doublereal *b, integer *ldb, integer *info);

extern void F77NAME(dsytrs)(char *uplo, integer *N, integer *nrhs, doublereal *A, integer *lda,
                            integer *ipiv, doublereal *b, integer *ldb, integer *info);

extern void F77NAME(dgttrs)(char *trans, integer *N, integer *nrhs, doublereal *DL,
                            doublereal *D, doublereal *DU, doublereal *DU2, integer *ipiv, doublereal *b,
                            integer *ldb, integer *info);

extern void F77NAME(dpotrs)(char *UPLO, integer *N, integer *nrhs, doublereal *A, integer *LDA,
                            doublereal *b, integer *ldb, integer *info);

extern void F77NAME(dpttrs)(integer *N, integer *nrhs, doublereal *D, doublereal *E,
                            doublereal *b, integer *ldb, integer *info);

extern void F77NAME(dpbtrs)(char *UPLO, integer *N, integer *KD, integer *nrhs, doublereal *AB,
                            integer *LDAB, doublereal *b, integer *ldb, integer *info);

// ********************* Eigen Solve Routines ***************************

extern void F77NAME(dsyev)(char *jobz, char *uplo, integer *N, doublereal *S,
                           integer *lda, doublereal *eig, doublereal *work, integer *lwork, integer *info);

// ********************* Eigenvalue/Singular Value Decomposition Drivers

extern void F77NAME(dsyevd)(char *jobz, char *uplo, integer *n, double *a, integer *lda, double *w, integer *info);
extern void F77NAME(dsyevr)(char *jobz, char *range, char *uplo, integer *n, double *a, integer *lda, double *vl, double *vu, integer *il, integer *iu, double *abstol, integer *m, double *w, double *z, integer *ldz, integer *isuppz, integer *info);

extern void F77NAME(dgesvd)(char *jobu, char *jobvt, integer *m, integer *n, double *a, integer *lda, double *sing, double *u, integer *ldu, double *vt, integer *ldvt, double *work, integer *lwork, integer *info);
extern void F77NAME(dgesdd)(char *jobz, integer *m, integer *n, double *a, integer *lda, double *s, double *u, integer *ldu, double *vt, integer *ldvt, double *work, integer *lwork, integer *iwork, integer *info);

// *******************************


#endif

#endif // BLASLAPACK_H


