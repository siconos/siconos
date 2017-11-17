/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#ifndef QP_SOLVERS_H
#define QP_SOLVERS_H

/*! \page QPSolvers Quadratic Programming problems (QP)
  \section qpIntro The problem

  Minimize: \n
  \f[
  \frac{1}{2} x' C x + d' x
  \f]

  subject to:\n
  \f{eqnarray*}
  A(j)*x  +  b(j)   =  0 & , & j=1,...,me \\
  A(j)*x  +  b(j)  >=  0 & , & j=me+1,...,m \\
  xl  <=  x  <=  xu
  \f}

  \section qpSolversList Available solvers

  The qp pack is not yet implemented.
  The only available function is ql0001() (fortran subroutine)

  (see the functions/solvers list in QP_Solvers.h)

*/

/*!\file QP_Solvers.h
  Subroutines for the resolution of QP problems.\n

  \author Franck Perignon
*/

#ifdef __cplusplus
extern "C"
{
#endif

  /**
    solution of quadratic programming problems
    ql0001 solves the quadratic programming problem

    minimize        .5*x'*C*x + d'*x
    subject to      A(j)*x  +  b(j)   =  0  ,  j=1,...,me
    A(j)*x  +  b(j)  >=  0  ,  j=me+1,...,m
    xl  <=  x  <=  xu

    Here C must be an n by n symmetric and positive matrix, d an n-dimensional
    vector, A an m by n matrix and b an m-dimensional vector. The above
    situation is indicated by iwar(1)=1. Alternatively, i.e. if iwar(1)=0,
    the objective function matrix can also be provided in factorized form.
    In this case, C is an upper triangular matrix.

    The subroutine reorganizes some data so that the problem can be solved
    by a modification of an algorithm proposed by Powell (1983).


    usage:

    call ql0001(m,me,mmax,n,nmax,mnn,c,d,a,b,xl,xu,x,u,iout,ifail,
    iprint,war,lwar,iwar,liwar)

    Definition of the parameters:

    \param m :        total number of constraints.
    \param me :       number of equality constraints.
    \param mmax :     row dimension of a. mmax must be at least one and greater
    than m.
    \param n :        number of variables.
    \param nmax :     row dimension of C. nmax must be greater or equal to n.
    \param mnn :      must be equal to m + n + n.
    \param c(nmax,nmax): objective function matrix which should be symmetric and
    positive definite. If iwar(1) = 0, c is supposed to be the
    choleskey-factor of another matrix, i.e. c is upper
    triangular.
    \param d(nmax) :  contains the constant vector of the objective function.
    \param a(mmax,nmax): contains the data matrix of the linear constraints.
    \param b(mmax) :  contains the constant data of the linear constraints.
    \param xl(n),xu(n): contain the lower and upper bounds for the variables.
    \param x(n) :     on return, x contains the optimal solution vector.
    \param u(mnn) :   on return, u contains the lagrange multipliers. The first
    m positions are reserved for the multipliers of the m
    linear constraints and the subsequent ones for the
    multipliers of the lower and upper bounds. On successful
    termination, all values of u with respect to inequalities
    and bounds should be greater or equal to zero.
    \param iout :     integer indicating the desired output unit number, i.e.
    all write-statements start with 'write(iout,... '.
    \param ifail :    shows the termination reason.
    ifail = 0 :   successful return.
    ifail = 1 :   too many iterations (more than 40*(n+m)).
    ifail = 2 :   accuracy insufficient to satisfy convergence
    criterion.
    ifail = 5 :   length of a working array is too short.
    ifail > 10 :  the constraints are inconsistent.
    \param iprint :   output control.
    iprint = 0 :  no output of ql0001.
    iprint > 0 :  brief output in error cases.
    \param war(lwar) : real working array. the length lwar should be grater than
    3*nmax*nmax/2 + 10*nmax + 2*mmax.
    \param iwar(liwar): integer working array. the length liwar should be at
    least n.
    if iwar(1)=1 initially, then the cholesky decomposition
    which is required by the dual algorithm to get the first
    unconstrained minimum of the objective function, is
    performed internally. otherwise, i.e. if iwar(1)=0, then
    it is assumed that the user provides the initial fac-
    torization by himself and stores it in the upper trian-
    gular part of the array c.
    a named common-block  /cmache/eps   must be provided by the user,
    where \param eps defines a guess for the underlying machine precision.

    \author (c): k. schittkowski,
    mathematisches institut,
    universitaet bayreuth,
    95440 bayreuth,
    germany, f.r.

    y version:    1.5  (june, 1991)

  */

  void ql0001_(int *m , int *me , int *mmax , int *n , int *nmax , int *mnn ,
                       double *c , double *d , double *a , double *b , double *xl , double *xu ,
                       double *x , double *u , int *iout , int *ifail , int *iprint , double *war ,
                       int *lwar , int *iwar , int *liwar , double *eps);

#ifdef __cplusplus
}
#endif

#endif
