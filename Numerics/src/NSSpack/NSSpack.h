/* Siconos-Numerics version 2.1.0, Copyright INRIA 2005-2006.
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
#ifndef SOLVERPACK_H
#define SOLVERPACK_H

#include "blaslapack.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*!\file solverpack.h
  \author Nineb Sheherazade and Dubois Frederic.
  Last Modifications : Mathieu Renouf , Pascal Denoyelle
 */



/*!\struct method_pr
   \brief A type definition for a structure method_pr.

  \param name       name of the solver.
  \param itermax    maximum number of iterations.
  \param tol        convergence criteria value.
  \param k_latin    latin coefficient
  \param a          upper bound.
  \param b          lower bound.
  \param chat       output boolean ( 0 = no output log ).
  \param normType   name norm (not yet available).

  \param iter       final number of iterations.
  \param err        final value of error criteria.
*/

typedef struct
{

  char     name[64];
  int      itermax;
  double   tol;
  double   k_latin;
  double   *a;
  double   *b;
  int      chat;
  char     normType[64];
  int      iter;
  double   err;

} method_pr;

/*!\struct method_dr

  \brief A type definition for a structure method_dr.

  \param name       name of the solver.
  \param itermax    maximum number of iterations.
  \param tol        convergence criteria value.
  \param k_latin    latin coefficient
  \param a          upper bound
  \param b          lower bound
  \param chat       output boolean ( 0 = no output log ).
  \param normType   name norm (not yet available).

  \param iter       final number of iterations
  \param err        final value of error criteria
*/

typedef struct
{

  char     name[64];
  int      itermax;
  double   tol;
  double   k_latin;
  double   *a;
  double   *b;
  int      chat;
  char     normType[64];
  int      iter;
  double   err;

} method_dr;

/*!\struct method_lcp

  \brief A type definition for a structure method_lcp.

  \param name       name of the solver.
  \param itermax    maximum number of iterations.
  \param tol        convergence criteria value.
  \param k_latin    latin coefficient.
  \param relax      relaxation coefficient.
  \param rho        regularization coefficient
  \param chat       output boolean ( 0 = no output log ).
  \param normType   name norm (not yet available).

  \param iter       final number of iterations.
  \param err        final value of error criteria.
*/

typedef struct
{

  char   name[64];
  int    itermax;
  double tol;
  double k_latin;
  double relax;
  double rho;
  int    chat;
  char   normType[64];
  int    iter;
  double err;

} method_lcp;

/*!\struct method_pfc_2D

  \brief A type definition for a structure method_pfc_2D.

  \param name       name of the solver.
  \param itermax    maximum number of iterations.
  \param mu         friction coefficient.
  \param tol        convergence criteria value.
  \param k_latin    search direction of the latin metod.
  \param chat       output boolean ( 0 = no output log ).
  \param normType   name norm (not yet available).
  \param iter       final number of iterations.
  \param err        final value of error criteria
*/

typedef struct
{

  char   name[64];
  int    itermax;
  double tol;
  double mu;
  double k_latin;
  int    chat;
  char   normType[64];
  int    iter;
  double err;

} method_pfc_2D;

/*!\struct method_pfc_3D

  \brief A type definition for a structure method_pfc_3D.

  \param name       name of the solver.
  \param itermax    maximum number of iterations.
  \param mu         friction coefficient.
  \param tol        convergence criteria value.
  \param k_latin    search direction of the latin metod.
  \param chat       output boolean ( 0 = no output log ).
  \param normType   name norm (not yet available).
  \param iter       final number of iterations.
  \param err        final value of error criteria
*/
typedef struct
{

  char   name[64];
  int    itermax;
  double tol;
  double mu;
  double k_latin;
  int    chat;
  char   normType[64];
  int    iter;
  double err;

} method_pfc_3D;

/*!\struct method_dfc_2D

  \brief A type definition for a structure method_dfc_2D

  \param name       name of the solver.
  \param itermax    maximum number of iterations.
  \param normType   name norm (not yet available).
  \param tol        convergence criteria value.
  \param mu         friction coefficient.
  \param k_latin    latin coefficient
  \param J1         gap in normal contact direction.
  \param ddl_n      the contact in normal direction dof (not prescribed),
  \param ddl_tt     the contact in tangential direction dof (not prescribed)

  \param ddl_d      the prescribed dof.
  \param dim_tt     the dimension of the vector ddl_tt.
  \param dim_d      the dimension of the vector ddl_d.
  \param chat       output boolean ( 0 = no output log ).
  \param iter       final number of iteration
  \param err        final value of error criteria
*/

typedef struct
{

  char   name[64];
  int    itermax;
  char   normType[64];
  double tol;
  double mu;
  double k_latin;

  double *J1;
  int    *ddl_n;
  int    *ddl_tt;
  int    *ddl_d;
  int    dim_tt;
  int    dim_d;

  int   chat;
  int    iter;
  double err;

} method_dfc_2D;

/*!\union method

  \brief A type definition for a union method.

  \param method_pr     : pr is a method_pr structure .
  \param method_dr     : dr is a method_dr structure .
  \param method_lcp    : lcp is a method_lcp structure .
  \param method_pfc_2D : pfc_2D is a method_pfc_2D structure .
  \param method_dfc_2D : dfc_2D is a method_dfc_2D structure .
  \param method_qp     : qp is a method_qp structure (not yet available).


*/

typedef union
{

  method_pr  pr;
  method_dr  dr;
  method_lcp lcp;
  method_pfc_2D pfc_2D;
  method_pfc_3D pfc_3D;
  method_dfc_2D dfc_2D;

  /*!
   * \todo method_qp does not exist
   */

} method;

/*!\struct SparseBlockStructuredMatrix

    \brief To store sparse block matrices with square diagonal blocks

    \param nbblocks         : the total number of non null blocks
    \param **block          : *block contains the double values of one block in Fortran storage (column by column)
                              **block is the list of non null blocks
    \param size             : the number of blocks along a row (or column)
    \param *blocksize       : the list of the sizes of diagonal (square) blocks
    \param *RowIndex        : the list of *block row indices (first row = 0)
    \param *ColumnIndex     : the list of *block column indices (first column = 0)
*/

typedef struct
{
  int nbblocks;
  double **block;
  int size;
  int *blocksize;
  int *RowIndex;
  int *ColumnIndex;
} SparseBlockStructuredMatrix;

/*
 * header for C++ compiling / and C compiling
 */

#ifdef __cplusplus

/* body of header */

/*  **************** LCP ********************       */

extern "C" int lcp_solver(double *vec, double *q , int *n , method *pt , double *z , double *w);

/*
extern "C" int lcp_solver_block( int *inb , int *iid , double *vec, double *q , int *nn , int *nb , method *pt , double *z ,
         double *w , int *it_end , int *itt_end , double *res );
*/
extern "C" int lcp_solver_block(SparseBlockStructuredMatrix *blmat, double *q, method *pt , double *z , double *w , int *it_end ,
                                int *itt_end , double *res);

extern  "C" void lcp_qp(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                        int *iparamLCP , double *dparamLCP);

extern  "C" void lcp_cpg(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                         int *iparamLCP , double *dparamLCP);

extern  "C" void lcp_pgs(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                         int *iparamLCP , double *dparamLCP);

extern  "C" void lcp_rpgs(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                          int *iparamLCP , double *dparamLCP);

extern  "C" void lcp_psor(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                          int *iparamLCP , double *dparamLCP);

extern  "C" void lcp_nsqp(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                          int *iparamLCP , double *dparamLCP);

extern  "C" void lcp_latin(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                           int *iparamLCP , double *dparamLCP);


extern  "C" void lcp_latin_w(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                             int *iparamLCP , double *dparamLCP);


extern  "C" void lcp_lexicolemke(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                                 int *iparamLCP , double *dparamLCP);

extern  "C" void lcp_newton_min(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                                int *iparamLCP , double *dparamLCP);

extern  "C" void lcp_newton_FB(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                               int *iparamLCP , double *dparamLCP);

extern  "C" int filter_result_LCP(int n, double *vec , double *q , double *z , double tol, int chat, double *w);

extern  "C" int lcp_compute_error(int n, double *vec , double *q , double *z , int chat, double *w, double error);

extern  "C" int filter_result_LCP_block(SparseBlockStructuredMatrix *blmat, double *q , double *z , double tol, int chat, double *w);

extern  "C" void freeSpBlMat(SparseBlockStructuredMatrix *blmat);

/* ******************************************* */

extern "C" int dr_solver(double* , double* , int* , method* , double* , double*);

extern "C" void dr_latin(double * , double *, int *, double * , double *, double *, int *, double *, int *, double* , double* , int *, double *, int *)  ;

extern "C" void dr_nlgs(double *vec , double *q , int *nn , double *a , double *b , int *itermax , double *tol , int *chat,
                        double *z , double *w , int *it_end , double *res , int *info);

/* ******************************************* */

extern "C" int dfc_2D_solver(double* , double* , int* , method* , double* , double*);

extern "C" void dfc_2D_latin(double* , double* , int* , double* , double* , int* , double* , int *, double* , double* , int* , double* , int*);


/* ******************************************* */
extern "C" int pfc_2D_solver(double *vec , double *q , int *n , method *pt , double *z , double *w);

extern "C" void pfc_2D_cpg(int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP);

extern "C" void pfc_2D_nlgs(int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP);

extern "C" void pfc_2D_latin(int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP);


/* ******************************************* */
extern "C" int pfc_3D_solver(double *vec , double *q , int *n , method *pt , double *z , double *w);

extern "C" void pfc_3D_cpg(int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP);

extern "C" void pfc_3D_nlgs(int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP);

extern "C" void pfc_3D_nlgsnewton(int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP);

extern "C" void pfc_3D_newtonfunction(int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP);

/* ******************************************* */
extern "C" int pr_solver(double* , double* , int* , method* , double* , double*);

extern "C" void pr_latin(double* , double* , int* , double* , double* , double* , int* ,
                         double* , int *, double* , double* , int* , double* , int*);

extern "C" void pr_nlgs(double* , double* , int* , double* , double* , int* , double* , int*, double* , double* , int* , double* , int *);


/* ******************************************* */

extern "C" void dfc_2D2lcp(int *, double *, double *, double *, int *, int *, int * , int *, int *, double * , double *, double *);





extern "C" void lcp2dfc_2D(int *, double *, double *, double *, double *, double *,  int *, int *,
                           int *, int *, int *,  double *, double *);



extern "C" void dfc_2Dcond_2D(int *, double *, double *, double *, int *, int *, int * , int *, int *, double * , double *, double *);





extern "C" void cond_2D2dfc_2D(int *, double *, double *, double *, double *, double *,  int *, int *,
                               int *, int *, int *,  double *, double *);





/* ******************************************* */

/*extern "C" void pfc_2D_projc( int nc , double mu , double *z , double *p , int *status );

extern "C" void pfc_2D_projf( int n , double *ww , double *zz , double *rr , double *pp , int *status )
*/


#endif

#ifndef __cplusplus


/**@defgroup group1 LCP (Linear Complementary Problem)
   @{
*/

/** \fn int extern lcp_solver( double *vec , double *q ,int *n , method *pt , double *z , double *w , int *it_end , double *res )
 *  \brief lcp_solver.c is a generic interface allowing the call of one of the @ref lcp solvers.
*/

/** @brief

 lcp_solver.c is a generic interface allowing the call of one of the @ref lcp solvers.
*/

extern int lcp_solver(double *vec, double *q , int *n , method *pt , double *z , double *w);

/**@}*/

/**@page lcp

The C routines that solve LCP:



lcp_pgs.c

lcp_rpgs.c

lcp_latin.c

lcp_latin_w.c

lcp_cpg.c

lcp_lexicolemke.c

lcp_qp.c

lcp_nsqp.c

lcp_newton_min.c
*/

/**@defgroup group2 Block LCP (Linear Complementary Problem)
   @{
*/

/** \fn extern int lcp_solver_block(SparseBlockStructuredMatrix *blmat, double *q, method *pt , double *z , double *w , int *it_end , int *itt_end ,double *res );

 * \brief lcp_solver_block() is a generic interface for block matrices allowing the call of one of the @ref lcp solvers.
*/

/** @brief
 lcp_solver_block() is a generic interface for block matrices allowing the call of one of the @ref lcp solvers.

*/
/*
extern int lcp_solver_block( int *inb , int *iid , double *vec, double *q , int *nn , int *nb , method *pt , double *z ,double *w , int *it_end , int *itt_end , double *res );
*/

extern int lcp_solver_block(SparseBlockStructuredMatrix *blmat, double *q, method *pt , double *z , double *w , int *it_end ,
                            int *itt_end , double *res);


/**@}*/




/**@defgroup group3 PR (Primal Relay)
   @{
*/

/** \fn int extern  pr_solver ( double* , double* , int* , method* , double* , double* )

 * \brief pr_solver.c is a generic interface allowing the call of one of the @ref pr solvers.
 */

/** @brief
pr_solver.c is a generic interface allowing the call of one of the @ref pr solvers.

 */
extern int pr_solver(double* , double* , int* , method* , double* , double*);

/**@}*/

/** @page pr

  The C routines that solve PR:

  pr_latin.c

  pr_nlgs.c

*/

/**@defgroup group4 DR (Dual Relay)
 @{
 */

/** \fn int extern  dr_solver( double* , double* , int* , method* , double* , double* )

 * \brief dr_solver.c is a generic interface allowing the call of one of the @ref dr solvers.

 */

/** @brief
dr_solver.c is a generic interface allowing the call of one of the @ref dr solvers.
*/

extern int dr_solver(double* , double* , int* , method* , double* , double*);
/**@}*/


/**@page dr

  The C routines that solve DR:

  dr_latin.c

  dr_nlgs.c
*/

/**@defgroup group5 2D PFC (Two-dimensional Primal Frictional Contact)
   @{
*/

/** \fn int extern  pfc_2D_solver( double *vec , double *q , int *n , method *pt , double *z , double *w )

 * \brief pfc_2D_solver.c is a generic interface allowing the call of one of the @ref pfc_2D solvers.

 */
/** @brief
pfc_2D_solver.c is a generic interface allowing the call of one of the @ref pfc_2D solvers.
*/
extern int pfc_2D_solver(double *vec , double *q , int *n , method *pt , double *z , double *w);

/**@}*/

/**@page pfc_2D

  The C routines that solve PFC:

  pfc_2D_latin.c

  pfc_2D_nlgs.c

  pfc_2D_cpg.c
*/

/**@defgroup group6 3D PFC (Three-dimensional Primal Frictional Contact)
   @{
*/

/** \fn int extern  pfc_3D_solver( double *vec , double *q , int *n , method *pt , double *z , double *w )

 * \brief pfc_3D_solver() is a generic interface allowing the call of one of the @ref pfc_3D solvers.

 */
/** @brief
pfc_3D_solver() is a generic interface allowing the call of one of the @ref  pfc_3D solvers.
*/
extern int pfc_3D_solver(double *vec , double *q , int *n , method *pt , double *z , double *w);

/**@}*/

/**@page pfc_3D

  The C routines that solve 3D PFC:

  pfc_3D_nlgs.c

  pfc_3D_nlgsnewton.c

  pfc_3D_newtonfunction.c

  pfc_3D_cpg.c

*/

/**@defgroup group7 2D DFC (Two-Dimensional Dual Frictional Contact)
   @{
 */

/** \fn int extern dfc_2D_solver( double *vec , double *q , int *n , method *pt , double *z , double *w )

 * \brief dfc_2D_solver() is a generic interface allowing the call of one of the @ref dfc solvers.

 */
/** @brief
dfc_2D_solver() is a generic interface allowing the call of one of the @ref dfc solvers.
*/
extern int dfc_2D_solver(double *vec , double *q , int *n , method *pt , double *z , double *w);

/**@}*/

/**@page dfc_2D

  The C routines that solve DFC:

  dfc_2D_latin.c

 */



/*!
 * \todo solve_qp does not exist
 */

/*  *********************************** LCP ******************************* */
extern void lcp_qp(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                   int *iparamLCP , double *dparamLCP);

extern void lcp_cpg(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                    int *iparamLCP , double *dparamLCP);

extern void lcp_pgs(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                    int *iparamLCP , double *dparamLCP);

extern void lcp_rpgs(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                     int *iparamLCP , double *dparamLCP);

extern void lcp_nsqp(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                     int *iparamLCP , double *dparamLCP);

extern void lcp_psor(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                     int *iparamLCP , double *dparamLCP);

extern void lcp_latin(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                      int *iparamLCP , double *dparamLCP);

extern void lcp_latin_w(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                        int *iparamLCP , double *dparamLCP);


extern void lcp_lexicolemke(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                            int *iparamLCP , double *dparamLCP);

extern void lcp_newton_min(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                           int *iparamLCP , double *dparamLCP);

extern void lcp_newton_FB(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                          int *iparamLCP , double *dparamLCP);

extern  int filter_result_LCP(int n, double *vec , double *q , double *z , double tol, int chat, double *w);

extern  int filter_result_LCP_block(SparseBlockStructuredMatrix *blmat, double *q , double *z , double tol, int chat, double *w);

extern  void freeSpBlMat(SparseBlockStructuredMatrix *blmat);

/* **************************** PR **************************************** */

extern void pr_latin(double* , double* , int* , double* , double* , double* , int* ,
                     double* , int *, double* , double* , int* , double* , int*);

extern void pr_nlgs(double* , double* , int* , double* , double* , int* , double* , int*, double* , double* , int* , double* , int *);

/* ********************************** DR ********************************** */

extern void dr_latin(double *, double *, int *, double * , double *, double *, int *, double *, int *, double*, double *, int *, double *, int *)  ;

extern void dr_nlgs(double *vec , double *q , int *nn , double *a , double *b , int *itermax , double *tol , int *chat,
                    double *z , double *w , int *it_end , double *res , int *info);

/* ********************************** PFC 2D *************************** */

extern void pfc_2D_cpg(int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP);

extern void pfc_2D_nlgs(int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP);

extern void pfc_2D_latin(int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP);

/*  ******************************* PFC 3D **************************** */

extern void pfc_3D_cpg(int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP);

extern void pfc_3D_nlgs(int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP);

extern void pfc_3D_nlgsnewton(int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP);

extern void pfc_3D_newtonfunction(int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP);

/* ************************* DFC 2D **************************************** */

extern void dfc_2D_latin(double* , double* , int* , double* , double* , int* , double* , int *, double* , double* , int* , double* , int*);

/* *****************************LCP SWITCH DFC 2D *********************** */



extern  void dfc_2D2lcp(int *, double *, double *, double *, int *, int *, int * , int *, int *, double * , double *, double *);





extern  void lcp2dfc_2D(int *, double *, double *, double *, double *, double *,  int *, int *,
                        int *, int *, int *,  double *, double *);

/* *****************************COND SWITCH DFC 2D ********************* */

extern  void dfc_2D2cond_2D(int *, double *, double *, double *, int *, int *, int * , int *, int *, double * , double *, double *);





extern  void cond_2D2dfc_2D(int *, double *, double *, double *, double *, double *,  int *, int *,
                            int *, int *, int *,  double *, double *);


/* ***************************** **************** ************************* */

/* extern void pfc_2D_projc( int nc , double mu , double *z , double *p , int *status );

extern void pfc_2D_projf( int n , double *ww , double *zz , double *rr , double *pp , int *status ) */

/* *****************************LCP SWITCH DFC 2D ********************** */

#endif

#endif /* SOLVERPACK_H */
