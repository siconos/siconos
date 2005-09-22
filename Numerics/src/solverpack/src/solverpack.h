#ifndef SOLVERPACK_H
#define SOLVERPACK_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*!\file SiconosNumerics.h
 * \author Nineb Sheherazade and Dubois Frederic.
 * Last Modifications : Mathieu Renouf
 *
 * !\struct method_pr
 *
 * \brief A type definition for a structure method_pr.
 *
 * \param name       name of the solver.
 * \param itermax    maximum number of iterations.
 * \param tol        convergence criteria value.
 * \param k_latin    latin coefficient
 * \param a          lower bound
 * \param a          upper bound
 * \param normType   name norm
 *
 * \param iter       final number of iteration
 * \param err        final value of error criteria
 *
 */

typedef struct
{

  char   name[64];
  int    itermax;
  double tol;
  double k_latin;
  double *a;
  double *b;
  int chat;
  char   normType[64];
  int    iter;
  double err;

} method_pr;

/*!\struct method_dr
 *
 * \brief A type definition for a structure method_dr.
 *
 * \param name       name of the solver.
 * \param itermax    maximum number of iterations.
 * \param tol        convergence criteria value.
 * \param k_latin    latin coefficient
 * \param a          lower bound
 * \param a          upper bound
 * \param normType   name norm
 *
 * \param iter       final number of iteration
 * \param err        final value of error criteria
 *
 */

typedef struct
{

  char   name[64];
  int    itermax;
  double tol;
  double k_latin;
  double *a;
  double *b;
  char   normType[64];
  int    iter;
  double err;

} method_dr;

/*!\struct method_lcp
 *
 * \brief A type definition for a structure method_lcp.
 *
 * \param name       name of the solver.
 * \param itermax    maximum number of iterations.
 * \param tol        convergence criteria value.
 * \param k_latin    latin coefficient
 * \param relax      relaxation coefficient
 * \param iout       output boolean ( 0 = no output log )
 * \param normType   name norm
 *
 * \param iter       final number of iteration
 * \param err        final value of error criteria
 *
 */

typedef struct
{

  char   name[64];
  int    itermax;
  double tol;
  double k_latin;
  double relax;
  int    iout;
  char   normType[64];
  int    iter;
  double err;

} method_lcp;

/*!\struct method_pfc_2D
 *
 * \brief A type definition for a structure method_pfc_2D.
 *
 * \param name       name of the solver.
 * \param itermax    maximum number of iterations.
 * \param mu         friction coefficient.
 * \param tol        convergence criteria value.
 * \param k_latin    latin coefficient
 * \param normType   name norm
 *
 * \param iter       final number of iteration
 * \param err        final value of error criteria
 *
 */

typedef struct
{

  char   name[64];
  int    itermax;
  double tol;
  double mu;
  double k_latin;
  int    iout;
  char   normType[64];
  int    iter;
  double err;

} method_pfc_2D;

/*!\struct method_dfc_2D
 *
 * \brief A type definition for a structure method_dfc_2D
 *
 * \param name       name of the solver.
 * \param itermax    maximum number of iterations.
 * \param normType   name norm
 * \param tol        convergence criteria value.
 * \param mu         friction coefficient.
 * \param k_latin    latin coefficient
 * \param
 * \param
 * \param
 * \param
 * \param
 * \param
 * \param
 * \param
 *
 * \param iter       final number of iteration
 * \param err        final value of error criteria
 *
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
  int    *ddl_i;
  int    *ddl_n;
  int    *ddl_tt;
  int    *ddl_c;
  int    dim_i;
  int    dim_c;
  int    dim_tt;
  int    dim_n;
  int    iter;
  double err;

} method_dfc_2D;

/*!\union method
 *
 * \brief A type definition for a union method.
 *
 * \param method_pr     : pr is a method_pr structure .
 * \param method_dr     : dr is a method_dr structure .
 * \param method_lcp    : lcp is a method_lpc structure .
 * \param method_pfc_2D : pfc_2D is a method_pfc_2D structure .
 * \param method_dfc_2D : dfc_2D is a method_dfc_2D structure .
 * \param method_qp     : qp is a method_qp structure .
 *
 * \param iter       final number of iteration
 * \param err        final value of error criteria
 *
 */

typedef union
{

  method_pr pr;
  method_dr dr;
  method_lcp lcp;
  method_pfc_2D pfc_2D;
  method_dfc_2D dfc_2D;

  /*!
   * \todo method_qp does not exist
   */

} method;


/*
 * header for C++ compiling / and C compiling
 */

#ifdef __cplusplus

//extern "C" {

/* body of header */

/**************** LCP *********************/

extern "C" int lcp_solver(double *vec, double *q , int *n , method *pt , double *z , double *w);

extern "C" int lcp_solver_block(int *inb , int *iid , double *vec, double *q , int *nn , int *nb , method *pt , double *z ,
                                double *w , int *it_end , int *itt_end , double *res);


extern "C" void lcp_lemke(double *vec, double *qqq, int *nn, int *itermax, double *zlem,
                          double *wlem, int *it_end, double *res, int *info);

extern  "C" void lcp_qp(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                        int *iparamLCP , double *dparamLCP);

extern  "C" void lcp_cpg(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                         int *iparamLCP , double *dparamLCP);

extern  "C" void lcp_nlgs(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                          int *iparamLCP , double *dparamLCP);

extern  "C" void lcp_nsqp(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                          int *iparamLCP , double *dparamLCP);

extern  "C" void lcp_latin(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                           int *iparamLCP , double *dparamLCP);

extern  "C" void lcp_lexicolemke(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                                 int *iparamLCP , double *dparamLCP);

extern  "C" void lcp_newton_min(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                                int *iparamLCP , double *dparamLCP);

/********************************************/

extern "C" int dr_solver(double* , double* , int* , method* , double* , double*);

extern "C" void dr_latin(double * , double *, int *, double * , double *, double *, int *, double *, double* , double* , int *, double *, int *)  ;

extern "C" void dr_nlgs(double *vec , double *q , int *nn , double *a , double *b , int *itermax , double *tol ,
                        double *z , double *w , int *it_end , double *res , int *info);

/********************************************/

extern "C" int dfc_2D_solver(double* , double* , int* , method* , double* , double*);

extern "C" void dfc_2D_latin(double* , double* , int* , double* , double* , int* , double* , double* , double* , int* , double* , int*);

/********************************************/

extern "C" int pfc_2D_solver(double *vec , double *q , int *n , method *pt , double *z , double *w);

extern "C" void pfc_2D_cpg(int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP);

extern "C" void pfc_2D_nlgs(int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP);

extern "C" void pfc_2D_latin(int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP);

/********************************************/

extern "C" int pr_solver(double* , double* , int* , method* , double* , double*);

extern "C" void pr_latin(double* , double* , int* , double* , double* , double* , int* ,
                         double* , int *, double* , double* , int* , double* , int*);

extern "C" void pr_nlgs(double* , double* , int* , double* , double* , int* , double* , int*, double* , double* , int* , double* , int *);

/********************************************/

extern "C" void dfc_2D2lcp(int *, double *, method *, double *, int *, int *, int * , int *, int *, int *,
                           int *, int *, double * , double * , int *, double *, double *);

extern "C" void lcp2dfc_2D(int *, double *, double *, method *, double *, double *, int *, double *, int *,
                           int *, int *, int *,  int *, int *, int *, double *, double *);

/****************************** **************** ************************************/

//extern "C" void pfc_2D_projc( int nc , double mu , double *z , double *p , int *status );

//extern "C" void pfc_2D_projf( int n , double *ww , double *zz , double *rr , double *pp , int *status )

#endif

#ifndef __cplusplus

/*extern { */

/**@defgroup group1 LCP (Linear Complementary Problem)
 * @{
 *
 * \fn int extern lcp_solver( double *vec , double *q ,int *n , method *pt , double *z , double *w , int *it_end , double *res )
 * \brief lcp_solver.c is a generic interface allowing the call of one of the @ref lcp solvers.
 *
 */

extern int lcp_solver(double *vec, double *q , int *n , method *pt , double *z , double *w);

/*
 * @}
 *
 * @page lcp
 *
 * The C routines that solve LCP:
 *
 * lcp_nlgs.c
 *
 * lcp_cpg.c
 *
 * lcp_latin.c
 *
 * lcp_lemke.c
 *
 * lcp_lexicolemke.c
 *
 * lcp_qp.c
 *
 * lcp_qpnonsym.c
 *
 * lcp_newtonmin.c
 *
 *
 * @defgroup group2 Block LCP (Linear Complementary Problem)
 * @{
 *
 * \fn int extern lcp_solver_block( int *inb , int *iid , double *vec, double *q ,\n
 *                                  int *nn , int *nb , method *pt , double *z ,\n
 *                  double *w , int *it_end , int *itt_end , double *res )
 *
 * \brief lcp_solver_block.c is a generic interface for block matrices allowing the call of one of the @ref lcp solvers.
 *
 */

extern int lcp_solver_block(int *inb , int *iid , double *vec, double *q , int *nn , int *nb , method *pt , double *z ,
                            double *w , int *it_end , int *itt_end , double *res);

/*!
 * @}
 *
 * @defgroup group3 PR (Primal Relay)
 * @{
 *
 * \fn int extern  pr_solver ( double* , double* , int* , method* , double* , double* )
 *
 * \brief pr_solver() is a generic interface allowing the call of one of the @ref pr solvers.
 */

extern int pr_solver(double* , double* , int* , method* , double* , double*);

/*!
 * @}
 *
 * @page pr
 *
 * The C routines that solve PR:
 *
 * pr_latin.c
 *
 * pr_nlgs.c
 *
 *
 * @defgroup group4 DR (Dual Relay)
 * @{
 *
 * \fn int extern  dr_solver( double* , double* , int* , method* , double* , double* )
 *
 * \brief dr_solver() is a generic interface allowing the call of one of the @ref dr solvers.
 *
 */

extern int dr_solver(double* , double* , int* , method* , double* , double*);

/*!
 * @}
 *
 * @page dr
 *
 * The C routines that solve DR:
 *
 * dr_latin.c
 *
 * dr_nlgs.c
 *
 * @defgroup group5 PFC (Primal Frictional Contact)
 * @{
 *
 * \fn int extern  pfc_2D_solver( double *vec , double *q , int *n , method *pt , double *z , double *w , int *it_end , double *res )
 *
 * \brief pfc_2D_solver() is a generic interface allowing the call of one of the @ref 2D pfc solvers.
 *
 */

extern int pfc_2D_solver(double *vec , double *q , int *n , method *pt , double *z , double *w);

/*!
 * @}
 *
 *
 * @page pfc
 *
 * The C routines that solve PFC:
 *
 * pfc_2D_latin.c
 *
 * pfc_2D_nlgs.c
 *
 * pfc_2D_cpg.c
 *
 * @defgroup group6 DFC (Dual Frictional Contact)
 * @{
 *
 * \fn int extern dfc_2D_solver( double *vec , double *q , int *n , method *pt , double *z , double *w )
 *
 * \brief dfc_2D_solver() is a generic interface allowing the call of one of the @ref dfc solvers.
 *
 */

extern int dfc_2D_solver(double *vec , double *q , int *n , method *pt , double *z , double *w);

/*!
 * @}
 *
 * @page dfc
 *
 * The C routines that solve DFC:
 *
 * dfc_2D_latin.c
 *
 */

/*!
 * \todo solve_qp does not exist
 */

/*********************************** LCP *****************************************/

extern void lcp_lemke(double *vec, double *qqq, int *nn, int *itermax, double *zlem,
                      double *wlem, int *it_end, double *res, int *info);


extern void lcp_qp(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                   int *iparamLCP , double *dparamLCP);

extern void lcp_cpg(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                    int *iparamLCP , double *dparamLCP);

extern void lcp_nlgs(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                     int *iparamLCP , double *dparamLCP);

extern void lcp_nsqp(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                     int *iparamLCP , double *dparamLCP);

extern void lcp_latin(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                      int *iparamLCP , double *dparamLCP);

extern void lcp_lexicolemke(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                            int *iparamLCP , double *dparamLCP);

extern void lcp_newton_min(int *nn , double *vec , double *q , double *z , double *w , int *info ,
                           int *iparamLCP , double *dparamLCP);


/*********************************** PR *****************************************/

extern void pr_latin(double* , double* , int* , double* , double* , double* , int* ,
                     double* , int *, double* , double* , int* , double* , int*);

extern void pr_nlgs(double* , double* , int* , double* , double* , int* , double* , int*, double* , double* , int* , double* , int *);

/*********************************** DR *****************************************/

extern void dr_latin(double *, double *, int *, double * , double *, double *, int *, double *, double*, double *, int *, double *, int *)  ;

extern void dr_nlgs(double *vec , double *q , int *nn , double *a , double *b , int *itermax , double *tol ,
                    double *z , double *w , int *it_end , double *res , int *info);

/*********************************** PFC 2D *****************************************/

extern void pfc_2D_cpg(int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP);

extern void pfc_2D_nlgs(int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP);

extern void pfc_2D_latin(int *nn , double *vec , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP);

/*********************************** DFC 2D *****************************************/

extern void dfc_2D_latin(double* , double* , int* , double* , double* , int* , double* , double* , double* , int* , double* , int*);

/******************************LCP SWITCH DFC 2D ************************************/


extern void dfc_2D2lcp(int *, double *, method *, double *, int *, int *, int * , int *, int *, int *, int *,
                       int *, double * , double * , int *, double *, double *);

extern void lcp2dfc_2D(int *, double *, double *, method *, double *, double *, int *, double *, int *, int *,
                       int *, int *,  int *, int *, int *, double *, double *);

/****************************** **************** ************************************/

//extern void pfc_2D_projc( int nc , double mu , double *z , double *p , int *status );

//extern void pfc_2D_projf( int n , double *ww , double *zz , double *rr , double *pp , int *status )

/******************************LCP SWITCH DFC 2D ************************************/

#endif

#endif // SOLVERPACK_H
