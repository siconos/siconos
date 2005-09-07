#ifndef SOLVERPACK_H
#define SOLVERPACK_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*!\file SiconosNumerics.h
 *   \author Nineb Sheherazade and Dubois Frederic.
 */

/*!\struct method_rp
  \brief A type definition for a structure method_rp.

  \param char * : nom_method A pointer over characters, the name of the solver.
  \param int : itermax The maximum iterations.
  \param char* : the kind of norm used for this method
  \param double : tol The tolerance value.
  \param double * :a A pointer over doubles, the bound.
  \param double : k_latin The latin coefficient.
*/

typedef struct
{
  char nom_method[64];
  int itermax;
  char normType[64];
  double tol;
  double *a;
  double *b;
  double k_latin;
} method_rp;


/*!\struct method_rd
   \brief A type definition for a structure method_rd.

   \param char* : nom_method A pointer over characters, the name of the solver.
   \param int : itermax The maximum iterations
   \param char* : the kind of norm used for this method
   \param double : tol The tolerance value.
   \param double* : a A pointer over doubles, the bound.
   \param double : k_latin The latin coefficient
*/

typedef struct
{
  char nom_method[64];
  int itermax;
  char normType[64];
  double tol;
  double *a;
  double *b;
  double k_latin;
} method_rd;

/*!\struct method_lcp
 *  \brief A type definition for a structure method_lcp.
 *
 * \param char*  : nom_method A pointer over characters, the name of the solver.
 * \param int    : itermax The maximum iterations.
 * \param double : tol The tolerance value.
 * \param double : k_latin The latin coefficient
 * \param char*  : the kind of norm used for this method
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

} method_lcp;


/*!\struct method_cfp
 * \brief A type definition for a structure method_cfp.
 *
 * \param char* : nom_method A pointer over characters, the name of the solver.
 * \param int : itermax The maximum iterations.
 * \param char* : the kind of norm used for this method
 * \param double : tol The tolerance value.
 * \param double : mu The friction coefficient.
 * \param double : k_latin The latin coefficient.
 *
 */

typedef struct
{

  char nom_method[64];
  int itermax;
  char normType[64];
  double tol;
  double mu;
  double k_latin;

} method_cfp;

/*!\struct method_cfd
 * \brief A type definition for a structure method_cfd.
 *
 * \param char* : nom_method A pointer over characters, the name of the solver.
 * \param int : itermax The maximum iterations.
 * \param char* : the kind of norm used for this method
 * \param double : tol The tolerance value.
 * \param double : mu The friction coefficient.
 * \param double : k_latin The latin coefficient
 */

typedef struct
{

  char nom_method[64];
  int itermax;
  char normType[64];
  double tol;
  double mu;
  double k_latin;
  double *J1;
  int *ddl_i;
  int *ddl_n;
  int *ddl_tt;
  int *ddl_c;
  int dim_i;
  int dim_c;
  int dim_tt;
  int dim_n;
} method_cfd;

/*!\union method
 * \brief A type definition for a union method.
 *
 *
 * \param method_rp  : rp is a method_rp structure .
 * \param method_rd  : rd is a method_rd structure .
 * \param method_lcp : lcp is a method_lpc structure .
 * \param method_cfp : cfp is a method_cfp structure .
 * \param method_cfd : cfd is a method_cfd structure .
 * \param method_qp  : qp is a method_qp structure .
 *
 */

typedef union
{

  method_rp rp;
  method_rd rd;
  method_lcp lcp;
  method_cfp cfp;
  method_cfd cfd;
  /* \todo method_qp does not exist
   */
  /* method_qp qp; */

} method;

/*
 * header for C++ compiling / and C compiling
 */

#ifdef __cplusplus

//extern "C" {

/* body of header */

/**************** LCP *********************/

extern "C" int lcp_solver(double *vec, double *q , int *n , method *pt , double *z , double *w , int *it_end , double *res);

extern "C" void lcp_lemke(double *vec, double *qqq, int *nn, int *itermax, double *zlem,
                          double *wlem, int *it_end, double *res, int *info);

extern "C" void lcp_nlgs(double *vec , double *q , int *nn , int *itermax , double *tol , double *omega , int *ispeak , double *z ,
                         double *w , int *it_end , double *res , int *info);

extern "C" void lcp_cpg(double *vec , double *q , int *nn , int *itermax , double *tol , int *ispeak , double *z ,
                        double *w , int *it_end , double * res , int *info);

extern "C" void lcp_latin(double*, double*, int *, double *, int *, double *, double [], double[], int *, double *, int *);

extern "C" void lcp_qp(double *vec , double *qq , int *nn , double *tol , double *z , double *w , int *info);

extern "C" void lcp_nsqp(double *vec , double *qq , int *nn , double *tol , double *z , double *w , int *info);

extern "C" void lcp_lexicolemke(double *vec, double *q , int *nn, int *itermax , int *ispeak , double *zlem ,
                                double *wlem, int *it_end , int *info);

/********************************************/

extern "C" int solve_rp(double* , double* , int* , method* , double* , double*);

extern "C" int solve_rd(double* , double* , int* , method* , double* , double*);

extern "C" int solve_cfp(double* , double* , int* , method* , double* , double*);

extern "C" int solve_cfd(double* , double* , int* , method* , double* , double*);

extern "C" void rp_latin(double [], double *, int *, double *, double [], double [], int *,
                         double *, double [], double [], int *, double *, int *);


extern "C" void rp_gsnl(double [], double *, int *, double *, double *, int *, double *, double[], double [], int *, double *, int *)  ;

extern "C" void rd_latin(double [], double *, int *, double * , double *, int *, double *, double[], double [], int *, double *, int *)  ;

extern "C" void cfp_gsnl(double [], double [], int *, double *, int *, double *, double [], double[], int *, double *, int*);

extern "C" void cfp_gcp(double [], double [], int *, double *, int *, double *, double [], double [], int *, double *, int *);

extern "C" void cfp_latin(double [], double [], int *, double *, double *, int *, double *, double [], double [], int *, double *, int*);

extern "C" void cfd_latin(double[], double[], int *, double*, double*, int *, double*, double[], double[], int *, double*, int*);

extern "C" int cfd_lcp(int *, double *, method *, double *, int *, int *, int * , int *, int *, int *,
                       int *, int *, double * , double * , int *, double *, double *);

extern "C" int lcp_cfd(int *, double *, double *, method *, double *, double *, int *, double *, int *,
                       int *, int *, int *,  int *, int *, int *, double *, double *);

#endif

#ifndef __cplusplus

/*extern { */
/**@defgroup group1 LCP (Linear Complementary Problem)
 * @{
 *
 * \fn int extern solve_lcp (double*,double*,int *, method *,double [],double [])
 * \brief solve_lcp.c is a generic interface allowing the call of one of the @ref lcp solvers.
 * @brief
 *
 * solve_lcp.c is a generic interface allowing the call of one of the @ref lcp solvers.
 */

extern int lcp_solver(double *vec, double *q , int *n , method *pt , double *z , double *w , int *it_end , double *res);

/**@}*/

/**@page lcp

The C routines that solve LCP:

gsnl_lcp.c

gcp_lcp.c

latin_lcp.c

lemke_lcp.c

lexicolemke_lcp.c

qp_lcp.c

qpnonsym_lcp.c

*/


/**@defgroup group2 PR (Primal Relay)
   @{
*/


/** \fn int extern  solve_rp (double*,double*,int *,method *, double [],double [])
 *  \brief solve_rp() is a generic interface allowing the call of one of the @ref pr solvers.
 */
/** @brief

    solve_rp() is a generic interface allowing the call of one of the @ref pr solvers.
*/
extern int solve_rp(double*, double*, int *, method *, double [], double []);

/**@}*/


/**@page pr

The C routines that solve PR:

rp_latin.c

rp_gsnl.c

*/



/**@defgroup group3 DR (Dual Relay)
   @{
*/


/** \fn int extern  solve_rd (double*,double*,int *,method *, double [],double [])
 *  \brief solve_rd() is a generic interface allowing the call of one of the @ref dr solvers.
 */
/** @brief

    solve_rd() is a generic interface allowing the call of one of the @ref dr solvers.
*/
extern int solve_rd(double*, double*, int *, method *, double [], double []);

/**@}*/


/**@page dr

The C routines that solve DR:

rd_latin.c

rd_gsnl.c

*/




/**@defgroup group4 PFC (Primal Frictional Contact)
   @{
*/


/** \fn int extern  solve_cfp (double*,double*,int *,method *, double [],double [])
 *  \brief solve_cfp() is a generic interface allowing the call of one of the @ref pfc solvers.
 */
/** @brief

    solve_cfp() is a generic interface allowing the call of one of the @ref pfc solvers.
*/
extern int solve_cfp(double*, double*, int *, method *, double [], double []);

/**@}*/


/**@page pfc

The C routines that solve PFC:

cfp_latin.c

cfp_gsnl.c

cfp_gcp.c

*/


/**@defgroup group5 DFC (Dual Frictional Contact)
   @{
*/


/** \fn int extern solve_cfd (double*,double*,int *, method *,double [],double [])
 *  \brief solve_cfd() is a generic interface allowing the call of one of the @ref dfc solvers.
 */
/** @brief

    solve_cfd() is a generic interface allowing the call of one of the @ref dfc solvers.
*/
extern int solve_cfd(double*, double*, int *, method *, double [], double []);
//extern double solve_cfd2 (double*, double*,int *, method *,double [],double []);
//extern double solve_cfd3 (double*, double*,int *, method *,double [],double []);


/**@}*/


/**@page dfc

The C routines that solve DFC:

cfd_latin.c

*/

/* \todo solve_qp does not exist
*/

/*********************************** LCP *****************************************/

extern void lcp_lemke(double *vec, double *qqq, int *nn, int *itermax, double *zlem,
                      double *wlem, int *it_end, double *res, int *info);

extern void lcp_nlgs(double *vec , double *q , int *nn , int *itermax , double *tol , double *omega , int *ispeak , double *z ,
                     double *w , int *it_end , double *res , int *info);

extern void lcp_cpg(double *vec , double *q , int *nn , int *itermax , double *tol , int *ispeak , double *z ,
                    double *w , int *it_end , double * res , int *info);

extern void lcp_latin(double*, double*, int *, double *, int *, double *, double [], double[], int *, double *, int *);

extern void lcp_qp(double *vec , double *qq , int *nn , double *tol , double *z , double *w , int *info);

extern void lcp_nsqp(double *vec , double *qq , int *nn , double *tol , double *z , double *w , int *info);

extern void lcp_lexicolemke(double *vec, double *q , int *nn, int *itermax , int *ispeak , double *zlem ,
                            double *wlem, int *it_end , int *info);

/*********************************** LCP *****************************************/

/* RP */

extern void rp_latin(double [], double *, int *, double *, double [], double [], int *, double *, double [], double [], int *, double *, int *);

extern int rp_gsnl(double [], double *, int *, double *, double *, int *, double *, double[], double [], int *, double *, int *)  ;

/* RD */

extern int rd_latin(double [], double *, int *, double * , double *, double *, int *, double *, double[], double [], int *, double *, int *)  ;

/* CFP */

extern cfp_gsnl(double [], double [], int *, double *, int *, double *, double [], double[], int *, double *, int*);

extern cfp_gcp(double [], double [], int *, double *, int *, double *, double [], double [], int *, double *, int *);

extern cfp_latin(double [], double [], int *, double *, double *, int *, double *, double [], double [], int *, double *, int*);

/* CFD */

extern  cfd_latin(double[], double[], int *, double*, double*, int *, double*, double[], double[], int *, double*, int*);

/* LCP SWITCH CFD */

extern int cfd_lcp(int *, double *, method *, double *, int *, int *, int * , int *, int *, int *, int *, int *, double * , double * , int *, double *, double *);

extern int lcp_cfd(int *, double *, double *, method *, double *, double *, int *, double *, int *, int *, int *, int *,  int *, int *, int *, double *, double *);


//extern cfp_gsnl (double [],double [],int *,double *,double *,double *, double [],double[],double *,double *,double*);

//extern int cfd_lcp26 (int *,double *, method *,double *,int *, int *,int * ,int *,int *, int *, int *,int *,double * , double * ,int *, double *,double *);

//extern void cfd_lcpfin (int *,double *, method *,double *,int *, int *,int * ,int *,int *, int *, int *,int *,double * , double * ,int *, double *,double *);

//extern lcp_cfd (int *,double[],double[],double[],double[]);

//extern cfp_gcp (double [],double [],int *,double *,double *,double *, double [],double [],double*,double *,double*);

//extern cfp_latin (double [],double [],int *,double *,double *,double *,double *,double [],double [], double *,double *,double*);

//extern cfd_latin2 (double[],double[],int *,double*,double*,double *,double*,double[],double[], double *,double*,double*);

//extern void lcp_cfd (int *,double *,double *,method *,double *, double *, double *,int *, int *, int *, int *, double *,double *);

//extern solve_qp (double [][maxcols],double [], int *, method *, double [],double []);

//extern int rp_latin(double [],double *,int *, double * ,double *,int *, double *,double[],double [],int *,double *,int *)  ;

//extern cfd_latin (double[],double[],int *,double*,double*,int *,double*,double[],double[], int *,double*,int*);

//extern cfd_lcp (int *,double*,double[],double[],double[],double[]);

//extern void cfd_lcp (int *,double *, method *,double *,int *,int *,int *,int *,double *, double *, double *,double *);

//extern void cfd_lcp (int *,double *, method *,double *,int *, int * ,int *,int * ,int *, int *, int *,int *,double *, double *, int *, double *,double *);

//lemke_lcp2 (double [],double [],int *,double *,double[],double[],double *, double *,double *);


//gsnl_lcp (double [],double [],int *,double *,double *,double [],double [],double*,double *,double *);

//extern gcp_lcp2 (double [],double [],int *,int *,double *,double [],double[],int *,double *,int *);

//extern int lcp_cfd26 (int *,double *,double *,method *,double *, double *, int *, double *,int *, int *, int *, int *,  int *, int *, int *, double *,double *);
//} /* closing brace for extern */
#endif

#endif // SOLVERPACK_H
