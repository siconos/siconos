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
#ifndef PFC_3D_SOLVERS_H
#define PFC_3D_SOLVERS_H

/*!\file pfc_3D_solvers.h
  Typedef and functions declarations related to primal frictional contact problems.
  \author Nineb Sheherazade and Dubois Frederic.
  Last Modifications : Mathieu Renouf , Pascal Denoyelle, Franck Perignon
*/

#include "NSSTools.h"

/**@defgroup group6 3D PFC (Three-dimensional Primal Frictional Contact)
   @{
*/

/** \fn int extern  pfc_3D_solver( double *vec , double *q , int *n , method *pt , double *z , double *w )

* \brief pfc_3D_solver() is a generic interface allowing the call of one of the @ref pfc_3D solvers.

*/
/** @brief
    pfc_3D_solver() is a generic interface allowing the call of one of the @ref  pfc_3D solvers.
*/
/**@}*/

/**@page pfc_3D

The C routines that solve 3D PFC:

pfc_3D_nlgs.c

pfc_3D_nlgsnewton.c

pfc_3D_cpg.c

pfc_3D_nsgs.c

pfc_3D_newton.c

pfc_3D_projection.c

*/

/*!\struct method_pfc_3D

\brief A type definition for a structure method_pfc_3D.

\param name       name of the solver.
\param itermax    maximum number of iterations.
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
  int    local_formulation;
  int    local_solver;
  int    local_itermax;
  double local_tol;
  double k_latin;
  int    chat;
  char   normType[64];
  int    iter;
  double err;
  int    local_iter;
  double local_err;

} method_pfc_3D;

/** Pointer to function type, used in pfc3D */
typedef void (*pfc3D_fPtr)(int, double*, double*, double*, double*, double*, double*, double*, double*, double, double, double);

/** Pointer to function type, used in pfc3D */
typedef void (*Linesearch_function)(int, double*, double*, double*, double*, double*, double*, double*, double*, double, double, double, double);

/** Local solver pointer to function for pfc_3D problems. */
typedef void (*PFC3D_local_solver)(int, double*, double*, double*, double*, double, pfc3D_fPtr*, pfc3D_fPtr*, double*, double*, double*, int*, double*);

#ifdef __cplusplus
extern "C" {
#endif

  /************ PFC_3D Related Functions and Definitions ******************************* */

  /** conjugate-projected gradient pfc3D solver */
  void pfc_3D_cpg(int, double *vec , double *q , double *z , double *w , double *mu, int *info , int *iparamLCP , double *dparamLCP);

  /** Non linear Gauss-Seidel pfc3D solver */
  void pfc_3D_nlgs(int, double *vec , double *q , double *z , double *w , double* mu, int *info , int *iparamLCP , double *dparamLCP);

  /** NLGS-Newton pfc3D solver */
  void pfc_3D_nlgsnewton(int, double *vec , double *q , double *z , double *w , double *mu, int *info , int *iparamLCP , double *dparamLCP);

  /** NSGS-Newton pfc3D solver */
  void pfc_3D_nsgs(int, double *vec , double *q , double *z , double *w , double *mu,  int *info , int *iparamLCP , double *dparamLCP);

  /** NSGS-Newton pfc3D solver new (temporary) version */
  void pfc_3D_nsgs_new(int, double *M , double *q , double *z , double *w , double *mu, int *info , int *iparamLCP , double *dparamLCP);

  /** Newton pfc3D solver */
  void pfc_3D_newton(int n , double *C , double *b , double *zz , double *ww , double mu , pfc3D_fPtr* Compute_G, pfc3D_fPtr* Compute_JacG, double *param1, double *param2, double *param3, int *iparam_local , double *dparam_local);

  /** projection-type pfc3D solver */
  void pfc_3D_projection(int n , double *C , double *b , double *zz , double *ww , double mu , pfc3D_fPtr* Compute_G, pfc3D_fPtr* Compute_JacG, double *param1, double *param2, double *param3, int *iparam_local , double *dparam_local);

  /** NSGS pfc3D solver, sparse block storage for M. */
  void pfc_3D_nsgs_block(int, SparseBlockStructuredMatrix*, double*, double *, double*, double*, int*, int*, double*);

  /** NSGS pfc3D solver, sparse block storage for M. */
  void pfc_3D_nsgs_block2(int, SparseBlockStructuredMatrix*, double*, double *, double*, double*, int*, int*, double*);

  /** Function to check convergence after PFC computation */
  int filter_result_pfc_block(int, int, double*, SparseBlockStructuredMatrix *, double* , double*, pfc3D_fPtr*, double*, int*, double*);

  /** Function to check convergence after PFC computation */
  int filter_result_pfc_block2(int, int, double*, SparseBlockStructuredMatrix *, double* , double*, pfc3D_fPtr*, double*, int*, double*, double*);

  /** Function to check convergence after PFC computation */
  int filter_result_pfc(int, int, double*, double*, double* , double*, pfc3D_fPtr*, double*, int*, double*, double**);

  /** Compute error function */
  void NCP_compute_error(int, double*, double*, double*, int, double*, double*);

  /** Compute error function */
  void NCP_block_compute_error(int, SparseBlockStructuredMatrix*, double*, double*, int, double*, double*);

#ifdef __cplusplus
}
#endif

#endif
