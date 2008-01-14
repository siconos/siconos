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

/*! \page pfc3DSolvers Primal friction-contact problems (3-dimensional) - OUT OF DATE

pfc3D directory handles all files related to  Primal friction-contact problems.

It has been (partially) replaced with FrictionContact3D.

Thus all pfc_3D_* files are to be removed when FrictionContact3D will be fully operational.

primal resolution of contact problems with friction.\n
*
* Try \f$(z,w)\f$ such that:\n
*  \f$
*   \left\lbrace
*    \begin{array}{l}
*     M z + q = w\\
*     0 \le z_n \perp w_n \ge 0\\
*     -w_t \in \partial\psi_{[-\mu z_n, \mu z_n]}(z_t)\\
*    \end{array}
*   \right.
*  \f$
*
* here M is an (n x n) matrix, q, z and w n-vector.
*/

/*!\file pfc_3D_Solvers.h
  Typedef and functions declarations related to primal frictional contact problems.
  \author Nineb Sheherazade and Dubois Frederic.
  Last Modifications : Mathieu Renouf , Pascal Denoyelle, Franck Perignon
*/

#include "NSSTools.h"
#include "pfc_3D_Alart_Curnier.h"
#include "pfc_3D_Fischer_Burmeister.h"

/** A type definition for a structure method_pfc_3D.
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

  /** pfc_3D_cpg is a specific cpg (conjugated projected gradient) for primal contact problem with friction.\n
   * Ref: Renouf, M. and Alart, P. "" Comp. Method Appl. Mech. Engrg. (2004).
   *
   * Generic pfc_3D parameters:\n
   * \param nc      Unchanged parameter which represents the number of contatcs. The dimension of the system is n = 3*nc.
   * \param vec     Unchanged parameter which contains the components of the matrix with a fortran storage.
   * \param q       Unchanged parameter which contains the components of the right hand side vector.
   * \param z       Modified parameter which contains the initial solution and returns the solution of the problem.
   * \param w       Modified parameter which returns the solution of the problem.
   * \param mu   the list of friction coefficients. mu[i] corresponds to contact number i.
   * \param info    Modified parameter which returns the termination value\n
   *                0 - convergence\n
   *                1 - iter = itermax\n
   *
   * Specific CPG parameters:\n
   *
   * \param iparamLCP[0] = itermax Input unchanged parameter which represents the maximum number of iterations allowed.
   * \param iparamLCP[1] = ispeak  Input unchanged parameter which represents the output log identifiant\n
   *                       0 - no output\n
   *                       0 < active screen output\n
   * \param iparamLCP[2] = it_end  Output modified parameter which returns the number of iterations performed by the algorithm.
   *
   * \param dparamLCP[0] = tol     Input unchanged parameter which represents the tolerance required.
   * \param dparamLCP[1] = res     Output modified parameter which returns the final error value.
   * \author Nineb Sheherazade & Mathieu Renouf.
   *
   */
  void pfc_3D_cpg(int, double *vec , double *q , double *z , double *w , double *mu, int *info , int *iparamLCP , double *dparamLCP);

  /** Non linear Gauss-Seidel pfc3D solver
   * \param nc      Unchanged parameter which represents the number of contacts. The dimension of the system is then n=3*nc.
   * \param vec     Unchanged parameter which contains the components of the matrix with a fortran storage.
   * \param q       Unchanged parameter which contains the components of the right hand side vector.
   * \param z       Modified parameter which contains the initial solution and returns the solution of the problem.
   * \param w       Modified parameter which returns the solution of the problem.
   * \param mu   the list of friction coefficients. mu[i] corresponds to contact number i.
   * \param info    Modified parameter which returns the termination value\n
   *                0 - convergence\n
   *                1 - iter = itermax\n
   *                2 - negative diagonal term\n
   *
   * Specific NLGS parameters:\n
   *
   * \param iparamLCP[0] = itermax Input unchanged parameter which represents the maximum number of iterations allowed.
   * \param iparamLCP[1] = ispeak  Input unchanged parameter which represents the output log identifiant\n
   *                       0 - no output\n
   *                       0 < active screen output\n
   * \param iparamLCP[2] = it_end  Output modified parameter which returns the number of iterations performed by the algorithm.
   *
   * \param dparamLCP[0] = tol     Input unchanged parameter which represents the tolerance required.
   * \param dparamLCP[1] = res     Output modified parameter which returns the final error value.
   *
   *
   * \author Mathieu Renouf & Nineb Sheherazade (Houari Khenous last modification (08/10/2007)).
   *
   */
  void pfc_3D_nlgs(int, double *vec , double *q , double *z , double *w , double* mu, int *info , int *iparamLCP , double *dparamLCP);

  /** NLGS-Newton pfc3D solver
   * Generic pfc_3D parameters:\n
   *
   * \param nc      Unchanged parameter which represents the number of contacts. The dimension of the system is 3*nc.
   * \param vec     Unchanged parameter which contains the components of the matrix with a fortran storage.
   * \param q       Unchanged parameter which contains the components of the right hand side vector.
   * \param z       Modified parameter which contains the initial solution and returns the solution of the problem.
   * \param w       Modified parameter which returns the solution of the problem.
   * \param mu   the list of friction coefficients. mu[i] corresponds to contact number i.
   * \param info    Modified parameter which returns the termination value\n
   *                0 - convergence\n
   *                1 - iter = itermax\n
   *                2 - negative diagonal term\n
   *
   * Specific NLGS parameters:\n
   *
   * \param iparamLCP[0] = itermax Input unchanged parameter which represents the maximum number of iterations allowed.
   * \param iparamLCP[1] = ispeak  Input unchanged parameter which represents the output log identifiant\n
   *                       0 - no output\n
   *                       0 < active screen output\n
   * \param iparamLCP[2] = it_end  Output modified parameter which returns the number of iterations performed by the algorithm.
   *
   * \param dparamLCP[0] = tol     Input unchanged parameter which represents the tolerance required.
   * \param dparamLCP[1] = res     Output modified parameter which returns the final error value.
   *
   *
   * \author Mathieu Renouf & Nineb Sheherazade .
   *
   */
  void pfc_3D_nlgsnewton(int, double *vec , double *q , double *z , double *w , double *mu, int *info , int *iparamLCP , double *dparamLCP);

  /** NSGS-Newton pfc3D solver
   *
   * Generic pfc_3D parameters:\n
   *
   * \param nn      Unchanged parameter which represents the number of contacts. The dimension of the system is 3*nc.
   * \param vec     Unchanged parameter which contains the components of the matrix with a fortran storage.
   * \param q       Unchanged parameter which contains the components of the right hand side vector.
   * \param z       Modified parameter which contains the initial solution and returns the solution of the problem.
   * \param w       Modified parameter which returns the solution of the problem.
   * \param info    Modified parameter which returns the termination value\n
   *                0 - convergence\n
   *                1 - iter = itermax\n
   *                2 - negative diagonal term\n
   *
   * Specific NLGS parameters:\n
   *
   * \param iparamLCP[0] = itermax Input unchanged parameter which represents the maximum number of iterations allowed.
   * \param iparamLCP[1] = ispeak  Input unchanged parameter which represents the output log identifiant\n
   *                       0 - no output\n
   *                       0 < active screen output\n
   * \param iparamLCP[2] = it_end  Output modified parameter which returns the number of iterations performed by the algorithm.
   * \param iparamLCP[5] = local_formulation Formulaton of frictional contact problem (Alart-Curnier: 0 and Fischer-Burmeister: 1)
   * \param iparamLCP[6] = local_solver      Solution method used to solve the frictional contact problem (projection: 0 and Newton: 1)
   *  [3] and [4] are useless.
   * \param dparamLCP[0] = tol     Input unchanged parameter which represents the tolerance required.
   * \param dparamLCP[1] = res     Output modified parameter which returns the final error value.
   *
   *
   * Note that the problem is solved with the Non-Smooth Gauss-Seidel method
   *
   * Important: The projection method is done only with Alart-Curnier formulation
   *
   * In progress: The Fischer-Burmeister formulation is nore ready to be used and is in progress
   *
   * This formulation is based on the idea of C. Glocker: Formulation of spatial contact situations
   * in rigid multibodiy systems, Comp. Meth. Appl. Mech. Engrg, 177 (1999), 199-214.
   *
   * Last modification 27/11/2007, H. Khenous
   *
   */
  void pfc_3D_nsgs(int, double *vec , double *q , double *z , double *w , double *mu,  int *info , int *iparamLCP , double *dparamLCP);

  /** NSGS-Newton pfc3D solver new (temporary) version
   *
   * Generic pfc_3D parameters:\n
   *
   * \param nc      Unchanged parameter which represents the number of contacts. The dimension of the system is 3*nc.
   * \param M       Unchanged parameter which contains the components of the matrix with a fortran storage (column major).
   * \param q       Unchanged parameter which contains the components of the right hand side vector.
   * \param z       Modified parameter which contains the initial solution and returns the solution of the problem.
   * \param w       Modified parameter which returns the solution of the problem.
   * \param mu   the list of friction coefficients. mu[i] corresponds to contact number i.
   * \param info    Modified parameter which returns the termination value\n
   *                0 - convergence\n
   *                1 - iter = itermax, ie the simulation reached the maximum number of iterations allowed\n
   *                2 - negative diagonal term(s) in M.\n
   *
   * Specific NSGS parameters:\n
   *
   * \param iparam[0] = itermax Input unchanged parameter which represents the maximum number of iterations allowed.
   * \param iparam[1] = ispeak  Input unchanged parameter which represents the output log identifiant\n
   *                       0 - no output\n
   *                       1 - active screen output\n
   * \param iparam[2] = it_end  Output modified parameter which returns the number of iterations performed by the algorithm.
   * \param iparam[3] = local iter_max
   * \param iparam[4] = iter local (output)
   * \param iparam[5] = local formulation (0: Alart-Curnier, 1: Fischer-Burmeister)
   * \param iparam[6] = local solver (0: projection, 1: Newton). Projection only for AC case.
   *
   * \param dparam[0] = tol     Input unchanged parameter which represents the tolerance required.
   * \param dparam[1] = error   Output modified parameter which returns the final error value.
   * \param dparam[2] = local tolerance
   * \param dparam[3] = Output modified parameter which returns the local error
   *
   *
   * \author Houari Khenous and Franck Perignon - Creation: 12/11/2007 - Last modification 14/11/2007.
   *
   *
   */
  void pfc_3D_nsgs_new(int, double *M , double *q , double *z , double *w , double *mu, int *info , int *iparamLCP , double *dparamLCP);

  /** Newton pfc3D solver
   *
   *
   * The method is done as follows:
   * 1) Take x0;
   * 2) Compute G(x{k}), Jac_G(x{k});
   * 3) if || G(x{k}) || < tolerance then stop
   * 4) else Solve Jac_G(x{k}) d = -G(x{k});
   *         Do linesearch to obtain a decent direction D depending on d
   *         x{k+1} = x{k} + D
   *         Compute G(x{k+1})
   *         if || G(x{k+1}) || < tolerance then stop
   *         else go to 2)
   *
   * \author houari khenous last modification 27/11/2007.
   *
   */
  void pfc_3D_newton(int n , double *C , double *b , double *zz , double *ww , double mu , pfc3D_fPtr* Compute_G, pfc3D_fPtr* Compute_JacG, double *param1, double *param2, double *param3, int *iparam_local , double *dparam_local);

  /** projection-type pfc3D solver
   *
   * \param nn      Unchanged parameter which represents the dimension of the system.
   * \param vec     Unchanged parameter which contains the components of the matrix with a fortran storage.
   * \param q       Unchanged parameter which contains the components of the right hand side vector.
   * \param z       Modified parameter which contains the initial solution and returns the solution of the problem.
   * \param w       Modified parameter which returns the solution of the problem.
   * \param coef    Unchanged parameter which represents the friction coefficient
   *
   * \author houari khenous last modification 08/10/2007 .
   *
   */
  void pfc_3D_projection(int n , double *C , double *b , double *zz , double *ww , double mu , pfc3D_fPtr* Compute_G, pfc3D_fPtr* Compute_JacG, double *param1, double *param2, double *param3, int *iparam_local , double *dparam_local);

  /** NSGS pfc3D solver, sparse block storage for M.
   *
   * Generic pfc_3D parameters:\n
   *
   * \param nc      Unchanged parameter which represents the number of contacts. The dimension of the system is 3*nc.
   * \param M       Unchanged parameter which contains the matrix M, saved as a list of non-null blocks (sparse block structure)
   * \param q       Unchanged parameter which contains the components of the right hand side vector.
   * \param z       Modified parameter which contains the initial solution and returns the solution of the problem.
   * \param w       Modified parameter which returns the solution of the problem.
   * \param mu   the list of friction coefficients. mu[i] corresponds to contact number i.
   * \param info    Modified parameter which returns the termination value\n
   *                0 - convergence\n
   *                1 - iter = itermax, ie the simulation reached the maximum number of iterations allowed\n
   *                2 - negative diagonal term(s) in M.\n
   *
   * Specific NSGS parameters:\n
   *
   * \param iparam[0] = itermax Input unchanged parameter which represents the maximum number of iterations allowed.
   * \param iparam[1] = ispeak  Input unchanged parameter which represents the output log identifiant\n
   *                       0 - no output\n
   *                       1 - active screen output\n
   * \param iparam[2] = it_end  Output modified parameter which returns the number of iterations performed by the algorithm.
   * \param iparam[3] = local iter_max
   * \param iparam[4] = iter local (output)
   * \param iparam[5] = local formulation (0: Alart-Curnier, 1: Fischer-Burmeister)
   * \param iparam[6] = local solver (0: projection, 1: Newton). Projection only for AC case.
   *
   * \param dparam[0] = tol     Input unchanged parameter which represents the tolerance required.
   * \param dparam[1] = error   Output modified parameter which returns the final error value.
   * \param dparam[2] = local tolerance
   * \param dparam[3] = Output modified parameter which returns the local error
   *
   * \author Houari Khenous and Franck Perignon - Creation: 12/11/2007 - Last modification 03/12/2007.
   *
   *
   */
  void pfc_3D_nsgs_block(int, SparseBlockStructuredMatrix*, double*, double *, double*, double*, int*, int*, double*);

  /** pfc_3D_projc is a specific projection operator related to CPG (conjugated projected gradient) algorithm
   *              for primal contact problem with friction.\n
   *
   * Ref: Renouf, M. and Alart, P. "" Comp. Method Appl. Mech. Engrg. (2004).
   *
   * \param n       Unchanged parameter which represents the half dimension of the system.
   * \param mu      Unchanged parameter which represents the friction coefficients list
   * \param z       Modified parameter which retruns the corrected iterate.
   * \param p       Unchanged parameter which contains the components of the descent direction.
   * \param status  Unchanged parameter which contains the vector status
   *
   * \author Mathieu Renouf.
   *
   */
  void pfc_3D_projc(int nc , double* mu , double *z , double *p , int *status);

  /** pfc_3D_projf is a specific projection operator related to CPG (conjugated projected gradient) algorithm
   *              for primal contact problem with friction.\n
   *
   * Ref: Renouf, M. and Alart, P. "" Comp. Method Appl. Mech. Engrg. (2004).
   *
   * \param n       Unchanged parameter which represents the half dimension of the system.
   * \param ww      Modified parameter which returns the projected residue.
   * \param zz      Modified parameter which retruns the projected descent direction.
   * \param rr      Unchanged parameter which contains the components of the residue vector.
   * \param pp      Unchanged parameter which contains the components of the descent direction.
   * \param status  Unchanged parameter which contains the vector status
   *
   * \author Mathieu Renouf.
   *
   */
  void pfc_3D_projf(int nc , double *ww , double *zz , double *rr , double *pp , int *status);

  /*   /\** NSGS pfc3D solver, sparse block storage for M. *\/ */
  /*   void pfc_3D_nsgs_block2(int, SparseBlockStructuredMatrix*, double*, double *, double*, double*, int*, int*, double*); */

  /** Function to check convergence after PFC computation */
  int filter_result_pfc_block(int, int, double*, SparseBlockStructuredMatrix *, double* , double*, pfc3D_fPtr*, double*, int*, double*);

  /** Function to check convergence after PFC computation */
  int filter_result_pfc_block2(int, int, double*, SparseBlockStructuredMatrix *, double* , double*, pfc3D_fPtr*, double*, int*, double*, double*);

  /** Function to check convergence after PFC computation */
  int filter_result_pfc(int, int, double*, double*, double* , double*, pfc3D_fPtr*, double*, int*, double*, double**);

#ifdef __cplusplus
}
#endif

#endif
