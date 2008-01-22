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
/*! \page NSSPACK Non-Smooth Solvers Pack

Numerics package proposes a set of non-smooth solvers dedicated to some specific formulations.\n

For each type of problem, a generic interface function is provided (see NSSpack.h for the complete list).\n
The main arguments of this function are the vectors and matrices used to formalize the problem, \n
the unknowns and a structure method, used to defined the solver type and its parameters. \n
To get more details on this, check for each type of problem in \ref NSSpackContents below.

\section NSSpackContents Contents
\subpage LCPSolvers \n
\subpage MLCPSolvers \n
\subpage fc2DSolvers \n
\subpage fc3DSolvers \n
\subpage RelaySolvers\n
\subpage QPSolvers\n

Other functions and useful tools related to NSSpack are listed in NSSTools.h.

*/

/*!\file NSSpack.h
 * \brief This file provides all generic functions, interfaces to the different solvers available in Numerics.
 *  \author Nineb Sheherazade and Dubois Frederic.
 *  Last Modifications : Mathieu Renouf , Pascal Denoyelle, Franck Perignon
 *  \todo solve_qp does not exist
 *
 */
#ifndef NSSPACK_H
#define NSSPACK_H

#include "blaslapack.h"
#include "Relay_Solvers.h"
#include "LCP_Solvers.h"
#include "MLCP_Solvers.h"
#include "FrictionContact2D_Solvers.h"
#include "FrictionContact3D_Solvers.h"
#include "pfc_3D_Solvers.h"
#include "NonSmoothNewton.h"


/** Union of specific methods (one for each type of problem)
    \param method_pr     : pr is a method_pr structure .
    \param method_dr     : dr is a method_dr structure .
    \param method_lcp    : lcp is a method_lcp structure .
    \param method_mlcp    : mlcp is a method_mlcp structure .
    \param method_pfc_2D : pfc_2D is a method_pfc_2D structure .
    \param method_dfc_2D : dfc_2D is a method_dfc_2D structure .
*/
typedef union
{
  method_pr  pr;
  method_dr  dr;
  method_lcp lcp;
  method_mlcp mlcp;
  method_pfc_2D pfc_2D;
  method_pfc_3D pfc_3D;
  method_dfc_2D dfc_2D;
  method_frictionContact_3D fc3D_parameters;
} method;

#ifdef __cplusplus
extern "C" {
#endif

  /** General interface to solver for LCP problems
      \param[in] , a (\f$n \times n\f$)-vector of doubles which contains the components of the LCP matrix with a Fortran storage.
      \param[in] , a n-vector of doubles which contains the components of the constant right hand side vector.
      \param[in] , an integer which represents the dimension of the LCP problem.
      \param[in] , a union containing the LCP structure.
      \param[in-out] , a n-vector of doubles which contains the solution of the problem.
      \param[in-out] , a n-vector of doubles which contains the complementary solution of the problem.
      \return integer
      - 0 : successful\n
      - >0 : otherwise see each solver for more information about the log info
      \author Nineb Sheherazade & Mathieu Renouf
  */
  int lcp_solver(double*, double*, int*, method*, double*, double*);

  /** General interface to solver for LCP problems, with M given as a list of blocks */
  int lcp_solver_block(SparseBlockStructuredMatrix *, double *, method *, double *, double *, int *, int *, double*);

  /** General interface to solver for MLCP problems
      \param[in] , a (\f$n \times n\f$)-vector of doubles which contains the components of the "A" MLCP matrix with a Fortran storage.
      \param[in] , a (\f$m \times m\f$)-vector of doubles which contains the components of the "B" MLCP matrix with a Fortran storage.
      \param[in] , a (\f$n \times m\f$)-vector of doubles which contains the components of the "C" MLCP matrix with a Fortran storage.
      \param[in] , a (\f$m \times n\f$)-vector of doubles which contains the components of the "D" MLCP matrix with a Fortran storage.
      \param[in] , a n-vector of doubles which contains the components of the constant right hand side vector.
      \param[in] , a m-vector of doubles which contains the components of the constant right hand side vector.
      \param[in] , an integer which represents one the dimension of the MLCP problem.
      \param[in] , an integer which represents one the dimension of the MLCP problem.
      \param[in] , a union containing the MLCP structure.
      \n \n
      \param[out] , a n-vector of doubles which contains the solution of the problem.
      \param[out] , a m-vector of doubles which contains the solution of the problem.
      \param[out] , a m-vector of doubles which contains the complementary solution of the problem.
      \return integer
      - 0 : successful\n
      - >0 : otherwise (see specific solvers for more information about the log info)
      \author V. Acary
      \todo Sizing the regularization paramter and apply it only on null diagnal term
  */
  int mlcp_solver(double*, double*, double*, double*, double*, double*, int*, int*, method*, double*, double*, double*);

  /** General interface to solver for pfc 3D problems */
  int pfc_3D_solver(int, double*, double*, method*, double*, double*, double*);

  /** General interface to solvers for primal friction 3D problem, with a sparse block storage for M.
   *   \param n    Dimension of the system. Then, the number of contacts is n/3.
   *  \param M    components of the double matrix with a fortran allocation.
   *  \param q    the components of the second member of the system.
   *  \param pt   structure
   *  \param z    the solution of the problem.
   *  \param w    the complementarity solution of the problem.
   *  \param mu   the list of friction coefficients. mu[i] corresponds to contact number i.
   *  \param iter On return, an integer, the final number of block iterations
   *  \param local_iter On return, an integer, the total number of local iterations
   *  \param error   On return, a double, the final value of error criteria.
   *
   *  \return     result (0 is successful otherwise 1).
   *
   *  \author Franck Perignon.
   *
   * The available solvers are (with the corresponding keywords):
   *
   *   - NSGS for non-smooth Gauss Seidel
   *
   */
  int pfc_3D_solver_block(int, SparseBlockStructuredMatrix*, double*, method*, double*, double*, double*);

  /** General interface to solver for pfc 2D problems
   *  \param[in] , a (n \f$\times\f$n)-vector of doubles which contains the components of the double matrix with a fortran allocation.
   *  \param[in] , a n-vector of doubles containing the components of the second member of the system.
   *  \param[in] , an integer, the dimension of the second member.
   *  \param[in] , a union (::method) containing the PFC_2D structure.
   *  \param[out] , a n-vector of doubles containing the solution of the problem.
   *  \param[out] , a n-vector of doubles containing the solution of the problem.
   *  \return     integer
   *                       - 0: successful,
   *                       - otherwise (see specific solvers for more information about the
   *                           termination reason).
   * \author Nineb Sheherazade.
   */
  int pfc_2D_solver(int, double*, double*, method*, double*, double*, double*);

  /** General interface to solver for dual friction-contact problems
      \param[in] , the stiffness, a vector of double (in which the components
      of the matrix have a Fortran storage),
      \param[in] , the right hand side, a vector of double,
      \param[in] , the dimension of the DFC_2D problem, an integer,
      \param[in] , the union (::method) containing the DFC_2D structure, with the following parameters:\n
      - char   name:      the name of the solver we want to use (on enter),
      - int    itermax:   the maximum number of iteration required (on enter)
      - double tol:       the tolerance required (on enter)
      - double  mu:       the friction coefficient (on enter)
      - int    *ddl_n:    contact in normal direction dof (not prescribed) (on enter)
      - int    *ddl_tt:   contact in tangential direction dof (not prescribed) (on enter)
      - int    *ddl_d:    prescribed dof (on enter)
      - int    dim_tt:    dimension of ddl_tt (= dimension of ddl_n) (on enter)
      - int    dim_d:     dimension of ddl_d (on enter)
      - double *J1:       gap in normal contact direction (on enter)
      - int    chat:      an integer that can make on or off the chattering (0=off, >0=on)(on enter)
      - int    iter:      the number of iteration carry out (on return)
      - double err:       the residue (on return) \n\n
      This problem can be solved thanks to dfc_2D solvers or thanks to lcp solvers after:\n
      - either a condensation makes thanks to dfc_2D2cond_2D.c and cond_2D2dfc_2D.c,
      - or a new formulation of this problem in the LCP form due to the dfc_2D2lcp.c and lcp2dfc_2D.c routines.
      \param[out] , the solution of the problem, vector of double.
      \param[out] , the solution of the problem, vector of double.
      \return integer     0 - successful\n
      0 >  - otherwise (see specific solvers for more information about the log info)
      \author Nineb Sheherazade
  */
  int dfc_2D_solver(double* , double* , int* , method* , double* , double*, double*);

  /** General interface to solver for dual-relay problems
   * \param[in] , (nn \f$\times\f$nn)-vector of doubles containing the components of the double matrix with a fortran90 allocation.
   * \param[in] , a nn-vector of doubles containing the components of the second member of the system.
   * \param[in] , an integer, the dimension of the second member.
   * \param[in] , a union (::method) containing the DR structure.\n \n
   * \param[out] , a nn-vector of doubles, the solution of the problem.
   * \param[out] , a nn-vector of doubles, the solution of the problem.
   * \author Nineb Sheherazade.
   */
  int dr_solver(double* , double* , int* , method* , double* , double*);

  /** General interface to solver for primal-relay problems
   *  \param[in] , a (\f$nn \times nn\f$)-vector of doubles containing the components of the  matrix with a fortran90 allocation.
   *  \param[in] , a nn-vector of doubles containing the components of the second member of the system.
   *  \param[in] , an integer which represents the dimension of the second member.
   *  \param[in] , a union (::method) containing the PR structure.
   *  \param[out] , a nn-vector of doubles which contains the solution of the problem.
   *  \param[out] , a nn-vector of doubles which contains the solution of the problem.
   *   \return integer : the termination reason\n
   *                    - 0 : successful\n
   *                    - otherwise : see specific solvers for more information about the log info.
   *   \author Nineb Sheherazade.
   */
  int pr_solver(double* , double* , int* , method* , double* , double*);

  /** General interface to solvers for friction-contact 3D problem
      \param[in] , number of contacts (dim of the problem n = 3*nc)
      \param[in] , M global matrix (n*n)
      \param[in] , q global vector (n)
      \param[in] , method
      \param[in-out] , reaction global vector (n)
      \param[in-out] , velocity global vector (n)
      \param[in] , mu vector of the friction coefficients (size nc)
      \return result (0 if successful otherwise 1).
  */
  int frictionContact3D_solver(FrictionContact3D_Problem* problem, double *reaction , double *velocity, Solver_Options* options, Numerics_Options* global_options);
  //int frictionContact3D_solver(int, double*, double*, method*, double*, double*, double*);

  /** General interface to solvers for friction-contact 3D problem with sparse-block storage for M
      \param[in] , number of contacts (dim of the problem n = 3*nc)
      \param[in] , M global matrix (n*n)
      \param[in] , q global vector (n)
      \param[in] , method
      \param[in-out] , reaction global vector (n)
      \param[in-out] , velocity global vector (n)
      \param[in] , mu vector of the friction coefficients (size nc)
      \return result (0 if successful otherwise 1).
  */
  int frictionContact3D_solver_SBS(int, SparseBlockStructuredMatrix*, double*, method*, double*, double*, double*);

#ifdef __cplusplus
}
#endif

#endif /* NSSPACK_H */
