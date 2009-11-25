/* Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
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
/*! \page NSDrivers  Non-Smooth Solvers

Numerics package proposes a set of non-smooth solvers dedicated to some specific formulations for a non-smooth problem.\n

For each type of problem, a generic interface function (driver) is provided. \n
The main arguments of this driver are:
 - a structure of type XXX_Problem (XXX being the formulation type: LinearComplementarity_Problem, FrictionContact_Problem ...), which holds the vectors and matrices used to formalize the problem, \n
 - the unknowns \n
 - a structure of type Solver_Options, used to defined the solver type and its parameters (see \ref NumericsSolver). \n
 - a Numerics_Options structure, used to define global options (verbose mode ...)

To get more details on each formulation, check for each type of problem in \ref NSSpackContents below.

All the drivers interfaces are defined in the file NonSmoothDrivers.h .\n
Each type of formulation is defined in a structure in the file XXX_Problem.h .

\section NSSpackContents Contents
\subpage LCProblem \n\n
\subpage MLCProblem \n\n
\subpage NCProblem \n\n
\subpage fcProblem \n\n
\subpage RelayProblem \n\n
\subpage QPSolvers\n\n

Moreover, details on matrix storage in Numerics can be found in:
\subpage NumericsMatrixPage \n

Other functions and useful tools related to NonSmoothSolvers are listed in NSSTools.h.

*/

/*!\file NonSmoothDrivers.h
  \brief This file provides all generic functions (drivers), interfaces to the different formulations for Non-Smooth Problems available in Numerics.
  \author Nineb Sheherazade and Dubois Frederic.
  Last Modifications : Mathieu Renouf , Pascal Denoyelle, Franck Perignon
  \todo solve_qp does not exist

  Note: all pfc_3D functions, structures ... are out of date and will be soon removed.
  Use FrictionContact3D tools.
*/
#ifndef NonSmoothSolvers_H
#define NonSmoothSolvers_H

#include "Relay_Solvers.h"
#include "LCP_Solvers.h"
#include "MLCP_Solvers.h"
#include "LinearSystem_Problem.h"
#include "FrictionContact2D_Solvers.h"
#include "FrictionContact3D_Solvers.h"
#include "PrimalFrictionContact3D_Solvers.h"
#include "pfc_3D_Solvers.h"
#include "NonSmoothNewton.h"


/** Union of specific methods (one for each type of problem)
    Deprecated.
*/
typedef union
{
  method_pfc_3D pfc_3D;
} method;

#ifdef __cplusplus
extern "C" {
#endif

  /** General interface to solvers for Linear Complementarity Problems
      \param[in] problem the LinearComplementarity_Problem structure which handles the problem (M,q)
      \param[in,out] z a n-vector of doubles which contains the solution of the problem.
      \param[in,out] w a n-vector of doubles which contains the solution of the problem.
      \param[in,out] options structure used to define the solver(s) and their parameters
      \param[in] numberOfSolvers size of the vector options, ie number of solvers \n
      (if numberOfSolvers>1, options[0] represent the global (block) solver and options[i>0] the local ones).
      \param[in] general options for Numerics (verbose mode ...)
      \return info termination value
      - 0 : successful\n
      - >0 : otherwise see each solver for more information about the log info
      \author Franck Perignon
  */
  int lcp_driver(LinearComplementarity_Problem* problem, double *z , double *w, Solver_Options* options, int numberOfSolvers, Numerics_Options* global_options);

  /** General interface to solver for MLCP problems
      \param[in] problem the MixedLinearComplementarity_Problem structure which handles the problem (M,q)
      \param[in,out] z a m+n-vector of doubles which contains the solution of the problem.
      \param[in,out] w a m+n-vector of doubles which contains the solution of the problem.
      \param[in,out] options structure used to define the solver(s) and their parameters
      \param[in] general options for Numerics (verbose mode ...)
      \return info termination value
      - 0 : successful\n
      - >0 : otherwise see each solver for more information about the log info
      \todo Sizing the regularization parameter and apply it only on null diagnal term
      \author Vincent Acary
  */
  int mlcp_driver(MixedLinearComplementarity_Problem* problem, double *z, double *w, Solver_Options* options, Numerics_Options* global_options);
  /** General interface to solver for linear system
      \param[in] problem the LinearSystem_Problem structure which handles the problem (M,q)
      \param[in,out] z a n-vector of doubles solution of the problem.
      \param[out] w a n-vector of doubles which contains zeros.
      \param[in,out] options structure used to define the solver(s) and their parameters
      \param[in] general options for Numerics (verbose mode ...)
      \return info termination value
      - 0 : successful\n
      - >0 : otherwise see each solver for more information about the log info
      \author Vincent Acary
  */
  int LinearSystem_driver(LinearSystem_Problem* problem, double *z , double *w, Solver_Options* options);

  /** General interface to solver for pfc 3D problems - Out of date (will be removed soon) */
  int pfc_3D_driver(int, double*, double*, method*, double*, double*, double*);

  /** General interface to solvers for primal friction 3D problem, with a sparse block storage for M - Out of date (will be removed soon)
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
  int pfc_3D_driver_block(int, SparseBlockStructuredMatrix*, double*, method*, double*, double*, double*);

  /** General interface to solvers for primal friction-contact 2D problem
      \param[in] , the structure which handles the Friction-Contact problem
      \param[in-out] , reaction global vector (n)
      \param[in-out] , velocity global vector (n)
      \param[in,out] options structure used to define the solver(s) and their parameters
      \param[in] general options for Numerics (verbose mode ...)
      \return result (0 if successful otherwise 1).
      \author Nineb Sheherazade.
  */
  int pfc_2D_driver(FrictionContact_Problem* problem, double *reaction , double *velocity, Solver_Options* options, Numerics_Options* global_options);

  /** General interface to solvers for friction-contact 2D problem
      \param[in] , the structure which handles the Friction-Contact problem
      \param[in-out] , reaction global vector (n)
      \param[in-out] , velocity global vector (n)
      \param[in,out] options structure used to define the solver(s) and their parameters
      \param[in] general options for Numerics (verbose mode ...)
      \return result (0 if successful otherwise 1).
  */
  int frictionContact2D_driver(FrictionContact_Problem* problem, double *reaction , double *velocity, Solver_Options* options, Numerics_Options* global_options);

  /** General interface to solver for dual friction-contact problems
      \param[in] , the structure which handles the Friction-Contact problem
      \param[in-out] , reaction global vector (n)
      \param[in-out] , velocity global vector (n)
      \param[in,out] options structure used to define the solver(s) and their parameters
      \param[in] general options for Numerics (verbose mode ...)
      \param[in,out] iparamDFC vector of dfc specific parameters (int):
         - [0] ddl_n:    contact in normal direction dof (not prescribed) (on enter)
         - [1] ddl_tt:   contact in tangential direction dof (not prescribed) (on enter)
         - [2] ddl_d:    prescribed dof (on enter)
         - [3] dim_tt:    dimension of ddl_tt (= dimension of ddl_n) (on enter)
         - [4] dim_d:     dimension of ddl_d (on enter)
      \param[in] J1 gap in normal contact direction
      \return result (0 if successful otherwise 1).
      This problem can be solved thanks to dfc_2D solvers or thanks to lcp solvers after:\n
      - either a condensation makes thanks to dfc_2D2cond_2D.c and cond_2D2dfc_2D.c,
      - or a new formulation of this problem in the LCP form due to the dfc_2D2lcp.c and lcp2dfc_2D.c routines.
      \author Nineb Sheherazade
  */
  int dfc_2D_driver(FrictionContact_Problem* problem, double *reaction , double *velocity, Solver_Options* options, Numerics_Options* global_options, int* iparamDFC, double* J1);

  /** General interface to solver for dual-relay problems
      \param[in] problem the Relay_Problem structure which handles the problem (M,q)
      \param[in,out] z a n-vector of doubles which contains the solution of the problem.
      \param[in,out] w a n-vector of doubles which contains the solution of the problem.
      \param[in,out] options structure used to define the solver(s) and their parameters
      \return info termination value
      - 0 : successful\n
      - >0 : otherwise see each solver for more information about the log info
   * \author Nineb Sheherazade.
   */
  int dr_driver(Relay_Problem* problem, double *z , double *w, Solver_Options* options, Numerics_Options* global_options);

  /** General interface to solver for primal-relay problems
      \param[in] problem the Relay_Problem structure which handles the problem (M,q)
      \param[in,out] z a n-vector of doubles which contains the solution of the problem.
      \param[in,out] w a n-vector of doubles which contains the solution of the problem.
      \param[in,out] options structure used to define the solver(s) and their parameters
      \return info termination value
      - 0 : successful\n
      - >0 : otherwise see each solver for more information about the log info
   * \author Nineb Sheherazade.
   */
  int relay_driver(Relay_Problem* problem, double *z , double *w, Solver_Options* options, int numberOfSolvers,  Numerics_Options* global_options);

  /** General interface to solvers for friction-contact 3D problem
      \param[in] , the structure which handles the Friction-Contact problem
      \param[in-out] , reaction global vector (n)
      \param[in-out] , velocity global vector (n)
      \param[in,out] options structure used to define the solver(s) and their parameters
      \param[in] numberOfSolvers size of the vector options, ie number of solvers \n
      \param[in] general options for Numerics (verbose mode ...)
      \return result (0 if successful otherwise 1).
  */
  int frictionContact3D_driver(FrictionContact_Problem* problem, double *reaction , double *velocity, Solver_Options* options, Numerics_Options* global_options);

  int primalFrictionContact3D_driver(PrimalFrictionContact_Problem* problem, double *reaction , double *velocity, double* globalVelocity, Solver_Options* options, Numerics_Options* global_options);

#ifdef __cplusplus
}
#endif

#endif
