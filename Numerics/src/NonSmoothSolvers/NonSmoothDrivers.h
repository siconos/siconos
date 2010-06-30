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
/*! \page NSDrivers  Non-Smooth Solvers

Numerics package proposes a set of non-smooth solvers dedicated to some specific formulations for a non-smooth problem.\n

For each type of problem, a generic interface function (driver) is provided. \n
The main arguments of this driver are:
 - a structure of type XXXProblem (XXX being the formulation type: LinearComplementarityProblem, FrictionContactProblem ...), which holds the vectors and matrices used to formalize the problem, \n
 - the unknowns \n
 - a structure of type SolverOptions, used to defined the solver type and its parameters (see \ref NumericsSolver). \n
 - a NumericsOptions structure, used to define global options (verbose mode ...)

To get more details on each formulation, check for each type of problem in \ref NSSpackContents below.

All the drivers interfaces are defined in the file NonSmoothDrivers.h .\n
Each type of formulation is defined in a structure in the file XXXProblem.h .

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

  Use FrictionContact3D tools.
*/
#ifndef NonSmoothSolvers_H
#define NonSmoothSolvers_H
#include "mlcp_cst.h"
#include "lcp_cst.h"
#include "Relay_Solvers.h"
#include "LCP_Solvers.h"
#include "MLCP_Solvers.h"
#include "LinearSystemProblem.h"
#include "FrictionContact2D_Solvers.h"
#include "FrictionContact3D_Solvers.h"
#include "PrimalFrictionContact3D_Solvers.h"
#include "GenericMechanical_Solvers.h"

#include "NonSmoothNewton.h"

/** Union of specific methods (one for each type of problem)
    Deprecated. 98
*/

#ifdef __cplusplus
extern "C"
{
#endif


  /** General interface to solver for MLCP problems
      \param[in] problem the MixedLinearComplementarityProblem structure which handles the problem (M,q)
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
  int mlcp_driver(MixedLinearComplementarityProblem* problem, double *z, double *w, SolverOptions* options, NumericsOptions* global_options);
  /** General interface to solver for linear system
      \param[in] problem the LinearSystemProblem structure which handles the problem (M,q)
      \param[in,out] z a n-vector of doubles solution of the problem.
      \param[out] w a n-vector of doubles which contains zeros.
      \param[in,out] options structure used to define the solver(s) and their parameters
      \param[in] general options for Numerics (verbose mode ...)
      \return info termination value
      - 0 : successful\n
      - >0 : otherwise see each solver for more information about the log info
      \author Vincent Acary
  */
  int LinearSystem_driver(LinearSystemProblem* problem, double *z , double *w, SolverOptions* options);

  /** General interface to solvers for friction-contact 2D problem
      \param[in] , the structure which handles the Friction-Contact problem
      \param[in-out] , reaction global vector (n)
      \param[in-out] , velocity global vector (n)
      \param[in,out] options structure used to define the solver(s) and their parameters
      \param[in] general options for Numerics (verbose mode ...)
      \return result (0 if successful otherwise 1).
  */
  int frictionContact2D_driver(FrictionContactProblem* problem, double *reaction , double *velocity, SolverOptions* options, NumericsOptions* global_options);






#ifdef __cplusplus
}
#endif

#endif
