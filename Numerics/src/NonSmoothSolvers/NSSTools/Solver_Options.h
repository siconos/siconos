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

#ifndef Solver_Options_H
#define Solver_Options_H

/*!\file Solver_Options.h
  \brief Structure used to send options (name, parameters and so on) to a specific solver-driver (mainly from Kernel to Numerics).
  \author Franck Perignon
*/

/*! \page NumericsSolver Solvers definition in Numerics

To define a non-smooth problem in Numerics, the structure Solver_Options is used. It handles the name of the solver and its input-output parameters.\n
Solver_Options main components are:
 - a name
 - two lists of input-output parameters (int: iparam, double: dparam) and their sizes

Check each type of formulation of the problem to find which solvers are available and what are the required parameters. \n
See for example:
 - \ref LCPSolvers
 - \ref FrictionContact3DSolvers

As an example, consider \ref LCProblem : \n
M is a NumericsMatrix and can be saved as a double* or as a SparseBlockStructuredMatrix.\n
One needs to define a Solver_Options, say "options", by choosing one solver among those given in \ref LCPSolvers and set:
\code
Solver_Options options;
strcpy(options.solverName,"PGS");
int iparam[2] ={maxIter, 0};
double dparam[2] = {tolerance,0.0};
options.iSize = 2;
options.dSize = 2;
options.iparam = iparam;
options.dparam = dparam;
options.isSet = 1;
\endcode
And then call the driver:
\code
int info = lcp_driver(myProblem, z,w, &options, &global_options);
\endcode
which will result in the resolution of the LCP defined in myProblem thanks to a PGS solver.

On the other side if M is saved as a SparseBlockStructuredMatrix, with N rows of blocks, one needs to used a \n
"block-solver" with possibly one or more specific local solver dedicated to each local problem.\n
In that case options must be a vector of Solver_Options, with:\n
 - options[0] the definition for the global "block" solver
 - options[i], i>0, the solvers used for each local problem.

\bf Example with a LCP:
\code
// First define a vector of options
Solver_Options options[3];

// The global solver:
strcpy(options[0].solverName,"GaussSeidel_SBM");
int iparam[2] ={maxIter, 0};
double dparam[2] = {tolerance,0.0};
options[0].iSize = 2;
options[0].dSize = 2;
options[0].iparam = iparam;
options[0].dparam = dparam;
options[0].isSet = 1;

// The local solvers:
strcpy(options[1].solverName,"PGS");
int iparam[2] ={maxIter, 0};
double dparam[2] = {tolerance,0.0};
options[1].iSize = 2;
options[1].dSize = 2;
options[1].iparam = iparam;
options[1].dparam = dparam;
options[1].isSet = 1;
strcpy(options[2].solverName,"Lemke");
int iparam[2] ={maxIter,0};
double dparam[2] = {tolerance,0.0};
options[2].iSize = 2;
options[2].dSize = 2;
options[2].iparam = iparam;
options[2].dparam = dparam;
options[2].isSet = 1;
\endcode
The call of the driver remains the same:
\code
int info = lcp_driver(myProblem, z,w, options, &global_options);
\endcode

In this case, if the matrix M has N rows of blocks, the global problem will be solved thanks to the Gauss-Seidel block solver, \n
with the first local problem (first row) solved thanks to a PGS and the others with a Lemke. \n
Note that options[i+1] is used for row i of M, while i<nbSolvers-1 and options[nbSolvers-1] for row i when i>=nbSolvers; nbSolvers being the size of options.


*/

#include "Numerics_Options.h"

/** Structure used to send options (name, parameters and so on) to a specific solver-driver (mainly from Kernel to Numerics).
    \param isSet, int equal to false(0) if the parameters below have not been set (ie need to read default values) else true(1)
    \param solverName name of the solver
    \param iSize size of vectors iparam \n
    \param iparam a list of int parameters (depends on each solver, see solver doc.)
    \param dSize size of vector dparam \n
    \param dparam a list of double parameters (depends on each solver, see solver doc.)
    \param storageType int to check storage type (0: double*, 1: SparseBlockStructuredMatrix)
    \param filterOn 1 to check solution validity after the driver call, else 0. Default = 1. (For example if \n
    filterOn = 1 for a LCP, lcp_compute_error() will be called at the end of the process)
    \param localSolver a list of local solvers, used only when the problem is written with SparseBlockStructuredMatrix .
*/
typedef struct
{
  int isSet;
  char solverName[64];
  int iSize;
  int * iparam;
  int dSize;
  double * dparam;
  int filterOn;
} Solver_Options;

#ifdef __cplusplus
extern "C" {
#endif

  /** Read default parameters values for a solver and save them in a Solver_Options structure
      \param[in] driverName, type of the considered problem \n
      0: LCP\n
      1: MLCP\n
      2: FrictionContact2D\n
      3: FrictionContact3D\n
      4: Relay\n
      5: QP\n
      6: NCP
      \param[out] options, structure used to save the parameters
   */
  void readSolverOptions(int, Solver_Options*);

  /** screen display of solver parameters
      \param options, the structure to be diplayed
  */
  void printSolverOptions(Solver_Options*);

#ifdef __cplusplus
}
#endif



#endif
