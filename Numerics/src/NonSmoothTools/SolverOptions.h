/* Siconos-Numerics, Copyright INRIA 2005-2011.
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



#ifndef SolverOptions_H
#define SolverOptions_H

/*!\file SolverOptions.h
  \brief Structure used to send options (name, parameters and so on) to a specific solver-driver (mainly from Kernel to Numerics).
  \author Franck Perignon
*/

/*! \page NumericsSolver Solvers definition in Numerics

To define a non-smooth problem in Numerics, the structure SolverOptions is used. It handles the name of the solver and its input-output parameters.\n
SolverOptions main components are:
 - a name
 - two lists of input-output parameters (int: iparam, double: dparam) and their sizes

Check each type of formulation of the problem to find which solvers are available and what are the required parameters. \n
See for example:
 - \ref LCPSolvers
 - \ref FC3DSolvers

As an example, consider \ref LCProblem : \n
M is a NumericsMatrix and can be saved as a double* or as a SparseBlockStructuredMatrix.\n
One needs to define a SolverOptions, say "options", by choosing one solver among those given in \ref LCPSolvers and set:
\code
int nbSolvers = 1;
SolverOptions options;
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
int info = lcp_driver(myProblem, z,w, &options, nbSolvers, &global_options);
\endcode
which will result in the resolution of the LCP defined in myProblem thanks to a PGS solver.

On the other side if M is saved as a SparseBlockStructuredMatrix, with N rows of blocks, one needs to used a \n
"block-solver" with possibly one or more specific local solver dedicated to each local problem.\n
In that case options must be a vector of SolverOptions, with:\n
 - options[0] the definition for the global "block" solver
 - options[i], i>0, the solvers used for each local problem.


Example with a LCP:
\code
// First define a vector of options
int nbSolvers = 3;
SolverOptions options[nbSolvers];

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
int info = lcp_driver(myProblem, z,w, options,nbSolvers, &global_options);
\endcode

In this case, if the matrix M has N rows of blocks, the global problem will be solved thanks to the Gauss-Seidel block solver, \n
with the first local problem (first row) solved thanks to a PGS and the others with a Lemke. \n
Note that options[i+1] is used for row i of M, while i<nbSolvers-1 and options[nbSolvers-1] for row i when i>=nbSolvers.


*/

#include "NumericsOptions.h"

/** \struct  SolverOptions SolverOptions.h
 Structure used to send options (name, parameters and so on) to a specific solver-driver (mainly from Kernel to Numerics).
    \param isSet int equal to false(0) if the parameters below have not been set (ie need to read default values) else true(1)
    \param solverId Id of the solver (see )
    \param numberOfSolvers : the number of internal or local solvers used by the solver
    \param iSize size of vectors iparam \n
    \param iparam a list of int parameters (depends on each solver, see solver doc.)
    \param dSize size of vector dparam \n
    \param dparam a list of double parameters (depends on each solver, see solver doc.)
    \param filterOn 1 to check solution validity after the driver call, else 0. Default = 1. (For example if \n
    filterOn = 1 for a LCP, lcp_compute_error() will be called at the end of the process)
    \param dWork is a pointer on a working memory zone (for doubles) reserved for the solver .
    \param iWork is a pointer on a working memory zone (for integers) reserved for the solver .
*/
typedef struct _SolverOptions
{
  int solverId;
  int isSet;
  /*char solverName[64] ;*/
  int iSize;
  int * iparam;
  int dSize;
  double * dparam;
  int filterOn;
  double * dWork;
  int * iWork;
  int numberOfInternalSolvers;
  struct _SolverOptions * internalSolvers;
} SolverOptions;

enum SICONOS_NUMERICS_PROBLEM_TYPE
{
  SICONOS_NUMERICS_PROBLEM_LCP = 0,
  SICONOS_NUMERICS_PROBLEM_MLCP = 1,
  SICONOS_NUMERICS_PROBLEM_EQUALITY = 2,
  SICONOS_NUMERICS_PROBLEM_FC2D = 3,
  SICONOS_NUMERICS_PROBLEM_FC3D = 4,
  SICONOS_NUMERICS_PROBLEM_NCP = 5,
  SICONOS_NUMERICS_PROBLEM_MCP = 6
};

extern char * SICONOS_NUMERICS_PROBLEM_LCP_STR;
extern char * SICONOS_NUMERICS_PROBLEM_MLCP_STR;
extern char * SICONOS_NUMERICS_PROBLEM_NCP_STR;
extern char * SICONOS_NUMERICS_PROBLEM_MCP_STR;
extern char * SICONOS_NUMERICS_PROBLEM_EQUALITY_STR;
extern char * SICONOS_NUMERICS_PROBLEM_FC2D_STR;
extern char * SICONOS_NUMERICS_PROBLEM_FC3D_STR;

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Read default parameters values for a solver and save them in a SolverOptions structure
      \param[in] int, type of the considered problem.
      Only the following solvers ids are allowed :\n
      0: LCP\n
      1: MLCP\n
      2: FrictionContact2D\n
      3: FrictionContact3D\n
      \param[out] options structure used to save the parameters
   */
  void readSolverOptions(int driverType, SolverOptions* options);

  /** screen display of solver parameters
      \param options the structure to be displayed
  */
  void printSolverOptions(SolverOptions* options);

  /** delete the solver parameters :
      delete iparam and dparam;
      \param options the structure to be destroyed
  */
  void deleteSolverOptions(SolverOptions * options);

  /* Free the working memory (options->iWork and options->dWork)
     \param[in] options structure used to define the solver(s) and their parameters
  */
  void free_working_memory(SolverOptions* options);

  int nameToId(char * pName);
  char * idToName(int Id);
  char * idProblemToChar(int id);
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif



#endif
