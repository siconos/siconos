/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

/*!\file NonSmoothDrivers.h
  \brief This file provides all generic functions (drivers), interfaces to the different formulations for Non-Smooth Problems available in Numerics.
  \todo solve_qp does not exist

  Use fc3d tools.
*/
#ifndef NonSmoothSolvers_H
#define NonSmoothSolvers_H

#include "SiconosConfig.h"

/* #include "mlcp_cst.h" */
/* #include "MCP_cst.h" */
/* #include "NCP_cst.h" */
/* #include "lcp_cst.h" */
/* #include "relay_cst.h" */
/* #include "Friction_cst.h" */
/* #include "VI_cst.h" */
/* #include "AVI_cst.h" */
/* #include "SOCLCP_cst.h" */
//#include "VariationalInequality.h"
//#include "VariationalInequality_Solvers.h"
//#include "SecondOrderConeLinearComplementarityProblem.h"
/* #include "SOCLCP_Solvers.h" */
/* #include "Relay_Solvers.h" */
/* #include "LCP_Solvers.h" */
/* #include "AVI_Solvers.h" */
/* #include "MLCP_Solvers.h" */
/* #include "NCP_Solvers.h" */
/* #include "MCP_Solvers.h" */
//#include "SolverOptions.h"
#include "NumericsFwd.h"
//#include "MixedComplementarityProblem.h"
/* #include "fc2d_Solvers.h" */
/* #include "fc3d_Solvers.h" */
/* #include "gfc3d_Solvers.h" */
//#include "GenericMechanical_Solvers.h"

//#include "NonSmoothNewton.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** General interface to solvers for Linear Complementarity Problems
    \param[in] problem the LinearComplementarityProblem structure which handles the problem (M,q)
    \param[in,out] z a n-vector of doubles which contains the solution of the problem.
    \param[in,out] w a n-vector of doubles which contains the solution of the problem.
    \param[in,out] options structure used to define the solver(s) and their parameters
    \return info termination value
    - 0 : successful
    - >0 : otherwise see each solver for more information about the log info
  */
  int linearComplementarity_driver(LinearComplementarityProblem* problem, double *z , double *w, SolverOptions* options);

  /** General interface to solver for MLCP problems
      \param[in] problem the MixedLinearComplementarityProblem structure which handles the problem (M,q)
      \param[in,out] z a m+n-vector of doubles which contains the solution of the problem.
      \param[in,out] w a m+n-vector of doubles which contains the solution of the problem.
      \param[in,out] options structure used to define the solver(s) and their parameters
      \return info termination value
      - 0 : successful
      - >0 : otherwise see each solver for more information about the log info
      \todo Sizing the regularization parameter and apply it only on null diagnal term
  */
  int mlcp_driver(MixedLinearComplementarityProblem* problem, double *z, double *w, SolverOptions* options);

  /** General interface to solvers for friction-contact 2D problem
   *  \param[in] problem the structure which handles the Friction-Contact problem
   *  \param[in,out] reaction global vector (n)
   *  \param[in,out] velocity global vector (n)
   *  \param[in,out] options structure used to define the solver(s) and their parameters
   *  \return result (0 if successful otherwise 1).
   */
  int fc2d_driver(FrictionContactProblem* problem, double *reaction , double *velocity, SolverOptions* options);


  /** General interface to solvers for friction-contact 3D problem
   *  \param[in] problem the structure which handles the Friction-Contact problem
   *  \param[in,out] reaction global vector (n)
   *  \param[in,out] velocity global vector (n)
   *  \param[in,out] options structure used to define the solver(s) and their parameters
   *  \return result (0 if successful otherwise 1).
   */
  int fc3d_driver(FrictionContactProblem* problem, double *reaction , double *velocity, SolverOptions* options);

  /** General interface to solvers for friction-contact 3D problem
   *  \param[in] problem the structure which handles the Friction-Contact problem
   *  \param[in,out] reaction global vector (n)
   *  \param[in,out] velocity global vector (n)
   *  \param[in,out] options structure used to define the solver(s) and their parameters
   *  \return result (0 if successful otherwise 1).
   */
  int rolling_fc3d_driver(RollingFrictionContactProblem* problem, double *reaction , double *velocity, SolverOptions* options);

  /** General interface to solvers for global friction-contact 3D problem
    \param[in] problem the structure which handles the Friction-Contact problem
    \param[in,out] reaction global vector (n)
    \param[in,out] velocity global vector (n)
    \param[in,out] globalVelocity global vector
    \param[in,out] options structure used to define the solver(s) and their parameters
    \return result (0 if successful otherwise 1).
  */
  int gfc3d_driver(GlobalFrictionContactProblem* problem, double *reaction ,
                                     double *velocity, double* globalVelocity,
                                     SolverOptions* options);

  /** General interface to solvers for friction-contact 3D problem
   *  \param[in] problem the structure which handles the Friction-Contact problem
   *  \param[in,out] x global vector (n)
   *  \param[in,out] w global vector (n)
   *  \param[in,out] options structure used to define the solver(s) and their parameters
   *  \return result (0 if successful otherwise 1).
   */
  int variationalInequality_driver(VariationalInequality* problem, double *x , double *w, SolverOptions* options);

  /** General interface to solvers for Affine Variational Inequalities (AVI)
    \param[in] problem the AffineVariationalInequalities structure which handles the problem (M,q)
    \param[in,out] z a n-vector of doubles which contains the solution of the problem.
    \param[in,out] w a n-vector of doubles which contains the solution of the problem.
    \param[in,out] options structure used to define the solver(s) and their parameters
    \return info termination value
    - 0 : successful
    - >0 : otherwise see each solver for more information about the log info
  */
  int avi_driver(AffineVariationalInequalities* problem, double* z, double* w, SolverOptions* options);

  /** General interface to solver for MCP problems
      \param[in] problem the MixedComplementarityProblem structure which handles the problem
      \param[in,out] z a m+n-vector of doubles which contains the solution of the problem.
      \param[in,out] w a m+n-vector of doubles which contains the solution of the problem.
      \param[in,out] options structure used to define the solver(s) and its(their) parameters
      \return info termination value  0 : successful, else error.
  */
  int mcp_driver(MixedComplementarityProblem* problem, double *z, double *w, SolverOptions* options);

  /** General interface to solver for MCP problems -- new version
      \param[in] problem the MixedComplementarityProblem2 structure which handles the problem
      \param[in,out] z a m+n-vector of doubles which contains the solution of the problem.
      \param[in,out] w a m+n-vector of doubles which contains the solution of the problem.
      \param[in,out] options structure used to define the solver(s) and its(their) parameters
      \return info termination value  0 : successful, else error.
  */
  int mcp_driver2(MixedComplementarityProblem2* problem, double *z, double *w, SolverOptions* options);

  /** General interface to solver for NCP problems
      \param[in] problem the NonlinearComplementarityProblem structure which handles the problem
      \param[in,out] z a n-vector of doubles which contains the solution of the problem.
      \param[in,out] F a n-vector of doubles which contains value of the function evaluated at the solution of the problem.
      \param[in,out] options structure used to define the solver(s) and its(their) parameters
      \return info termination value  0 : successful, else error
  */
  int ncp_driver(NonlinearComplementarityProblem* problem, double *z , double *F, SolverOptions* options);

  /** General interface to solvers for SOCLCP problem
      \param[in] problem the structure which handles the Friction-Contact problem
      \param[in,out] r global vector (n)
      \param[in,out] v global vector (n)
      \param[in,out] options structure used to define the solver(s) and their parameters
      \return result (0 if successful otherwise 1).
  */
  int soclcp_driver(SecondOrderConeLinearComplementarityProblem* problem, double *r , double *v, SolverOptions* options);

  /** LMGC interface to solvers for friction-contact 3D problem
   *  \param[in,out] reaction global vector (nc*3)
   *  \param[in,out] velocity global vector (nc*3)
   *  \param[in] q global vector (nc*3)
   *  \param[in] mu global vector (nc)
   *  \param[in] W the block matrix in coordinate format
   *  \param[in] row block row indices
   *  \param[in] column block column indices
   *  \param[in] nc number of contacts
   *  \param[in] nb number of blocks
   *  \param[in] solver_id id an int to be mapped to actual solver in Numerics
   *  \param[in] tolerance threshold used to validate the solution: if the error is less than this value, the solution is accepted
   *  \param[in] itermax the maximum number of iteration
   *  \param[in] verbose level 0 : nothing, 1: mid level 2: high level
   *  \param[in] outputFile outputFile option 0 : nothing 1 : dat file 2: FCLIB HDF5 file if FCLIB is found
   *  \param[in] freq_output
   *  \param[in] ndof the numbe of dof in the dynamical systems involved in contact (for output in file.)
   *  \return result (0 if successful otherwise 1).
   *
   */
  int fc3d_LmgcDriver(double *reaction,
                                   double *velocity,
                                   double *q,
                                   double *mu,
                                   double* W,
                                   unsigned int *row,
                                   unsigned int *column,
                                   unsigned int nc,
                                   unsigned int nb,
                                   int solver_id,
                                   double tolerance,
                                   int itermax,
                                   int verbose,
                                   int outputFile,
                                   int freq_output,
                                   int ndof);

  /** LMGC interface to solvers for global friction-contact 3D problem
   *  \param[in,out] reaction global vector (nc*3)
   *  \param[in,out] velocity global vector (nc*3)
   *  \param[in,out] globalVelocity global velocity vector (n)
   *  \param[in] q global vector (n)
   *  \param[in] b global vector (nc*3)
   *  \param[in] mu global vector (nc)
   *  \param[in] Mdata the sparse matrix in coordinate format
   *  \param[in] nzM number of non zeros in Mdata
   *  \param[in] rowM  row indices of M
   *  \param[in] colM  column indices of M
   *  \param[in] Hdata the sparse matrix in coordinate format
   *  \param[in] nzH number of non zeros in Hdata
   *  \param[in] rowH  row indices of H
   *  \param[in] colH  column indices of H
   *  \param[in] n size of global velocity
   *  \param[in] nc number of contacts
   *  \param[in] solver_id id an int to be mapped to actual solver in Numerics
   *  \param[in] isize sive of integer parameters array
   *  \param[in] iparam integer parameters array
   *  \param[in] dsize sive of double parameters array
   *  \param[in] dparam double parameters array
   *  \param[in] verbose level 0 : nothing, 1: mid level 2: high level
   *  \param[in] outputFile outputFile option 0 : nothing 1 : C file , 1 :  dat file 3: FCLIB HDF5 file if FCLIB is found
   *  \param[in] freq_output
   *  \return result (0 if successful otherwise 1).
   *
   */
  int gfc3d_LmgcDriver(double *reaction,
                                         double *velocity,
                                         double *globalVelocity,
                                         double *q,
                                         double *b,
                                         double *mu,
                                         double *Mdata,
                                         unsigned int nzM,
                                         unsigned int *rowM,
                                         unsigned int *colM,
                                         double* Hdata,
                                         unsigned int nzH,
                                         unsigned int *rowH,
                                         unsigned int *colH,
                                         unsigned int n,
                                         unsigned int nc,
                                         int solver_id,
                                         int isize,
                                         int *iparam,
                                         int dsize,
                                         double *dparam,
                                         int verbose,
                                         int outputFile,
                                         int freq_output);

  /** General interface to solver for relay problems
      \param[in] problem the RelayProblem structure which handles the problem (M,q)
      \param[in,out] z a n-vector of doubles which contains the solution of the problem.
      \param[in,out] w a n-vector of doubles which contains the solution of the problem.
      \param[in,out] options structure used to define the solver(s) and its (their) parameters
      \return info termination value
      - 0 : successful
      - >0 : otherwise see each solver for more information about the log info
   */
  int relay_driver(RelayProblem* problem, double *z , double *w, SolverOptions* options);


  
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
