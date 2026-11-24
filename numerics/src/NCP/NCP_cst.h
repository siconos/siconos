/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef NCP_CST_H
#define NCP_CST_H

/*!\file NCP_cst.h
  \brief Constants to define the list of available NCP solvers. See the solver list \ref
  ncpSolversList
*/

/**\enum NCP_SOLVER
   Each SICONOS_NCP_XXX refers to number of the solver XXX for NCP. See the solver list \ref
   ncpSolversList
 */
enum NCP_SOLVER {
  SICONOS_NCP_NEWTON_FB_FBLSA = 900,
  SICONOS_NCP_NEWTON_MIN_FBLSA = 901,
  SICONOS_NCP_PATHSEARCH = 902,
  SICONOS_NCP_PATH = 903
};

#endif
