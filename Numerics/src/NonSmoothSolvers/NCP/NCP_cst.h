/* Siconos-Numerics, Copyright INRIA 2005-2014
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

#ifndef NCP_CST_H
#define NCP_CST_H

/*!\file NCP_cst.h
  \brief Constants to define the list of available NCP solvers. See the solver list \ref ncpSolversList
*/

/**\enum NCP_SOLVER
   Each SICONOS_NCP_XXX refers to number of the solver XXX for NCP. See the solver list \ref ncpSolversList
 */
enum NCP_SOLVER
{
  SICONOS_NCP_NEWTON_FBLSA = 900,
  SICONOS_NCP_NEWTON_MINFBLSA = 901,
  SICONOS_NCP_PATHSEARCH = 902,
  SICONOS_NCP_PATH = 903
};


extern char SICONOS_NCP_NEWTON_FBLSA_STR[];
extern char SICONOS_NCP_NEWTON_MINFBLSA_STR[];
extern char SICONOS_NCP_PATHSEARCH_STR[];
extern char SICONOS_NCP_PATH_STR[];

#endif
