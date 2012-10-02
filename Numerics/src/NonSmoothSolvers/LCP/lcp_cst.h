/* Siconos-Numerics, Copyright INRIA 2005-2012.
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

#ifndef LCP_CST_H
#define LCP_CST_H
/*!\file lcp_cst.h
  \brief Constants to define the list of available LCP solvers. See the solver list \ref lcpSolversList
*/
/**\enum LCP_SOLVER
   \brief Each SICONOS_LCP_XXX refers to number of the solver XXX for LCP. See the solver list \ref lcpSolversList
 */
enum LCP_SOLVER
{
  SICONOS_LCP_LEMKE = 200,
  SICONOS_LCP_NSGS_SBM = 201,
  SICONOS_LCP_PGS = 202,
  SICONOS_LCP_CPG = 203,
  SICONOS_LCP_LATIN = 204,
  SICONOS_LCP_LATIN_W = 205,
  SICONOS_LCP_QP = 206,
  SICONOS_LCP_NSQP = 207,
  SICONOS_LCP_NEWTONMIN = 208,
  SICONOS_LCP_NEWTONFB = 209,
  SICONOS_LCP_PSOR = 210,
  SICONOS_LCP_RPGS = 211,
  SICONOS_LCP_PATH = 212,
  SICONOS_LCP_ENUM = 213
};



extern char *  SICONOS_LCP_LEMKE_STR;
extern char *  SICONOS_LCP_NSGS_SBM_STR;
extern char *  SICONOS_LCP_PGS_STR;
extern char *  SICONOS_LCP_CPG_STR;
extern char *  SICONOS_LCP_LATIN_STR;
extern char *  SICONOS_LCP_LATIN_W_STR;
extern char *  SICONOS_LCP_QP_STR;
extern char *  SICONOS_LCP_NSQP_STR;
extern char *  SICONOS_LCP_NEWTONMIN_STR;
extern char *  SICONOS_LCP_NEWTONFB_STR;
extern char *  SICONOS_LCP_PSOR_STR;
extern char *  SICONOS_LCP_RPGS_STR;
extern char *  SICONOS_LCP_PATH_STR;
extern char *  SICONOS_LCP_ENUM_STR;

#endif
