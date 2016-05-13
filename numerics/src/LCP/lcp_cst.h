/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.

 * Copyright 2016 INRIA.

 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at

 * http://www.apache.org/licenses/LICENSE-2.0

 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#ifndef LCP_CST_H
#define LCP_CST_H
/*!\file lcp_cst.h
  \brief Constants to define the list of available LCP solvers. See the solver list \ref lcpSolversList
*/
/**\enum LCP_SOLVER
   Each SICONOS_LCP_XXX refers to number of the solver XXX for LCP. See the solver list \ref lcpSolversList
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
  SICONOS_LCP_NEWTON_FBLSA = 209,
  SICONOS_LCP_PSOR = 210,
  SICONOS_LCP_RPGS = 211,
  SICONOS_LCP_PATH = 212,
  SICONOS_LCP_ENUM = 213,
  SICONOS_LCP_AVI_CAOFERRIS = 214,
  SICONOS_LCP_PIVOT = 215,
  SICONOS_LCP_BARD = 216,
  SICONOS_LCP_MURTY = 217,
  SICONOS_LCP_NEWTON_MINFBLSA = 218,
  SICONOS_LCP_PATHSEARCH = 219,
  SICONOS_LCP_PIVOT_LUMOD = 220,
  SICONOS_LCP_GAMS = 221
};

enum LCP_PIVOT
{
  SICONOS_LCP_PIVOT_BARD = 1,
  SICONOS_LCP_PIVOT_LEAST_INDEX = 2,
  SICONOS_LCP_PIVOT_LEMKE = 3,
  SICONOS_LCP_PIVOT_PATHSEARCH = 4
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
extern char *  SICONOS_LCP_NEWTON_FBLSA_STR;
extern char *  SICONOS_LCP_NEWTON_MINFBLSA_STR;
extern char *  SICONOS_LCP_PSOR_STR;
extern char *  SICONOS_LCP_RPGS_STR;
extern char *  SICONOS_LCP_PATH_STR;
extern char *  SICONOS_LCP_ENUM_STR;
extern char *  SICONOS_LCP_AVI_CAOFERRIS_STR;
extern char *  SICONOS_LCP_PIVOT_STR;
extern char *  SICONOS_LCP_BARD_STR;
extern char *  SICONOS_LCP_MURTY_STR;
extern char *  SICONOS_LCP_PATHSEARCH_STR;
extern char *  SICONOS_LCP_PIVOT_LUMOD_STR;
extern char *  SICONOS_LCP_GAMS_STR;
#endif
