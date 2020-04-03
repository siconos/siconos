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

#ifndef LCP_CST_H
#define LCP_CST_H
/*!\file lcp_cst.h
  \brief Constants to define the list of available LCP solvers.

  \rst Check the detailed documentation in :ref:`lcp_solvers`. \endrst
*/

/**\enum LCP_SOLVER
   Each SICONOS_LCP_XXX refers to number of the solver XXX for LCP.
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
  SICONOS_LCP_NEWTON_FB_FBLSA = 209,
  SICONOS_LCP_PSOR = 210,
  SICONOS_LCP_RPGS = 211,
  SICONOS_LCP_PATH = 212,
  SICONOS_LCP_ENUM = 213,
  SICONOS_LCP_AVI_CAOFERRIS = 214,
  SICONOS_LCP_PIVOT = 215,
  SICONOS_LCP_BARD = 216,
  SICONOS_LCP_MURTY = 217,
  SICONOS_LCP_NEWTON_MIN_FBLSA = 218,
  SICONOS_LCP_PATHSEARCH = 219,
  SICONOS_LCP_PIVOT_LUMOD = 220,
  SICONOS_LCP_GAMS = 221,
  SICONOS_LCP_CONVEXQP_PG = 222
};


enum SICONOS_LCP_IPARAM
  {
   /** index in iparam to store the sum of local solver iterations number */
   SICONOS_LCP_IPARAM_NSGS_ITERATIONS_SUM =4,
   /** index in iparam to store type of pivoting methods */
   SICONOS_LCP_IPARAM_PIVOTING_METHOD_TYPE =5,
   /** index in iparam to skip trivial solution */
   SICONOS_LCP_IPARAM_SKIP_TRIVIAL =6,
   /** index in iparam to store the number of solutions */
   SICONOS_LCP_IPARAM_ENUM_NUMBER_OF_SOLUTIONS =7,
  /** index in iparam to store  the current enum */
   SICONOS_LCP_IPARAM_ENUM_CURRENT_ENUM =8,
   /** index in iparam to store the seed for starting enum*/
   SICONOS_LCP_IPARAM_ENUM_SEED =9,
   /** index in iparam to store the first seed for enum */
   SICONOS_LCP_IPARAM_ENUM_USE_DGELS =10,
   /** index in iparam to store to activate multiple solutions search */
   SICONOS_LCP_IPARAM_ENUM_MULTIPLE_SOLUTIONS =11,
   /** **/
   
  };

enum SICONOS_LCP_DPARAM
{
  /** index in dparam to store to the relaxation or regularization parameter */
  SICONOS_LCP_DPARAM_RHO =3,
  /** index in dparam to store the sum of local error values */
  SICONOS_LCP_DPARAM_NSGS_LOCAL_ERROR_SUM =4,
  /** index in dparam to store the latin parameter */
  SICONOS_LCP_DPARAM_LATIN_PARAMETER =12
};

enum SICONOS_LCP_SKIP_TRIVIAL
{
  SICONOS_LCP_SKIP_TRIVIAL_NO=0,
  SICONOS_LCP_SKIP_TRIVIAL_YES=1
};


/** Allowed values for iparam[SICONOS_LCP_IPARAM_PIVOTING_METHOD_TYPE] */
enum SICONOS_LCP_PIVOT_TYPE
{
  SICONOS_LCP_PIVOT_BARD = 1,
  SICONOS_LCP_PIVOT_LEAST_INDEX = 2,
  SICONOS_LCP_PIVOT_LEMKE = 3,
  SICONOS_LCP_PIVOT_PATHSEARCH = 4
};

extern const char* const   SICONOS_LCP_LEMKE_STR;
extern const char* const   SICONOS_LCP_NSGS_SBM_STR;
extern const char* const   SICONOS_LCP_PGS_STR;
extern const char* const   SICONOS_LCP_CPG_STR;
extern const char* const   SICONOS_LCP_LATIN_STR;
extern const char* const   SICONOS_LCP_LATIN_W_STR;
extern const char* const   SICONOS_LCP_QP_STR;
extern const char* const   SICONOS_LCP_NSQP_STR;
extern const char* const   SICONOS_LCP_NEWTONMIN_STR;
extern const char* const   SICONOS_LCP_NEWTON_FB_FBLSA_STR;
extern const char* const   SICONOS_LCP_NEWTON_MIN_FBLSA_STR;
extern const char* const   SICONOS_LCP_PSOR_STR;
extern const char* const   SICONOS_LCP_RPGS_STR;
extern const char* const   SICONOS_LCP_PATH_STR;
extern const char* const   SICONOS_LCP_ENUM_STR;
extern const char* const   SICONOS_LCP_AVI_CAOFERRIS_STR;
extern const char* const   SICONOS_LCP_PIVOT_STR;
extern const char* const   SICONOS_LCP_BARD_STR;
extern const char* const   SICONOS_LCP_MURTY_STR;
extern const char* const   SICONOS_LCP_PATHSEARCH_STR;
extern const char* const   SICONOS_LCP_PIVOT_LUMOD_STR;
extern const char* const   SICONOS_LCP_GAMS_STR;
extern const char* const   SICONOS_LCP_CONVEXQP_PG_STR;
#endif
