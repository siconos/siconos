/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

#ifndef VI_CST_H
#define VI_CST_H
/*!\file VI_cst.h
  \brief Enum and constants relative to Variational Inequality problems and solvers.
*/

/** \enum VI_SOLVER
    Available ids for VI solvers.
*/
enum VI_SOLVER
  {
   SICONOS_VI_EG = 1000,
   SICONOS_VI_FPP = 1001,
   SICONOS_VI_HP = 1002,
   SICONOS_VI_BOX_QI = 1020,
   SICONOS_VI_BOX_AVI_LSA = 1021,
   SICONOS_VI_BOX_PATH = 1022
  };

enum SICONOS_VI_IPARAM
  {
   SICONOS_VI_IPARAM_DECREASE_RHO = 3,
   /** index in iparam to store the linesearch method */
   SICONOS_VI_IPARAM_LINESEARCH_METHOD = 6,
   /** index in iparam to store the error evaluation method */
   SICONOS_VI_IPARAM_ERROR_EVALUATION = 7,
   /** index in iparam to store the frequency of error evaluation method */
   SICONOS_VI_IPARAM_ERROR_EVALUATION_FREQUENCY = 8,
   /** index for max number of iterations allowed in line search*/
   SICONOS_VI_IPARAM_LS_MAX_ITER = 9,
   /** activate the update in the loop (0:false default choice) */
   SICONOS_VI_IPARAM_ACTIVATE_UPDATE = 10,
  };

/** \enum allowed values for iparam[SICONOS_VI_IPARAM_LINESEARCH_METHOD]   
*/
enum SICONOS_VI_LINESEARCH_METHOD
  {
   /**  : Armijo rule with Khotbotov ratio (default) */
   SICONOS_VI_LS_ARMIJO = 0,
   /** : Armijo rule with Solodov.Tseng ratio */
   SICONOS_VI_LS_SOLODOV = 1,
   /**  Armijo rule with Han.Sun ratio */
   SICONOS_VI_LS_HANSUN = 2,
  };
  
enum SICONOS_VI_DPARAM
{
  /** index in dparam to store the initial value of rho */
  SICONOS_VI_DPARAM_RHO = 3,
  /** index in dparam to store the tau coeff of line-search */
  SICONOS_VI_DPARAM_LS_TAU = 4,
  /** index in dparam to store the tauinv coeff of line-search */
  SICONOS_VI_DPARAM_LS_TAUINV = 5,
  /** index in dparam to store the L coeff of line-search */
  SICONOS_VI_DPARAM_LS_L = 6,
  /** index in dparam to store the LMIN coeff of line-search */
  SICONOS_VI_DPARAM_LS_LMIN = 7,
  /** index in dparam to store the sigma coeff (HP) */
  SICONOS_VI_DPARAM_SIGMA = 8,
};



/** \enum allowed values for iparam[SICONOS_VI_IPARAM_ERROR_EVALUATION]   
*/
enum SICONOS_VI_ERROR_EVALUATION_ENUM
{
  SICONOS_VI_ERROR_EVALUATION_FULL = 0,
  SICONOS_VI_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL = 1,
  SICONOS_VI_ERROR_EVALUATION_LIGHT = 2,
  SICONOS_VI_ERROR_EVALUATION_ADAPTIVE =3,
};



extern const char* const   SICONOS_VI_EG_STR ;
extern const char* const   SICONOS_VI_FPP_STR ;
extern const char* const   SICONOS_VI_HP_STR ;
extern const char* const   SICONOS_VI_BOX_QI_STR ;
extern const char* const   SICONOS_VI_BOX_AVI_LSA_STR ;
extern const char* const   SICONOS_VI_BOX_PATH_STR ;

#endif
