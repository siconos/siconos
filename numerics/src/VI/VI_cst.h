/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
/*!\file VI_cst.h */

/** \enum VI_SOLVER VI_cst.h
 * Enum that allows one to encode the list of solvers in a proper to avoid mispelling
 * with const char* const  variables
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
  /** index in iparam to store the error evaluation method */
  SICONOS_VI_IPARAM_LINESEARCH_METHOD = 6,
  /** index in iparam to store the error evaluation method */
  SICONOS_VI_IPARAM_ERROR_EVALUATION = 7,
  /** index in iparam to store the frequency of error evaluation method */
  SICONOS_VI_IPARAM_ERROR_EVALUATION_FREQUENCY = 8
};

enum SICONOS_VI_EG_DPARAM
{
  /** index in iparam to store the initial value of rho */
  SICONOS_VI_EG_DPARAM_RHO = 3,
  /** index in iparam to store the tau coeff of line-search */
  SICONOS_VI_EG_DPARAM_LS_TAU = 4,
  /** index in iparam to store the tauinv coeff of line-search */
  SICONOS_VI_EG_DPARAM_LS_TAUINV = 5,
  /** index in iparam to store the L coeff of line-search */
  SICONOS_VI_EG_DPARAM_LS_L = 6,
  /** index in iparam to store the LMIN coeff of line-search */
  SICONOS_VI_EG_DPARAM_LS_LMIN = 7,
};



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
