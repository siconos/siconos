/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

#ifndef CONVEXQP_CST_H
#define CONVEXQP_CST_H

/*!\file ConvexQP_cst.h 
   \brief Enum and constants relative to QP solvers.
*/


/** \enum CONVEXQP_SOLVER 
    List of available ids for QP solvers.
*/
enum CONVEXQP_SOLVER
{
 /** convex QP, projected gradient */
  SICONOS_CONVEXQP_PG = 1200,
  /** convex QP, VI solver, fixed point projection */
  SICONOS_CONVEXQP_VI_FPP= 1201,
  /** convex QP, VI solver, extra-gradient */
  SICONOS_CONVEXQP_VI_EG= 1202,
  /** convex QP, alternating direction method of multipliers (ADMM) */
  SICONOS_CONVEXQP_ADMM= 1203
};

extern const char* const   SICONOS_CONVEXQP_PG_STR;
extern const char* const   SICONOS_CONVEXQP_VI_FPP_STR;
extern const char* const   SICONOS_CONVEXQP_VI_EG_STR;
extern const char* const   SICONOS_CONVEXQP_ADMM_STR;


/** iparam indices specific to QP solvers. */
enum SICONOS_CONVEXQP_PGOC_IPARAM_ENUM
{
  /** index in iparam to store the maximum number of iterations */
  SICONOS_CONVEXQP_PGOC_LINESEARCH_MAX_ITER = 10
};

/** dparam indices specific to QP solvers. */
enum SICONOS_CONVEXQP_PGOC_DPARAM_ENUM
{
  /** index in dparam to store the rho value for projection formulation */
  SICONOS_CONVEXQP_PGOC_RHO = 3,
  /** index in dparam to store the minrho value for projection formulation */
  SICONOS_CONVEXQP_PGOC_RHOMIN = 4,
  /** index in dparam to store the mu value for line search algo */
  SICONOS_CONVEXQP_PGOC_LINESEARCH_MU = 5,
  /** index in dparam to store the tau value for line search algo */
  SICONOS_CONVEXQP_PGOC_LINESEARCH_TAU  = 6 
};

enum SICONOS_CONVEXQP_ADMM_IPARAM_ENUM
{
  /** index in iparam to store the strategy for computing rho */
  SICONOS_CONVEXQP_ADMM_IPARAM_RHO_STRATEGY = 9,
  /** index in iparam to store the acceleration paramter */
  SICONOS_CONVEXQP_ADMM_IPARAM_ACCELERATION= 10

};

enum SICONOS_CONVEXQP_ADMM_DPARAM_ENUM
{
  /** index in dparam to store the rho value for projection formulation */
  SICONOS_CONVEXQP_ADMM_RHO = 3,
  /** index in dparam to store the rho value for projection formulation */
  SICONOS_CONVEXQP_ADMM_RESTART_ETA = 4,
  /** index in dparam to store the tau value for the balancing residual technique */
  SICONOS_CONVEXQP_ADMM_BALANCING_RESIDUAL_TAU = 5,
  /** index in dparam to store the phi value for the balancing residual technique */
  SICONOS_CONVEXQP_ADMM_BALANCING_RESIDUAL_PHI = 6

};

enum SICONOS_CONVEXQP_ADMM_ACCELERATION_ENUM
{
  SICONOS_CONVEXQP_ADMM_NO_ACCELERATION= 0,
  SICONOS_CONVEXQP_ADMM_ACCELERATION= 1,
  SICONOS_CONVEXQP_ADMM_ACCELERATION_AND_RESTART= 2
};

enum SICONOS_CONVEXQP_RHO_STRATEGY_ENUM
{
  /** A constant value given in dparam[CONVEXQP_RHO_RHO] is used */
  SICONOS_CONVEXQP_RHO_STRATEGY_CONSTANT = 0,
  /** A computed value stored in dparam[SICONOS_CONVEXQP_NSN_RHO] is used */
  SICONOS_CONVEXQP_ADMM_RHO_STRATEGY_NORM_INF =1,
  /** An adaptive strategy for rho is used */
  SICONOS_CONVEXQP_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING =2,
  /** An adaptive strategy for rho is used */
  SICONOS_CONVEXQP_ADMM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING =3
};




#endif
