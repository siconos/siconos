#ifndef CONVEXQP_CST_H
#define CONVEXQP_CST_H
/** \file ConvexQP_cst.h */


/** \enum CONVEXQP_SOLVER encode the list of solvers as integers, to avoid mispelling
 * with const char* const  variables
 */
enum CONVEXQP_SOLVER
{
  SICONOS_CONVEXQP_PG = 1200,
  SICONOS_CONVEXQP_VI_FPP= 1201,
  SICONOS_CONVEXQP_VI_EG= 1202,
  SICONOS_CONVEXQP_ADMM= 1203
};

extern const char* const   SICONOS_CONVEXQP_PG_STR;
extern const char* const   SICONOS_CONVEXQP_VI_FPP_STR;
extern const char* const   SICONOS_CONVEXQP_VI_EG_STR;
extern const char* const   SICONOS_CONVEXQP_ADMM_STR;

enum SICONOS_CONVEXQP_PGOC_IPARAM_ENUM
{
  /** index in iparam to store the maximum number of iterations */
  SICONOS_CONVEXQP_PGOC_LINESEARCH_MAXITER = 10
};
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
