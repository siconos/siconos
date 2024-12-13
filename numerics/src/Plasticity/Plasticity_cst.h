#ifndef PLASTICITY_CST_H
#define PLASTICITY_CST_H
/** \file PLASTICITY_cst.h */
/** \enum PLASTICITY_SOLVER encode the list of solvers as integers, to avoid mispelling
 * with const char* const  variables
 */
enum PLASTICITY_SOLVER {
  /** 2D Mohr Coulomb solvers */
  MOHR_COULOMB_2D_NSGS = 20000,

  /* 2D Mohr Coulomb  solvers for one cone (used mainly inside NSGS solvers) */

  /** Non-smooth Newton Alart-Curnier, 'direct', one contact solver */
  MOHR_COULOMB_2D_ONECONE_NSN = 20050,
  /** Non-smooth Newton Alart-Curnier, 'damped', one contact solver */
  MOHR_COULOMB_2D_ONECONE_NSN_GP = 20051,
  /** Projection on cone, one contact solver */
  MOHR_COULOMB_2D_ONECONE_ProjectionOnCone = 20052,
  /** Projection on cone, one contact solver */
  MOHR_COULOMB_2D_ONECONE_ProjectionOnConeWithLocalIteration = 20053,
  /** Non-smooth Newton Alart-Curnier, 'damped' and hybrid with projection, one contact solver
   */
  MOHR_COULOMB_2D_ONECONE_NSN_GP_HYBRID = 20065
};

extern const char* const MOHR_COULOMB_2D__NSGS_STR;

enum PLASTICITY_IPARAM {
  /** index in iparam to store the error strategy for the internal solver */
  PLASTICITY_IPARAM_INTERNAL_ERROR_STRATEGY = 2,
  /** index in iparam to store the rescaling  */
  PLASTICITY_IPARAM_RESCALING = 3,
  /** index in iparam to store the rescaling  */
  PLASTICITY_IPARAM_RESCALING_CONE = 4,
  /** current contact number (example of use: one contact solvers) **/
  PLASTICITY_CURRENT_CONE_NUMBER = 5,
  /** index in iparam to store the error evaluation method */
  PLASTICITY_IPARAM_ERROR_EVALUATION = 7,
  /** index in iparam to store the frequency of error evaluation method */
  PLASTICITY_IPARAM_ERROR_EVALUATION_FREQUENCY = 8,
  PLASTICITY_NUMBER_OF_CONTACTS = 17,
};

enum PLASTICITY_INTERNAL_ERROR_STRATEGY {
  PLASTICITY_INTERNAL_ERROR_STRATEGY_ADAPTIVE = 0,
  PLASTICITY_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE = 1,
  PLASTICITY_INTERNAL_ERROR_STRATEGY_ADAPTIVE_N_CONTACT = 2
};

enum PLASTICITY_RESCALING_ENUM {
  PLASTICITY_RESCALING_NO = 0,
  PLASTICITY_RESCALING_SCALAR = 1,
  PLASTICITY_RESCALING_BALANCING_M = 2,
  PLASTICITY_RESCALING_BALANCING_MH = 3,
  PLASTICITY_RESCALING_BALANCING_MHHT = 4
};

enum PLASTICITY_RESCALING_CONE_ENUM {
  PLASTICITY_RESCALING_CONE_NO = 0,
  PLASTICITY_RESCALING_CONE_YES = 1
};

enum PLASTICITY_DPARAM {
  /** index in dparam to store the internal solver error ratio*/
  PLASTICITY_DPARAM_INTERNAL_ERROR_RATIO = 2
};

enum PLASTICITY_NSGS_IPARAM {
  /** index in iparam to store the relaxation strategy */
  PLASTICITY_NSGS_RELAXATION = 4,
  /** index in iparam to store the shuffle strategy */
  PLASTICITY_NSGS_SHUFFLE = 5,
  /** index in iparam to store the shuffle seed */
  PLASTICITY_NSGS_SHUFFLE_SEED = 6,
  /** index in iparam to store the  */
  PLASTICITY_NSGS_FREEZING_CONTACT = 19,
  /** index in iparam to store the  */
  PLASTICITY_NSGS_FILTER_LOCAL_SOLUTION = 14,
};
enum PLASTICITY_NSGS_DPARAM {
  /** index in dparam to store the relaxation strategy */
  PLASTICITY_NSGS_RELAXATION_VALUE = 8,
};

enum PLASTICITY_NSGS_LOCALSOLVER_IPARAM {
  PLASTICITY_NSGS_LOCALSOLVER_IPARAM_USE_TRIVIAL_SOLUTION = 10
};

enum PLASTICITY_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION {
  PLASTICITY_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_FALSE = 0,
  PLASTICITY_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_TRUE = 1
};

enum PLASTICITY_NSGS_ERROR_EVALUATION_ENUM {
  /** Evaluation of the error with the expensive function fc3d_compute_error **/
  PLASTICITY_NSGS_ERROR_EVALUATION_FULL = 0,
  /** Evaluation of the error with the cheap incremental variation **/
  PLASTICITY_NSGS_ERROR_EVALUATION_LIGHT = 1,
  /** Evaluation of the error with the cheap incremental variation but we modify
      the incremental toleranve to reach the requred accuracy **/
  PLASTICITY_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL = 2,
  /** Evaluation of the error with the expensive function fc3d_compute_error and
      an adaptive frequency for calling the error function  **/
  PLASTICITY_NSGS_ERROR_EVALUATION_ADAPTIVE = 3,
};
enum PLASTICITY_NSGS_SHUFFLE_ENUM {
  PLASTICITY_NSGS_SHUFFLE_FALSE = 0,
  PLASTICITY_NSGS_SHUFFLE_TRUE = 1,
  PLASTICITY_NSGS_SHUFFLE_TRUE_EACH_LOOP = 2
};

enum PLASTICITY_NSGS_RELAXATION_ENUM {
  PLASTICITY_NSGS_RELAXATION_FALSE,
  PLASTICITY_NSGS_RELAXATION_TRUE
};
enum PLASTICITY_NSGS_FILTER_LOCAL_SOLUTION_ENUM {

  PLASTICITY_NSGS_FILTER_LOCAL_SOLUTION_FALSE = 0,
  PLASTICITY_NSGS_FILTER_LOCAL_SOLUTION_TRUE = 1
};

enum PLASTICITY_NSN_IPARAM {
  /** index in iparam to store the strategy for computing rho */
  PLASTICITY_NSN_RHO_STRATEGY = 9,
  /** index in iparam to store the formulation */
  PLASTICITY_NSN_FORMULATION = 10,
  /** index in iparam to store the line-search */
  PLASTICITY_NSN_LINESEARCH = 11,
  /** index in iparam to store the maximum number of iterations */
  PLASTICITY_NSN_LINESEARCH_MAX_ITER = 12,
  /** index in iparam to set the linear solver used at each Newton iteration
   cs_lusol or mumps */
  PLASTICITY_NSN_LINEAR_SOLVER = 13,
  /** index in iparam to store the strategy for the hybrid solver */
  PLASTICITY_NSN_HYBRID_STRATEGY = 14,
  /** index in iparam to store the maximum number of loop for the hybrid solver */
  PLASTICITY_NSN_HYBRID_MAX_LOOP = 15,
  /** index in iparam to store the maximum number of iterations for the projection solver */
  PLASTICITY_NSN_HYBRID_MAX_ITER = 16,
  /** index in iparam used to check if memory allocation has already be done (if true/1) or not
     (if 0/false) for internal work array. */
  PLASTICITY_NSN_MEMORY_ALLOCATED = 17,
  /** index in iparam to store the boolean to know if allocation of dwork is needed */
  PLASTICITY_NSN_MPI_COM = 18

};

enum PLASTICITY_NSN_LINEAR_SOLVER {
  PLASTICITY_NSN_USE_CSLUSOL = 0,
  PLASTICITY_NSN_USE_MUMPS = 1
};

enum PLASTICITY_NSN_DPARAM {
  /** index in dparam to store the rho value for projection formulation */
  PLASTICITY_NSN_RHO = 3,
};

enum PLASTICITY_NSN_RHO_STRATEGY_ENUM {
  /** A constant value given in dparam[PLASTICITY_NSN_RHO] is used */
  PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_CONSTANT = 0,
  /** A computed value stored in dparam[PLASTICITY_NSN_RHO] is used */
  PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM = 1,
  /** A computed value stored in dparam[PLASTICITY_NSN_RHO] is used */
  PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM_COND = 2,
  /** An adaptive strategy for rho is used */
  PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM = 3,
  /** An adaptive strategy for rho is used */
  PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_ADAPTIVE = 4,
};

enum PLASTICITY_NSN_FORMULATION_ENUM {
  PLASTICITY_NSN_FORMULATION_ALARTCURNIER_STD = 0,
  PLASTICITY_NSN_FORMULATION_NATURALMAP = 1,
  PLASTICITY_NSN_FORMULATION_NULL = 2,
};

enum PLASTICITY_NSN_LINESEARCH_ENUM {
  PLASTICITY_NSN_LINESEARCH_GOLDSTEINPRICE = 0,
  PLASTICITY_NSN_LINESEARCH_ARMIJO = 1,
  PLASTICITY_NSN_LINESEARCH_NO = -1,
};

enum PLASTICITY_NSN_HYBRID_ENUM {
  /** No strategy for the hybrid solver */
  PLASTICITY_NSN_HYBRID_STRATEGY_NO = 0,
  /** Loop PLI-NSN strategy for the hybrid solver */
  PLASTICITY_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP = 1,
  /** NSN and after Loop PLI-NSN strategy for the hybrid solver */
  PLASTICITY_NSN_HYBRID_STRATEGY_NSN_AND_PLI_NSN_LOOP = 2,
  /** VI_EG preconditionning to NSN strategy for the hybrid solver */
  PLASTICITY_NSN_HYBRID_STRATEGY_VI_EG_NSN = 3
};

#endif
