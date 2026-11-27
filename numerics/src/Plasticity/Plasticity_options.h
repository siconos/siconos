// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef PLASTICITY_CST_H
#define PLASTICITY_CST_H

#include "natural_map_error.h"  // for the NATURAL_MAP_ERROR_* scaling options
/** \file PLASTICITY_cst.h */
/** \enum PLASTICITY_SOLVER encode the list of solvers as integers, to avoid mispelling
 * with const char* const  variables
 */
enum PLASTICITY_SOLVER {
  /** 2D Mohr Coulomb solvers */
  PLASTICITY_2D_NSGS = 20000,
  /** 2D Mohr Coulomb solvers using generic NSGS framework */
  PLASTICITY_2D_NSGS_GENERIC = 20001,
  /** Global semi-smooth Newton method on the natural map of the whole problem (sparse LU of
      the generalized Jacobian, Armijo line search, Levenberg-Marquardt fallback) */
  PLASTICITY_2D_NSN = 20002,

  /* 2D Mohr Coulomb  solvers for one cone (used mainly inside NSGS solvers) */

  /** Non-smooth Newton Alart-Curnier, 'direct', one cone solver */
  PLASTICITY_2D_ONECONE_NSN = 20050,
  /** Non-smooth Newton Alart-Curnier, 'damped', one cone solver */
  PLASTICITY_2D_ONECONE_NSN_GP = 20051,
  /** Projection on cone, one cone solver */
  PLASTICITY_2D_ONECONE_ProjectionOnCone = 20052,
  /** Projection on cone, one cone solver */
  PLASTICITY_2D_ONECONE_ProjectionOnConeWithLocalIteration = 20053,
  /** Non-smooth Newton Alart-Curnier, 'damped' and hybrid with projection, one cone solver
   */
  PLASTICITY_2D_ONECONE_NSN_GP_HYBRID = 20065,

  /** Von Mises (J2), one element: return mapping (exact solution, closest-point projection
      of the trial stress -M^{-1} q on the yield set in the M-norm) */
  PLASTICITY_2D_ONECONE_VONMISES_RETURN_MAPPING = 20070,
  /** Von Mises (J2), one element: fixed point on the natural map
      sigma <- P_C(sigma - rho (M sigma + q)) (rho strategy: PLASTICITY_NSN_RHO_STRATEGY) */
  PLASTICITY_2D_ONECONE_VONMISES_FP = 20071,
  /** Von Mises (J2), one element: semi-smooth Newton on the natural map
      sigma - P_C(sigma - rho (M sigma + q)) = 0, with line search */
  PLASTICITY_2D_ONECONE_VONMISES_NSN = 20072,
  /** Drucker-Prager, one cone: return mapping (exact solution by enumeration: elastic,
      apex, sliding with a quartic equation), any M and dilatancy */
  PLASTICITY_2D_ONECONE_DRUCKERPRAGER_RETURN_MAPPING = 20080,
  /** Drucker-Prager, one cone: semi-smooth Newton on the natural map with line search and
      fallback to the return mapping */
  PLASTICITY_2D_ONECONE_DRUCKERPRAGER_NSN = 20081,
  /** deprecated name of PLASTICITY_2D_ONECONE_VONMISES_RETURN_MAPPING */
  PLASTICITY_2D_ONECONE_VONMISES = PLASTICITY_2D_ONECONE_VONMISES_RETURN_MAPPING
};

enum PLASTICITY_IPARAM {
  /** index in iparam to store the error strategy for the internal solver */
  PLASTICITY_IPARAM_INTERNAL_ERROR_STRATEGY = 2,
  /** index in iparam to store the rescaling  */
  PLASTICITY_IPARAM_RESCALING = 3,
  /** index in iparam to store the rescaling  */
  PLASTICITY_IPARAM_RESCALING_CONE = 4,
  /** current cone number (example of use: one cone solvers) **/
  PLASTICITY_CURRENT_CONE_NUMBER = 5,
  /** index in iparam to store the error evaluation method */
  PLASTICITY_IPARAM_ERROR_EVALUATION = 7,
  /** index in iparam to store the frequency of error evaluation method */
  PLASTICITY_IPARAM_ERROR_EVALUATION_FREQUENCY = 8,
  PLASTICITY_NUMBER_OF_CONES = 17,
  /** Norm used to scale the error of plasticity_2d_compute_error (see PLASTICITY_SCALING;
      same index and values as SICONOS_FRICTION_3D_SCALING_NORM) */
  PLASTICITY_SCALING_NORM = NATURAL_MAP_ERROR_IPARAM_NORM,
  /** Scope of the scaling: per element or whole vector (see PLASTICITY_SCALING_SCOPE) */
  PLASTICITY_SCALING_SCOPE = NATURAL_MAP_ERROR_IPARAM_SCOPE,
  /** Reduction of the two norms in the scaling (see PLASTICITY_SCALING_REDUCTION) */
  PLASTICITY_SCALING_REDUCTION = NATURAL_MAP_ERROR_IPARAM_REDUCTION,
  /** rho of the natural map of the error, balancing the stress and the strain rate (see
      PLASTICITY_SCALING_RHO_ENUM) */
  PLASTICITY_SCALING_RHO = 27,
};

/** rho of the natural map r - P_K(r - rho v) used by plasticity_2d_compute_error: with
    rho ~ 1/||M||, rho v is a stress and the scaling of the error compares quantities of the
    same magnitude. */
enum PLASTICITY_SCALING_RHO_ENUM {
  /** rho = 1 (as in fc3d) */
  PLASTICITY_SCALING_RHO_ONE = 0,
  /** rho_e = 1/||M_ee||_F for each element (diagonal block in the variables of the solvers):
      the error does not depend on the units of the strain rate (default of the Plasticity
      NSGS) */
  PLASTICITY_SCALING_RHO_ELEMENT = 1,
  /** rho = 1/max_e ||M_ee||_F for all the elements */
  PLASTICITY_SCALING_RHO_GLOBAL = 2,
};

/** Norm for the scaling of the error (plasticity_2d_compute_error), as in fc3d. */
enum PLASTICITY_SCALING {
  /** Euclidean norm (default) */
  PLASTICITY_SCALING_EUCLIDEAN = NATURAL_MAP_ERROR_NORM_EUCLIDEAN,
  /** no scaling: absolute error */
  PLASTICITY_SCALING_NO = NATURAL_MAP_ERROR_NORM_NO,
  /** cone norm |x_N| + ||x_T|| */
  PLASTICITY_SCALING_CONE = NATURAL_MAP_ERROR_NORM_CONE,
  /** cone norm |x_N| + eta ||x_T|| (Von Mises: eta = 1) */
  PLASTICITY_SCALING_CONE_ETA = NATURAL_MAP_ERROR_NORM_CONE_MU,
  /** sqrt(x_N^2 + alpha ||x_T||^2), alpha = dparam[PLASTICITY_SCALING_ALPHA] */
  PLASTICITY_SCALING_CONE_WEIGHTED = NATURAL_MAP_ERROR_NORM_CONE_WEIGHTED,
};

/** Scope of the scaling of the error. */
enum PLASTICITY_SCALING_SCOPE {
  /** one scale per element (default) */
  PLASTICITY_SCALING_SCOPE_PER_CONE = NATURAL_MAP_ERROR_SCOPE_PER_BLOCK,
  /** one scale for the whole vectors */
  PLASTICITY_SCALING_SCOPE_PER_VECTOR = NATURAL_MAP_ERROR_SCOPE_PER_VECTOR,
};

/** Combination of the (safeguarded) norms of the stress and of the modified strain rate. */
enum PLASTICITY_SCALING_REDUCTION {
  /** min(max(1, |r|), max(1, |v|)) (named MAX as in fc3d, its default) */
  PLASTICITY_SCALING_REDUCTION_MAX = NATURAL_MAP_ERROR_REDUCTION_MAX,
  /** geometric mean of max(1, |r|) and max(1, |v|) */
  PLASTICITY_SCALING_REDUCTION_MEAN = NATURAL_MAP_ERROR_REDUCTION_MEAN,
  /** |r| + |v| (safeguarded) (default of the Plasticity NSGS) */
  PLASTICITY_SCALING_REDUCTION_MIXED = NATURAL_MAP_ERROR_REDUCTION_MIXED,
};

enum PLASTICITY_INTERNAL_ERROR_STRATEGY {
  PLASTICITY_INTERNAL_ERROR_STRATEGY_ADAPTIVE = 0,
  PLASTICITY_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE = 1,
  PLASTICITY_INTERNAL_ERROR_STRATEGY_ADAPTIVE_N_CONE = 2
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
  PLASTICITY_DPARAM_INTERNAL_ERROR_RATIO = 2,
  /** weight alpha of PLASTICITY_SCALING_CONE_WEIGHTED (same index as
      SICONOS_FRICTION_3D_SCALING_ALPHA) */
  PLASTICITY_SCALING_ALPHA = NATURAL_MAP_ERROR_DPARAM_ALPHA
};

enum PLASTICITY_NSGS_IPARAM {
  /** index in iparam to store the relaxation strategy */
  PLASTICITY_NSGS_RELAXATION = 4,
  /** index in iparam to store the shuffle strategy */
  PLASTICITY_NSGS_SHUFFLE = 5,
  /** index in iparam to store the shuffle seed */
  PLASTICITY_NSGS_SHUFFLE_SEED = 6,
  /** index in iparam to store the  */
  PLASTICITY_NSGS_FREEZING_CONE = 19,
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
  /** PLASTICITY_2D_NSN: coefficient c of the Levenberg-Marquardt damping
      lambda = c ||Phi|| (default 1) */
  PLASTICITY_NSN_LM_COEF = 28,
  /** PLASTICITY_2D_NSN: factor c of rho_e = c ||M_ee^{-1}||_2 (or c/||M_ee||_2, see
      PLASTICITY_NSN_RHO_STRATEGY) for the full Newton steps, default 10 (the merit phase
      uses c = 1) */
  PLASTICITY_NSN_RHO_FACTOR = 29,
};

/** PLASTICITY_2D_NSN: index in iparam of the Levenberg-Marquardt fallback of the natural-map
    merit (1: used when the iteration stalls or the Newton system is singular, default;
    0: never) */
enum PLASTICITY_NSN_GLOBAL_IPARAM {
  PLASTICITY_NSN_LM_FALLBACK = 28,
  /** PLASTICITY_2D_NSN: index in iparam of the globalization strategy
      (PLASTICITY_NSN_GLOBALIZATION_ENUM) */
  PLASTICITY_NSN_GLOBALIZATION = 29
};

/** PLASTICITY_2D_NSN globalization strategies */
enum PLASTICITY_NSN_GLOBALIZATION_ENUM {
  /** full Newton steps; after PLASTICITY_NSN_WATCHDOG_ITERATIONS of them without
      convergence (or a non-finite iterate / singular system), the merit function from the
      best iterate (default) */
  PLASTICITY_NSN_GLOBALIZATION_WATCHDOG = 0,
  /** merit function only: QP objective with projected arc search for convex QP problems,
      1/2 ||Phi||^2 with Armijo and Levenberg-Marquardt otherwise */
  PLASTICITY_NSN_GLOBALIZATION_MERIT = 1,
  /** full Newton steps only (no globalization; returns the best iterate) */
  PLASTICITY_NSN_GLOBALIZATION_NONE = 2
};

/** PLASTICITY_2D_NSN: number of full Newton steps before the watchdog switches to the merit
    function */
#define PLASTICITY_NSN_WATCHDOG_ITERATIONS 50

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
  /** PLASTICITY_2D_NSN: rho_e = ||M_ee^{-1}||_2 = 1/lambda_min(M_ee) for each element
      (default of PLASTICITY_2D_NSN; SPECTRAL_NORM: rho_e = 1/||M_ee||_2) */
  PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_INVERSE_SPECTRAL_NORM = 5,
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
