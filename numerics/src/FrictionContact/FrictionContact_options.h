/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

/*!\file FrictionContact_options.h
 * \brief Solver options and constants for Friction Contact problems
 *
 * This header defines:
 * - Solver IDs (enum FRICTION_SOLVER): All available solvers for friction contact problems
 * - Integer parameter indices (iparam): Configuration options stored as integers
 * - Double parameter indices (dparam): Configuration options stored as doubles
 * - Enumerated values: Valid values for various configuration options
 *
 * The friction contact solvers are organized into several categories:
 * - 2D Frictional Contact solvers
 * - 3D Frictional Contact solvers (local formulation)
 * - One-contact solvers (used as local solvers within NSGS)
 * - Global formulation solvers
 * - Rolling friction solvers
 *
 * NSGS-specific options have been moved to NonSmoothSolvers/NonSmoothGaussSeidel_options.h
 * but backward compatibility macros are provided.
 */

#ifndef FRICTION_CST_H
#define FRICTION_CST_H

/* ===========================================================================
 * Solver IDs (enum FRICTION_SOLVER)
 * ===========================================================================
 * Each solver is assigned a unique integer ID to avoid string-based errors.
 * The IDs are organized by solver category.
 */
/** \enum FRICTION_SOLVER
 * \brief Unique identifiers for all friction contact solvers
 *
 * Solvers are grouped by problem formulation (2D, 3D, global, rolling)
 * and solver family (NSGS, Proximal, Newton, etc.)
 */
enum FRICTION_SOLVER {
  /* -----------------------------------------------------------------------
   * 2D Frictional Contact solvers (IDs 400-499)
   * ----------------------------------------------------------------------- */
  /** Non-smooth Gauss-Seidel for 2D friction contact */
  SICONOS_FRICTION_2D_NSGS = 400,
  SICONOS_FRICTION_2D_CPG = 402,
  SICONOS_FRICTION_2D_LEMKE = 404,
  /** Enumerative solver for 2D friction contact */
  SICONOS_FRICTION_2D_ENUM = 405,

  /* -----------------------------------------------------------------------
   * 3D Frictional Contact solvers - Local formulation (IDs 500-549)
   * -----------------------------------------------------------------------
   * These solvers work on the local (contact-level) formulation of the
   * 3D friction contact problem. They typically use a block decomposition
   * where each contact is solved individually (e.g., NSGS).
   */
  SICONOS_FRICTION_3D_NSGS = 500,
  /** Non-smooth Gauss Seidel-velocity, local formulation */
  SICONOS_FRICTION_3D_NSGSV = 501,
  /** proximal, local formulation */
  SICONOS_FRICTION_3D_PROX = 502,
  /** Tresca, fixed point, local formulation */
  SICONOS_FRICTION_3D_TFP = 503,
  /** Non-smooth Newton Alart-Curnier, local formulation */
  SICONOS_FRICTION_3D_NSN_AC = 504,
  /** De Saxce fixed point, local formulation */
  SICONOS_FRICTION_3D_DSFP = 505,
  /** VI formulation, fixed point projection, local formulation */
  SICONOS_FRICTION_3D_VI_FPP = 506,
  /** VI formulation, Extra-gradient, local formulation */
  SICONOS_FRICTION_3D_VI_EG = 507,
  /** Hyperplane projection, local formulation */
  SICONOS_FRICTION_3D_HP = 508,
  /** Fischer Burmeister fixed point, local formulation */
  SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint = 510,
  /** Fixed point projection, local formulation */
  SICONOS_FRICTION_3D_FPP = 511,
  /** Extra-gradient, local formulation */
  SICONOS_FRICTION_3D_EG = 512,
  /** Non-smooth Newton Fischer Burmeister, local formulation */
  SICONOS_FRICTION_3D_NSN_FB = 513,
  /** GAMS/Path (Ferris), local formulation */
  SICONOS_FRICTION_3D_GAMS_PATH = 514,
  /** VI formulation, GAMS/Path (Ferris), local formulation */
  SICONOS_FRICTION_3D_GAMS_PATHVI = 515,
  /** Alart-Curnier fixed point, local formulation */
  SICONOS_FRICTION_3D_ACLMFP = 516,
  /** Second-order Cone LCP, local formulation */
  SICONOS_FRICTION_3D_SOCLCP = 517,
  /** GAMS/PATH (Ferris) LCP, local formulation */
  SICONOS_FRICTION_3D_GAMS_LCP_PATH = 518,
  /** VI formulation, GAMS/PATH (Ferris) LCP, local formulation */
  SICONOS_FRICTION_3D_GAMS_LCP_PATHVI = 519,
  /** Non-smooth Newton, natural map, local formulation */
  SICONOS_FRICTION_3D_NSN_NM = 520,
  /** Non-smooth Newton, Alart-Curnier, new implementation */
  SICONOS_FRICTION_3D_NSN_AC_NEW = 521,
  /** Panagiotopoulos fixed point, local formulation */
  SICONOS_FRICTION_3D_PFP = 522,
  /** ADMM local formulation */
  SICONOS_FRICTION_3D_ADMM = 523,
  /** IPM local formulation */
  SICONOS_FRICTION_3D_IPM = 524,
  /** IPM-SNM local formulation */
  SICONOS_FRICTION_3D_IPM_SNM = 525,

  /* -----------------------------------------------------------------------
   * One-contact solvers (IDs 550-599)
   * -----------------------------------------------------------------------
   * These solvers are designed to solve a single contact problem.
   * They are primarily used as local solvers within NSGS-type methods.
   */

  /** Non-smooth Newton, Alart-Curnier 'direct', one contact solver */
  SICONOS_ONECONE_NSN = 550,
  /** Non-smooth Newton Alart-Curnier, 'damped', one contact solver */
  SICONOS_ONECONE_NSN_GP = 551,
  /** Projection on cone, one contact solver */
  SICONOS_ONECONE_ProjectionOnCone = 552,
  /** Projection on cone, one contact solver */
  SICONOS_ONECONE_ProjectionOnConeWithLocalIteration = 553,
  SICONOS_ONECONE_ProjectionOnConeWithRegularization = 554,
  SICONOS_ONECONE_ProjectionOnConeWithDiagonalization = 555,
  SICONOS_ONECONE_ProjectionOnCone_velocity = 558,

  /** Fischer Burmeister/Path, Glocker formulation, one contact solver */
  SICONOS_FRICTION_3D_NCPGlockerFBPATH = 556,
  /** Newton/Fischer Burmeister, Glocker formulation, one contact solver */
  SICONOS_FRICTION_3D_NCPGlockerFBNewton = 561,
  SICONOS_FRICTION_3D_ONECONTACT_QUARTIC = 562,
  SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU = 563,
  SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinder = 557,
  SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration = 564,
  SICONOS_ONECONE_NSN_GP_HYBRID = 565,
  SICONOS_FRICTION_3D_VI_FPP_Cylinder = 566,
  SICONOS_FRICTION_3D_CONVEXQP_PG_CYLINDER = 567,

  /* -----------------------------------------------------------------------
   * Global formulation solvers (IDs 600-699)
   * -----------------------------------------------------------------------
   * These solvers work on the global (aggregated) formulation of the
   * friction contact problem, treating all contacts simultaneously.
   * The "_WR" suffix indicates wrappers for reduced-space formulations.
   */

  /** Non-smooth Gauss-Seidel, global formulation with wrapper */
  SICONOS_GLOBAL_FRICTION_3D_NSGS_WR = 600,
  SICONOS_GLOBAL_FRICTION_3D_NSGSV_WR = 601,
  SICONOS_GLOBAL_FRICTION_3D_PROX_WR = 602,
  SICONOS_GLOBAL_FRICTION_3D_DSFP_WR = 603,
  SICONOS_GLOBAL_FRICTION_3D_TFP_WR = 604,
  SICONOS_GLOBAL_FRICTION_3D_NSGS = 605,
  SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR = 606,
  SICONOS_GLOBAL_FRICTION_3D_NSN_AC = 607,
  SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH = 608,
  SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI = 609,
  /** VI formulation, Fixed Point Projection, local formulation */
  SICONOS_GLOBAL_FRICTION_3D_VI_FPP = 610,
  /** VI formulation, Extra-gradient, local formulation */
  SICONOS_GLOBAL_FRICTION_3D_VI_EG = 611,
  SICONOS_GLOBAL_FRICTION_3D_ACLMFP = 612,
  SICONOS_GLOBAL_FRICTION_3D_ADMM = 613,
  SICONOS_GLOBAL_FRICTION_3D_ADMM_WR = 614,
  SICONOS_GLOBAL_FRICTION_3D_IPM = 615,
  SICONOS_GLOBAL_FRICTION_3D_IPM_WR = 616,
  SICONOS_GLOBAL_FRICTION_3D_IPM_SNM = 617,
  SICONOS_GLOBAL_FRICTION_3D_IPM_SNM_WR = 618,
  /** IPM-SNM with proximal regularization */
  SICONOS_GLOBAL_FRICTION_3D_IPM_SNM_PROX = 620,

}; /* end of enum FRICTION_SOLVER */

/* ===========================================================================
 * Rolling Friction Solver Options
 * ===========================================================================
 * Rolling friction solvers are defined in RollingFrictionContact_options.h
 * Include it for rolling friction specific solver IDs and options.
 */
#include "RollingFrictionContact/RollingFrictionContact_options.h"

/* ===========================================================================
 * Solver Name Strings
 * ===========================================================================
 * These are human-readable string representations of each solver ID,
 * used for logging, debugging, and output purposes.
 */

/** String name for 2D NSGS solver */
extern const char* const SICONOS_FRICTION_2D_NSGS_STR;
extern const char* const SICONOS_FRICTION_2D_CPG_STR;
extern const char* const SICONOS_FRICTION_2D_LEMKE_STR;
extern const char* const SICONOS_FRICTION_2D_ENUM_STR;

extern const char* const SICONOS_FRICTION_3D_NSGS_STR;
extern const char* const SICONOS_FRICTION_3D_NSGSV_STR;
extern const char* const SICONOS_FRICTION_3D_PROX_STR;
extern const char* const SICONOS_FRICTION_3D_TFP_STR;
extern const char* const SICONOS_FRICTION_3D_PFP_STR;
extern const char* const SICONOS_FRICTION_3D_NSN_AC_STR;
extern const char* const SICONOS_FRICTION_3D_NSN_AC_NEW_STR;
extern const char* const SICONOS_FRICTION_3D_NSN_FB_STR;
extern const char* const SICONOS_FRICTION_3D_NSN_NM_STR;
extern const char* const SICONOS_FRICTION_3D_DSFP_STR;
extern const char* const SICONOS_FRICTION_3D_VI_EG_STR;
extern const char* const SICONOS_FRICTION_3D_VI_FPP_STR;
extern const char* const SICONOS_FRICTION_3D_EG_STR;
extern const char* const SICONOS_FRICTION_3D_FPP_STR;
extern const char* const SICONOS_FRICTION_3D_HP_STR;
extern const char* const SICONOS_FRICTION_3D_IPM_STR;
extern const char* const SICONOS_FRICTION_3D_IPM_SNM_STR;
extern const char* const SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint_STR;
extern const char* const SICONOS_ONECONE_NSN_STR;
extern const char* const SICONOS_ONECONE_NSN_GP_STR;
extern const char* const SICONOS_ONECONE_NSN_GP_HYBRID_STR;
extern const char* const SICONOS_FRICTION_3D_NCPGlockerFBNewton_STR;
extern const char* const
    SICONOS_ONECONE_ProjectionOnConeWithDiagonalization_STR;
extern const char* const SICONOS_ONECONE_ProjectionOnCone_STR;
extern const char* const SICONOS_ONECONE_ProjectionOnConeWithLocalIteration_STR;
extern const char* const SICONOS_ONECONE_ProjectionOnConeWithRegularization_STR;
extern const char* const SICONOS_FRICTION_3D_NCPGlockerFBPATH_STR;
extern const char* const SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinder_STR;
extern const char* const
    SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration_STR;
extern const char* const SICONOS_ONECONE_ProjectionOnCone_velocity_STR;
extern const char* const SICONOS_FRICTION_3D_CONVEXQP_PG_CYLINDER_STR;
extern const char* const SICONOS_FRICTION_3D_VI_FPP_Cylinder_STR;
extern const char* const SICONOS_FRICTION_3D_GAMS_PATH_STR;
extern const char* const SICONOS_FRICTION_3D_GAMS_PATHVI_STR;
extern const char* const SICONOS_FRICTION_3D_GAMS_LCP_PATH_STR;
extern const char* const SICONOS_FRICTION_3D_GAMS_LCP_PATHVI_STR;
extern const char* const SICONOS_FRICTION_3D_SOCLCP_STR;
extern const char* const SICONOS_FRICTION_3D_ACLMFP_STR;
extern const char* const SICONOS_FRICTION_3D_ADMM_STR;
extern const char* const SICONOS_GLOBAL_FRICTION_3D_IPM_STR;
extern const char* const SICONOS_GLOBAL_FRICTION_3D_IPM_WR_STR;
extern const char* const SICONOS_GLOBAL_FRICTION_3D_IPM_SNM_STR;
extern const char* const SICONOS_GLOBAL_FRICTION_3D_IPM_SNM_WR_STR;
extern const char* const SICONOS_GLOBAL_FRICTION_3D_IPM_SNM_PROX_STR;
extern const char* const SICONOS_GLOBAL_FRICTION_3D_NSGS_WR_STR;
extern const char* const SICONOS_GLOBAL_FRICTION_3D_NSGSV_WR_STR;
extern const char* const SICONOS_GLOBAL_FRICTION_3D_PROX_WR_STR;
extern const char* const SICONOS_GLOBAL_FRICTION_3D_DSFP_WR_STR;
extern const char* const SICONOS_GLOBAL_FRICTION_3D_TFP_WR_STR;
extern const char* const SICONOS_GLOBAL_FRICTION_3D_NSGS_STR;
extern const char* const SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR_STR;
extern const char* const SICONOS_GLOBAL_FRICTION_3D_NSN_AC_STR;
extern const char* const SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH_STR;
extern const char* const SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI_STR;
extern const char* const SICONOS_GLOBAL_FRICTION_3D_VI_FPP_STR;
extern const char* const SICONOS_GLOBAL_FRICTION_3D_VI_EG_STR;
extern const char* const SICONOS_GLOBAL_FRICTION_3D_ACLMFP_STR;
extern const char* const SICONOS_GLOBAL_FRICTION_3D_ADMM_STR;
extern const char* const SICONOS_GLOBAL_FRICTION_3D_ADMM_WR_STR;
extern const char* const SICONOS_GLOBAL_FRICTION_3D_IPM_STR;

extern const char* const SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_STR;
extern const char* const SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU_STR;

/* Rolling friction solver strings are declared in RollingFrictionContact_options.h */

/* ===========================================================================
 * Integer Parameter Indices (iparam)
 * ===========================================================================
 * These enums define indices into the iparam array of SolverOptions.
 * Use these to configure integer-valued solver options.
 *
 * Example:
 *   options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] =
 *       SICONOS_NSGS_ERROR_EVALUATION_LIGHT;
 */

/** Indices for general friction contact 3D solver parameters */
enum SICONOS_FRICTION_3D_IPARAM {

  /** Error strategy for the internal/local solver (see SICONOS_NSGS_INTERNAL_ERROR_STRATEGY) */
  SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY = 2,
  /** Rescaling strategy for the problem (see SICONOS_FRICTION_3D_RESCALING_ENUM) */
  SICONOS_FRICTION_3D_IPARAM_RESCALING = 3,
  /** Rescaling strategy for friction cones (see SICONOS_FRICTION_3D_RESCALING_CONE_ENUM) */
  SICONOS_FRICTION_3D_IPARAM_RESCALING_CONE = 4,
  /** Current contact/block number (used by one-contact solvers) */
  SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER = 5,
  /** Error evaluation method (see SICONOS_NSGS_ERROR_EVALUATION) - now in NonSmoothGaussSeidel_options.h */
  SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION = 7,
  /** Frequency of error evaluation (0 = every iteration) */
  SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY = 8,
  /** Number of contacts in the problem */
  SICONOS_FRICTION_3D_NUMBER_OF_CONTACTS = 17,
};

/** Matrix rescaling strategies for friction contact problems */
enum SICONOS_FRICTION_3D_RESCALING {
  /** No rescaling applied */
  SICONOS_FRICTION_3D_RESCALING_NO = 0,
  /** Scalar rescaling (simple scaling factor) */
  SICONOS_FRICTION_3D_RESCALING_SCALAR = 1,
  /** Matrix balancing on M */
  SICONOS_FRICTION_3D_RESCALING_BALANCING_M = 2,
  /** Matrix balancing on M and H */
  SICONOS_FRICTION_3D_RESCALING_BALANCING_MH = 3,
  /** Matrix balancing on M, H, and H^T */
  SICONOS_FRICTION_3D_RESCALING_BALANCING_MHHT = 4
};

/** Friction cone rescaling options */
enum SICONOS_FRICTION_3D_RESCALING_CONE {
  /** No cone rescaling */
  SICONOS_FRICTION_3D_RESCALING_CONE_NO = 0,
  /** Apply cone rescaling */
  SICONOS_FRICTION_3D_RESCALING_CONE_YES = 1
};

/* ===========================================================================
 * Double Parameter Indices (dparam)
 * ===========================================================================
 * These enums define indices into the dparam array of SolverOptions.
 * Use these to configure floating-point solver options like tolerances
 * and numerical parameters.
 */

/** Indices for friction contact 3D double parameters */
enum SICONOS_FRICTION_3D_DPARAM {
  /** Ratio for internal solver error tolerance (relative to outer tolerance) */
  SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO = 2
};

/* ===========================================================================
 * NSGS Options
 * ===========================================================================
 * NSGS-specific options are now defined in NonSmoothGaussSeidel_options.h.
 * The include below provides backward compatibility macros so existing code
 * using the old FRICTION_3D_NSGS_* names continues to work.
 */
#include "NonSmoothSolvers/NonSmoothGaussSeidel_options.h"

/* ===========================================================================
 * Non-Smooth Newton (NSN) Solver Options
 * ===========================================================================
 * Parameters for Alart-Curnier, Fischer-Burmeister, and Natural Map
 * Newton-type solvers. These are nonsmooth Newton methods for solving
 * the friction contact problem via complementarity formulations.
 */

/** Indices for NSN solver integer parameters */
enum SICONOS_FRICTION_3D_NSN_IPARAM {
  /** Strategy for computing the rho parameter (see SICONOS_FRICTION_3D_NSN_RHO_STRATEGY_ENUM) */
  SICONOS_FRICTION_3D_NSN_RHO_STRATEGY = 9,
  /** Complementarity formulation (see SICONOS_FRICTION_3D_NSN_FORMULATION_ENUM) */
  SICONOS_FRICTION_3D_NSN_FORMULATION = 10,
  /** Line-search strategy (see SICONOS_FRICTION_3D_NSN_LINESEARCH_ENUM) */
  SICONOS_FRICTION_3D_NSN_LINESEARCH = 11,
  /** Maximum iterations for line-search */
  SICONOS_FRICTION_3D_NSN_LINESEARCH_MAX_ITER = 12,
  /** Linear solver for Newton systems (see SICONOS_FC3D_NSN_LINEAR_SOLVER) */
  SICONOS_FRICTION_3D_NSN_LINEAR_SOLVER = 13,
  /** Hybrid solver strategy (see SICONOS_FRICTION_3D_NSN_HYBRID_ENUM) */
  SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY = 14,
  /** Maximum outer loops for hybrid solver */
  SICONOS_FRICTION_3D_NSN_HYBRID_MAX_LOOP = 15,
  /** Maximum iterations for projection solver in hybrid mode */
  SICONOS_FRICTION_3D_NSN_HYBRID_MAX_ITER = 16,
  /** Flag: internal work array allocated (1=yes, 0=no) */
  SICONOS_FRICTION_3D_NSN_MEMORY_ALLOCATED = 17,
  /** MPI communicator for distributed solving */
  SICONOS_FRICTION_3D_NSN_MPI_COM = 18

};

enum SICONOS_FC3D_NSN_LINEAR_SOLVER {
  SICONOS_FRICTION_3D_NSN_USE_CSLUSOL = 0,
  SICONOS_FRICTION_3D_NSN_USE_MUMPS = 1
};

enum SICONOS_FRICTION_3D_NSN_DPARAM {
  /** index in dparam to store the rho value for projection formulation */
  SICONOS_FRICTION_3D_NSN_RHO = 3,
};

enum SICONOS_FRICTION_3D_NSN_RHO_STRATEGY {
  /** A constant value given in dparam[SICONOS_FRICTION_3D_NSN_RHO] is used */
  SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_CONSTANT = 0,
  /** A computed value stored in dparam[SICONOS_FRICTION_3D_NSN_RHO] is used */
  SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM = 1,
  /** A computed value stored in dparam[SICONOS_FRICTION_3D_NSN_RHO] is used */
  SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM_COND = 2,
  /** An adaptive strategy for rho is used */
  SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM = 3,
  /** An adaptive strategy for rho is used */
  SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_ADAPTIVE = 4,
};

enum SICONOS_FRICTION_3D_NSN_FORMULATION {
  SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD = 0,
  SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD = 1,
  SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_GENERATED = 2,
  SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_GENERATED = 3,
  SICONOS_FRICTION_3D_NSN_FORMULATION_NATURALMAP = 4,
  SICONOS_FRICTION_3D_NSN_FORMULATION_NULL = 5,
};

enum SICONOS_FRICTION_3D_NSN_LINESEARCH {
  SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE = 0,
  SICONOS_FRICTION_3D_NSN_LINESEARCH_ARMIJO = 1,
  SICONOS_FRICTION_3D_NSN_LINESEARCH_NO = -1,
};

enum SICONOS_FRICTION_3D_NSN_HYBRID {
  /** No strategy for the hybrid solver */
  SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NO = 0,
  /** Loop PLI-NSN strategy for the hybrid solver */
  SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP = 1,
  /** NSN and after Loop PLI-NSN strategy for the hybrid solver */
  SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NSN_AND_PLI_NSN_LOOP = 2,
  /** VI_EG preconditionning to NSN strategy for the hybrid solver */
  SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_VI_EG_NSN = 3
};

/* ===========================================================================
 * Proximal Solver Options
 * ===========================================================================
 * Parameters for the proximal point algorithm and regularization methods
 * for friction contact problems.
 */

/** Indices for Proximal solver integer parameters */
enum SICONOS_FRICTION_3D_PROXIMAL_IPARAM {
  /** index in iparam to store the error strategy for the internal solver */
  SICONOS_FRICTION_3D_FP_ERROR_STRATEGY = 2,
  /** index in iparam to store the relaxation strategy*/
  SICONOS_FRICTION_3D_PROXIMAL_IPARAM_CUMULATIVE_ITER_DONE = 6,
  /** index in iparam to store the relaxation strategy*/
  SICONOS_FRICTION_3D_PROXIMAL_IPARAM_RELAXATION = 8,
  /** index in iparam to store the proximal strategy*/
  SICONOS_FRICTION_3D_PROXIMAL_IPARAM_STRATEGY = 9,
};

enum SICONOS_FRICTION_3D_PROXIMAL_DPARAM {
  /** index in dparam to store the parameter alpha*/
  SICONOS_FRICTION_3D_PROXIMAL_DPARAM_ALPHA = 3,
  SICONOS_FRICTION_3D_PROXIMAL_DPARAM_SIGMA = 4,
  SICONOS_FRICTION_3D_PROXIMAL_DPARAM_NU = 5,
  SICONOS_FRICTION_3D_PROXIMAL_DPARAM_RELAXATION = 8,

};

enum SICONOS_FRICTION_3D_PROXIMAL {
  /** Proximal algorithm */
  SICONOS_FRICTION_3D_PROXIMAL_PROX = 0,

  /** Regularization algorithm */
  SICONOS_FRICTION_3D_PROXIMAL_REGULARIZATION = 1

};

/* ===========================================================================
 * ADMM Solver Options (Alternating Direction Method of Multipliers)
 * ===========================================================================
 * Parameters for the ADMM solver, a first-order optimization method that
 * splits the problem into simpler subproblems.
 */

/** Indices for ADMM solver integer parameters */
enum SICONOS_FRICTION_3D_ADMM_IPARAM {
  /** index in iparam to store the strategy for computing rho */
  SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY = 9,
  /** index in iparam to store the strategy for computing rho */
  SICONOS_FRICTION_3D_ADMM_IPARAM_INITIAL_RHO = 10,
  /** index in iparam to store the acceleration parameter */
  SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION = 11,
  /** index in iparam to store the symmetry parameter */
  SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY = 12,
  /** index in iparam to store the sparse storage parameter */
  SICONOS_FRICTION_3D_ADMM_IPARAM_SPARSE_STORAGE = 13,
  /** index in iparam to get problem info */
  SICONOS_FRICTION_3D_ADMM_IPARAM_GET_PROBLEM_INFO = 14,
  SICONOS_FRICTION_3D_ADMM_IPARAM_UPDATE_S = 15,
  SICONOS_FRICTION_3D_ADMM_IPARAM_FULL_H = 17,
  /** index in iparam to print on screen the output in the same style as in IPM solver */
  SICONOS_FRICTION_3D_ADMM_PRINTING_LIKE_IPM = 18
};

enum SICONOS_FRICTION_3D_ADMM_DPARAM {
  /** index in dparam to store the rho value for projection formulation */
  SICONOS_FRICTION_3D_ADMM_RHO = 3,
  /** index in dparam to store the eta value for the restarting criteria */
  SICONOS_FRICTION_3D_ADMM_RESTART_ETA = 4,
  /** index in dparam to store the tau value for the balancing residual technique */
  SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_TAU = 5,
  /** index in dparam to store the phi value for the balancing residual technique */
  SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_PHI = 6
};

enum SICONOS_FRICTION_3D_ADMM_ACCELERATION {
  SICONOS_FRICTION_3D_ADMM_NO_ACCELERATION = 0,
  SICONOS_FRICTION_3D_ADMM_ACCELERATION = 1,
  SICONOS_FRICTION_3D_ADMM_ACCELERATION_AND_RESTART = 2
};

enum SICONOS_FRICTION_3D_ADMM_SYMMETRY {
  /* default choice. We check symmetry of the problem (Matrix M)
   * if the problem is not symmetric, we called an asymmetric
   * version of the algo is possible */
  SICONOS_FRICTION_3D_ADMM_CHECK_SYMMETRY = 0,
  /* The symmetric version of the algorithm is used even if
   *  the system is not symmetric using the LU solver */
  SICONOS_FRICTION_3D_ADMM_FORCED_SYMMETRY = 1,
  /* The asymmetric version of the algorithm is used even if
   *  the system is symmetric */
  SICONOS_FRICTION_3D_ADMM_FORCED_ASYMMETRY = 2,
  /* The symmetric version of the algorithm is used and the matrix
   *is systematically symmetrized*/
  SICONOS_FRICTION_3D_ADMM_SYMMETRIZE = 3,
  /* The symmetric version of the algorithm is used and we assume
   *  that the data are symmetric */
  SICONOS_FRICTION_3D_ADMM_ASSUME_SYMMETRY = 4
};

enum SICONOS_FRICTION_3D_ADMM_STORAGE {
  SICONOS_FRICTION_3D_ADMM_KEEP_STORAGE = 0,
  SICONOS_FRICTION_3D_ADMM_FORCED_SPARSE_STORAGE = 1
};

enum SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO {
  SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO_NO = 0,
  SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO_YES = 1
};

enum SICONOS_FRICTION_3D_ADMM_UPDATE_S {
  SICONOS_FRICTION_3D_ADMM_UPDATE_S_YES = 0,
  SICONOS_FRICTION_3D_ADMM_UPDATE_S_NO = 1
};

enum SICONOS_FRICTION_3D_ADMM_FULL_H {
  SICONOS_FRICTION_3D_ADMM_FULL_H_NO = 0,
  SICONOS_FRICTION_3D_ADMM_FULL_H_YES = 1
};

enum SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY {
  /** A constant value given in dparam[SICONOS_FRICTION_3D_NSN_RHO] is used */
  SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_CONSTANT = 0,
  /** An adaptive strategy for rho is used */
  SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING = 2,
  /** An adaptive strategy for rho is used */
  SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING = 3
};

enum SICONOS_FRICTION_3D_ADMM_INITIAL_RHO {
  /** A constant value given in dparam[SICONOS_FRICTION_3D_NSN_RHO] is used */
  SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_GIVEN = 0,
  /** A computed value stored in dparam[SICONOS_FRICTION_3D_NSN_RHO] is used */
  SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_NORM_INF = 1,
  /** An adaptive strategy for rho is used */
  SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_EIGENVALUES = 2
};

enum SICONOS_FRICTION_3D_ADMM_PRINTING_LIKE_IPM {
  SICONOS_FRICTION_3D_ADMM_PRINTING_LIKE_IPM_FALSE = 0,
  SICONOS_FRICTION_3D_ADMM_PRINTING_LIKE_IPM_TRUE = 1
};

/* ===========================================================================
 * IPM Solver Options (Interior Point Method)
 * ===========================================================================
 * Parameters for the interior point method, a second-order optimization
 * approach that follows the central path to the solution.
 */

/** Indices for IPM solver integer parameters */
enum SICONOS_FRICTION_3D_IPM_IPARAM {
  /** index in iparam to use NT scaling technique */
  SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING = 11,
  /** index in iparam to store the sparse storage parameter */
  SICONOS_FRICTION_3D_IPM_IPARAM_SPARSE_STORAGE = 12,
  /** index in iparam to get problem info */
  SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO = 13,
  /** index in iparam to print iterates (including problem data) into a Matlab file */
  SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE = 14,
  /** index in iparam to use reduce the linear system */
  SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM = 15,
  /** index in iparam to finish the solution without scaling */
  SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING = 16,
  /** index in iparam to update the vector w for solving nonconvex problem */
  SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S = 17,
  /** index in iparam to use Qp or F formula for computing Nesterov-Todd scaling **/
  SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD = 18,
  /** index in iparam to use a reduced symmetric system with Qp2 or QpH inside the matrix **/
  SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM_METHOD = 19,
  /** index in iparam to choose the linear system formulation **/
  SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM = 10,
  /** index in iparam to perform iterative refinement with MA57 (LDT) **/
  SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT = 9,
  /** index in iparam to perform Cholesky factorization **/
  SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY = 8,
  /** index in iparam to print iterates (including problem data) into a Python file */
  SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_PYTHON_FILE = 7,
  /** index in iparam to perform Mehrotra scheme **/
  SICONOS_FRICTION_3D_IPM_IPARAM_MEHROTRA = 6
};

enum SICONOS_FRICTION_3D_IPM_DPARAM {
  /** index in dparam to store the parameter for computation the power of sigma */
  SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1 = 7,
  SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2 = 8,
  SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3 = 9,

  /** index in dparam to store the parameter for computation the safity coefficient of step
     length */
  SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1 = 10,
  SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2 = 11,
};

enum SICONOS_FRICTION_3D_IPM_STORAGE {
  SICONOS_FRICTION_3D_IPM_KEEP_STORAGE = 0,
  SICONOS_FRICTION_3D_IPM_FORCED_SPARSE_STORAGE = 1
};

enum SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO {
  SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_NO = 0,
  SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_YES = 1
};

enum SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_METHOD {
  SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP = 0,
  SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F = 1
};

enum SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM_METHOD {
  SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM_WITH_QP2 = 0,
  SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM_WITH_QPH = 1
};
enum SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM {
  SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL = 0,
  SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QP2 = 1,
  SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QPH = 2,
  SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QP2 = 3,
  SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QPH = 4,
  SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_QPH = 5,
  SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_JQJ = 6,
  SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_JQJ = 7,
  SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_invPH = 8,
  SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_JQinv = 9,
  SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_NOSCAL = 10,
  SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_QP2 = 11,
};

enum SICONOS_FRICTION_3D_IPM_IPARAM_MEHROTRA {
  SICONOS_FRICTION_3D_IPM_IPARAM_MEHROTRA_NO = 0,
  SICONOS_FRICTION_3D_IPM_IPARAM_MEHROTRA_YES = 1
};

enum SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING {
  SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING_NO = 0,
  SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING_YES = 1
};

enum SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT {
  SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_NO = 0,
  SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES = 1,
  SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_AFTER = 2
};

enum SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY {
  SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_NO = 0,
  SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_YES = 1
};

enum SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S {
  SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S_EXTERNAL = 1,
  SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S_INTERNAL = 2,
};

/* Backward compatibility macros for renamed enums */
#define SICONOS_FRICTION_3D_RESCALING_ENUM SICONOS_FRICTION_3D_RESCALING
#define SICONOS_FRICTION_3D_RESCALING_CONE_ENUM SICONOS_FRICTION_3D_RESCALING_CONE
#define SICONOS_FRICTION_3D_NSN_RHO_STRATEGY_ENUM SICONOS_FRICTION_3D_NSN_RHO_STRATEGY
#define SICONOS_FRICTION_3D_NSN_FORMULATION_ENUM SICONOS_FRICTION_3D_NSN_FORMULATION
#define SICONOS_FRICTION_3D_NSN_LINESEARCH_ENUM SICONOS_FRICTION_3D_NSN_LINESEARCH
#define SICONOS_FRICTION_3D_NSN_HYBRID_ENUM SICONOS_FRICTION_3D_NSN_HYBRID
#define SICONOS_FRICTION_3D_ADMM_IPARAM_ENUM SICONOS_FRICTION_3D_ADMM_IPARAM
#define SICONOS_FRICTION_3D_ADMM_DPARAM_ENUM SICONOS_FRICTION_3D_ADMM_DPARAM
#define SICONOS_FRICTION_3D_ADMM_ACCELERATION_ENUM SICONOS_FRICTION_3D_ADMM_ACCELERATION
#define SICONOS_FRICTION_3D_ADMM_SYMMETRY_ENUM SICONOS_FRICTION_3D_ADMM_SYMMETRY
#define SICONOS_FRICTION_3D_ADMM_STORAGE_ENUM SICONOS_FRICTION_3D_ADMM_STORAGE
#define SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO_ENUM SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO
#define SICONOS_FRICTION_3D_ADMM_UPDATE_S_ENUM SICONOS_FRICTION_3D_ADMM_UPDATE_S
#define SICONOS_FRICTION_3D_ADMM_FULL_H_ENUM SICONOS_FRICTION_3D_ADMM_FULL_H
#define SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_ENUM SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY
#define SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_ENUM SICONOS_FRICTION_3D_ADMM_INITIAL_RHO
#define SICONOS_FRICTION_3D_ADMM_PRINTING_LIKE_IPM_ENUM SICONOS_FRICTION_3D_ADMM_PRINTING_LIKE_IPM
#define SICONOS_FRICTION_3D_IPM_IPARAM_ENUM SICONOS_FRICTION_3D_IPM_IPARAM
#define SICONOS_FRICTION_3D_IPM_DPARAM_ENUM SICONOS_FRICTION_3D_IPM_DPARAM
#define SICONOS_FRICTION_3D_IPM_STORAGE_ENUM SICONOS_FRICTION_3D_IPM_STORAGE
#define SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_ENUM SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO
#define SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_METHOD_ENUM SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_METHOD
#define SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM_METHOD_ENUM SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM_METHOD
#define SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM_ENUM SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM
#define SICONOS_FRICTION_3D_IPM_IPARAM_MEHROTRA_ENUM SICONOS_FRICTION_3D_IPM_IPARAM_MEHROTRA
#define SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING_ENUM SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING
#define SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_ENUM SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT
#define SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY_ENUM SICONOS_FRICTION_3D_IPM_IPARAM_CHOLESKY
#define SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S_ENUM SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S

#endif