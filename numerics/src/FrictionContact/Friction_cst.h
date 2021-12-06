#ifndef FRICTION_CST_H
#define FRICTION_CST_H
/** \file Friction_cst.h */
/** \enum FRICTION_SOLVER encode the list of solvers as integers, to avoid mispelling
 * with const char* const  variables
 */
enum FRICTION_SOLVER
{
  /** 2D Frictional Contact solvers */
  SICONOS_FRICTION_2D_NSGS = 400,
  SICONOS_FRICTION_2D_CPG = 402,
  SICONOS_FRICTION_2D_LEMKE = 404,
  SICONOS_FRICTION_2D_ENUM = 405,

  /* 3D frictional contact solvers on local formulation */

  /** Non-smooth Gauss Seidel, local formulation */
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
  SICONOS_FRICTION_3D_NSN_AC_TEST = 521,
  /** Panagiotopoulos, fixed point, local formulation */
  SICONOS_FRICTION_3D_PFP = 522,
  /** ADMM local formulation */
  SICONOS_FRICTION_3D_ADMM = 523,

  /* 3D Frictional Contact solvers for one contact (used mainly inside NSGS solvers) */

  /** Non-smooth Newton Alart-Curnier, 'direct', one contact solver */
  SICONOS_FRICTION_3D_ONECONTACT_NSN= 550,
  /** Non-smooth Newton Alart-Curnier, 'damped', one contact solver */
  SICONOS_FRICTION_3D_ONECONTACT_NSN_GP = 551,
  /** Projection on cone, one contact solver */
  SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone = 552,
  /** Projection on cone, one contact solver */
  SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration = 553,
  SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization = 554,
  SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithDiagonalization = 555,
  SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone_velocity = 558,

  /** Fischer Burmeister/Path, Glocker formulation, one contact solver */
  SICONOS_FRICTION_3D_NCPGlockerFBPATH = 556,
  /** Newton/Fischer Burmeister, Glocker formulation, one contact solver */
  SICONOS_FRICTION_3D_NCPGlockerFBNewton = 561,
  SICONOS_FRICTION_3D_ONECONTACT_QUARTIC = 562,
  SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU = 563,
  SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinder = 557,
  SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration = 564,
  SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID = 565,
  SICONOS_FRICTION_3D_VI_FPP_Cylinder = 566,
  SICONOS_FRICTION_3D_CONVEXQP_PG_CYLINDER = 567,

  /** 3D Frictional contact local solvers on global formulation */
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


  /** Non-smooth Gauss Seidel, local formulation */
  SICONOS_ROLLING_FRICTION_3D_NSGS = 3000,
  SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnCone= 3001,
  SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration = 3002,
  SICONOS_ROLLING_FRICTION_3D_ADMM = 3003,

  /** Non-smooth Gauss Seidel, local formulation */
  SICONOS_ROLLING_FRICTION_2D_NSGS = 4000,
  SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnCone= 4001,
  SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnConeWithLocalIteration = 4002,

  /** Non-smooth Gauss Seidel, global formulation */
  SICONOS_GLOBAL_ROLLING_FRICTION_3D_NSGS_WR = 5000

};



extern const char* const   SICONOS_FRICTION_2D_NSGS_STR ;
extern const char* const   SICONOS_FRICTION_2D_CPG_STR ;
extern const char* const   SICONOS_FRICTION_2D_LEMKE_STR ;
extern const char* const   SICONOS_FRICTION_2D_ENUM_STR ;

extern const char* const   SICONOS_FRICTION_3D_NSGS_STR ;
extern const char* const   SICONOS_FRICTION_3D_NSGSV_STR ;
extern const char* const   SICONOS_FRICTION_3D_PROX_STR;
extern const char* const   SICONOS_FRICTION_3D_TFP_STR ;
extern const char* const   SICONOS_FRICTION_3D_PFP_STR ;
extern const char* const   SICONOS_FRICTION_3D_NSN_AC_STR ;
extern const char* const   SICONOS_FRICTION_3D_NSN_AC_TEST_STR ;
extern const char* const   SICONOS_FRICTION_3D_NSN_FB_STR ;
extern const char* const   SICONOS_FRICTION_3D_NSN_NM_STR ;
extern const char* const   SICONOS_FRICTION_3D_DSFP_STR ;
extern const char* const   SICONOS_FRICTION_3D_VI_EG_STR ;
extern const char* const   SICONOS_FRICTION_3D_VI_FPP_STR ;
extern const char* const   SICONOS_FRICTION_3D_EG_STR ;
extern const char* const   SICONOS_FRICTION_3D_FPP_STR ;
extern const char* const   SICONOS_FRICTION_3D_HP_STR ;
extern const char* const   SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint_STR;
extern const char* const   SICONOS_FRICTION_3D_ONECONTACT_NSN_STR;
extern const char* const   SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_STR;
extern const char* const   SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID_STR;
extern const char* const   SICONOS_FRICTION_3D_NCPGlockerFBNewton_STR;
extern const char* const   SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithDiagonalization_STR;
extern const char* const   SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone_STR;
extern const char* const   SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration_STR;
extern const char* const   SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization_STR;
extern const char* const   SICONOS_FRICTION_3D_NCPGlockerFBPATH_STR;
extern const char* const   SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinder_STR;
extern const char* const   SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration_STR;
extern const char* const   SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone_velocity_STR;
extern const char* const   SICONOS_FRICTION_3D_CONVEXQP_PG_CYLINDER_STR;
extern const char* const   SICONOS_FRICTION_3D_VI_FPP_Cylinder_STR;
extern const char* const   SICONOS_FRICTION_3D_GAMS_PATH_STR;
extern const char* const   SICONOS_FRICTION_3D_GAMS_PATHVI_STR;
extern const char* const   SICONOS_FRICTION_3D_GAMS_LCP_PATH_STR;
extern const char* const   SICONOS_FRICTION_3D_GAMS_LCP_PATHVI_STR;
extern const char* const   SICONOS_FRICTION_3D_SOCLCP_STR;
extern const char* const   SICONOS_FRICTION_3D_ACLMFP_STR;
extern const char* const   SICONOS_FRICTION_3D_ADMM_STR;
extern const char* const   SICONOS_GLOBAL_FRICTION_3D_NSGS_WR_STR ;
extern const char* const   SICONOS_GLOBAL_FRICTION_3D_NSGSV_WR_STR ;
extern const char* const   SICONOS_GLOBAL_FRICTION_3D_PROX_WR_STR ;
extern const char* const   SICONOS_GLOBAL_FRICTION_3D_DSFP_WR_STR ;
extern const char* const   SICONOS_GLOBAL_FRICTION_3D_TFP_WR_STR ;
extern const char* const   SICONOS_GLOBAL_FRICTION_3D_NSGS_STR ;
extern const char* const   SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR_STR ;
extern const char* const   SICONOS_GLOBAL_FRICTION_3D_NSN_AC_STR;
extern const char* const   SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH_STR;
extern const char* const   SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI_STR;
extern const char* const   SICONOS_GLOBAL_FRICTION_3D_VI_FPP_STR;
extern const char* const   SICONOS_GLOBAL_FRICTION_3D_VI_EG_STR;
extern const char* const   SICONOS_GLOBAL_FRICTION_3D_ACLMFP_STR;
extern const char* const   SICONOS_GLOBAL_FRICTION_3D_ADMM_STR;
extern const char* const   SICONOS_GLOBAL_FRICTION_3D_ADMM_WR_STR;
extern const char* const   SICONOS_GLOBAL_FRICTION_3D_IPM_STR;
extern const char* const   SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_STR ;
extern const char* const   SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU_STR ;


extern const char* const   SICONOS_ROLLING_FRICTION_3D_NSGS_STR ;
extern const char* const   SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnCone_STR;
extern const char* const   SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration_STR;
extern const char* const   SICONOS_ROLLING_FRICTION_3D_ADMM_STR ;

extern const char* const   SICONOS_ROLLING_FRICTION_2D_NSGS_STR ;
extern const char* const   SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnCone_STR;
extern const char* const   SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnConeWithLocalIteration_STR;

extern const char* const   SICONOS_GLOBAL_ROLLING_FRICTION_3D_NSGS_WR_STR ;

enum SICONOS_FRICTION_3D_IPARAM
{
  /** index in iparam to store the error strategy for the internal solver */
  SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY =2,
  /** index in iparam to store the rescaling  */
  SICONOS_FRICTION_3D_IPARAM_RESCALING =3,
  /** index in iparam to store the rescaling  */
  SICONOS_FRICTION_3D_IPARAM_RESCALING_CONE =4,
  /** current contact number (example of use: one contact solvers) **/
  SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER =5,
  /** index in iparam to store the error evaluation method */
  SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION = 7,
  /** index in iparam to store the frequency of error evaluation method */
  SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY = 8,
  SICONOS_FRICTION_3D_NUMBER_OF_CONTACTS = 17,
};

enum SICONOS_FRICTION_INTERNAL_ERROR_STRATEGY
{
  SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE =0,
  SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE =1,
  SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE_N_CONTACT =2
};

enum SICONOS_FRICTION_3D_RESCALING_ENUM
{
  SICONOS_FRICTION_3D_RESCALING_NO =0,
  SICONOS_FRICTION_3D_RESCALING_SCALAR=1,
  SICONOS_FRICTION_3D_RESCALING_BALANCING_M=2,
  SICONOS_FRICTION_3D_RESCALING_BALANCING_MH=3,
  SICONOS_FRICTION_3D_RESCALING_BALANCING_MHHT=4
};

enum SICONOS_FRICTION_3D_RESCALING_CONE_ENUM
{
  SICONOS_FRICTION_3D_RESCALING_CONE_NO =0,
  SICONOS_FRICTION_3D_RESCALING_CONE_YES=1
};

enum SICONOS_FRICTION_3D_DPARAM
{
  /** index in dparam to store the internal solver error ratio*/
  SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO = 2
};


enum SICONOS_FRICTION_3D_NSGS_IPARAM
{
  /** index in iparam to store the relaxation strategy */
  SICONOS_FRICTION_3D_NSGS_RELAXATION=4,
  /** index in iparam to store the shuffle strategy */
  SICONOS_FRICTION_3D_NSGS_SHUFFLE=5,
  /** index in iparam to store the shuffle seed */
  SICONOS_FRICTION_3D_NSGS_SHUFFLE_SEED=6,
  /** index in iparam to store the  */
  SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT =19,
  /** index in iparam to store the  */
  SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION =14,
};
enum SICONOS_FRICTION_3D_NSGS_DPARAM
{
  /** index in dparam to store the relaxation strategy */
  SICONOS_FRICTION_3D_NSGS_RELAXATION_VALUE=8,
};


enum SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_IPARAM
{
  SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_IPARAM_USE_TRIVIAL_SOLUTION=10
};

enum SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION
{
  SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_FALSE=0,
  SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_TRUE=1
};

enum SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ENUM
{
  /** Evaluation of the error with the expensive function fc3d_compute_error **/
  SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL = 0,
  /** Evaluation of the error with the cheap incremental variation **/
  SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT = 1,
  /** Evaluation of the error with the cheap incremental variation but we modify
      the incremental toleranve to reach the requred accuracy **/
  SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL = 2,
  /** Evaluation of the error with the expensive function fc3d_compute_error and
      an adaptive frequency for calling the error function  **/
  SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE =3,
};
enum SICONOS_FRICTION_3D_NSGS_SHUFFLE_ENUM
{
  SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE=0,
  SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE=1,
  SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP=2
};

enum SICONOS_FRICTION_3D_NSGS_RELAXATION_ENUM
{
  SICONOS_FRICTION_3D_NSGS_RELAXATION_FALSE,
  SICONOS_FRICTION_3D_NSGS_RELAXATION_TRUE
};
enum SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_ENUM
{

  SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_FALSE =0,
  SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_TRUE =1
};

enum SICONOS_FRICTION_3D_NSN_IPARAM
{
  /** index in iparam to store the strategy for computing rho */
  SICONOS_FRICTION_3D_NSN_RHO_STRATEGY = 9,
  /** index in iparam to store the formulation */
  SICONOS_FRICTION_3D_NSN_FORMULATION = 10,
  /** index in iparam to store the line-search */
  SICONOS_FRICTION_3D_NSN_LINESEARCH = 11,
  /** index in iparam to store the maximum number of iterations */
  SICONOS_FRICTION_3D_NSN_LINESEARCH_MAX_ITER = 12,
  /** index in iparam to set the linear solver used at each Newton iteration
   cs_lusol or mumps */
  SICONOS_FRICTION_3D_NSN_LINEAR_SOLVER = 13,
  /** index in iparam to store the strategy for the hybrid solver */
  SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY = 14,
  /** index in iparam to store the maximum number of loop for the hybrid solver */
  SICONOS_FRICTION_3D_NSN_HYBRID_MAX_LOOP = 15,
  /** index in iparam to store the maximum number of iterations for the projection solver */
  SICONOS_FRICTION_3D_NSN_HYBRID_MAX_ITER = 16,
  /** index in iparam used to check if memory allocation has already be done (if true/1) or not (if 0/false) for internal work array. */
  SICONOS_FRICTION_3D_NSN_MEMORY_ALLOCATED= 17,
  /** index in iparam to store the boolean to know if allocation of dwork is needed */
  SICONOS_FRICTION_3D_NSN_MPI_COM= 18

};

enum SICONOS_FC3D_NSN_LINEAR_SOLVER
  {
   SICONOS_FRICTION_3D_NSN_USE_CSLUSOL = 0,
   SICONOS_FRICTION_3D_NSN_USE_MUMPS = 1
  };

enum SICONOS_FRICTION_3D_NSN_DPARAM
{
  /** index in dparam to store the rho value for projection formulation */
  SICONOS_FRICTION_3D_NSN_RHO = 3,
};


enum SICONOS_FRICTION_3D_NSN_RHO_STRATEGY_ENUM
{
  /** A constant value given in dparam[SICONOS_FRICTION_3D_NSN_RHO] is used */
  SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_CONSTANT = 0,
  /** A computed value stored in dparam[SICONOS_FRICTION_3D_NSN_RHO] is used */
  SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM =1,
  /** A computed value stored in dparam[SICONOS_FRICTION_3D_NSN_RHO] is used */
  SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM_COND = 2,
  /** An adaptive strategy for rho is used */
  SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM =3,
  /** An adaptive strategy for rho is used */
  SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_ADAPTIVE =4,
};


enum SICONOS_FRICTION_3D_NSN_FORMULATION_ENUM
{
  SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD =0,
  SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD = 1,
  SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_GENERATED =2,
  SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_GENERATED =3,
  SICONOS_FRICTION_3D_NSN_FORMULATION_NULL = 4 ,
};


enum SICONOS_FRICTION_3D_NSN_LINESEARCH_ENUM
{
  SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE = 0 ,
  SICONOS_FRICTION_3D_NSN_LINESEARCH_ARMIJO = 1,
  SICONOS_FRICTION_3D_NSN_LINESEARCH_NO=-1,
};

enum SICONOS_FRICTION_3D_NSN_HYBRID_ENUM
{
  /** No strategy for the hybrid solver */
  SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NO = 0,
  /** Loop PLI-NSN strategy for the hybrid solver */
  SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP = 1,
  /** NSN and after Loop PLI-NSN strategy for the hybrid solver */
  SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NSN_AND_PLI_NSN_LOOP = 2,
  /** VI_EG preconditionning to NSN strategy for the hybrid solver */
  SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_VI_EG_NSN =3
};



enum SICONOS_FRICTION_3D_PROXIMAL_IPARAM
{
  /** index in iparam to store the error strategy for the internal solver */
  SICONOS_FRICTION_3D_FP_ERROR_STRATEGY = 2,
  /** index in iparam to store the relaxation strategy*/
  SICONOS_FRICTION_3D_PROXIMAL_IPARAM_CUMULATIVE_ITER_DONE = 6,
  /** index in iparam to store the relaxation strategy*/
  SICONOS_FRICTION_3D_PROXIMAL_IPARAM_RELAXATION = 8,
  /** index in iparam to store the proximal strategy*/
  SICONOS_FRICTION_3D_PROXIMAL_IPARAM_STRATEGY = 9,
};


enum SICONOS_FRICTION_3D_PROXIMAL_DPARAM
{
  /** index in dparam to store the parameter alpha*/
  SICONOS_FRICTION_3D_PROXIMAL_DPARAM_ALPHA =3,
  SICONOS_FRICTION_3D_PROXIMAL_DPARAM_SIGMA =4,
  SICONOS_FRICTION_3D_PROXIMAL_DPARAM_NU =5,
  SICONOS_FRICTION_3D_PROXIMAL_DPARAM_RELAXATION =8,

};

enum SICONOS_FRICTION_3D_PROXIMAL
{
  /** Proximal algorithm */
  SICONOS_FRICTION_3D_PROXIMAL_PROX = 0,

  /** Regularization algorithm */
  SICONOS_FRICTION_3D_PROXIMAL_REGULARIZATION = 1

};

enum SICONOS_FRICTION_3D_ADMM_IPARAM_ENUM
{
  /** index in iparam to store the strategy for computing rho */
  SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY = 9,
  /** index in iparam to store the strategy for computing rho */
  SICONOS_FRICTION_3D_ADMM_IPARAM_INITIAL_RHO = 10,
  /** index in iparam to store the acceleration parameter */
  SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION= 11,
  /** index in iparam to store the symmetry parameter */
  SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY= 12,
  /** index in iparam to store the sparse storage parameter */
  SICONOS_FRICTION_3D_ADMM_IPARAM_SPARSE_STORAGE= 13,
  /** index in iparam to get problem info */
  SICONOS_FRICTION_3D_ADMM_IPARAM_GET_PROBLEM_INFO= 14,
  SICONOS_FRICTION_3D_ADMM_IPARAM_UPDATE_S= 15,
  SICONOS_FRICTION_3D_ADMM_IPARAM_FULL_H= 17
};

enum SICONOS_FRICTION_3D_ADMM_DPARAM_ENUM
{
  /** index in dparam to store the rho value for projection formulation */
  SICONOS_FRICTION_3D_ADMM_RHO = 3,
  /** index in dparam to store the eta value for the restarting criteria */
  SICONOS_FRICTION_3D_ADMM_RESTART_ETA = 4,
  /** index in dparam to store the tau value for the balancing residual technique */
  SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_TAU = 5,
  /** index in dparam to store the phi value for the balancing residual technique */
  SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_PHI = 6
};

enum SICONOS_FRICTION_3D_ADMM_ACCELERATION_ENUM
{
  SICONOS_FRICTION_3D_ADMM_NO_ACCELERATION= 0,
  SICONOS_FRICTION_3D_ADMM_ACCELERATION= 1,
  SICONOS_FRICTION_3D_ADMM_ACCELERATION_AND_RESTART= 2
};

enum SICONOS_FRICTION_3D_ADMM_SYMMETRY_ENUM
{
  /* default choice. We check symmetry of the problem (Matrix M)
   * if the problem is not symmetric, we called an asymmetric
   * version of the algo is possible */
  SICONOS_FRICTION_3D_ADMM_CHECK_SYMMETRY= 0,
  /* The symmetric version of the algorithm is used even if
   *  the system is not symmetric using the LU solver */
  SICONOS_FRICTION_3D_ADMM_FORCED_SYMMETRY= 1,
  /* The asymmetric version of the algorithm is used even if
   *  the system is symmetric */
  SICONOS_FRICTION_3D_ADMM_FORCED_ASYMMETRY= 2,
  /* The symmetric version of the algorithm is used and the matrix
   *is systematically symmetrized*/
  SICONOS_FRICTION_3D_ADMM_SYMMETRIZE= 3,
  /* The symmetric version of the algorithm is used and we assume
   *  that the data are symmetric */
  SICONOS_FRICTION_3D_ADMM_ASSUME_SYMMETRY= 4
};

enum SICONOS_FRICTION_3D_ADMM_STORAGE_ENUM
{
  SICONOS_FRICTION_3D_ADMM_KEEP_STORAGE= 0,
  SICONOS_FRICTION_3D_ADMM_FORCED_SPARSE_STORAGE= 1
};

enum SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO_ENUM
{
  SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO_NO= 0,
  SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO_YES= 1
};

enum SICONOS_FRICTION_3D_ADMM_UPDATE_S_ENUM
{
  SICONOS_FRICTION_3D_ADMM_UPDATE_S_YES= 0,
  SICONOS_FRICTION_3D_ADMM_UPDATE_S_NO= 1
};

enum SICONOS_FRICTION_3D_ADMM_FULL_H_ENUM
{
  SICONOS_FRICTION_3D_ADMM_FULL_H_NO= 0,
  SICONOS_FRICTION_3D_ADMM_FULL_H_YES= 1
};

enum SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_ENUM
{
  /** A constant value given in dparam[SICONOS_FRICTION_3D_NSN_RHO] is used */
  SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_CONSTANT = 0,
  /** An adaptive strategy for rho is used */
  SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING =2,
  /** An adaptive strategy for rho is used */
  SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING =3
};

enum SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_ENUM
{
  /** A constant value given in dparam[SICONOS_FRICTION_3D_NSN_RHO] is used */
  SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_GIVEN = 0,
  /** A computed value stored in dparam[SICONOS_FRICTION_3D_NSN_RHO] is used */
  SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_NORM_INF =1,
  /** An adaptive strategy for rho is used */
  SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_EIGENVALUES =2
};

enum SICONOS_FRICTION_3D_IPM_IPARAM_ENUM
{
  SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING= 11,
  /** index in iparam to store the sparse storage parameter */
  SICONOS_FRICTION_3D_IPM_IPARAM_SPARSE_STORAGE= 12,
  /** index in iparam to get problem info */
  SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO= 13,
  /** index in iparam to print iterates (including problem data) into a Matlab file */
  SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE= 14,
  /** index in iparam to use reduce the linear system */
  SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM= 15,
  /** index in iparam to finish the solution without scaling */
  SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING= 16,
  /** index in iparam to update the vector w for solving nonconvex problem */
  SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S = 17,
  /** index in iparam to use Qp or F formula for computing Nesterov-Todd scaling 
  SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD = 18
};

enum SICONOS_FRICTION_3D_IPM_DPARAM_ENUM
{
  /** index in dparam to store the parameter for computation the power of sigma */
  SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1 = 7,
  SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2 = 8,
  SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3 = 9,

  /** index in dparam to store the parameter for computation the safity coefficient of step length */
  SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1 = 10,
  SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2 = 11,
};

enum SICONOS_FRICTION_3D_IPM_STORAGE_ENUM
{
  SICONOS_FRICTION_3D_IPM_KEEP_STORAGE= 0,
  SICONOS_FRICTION_3D_IPM_FORCED_SPARSE_STORAGE= 1
};

enum SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_ENUM
{
  SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_NO= 0,
  SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_YES= 1
};

enum SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_METHOD_ENUM
  {
    SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP = 0,
    SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F = 1
  };

#endif
