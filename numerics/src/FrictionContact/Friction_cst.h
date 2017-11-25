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
  SICONOS_FRICTION_2D_PGS = 401,
  SICONOS_FRICTION_2D_CPG = 402,
  SICONOS_FRICTION_2D_LATIN = 403,
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

  /** De Saxce fixed point, one contact solver */
  SICONOS_FRICTION_3D_DeSaxceFixedPoint = 560,
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
  SICONOS_FRICTION_3D_ConvexQP_PG_Cylinder = 567,

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
  SICONOS_GLOBAL_FRICTION_3D_VI_EG = 611
};



extern const char* const   SICONOS_FRICTION_2D_NSGS_STR ;
extern const char* const   SICONOS_FRICTION_2D_PGS_STR ;
extern const char* const   SICONOS_FRICTION_2D_CPG_STR ;
extern const char* const   SICONOS_FRICTION_2D_LATIN_STR ;
extern const char* const   SICONOS_FRICTION_2D_LEMKE_STR ;
extern const char* const   SICONOS_FRICTION_2D_ENUM_STR ;
extern const char* const   SICONOS_FRICTION_3D_NSGS_STR ;
extern const char* const   SICONOS_FRICTION_3D_NSGSV_STR ;
extern const char* const   SICONOS_FRICTION_3D_PROX_STR;
extern const char* const   SICONOS_FRICTION_3D_TFP_STR ;
extern const char* const   SICONOS_FRICTION_3D_PFP_STR ;
extern const char* const   SICONOS_FRICTION_3D_NSN_AC_STR ;
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
extern const char* const   SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone_velocity_STR;
extern const char* const   SICONOS_FRICTION_3D_ConvexQP_PG_Cylinder_STR;
extern const char* const   SICONOS_FRICTION_3D_VI_FPP_Cylinder_STR;
extern const char* const   SICONOS_FRICTION_3D_DeSaxceFixedPoint_STR;
extern const char* const   SICONOS_FRICTION_3D_GAMS_PATH_STR;
extern const char* const   SICONOS_FRICTION_3D_GAMS_PATHVI_STR;
extern const char* const   SICONOS_FRICTION_3D_GAMS_LCP_PATH_STR;
extern const char* const   SICONOS_FRICTION_3D_GAMS_LCP_PATHVI_STR;
extern const char* const   SICONOS_FRICTION_3D_SOCLCP_STR;
extern const char* const   SICONOS_FRICTION_3D_ACLMFP_STR;
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
extern const char* const   SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_STR ;
extern const char* const   SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU_STR ;


enum SICONOS_FRICTION_3D_IPARAM
{
  /** index in iparam to store the error strategy for the internal solver */
  SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY =2
};

enum SICONOS_FRICTION_INTERNAL_ERROR_STRATEGY
{
  SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE =0,
  SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE =1,
  SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE_N_CONTACT =2
};

enum SICONOS_FRICTION_3D_DPARAM
{
  /** index in iparam to store the internal solver error ratio*/
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
  /** index in iparam to store the error evaluation method */
  SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION = 7,
  /** index in iparam to store the frequency of error evaluation method */
  SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FREQUENCY = 8,
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
  SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_CONTACTNUMBER = 4,
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
      an adaptive frequncy for calling the error function  **/
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
  SICONOS_FRICTION_3D_NSN_LINESEARCH_MAXITER = 12,
  /** index in iparam to store the strategy for the hybrid solver */
  SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY = 14,
  /** index in iparam to store the maximum number of loop for the hybrid solver */
  SICONOS_FRICTION_3D_NSN_HYBRID_MAX_LOOP = 15,
  /** index in iparam to store the maximum number of iterations for the projection solver */
  SICONOS_FRICTION_3D_NSN_HYBRID_MAX_ITER = 16
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
  /** Loop NSN-PLI strategy for the hybrid solver */
  SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NO = 0,
  /** Loop NSN-PLI strategy for the hybrid solver */
  SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP = 1,
  /** Loop PLI-NSN-PLI strategy for the hybrid solver */
  SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NSN_PLI_LOOP = 2,
  /** NSN and after Loop NSN-PLI strategy for the hybrid solver */
  SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NSN_AND_NSN_PLI_LOOP = 3,
  SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_VI_EG_NSN =4,
};



enum SICONOS_FRICTION_3D_PROXIMAL_IPARAM
{
  /** index in iparam to store the error strategy for the internal solver */
  SICONOS_FRICTION_3D_FP_ERROR_STRATEGY =2,
  /** index in iparam to store the relaxation strategy*/
  SICONOS_FRICTION_3D_PROXIMAL_IPARAM_CUMULATIVE_ITER_DONE =6,
  /** index in iparam to store the relaxation strategy*/
  SICONOS_FRICTION_3D_PROXIMAL_IPARAM_RELAXATION =8,
  /** index in iparam to store the proximal strategy*/
  SICONOS_FRICTION_3D_PROXIMAL_IPARAM_STRATEGY =9,
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






#endif
