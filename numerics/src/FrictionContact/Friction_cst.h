#ifndef FRICTION_CST_H
#define FRICTION_CST_H
/** \file Friction_cst.h */
/** \enum FRICTION_SOLVER encode the list of solvers as integers, to avoid mispelling
 * with char * variables
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

  /* 3D Frictional Contact solvers for one contact (used mainly inside NSGS solvers) */
  
  /** Non-smooth Newton Alart-Curnier, 'direct', one contact solver */
  SICONOS_FRICTION_3D_ONECONTACT_NSN_AC= 550,
  /** Non-smooth Newton Alart-Curnier, 'damped', one contact solver */
  SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP = 551,
  /** Projection on cone, one contact solver */
  SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone = 552,
  /** Projection on cone, one contact solver */
  SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration = 553,
  SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization = 554,
  SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithDiagonalization = 555,
  SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone_velocity = 558,
  SICONOS_FRICTION_3D_PGoC = 559,
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
  SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP_P = 565,
  
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
  SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI = 609

};



extern char *  SICONOS_FRICTION_2D_NSGS_STR ;
extern char *  SICONOS_FRICTION_2D_PGS_STR ;
extern char *  SICONOS_FRICTION_2D_CPG_STR ;
extern char *  SICONOS_FRICTION_2D_LATIN_STR ;
extern char *  SICONOS_FRICTION_2D_LEMKE_STR ;
extern char *  SICONOS_FRICTION_2D_ENUM_STR ;
extern char *  SICONOS_FRICTION_3D_NSGS_STR ;
extern char *  SICONOS_FRICTION_3D_NSGSV_STR ;
extern char *  SICONOS_FRICTION_3D_PROX_STR;
extern char *  SICONOS_FRICTION_3D_TFP_STR ;
extern char *  SICONOS_FRICTION_3D_NSN_AC_STR ;
extern char *  SICONOS_FRICTION_3D_NSN_FB_STR ;
extern char *  SICONOS_FRICTION_3D_NSN_NM_STR ;
extern char *  SICONOS_FRICTION_3D_DSFP_STR ;
extern char *  SICONOS_FRICTION_3D_VI_EG_STR ;
extern char *  SICONOS_FRICTION_3D_VI_FPP_STR ;
extern char *  SICONOS_FRICTION_3D_EG_STR ;
extern char *  SICONOS_FRICTION_3D_FPP_STR ;
extern char *  SICONOS_FRICTION_3D_HP_STR ;
extern char *  SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint_STR;
extern char *  SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_STR;
extern char *  SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP_STR;
extern char *  SICONOS_FRICTION_3D_NCPGlockerFBNewton_STR;
extern char *  SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithDiagonalization_STR;
extern char *  SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone_STR;
extern char *  SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration_STR;
extern char *  SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization_STR;
extern char *  SICONOS_FRICTION_3D_NCPGlockerFBPATH_STR;
extern char *  SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinder_STR;
extern char *  SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone_velocity_STR;
extern char *  SICONOS_FRICTION_3D_PGoC_STR;
extern char *  SICONOS_FRICTION_3D_DeSaxceFixedPoint_STR;
extern char *  SICONOS_FRICTION_3D_GAMS_PATH_STR;
extern char *  SICONOS_FRICTION_3D_GAMS_PATHVI_STR;
extern char *  SICONOS_FRICTION_3D_GAMS_LCP_PATH_STR;
extern char *  SICONOS_FRICTION_3D_GAMS_LCP_PATHVI_STR;
extern char *  SICONOS_GLOBAL_FRICTION_3D_NSGS_WR_STR ;
extern char *  SICONOS_GLOBAL_FRICTION_3D_NSGSV_WR_STR ;
extern char *  SICONOS_GLOBAL_FRICTION_3D_PROX_WR_STR ;
extern char *  SICONOS_GLOBAL_FRICTION_3D_DSFP_WR_STR ;
extern char *  SICONOS_GLOBAL_FRICTION_3D_TFP_WR_STR ;
extern char *  SICONOS_GLOBAL_FRICTION_3D_NSGS_STR ;
extern char *  SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR_STR ;
extern char *  SICONOS_GLOBAL_FRICTION_3D_NSN_AC_STR;
extern char *  SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH_STR;
extern char *  SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI_STR;
extern char *  SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_STR ;
extern char *  SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU_STR ;

#endif
