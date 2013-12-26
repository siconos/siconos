#ifndef VI_CST_H
#define VI_CST_H
/** \file VI_cst.h */
/** \enum VI_SOLVER VI_cst.h
 * Enum that allows one to encode the list of solvers in a proper to avoid mispelling
 * with char * variables
 */
enum VI_SOLVER
{
  SICONOS_VI_EG = 1000,
  SICONOS_VI_FPP = 1001
};

extern char *  SICONOS_VI_EG_STR ;
extern char *  SICONOS_VI_FPP_STR ;

#endif
