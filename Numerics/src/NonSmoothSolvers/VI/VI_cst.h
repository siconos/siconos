#ifndef FRICTION_CST_H
#define FRICTION_CST_H
/** \file VI_cst.h */
/** \enum VI_SOLVER Friction_cst.h
 * Enum that allows one to encode the list of solvers in a proper to avoid mispelling
 * with char * variables
 */
enum VI_SOLVER
{
  SICONOS_VI_EG = 1000
};

extern char *  SICONOS_VI_EG_STR ;

#endif
