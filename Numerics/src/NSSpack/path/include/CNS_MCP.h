/*****************************************************************************/
/* CNS_MCP.h                                                                 */
/*                                                                           */
/* DESCRIPTION                                                               */
/*   The constrained nonlinear system is to find an x_star such that:        */
/*     F(x_star) = 0                                                         */
/*     lb <= x_start <= ub                                                   */
/*   where F is a given function from R^n to R^n, and lb and ub are          */
/*   prescribed lower and upper bounds.                                      */
/*                                                                           */
/*   Instead of directly solving this problem, we can convert it into an     */
/*   equivalent mixed complementarity problem to find and x and y such that  */
/*        -y  perp l <= x <= u                                               */
/*       F(x) perp      y                                                    */
/*                                                                           */
/*   This header file defines such a conversion mechanism.  The code itself  */
/*   only adds extra variables for the x's that have at least one finite     */
/*   bound.                                                                  */
/*****************************************************************************/

#ifndef CNS_MCP_H
#define CNS_MCP_H

#include "Types.h"
#include "CNS_Interface.h"
#include "MCP_Interface.h"

/*****************************************************************************/
/* Converson routines.                                                       */
/*****************************************************************************/
/*                                                                           */
/* - CNStoMCP   - reformulate the specified CNS as a complementarity         */
/*                problem.  The new complementarity problem structure is     */
/*                returned.  NOTE: the CNS structure be created and an       */
/*                interface set before calling the conversion routine.       */
/*                                                                           */
/* - CNSfromMCP - takes as input the CNS and a pointed to the MCP structure  */
/*                allocated by the CNStoMCP routine.  This function takes    */
/*                the solution contained in the MCP and communicates is to   */
/*                the CNS.  The MCP is destroyed on exit.                    */
/*                                                                           */
/*****************************************************************************/

FUN_DECL(CNS_Termination) CNStoMCP(CNS *c, MCP **m, Int *m_size, Int *m_nnz);
FUN_DECL(Void)            CNSfromMCP(CNS *c, MCP **m);

#endif
