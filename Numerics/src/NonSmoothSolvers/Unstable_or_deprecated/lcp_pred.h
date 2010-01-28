/* Siconos-Numerics, Copyright INRIA 2005-2010.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

/*!\file lcp_pred.h
  \brief Header for lcp solvers with extract-predict mechanism (Unstable, under devel.)

*/
#ifndef LCP_pred_H
#define LCP_pred_H

#include "LCP_Solvers.h"

/** structure used to define the lcp problem - Deprecated.
*/
typedef struct
{
} method_lcp;


#ifdef __cplusplus
extern "C"
{
#endif

  /** Solver with extract-predict mechanism - Unstable version (under development) */
  int lcp_solver_pred(double *vec, double *q , int *n , method_lcp *pt , double *z , double *w ,
                      int firsttime, int *soltype , int *indic , int *indicop , double *submatlcp , double *submatlcpop ,
                      int *ipiv , int *sizesublcp , int *sizesublcpop ,
                      double *subq , double *bufz , double *newz , double *workspace);

  /** Solver with extract-predict mechanism - Unstable version (under development) */
  int lcp_solver_block_pred_vec(SparseBlockStructuredMatrix *blmat, SparseBlockStructuredMatrixPred *blmatpred, int nbmethod,
                                int maxiterglob, double tolglob,
                                double *q, method_lcp **pt , double *z , double *w , int *it_end , int *itt_end , double *res);

#ifdef __cplusplus
}
#endif

#endif
