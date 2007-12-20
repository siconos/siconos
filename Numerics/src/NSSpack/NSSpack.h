/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#ifndef NSSPACK_H
#define NSSPACK_H

#include "blaslapack.h"
#include "lcp_solvers.h"
#include "mlcp_solvers.h"
#include "pfc_3D_solvers.h"
#include "pfc_2D_solvers.h"
#include "dfc_solvers.h"
#include "dr_solvers.h"
#include "pr_solvers.h"


/*!\file solverpack.h
  \author Nineb Sheherazade and Dubois Frederic.
  Last Modifications : Mathieu Renouf , Pascal Denoyelle

  \todo solve_qp does not exist

*/

/*!\union method

  \brief A type definition for a union method.

  \param method_pr     : pr is a method_pr structure .
  \param method_dr     : dr is a method_dr structure .
  \param method_lcp    : lcp is a method_lcp structure .
  \param method_mlcp    : mlcp is a method_mlcp structure .
  \param method_pfc_2D : pfc_2D is a method_pfc_2D structure .
  \param method_dfc_2D : dfc_2D is a method_dfc_2D structure .
  \param method_qp     : qp is a method_qp structure (not yet available).


*/
typedef union
{

  method_pr  pr;
  method_dr  dr;
  method_lcp lcp;
  method_mlcp mlcp;
  method_pfc_2D pfc_2D;
  method_pfc_3D pfc_3D;
  method_dfc_2D dfc_2D;

} method;

/*
 * header for C++ compiling / and C compiling
 */

#ifdef __cplusplus
extern "C" {
#endif

  /** General interface to solver for lcp problems */
  int lcp_solver(double *M, double *q , int *n , method *pt , double *z , double *w);

  /** General interface to solver for lcp problems, with M given as a list of blocks */
  int lcp_solver_block(SparseBlockStructuredMatrix *blmat, double *q, method *pt , double *z , double *w , int *it_end ,
                       int *itt_end , double *res);

  /** General interface to solver for lcp problems */
  int mlcp_solver(double *A , double *B , double *C , double *D , double *a, double *b, int *n , int* m, method *pt ,  double *u, double *v, double *w);

  /** General interface to solver for pfc 3D problems */
  int pfc_3D_solver(int, double*, double*, method*, double*, double*, double*);

  /** General interface to solvers for primal friction 3D problem, with a sparse block storage for M. */
  int pfc_3D_solver_block(int, SparseBlockStructuredMatrix*, double*, method*, double*, double*, double*);

  /** General interface to solver for pfc 2D problems */
  int pfc_2D_solver(int, double*, double*, method*, double*, double*, double*);

  /** General interface to solver for dual-relay problems */
  int dr_solver(double* , double* , int* , method* , double* , double*);

  /** General interface to solver for primal-relay problems */
  int pr_solver(double* , double* , int* , method* , double* , double*);

  /** General interface to solver for dual friction-contact problems */
  int dfc_2D_solver(double* , double* , int* , method* , double* , double*, double*);

#ifdef __cplusplus
}
#endif

#endif /* NSSPACK_H */
