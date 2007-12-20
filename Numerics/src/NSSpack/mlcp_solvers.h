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
#ifndef MLCP_SOLVERS_H
#define MLCP_SOLVERS_H

/*!\file mlcp_solvers.h
  \author Vincent Acary
*/

/**@defgroup group1 MLCP (Mixed Linear Complementary Problem)
   @{
*/

/** \fn int extern mlcp_solver( double *vec , double *q ,int *n , method *pt , double *z , double *w , int *it_end , double *res )
 *  \brief mlcp_solver.c is a generic interface allowing the call of one of the @ref mlcp solvers.
 */
/**@}*/

/**@page mlcp

The C routines that solve MLCP:

mlcp_pgs.c

mlcp_rpgs.c
*/

/**@defgroup group2 Block MLCP (Mixed Linear Complementary Problem)
   @{
*/

/*!\struct method_mlcp

\brief A type definition for a structure method_mlcp.

\param name       name of the solver.
\param itermax    maximum number of iterations.
\param tol        convergence criteria value.
\param relax      relaxation coefficient.
\param rho        regularization coefficient
\param chat       output boolean ( 0 = no output log ).


\param iter       final number of iterations.
\param err        final value of error criteria.
*/

typedef struct
{

  char   name[64];
  int    itermax;
  double tol;
  double relax;
  double rho;
  int    chat;
  int    iter;
  double err;

} method_mlcp;

#ifdef __cplusplus
extern "C" {
#endif

  void mlcp_pgs(int* nn , int* mm, double *A , double *B , double *C , double *D , double *a , double *b, double *u, double *v, double *w , int *info , int *iparamMLCP , double *dparamMLCP);

  void mlcp_rpgs(int* nn , int* mm, double *A , double *B , double *C , double *D , double *a , double *b, double *u, double *v, double *w , int *info , int *iparamMLCP , double *dparamMLCP);

  void mlcp_psor(int* nn , int* mm, double *A , double *B , double *C , double *D , double *a , double *b, double *u, double *v, double *w , int *info , int *iparamMLCP , double *dparamMLCP);

  void mlcp_rpsor(int* nn , int* mm, double *A , double *B , double *C , double *D , double *a , double *b, double *u, double *v, double *w , int *info , int *iparamMLCP , double *dparamMLCP);

  void mlcp_path(int* nn , int* mm, double *A , double *B , double *C , double *D , double *a , double *b, double *u, double *v, double *w , int *info , int *iparamMLCP , double *dparamMLCP);

  int mlcp_filter_result(int* n, int* mm, double *A , double *B , double *C , double *D , double *a , double *b, double *u, double *v,  double tol, int chat, double *w);

  int mlcp_compute_error(int* n, int* mm,  double *A , double *B , double *C , double *D , double *a , double *b, double *u, double *v,  int chat, double *w, double * error);

#ifdef __cplusplus
}
#endif

#endif
