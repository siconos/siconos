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
#ifndef PFC_2D_SOLVERS_H
#define PFC_2D_SOLVERS_H

/*!\file pfc_2D_solvers.h
  \author Nineb Sheherazade and Dubois Frederic.
  Last Modifications : Mathieu Renouf , Pascal Denoyelle, Franck Perignon
*/

/**@defgroup group5 2D PFC (Two-dimensional Primal Frictional Contact)
   @{
*/

/** \fn int extern  pfc_2D_solver( double *vec , double *q , int *n , method *pt , double *z , double *w )

* \brief pfc_2D_solver.c is a generic interface allowing the call of one of the @ref pfc_2D solvers.

*/

/**@page pfc_2D

The C routines that solve PFC:

pfc_2D_latin.c

pfc_2D_nlgs.c

pfc_2D_cpg.c
*/

/*!\struct method_pfc_2D

\brief A type definition for a structure method_pfc_2D.

\param name       name of the solver.
\param itermax    maximum number of iterations.
\param tol        convergence criteria value.
\param k_latin    search direction of the latin metod.
\param chat       output boolean ( 0 = no output log ).
\param normType   name norm (not yet available).
\param iter       final number of iterations.
\param err        final value of error criteria
*/

typedef struct
{

  char   name[64];
  int    itermax;
  double tol;
  double k_latin;
  int    chat;
  char   normType[64];
  int    iter;
  double err;

} method_pfc_2D;

#ifdef __cplusplus
extern "C" {
#endif

  void pfc_2D_cpg(int, double *vec , double *q , double *z , double *w , double *mu , int *info , int *iparamLCP , double *dparamLCP);

  void pfc_2D_nlgs(int, double *vec , double *q , double *z , double *w, double *mu  , int *info , int *iparamLCP , double *dparamLCP);

  void pfc_2D_latin(int, double *vec , double *q , double *z , double *w , double *mu , int *info , int *iparamLCP , double *dparamLCP);

  /*void pfc_2D_projc( int nc , double mu , double *z , double *p , int *status );
    void pfc_2D_projf( int n , double *ww , double *zz , double *rr , double *pp , int *status )
  */

#ifdef __cplusplus
}
#endif

#endif
