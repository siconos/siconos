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
#ifndef DRSOLVERS_H
#define DRSOLVERS_H

/*!\file dr_solvers.h
  \author Nineb Sheherazade and Dubois Frederic.
  Last Modifications : Mathieu Renouf , Pascal Denoyelle, Franck Perignon
*/

/**@defgroup group4 DR (Dual Relay)
   @{
*/

/** \fn int extern  dr_solver( double* , double* , int* , method* , double* , double* )

* \brief dr_solver.c is a generic interface allowing the call of one of the @ref dr solvers.

*/

/**@page dr

The C routines that solve DR:

dr_latin.c

dr_nlgs.c
*/

/*!\struct method_dr

\brief A type definition for a structure method_dr.

\param name       name of the solver.
\param itermax    maximum number of iterations.
\param tol        convergence criteria value.
\param k_latin    latin coefficient
\param a          upper bound
\param b          lower bound
\param chat       output boolean ( 0 = no output log ).
\param normType   name norm (not yet available).

\param iter       final number of iterations
\param err        final value of error criteria
*/

typedef struct
{

  char     name[64];
  int      itermax;
  double   tol;
  double   k_latin;
  double   *a;
  double   *b;
  int      chat;
  char     normType[64];
  int      iter;
  double   err;

} method_dr;

#ifdef __cplusplus
extern "C" {
#endif

  void dr_latin(double * , double *, int *, double * , double *, double *, int *, double *, int *, double* , double* , int *, double *, int *)  ;

  void dr_nlgs(double *vec , double *q , int *nn , double *a , double *b , int *itermax , double *tol , int *chat,
               double *z , double *w , int *it_end , double *res , int *info);

#ifdef __cplusplus
}
#endif

#endif
