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
#ifndef PRSOLVERS_H
#define PRSOLVERS_H

/*!\file prsolvers_h
  \author Nineb Sheherazade and Dubois Frederic.
  Last Modifications : Mathieu Renouf , Pascal Denoyelle, Franck Perignon
 */

/**@defgroup group3 PR (Primal Relay)
   @{
*/

/** \fn int extern  pr_solver ( double* , double* , int* , method* , double* , double* )

 * \brief pr_solver.c is a generic interface allowing the call of one of the @ref pr solvers.
 */

/** @page pr

  The C routines that solve PR:

  pr_latin.c

  pr_nlgs.c

*/

/*!\struct method_pr
   \brief A type definition for a structure method_pr.

  \param name       name of the solver.
  \param itermax    maximum number of iterations.
  \param tol        convergence criteria value.
  \param k_latin    latin coefficient
  \param a          upper bound.
  \param b          lower bound.
  \param chat       output boolean ( 0 = no output log ).
  \param normType   name norm (not yet available).

  \param iter       final number of iterations.
  \param err        final value of error criteria.
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

} method_pr;

#ifdef __cplusplus
extern "C" {
#endif

  void pr_latin(double* , double* , int* , double* , double* , double* , int* , double* , int *, double* , double* , int* , double* , int*);

  void pr_nlgs(double* , double* , int* , double* , double* , int* , double* , int*, double* , double* , int* , double* , int *);

#ifdef __cplusplus
}
#endif

#endif
