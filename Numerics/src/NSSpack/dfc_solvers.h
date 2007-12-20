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
#ifndef DFCSOLVERS_H
#define DFCSOLVERS_H

/*!\file dfc_solvers.h
  \author Nineb Sheherazade and Dubois Frederic.
  Last Modifications : Mathieu Renouf , Pascal Denoyelle, Franck Perignon
*/

/**@defgroup group7 2D DFC (Two-Dimensional Dual Frictional Contact)
   @{
*/

/** \fn int extern dfc_2D_solver( double *vec , double *q , int *n , method *pt , double *z , double *w, double *mu )

* \brief dfc_2D_solver() is a generic interface allowing the call of one of the @ref dfc solvers.

*/
/**@page dfc_2D

The C routines that solve DFC:

dfc_2D_latin.c

*/

/*!\struct method_dfc_2D

\brief A type definition for a structure method_dfc_2D

\param name       name of the solver.
\param itermax    maximum number of iterations.
\param normType   name norm (not yet available).
\param tol        convergence criteria value.
\param k_latin    latin coefficient
\param J1         gap in normal contact direction.
\param ddl_n      the contact in normal direction dof (not prescribed),
\param ddl_tt     the contact in tangential direction dof (not prescribed)

\param ddl_d      the prescribed dof.
\param dim_tt     the dimension of the vector ddl_tt.
\param dim_d      the dimension of the vector ddl_d.
\param chat       output boolean ( 0 = no output log ).
\param iter       final number of iteration
\param err        final value of error criteria
*/

typedef struct
{

  char   name[64];
  int    itermax;
  char   normType[64];
  double tol;
  double k_latin;

  double *J1;
  int    *ddl_n;
  int    *ddl_tt;
  int    *ddl_d;
  int    dim_tt;
  int    dim_d;

  int   chat;
  int    iter;
  double err;

} method_dfc_2D;

#ifdef __cplusplus
extern "C" {
#endif

  void dfc_2D_latin(double* , double* , int* , double* , double* , int* , double* , int *, double* , double* , int* , double* , int*);

  void dfc_2D2lcp(int *, double *, double *, double *, int *, int *, int * , int *, int *, double * , double *, double *);

  void lcp2dfc_2D(int *, double *, double *, double *, double *, double *,  int *, int *,
                  int *, int *, int *,  double *, double *);

  void dfc_2Dcond_2D(int *, double *, double *, double *, int *, int *, int * , int *, int *, double * , double *, double *);

  void cond_2D2dfc_2D(int *, double *, double *, double *, double *, double *,  int *, int *,
                      int *, int *, int *,  double *, double *);

#ifdef __cplusplus
}
#endif

#endif
