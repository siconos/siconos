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
#ifndef LCP_SOLVERS_H
#define LCP_SOLVERS_H

/*!\file lcp_solvers.h
  \author Nineb Sheherazade and Dubois Frederic.
  Last Modifications : Mathieu Renouf , Pascal Denoyelle, Franck Perignon
*/

/**@defgroup group1 LCP (Linear Complementary Problem)
   @{
*/

/** \fn int extern lcp_solver( double *vec , double *q ,int *n , method *pt , double *z , double *w , int *it_end , double *res )
 *  \brief lcp_solver.c is a generic interface allowing the call of one of the @ref lcp solvers.
 */
/**@}*/

/**@page lcp

The C routines that solve LCP:



lcp_pgs.c

lcp_rpgs.c

lcp_latin.c

lcp_latin_w.c

lcp_cpg.c

lcp_lexicolemke.c

lcp_qp.c

lcp_nsqp.c

lcp_newton_min.c
*/

/**@defgroup group2 Block LCP (Linear Complementary Problem)
   @{
*/

/*!\struct method_lcp

\brief A type definition for a structure method_lcp.

\param name       name of the solver.
\param itermax    maximum number of iterations.
\param tol        convergence criteria value.
\param k_latin    latin coefficient.
\param relax      relaxation coefficient.
\param rho        regularization coefficient
\param chat       output boolean ( 0 = no output log ).
\param normType   name norm (not yet available).

\param iter       final number of iterations.
\param err        final value of error criteria.
*/

typedef struct
{

  char   name[64];
  int    itermax;
  double tol;
  double k_latin;
  double relax;
  double rho;
  int    chat;
  char   normType[64];
  int    iter;
  double err;

} method_lcp;

#ifdef __cplusplus
extern "C" {
#endif

  void lcp_qp(int *nn , double *M , double *q , double *z , double *w , int *info ,
              int *iparamLCP , double *dparamLCP);

  void lcp_cpg(int *nn , double *M , double *q , double *z , double *w , int *info ,
               int *iparamLCP , double *dparamLCP);

  void lcp_pgs(int *nn , double *M , double *q , double *z , double *w , int *info ,
               int *iparamLCP , double *dparamLCP);

  void lcp_rpgs(int *nn , double *M , double *q , double *z , double *w , int *info ,
                int *iparamLCP , double *dparamLCP);

  void lcp_psor(int *nn , double *M , double *q , double *z , double *w , int *info ,
                int *iparamLCP , double *dparamLCP);

  void lcp_nsqp(int *nn , double *M , double *q , double *z , double *w , int *info ,
                int *iparamLCP , double *dparamLCP);

  void lcp_latin(int *nn , double *M , double *q , double *z , double *w , int *info ,
                 int *iparamLCP , double *dparamLCP);

  void lcp_latin_w(int *nn , double *M , double *q , double *z , double *w , int *info ,
                   int *iparamLCP , double *dparamLCP);

  void lcp_lexicolemke(int *nn , double *M , double *q , double *z , double *w , int *info ,
                       int *iparamLCP , double *dparamLCP);

  void lcp_newton_min(int *nn , double *M , double *q , double *z , double *w , int *info ,
                      int *iparamLCP , double *dparamLCP);

  void lcp_newton_FB(int *nn , double *M , double *q , double *z , double *w , int *info ,
                     int *iparamLCP , double *dparamLCP);

  int filter_result_LCP(int n, double *M , double *q , double *z , double tol, int chat, double *w);

  int lcp_compute_error(int n, double *M , double *q , double *z , int chat, double *w, double * error);

#ifdef __cplusplus
}
#endif

#endif
