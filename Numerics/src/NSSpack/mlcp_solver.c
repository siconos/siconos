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
/*!\file mlcp_solver.c

  This subroutine allows the resolution of MLCP (Mixed Linear Complementary Problem).\n
  Try \f$(u,v,w)\f$ such that:\n

  \f$
   \left\lbrace
    \begin{array}{l}
    A u + Cv +a =0\\
    D u + Bv +b = w
    0 \le v \perp  w \ge 0\\
    \end{array}
   \right.
  \f$

  A is an (\f$ n \times n\f$ ) matrix, B is an (\f$ m \times m\f$ ) matrix,  C is an (\f$ n \times m\f$ ) matrix, D is an (\f$ m \times n\f$ ) matrix,    a and u is an (\f$ n \f$ ) vectors b,v and w is an (\f$ m \f$ ) vectors.
  This system of equalities and inequalities
  is solved thanks to @ref mlcp solvers. The routine's call is due to the function mlcp_solver.c.

 !\fn int mlcp_solver(double *A , double *B , double *C , double *D , double *a  double *b, int *n , int* m, method *pt ,  double *u, double *v, double *w  );


  mlcp_solver is a generic interface allowing the call of one of the MLCP solvers.

  \param A            On enter, a (\f$n \times n\f$)-vector of doubles which contains the components of the "A" MLCP matrix with a Fortran storage.
  \param B            On enter, a (\f$m \times m\f$)-vector of doubles which contains the components of the "B" MLCP matrix with a Fortran storage.
  \param C            On enter, a (\f$n \times m\f$)-vector of doubles which contains the components of the "C" MLCP matrix with a Fortran storage.
  \param D            On enter, a (\f$m \times n\f$)-vector of doubles which contains the components of the "D" MLCP matrix with a Fortran storage.
  \param a            On enter, a n-vector of doubles which contains the components of the constant right hand side vector.
  \param b            On enter, a m-vector of doubles which contains the components of the constant right hand side vector.
  \param n           On enter, an integer which represents one the dimension of the MLCP problem.
  \param m           On enter, an integer which represents one the dimension of the MLCP problem.
  \param pt           On enter, a union containing the MLCP structure.
  \n \n
  \param u            On return, a n-vector of doubles which contains the solution of the problem.
  \param v            On return, a m-vector of doubles which contains the solution of the problem.
  \param w            On return, a m-vector of doubles which contains the complementary solution of the problem.



  \return integer
                   - 0 : successful\n
                   - >0 : otherwise (see specific solvers for more information about the log info)

  \author V. Acary
  \todo Sizing the regularization paramter and apply it only on null diagnal term
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#ifndef MEXFLAG
#include "NSSpack.h"
#endif


int mlcp_solver(double *A , double *B , double *C , double *D , double *a , double *b, int *n , int* m, method *pt ,  double *u, double *v, double *w)
{

  const char  mlcpkey1[10] = "PGS", mlcpkey2[10] = "RPGS" ;


  int i, j, info = 1;

  int     iparamMLCP[5];
  double  dparamMLCP[5];

  for (i = 0 ; i < 5 ; ++i) iparamMLCP[i] = 0;
  for (i = 0 ; i < 5 ; ++i) dparamMLCP[i] = 0.0;


  if (strcmp(pt->mlcp.name , mlcpkey1) == 0)
  {

    iparamMLCP[0] = pt->mlcp.itermax;
    iparamMLCP[1] = pt->mlcp.chat;
    dparamMLCP[0] = pt->mlcp.tol;
    /* dparamMLCP[1] = pt->mlcp.relax;*/


    mlcp_pgs(n , m, A , B , C , D , a  , b, u, v, w , &info , iparamMLCP , dparamMLCP);

    pt->mlcp.iter = iparamMLCP[2];
    pt->mlcp.err  = dparamMLCP[2];

  }
  else if (strcmp(pt->mlcp.name , mlcpkey2) == 0)
  {

    iparamMLCP[0] = pt->mlcp.itermax;
    iparamMLCP[1] = pt->mlcp.chat;
    dparamMLCP[0] = pt->mlcp.tol;
    dparamMLCP[1] = pt->mlcp.rho;
    /* dparamMLCP[1] = pt->mlcp.relax;*/


    mlcp_rpgs(n , m, A , B , C , D , a  , b, u, v, w , &info , iparamMLCP , dparamMLCP);

    pt->mlcp.iter = iparamMLCP[2];
    pt->mlcp.err  = dparamMLCP[2];

  }

  else printf("Warning : Unknown solver : %s\n", pt->mlcp.name);

  /* Checking validity of z found  */

  /*  if (info == 0) info = filter_result_MLCP(*n,vec,q,z,pt->mlcp.tol,pt->mlcp.chat,w);*/

  //   info= mlcp_filter_result(n,  m, A ,B , C ,D ,a ,b,u, v,  pt->mlcp.tol, pt->mlcp.chat, w);

  return info;

}
