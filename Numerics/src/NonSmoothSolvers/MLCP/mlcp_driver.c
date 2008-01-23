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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#ifndef MEXFLAG
#include "NonSmoothDrivers.h"
#endif


int mlcp_driver(double *A , double *B , double *C , double *D , double *a , double *b, int *n , int* m, method *pt ,  double *u, double *v, double *w)
{

  const char  mlcpkey1[10] = "PGS", mlcpkey2[10] = "RPGS", mlcpkey3[10] = "PSOR" , mlcpkey4[10] = "RPSOR", mlcpkey5[10] = "PATH" ;


  int i, info = 1;

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
  else if (strcmp(pt->mlcp.name , mlcpkey3) == 0)
  {
    iparamMLCP[0] = pt->mlcp.itermax;
    iparamMLCP[1] = pt->mlcp.chat;
    dparamMLCP[0] = pt->mlcp.tol;
    dparamMLCP[1] = pt->mlcp.rho;
    dparamMLCP[2] = pt->mlcp.relax;


    mlcp_psor(n , m, A , B , C , D , a  , b, u, v, w , &info , iparamMLCP , dparamMLCP);

    pt->mlcp.iter = iparamMLCP[2];
    pt->mlcp.err  = dparamMLCP[2];


  }
  else if (strcmp(pt->mlcp.name , mlcpkey4) == 0)
  {
    iparamMLCP[0] = pt->mlcp.itermax;
    iparamMLCP[1] = pt->mlcp.chat;
    dparamMLCP[0] = pt->mlcp.tol;
    dparamMLCP[1] = pt->mlcp.rho;
    dparamMLCP[2] = pt->mlcp.relax;


    /*      mlcp_rpsor( n ,m, A , B , C , D , a  , b, u, v, w , &info , iparamMLCP , dparamMLCP ); */

    pt->mlcp.iter = iparamMLCP[2];
    pt->mlcp.err  = dparamMLCP[2];


  }
  else if (strcmp(pt->mlcp.name , mlcpkey5) == 0)
  {
    iparamMLCP[0] = pt->mlcp.itermax;
    iparamMLCP[1] = pt->mlcp.chat;
    dparamMLCP[0] = pt->mlcp.tol;
    dparamMLCP[1] = pt->mlcp.rho;
    dparamMLCP[2] = pt->mlcp.relax;


    mlcp_path(n , m, A , B , C , D , a  , b, u, v, w , &info , iparamMLCP , dparamMLCP);

    pt->mlcp.iter = iparamMLCP[2];
    pt->mlcp.err  = dparamMLCP[2];


  }
  else printf("Warning : Unknown driver : %s\n", pt->mlcp.name);

  /* Checking validity of z found  */

  /*  if (info == 0) info = filter_result_MLCP(*n,vec,q,z,pt->mlcp.tol,pt->mlcp.chat,w);*/

  //   info= mlcp_filter_result(n,  m, A ,B , C ,D ,a ,b,u, v,  pt->mlcp.tol, pt->mlcp.chat, w);

  return info;

}
