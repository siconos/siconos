/*!\file solve_cfd.c

  This subroutine allows the dual resolution of contact problems with friction.

  Try \f$(z,w)\f$ such that:
\f$
\left\lbrace
\begin{array}{l}
M z- w=q\\
0 \le z_n \perp w_n \ge 0\\
-z_t \in \partial\psi_{[-\mu w_n, \mu w_n]}(w_t)\\
\end{array}
\right.
\f$

 here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional  vector and w an n-dimensional vector.

 This system of equations and inequalities is solved thanks to @ref dfc solvers or thanks to @ref lcp routines after a new formulation of this problem in the LCP form due to the cfd_lcp.c and lcp_cfd.c routines.
 The routine's call is due to the function solve_cfd.c.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SiconosNumerics.h"


/*!\fn  int solve_cfd (double *vec,double *q,int *n,methode *pt,double z[],double w[])

   solve_cfd is a generic interface allowing the call of one of the DFC solvers.
   \param vec On enter a double vector containing the components of the double matrix with a fortran90 allocation.
   \param q On enter a pointer over doubles containing the components of the double vector.
   \param n On enter a pointer over integers, the dimension of the second member.
   \param pt On enter a pointer over a structure (::methode).
   \param z On return double vector, the solution of the problem.
   \param w On return double vector, the solution of the problem.

   \return On return a pointer over integers, the termination reason (0 is successful otherwise 1).

   \author Nineb Sheherazade.
*/
int solve_cfd(double *vec, double *q, int *n, methode *pt, double z[], double w[])
{
  int info = -1, choix, it_end, fail, i, j;
  const char mot1[10] = "Latin", mot2[10] = "Lemke", mot3[10] = "Gsnl", mot4[10] = "Gcp";
  double res, Mij, *vtel;
  int nn = *n;
  int nc = nn / 2;
  int nc3 = 3 * nc;
  double *Mtel, *Mteltoto;
  double qtel[3 * nc], ztel[3 * nc], wtel[3 * nc], qteltoto[3 * nc];

  vtel = (double*)malloc(3 * nc * 3 * nc * sizeof(double));
  Mtel = (double*)malloc(3 * nc * 3 * nc * sizeof(double));
  Mteltoto = (double*)malloc(3 * nc * 3 * nc * sizeof(double));

  if (strcmp(pt->cfd.nom_method, mot1) == 0)
    cfd_latin(vec, q, &nn, & pt->cfd.k_latin, & pt->cfd.mu, & pt->cfd.itermax, & pt->cfd.tol, z, w, & it_end, &res, &info);
  else if (strcmp(pt->cfd.nom_method, mot2) == 0)
  {
    //    cfd_lcp (&nc,& pt->cfd.mu,vec,q,Mtel,qtel);
    //    //// valeurs du tableau dans vtel (compatibilite allocation memoire f90)///
    //    for (i=0;i<3*nc;i++)
    //      for (j=0;j<3*nc;j++)
    //        vtel[j*3*nc+i]= Mtel[i][j];
    //
    //    lemke_lcp_(vtel,qtel,&nc3,& pt->cfd.itermax,ztel,wtel,&it_end,&res,&info);
    //    lcp_cfd (&nc,ztel,wtel,z,w);
  }
  else if (strcmp(pt->cfd.nom_method, mot3) == 0)
  {
    cfd_lcp(&nc, & pt->cfd.mu, vec, q, Mtel, qtel);

    //// valeurs du tableau dans vtel (compatibilite allocation memoire f90/C)///
    for (i = 0; i < 3 * nc; i++)
      for (j = 0; j < 3 * nc; j++)
        vtel[3 * nc * i + j] = Mtel[i + 3 * nc * j];

    gsnl_lcp(vtel, qtel, &nc3, & pt->cfd.itermax, & pt->cfd.tol, ztel, wtel, &it_end, &res, &info);
    lcp_cfd(&nc, ztel, wtel, z, w);
  }
  else if (strcmp(pt->cfd.nom_method, mot4) == 0)
  {
    cfd_lcp(&nc, & pt->cfd.mu, vec, q, Mtel, qtel);

    //// valeurs du tableau dans vtel (compatibilite allocation memoire f90/C)///
    for (i = 0; i < 3 * nc; i++)
      for (j = 0; j < 3 * nc; j++)
        vtel[3 * nc * i + j] = Mtel[i + 3 * nc * j];

    gcp_lcp(vtel, qtel, &nc3, & pt->cfd.itermax, & pt->cfd.tol, ztel, wtel, &it_end, &res, &info);
    lcp_cfd(&nc, ztel, wtel, z, w);
  }
  else printf("Warning : Unknown solving method : %s\n", pt->cfd.nom_method);

  free(vtel);
  free(Mteltoto);
  free(Mtel);
  return info;
}
