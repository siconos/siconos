
/*!\file lcp_cfd.c

   This file allows to give the solution of the contact problem with friction given.

*/



/*!\fn lcp_cfd (int *dim_nn,double *ztel,double *wtel,double *z,double *w)

   lcp_cfd subroutine allows to give the solution of the contact problem with friction given.
    \sa cfd_lcp subroutine.

   \param dim_nn On enter a pointer over integers, the number of normal variables after the condensation of a contact friction problem.
   \param ztel On enter a pointer over doubles, the solution given by an LCP ztel(3*dim_nn).
   \param wtel On enter a pointer over doubles, the solution given by an LCP wtel(3*dim_nn).
   \param z On return a pointer over doubles, the solution of the contact friction problem z(2*dim_nn).
   \param w On return a pointer over doubles, the solution of the contact friction problem w(3*dim_n).

   \author Nineb Sheherazade.


*/
void lcp_cfd(int *dim_nn, double *ztel, double *wtel, double *z, double *w)
{
  int dim_n = *dim_nn, i, j;
  double a1[dim_n], a2[dim_n], b1[dim_n], b2[dim_n], Ut[dim_n], Ft[dim_n];


  for (i = 0; i < dim_n; i++)
  {
    z[i] = ztel[i];
    b1[i] = ztel[dim_n + i];
    b2[i] = ztel[2 * dim_n + i];
    z[dim_n + i] = b2[i] - b1[i];
    w[i] = wtel[i];
    a1[i] = wtel[dim_n + i];
    a2[i] = wtel[2 * dim_n + i];
    w[dim_n + i] = 0.5 * (a2[i] - a1[i]);
  }
}
