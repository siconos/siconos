#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>



/*!\file lemke_lcp.c


   This subroutine allows the resolution of LCP (Linear Complementary Problem).
   Try \f$(z,w)\f$ such that:

\f$
\left\lbrace
\begin{array}{l}
M z- w=q\\
0 \le z \perp w \ge 0\\
\end{array}
\right.
\f$

  here M is an n by n  matrix, q an n-dimensional vector, w an n-dimensional  vector and z an n-dimensional vector.
*/

double ddot_(int*, double [], int*, double[], int*);
/*!\fn  lemke_lcp(double vec[], double *qq,int *nn, int *itermax, double *z, double *w, int *it_end, double *res, int *info )


   lemke_lcp is a direct solver for LCP.


   \param vec On enter a double vector containing the components of the double matrix with a fortran90 allocation.
   \param qq On enter a pointer over doubles containing the components of the double vector.
   \param nn On enter a pointer over integers, the dimension of the second member.
   \param itermax On enter a pointer over integers, the maximum iterations required.
   \param it_end On enter a pointer over integers, the number of iterations carried out.
   \param res On return a pointer over doubles, the error value.
   \param z On return double vector, the solution of the problem.
   \param w On return double vector, the solution of the problem.
   \param info On return a pointer over integers, the termination reason (0 is successful otherwise 1).
   \author Nineb Sheherazade.
*/

lemke_lcp(double *vec, double *qqq, int *nn, int *itermax, double *zlem, double *wlem, int *it_end, double *res, int *info)
{


  int n = *nn, i, j, k, pos_z0, pos_pivot, pos_q, incx, incy, pos,  itt = *itermax, s[n], *ss, kk, test_out;
  int *basic, iter, pos_out0, pos_out, *posA, dim_posA, dimposAA, *posAA, nonbasic, compt;
  double  alpha, beta, qs, z0, *tempo, *rtmp, *sol, pdt, pivot;
  double M[n][n], val, aa, bb, cc, WI[n][n], W[n][n], MM[n][n], q[n], qq[n];
  double z[n], w[n];
  double A[n][2 * n + 2], E[n][2 * n + 2], EE[n][2 * n + 2], AA[n][2 * n + 2], num, den, resi, minw, minz, *Mz, m0, mm0;
  char trans = 'T';
  int trouve;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      M[i][j] = vec[i * n + j];
    }

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      WI[i][j] = 0.;
    }


  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      W[i][j] = 0.;
    }

  for (i = 0; i < n; i++)
  {
    WI[i][i] = sqrt(fabs(M[i][i])) ;
  }


  for (i = 0; i < n; i++)
  {
    W[i][i] = 1 / sqrt(fabs(M[i][i])) ;
  }


  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      m0 = 0.0;
      for (kk = 0; kk < n; kk++)
      {
        MM[i][j] = W[i][kk] * M[kk][j] + m0;
        m0 = MM[i][j];
      }
    }
  }


  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      m0 = 0.0;
      for (kk = 0; kk < n; kk++)
      {
        M[i][j] = MM[i][kk] * W[kk][j] + m0;
        m0 = M[i][j];
      }
    }
  }


  for (i = 0; i < n; i++)
  {
    m0 = 0.0;
    for (j = 0; j < n; j++)
    {
      qq[i] = W[i][j] * qqq[j] + m0;
      m0 = qq[i];
    }

  }



  for (i = 0; i < n; i++)
  {
    q[i] = -qq[i];
  }




  qs = q[0];
  for (i = 1; i < n; i++)
    if (q[i] < qs) qs = q[i];

  if (qs >= 0)
  {
    for (i = 0; i < n; i++)
    {
      z[i] = 0.0;
      w[i] = q[i];
      z0 = 0.0;
    }
    *info = 0;
  }
  else
  {

    for (i = 0 ; i < n; i++)
      for (j = 0; j < 2 * n + 2; j++)
        if (i != j) A[i][j] = 0.;
        else A[i][j] = 1.0;

    for (i = 0 ; i < n; i++)
    {
      A[i][n] = -1.;
      for (j = 0; j < n; j++)
      {
        A[i][n + 1 + j] = -M[i][j];
      }
      A[i][2 * n + 1] = q[i];
    }


    pos_z0 = n;
    pos_q = 2 * n + 1;

    for (i = 0 ; i < n; i++)
      for (j = 0; j < 2 * n + 2; j++)
        if (A[i][j] > 1.e-16) E[i][j] = 16 * (1.e-16) * A[i][j];
        else  E[i][j] = -16 * (1.e-16) * A[i][j];

    basic = (int *)malloc(n * sizeof(int));
    for (i = 0 ; i < n; i++)
    {
      basic[i] = i;
      s[i] = i;
    }
    pos_pivot = pos_z0;
    iter = 0;

    trouve = 0;

    tempo = (double *)malloc(n * sizeof(double));



    while (iter < itt && !trouve)
    {
      iter = iter + 1;


      for (i = 0; i < n; i++)
      {
        tempo[i] = A[i][pos_q] / A[i][pos_pivot];
      }


      if (pos_pivot == pos_z0)
      {
        val = tempo[0];
        pos_out = 0;
        for (i = 0; i < n; i++)
        {
          if (tempo[i] > val)
          {
            val = tempo[i];
            pos_out = i;
          }
        }

        pos_out0 = pos_out;
      }
      else
      {
        dim_posA = 0;
        for (i = 0; i < n; i++)
        {
          if (A[i][pos_pivot] > 0)
            dim_posA = dim_posA + 1;
        }

        if (dim_posA == 0)
        {
          printf("there is a matter\n");
          /*  exit(2);*/

          break;
        }

        posA = (int *)malloc(dim_posA * sizeof(int));

        posAA = (int *)malloc(dim_posA * sizeof(int));

        compt = 0;
        for (i = 0; i < n; i++)
        {
          if (A[i][pos_pivot] > 0)
          {
            posA[compt] = i;
            compt = compt + 1;
          }
        }

        rtmp = (double *)malloc((dim_posA) * sizeof(double));
        ss = (int *)malloc((dim_posA) * sizeof(int));

        for (i = 0; i  < dim_posA; i++)
        {
          k = posA[i];
          rtmp[i] = tempo[k];
          ss[i] = s[k];
        }


        free(posA);
        val = -rtmp[0];
        pos = 0;

        for (i = 0; i < dim_posA; i++)
        {
          if (val < -rtmp[i])
          {
            val = -rtmp[i];
            pos = i;
          }
        }

        pos_out = pos;
        pos_out = ss[pos_out];
        free(ss);
        free(rtmp);
      }




      pivot = A[pos_out][pos_pivot];
      if (fabs(pivot) < 1.e-16)
      {
        printf("Error: nul pivot\n");
      }


      for (j = 0; j < 2 * n + 2; j++)
      {
        E[pos_out][j] = E[pos_out][j] / fabs(pivot);
        A[pos_out][j] = A[pos_out][j] / pivot;
      }

      A[pos_out][pos_pivot] = 1.;
      E[pos_out][pos_pivot] = 0.;

      for (i = 0 ; i < n; i++)
        for (j = 0 ; j < 2 * n + 2; j++)
        {
          AA[i][j] = A[i][j];
          EE[i][j] = E[i][j];
        }

      for (i = 0; i < n; i++)
      {
        if (i != pos_out)
        {
          for (j = 0; j < (2 * n + 2); j++)
          {
            aa = fabs(A[i][pos_pivot]) * E[pos_out][j];
            bb = E[i][pos_pivot] * fabs(A[pos_out][j]);
            cc = E[i][pos_pivot] * E[pos_out][j];
            EE[i][j] = E[i][j] + aa + bb + cc;
            EE[i][j] = 0.8 * EE[i][j];
          }

          for (j = 0; j < (2 * n + 2); j++)
          {
            AA[i][j] = A[i][j] - A[i][pos_pivot] * A[pos_out][j];
          }

          if (AA[i][pos_q] < (-EE[i][pos_q]))
          {
            printf("non positive AA %g\n", AA[i][pos_q]);
            /*exit(3);   */
            break;
          }

          AA[i][pos_pivot] = 0.;
          EE[i][pos_pivot] = 0.;

        }
      }

      for (i = 0 ; i < n; i++)
        for (j = 0 ; j < 2 * n + 2; j++)
        {
          A[i][j] = AA[i][j];
          E[i][j] = EE[i][j];
        }

      for (i = 0; i < n; i++)
        for (j = 0; j < 2 * n + 2; j++)
          E[i][j] = E[i][j] + (1.e-16) * fabs(A[i][j]);

      for (i = 0; i < n; i++)
      {
        A[i][pos_q] = fabs(A[i][pos_q]);
      }

      /*      printf("we enter %d in the base in the place of %d \n",pos_pivot, pos_out);*/

      nonbasic = basic[pos_out];
      basic[pos_out] = pos_pivot;

      if (nonbasic >= n)
      {
        pos_pivot = nonbasic - 1 - (n);
      }
      else
      {
        pos_pivot = nonbasic + (n) + 1;
      }

      trouve = (A[pos_out0][pos_q] < E[pos_out0][pos_q] || (nonbasic == pos_z0));



      /* Solution */


      sol = (double *)malloc((2 * n + 1) * sizeof(double));
      for (i = 0; i < 2 * n + 1 ; i++)
      {
        sol[i] = 0.;
      }


      for (i = 0; i < n; i++)
      {

        j = basic[i];
        sol[j] = A[i][pos_q];
      }

      for (i = 0; i < n; i++)
      {
        w[i] = sol[i];
        z[i] = sol[n + 1 + i];
      }
      z0 = sol[n];

      *it_end = iter;


      if (trouve == 0)
      {
        *info = 1;
      }
      else
      {
        *info = 0;
        break;
      }



    }


    /* verification */

    Mz = (double *)malloc(n * sizeof(double));
    incx = 1;
    incy = 1;
    alpha = -1.;
    beta = 1.;
    dcopy_(&n, w, &incx, Mz, &incy);
    dgemv_(&trans, &n, &n, &alpha, M, &n, z, &incx, &beta, Mz, &incy);
    beta = -1.;
    daxpy_(&n, &beta, q, &incx, Mz, &incy);
    num = ddot_(&n, Mz, &incx, Mz, &incy);
    den = ddot_(&n, q, &incx, q, &incy);
    resi = sqrt(num) / sqrt(den);
    *res = resi;
    minw = w[0];
    minz = z[0];
    for (i = 1; i < n; i++)
    {
      if (minw > w[i]) minw = w[i];
      if (minz > z[i]) minz = z[i];
    }


    pdt = ddot_(&n, w, &incx, z, &incy);
    printf("equilibrium %g min(w) %g min(z) %g w'z %g z0 %g\n", resi, minw, minz, pdt, z0);

    for (i = 0; i < n; i++)
    {
      m0 = 0.0;
      mm0 = 0.0;
      for (j = 0; j < n; j++)
      {
        zlem[i] = W[i][j] * z[j] + m0;
        m0 = zlem[i];
        wlem[i] = WI[i][j] * w[j] + mm0;
        mm0 = wlem[i];
      }
    }



    free(Mz);
    free(sol);
    free(posAA);
    free(tempo);
    free(basic);


  }

}
