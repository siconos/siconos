// This main file allows the primal resolution of relay problems:
//  try (z,w) such that:
//
//
//
//                    M z  + q = w                     (1)
//                    -w in diff(PSI) (z)              (2)
//                               [-b,a]
//
//
//  here M is an n by n  matrix, q an n-dimensional vector, z,w,b,a an n-dimensional vectors.
//
//  This system of equations and inequalities is solved thanks to rp_subroutine:
//        pr_gsnl  ( M, q, n, a, b, itermax, tol, chat, z, w, it_end, res, info)
//        pr_latin ( M, q, n, k_latin, a, b, itermax, tol, chat, z, w, it_end, res, info)
//
//
//  where _ itermax    is the maximum iterations required, it's an integer
//        _ res        is the residue, it's a float (positive float)
//        _ it_end     is the number of iterations carried out, it's an integer
//        _ tol        is the tolerance value desired, it's a positive float (if you make it non positive then you force the convergence)
//        _ k_latin    is the parameter research of the latin, it's a float (strictly positive)
//        _ a          is the upper bound, it's a vector of floats
//        _ b          is the down bound, it's a vector of floats
//        _ chat       is the output log identifiant
//        _ z and w    are the solutions of the problem
//        _ info       shows the termination reason,0 is successful (otherwise see the termination reason of the solver used), it's an integer.
//
//
//    For more information about the methods see the chapter 4 of the Siconos manual theory.
//
//
//
//  The subroutine's call is due to the function solve_rp:
//
//  int pr_solver ( double *M, double *q, int n, methode *pt, double *z, double *w )
//
//  where M        is an n by n matrix,
//        q        an n-dimensional vector,
//        n        is the row dimension of M,
//        pt       a pointer other a union ( methode see "solverpack.h").
//        z and w are n-dimensional vectors solution.
//
//  methode is a variable with a union type; in this union you find the structure (method_rp) that gives to the function
// rp_solver, the name and the parameters ( itermax, tol, k_latin, a, b, chat) of the method we want to use.
//  This function return an interger:  0 for successful
//
//
///////////////////////////////////////////////////////////////////////////////


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "solverpack.h"

#define CHAT


void test_mmc(void)
{
  FILE *f1, *f2, *f3, *f4;
  int i, j, nl, nc, nll, n = 15;
  int info1 = -1, info2 = -1;
  int nonsymmetric;

  double *q, *z1, *w1, *z2, *w2, *vec, *a, *b;
  double qi, Mij;


  char val[50], vall[50];
  method meth_pr1, meth_pr2;


  printf("\n* *** ******************** *** * \n");
  printf("* ***        TEST MMC      *** * \n");
  printf("* *** ******************** *** * \n");



  /*           Allocations                     */

  q   = (double*)malloc(n * sizeof(double));
  z1  = (double*)malloc(n * sizeof(double));
  w1  = (double*)malloc(n * sizeof(double));
  z2  = (double*)malloc(n * sizeof(double));
  w2  = (double*)malloc(n * sizeof(double));
  a   = (double*) malloc(n * sizeof(double));
  b   = (double*)malloc(n * sizeof(double));
  vec = (double*) malloc(n * n * sizeof(double));


  /*     Data loading of M , q , a and b           */

  if ((f1 = fopen("DATA/M_relay1.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }


  if ((f2 = fopen("DATA/q_relay1.dat", "r")) == NULL)
  {
    perror("fopen 2");
    exit(2);
  }


  if ((f3 = fopen("DATA/a_relay1.dat", "r")) == NULL)
  {
    perror("fopen 5");
    exit(5);
  }


  if ((f4 = fopen("DATA/b_relay1.dat", "r")) == NULL)
  {
    perror("fopen 6");
    exit(6);
  }



  while (!feof(f1))
  {
    fscanf(f1, "%d", &nl);
    fscanf(f1, "%d", &nc);
    fscanf(f1, "%s", val);
    Mij = atof(val);


    vec[(nc - 1)*n + nl - 1 ] = Mij;

  }




  while (!feof(f2))
  {
    fscanf(f2, "%d", &nll);
    fscanf(f2, "%s", vall);
    qi = atof(vall);
    q[nll - 1] = -qi;

  }


  while (!feof(f3))
  {
    fscanf(f3, "%d", &nll);
    fscanf(f3, "%s", vall);
    qi = atof(vall);
    a[nll - 1] = qi;
  }


  while (!feof(f4))
  {
    fscanf(f4, "%d", &nll);
    fscanf(f4, "%s", vall);
    qi = atof(vall);
    b[nll - 1] = qi;
  }


  fclose(f2);
  fclose(f4);
  fclose(f1);
  fclose(f3);



  /*           End od the data loading        */


  nonsymmetric = 0;


  /*           Is M symmetric ?                */

  for (i = 0 ; i < n ; ++i)
  {
    for (j = 0 ; j < i ; ++j)
    {
      if (abs(vec[i * n + j] - vec[j * n + i]) > 1e-16)
      {
        nonsymmetric = 1;
        break;
      }
    }
  }



#ifdef CHAT
  if (nonsymmetric) printf("\n !! WARNING !!\n M is a non symmetric matrix \n \n");
  else printf(" \n M is a symmetric matrix \n \n");
#endif





  /*           Initialization                 */

  strcpy(meth_pr1.pr.name, "NLGS");
  meth_pr1.pr.itermax  =  30000;
  meth_pr1.pr.tol      =  0.0000001;
  meth_pr1.pr.chat  =  0;


  strcpy(meth_pr2.pr.name, "Latin");
  meth_pr2.pr.itermax  =  1500;
  meth_pr2.pr.tol      =  0.0001;
  meth_pr2.pr.k_latin  =  0.003;
  meth_pr2.pr.chat  =  0;




  meth_pr1.pr.a = (double*)malloc(n * sizeof(double));
  meth_pr1.pr.b = (double*)malloc(n * sizeof(double));
  meth_pr2.pr.a = (double*)malloc(n * sizeof(double));
  meth_pr2.pr.b = (double*)malloc(n * sizeof(double));

  for (i = 0; i <= n - 1; i++)
  {

    meth_pr1.pr.a[i] =  a[i];
    meth_pr1.pr.b[i] = -b[i];
    meth_pr2.pr.a[i] =  a[i];
    meth_pr2.pr.b[i] = -b[i];


  }




#ifdef CHAT
  printf("**** NLGS TEST ****\n \n");
#endif

  info1 = pr_solver(vec, q, &n, &meth_pr1, z1, w1);

#ifdef CHAT
  printf("\n**** LATIN TEST ***\n \n");
#endif

  info2 = pr_solver(vec, q, &n, &meth_pr2, z2, w2);


#ifdef CHAT
  printf(" *** ************************************** ***\n");
  for (i = 0 ; i < n ; ++i)
    printf(" NLGS RESULT : %14.7g   LATIN RESULT : % 14.7g \n" , z1[i] , z2[i]);

  printf(" -- SOLVEUR ------ ITER ----- ERR\n");
  printf(" \n NLGS RESULT : %d   % 14.7g \n" , meth_pr1.pr.iter , meth_pr1.pr.err);
  printf(" LATIN RESULT : %d   % 14.7g \n" , meth_pr2.pr.iter , meth_pr2.pr.err);





  printf(" *** ************************************** ***\n");
#endif


  free(vec);
  free(q);
  free(z1);
  free(w1);
  free(z2);
  free(w2);
  free(a);
  free(b);
  free(meth_pr1.pr.a);
  free(meth_pr1.pr.b);
  free(meth_pr2.pr.a);
  free(meth_pr2.pr.b);
}


int main(void)
{
  test_mmc();
  return 1;
}


