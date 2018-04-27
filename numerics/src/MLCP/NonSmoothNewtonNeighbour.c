/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#include "NonSmoothNewton.h"
#include "NonSmoothNewtonNeighbour.h"
//#include "MixedLinearComplementarityProblem.h"
#include <string.h>
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "SiconosLapack.h"
#include "mlcp_enum_tool.h"
#include "numerics_verbose.h"



static  int sN ;
static  int sN2 ;

static  double * sphi_z ;
static  double * sdir_descent ;
static  double * sphi_zaux ;
static  double *sjacobianPhi_z ;
static  double *sjacobianPhi_zaux ;
static  double *sgrad_psi_z ;
static  double *sgrad_psi_zaux ;
static  double *sPrevDirDescent;
static  double *szaux ;
static  double *szzaux ;
static  double *sz2 ;
static  lapack_int* sipiv ;
static  int* sW2V;

static int scmp = 0;

static int sPlotMerit = 1;
static char fileName[64];
/* static char fileId[16]; */

static double* sZsol = 0;

static NewtonFunctionPtr* sFphi;
static NewtonFunctionPtr* sFjacobianPhi;


static void plotMerit(double *z, double psi_k, double descentCondition);
//static void plotMeritToZsol(double *z);
static int linesearch2_Armijo(int n, double *z, double psi_k, double descentCondition);
//static int lineSearch_Wolfe(double *z, double qp_0);
//static int NonMonotomnelineSearch(double *z, double Rk);


/************************************************************************/
/*useful for debug*/
void NSNN_thisIsTheSolution(int n, double * z)
{
  sZsol = (double *)malloc(n * sizeof(double));
  memcpy(sZsol, z, n * sizeof(double));
}
void  NSNN_reset()
{
  if (sZsol)
    free(sZsol);
  sZsol = 0;
}
void plotMerit(double *z, double psi_k, double descentCondition)
{
  int incx = 1, incy = 1;
  double q_0, q_tk, qp_tk, merit_k;
  /* double tmin = 1e-12; */
  double tk = 1, aux;
  double m1 = 1e-4;
  double Nstep = 0;
  int i = 0;

  FILE *fp;

  (*sFphi)(sN, z, sphi_z, 0);
  aux = cblas_dnrm2(sN, sphi_z, 1);
  /* Computes merit function */
  aux = 0.5 * aux * aux;
  printf("plot psi_z %e\n", aux);


  if (!sPlotMerit)
    return;

  if (sPlotMerit)
  {
    /*    sPlotMerit=0;*/
    strcpy(fileName, "outputLS");


    (*sFphi)(sN, z, sphi_z, 0);
    q_0 =  cblas_dnrm2(sN, sphi_z , incx);
    q_0 = 0.5 * q_0 * q_0;

    fp = fopen(fileName, "w");

    /*    sPlotMerit=0;*/
    tk = 5e-7;
    aux = -tk;
    Nstep = 1e4;
    for (i = 0; i < 2 * Nstep; i++)
    {
      cblas_dcopy(sN, z, incx, sz2, incx);
      cblas_daxpy(sN , aux , sdir_descent , incx , sz2 , incy);
      (*sFphi)(sN, sz2, sphi_z, 0);
      q_tk =  cblas_dnrm2(sN, sphi_z , incx);
      q_tk = 0.5 * q_tk * q_tk;


      (*sFjacobianPhi)(sN, sz2, sjacobianPhi_z, 1);
      /* Computes the jacobian of the merit function, jacobian_psi = transpose(jacobianPhiMatrix).phiVector */
      cblas_dgemv(CblasColMajor,CblasTrans, sN, sN, 1.0, sjacobianPhi_z, sN, sphi_z, incx, 0.0, sgrad_psi_z, incx);
      qp_tk = cblas_ddot(sN, sgrad_psi_z, 1, sdir_descent, 1);

      merit_k = psi_k + m1 * aux * descentCondition;


      fprintf(fp, "%e %.16e %.16e %e\n", aux, q_tk, merit_k, qp_tk);
      if (i == Nstep - 1)
        aux = 0;
      else
        aux += tk / Nstep;
    }

    fclose(fp);
  }
}
/* void plotMeritToZsol(double *z) */
/* { */
/*   int incx = 1, incy = 1; */
/*   double q_0, q_tk; */
/*   /\*   double merit_k;  *\/ */
/*   /\*  double tmin = 1e-12; *\/ */
/*   double tk = 1; */
/*   /\*   double m1=0.5; *\/ */
/*   double aux; */
/*   int i = 0; */
/*   int ii; */
/*   if (!sPlotMerit || !sZsol) */
/*     return; */
/*   FILE *fp; */

/*   for (ii = 0; ii < sN; ii++) */
/*     szzaux[ii] = sZsol[ii] - z[ii]; */



/*   if (sPlotMerit) */
/*   { */
/*     /\*    sPlotMerit=0;*\/ */
/*     strcpy(fileName, "outputLSZsol"); */
/*     (*sFphi)(sN, z, sphi_z, 0); */
/*     q_0 =  cblas_dnrm2(sN, sphi_z , incx); */
/*     q_0 = 0.5 * q_0 * q_0; */

/*     fp = fopen(fileName, "w"); */

/*     /\*    sPlotMerit=0;*\/ */
/*     tk = 1; */
/*     aux = -tk; */
/*     for (i = 0; i < 2e3; i++) */
/*     { */
/*       cblas_dcopy(sN, z, incx, sz2, incx); */
/*       cblas_daxpy(sN , aux , szzaux , incx , sz2 , incy); */
/*       (*sFphi)(sN, sz2, sphi_z, 0); */
/*       q_tk =  cblas_dnrm2(sN, sphi_z , incx); */
/*       q_tk = 0.5 * q_tk * q_tk; */
/*       fprintf(fp, "%e %e\n", aux, q_tk); */
/*       aux += tk / 1e3; */
/*     } */

/*     fclose(fp); */
/*   } */
/* }*/

/************************************************************************/

/* Linesearch */
int linesearch2_Armijo(int n, double *z, double psi_k, double descentCondition)
{

  /* IN :
     psi_k (merit function for current iteration)
     jacobian_psi_k (jacobian of the merit function)
     dk: descent direction

     OUT: tk, z
  */

  double m1 = 0.1;
  double tk = 1;
  double tkl, tkr, tkaux;
  int incx = 1, incy = 1;
  double merit, merit_k;
  double tmin = 1e-14;
  double qp_tk;

  /*  cblas_dcopy(sN, z, incx,sz2,incx);*/

  /* z1 = z0 + dir */
  /*  cblas_daxpy(n , 1.0 , sdir_descent , incx , z , incy );*/

  tk = 3.25;


  while (tk > tmin)
  {

    /* Computes merit function = 1/2*norm(phi(z_{k+1}))^2 */
    cblas_dcopy(sN, z, incx, sz2, incx);
    cblas_daxpy(n , tk , sdir_descent , incx , sz2 , incy);


    (*sFphi)(n, sz2, sphi_z, 0);
    merit =  cblas_dnrm2(n, sphi_z , incx);
    merit = 0.5 * merit * merit;
    merit_k = psi_k + m1 * tk * descentCondition;
    if (merit < merit_k)
    {
      tkl = 0;
      tkr = tk;

      /*calcul merit'(tk)*/
      (*sFjacobianPhi)(sN, sz2, sjacobianPhi_z, 1);
      /* Computes the jacobian of the merit function, jacobian_psi = transpose(jacobianPhiMatrix).phiVector */
      cblas_dgemv(CblasColMajor,CblasTrans, sN, sN, 1.0, sjacobianPhi_z, sN, sphi_z, incx, 0.0, sgrad_psi_zaux, incx);
      qp_tk = cblas_ddot(sN, sgrad_psi_zaux, 1, sdir_descent, 1);

      if (qp_tk > 0)
      {
        while (fabs(tkl - tkr) > tmin)
        {
          tkaux = 0.5 * (tkl + tkr);
          cblas_dcopy(sN, z, incx, sz2, incx);
          cblas_daxpy(n , tkaux , sdir_descent , incx , sz2 , incy);
          /*calcul merit'(tk)*/
          (*sFphi)(n, sz2, sphi_z, 0);
          (*sFjacobianPhi)(sN, sz2, sjacobianPhi_z, 1);
          /* Computes the jacobian of the merit function, jacobian_psi = transpose(jacobianPhiMatrix).phiVector */
          cblas_dgemv(CblasColMajor,CblasTrans, sN, sN, 1.0, sjacobianPhi_z, sN, sphi_z, incx, 0.0, sgrad_psi_zaux, incx);
          qp_tk = cblas_ddot(sN, sgrad_psi_zaux, 1, sdir_descent, 1);
          if (qp_tk > 0)
          {
            tkr = tkaux;
          }
          else
          {
            tkl = tkaux;
          }
        }
      }

      /* printf("merit = %e, merit_k=%e,tk= %e,tkaux=%e \n",merit,merit_k,tk,tkaux);*/
      cblas_dcopy(sN, sz2, incx, z, incx);
      break;
    }
    tk = tk * 0.5;
  }
  if (tk <= tmin)
  {
    cblas_dcopy(sN, sz2, incx, z, incx);
    printf("NonSmoothNewton::linesearch2_Armijo warning, resulting tk=%e < tmin, linesearch stopped.\n", tk);
    return 0;

  }
  return 1;

}


/* int lineSearch_Wolfe(double *z, double qp_0) */
/* { */
/*   int incx = 1, incy = 1; */
/*   double q_0, q_tk, qp_tk; */
/*   double tmin = 1e-12; */
/*   int maxiter = 100; */
/*   int niter = 0; */
/*   double tk = 1; */
/*   double tg, td; */
/*   double m1 = 0.1; */
/*   double m2 = 0.9; */


/*   (*sFphi)(sN, z, sphi_z, 0); */
/*   q_0 =  cblas_dnrm2(sN, sphi_z , incx); */
/*   q_0 = 0.5 * q_0 * q_0; */

/*   tg = 0; */
/*   td = 10e5; */

/*   tk = (tg + td) / 2.0; */

/*   while (niter < maxiter || (td - tg) < tmin) */
/*   { */
/*     niter++; */
/*     /\*q_tk = 0.5*|| phi(z+tk*d) ||*\/ */
/*     cblas_dcopy(sN, z, incx, sz2, incx); */
/*     cblas_daxpy(sN , tk , sdir_descent , incx , sz2 , incy); */
/*     (*sFphi)(sN, sz2, sphi_z, 0); */
/*     q_tk =  cblas_dnrm2(sN, sphi_z , incx); */
/*     q_tk = 0.5 * q_tk * q_tk; */

/*     (*sFjacobianPhi)(sN, sz2, sjacobianPhi_z, 1); */
/*     /\* Computes the jacobian of the merit function, jacobian_psi = transpose(jacobianPhiMatrix).phiVector *\/ */
/*     cblas_dgemv(CblasColMajor,CblasTrans, sN, sN, 1.0, sjacobianPhi_z, sN, sphi_z, incx, 0.0, sgrad_psi_z, incx); */
/*     qp_tk = cblas_ddot(sN, sgrad_psi_z, 1, sdir_descent, 1); */
/*     if (qp_tk <  m2 * qp_0 && q_tk < q_0 + m1 * tk * qp_0) */
/*     { */
/*       /\*too small*\/ */
/*       if (niter == 1) */
/*         break; */
/*       tg = tk; */
/*       tk = (tg + td) / 2.0; */
/*       continue; */
/*     } */
/*     else if (q_tk > q_0 + m1 * tk * qp_0) */
/*     { */
/*       /\*too big*\/ */
/*       td = tk; */
/*       tk = (tg + td) / 2.0; */
/*       continue; */
/*     } */
/*     else */
/*       break; */
/*   } */

/*   cblas_dcopy(sN, sz2, incx, z, incx); */

/*   if ((td - tg) <= tmin) */
/*   { */
/*     printf("NonSmoothNewton2::lineSearchWolfe warning, resulting tk < tmin, linesearch stopped.\n"); */
/*     return 0; */
/*   } */
/*   return 1; */

/* } */

/* nt NonMonotomnelineSearch(double *z, double Rk) */
/* { */
/*   int incx = 1, incy = 1; */
/*   double q_0, q_tk; */
/*   /\*   double merit_k;  *\/ */
/*   double tmin = 1e-12; */
/*   double tmax = 1000; */
/*   double tk = 1; */
/*   /\*   double m1=0.5; *\/ */



/*   (*sFphi)(sN, z, sphi_z, 0); */
/*   q_0 =  cblas_dnrm2(sN, sphi_z , incx); */
/*   q_0 = 0.5 * q_0 * q_0; */


/*   while ((tmax - tmin) > 1e-1) */
/*   { */
/*     tk = (tmax + tmin) / 2; */
/*     /\*q_tk = 0.5*|| phi(z+tk*d) ||*\/ */
/*     cblas_dcopy(sN, z, incx, sz2, incx); */
/*     cblas_daxpy(sN , tk , sdir_descent , incx , sz2 , incy); */
/*     (*sFphi)(sN, sz2, sphi_z, 0); */
/*     q_tk =  cblas_dnrm2(sN, sphi_z , incx); */
/*     q_tk = 0.5 * q_tk * q_tk; */
/*     if (fabs(q_tk - q_0) < Rk) */
/*       tmin = tk; */
/*     else */
/*       tmax = tk; */
/*   } */
/*   printf("NonMonotomnelineSearch, tk = %e\n", tk); */
/*   cblas_dcopy(sN, sz2, incx, z, incx); */

/*   if (tk <= tmin) */
/*   { */
/*     printf("NonMonotomnelineSearch warning, resulting tk < tmin, linesearch stopped.\n"); */
/*     return 0; */
/*   } */
/*   return 1; */

/* } */


int nonSmoothNewtonNeigh_getNbIWork(int n, int m)
{
  return 2 * (n + m);

}
int nonSmoothNewtonNeigh_getNbDWork(int n, int m)
{
  return (11 + 2 * (n + m)) * (n + m) + 1;
}

double * nonSmoothNewtonNeighInitMemory(int n, double * dWork, int * iWork)
{

  if (dWork == NULL || iWork == NULL)
  {
    fprintf(stderr, "nonSmoothNewtonNeighInitMemory, memory allocation failed.\n");
    exit(EXIT_FAILURE);
  }
  sN = n;
  sN2 = n * n;
  sphi_z = dWork;//(double*)malloc(n*sizeof(*phi_z));
  sdir_descent = sphi_z + sN;//(double*)malloc(n*sizeof(double));
  sphi_zaux = sdir_descent + sN ; //(double*)malloc(n*sizeof(double));
  sjacobianPhi_z = sphi_zaux + sN; //(double*)malloc(n2*sizeof(*jacobianPhi_z));
  sjacobianPhi_zaux = sjacobianPhi_z + sN2;//(double*)malloc(n2*sizeof(double));
  sgrad_psi_z = sjacobianPhi_zaux + sN2;//(double*)malloc(n*sizeof(*jacobian_psi_z));
  sgrad_psi_zaux = sgrad_psi_z + sN;//(double*)malloc(n*sizeof(double));
  sPrevDirDescent = sgrad_psi_zaux + sN;//(double*)malloc((n)*sizeof(double));
  szaux = sPrevDirDescent + sN;//(double*)malloc(n*sizeof(double));
  szzaux = szaux + sN; //(double*)malloc(n*sizeof(double));
  sz2 = szzaux + sN;// size n


  sipiv = iWork;//(int *)malloc(n*sizeof(*ipiv));
  sW2V = sipiv + sN;

  return sz2 + sN;

}


int nonSmoothNewtonNeigh(int n, double* z, NewtonFunctionPtr* phi, NewtonFunctionPtr* jacobianPhi, int* iparam, double* dparam)
{


  int itermax = iparam[0]; // maximum number of iterations allowed
  int iterMaxWithSameZ = itermax / 4;
  int niter = 0; // current iteration number
  double tolerance = dparam[0];
  /*   double coef; */
  sFphi = phi;
  sFjacobianPhi = jacobianPhi;
  //  verbose=1;
  if (verbose > 0)
  {
    printf(" ============= Starting of Newton process =============\n");
    printf(" - tolerance: %14.7e\n - maximum number of iterations: %i\n", tolerance, itermax);
  }

  int incx = 1;
  /*   int n2 = n*n; */
  lapack_int infoDGESV;

  /** merit function and its jacobian */
  double psi_z;

  /** The algorithm is alg 4.1 of the paper of Kanzow and Kleinmichel, "A new class of semismooth Newton-type methods
      for nonlinear complementarity problems", in Computational Optimization and Applications, 11, 227-251 (1998).

      We try to keep the same notations
  */

  double rho = 1e-8;
  double descentCondition, criterion, norm_jacobian_psi_z, normPhi_z;
  double p = 2.1;
  double terminationCriterion = 1;
  double norm;
  int findNewZ, i, j, NbLookingForANewZ;
  /*   int naux=0; */
  double aux = 0;
  /*   double aux1=0; */
  int ii;
  int resls = 1;
  /*   char c; */
  /*  double * oldz; */
  /*  oldz=(double*)malloc(n*sizeof(double));*/

  NbLookingForANewZ = 0;

  /** Iterations ... */
  while ((niter < itermax) && (terminationCriterion > tolerance))
  {
    scmp++;
    ++niter;
    /** Computes phi and its jacobian */
    if (sZsol)
    {
      for (ii = 0; ii < sN; ii++)
        szzaux[ii] = sZsol[ii] - z[ii];
      printf("dist zzsol %.32e.\n", cblas_dnrm2(n, szzaux, 1));
    }

    (*sFphi)(n, z, sphi_z, 0);
    (*sFjacobianPhi)(n, z, sjacobianPhi_z, 1);
    /* Computes the jacobian of the merit function, jacobian_psi = transpose(jacobianPhiMatrix).phiVector */
    cblas_dgemv(CblasColMajor,CblasTrans, n, n, 1.0, sjacobianPhi_z, n, sphi_z, incx, 0.0, sgrad_psi_z, incx);
    norm_jacobian_psi_z = cblas_dnrm2(n, sgrad_psi_z, 1);

    /* Computes norm2(phi) */
    normPhi_z = cblas_dnrm2(n, sphi_z, 1);
    /* Computes merit function */
    psi_z = 0.5 * normPhi_z * normPhi_z;

    if (normPhi_z < tolerance)
    {
      /*it is the solution*/
      terminationCriterion = tolerance / 2.0;
      break;
    }

    if (verbose > 0)
    {
      printf("Non Smooth Newton, iteration number %i, norm grad psi= %14.7e , psi = %14.7e, normPhi = %e .\n", niter, norm_jacobian_psi_z, psi_z, normPhi_z);
      printf(" -----------\n");
    }

    NbLookingForANewZ++;

    if (niter > 2)
    {
      if (10 * norm_jacobian_psi_z < tolerance || !resls || NbLookingForANewZ > iterMaxWithSameZ)
      {
        NbLookingForANewZ = 0;
        resls = 1;
        /*   if (NbLookingForANewZ % 10 ==1 && 0){
          printf("Try NonMonotomnelineSearch\n");
          cblas_dcopy(n,sgrad_psi_z,1,sdir_descent,1);
          cblas_dscal( n , -1.0 ,sdir_descent,incx);
          NonMonotomnelineSearch( z,  phi, 10);
          continue;
        }
        */

        /* FOR DEBUG ONLY*/
        if (sZsol)
        {
          printf("begin plot prev dir\n");
          plotMerit(z, 0, 0);
          printf("end\n");
          /*     gets(&c);*/
          (*sFphi)(n, sZsol, szaux, 0);
          printf("value psi(zsol)=%e\n", cblas_dnrm2(n, szaux, 1));
          cblas_dcopy(n, sZsol, incx, szaux, incx);
          cblas_daxpy(n , -1 , z , 1 , szaux , 1);
          printf("dist to sol %e \n", cblas_dnrm2(n, szaux, 1));
          for (ii = 0; ii < n; ii++)
            sdir_descent[ii] = sZsol[ii] - z[ii];

          aux = norm;
          norm = 1;
          printf("begin plot zzsol dir\n");
          plotMerit(z, 0, 0);
          printf("end\n");
          /*     gets(&c);*/
          norm = aux;
        }

        printf("looking for a new Z...\n");
        /*may be a local minimal*/
        /*find a gradiant going out of this cul-de-sac.*/
        norm = n / 2;
        findNewZ = 0;
        for (j = 0; j < 20; j++)
        {

          for (i = 0; i < n; i++)
          {
            if (sZsol)
            {
              /* FOR DEBUG ONLY*/
              (*sFphi)(n, sZsol, sphi_zaux, 0);
              norm = cblas_dnrm2(n, sphi_zaux, 1);
              printf("Norm of the sol %e.\n", norm);

              for (ii = 0; ii < n; ii++)
                sdir_descent[ii] = sZsol[ii] - z[ii];
              norm = 1;
            }
            else
            {
              for (ii = 0; ii < n; ii++)
              {
                sdir_descent[ii] = 1.0 * rand();
              }
              cblas_dscal(n, 1 / cblas_dnrm2(n, sdir_descent, 1), sdir_descent, incx);
              cblas_dscal(n, norm, sdir_descent, incx);
            }
            cblas_dcopy(n, z, incx, szaux, incx);
            // cblas_dscal(n,0.0,zaux,incx);
            /* zaux = z + dir */
            cblas_daxpy(n , norm , sdir_descent , 1 , szaux , 1);
            /* Computes the jacobian of the merit function, jacobian_psi_zaux = transpose(jacobianPhi_zaux).phi_zaux */
            (*sFphi)(n, szaux, sphi_zaux, 0);
            (*sFjacobianPhi)(n, szaux, sjacobianPhi_zaux, 1);

            /* FOR DEBUG ONLY*/
            if (sZsol)
            {
              aux = cblas_dnrm2(n, sphi_zaux, 1);
              printf("Norm of the sol is now %e.\n", aux);
              for (ii = 0; ii < n; ii++)
                printf("zsol %e zaux %e \n", sZsol[ii], szaux[ii]);
            }


            cblas_dgemv(CblasColMajor, CblasTrans, n, n, 1.0, sjacobianPhi_zaux, n, sphi_zaux, incx, 0.0, sgrad_psi_zaux, incx);
            cblas_dcopy(n, szaux, 1, szzaux, 1);
            cblas_daxpy(n , -1 , z , incx , szzaux , incx);
            /*zzaux must be a descente direction.*/
            /*ie jacobian_psi_zaux.zzaux <0
            printf("jacobian_psi_zaux : \n");*/
            /*cblas_dcopy(n,sdir,incx,sdir_descent,incx);
            plotMerit(z, phi);*/


            aux = cblas_ddot(n, sgrad_psi_zaux, 1, szzaux, 1);
            /*       aux1 = cblas_dnrm2(n,szzaux,1);
            aux1 = cblas_dnrm2(n,sgrad_psi_zaux,1);*/
            aux = aux / (cblas_dnrm2(n, szzaux, 1) * cblas_dnrm2(n, sgrad_psi_zaux, 1));
            /*       printf("aux: %e\n",aux);*/
            if (aux < 0.1 * (j + 1))
            {
              //zaux is the new point.
              findNewZ = 1;
              cblas_dcopy(n, szaux, incx, z, incx);
              break;
            }
          }
          if (findNewZ)
            break;
          if (j == 10)
          {
            norm = n / 2;
          }
          else if (j > 10)
            norm = -2 * norm;
          else
            norm = -norm / 2.0;
        }
        if (! findNewZ)
        {
          printf("failed to find a new z\n");
          /* exit(1);*/
          continue;

        }
        else
          continue;
      }
    }

    /* Stops if the termination criterion is satisfied */
    terminationCriterion = norm_jacobian_psi_z;
    /*      if(terminationCriterion < tolerance){
    break;
    }*/

    /* Search direction calculation
    Find a solution dk of jacobianPhiMatrix.d = -phiVector.
    dk is saved in phiVector.
    */
    cblas_dscal(n , -1.0 , sphi_z, incx);
    DGESV(n, 1, sjacobianPhi_z, n, sipiv, sphi_z, n, &infoDGESV);
    if (infoDGESV)
    {
      printf("DGEV error %d.\n", infoDGESV);
    }
    cblas_dcopy(n, sphi_z, 1, sdir_descent, 1);
    criterion = cblas_dnrm2(n, sdir_descent, 1);
    /*      printf("norm dir descent %e\n",criterion);*/

    /*printf("begin plot descent dir\n");
    plotMerit(z, phi);
    printf("end\n");
          gets(&c);*/

    /*printf("begin plot zzsol dir\n");
    plotMeritToZsol(z,phi);
    printf("end\n");
          gets(&c);*/


    /*
    norm = cblas_dnrm2(n,sdir_descent,1);
    printf("norm desc %e \n",norm);
    cblas_dscal( n , 1/norm , sdir_descent, 1);
    */
    /* descentCondition = jacobian_psi.dk */
    descentCondition = cblas_ddot(n, sgrad_psi_z,  1,  sdir_descent, 1);

    /* Criterion to be satisfied: error < -rho*norm(dk)^p */
    criterion = -rho * pow(criterion, p);
    /*      printf("ddddddd %d\n",scmp);
    if (scmp>100){
    NM_dense_display(sjacobianPhi_z,n,n,n);
    exit(1);
    }*/

//    if ((infoDGESV != 0 || descentCondition > criterion) && 0)
//    {
//      printf("no a desc dir, get grad psy\n");
      /* dk = - jacobian_psi (remind that dk is saved in phi_z) */
//      cblas_dcopy(n, sgrad_psi_z, 1, sdir_descent, 1);
//      cblas_dscal(n , -1.0 , sdir_descent, incx);
      /*DEBUG ONLY*/
      /*printf("begin plot new descent dir\n");
      plotMerit(z);
      printf("end\n");
       gets(&c);*/
//    }
    /*      coef=fabs(norm_jacobian_psi_z*norm_jacobian_psi_z/descentCondition);
    if (coef <1){
    cblas_dscal(n,coef,sdir_descent,incx);
    printf("coef %e norm dir descent is now %e\n",coef,cblas_dnrm2(n,sdir_descent,1));
    }*/


    /* Step-3 Line search: computes z_k+1 */
    /*linesearch_Armijo(n,z,sdir_descent,psi_z, descentCondition, phi);*/
    /*            if (niter == 10){
    printf("begin plot new descent dir\n");
    plotMerit(z);
    printf("end\n");
     gets(&c);
    }*/
    /*      memcpy(oldz,z,n*sizeof(double));*/

    resls = linesearch2_Armijo(n, z, psi_z, descentCondition);
    if (!resls && niter > 1)
    {

      /* NM_dense_display(sjacobianPhi_z,n,n,n);
      printf("begin plot new descent dir\n");
      plotMerit(oldz,psi_z, descentCondition);
      printf("end\n");
      gets(&c);*/
    }


    /*      lineSearch_Wolfe(z, descentCondition, phi,jacobianPhi);*/
    /*      if (niter>3){
    printf("angle between prev dir %e.\n",acos(cblas_ddot(n, sdir_descent,  1,  sPrevDirDescent, 1)/(cblas_dnrm2(n,sdir_descent,1)*cblas_dnrm2(n,sPrevDirDescent,1))));
    }*/
    cblas_dcopy(n, sdir_descent, 1, sPrevDirDescent, 1);

    /*      for (j=20;j<32;j++){
    if (z[j]<0)
    z[j]=0;
    }*/

    /*      if( 1 || verbose>0)
    {
     printf("Non Smooth Newton, iteration number %i, error grad equal to %14.7e , psi value is %14.7e .\n",niter, terminationCriterion,psi_z);
       printf(" -----------\n");
       }*/
  }

  /* Total number of iterations */
  iparam[1] = niter;
  /* Final error */
  dparam[1] = terminationCriterion;

  /** Free memory*/

  if (verbose > 0)
  {
    if (dparam[1] > tolerance)
      printf("Non Smooth Newton warning: no convergence after %i iterations\n" , niter);

    else
      printf("Non Smooth Newton: convergence after %i iterations\n" , niter);
    printf(" The residue is : %e \n", dparam[1]);
  }

  /*  free(oldz);*/

  if (dparam[1] > tolerance)
    return 1;
  else return 0;
}
