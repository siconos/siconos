/* Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
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
#include "NonSmoothNewton.h"
#include "NonSmoothNewtonNeighbour.h"
#include "Numerics_Options.h"
#include "LA.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"

#include "mlcp_enum_tool.h"



static  int sN ;
static  int sN2 ;

static  double * sphi_z ;
static  double * sdir_descent ;
static  double * sphi_zaux ;
static  double *sjacobianPhi_z ;
static  double *sjacobianPhi_zaux ;
static  double *sjacobian_psi_z ;
static  double *sjacobian_psi_zaux ;
static  double *sdir;
static  double *szaux ;
static  double *szzaux ;
static  double *sz2 ;
static  int* sipiv ;
static  int* sW2V;

static int sPlotMerit = 1;
static char fileName[64];
static char fileId[16];

static double* sZsol = 0;
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
void plotMerit(double *z,  NewtonFunctionPtr* phi)
{
  int incx = 1, incy = 1;
  double q_0, q_tk, merit_k;
  double tmin = 1e-12;
  double tk = 1;
  double m1 = 0.5;
  int i = 0;

  FILE *fp;

  if (sPlotMerit)
  {
    /*    sPlotMerit=0;*/
    strcpy(fileName, "outputLS");


    (*phi)(sN, z, sphi_z, 0);
    q_0 =  DNRM2(sN, sphi_z , incx);
    q_0 = 0.5 * q_0 * q_0;

    fp = fopen(fileName, "w");

    /*    sPlotMerit=0;*/
    tk = 0;
    for (i = 0; i < 2e4; i++)
    {
      DCOPY(sN, z, incx, sz2, incx);
      DAXPY(sN , tk , sdir_descent , incx , sz2 , incy);
      (*phi)(sN, sz2, sphi_z, 0);
      q_tk =  DNRM2(sN, sphi_z , incx);
      q_tk = 0.5 * q_tk * q_tk;
      fprintf(fp, "%e %e\n", tk, q_tk);
      tk += 0.5e-4;
    }

    fclose(fp);
  }
}
/************************************************************************/

int lineSearch_Wolfe(double *z, double qp_0,  NewtonFunctionPtr* phi)
{
  int incx = 1, incy = 1;
  double q_0, q_tk, merit_k;
  double tmin = 1e-12;
  double tk = 1;
  double m1 = 0.5;
  int i = 0;


  (*phi)(sN, z, sphi_z, 0);
  q_0 =  DNRM2(sN, sphi_z , incx);
  q_0 = 0.5 * q_0 * q_0;


  tk = 1;
  while (tk > tmin)
  {
    /*q_tk = 0.5*|| phi(z+tk*d) ||*/
    DCOPY(sN, z, incx, sz2, incx);
    DAXPY(sN , tk , sdir_descent , incx , sz2 , incy);
    (*phi)(sN, sz2, sphi_z, 0);
    q_tk =  DNRM2(sN, sphi_z , incx);
    q_tk = 0.5 * q_tk * q_tk;
    if (q_tk < q_0 + m1 * tk * qp_0)
      break;
    tk = tk * 0.5;

  }

  DCOPY(sN, sz2, incx, z, incx);

  if (tk <= tmin)
  {
    printf("NonSmoothNewton2::lineSearchWolfe warning, resulting tk < tmin, linesearch stopped.\n");
    return 0;
  }
  return 1;

}
double * nonSmoothNewtonNeighInitMemory(int n, double * dWork, int * iWork)
{

  sN = n;
  sN2 = n * n;
  sphi_z = dWork;//(double*)malloc(n*sizeof(*phi_z));
  sdir_descent = sphi_z + sN;//(double*)malloc(n*sizeof(double));
  sphi_zaux = sdir_descent + sN ; //(double*)malloc(n*sizeof(double));
  sjacobianPhi_z = sphi_zaux + sN; //(double*)malloc(n2*sizeof(*jacobianPhi_z));
  sjacobianPhi_zaux = sjacobianPhi_z + sN2;//(double*)malloc(n2*sizeof(double));
  sjacobian_psi_z = sjacobianPhi_zaux + sN2;//(double*)malloc(n*sizeof(*jacobian_psi_z));
  sjacobian_psi_zaux = sjacobian_psi_z + sN;//(double*)malloc(n*sizeof(double));
  sdir = sjacobian_psi_zaux + sN;//(double*)malloc((n)*sizeof(double));
  szaux = sdir + sN;//(double*)malloc(n*sizeof(double));
  szzaux = szaux + sN; //(double*)malloc(n*sizeof(double));
  sz2 = szzaux + sN;// size n
  sipiv = iWork;//(int *)malloc(n*sizeof(*ipiv));
  sW2V = sipiv + sN;
  if (dWork == NULL || iWork == NULL)
  {
    fprintf(stderr, "NonSmoothNewton, memory allocation failed.\n");
    exit(EXIT_FAILURE);
  }

  return sz2 + sN;

}


int nonSmoothNewtonNeigh(int n, double* z, NewtonFunctionPtr* phi, NewtonFunctionPtr* jacobianPhi, int* iparam, double* dparam)
{


  int itermax = iparam[0]; // maximum number of iterations allowed
  int niter = 0; // current iteration number
  double tolerance = dparam[0];
  if (verbose > 0)
  {
    printf(" ============= Starting of Newton process =============\n");
    printf(" - tolerance: %14.7e\n - maximum number of iterations: %i\n", tolerance, itermax);
  }

  int incx = 1;
  int n2 = n * n;
  int infoDGESV;

  /** merit function and its jacobian */
  double psi_z, prev_psi_z;

  /** The algorithm is alg 4.1 of the paper of Kanzow and Kleinmichel, "A new class of semismooth Newton-type methods
      for nonlinear complementarity problems", in Computational Optimization and Applications, 11, 227-251 (1998).

      We try to keep the same notations
  */

  double rho = 1e-8;
  double descentCondition, criterion, norm_jacobian_psi_z, normPhi_z;
  double p = 2.1;
  double terminationCriterion = 1;
  double prev_norm_jacobian_psi_z = 0;
  double norm;
  int findNewZ, i, j;
  int lastN = 0;
  int naux = 0;
  double aux = 0, aux1 = 0;
  int ii;
  int useNewZ = 1;
  char c;

  prev_psi_z = 0;

  /** Iterations ... */
  while ((niter < itermax) && (terminationCriterion > tolerance))
  {
    ++niter;
    /** Computes phi and its jacobian */
    (*phi)(n, z, sphi_z, 0);
    (*jacobianPhi)(n, z, sjacobianPhi_z, 1);
    /* Computes the jacobian of the merit function, jacobian_psi = transpose(jacobianPhiMatrix).phiVector */
    DGEMV(LA_TRANS, n, n, 1.0, sjacobianPhi_z, n, sphi_z, incx, 0.0, sjacobian_psi_z, incx);
    norm_jacobian_psi_z = DNRM2(n, sjacobian_psi_z, 1);

    /* Computes norm2(phi) */
    normPhi_z = DNRM2(n, sphi_z, 1);
    /* Computes merit function */
    psi_z = 0.5 * normPhi_z * normPhi_z;

    if (verbose > 0)
    {
      printf("Non Smooth Newton, iteration number %i, error grad equal to %14.7e , psi value is %14.7e .\n", niter, norm_jacobian_psi_z, psi_z);
      printf(" -----------------------------------------------------------------------\n");
    }



    if (niter > 2 && useNewZ)
    {
      if (10 * (itermax)*fabs(prev_psi_z - psi_z) < psi_z)
      {
        /*    plotMerit(z, descentCondition, phi);*/

        if (sZsol && 0)
        {
          DCOPY(n, sZsol, incx, szaux, incx);
          DAXPY(n , -1 , z , 1 , szaux , 1);
          printf("dist to sol %e \n", DNRM2(n, szaux, 1));
          for (ii = 0; ii < n; ii++)
            sdir_descent[ii] = sZsol[ii] - z[ii];

          aux = norm;
          norm = 1;
          plotMerit(z, phi);
          printf("type1\n");
          gets(&c);
          printf("type2\n");
          norm = aux;

        }

        prev_norm_jacobian_psi_z = 0;
        prev_psi_z = 0;
        printf("looking for a new Z...\n");
        /*may be a local minimal*/
        /*find a gradiant going out of this cul-de-sac.*/
        norm = n / 2;
        findNewZ = 0;
        for (j = 0; j < 10; j++)
        {

          for (i = 0; i < n; i++)
          {
            if (sZsol && 0)
            {
              (*phi)(n, sZsol, sphi_z, 0);
              norm = DNRM2(n, sphi_z, 1);
              printf("Norm of the sol %e.\n", norm);

              for (ii = 0; ii < n; ii++)
                sdir[ii] = sZsol[ii] - z[ii];
              norm = 0.5;
            }
            else
            {
              for (ii = 0; ii < n; ii++)
              {
                sdir[ii] = 1.0 * rand();
              }
              DSCAL(n, 1 / DNRM2(n, sdir, 1), sdir, incx);
              DSCAL(n, norm, sdir, incx);
            }
            DCOPY(n, z, incx, szaux, incx);
            // DSCAL(n,0.0,zaux,incx);
            /* zaux = z + dir */
            DAXPY(n , norm , sdir , 1 , szaux , 1);
            /* Computes the jacobian of the merit function, jacobian_psi_zaux = transpose(jacobianPhi_zaux).phi_zaux */
            (*phi)(n, szaux, sphi_zaux, 0);
            (*jacobianPhi)(n, szaux, sjacobianPhi_zaux, 1);

            /*norm = DNRM2(n,sphi_zaux,1);
            printf("Norm of the sol is now %e.\n",norm);*/
            /*  for (ii=0;ii<n;ii++)
              printf("zsol %e zaux %e \n",sZsol[ii],szaux[ii]);*/





            DGEMV(LA_TRANS, n, n, 1.0, sjacobianPhi_zaux, n, sphi_zaux, incx, 0.0, sjacobian_psi_zaux, incx);
            DCOPY(n, szaux, 1, szzaux, 1);
            DAXPY(n , -1 , z , incx , szzaux , incx);
            /*zzaux must be a descente direction.*/
            /*ie jacobian_psi_zaux.zzaux <0
            printf("jacobian_psi_zaux : \n");*/
            /*DCOPY(n,sdir,incx,sdir_descent,incx);
            plotMerit(z, phi);*/


            aux = DDOT(n, sjacobian_psi_zaux, 1, szzaux, 1);
            /*        aux1 = DNRM2(n,szzaux,1);
            aux1 = DNRM2(n,sjacobian_psi_zaux,1);*/
            aux = aux / (DNRM2(n, szzaux, 1) * DNRM2(n, sjacobian_psi_zaux, 1));
            printf("aux: %e\n", aux);
            if (aux < 0.1 * (j + 1))
            {
              //zaux is the new point.
              findNewZ = 1;
              DCOPY(n, szaux, incx, z, incx);
              lastN = i + 1;
              break;
            }
          }
          if (findNewZ)
            break;
          norm = -2 * norm;
        }
        if (! findNewZ)
        {
          printf("failed to find a new z\n");
          useNewZ = 0;
          /* exit(1);*/
          continue;

        }
        else
          continue;
      }
    }
    prev_norm_jacobian_psi_z = norm_jacobian_psi_z;
    prev_psi_z = psi_z;

    /* Stops if the termination criterion is satisfied */
    terminationCriterion = norm_jacobian_psi_z;
    if (terminationCriterion < tolerance)
    {
      break;
    }
    if (normPhi_z < tolerance)
    {
      terminationCriterion = tolerance / 2.0;
      break;
    }

    /* Search direction calculation
    Find a solution dk of jacobianPhiMatrix.d = -phiVector.
    dk is saved in phiVector.
    */
    DSCAL(n , -1.0 , sphi_z, incx);
    DGESV(n, 1, sjacobianPhi_z, n, sipiv, sphi_z, n, infoDGESV);
    DCOPY(n, sphi_z, 1, sdir_descent, 1);
    /*     norm = DNRM2(n,sdir_descent,1);
    printf("norm desc %e \n",norm);
    DSCAL( n , 1/norm , sdir_descent, 1);
    */
    /* descentCondition = jacobian_psi.dk */
    descentCondition = DDOT(n, sjacobian_psi_z,  1,  sdir_descent, 1);

    /* Criterion to be satisfied: error < -rho*norm(dk)^p */
    criterion = DNRM2(n, sdir_descent, 1);
    criterion = -rho * pow(criterion, p);

    if (infoDGESV != 0 || descentCondition > criterion)
    {
      /* dk = - jacobian_psi (remind that dk is saved in phi_z) */
      DCOPY(n, sjacobian_psi_z, 1, sdir_descent, 1);
      DSCAL(n , -1.0 , sdir_descent, incx);
    }

    /* Step-3 Line search: computes z_k+1 */
    /*linesearch_Armijo(n,z,sphi_z,psi_z, descentCondition, phi);*/

    lineSearch_Wolfe(z, descentCondition, phi);
    /*      for (j=20;j<32;j++){
    if (z[j]<0)
    z[j]=0;
    }*/

    /*      if( 1 || verbose>0)
    {
    printf("Non Smooth Newton, iteration number %i, error grad equal to %14.7e , psi value is %14.7e .\n",niter, terminationCriterion,psi_z);
     printf(" -----------------------------------------------------------------------\n");
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

  if (dparam[1] > tolerance)
    return 1;
  else return 0;
}
