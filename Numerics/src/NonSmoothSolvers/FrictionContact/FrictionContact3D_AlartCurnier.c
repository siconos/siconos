/* Siconos-Numerics, Copyright INRIA 2005-2010.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/
#include "LA.h"
#include "FrictionContact3D_Solvers.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "Friction_cst.h"
#include "op3x3.h"

typedef void (*computeNonsmoothFunction)(double *, double * , double , double * , double *, double *, double *);

//#define VERBOSE_DEBUG
//#define AC_STD
//#define AC_Generated
#define AC_CKPS // default is Christensen & Pang version

// Set the function for computing F and its gradient
// \todo should nbe done in initialization
#ifdef AC_STD
computeNonsmoothFunction  Function = &(computeAlartCurnierSTD);
#endif
#ifdef AC_CKPS
computeNonsmoothFunction  Function = &(computeAlartCurnierCKPS);
#endif

// computeAlartCurnierCKPS == AC_Generated (but with transpose(A) why ?)

#ifdef AC_Generated
computeNonsmoothFunction  Function = &(frictionContact3D_localAlartCurnierFunctionGenerated);
#endif

// HandMade not done
#ifdef AC_HandMade
computeNonsmoothFunction  Function = &(frictionContact3D_localAlartCurnierFunctionHandMade);
#endif


#define OPTI_RHO

/*Static variables */

/* The global problem of size n= 3*nc, nc being the number of contacts, is locally saved in MGlobal and qGlobal */
/* mu corresponds to the vector of friction coefficients */
/* note that either MGlobal or MBGlobal is used, depending on the chosen storage */
/* static int n=0; */
/* static const NumericsMatrix* MGlobal = NULL; */
/* static const double* qGlobal = NULL; */
/* static const double* mu = NULL; */

/* Local problem operators */
/* static const int nLocal = 3; */
/* static double* MLocal; */
/* static int isMAllocatedIn = 0; /\* True if a malloc is done for MLocal, else false *\/ */
static double velocityLocal[3];
/* static double qLocal[3]; */
/* static double mu_i = 0.0; */

static FrictionContactProblem* localFC3D = NULL;
static FrictionContactProblem* globalFC3D = NULL;

static double an;
static double at;
static double projN;
static double projT;
static double projS;

void AC_fillMLocal(FrictionContactProblem * problem, FrictionContactProblem * localproblem, int contact)
{

  NumericsMatrix * MGlobal = problem->M;
  int n = 3 * problem->numberOfContacts;



  // Dense storage
  int storageType = MGlobal->storageType;
  if (storageType == 0)
  {
    int in = 3 * contact, it = in + 1, is = it + 1;
    int inc = n * in;
    double * MM = MGlobal->matrix0;
    double * MLocal =  localproblem->M->matrix0;
    /* The part of MM which corresponds to the current block is copied into MLocal */
    MLocal[0] = MM[inc + in];
    MLocal[1] = MM[inc + it];
    MLocal[2] = MM[inc + is];
    inc += n;
    MLocal[3] = MM[inc + in];
    MLocal[4] = MM[inc + it];
    MLocal[5] = MM[inc + is];
    inc += n;
    MLocal[6] = MM[inc + in];
    MLocal[7] = MM[inc + it];
    MLocal[8] = MM[inc + is];
  }
  else if (storageType == 1)
  {
    int diagPos = getDiagonalBlockPos(MGlobal->matrix1, contact);
    localproblem->M->matrix0 = MGlobal->matrix1->block[diagPos];
    /*     DCOPY(9, MGlobal->matrix1->block[diagPos], 1, localproblem->M->matrix0 , 1); */

  }
  else
    numericsError("FrictionContact3D_AlartCurnier:AC_fillMLocal() -", "unknown storage type for matrix M");

}

void frictionContact3D_AC_initialize(FrictionContactProblem* problem, FrictionContactProblem* localproblem)
{
  /*
    In initialize, these operators are "connected" to their corresponding static variables, that will be used to build local problem
    for each considered contact.
    Local problem is built during call to update (which depends on the storage type for M).
  */

  localFC3D = localproblem;
  globalFC3D = problem;



}

void frictionContact3D_AC_update(int contact, FrictionContactProblem* problem, FrictionContactProblem* localproblem, double * reaction, SolverOptions* options)
{
  /* Build a local problem for a specific contact
     reaction corresponds to the global vector (size n) of the global problem.
  */
  /* Call the update function which depends on the storage for MGlobal/MBGlobal */
  /* Build a local problem for a specific contact
   reaction corresponds to the global vector (size n) of the global problem.
  */

  /* The part of MGlobal which corresponds to the current block is copied into MLocal */
  AC_fillMLocal(problem, localproblem, contact);

  /****  Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products MLocal.reactionBlock,
     excluding the block corresponding to the current contact. ****/
  frictionContact3D_nsgs_computeqLocal(problem, localproblem, reaction, contact);

  /* Friction coefficient for current block*/
  localproblem->mu[0] = problem->mu[contact];


}

/* Compute function F(Reaction) */
void F_AC(int Fsize, double *reaction , double *F, int up2Date)
{

  int nLocal = 3;



  if (F == NULL)
  {
    fprintf(stderr, "FrictionContact3D_AlartCurnier::F_AC error:  F == NULL.\n");
    exit(EXIT_FAILURE);
  }
  if (Fsize != nLocal)
  {
    fprintf(stderr, "FrictionContact3D_AlartCurnier::F error, wrong block size.\n");
    exit(EXIT_FAILURE);
  }
  double * qLocal = localFC3D->q;
  double * MLocal = localFC3D->M->matrix0;
  double mu_i = localFC3D->mu[0];


  /* up2Date = 1 = true if jacobianF(n, reaction,jacobianF) has been called just before jacobianFMatrix(...). In that case the computation of
     velocityLocal is not required again.
  */
  if (up2Date == 0)
  {
    /* velocityLocal = M.reaction + qLocal */
    int incx = 1, incy = 1;
    DCOPY(Fsize, qLocal, incx, velocityLocal, incy);
    DGEMV(LA_NOTRANS, nLocal, nLocal, 1.0, MLocal, 3, reaction, incx, 1.0, velocityLocal, incy);
    /*   velocityLocal[0] = MLocal[0]*reaction[0] + MLocal[Fsize]*reaction[1] + MLocal[2*Fsize]*reaction[2] + qLocal[0]; */
    /*   velocityLocal[1] = MLocal[1]*reaction[0] + MLocal[Fsize+1]*reaction[1] + MLocal[2*Fsize+1]*reaction[2] + qLocal[1]; */
    /*   velocityLocal[2] = MLocal[2]*reaction[0] + MLocal[Fsize+2]*reaction[1] + MLocal[2*Fsize+2]*reaction[2] + qLocal[2]; */

    an = 1. / MLocal[0];
    double alpha = MLocal[Fsize + 1] + MLocal[2 * Fsize + 2];
    double det = MLocal[1 * Fsize + 1] * MLocal[2 * Fsize + 2] - MLocal[2 * Fsize + 1] * MLocal[1 * Fsize + 2];
    double beta = alpha * alpha - 4 * det;
    if (beta > 0.)
      beta = sqrt(beta);
    else
      beta = 0.;
    at = 2 * (alpha - beta) / ((alpha + beta) * (alpha + beta));

    /* Projection on [0, +infty[ and on D(0, mu*reactionn) */
    projN = reaction[0] - an * velocityLocal[0];
    projT = reaction[1] - at * velocityLocal[1];
    projS = reaction[2] - at * velocityLocal[2];
  }

  double num;

  double coef2 = mu_i * mu_i;
  if (projN > 0)
  {
    F[0] = velocityLocal[0];
  }
  else
  {
    F[0] = reaction[0] / an;
  }

  double mrn = projT * projT + projS * projS;
  if (mrn <= coef2 * reaction[0]*reaction[0])
  {
    F[1] = velocityLocal[1];
    F[2] = velocityLocal[2];
  }
  else
  {
    num  = mu_i / sqrt(mrn);
    F[1] = (reaction[1] - projT * reaction[0] * num) / at;
    F[2] = (reaction[2] - projS * reaction[0] * num) / at;
  }
  /*   for(int l = 0; l < 3 ; ++l) */
  /*     printf(" %lf", F[l]); */
  /*   printf("\n"); */

}

/* Compute Jacobian of function F */
void jacobianF_AC(int Fsize, double *reaction, double *jacobianFMatrix, int up2Date)
{
  int nLocal = 3;

  double * qLocal = localFC3D->q;
  double * MLocal = localFC3D->M->matrix0;
  double mu_i = localFC3D->mu[0];



  if (jacobianFMatrix == NULL)
  {
    fprintf(stderr, "FrictionContact3D_AlartCurnier::jacobianF_AC error: jacobianMatrix == NULL.\n");
    exit(EXIT_FAILURE);
  }
  if (Fsize != nLocal)
  {
    fprintf(stderr, "FrictionContact3D_AlartCurnier::jacobianF_AC error, wrong block size.\n");
    exit(EXIT_FAILURE);
  }





  /* Warning: input reaction is not necessary equal to the last computed value of reactionBlock */

  /* up2Date = 1 = true if F(n, reaction,F) has been called just before jacobianFMatrix(...). In that case the computation of
     velocityLocal is not required again.
  */
  if (up2Date == 0)
  {
    /* velocityLocal = M.reaction + qLocal */
    int incx = 1, incy = 1;
    DCOPY(Fsize, qLocal, incx, velocityLocal, incy);
    DGEMV(LA_NOTRANS, nLocal, nLocal, 1.0, MLocal, 3, reaction, incx, 1.0, velocityLocal, incy);

    an = 1. / MLocal[0];
    double alpha = MLocal[Fsize + 1] + MLocal[2 * Fsize + 2];
    double det = MLocal[1 * Fsize + 1] * MLocal[2 * Fsize + 2] - MLocal[2 * Fsize + 1] * MLocal[1 * Fsize + 2];
    double beta = alpha * alpha - 4 * det;
    if (beta > 0.)
      beta = sqrt(beta);
    else
      beta = 0.;
    at = 2 * (alpha - beta) / ((alpha + beta) * (alpha + beta));

    /* Projection on [0, +infty[ and on D(0, mu*reactionn) */
    projN = reaction[0] - an * velocityLocal[0];
    projT = reaction[1] - at * velocityLocal[1];
    projS = reaction[2] - at * velocityLocal[2];
  }

  double coef2 = mu_i * mu_i;

  int i, j;
  if (projN > 0)
  {
    for (j = 0; j < Fsize; ++j)
      jacobianFMatrix[j * Fsize] = MLocal[j * Fsize];
  }
  else
  {
    jacobianFMatrix[0] = 1.0 / an;
  }

  double mrn = projT * projT + projS * projS;
  double num, rcof, mrn3;
  double tmp;
  if (mrn <= coef2 * reaction[0]*reaction[0])
    for (i = 1; i < Fsize; ++i)
      for (j = 0; j < Fsize; ++j)
        jacobianFMatrix[j * Fsize + i] = MLocal[j * Fsize + i];
  else
  {
    num  = 1. / sqrt(mrn);
    mrn3 = 1. / sqrt(mrn) * sqrt(mrn) * sqrt(mrn);
    rcof = mu_i / at;
    tmp = at * mrn3 * (MLocal[1] * projT + MLocal[2] * projS);
    jacobianFMatrix[1] = -rcof * (num * projT + reaction[0] * projT * tmp);
    jacobianFMatrix[2] = -rcof * (num * projS + reaction[0] * projS * tmp);

    tmp = mrn3 * ((1 - at * MLocal[Fsize + 1]) * projT - at * MLocal[Fsize + 2] * projS);
    jacobianFMatrix[1 * Fsize + 1] = (1 - mu_i * reaction[0] * (num * (1 - at * MLocal[Fsize + 1]) - projT * tmp)) / at;
    jacobianFMatrix[1 * Fsize + 2] =  - rcof * reaction[0] * ((-num * at * MLocal[Fsize + 2]) - projS * tmp);

    tmp = mrn3 * ((1 - at * MLocal[2 * Fsize + 2]) * projS - at * MLocal[2 * Fsize + 1] * projT);
    jacobianFMatrix[2 * Fsize + 1] =  - rcof * reaction[0] * ((-num * at * MLocal[2 * Fsize + 1]) - projT * tmp);
    jacobianFMatrix[2 * Fsize + 2] = (1 - mu_i * reaction[0] * (num * (1 - at * MLocal[2 * Fsize + 2]) - projS * tmp)) / at;
  }
  /*   for(int l = 0; l < 9 ; ++l) */
  /*     printf(" %lf", jacobianFMatrix[l]); */
  /*   printf("\n"); */
}

void frictionContact3D_AC_post(int contact, double* reaction)
{
  /* This function is required in the interface but useless in Alart-Curnier case */
}

void computeFGlobal_AC(double* reaction, double* FGlobal)
{

  int numberOfContacts =  globalFC3D->numberOfContacts;

  int n = 3 * numberOfContacts;

  NumericsMatrix * MGlobal = globalFC3D->M;
  double * MLocal = localFC3D->M->matrix0;
  double * qLocal = localFC3D->q;
  double *mu = globalFC3D->mu;


  int contact, diagPos = 0, position;
  int in, it, is, inc, incx;
  double * reactionLocal;
  double alpha, det, beta, num, coef2, mrn;
  for (contact = 0; contact < numberOfContacts; ++contact)
  {
    position = 3 * contact;
    if (MGlobal->storageType == 1) /* Sparse storage */
    {
      /* The part of MGlobal which corresponds to the current block is copied into MLocal */
      diagPos = numberOfContacts * contact + contact;
      MLocal = MGlobal->matrix1->block[diagPos];
    }
    else if (MGlobal->storageType == 0)
    {
      in = 3 * contact;
      it = in + 1;
      is = it + 1;
      inc = n * in;
      double *MM = MGlobal->matrix0;
      /* The part of MM which corresponds to the current block is copied into MLocal */
      MLocal[0] = MM[inc + in];
      MLocal[1] = MM[inc + it];
      MLocal[2] = MM[inc + is];
      inc += n;
      MLocal[3] = MM[inc + in];
      MLocal[4] = MM[inc + it];
      MLocal[5] = MM[inc + is];
      inc += n;
      MLocal[6] = MM[inc + in];
      MLocal[7] = MM[inc + it];
      MLocal[8] = MM[inc + is];
    }

    reactionLocal = &reaction[3 * contact];
    incx = 3;
    velocityLocal[0] = DDOT(3 , MLocal , incx , reactionLocal , 1) + qLocal[0];
    velocityLocal[1] = DDOT(3 , MLocal , incx , reactionLocal , 1) + qLocal[1];
    velocityLocal[2] = DDOT(3 , MLocal , incx , reactionLocal , 1) + qLocal[2];
    an = 1. / MLocal[0];
    alpha = MLocal[4] + MLocal[8];
    det = MLocal[4] * MLocal[8] - MLocal[7] + MLocal[5];
    beta = alpha * alpha - 4 * det;
    at = 2 * (alpha - beta) / ((alpha + beta) * (alpha + beta));
    projN = reactionLocal[0] - an * velocityLocal[0];
    projT = reactionLocal[1] - at * velocityLocal[1];
    projS = reactionLocal[2] - at * velocityLocal[2];
    coef2 = mu[contact] * mu[contact];
    if (projN > 0)
    {
      FGlobal[position] = velocityLocal[0];
    }
    else
    {
      FGlobal[position] = reactionLocal[0] / an;
    }

    mrn = projT * projT + projS * projS;
    if (mrn <= coef2 * reactionLocal[0]*reactionLocal[0])
    {
      FGlobal[position + 1] = velocityLocal[1];
      FGlobal[position + 2] = velocityLocal[2];
    }
    else
    {
      num  = mu[contact] / sqrt(mrn);
      FGlobal[position + 1] = (reactionLocal[1] - projT * reactionLocal[0] * num) / at;
      FGlobal[position + 2] = (reactionLocal[2] - projS * reactionLocal[0] * num) / at;
    }
  }
}

void frictionContact3D_AC_free()
{

}

int frictionContact3D_AlartCurnierNewton_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the NSGS Solver\n");
  }

  options->solverId = SICONOS_FRICTION_3D_DampedAlartCurnierNewton;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  options->iWork = NULL;
  for (i = 0; i < 5; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }

  options->iparam[0] = 10;
  options->iparam[1] = 10;
  options->dparam[0] = 1e-16;
  return 0;
}

/* Alart & Curnier version (Radius = mu*max(0,RVN)) */
void computeAlartCurnierSTD(double R[3], double velocity[3], double mu, double rho[3], double F[3], double A[9], double B[9])
{
  double RVN, RVT, RVS;
  double RhoN = rho[0];
  double RhoT = rho[1];
  double Radius, RV, RV1, RV3, GammaTT, GammaTS, GammaST, GammaSS;

  RVN = R[0] - RhoN * velocity[0];
  RVT = R[1] - RhoT * velocity[1];
  RVS = R[2] - RhoT * velocity[2];
  RV = sqrt(RVT * RVT + RVS * RVS);
  //Radius = mu*R[0];


  if (RVN >= 0.0)
  {
#ifdef VERBOSE_DEBUG
    printf("Normal part in the cone\n");
#endif


    Radius = mu * RVN;
    F[0] = RhoN * (velocity[0]);
    if (A && B)
    {
      A[0 + 3 * 0] =  RhoN;
      B[0 + 3 * 0] = 0.0;
    }
  }
  else
  {
#ifdef VERBOSE_DEBUG
    printf("Normal part out the cone\n");
#endif
    Radius = 0.0;
    F[0] = R[0];
    if (A && B)
    {
      A[0 + 3 * 0] = 0.0;
      B[0 + 3 * 0] = 1.0;
    }
  }

  // Compute the value of the Alart--Curnier Function and its gradient for the tangential part


#ifdef VERBOSE_DEBUG
  printf("Radius=%le\n", Radius);
  printf("RV=%le\n", RV);
#endif
  if (RV <= Radius) // We are in the disk and Radius is positive
  {
#ifdef VERBOSE_DEBUG
    printf("We are in the disk\n");
#endif
    F[1] = RhoT * (velocity[1]);
    F[2] = RhoT * (velocity[2]);
    if (A && B)
    {
      A[1 + 3 * 1] = RhoT;
      A[1 + 3 * 2] = 0.0;
      A[2 + 3 * 1] = 0.0;
      A[2 + 3 * 2] = RhoT;
      B[1 + 3 * 0] = 0.0;
      B[1 + 3 * 1] = 0.0;
      B[1 + 3 * 2] = 0.0;
      B[2 + 3 * 0] = 0.0;
      B[2 + 3 * 1] = 0.0;
      B[2 + 3 * 2] = 0.0;
    }
  }
  else if (RV > Radius) // We are out the disk and Radius is postive
  {

    if (Radius > 0)
    {
#ifdef VERBOSE_DEBUG
      printf("We are out the disk and Radius is positive\n");
#endif
      RV1 = 1.0 / RV;
      F[1] = R[1] - Radius * RVT * RV1;
      F[2] = R[2] - Radius * RVS * RV1;
      if (A && B)
      {
        RV3 = RV1 * RV1 * RV1;
        GammaTT = RV1 - RVT * RVT * RV3;
        GammaTS =  - RVT * RVS * RV3;
        GammaST =  GammaTS;
        GammaSS = RV1 - RVS * RVS * RV3;

        A[0 + 3 * 1] = mu * RhoN * RVT * RV1;
        A[0 + 3 * 2] = mu * RhoN * RVS * RV1;


        A[1 + 3 * 1] = GammaTT * RhoT * Radius;

        A[1 + 3 * 2] = GammaTS * RhoT * Radius;
        A[2 + 3 * 1] = GammaST * RhoT * Radius;

        A[2 + 3 * 2] = GammaSS * RhoT * Radius;

        B[1 + 3 * 0] = -mu * RVT * RV1;

        B[1 + 3 * 1] = 1.0 - GammaTT * Radius ;
        B[1 + 3 * 2] = - GammaTS  * Radius ;

        B[2 + 3 * 0] = -mu * RVS * RV1;

        B[2 + 3 * 1] = - GammaST  * Radius;
        B[2 + 3 * 2] = 1.0 - GammaSS * Radius;
      }
    }
    else
#ifdef VERBOSE_DEBUG
      printf("We are out the disk and Radius is zero\n");
#endif
    {

      F[1] = R[1] ;
      F[2] = R[2] ;
      if (A && B)
      {
        A[1 + 3 * 1] = 0.0;
        A[1 + 3 * 2] = 0.0;
        A[2 + 3 * 1] = 0.0;
        A[2 + 3 * 2] = 0.0;

        B[1 + 3 * 0] = 0.0;
        B[1 + 3 * 1] = 1.0;
        B[1 + 3 * 2] = 0.0;
        B[2 + 3 * 0] = 0.0;
        B[2 + 3 * 1] = 0.0;
        B[2 + 3 * 2] = 1.0;
      }
    }

  }
  /*   else // We are out the disk and Radius is negative */
  /*     { */
  /* #ifdef VERBOSE_DEBUG */
  /*       printf("We are out the disk and Radius is negative\n"); */
  /* #endif */

  /*       /\*Version original *\/ */
  /*       F[1] = R[1] ; */
  /*       F[2] = R[2] ; */
  /*       if (A && B){ */
  /*    A[1+3*1]=0.0; */
  /*    A[1+3*2]=0.0; */
  /*    A[2+3*1]=0.0; */
  /*    A[2+3*2]=0.0; */

  /*    B[1+3*0]=0.0; */
  /*    B[1+3*1]=1.0; */
  /*    B[1+3*2]=0.0; */
  /*    B[2+3*0]=0.0; */
  /*    B[2+3*1]=0.0; */
  /*    B[2+3*2]=1.0;} */
  /*     } */

#ifdef VERBOSE_DEBUG
  printf("F[0] = %le\t", F[0]);
  printf("F[1] = %le\t", F[1]);
  printf("F[2] = %le\n", F[2]);

  if (A && B)
  {
    for (int l = 0; l < 3; l++)
    {
      for (int k = 0; k < 3; k++)
      {
        printf("A[%i+3*%i] = %le\t", l, k, A[l + 3 * k]);
      }
      printf("\n");
    }
    for (int l = 0; l < 3; l++)
    {
      for (int k = 0; k < 3; k++)
      {
        printf("B[%i+3*%i] = %le\t", l, k, B[l + 3 * k]);
      }
      printf("\n");
    }
    double diago = 0.0;
    for (int l = 0; l < 3; l++)
    {
      for (int k = 0; k < 3; k++)
      {
        if (k == l)  diago = 1.0;
        else diago = 0.0;
        printf("I-B[%i+3*%i] = %le\t", l, k, diago - B[l + 3 * k]);
      }
      printf("\n");
    }
  }
#endif



}


/* Christensen & Pang version (Radius = mu* R[0])*/
void computeAlartCurnierCKPS(double R[3], double velocity[3], double mu, double rho[3], double F[3], double A[9], double B[9])
{

  double RVN, RVT, RVS;
  double RhoN = rho[0];
  double RhoT = rho[1];
  double Radius, RV, RV1, RV3, GammaTT, GammaTS, GammaST, GammaSS;

  RVN = R[0] - RhoN * velocity[0];
  RVT = R[1] - RhoT * velocity[1];
  RVS = R[2] - RhoT * velocity[2];
  RV = sqrt(RVT * RVT + RVS * RVS);
  Radius = mu * R[0];

  // Compute the value of the Alart--Curnier Function and its gradient for the normal part

  if (RVN >= 0.0)
  {
    F[0] = RhoN * (velocity[0]);
    if (A && B)
    {
      A[0 + 3 * 0] =  RhoN;
      B[0 + 3 * 0] = 0.0;
    }
  }
  else
  {
    F[0] = R[0];
    if (A && B)
    {
      A[0 + 3 * 0] = 0.0;
      B[0 + 3 * 0] = 1.0;
    }
  }

  // Compute the value of the Alart--Curnier Function and its gradient for the tangential part


#ifdef VERBOSE_DEBUG
  printf("Radius=%le\n", Radius);
  printf("RV=%le\n", RV);
#endif
  if (RV < Radius || RV < 1e-20)  // We are in the disk
  {
#ifdef VERBOSE_DEBUG
    printf("We are in the disk \n");
#endif
    F[1] = RhoT * (velocity[1]);
    F[2] = RhoT * (velocity[2]);
    if (A && B)
    {
      A[1 + 3 * 1] = RhoT;
      A[1 + 3 * 2] = 0.0;
      A[2 + 3 * 1] = 0.0;
      A[2 + 3 * 2] = RhoT;
      B[1 + 3 * 0] = 0.0;
      B[1 + 3 * 1] = 0.0;
      B[1 + 3 * 2] = 0.0;
      B[2 + 3 * 0] = 0.0;
      B[2 + 3 * 1] = 0.0;
      B[2 + 3 * 2] = 0.0;
    }
  }
  else  // We are out the disk
  {
#ifdef VERBOSE_DEBUG
    printf("We are out the disk\n");
#endif
    /*        RV1 = 1.0/RV; */
    /*        F[1] = R[1] - Radius*RVT*RV1; */
    /*        F[2] = R[2] - Radius*RVS*RV1; */


    /*        RV3 = RV1*RV1*RV1; */
    /*        GammaTT = (RV - RVT*RVT)*RV3; */
    /*        GammaTS =  - RVT*RVS*RV3; */
    /*        GammaST =  GammaTS; */
    /*        GammaSS = (RV - RVS*RVS)*RV3; */


    RV1 = 1.0 / RV;
    F[1] = R[1] - Radius * RVT * RV1;
    F[2] = R[2] - Radius * RVS * RV1;
    if (A && B)
    {
      RV3 = RV1 * RV1 * RV1;
      GammaTT = RV1 - RVT * RVT * RV3;
      GammaTS =  - RVT * RVS * RV3;
      GammaST =  GammaTS;
      GammaSS = RV1 - RVS * RVS * RV3;

      // come back to r2146
      //        A[0+3*1]= mu * RhoN * RVT*RV1;
      //        A[0+3*2]= mu * RhoN * RVS*RV1;

      A[1 + 3 * 1] = GammaTT * RhoT * Radius;

      A[1 + 3 * 2] = GammaTS * RhoT * Radius;

      A[2 + 3 * 1] = GammaST * RhoT * Radius;
      A[2 + 3 * 2] = GammaSS * RhoT * Radius;

      B[1 + 3 * 0] = -mu * RVT * RV1;

      B[1 + 3 * 1] = 1.0 - GammaTT * Radius ;
      B[1 + 3 * 2] = - GammaTS * Radius ;

      B[2 + 3 * 0] = -mu * RVS * RV1;

      B[2 + 3 * 1] = - GammaST * Radius;
      B[2 + 3 * 2] = 1.0 - GammaSS * Radius;
    }
  }




#ifdef VERBOSE_DEBUG
  printf("F[0] = %le\n", F[0]);
  printf("F[1] = %le\n", F[1]);
  printf("F[2] = %le\n", F[2]);
  if (A && B)
  {
    for (int l = 0; l < 3; l++)
    {
      for (int k = 0; k < 3; k++)
      {
        printf("A[%i+3*%i] = %le\t", l, k, A[l + 3 * k]);
      }
      printf("\n");
    }
    for (int l = 0; l < 3; l++)
    {
      for (int k = 0; k < 3; k++)
      {
        printf("B[%i+3*%i] = %le\t", l, k, B[l + 3 * k]);
      }
      printf("\n");
    }
  }

#endif
}


void computerho(FrictionContactProblem* localproblem, double * rho)
{




  double * MLocal = localproblem->M->matrix0;

  assert(MLocal[0 + 0 * 3] > 0);

  double sw = MLocal[1 + 1 * 3] + MLocal[2 + 2 * 3];

  double dw = sw * sw - 4.0 * (MLocal[1 + 1 * 3] + MLocal[2 + 2 * 3] -  MLocal[2 + 1 * 3] + MLocal[1 + 2 * 3]);

  if (dw > 0.0) dw = sqrt(dw);
  else dw = 0.0;

  rho[0] = 1.0 / MLocal[0 + 0 * 3];
  rho[1] = 2.0 * (sw - dw) / ((sw + dw) * (sw + dw));
  rho[2] = 2.0 * (sw - dw) / ((sw + dw) * (sw + dw));

  assert(rho[0] > 0);
  assert(rho[1] > 0);
  assert(rho[2] > 0);




#ifdef VERBOSE_DEBUG
  printf("sw=%le = ", sw);
  printf("dw=%le = ", dw);
  printf("rho[0]=%le = ", rho[0]);
  printf("rho[1]=%le\n", rho[1]);
  printf("rho[2]=%le\n", rho[2]);
  printf("\n");
#endif
}




int LocalNonsmoothNewtonSolver(FrictionContactProblem* localproblem, double * R, int *iparam, double *dparam)
{
  double mu = localproblem->mu[0];
  double * qLocal = localproblem->q;

  double * MLocal = localproblem->M->matrix0;





  double Tol = dparam[0];
  double itermax = iparam[0];


  int i, j, k, inew;

  // store the increment
  double dR[3] = {0., 0., 0.};

  // store the value fo the function
  double F[3] = {0., 0., 0.};

  // Store the (sub)-gradient of the function
  double A[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
  double B[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};

  // Value of AW+B
  double AWplusB[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};

  // Compute values of Rho (should be here ?)
  double rho[3] = {1., 1., 1.};
#ifdef OPTI_RHO
  computerho(localproblem, rho);
#endif

  // compute the velocity
  double velocity[3] = {0., 0., 0.};

  for (i = 0; i < 3; i++) velocity[i] = MLocal[i + 0 * 3] * R[0] + qLocal[i]
                                          + MLocal[i + 1 * 3] * R[1] +
                                          + MLocal[i + 2 * 3] * R[2] ;

  for (inew = 0 ; inew < itermax ; ++inew) // Newton iteration
  {
    //Update function and gradient

    Function(R, velocity, mu, rho, F, A, B);

#ifndef NDEBUG
    double Fg[3] = {0., 0., 0.};
    double Ag[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
    double Bg[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
    double AWpB[9];


    assert(*rho > 0. && *(rho + 1) > 0. && *(rho + 2) > 0.);

#ifdef AC_STD
    frictionContact3D_localAlartCurnierFunctionGenerated(R, velocity, mu, rho, Fg, Ag, Bg);
#endif

#ifdef AC_CKPS
    frictionContact3D_localAlartCurnierCKPSFunctionGenerated(R, velocity, mu, rho, Fg, Ag, Bg);
#endif

    sub3x3(F, Fg);
    sub3x3(A, Ag);
    sub3x3(B, Bg);

    assert(hypot3(Fg) <= 1e-7);
    assert(hypot9(Ag) <= 1e-7);
    assert(hypot9(Bg) <= 1e-7);
    cpy3x3(A, Ag);
    cpy3x3(B, Bg);
    mm3x3(A, MLocal, AWpB);
    add3x3(B, AWpB);

#endif



    // compute -(A MLocal +B)
    for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
      {
        AWplusB[i + 3 * j] = 0.0;
        for (k = 0; k < 3; k++)
        {
          AWplusB[i + 3 * j] -= A[i + 3 * k] * MLocal[k + j * 3];
        }
        AWplusB[i + 3 * j] -= B[i + 3 * j];
      }
    }

#ifdef AC_STD
#ifndef NDEBUG
    scal3x3(-1., AWpB);
    sub3x3(AWplusB, AWpB);
    assert(hypot9(AWpB) <= 1e-7);
#endif
#endif

    // Solve the linear system
    solv3x3(AWplusB, dR, F);
    // upate iterates
    R[0] += dR[0];
    R[1] += dR[1];
    R[2] += dR[2];
    // compute new residue
    for (i = 0; i < 3; i++) velocity[i] = MLocal[i + 0 * 3] * R[0] + qLocal[i]
                                            + MLocal[i + 1 * 3] * R[1] +
                                            + MLocal[i + 2 * 3] * R[2] ;
    //      Function(R,velocity,mu,rho,F,NULL,NULL);
    //      dparam[1] = 0.5 *(F[0]*F[0]+F[1]*F[1]+F[2]*F[2])/(1.0 + sqrt(R[0]*R[0]+R[1]*R[1]+R[2]*R[2]) ) ; // improve with relative tolerance

    /*      dparam[2] =0.0;  */
    /*      FrictionContact3D_unitary_compute_and_add_error( R , velocity,mu, &(dparam[2])); */




    if (verbose > 1) printf("-----------------------------------    LocalNewtonSolver number of iteration = %i  error = %.10e \n", inew, dparam[1]);

    if (dparam[1] < Tol)
    {
      /*    printf("-----------------------------------    LocalNewtonSolver number of iteration = %i  error = %.10e \t error2 = %.10e \n",inew,dparam[1], dparam[2]); */

      return 0;
    }

  }// End of the Newton iteration

  /*  printf("-----------------------------------    LocalNewtonSolver number of iteration = %i  error = %.10e \t error2 = %.10e \n",inew,dparam[1], dparam[2]); */
  return 1;

}



int  LineSearchGP(FrictionContactProblem* localproblem,
                  computeNonsmoothFunction  Function,
                  double * t_opt,
                  double R[3],
                  double dR[3],
                  double *rho,
                  int LSitermax,
                  double * F,
                  double * A,
                  double * B,
                  double * velocity)
{
  double alpha = *t_opt;

  double inf = 1e20;

  double alphamin = 0.0;
  double alphamax = inf;

  double m1 = 0.1, m2 = 0.9;


  /*     // store the value fo the function */
  /*     double F[3]={0.,0.,0.}; */

  /*     // Store the (sub)-gradient of the function */
  /*     double A[9]={0.,0.,0.,0.,0.,0.,0.,0.,0.}; */
  /*     double B[9]={0.,0.,0.,0.,0.,0.,0.,0.,0.}; */

  /*     double velocity[3]={0.,0.,0.}; */

  double mu = localproblem->mu[0];
  double * qLocal = localproblem->q;
  double * MLocal = localproblem->M->matrix0;

  /*     for (int i=0; i<3; i++) velocity[i] = MLocal[i+0*3]*R[0] + qLocal[i] */
  /*          + MLocal[i+1*3]*R[1] + */
  /*          + MLocal[i+2*3]*R[2] ; */

  /*     Function(R,velocity,mu,rho,F,A,B); */


  // Computation of q(t) and q'(t) for t =0

  double q0 = 0.5 * DDOT(3 , F , 1 , F , 1);

  double tmp[3] = {0., 0., 0.};

  // Value of AW+B
  double AWplusB[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};

  // compute A MLocal +B
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      AWplusB[i + 3 * j] = 0.0;
      for (int k = 0; k < 3; k++)
      {
        AWplusB[i + 3 * j] += A[i + 3 * k] * MLocal[k + j * 3];
      }
      AWplusB[i + 3 * j] += B[i + 3 * j];
    }
  }

#ifdef VERBOSE_DEBUG
  for (int l = 0; l < 3; l++)
  {
    for (int k = 0; k < 3; k++)
    {
      printf("AWplusB[%i+3*%i] = %le\t", l, k, AWplusB[l + 3 * k]);
    }
    printf("\n");
  }
#endif

  for (int i = 0; i < 3; i++)
  {
    tmp[i] = 0.0;
    for (int j = 0; j < 3; j++)
    {
      tmp[i] += AWplusB[i + 3 * j] * dR[j]  ;
    }
  }




  double dqdt0 = 0.0;
  for (int i = 0; i < 3; i++)
  {
    dqdt0 += F[i] * tmp[i];
  }
#ifdef VERBOSE_DEBUG
  printf("q0 = %12.8e \n", q0);
  printf("dqdt0 = %12.8e \n", dqdt0);
  for (int i = 0; i < 3; i++)
  {
    printf("tmp[%i] = %12.8e \t", i, tmp[i]);
  }
  printf("\n");
  for (int i = 0; i < 3; i++)
  {
    printf("dR[%i] = %12.8e \t", i, dR[i]);
  }
  printf("\n");
#endif

  for (int iter = 0; iter < LSitermax; iter++)
  {

    for (int i = 0; i < 3; i++)  tmp[i] = R[i] + alpha * dR[i];

    for (int i = 0; i < 3; i++) velocity[i] = MLocal[i + 0 * 3] * tmp[0] + qLocal[i]
          + MLocal[i + 1 * 3] * tmp[1] +
          + MLocal[i + 2 * 3] * tmp[2] ;

    Function(tmp, velocity, mu, rho, F, NULL, NULL);

    double q  = 0.5 * DDOT(3 , F , 1 , F , 1);

    double slope = (q - q0) / alpha;

#ifdef VERBOSE_DEBUG
    printf("q = %12.8e \n", q);
    printf("slope = %12.8e \n", slope);
#endif


    int C1 = (slope >= m2 * dqdt0);
    int C2 = (slope <= m1 * dqdt0);

    if (C1 && C2)
    {
#ifdef VERBOSE_DEBUG
      printf("Sucess in LS: alpha = %12.8e\n", alpha);
#endif
      *t_opt = alpha;
      if (verbose > 1)
      {
        printf("-----------------------------------------    LineSearchGP success number of iteration = %i  alpha = %.10e \n", iter, alpha);
      }
      return 0;

    }
    else if (!C1)
    {
#ifdef VERBOSE_DEBUG
      printf("LS: alpha too small = %12.8e\t, slope =%12.8e\n", alpha, slope);
      printf(" m1*dqdt0 =%12.8e\t, m2*dqdt0 =%12.8e\n ", m1 * dqdt0 , m2 * dqdt0);
#endif
      //std::cout << "t = " << t << " is too small : slope = " << slope << ", m2*qp0 = " << m2*qp0 << std::endl;
      alphamin = alpha;
    }
    else     // not(C2)
    {
#ifdef VERBOSE_DEBUG
      printf("LS: alpha too big = %12.8e\t, slope =%12.8e\n", alpha, slope);
      printf(" m1*dqdt0 =%12.8e\t, m2*dqdt0 =%12.8e\n ", m1 * dqdt0 , m2 * dqdt0);
#endif
      //std::cout << "t = " << t << " is too big : slope = " << slope << ", m1*qp0 = " << m1*qp0 << std::endl;
      alphamax = alpha;
    }
    if (alpha < inf)
    {
      alpha = 0.5 * (alphamin + alphamax);
    }
    else
    {
      alpha = 10 * alpha;
    }


  }
  if (verbose > 1)
  {
    printf("-----------------------------------------    LineSearchGP failed max number of iteration reached  = %i  with alpha = %.10e \n", LSitermax, alpha);
  }
  *t_opt = alpha;
  return -1;
}



int DampedLocalNonsmoothNewtonSolver(FrictionContactProblem* localproblem, double * R, int *iparam, double *dparam)
{
  double mu = localproblem->mu[0];
  double * qLocal = localproblem->q;
  double * MLocal = localproblem->M->matrix0;

  double Tol = dparam[0];
  double itermax = iparam[0];
  double LSitermax = iparam[1];


  int i, j, k, inew;

  // store the value fo the function
  double F[3] = {0., 0., 0.};

  // Store the (sub)-gradient of the function
  double A[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
  double B[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};

  // store the serach direction
  double dR[3] = {0., 0., 0.};

  // path length
  double t = 1.;
  double t_opt = 1.;
  double t_init = 1.;
  int NumberofLSfailed = 0;

  // Value of AW+B
  double AWplusB[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};

  // Compute values of Rho (should be here ?)
  double rho[3] = {1., 1., 1.};
#ifdef OPTI_RHO
  computerho(localproblem, rho);
#endif

  // compute the velocity
  double velocity[3] = {0., 0., 0.};

  //cpy3(qLocal,velocity);
  //mvp3x3(MLocal,velocity)

  for (i = 0; i < 3; i++) velocity[i] = MLocal[i + 0 * 3] * R[0] + qLocal[i]
                                          + MLocal[i + 1 * 3] * R[1] +
                                          + MLocal[i + 2 * 3] * R[2] ;

  for (inew = 0 ; inew < itermax ; ++inew) // Newton iteration
  {
    //Update function and gradient
    Function(R, velocity, mu, rho, F, A, B);

    // compute -(A MLocal +B)
    for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
      {
        AWplusB[i + 3 * j] = 0.0;
        for (k = 0; k < 3; k++)
        {
          AWplusB[i + 3 * j] -= A[i + 3 * k] * MLocal[k + j * 3];
        }
        AWplusB[i + 3 * j] -= B[i + 3 * j];
      }
    }

    solv3x3(AWplusB, dR, F);

    // Perform Line Search

    t_opt = t_init;
    int infoLS = LineSearchGP(localproblem, Function, &t_opt, R, dR, rho, LSitermax, F, A, B, velocity);

    if (infoLS == 0) t = t_opt;
    else
    {
      NumberofLSfailed++;
      if (NumberofLSfailed > 5)
      {
        t = 100.0;
        if (verbose > 1) printf("-----------------------------------------  Max Number of LineSearchGP failed =%i Tilt point\n ", NumberofLSfailed);
        NumberofLSfailed = 0;
      }
    }

    // upate iterates
    R[0] = R[0] + t * dR[0];
    R[1] = R[1] + t * dR[1];
    R[2] = R[2] + t * dR[2];

    // compute new residue
    for (i = 0; i < 3; i++) velocity[i] = MLocal[i + 0 * 3] * R[0] + qLocal[i]
                                            + MLocal[i + 1 * 3] * R[1] +
                                            + MLocal[i + 2 * 3] * R[2] ;

    Function(R, velocity, mu, rho, F, NULL, NULL);
    dparam[1] = 0.5 * (F[0] * F[0] + F[1] * F[1] + F[2] * F[2]) / (1.0 + sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2])) ; // improve with relative tolerance

    if (verbose > 1) printf("-----------------------------------    LocalNewtonSolver number of iteration = %i  error = %.10e \n", inew, dparam[1]);
    if (dparam[1] < Tol) return 0;


  }// End of the Newton iteration


  return 1;

}
