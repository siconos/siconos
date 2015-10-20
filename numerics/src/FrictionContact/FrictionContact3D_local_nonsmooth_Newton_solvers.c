/* Siconos-Numerics, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
2 * the Free Software Foundation; either version 2 of the License, or
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

#include "NonSmoothNewton.h"
#include "NCP_Solvers.h"
#include "FrictionContact3D_Solvers.h"
#include <stdlib.h>
#include <stdio.h>
#include "FischerBurmeister.h"
#include "SiconosBlas.h"
#include "AlartCurnierGenerated.h"
#include "FrictionContact3D_GlockerFischerBurmeister_functions.h"
#include "op3x3.h"


static computeNonsmoothFunction  Function = NULL;
static NewtonFunctionPtr F = NULL;
static NewtonFunctionPtr jacobianF = NULL;
static UpdateSolverPtr updateSolver = NULL;
static PostSolverPtr postSolver = NULL;
static FreeSolverPtr freeSolver = NULL;

/* size of a block */
static int Fsize;

void frictionContact3D_AC_free()
{

}
void frictionContact3D_AC_post(int contact, double* reaction)
{
  /* This function is required in the interface but useless in Alart-Curnier case */
}

void frictionContact3D_local_nonsmooth_Newton_solvers_initialize(FrictionContactProblem* problem, FrictionContactProblem* localproblem,    SolverOptions * localsolver_options)
{

  /*
     Initialize solver (Connect F and its jacobian, set local size ...) according to the chosen formulation.
  */

  /* Alart-Curnier formulation */
  if (localsolver_options->solverId == SICONOS_FRICTION_3D_AlartCurnierNewton)
  {
    frictionContact3D_AC_initialize(problem, localproblem,localsolver_options);
    /* Fsize = 3; */
    /* F = &F_AC; */
    /* jacobianF = &jacobianF_AC; */
    /*     updateSolver = &frictionContact3D_AC_update; */
    postSolver = &frictionContact3D_AC_post;
    freeSolver = (FreeSolverPtr)&frictionContact3D_AC_free;

  }
  else if (localsolver_options->solverId == SICONOS_FRICTION_3D_DampedAlartCurnierNewton)
  {
    frictionContact3D_AC_initialize(problem, localproblem,localsolver_options);
    /* Fsize = 3; */
    /* F = &F_AC; */
    /* jacobianF = &jacobianF_AC; */
    /*     updateSolver = &frictionContact3D_AC_update; */
    postSolver = &frictionContact3D_AC_post;
    freeSolver = (FreeSolverPtr)&frictionContact3D_AC_free;

  }




  /* Glocker formulation - Fischer-Burmeister function used in Newton */
  else if (localsolver_options->solverId == SICONOS_FRICTION_3D_NCPGlockerFBNewton)
  {
    Fsize = 5;
    NCPGlocker_initialize(problem, localproblem);
    F = &F_GlockerFischerBurmeister;
    jacobianF = &jacobianF_GlockerFischerBurmeister;
    /*     updateSolver = &NCPGlocker_update; */
    postSolver = &NCPGlocker_post;
    freeSolver = (FreeSolverPtr)&NCPGlocker_free;
  }
  else
  {
    fprintf(stderr, "Numerics, FrictionContact3D_nsgs failed. Unknown formulation type.\n");
    exit(EXIT_FAILURE);
  }
}

int frictionContact3D_local_nonsmooth_Newton_solvers_solve(FrictionContactProblem* localproblem, double* reaction, SolverOptions * options)
{

  /*  (*updateSolver)(contact, reaction); */
  double * reactionBlock = reaction;

  int * iparam = options->iparam;
  double * dparam = options->dparam;

  int info;
  if (options->solverId == SICONOS_FRICTION_3D_AlartCurnierNewton)
  {
    info = LocalNonsmoothNewtonSolver(localproblem, reactionBlock, iparam, dparam);
  }
  else if (options->solverId == SICONOS_FRICTION_3D_DampedAlartCurnierNewton)
  {
    info = DampedLocalNonsmoothNewtonSolver(localproblem, reactionBlock, iparam, dparam);
  }
  else
  {
    info = nonSmoothDirectNewton(Fsize, reactionBlock, &F, &jacobianF, iparam, dparam);
  }
  if (info > 0)
  {
    if (verbose > 0)
    {
      printf("Numerics, frictionContact3D_local_nonsmooth_Newton_solvers_solve, warning. reached max. number of iterations without convergence. Error = %12.8e\n", dparam[1]);
      /* note : exit on failure should be done in DefaultCheckSolverOutput */
    }
  }
  return info;
  /*  (*postSolver)(contact,reaction); */
}

void frictionContact3D_local_nonsmooth_Newton_solvers_free(FrictionContactProblem* localproblem)
{
  F = NULL;
  jacobianF = NULL;
  updateSolver = NULL;
  postSolver = NULL;
  (*freeSolver)();
}

void frictionContact3D_local_nonsmooth_Newton_solvers_computeError(int n, double* velocity, double*reaction, double * error)
{
  /*   int numberOfContacts = n/3; */
  /*   int sizeGlobal = numberOfContacts*FSize; */
  /*   //  double * FGlobal = (double*)malloc(sizeGlobal*sizeof(*FGlobal));  */
  /*   (*computeFGlobal)(reaction,velocity); */
  /*   int i; */
  /*   double Fz; */
  /*   *error = 0; */
  /*   for(i=0;i<sizeGlobal;++i) */
  /*     { */
  /*       Fz = velocity[i]*reaction[i]; */
  /*       if(Fz>0) */
  /*  *error+=Fz; */
  /*       if(reaction[i]<0) */
  /*  *error+=reaction[i]; */
  /*       if(velocity[i]<0) */
  /*  *error+=velocity[i]; */
  /*     } */

  /*   // (*computeVelocity)(FGlobal); */

  /*   free(FGlobal); */

}


static void AC_fillMLocal(FrictionContactProblem * problem, FrictionContactProblem * localproblem, int contact)
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
    /*     cblas_dcopy(9, MGlobal->matrix1->block[diagPos], 1, localproblem->M->matrix0 , 1); */

  }
  else
    numericsError("FrictionContact3D_AlartCurnier:AC_fillMLocal() -", "unknown storage type for matrix M");

}

void frictionContact3D_AC_initialize(FrictionContactProblem* problem, FrictionContactProblem* localproblem, SolverOptions * options)
{
  /*
    In initialize, these operators are "connected" to their corresponding static variables, that will be used to build local problem
    for each considered contact.
    Local problem is built during call to update (which depends on the storage type for M).
  */

  /* localFC3D = localproblem; */
  /* globalFC3D = problem; */

  if (options->iparam[10] == 0 )
  {
    Function = &(computeAlartCurnierSTD);
  }
  else if (options->iparam[10] == 1 )
  {
    Function = &(computeAlartCurnierJeanMoreau);
  }
  else if (options->iparam[10] == 2 )
  {
    Function = &(frictionContact3D_AlartCurnierFunctionGenerated);
  }
  else if (options->iparam[10] == 3 )
  {
    Function = &frictionContact3D_AlartCurnierJeanMoreauFunctionGenerated;;
  }
  else if (options->iparam[10] == 4 )
  {
    Function = NULL;
  }


/* #ifdef AC_STD */
/* computeNonsmoothFunction  Function = &(computeAlartCurnierSTD); */
/* #endif */
/* #ifdef AC_JeanMoreau */
/* computeNonsmoothFunction  Function = &(computeAlartCurnierJeanMoreau); */
/* #endif */

/* // computeAlartCurnier[JeanMoreau] == AC_Generated */

/* #ifdef AC_Generated */
/* computeNonsmoothFunction  Function = &(frictionContact3D_AlartCurnierFunctionGenerated); */
/* #endif */

/* // HandMade not done */
/* #ifdef AC_HandMade */
/* computeNonsmoothFunction  Function = &(frictionContact3D_AlartCurnierFunctionHandMade); */
/* #endif */

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
  options->iSize = 14;
  options->dSize = 14;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  null_SolverOptions(options);
  for (i = 0; i < 14; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }

  options->iparam[0] = 10;

  options->dparam[0] = 1e-16;
  options->iparam[10] = 0;     /* 0 STD AlartCurnier, 1 JeanMoreau, 2 STD generated, 3 JeanMoreau generated */
  options->iparam[11] = 0;     /* 0 GoldsteinPrice line search, 1 FBLSA */
  options->iparam[12] = 10;   /* max iter line search */

  return 0;
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

/* #ifndef AC_Generated */
#ifndef NDEBUG
    double AWpB[9];
    if (iparam[10] != 3 && iparam[10] != 4)
    {
      double Fg[3] = {0., 0., 0.};
      double Ag[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
      double Bg[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};



      assert(*rho > 0. && *(rho + 1) > 0. && *(rho + 2) > 0.);

/* #ifdef AC_STD */
      if  (iparam[10] == 0 )
      {
        frictionContact3D_AlartCurnierFunctionGenerated(R, velocity, mu, rho, Fg, Ag, Bg);
      }

/* #endif */

/* #ifdef AC_JeanMoreau */
      if  (iparam[10] == 1 )
      {
        frictionContact3D_AlartCurnierJeanMoreauFunctionGenerated(R, velocity, mu, rho, Fg, Ag, Bg);
      }
/* #endif */

      sub3(F, Fg);
      sub3x3(A, Ag);
      sub3x3(B, Bg);

      assert(hypot3(Fg) <= 1e-7);
      assert(hypot9(Ag) <= 1e-7);
      assert(hypot9(Bg) <= 1e-7);
      cpy3x3(A, Ag);
      cpy3x3(B, Bg);
      mm3x3(A, MLocal, AWpB);
      add3x3(B, AWpB);
    }
#endif
/* #endif */


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

/* #ifdef AC_STD */

#ifndef NDEBUG
    if (iparam[10]==0)
    {
      scal3x3(-1., AWpB);
      sub3x3(AWplusB, AWpB);
      assert(hypot9(AWpB) <= 1e-7);
    }
#endif
/* #endif */

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
    Function(R, velocity, mu, rho, F, NULL, NULL);
    dparam[1] = 0.5 * (F[0] * F[0] + F[1] * F[1] + F[2] * F[2]) / (1.0 + sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2])) ; // improve with relative tolerance

    /*      dparam[2] =0.0;
            FrictionContact3D_unitary_compute_and_add_error( R , velocity,mu, &(dparam[2]));*/




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



static int LineSearchGP(FrictionContactProblem* localproblem,
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

  double q0 = 0.5 * cblas_ddot(3 , F , 1 , F , 1);

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

    double q  = 0.5 * cblas_ddot(3 , F , 1 , F , 1);

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
    else   // not(C2)
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
  int LSitermax = iparam[12];


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
