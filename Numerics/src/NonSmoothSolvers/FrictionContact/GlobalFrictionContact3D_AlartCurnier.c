/* Siconos-Numerics, Copyright INRIA 2005-2012.
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


/* /!\ work in progress */

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

#include "NumericsConfig.h"

#include "GlobalFrictionContact3D_AlartCurnier.h"
#include "GlobalFrictionContact3D_Solvers.h"
#include "GlobalFrictionContact3D_compute_error.h"
#include "AlartCurnierGenerated.h"
#include "FrictionContact3D_AlartCurnier.h"
#include "op3x3.h"
#include "SparseBlockMatrix.h"
#include "misc.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "SparseMatrix.h"
#include "FrictionContact3D_localAlartCurnier.h"

#include "GlobalFrictionContact3D_compute_error.h"
#include "SiconosBlas.h"

//#define DEBUG_MESSAGES 1
#include <debug.h>


#ifdef WITH_MUMPS
#include <mpi.h>
#include <dmumps_c.h>

#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654
#define ICNTL(I) icntl[(I)-1]
#define CNTL(I) cntl[(I)-1]
#endif

/* size of whole problem */
unsigned int sizeOfPsi(
  CSparseMatrix* M,
  CSparseMatrix* H)
{
  assert(M->n + H->n + H->m >= 0);
  return (unsigned)(M->n + H->n + H->m);
}

/* compute psi function */
void ACPsi(
  GlobalFrictionContactProblem* problem,
  AlartCurnierFun3x3Ptr computeACFun3x3,
  double *globalVelocity,
  double *reaction,
  double *velocity,
  double *rho,
  double *psi)
{

  assert(problem->H->size1 == problem->dimension * problem->numberOfContacts);

  unsigned int localProblemSize = problem->H->size1;

  unsigned int ACProblemSize = sizeOfPsi(NM_triplet(problem->M),
                                         NM_triplet(problem->H));

  unsigned int globalProblemSize = problem->M->size0;


  /* psi <-
       compute -problem->M * globalVelocity + problem->H * reaction + problem->q
     ... */
  cblas_dscal(ACProblemSize, 0., psi, 1);
  cblas_dcopy(globalProblemSize, problem->q, 1, psi, 1);
  NM_gemv(1., problem->H, reaction, 1, psi);
  NM_gemv(-1., problem->M, globalVelocity, 1, psi);


  /* psi + globalProblemSize <-
     compute -velocity + trans(problem->H) * globalVelocity + problem->b
   ... */
  cblas_daxpy(localProblemSize, -1., velocity, 1, psi + globalProblemSize, 1);
  cblas_daxpy(localProblemSize, 1, problem->b, 1, psi + globalProblemSize, 1);
  NM_tgemv(1., problem->H, globalVelocity, 1, psi + globalProblemSize);



  /* compute AC function */
  frictionContact3D_AlartCurnierFunction(localProblemSize,
                                         computeACFun3x3,
                                         reaction,
                                         velocity, problem->mu, rho,
                                         psi+globalProblemSize+
                                         problem->H->size1,
                                         NULL, NULL);

}





/* init memory for jacobian */
csi initACPsiJacobian(
  CSparseMatrix* M,
  CSparseMatrix* H,
  CSparseMatrix *A,
  CSparseMatrix *B,
  CSparseMatrix *J)
{
  /* only triplet matrix */
  assert(M->nz>=0);
  assert(H->nz>=0);
  assert(A->nz>=0);
  assert(B->nz>=0);

  /* M square */
  assert(M->m == M->n);

  /* A & B squares */
  assert(A->m == A->n);
  assert(B->m == B->n);

  assert(A->nz == B->nz);

  /* - M */
  for(int e = 0; e < M->nz; ++e)
  {
    DEBUG_PRINTF("e=%d, M->i[e]=%d, M->p[e]=%d, M->x[e]=%g\n", e, M->i[e], M->p[e], M->x[e]);
    CHECK_RETURN(cs_zentry(J, M->i[e], M->p[e], - M->x[e]));
  }

  /* H */
  assert(M->n == H->m);
  for(int e = 0; e < H->nz; ++e)
  {
    DEBUG_PRINTF("e=%d, H->i[e]=%d, H->p[e] + M->n + A->n=%d, H->x[e]=%g\n", 
                 e, H->i[e], H->p[e] + M->n + A->n , H->x[e]);
    CHECK_RETURN(cs_zentry(J, H->i[e], H->p[e] + M->n + A->n, H->x[e]));
  }

  /* Ht */
  for(int e = 0; e < H->nz; ++e)
  {
    CHECK_RETURN(cs_zentry(J, H->p[e] + M->m, H->i[e], H->x[e]));
  }

  /* -I */
  for(int e = 0; e < A->m; ++e)
  {
    CHECK_RETURN(cs_zentry(J, e + M->m, e + M->n, -1.));
  }

  /* keep A start indice for update */
  csi Astart = J->nz;

  /* A */
  for(int e = 0; e < A->nz; ++e)
  {
    CHECK_RETURN(cs_zentry(J, A->i[e] + M->m + H->n, A->p[e] + M->n, A->x[e]));
  }

  /* B */
  for(int e = 0; e < B->nz; ++e)
  {
    CHECK_RETURN(cs_zentry(J, B->i[e] + M->m + H->n, B->p[e] + M->n + A->n, B->x[e]));
  }

  return Astart;
}

/* update J with new A and B */
void updateACPsiJacobian(
  CSparseMatrix* M,
  CSparseMatrix* H,
  CSparseMatrix *A,
  CSparseMatrix *B,
  CSparseMatrix *J,
  csi Astart)
{
  /* only triplet matrix */
  assert(M->nz>=0);
  assert(H->nz>=0);
  assert(A->nz>=0);
  assert(B->nz>=0);

  /* M square */
  assert(M->m == M->n);

  /* A & B squares */
  assert(A->m == A->n);
  assert(B->m == B->n);


  assert(J->n == M->n + A->m + B->m);
  assert(J->m == M->m + A->n + B->n);

  assert(J->p);
  assert(J->i);
  assert(J->x);

  if(((Astart + A->nz + B->nz) > J->nzmax))
  {
    CHECK_RETURN(cs_sprealloc(J, Astart + A->nz + B->nz));
  }

  /* A */
  J->nz = Astart;

  for(int e = 0; e < A->nz; ++e)
  {
    if(fabs(A->x[e]) >= DBL_EPSILON)
    {

      J->i[J->nz] = A->i[e] + M->m + H->n;
      J->p[J->nz] = A->p[e] + M->n;
      J->x[J->nz] = A->x[e];
      J->nz++;

      assert(J->nz <= J->nzmax);

    }
  }

  /* B */
  for(int e = 0; e < B->nz; ++e)
  {
    if(fabs(B->x[e]) > DBL_EPSILON)
    {
      J->i[J->nz] = B->i[e] + M->m + H->n;
      J->p[J->nz] = B->p[e] + M->n + A->n;
      J->x[J->nz] = B->x[e];
      J->nz++;

      assert(J->nz <= J->nzmax);

    }
  }
}

/* 3D lists blocks => sparse matrix */
void init3x3DiagBlocks(int nc, double* P, CSparseMatrix* R)
{

  R->m = 3 * nc;
  R->n = 3 * nc;
  R->x = P;
  R->nz = 9*nc;
  R->nzmax = R->nz;

  for(int ib = 0; ib < nc; ++ib)
  {
    for(int j = 0; j < 3; ++j)
    {
      for(int i = 0; i < 3; ++i)
      {
        R->i[9*ib+i+3*j] = ib * 3 + i;
        R->p[9*ib+i+3*j] = ib * 3 + j;
      }
    }
  }
}

/*
*/
int _globalLineSearchSparseGP(
  GlobalFrictionContactProblem *problem,
  AlartCurnierFun3x3Ptr computeACFun3x3,
  double *solution,
  double *direction,
  double *globalVelocity,
  double *reaction,
  double *velocity,
  double *mu,
  double *rho,
  double *F,
  double *psi,
  CSparseMatrix *J,
  double *tmp,
  double alpha[1],
  unsigned int maxiter_ls)
{
  double inf = 1e10;
  double alphamin = 0.0;
  double alphamax = inf;

  double m1 = 0.01, m2 = 0.99;

  unsigned int ACProblemSize = sizeOfPsi(NM_triplet(problem->M),
                                         NM_triplet(problem->H));

  // Computation of q(t) and q'(t) for t =0

  double q0 = 0.5 * cblas_ddot(ACProblemSize, psi, 1, psi, 1);

  //  tmp <- J * direction
  cblas_dscal(ACProblemSize, 0., tmp, 1);
  cs_gaxpy(J, direction, tmp);

  double dqdt0 = cblas_ddot(ACProblemSize, psi, 1, tmp, 1);

  for(unsigned int iter = 0; iter < maxiter_ls; ++iter)
  {

    // tmp <- alpha*direction+solution
    cblas_dcopy(ACProblemSize, solution, 1, tmp, 1);
    cblas_daxpy(ACProblemSize, alpha[0], direction, 1, tmp, 1);

    ACPsi(
      problem,
      computeACFun3x3,
      tmp,  /* v */
      tmp+problem->M->size0+problem->H->size1, /* P */
      tmp+problem->M->size0, /* U */
      rho, psi);

    double q  = 0.5 * cblas_ddot(ACProblemSize, psi, 1, psi, 1);

    assert(q >= 0);

    double slope = (q - q0) / alpha[0];

    int C1 = (slope >= m2 * dqdt0);
    int C2 = (slope <= m1 * dqdt0);

    if(C1 && C2)
    {
      if(verbose > 0)
      {
        printf("global line search success.\n");
        printf("Number of iteration = %i  alpha = %.10e, q = %.10e\n",
               iter,
               alpha[0], q);
      }

      return 0;

    }
    else if(!C1)
    {
      alphamin = alpha[0];
    }
    else
    {
      // not(C2)
      alphamax = alpha[0];
    }

    if(alpha[0] < inf)
    {
      alpha[0] = 0.5 * (alphamin + alphamax);
    }
    else
    {
      alpha[0] = alphamin;
    }

  }
  if(verbose > 0)
  {
    printf("global line search unsuccessfull.\n");
    printf("max number of iteration reached  = %i  with alpha = %.10e \n",
           maxiter_ls, alpha[0]);
  }

  return -1;
}

int globalFrictionContact3D_AlartCurnier_setDefaultSolverOptions(
  SolverOptions* options)
{
  if(verbose > 0)
  {
    printf("Set the default solver options for the GLOBALAC Solver\n");
  }

  options->solverId = SICONOS_FRICTION_3D_GLOBAL_AC;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 14;
  options->dSize = 11;
  options->iparam = (int *) calloc(options->iSize, sizeof(int));
  options->dparam = (double *) calloc(options->dSize,  sizeof(double));
  options->dWork = NULL;
  null_SolverOptions(options);
  options->iparam[0] = 200;    /* input :  itermax */
  options->iparam[1] = 1;      /* output : #iter */
  options->iparam[2] = 0;      /* unused */
  options->iparam[3] = 100000; /* nzmax*/
  options->iparam[4] = 0;      /* unused */
  options->iparam[5] = 1;      /* unused */
  options->iparam[7] = 1;      /* erritermax */
  options->iparam[8] = -1;     /* mpi com fortran */

  options->iparam[9] = 0;      /* > 0 memory is allocated */

  options->iparam[10] = 0;     /* 0 STD AlartCurnier, 1 JeanMoreau, 2 STD generated, 3 JeanMoreau generated */

  options->dparam[0] = 1e-10;


  options->internalSolvers = NULL;

  return 0;
}

void globalFrictionContact3D_sparseGlobalAlartCurnierInit(
  SolverOptions *SO)
{
#ifdef WITH_MUMPS
  DMUMPS_STRUC_C* mumps_id = (DMUMPS_STRUC_C*)malloc(sizeof(DMUMPS_STRUC_C));

  // SO with void pointers ?
  SO->dparam[7] = (double) (long long) mumps_id;

  // Initialize a MUMPS instance. Use MPI_COMM_WORLD.
  mumps_id->job = JOB_INIT;
  mumps_id->par = 1;
  mumps_id->sym = 0;
  mumps_id->comm_fortran = USE_COMM_WORLD;
  dmumps_c(mumps_id);

  if(verbose > 1)
  {
    mumps_id->ICNTL(4) = 0;
    mumps_id->ICNTL(10) = 1;
    mumps_id->ICNTL(11) = 1;
  }
  else
  {
    mumps_id->ICNTL(1) = -1;
    mumps_id->ICNTL(2) = -1;
    mumps_id->ICNTL(3) = -1;
  }

  mumps_id->ICNTL(24) = 1; /* Null pivot row detection see also
                              CNTL(3) & CNTL(5) */

#endif
}

/* Alart & Curnier solver for sparse global problem */
void globalFrictionContact3D_AlartCurnier(
  GlobalFrictionContactProblem* problem,
  double *reaction,
  double *velocity,
  double *globalVelocity,
  int *info,
  SolverOptions* options)
{

  assert(problem);
  assert(reaction);
  assert(velocity);
  assert(info);
  assert(options);

  assert(problem->dimension == 3);

  assert(options->iparam);
  assert(options->dparam);

  assert(problem->q);
  assert(problem->mu);
  assert(problem->M);
  assert(problem->H);

  assert(!problem->M->matrix0);
//  assert(problem->M->matrix1);

  assert(!options->iparam[4]); // only host

  /* M is square */
  assert(problem->M->size0 == problem->M->size1);

  assert(problem->M->size0 == problem->H->size0);

  unsigned int iter = 0;
  unsigned int itermax = options->iparam[0];
  unsigned int erritermax = options->iparam[7];

  if (erritermax == 0)
  {
    /* output a warning here */
    erritermax = 1;
  }

#ifdef WITH_MUMPS
  DMUMPS_STRUC_C* mumps_id = (DMUMPS_STRUC_C*)(long) options->dparam[7];
  int mumps_com = options->iparam[8];

  if(mumps_com<0)
  {
    /* we suppose mpi init has not been done */
    int ierr MAYBE_UNUSED;
    int myid;
    /* c99 mandates that argv[argc] = NULL --xhub */
    /* Also where are we calling MPI_Finalize ?  */
    int argc = 0;
    char * argv0 = NULL;
    char ** argv = &argv0;
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    globalFrictionContact3D_sparseGlobalAlartCurnierInit(options);
    mumps_id = (DMUMPS_STRUC_C*)(long) options->dparam[7];
  }
  else /* an mpi communicator is passed through iparam[8] so mpi init
        * has been done */
  {
    globalFrictionContact3D_sparseGlobalAlartCurnierInit(options);
    mumps_id = (DMUMPS_STRUC_C*)(long) options->dparam[7];
    mumps_id->comm_fortran = mumps_com;
  }

  assert(mumps_id);
#endif

  assert(itermax > 0);
  assert(options->iparam[3] > 0);

  double tolerance = options->dparam[0];
  assert(tolerance > 0);

  /* sparse triplet storage */
  NM_triplet(problem->M);
  NM_triplet(problem->H);

  unsigned int ACProblemSize = sizeOfPsi(NM_triplet(problem->M),
                                         NM_triplet(problem->H));

  unsigned int globalProblemSize = (unsigned)NM_triplet(problem->M)->m;

  unsigned int localProblemSize = problem->H->size1;

  assert(localProblemSize == problem->numberOfContacts * problem->dimension);

  assert(globalProblemSize == problem->H->size0); /* size(velocity) ==
                                                   * Htrans*globalVelocity */


  AlartCurnierFun3x3Ptr computeACFun3x3 = NULL;

  switch (options->iparam[10])
  {
  case 0:
  {
    computeACFun3x3 = &computeAlartCurnierSTD;
    break;
  }
  case 1:
  {
    computeACFun3x3 = &computeAlartCurnierJeanMoreau;
    break;
  };
  case 2:
  {
    computeACFun3x3 = &frictionContact3D_AlartCurnierFunctionGenerated;
    break;
  }
  case 3:
  {
    computeACFun3x3 = &frictionContact3D_AlartCurnierJeanMoreauFunctionGenerated;
    break;
  }
  }

  if(options->iparam[9] == 0)
  {
    /* allocate memory */
    assert(options->dWork == NULL);
    assert(options->iWork == NULL);
    options->dWork = (double *) malloc(
                       (localProblemSize + /* F */
                        3 * localProblemSize + /* A */
                        3 * localProblemSize + /* B */
                        localProblemSize + /* rho */
                        ACProblemSize + /* psi */
                        ACProblemSize + /* rhs */
                        ACProblemSize + /* tmp2 */
                        ACProblemSize   /* solution */) *  sizeof(double));

    options->iWork = (int *) malloc(
                       (3 * localProblemSize + /* iA */
                        3 * localProblemSize + /* iB */
                        3 * localProblemSize + /* pA */
                        3 * localProblemSize)  /* pB */
                       * sizeof(int));

    options->iparam[9] = 1;

  }

  assert(options->dWork != NULL);
  assert(options->iWork != NULL);

  double *F = options->dWork;
  double *A = F +   localProblemSize;
  double *B = A +   3 * localProblemSize;
  double *rho = B + 3 * localProblemSize;

  double * psi = rho + localProblemSize;
  double * rhs = psi + ACProblemSize;
  double * tmp2 = rhs + ACProblemSize;
  double * solution = tmp2 + ACProblemSize;

  /* XXX big hack --xhub*/
  csi * iA = (csi *)options->iWork;
  csi * iB = iA + 3 * localProblemSize;
  csi * pA = iB + 3 * localProblemSize;
  csi * pB = pA + 3 * localProblemSize;

  CSparseMatrix A_;
  CSparseMatrix B_;
  CSparseMatrix *J;

  A_.p = pA;
  B_.p = pB;
  A_.i = iA;
  B_.i = iB;

  init3x3DiagBlocks(problem->numberOfContacts, A, &A_);
  init3x3DiagBlocks(problem->numberOfContacts, B, &B_);

  J = cs_spalloc(NM_triplet(problem->M)->n + A_.m + B_.m,
                 NM_triplet(problem->M)->n + A_.m + B_.m,
                 NM_triplet(problem->M)->nzmax + 2*NM_triplet(problem->H)->nzmax +
                 2*A_.n + A_.nzmax + B_.nzmax, 1, 1);

  assert(A_.n == problem->H->size1);
  assert(A_.nz == problem->numberOfContacts * 9);
  assert(B_.n == problem->H->size1);
  assert(B_.nz == problem->numberOfContacts * 9);

  frictionContact3D_AlartCurnierFunction(
    localProblemSize,
    computeACFun3x3,
    reaction, velocity,
    problem->mu, rho,
    F, A, B);

  csi Astart = initACPsiJacobian(NM_triplet(problem->M),
                                 NM_triplet(problem->H),
                                 &A_, &B_, J);

  assert(Astart > 0);

  assert(A_.m == A_.n);
  assert(B_.m == B_.n);

  assert(A_.m == problem->H->size1);

  // compute rho here
  for(unsigned int i = 0; i < localProblemSize; ++i) rho[i] = 1.;

  // direction
  for(unsigned int i = 0; i < ACProblemSize; ++i) rhs[i] = 0.;


#ifdef WITH_MUMPS
  mumps_id->n = (int) J->m;
  mumps_id->nz = (int) J->nz;
  mumps_id->irn = (int) J->i;
  mumps_id->jcn = (int) J->p;
  mumps_id->a = J->x;
  mumps_id->rhs = rhs;

  mumps_id->job = 6;
#endif

  info[0] = 1;

  /* update local velocity from global velocity */
  /* an assertion ? */
  cblas_dcopy(localProblemSize, problem->b, 1, velocity, 1);
  NM_tgemv(1., problem->H, globalVelocity, 1, velocity);

  while(iter++ < itermax)
  {

    /* compute psi */
    ACPsi(problem, computeACFun3x3, globalVelocity, reaction, velocity, rho, psi);

    /* compute A & B */
    frictionContact3D_AlartCurnierFunction(localProblemSize,
                                           computeACFun3x3,
                                           reaction, velocity,
                                           problem->mu, rho,
                                           F, A, B);
    /* update J */
    updateACPsiJacobian(NM_triplet(problem->M),
                        NM_triplet(problem->H),
                        &A_, &B_, J, Astart);

    /* rhs = -psi */
    cblas_dcopy(ACProblemSize, psi, 1, rhs, 1);
    cblas_dscal(ACProblemSize, -1., rhs, 1);

    /* get compress column storage for linear ops */
    CSparseMatrix* Jcsc = cs_compress(J);

    /* Solve: J X = -psi */

#ifdef WITH_MUMPS
    mumps_id->n = (int) J->m;
    mumps_id->nz = (int) J->nz;

    /* set fortran indices in place so we do not have to bother with
     * yet another heap allocation */
    for(int e = 0 ; e < J->nz ; ++e)
    {
      J->i[e]++;
      J->p[e]++;
    }

    dmumps_c(mumps_id);

    /* come back to C indices */
    for(int e = 0 ; e < J->nz ; ++e)
    {
      J->i[e]--;
      J->p[e]--;
    }

    assert(mumps_id->info[0] >= 0);

    if(mumps_id->info[0] > 0)
      /*if (verbose>0)*/
      printf("GLOBALAC: MUMPS warning : info(1)=%d, info(2)=%d\n",
             mumps_id->info[0], mumps_id->info[1]);


    if(verbose > 1)
    {

      printf("mumps : condition number %g\n", mumps_id->rinfog[9]);
      printf("mumps : component wise scaled residual %g\n",
             mumps_id->rinfog[6]);
      printf("mumps : \n");

    }
#else
    /* use csparse LU factorization */
    CHECK_RETURN(cs_lusol(1, Jcsc, rhs, DBL_EPSILON));

#endif

    /* line search */
    double alpha = 1;

    /* set current solution */
    for(unsigned int i = 0; i < globalProblemSize; ++i)
    {
      solution[i] = globalVelocity[i];
    }
    for(unsigned int i = 0; i < localProblemSize; ++i)
    {
      solution[i+globalProblemSize] = velocity[i];
      solution[i+globalProblemSize+localProblemSize] = reaction[i];
    }

    int info_ls = _globalLineSearchSparseGP(problem,
                                            computeACFun3x3,
                                            solution,
                                            rhs,
                                            globalVelocity,
                                            reaction, velocity,
                                            problem->mu, rho, F, psi, Jcsc,
                                            tmp2, &alpha, 100);


    cs_spfree(Jcsc);

    if(!info_ls)
    {
      cblas_daxpy(ACProblemSize, alpha, rhs, 1, solution, 1);
    }
    else
    {
      cblas_daxpy(ACProblemSize, 1, rhs, 1., solution, 1);
    }

    for(unsigned int e = 0 ; e < globalProblemSize; ++e)
    {
      globalVelocity[e] = solution[e];
    }

    for(unsigned int e = 0 ; e < localProblemSize; ++e)
    {
      velocity[e] = solution[e+globalProblemSize];
    }

    for(unsigned int e = 0; e < localProblemSize; ++e)
    {
      reaction[e] = solution[e+globalProblemSize+localProblemSize];
    }

    options->dparam[1] = INFINITY;

    if(!(iter % erritermax))
    {

      GlobalFrictionContact3D_compute_error(problem,
                                            reaction, velocity, globalVelocity,
                                            tolerance,
                                            &(options->dparam[1]));

    }

    if(verbose > 0)
      printf("GLOBALAC: iteration %d : error=%g\n", iter, options->dparam[1]);

    if(options->dparam[1] < tolerance)
    {
      info[0] = 0;
      break;
    }


  }

  if(verbose > 0)
  {
    if(!info[0])
      printf("GLOBALAC: convergence after %d iterations, error : %g\n",
             iter, options->dparam[1]);
    else
    {
      printf("GLOBALAC: no convergence after %d iterations, error : %g\n",
             iter, options->dparam[1]);
    }
  }

  options->iparam[1] = iter;

#ifdef DUMP_PROBLEM
  if(info[0])
  {
    static int file_counter = 0;
    char filename[64];
    printf("GLOBALAC: dumping problem\n");
    sprintf(filename, "GLOBALAC_failure%d.dat", file_counter++);
    FILE* file = fopen(filename, "w");
    frictionContact_printInFile(problem, file);
    fclose(file);
  }
#endif

  {
    cs_spfree(J);
  }
}



