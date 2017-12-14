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


/* /!\ work in progress */

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

#include "SiconosConfig.h"
#include "SparseMatrix_internal.h"
#include "gfc3d_nonsmooth_Newton_AlartCurnier.h"
#include "gfc3d_Solvers.h"
#include "gfc3d_compute_error.h"
#include "AlartCurnierGenerated.h"
#include "fc3d_AlartCurnier_functions.h"
#include "op3x3.h"
#include "SparseBlockMatrix.h"
#include "numerics_verbose.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "sanitizer.h"

#include "SparseMatrix.h"
#include "fc3d_nonsmooth_Newton_AlartCurnier.h"

#include "gfc3d_compute_error.h"
#include "SiconosBlas.h"
#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"
#include "NumericsVector.h"

/* #define DEBUG_MESSAGES 1 */
/* #define DEBUG_STDOUT */
#include <debug.h>

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
  unsigned int m = problem->H->size1;
  unsigned int n = problem->M->size0;
  unsigned int problem_size = n +  2*m ;

  cblas_dscal(problem_size, 0., psi, 1);

  /* -problem->M * globalVelocity + problem->H * reaction + problem->q ==> psi  */
  cblas_dscal(problem_size, 0., psi, 1);

  cblas_dcopy(n, problem->q, 1, psi, 1);
  NM_gemv(1., problem->H, reaction, 1, psi);
  NM_gemv(-1., problem->M, globalVelocity, 1, psi);

  /* -velocity + trans(problem->H) * globalVelocity + problem->b ==> psi + n */
  cblas_daxpy(m, -1., velocity, 1, psi + n, 1);
  cblas_daxpy(m, 1, problem->b, 1, psi + n, 1);
  NM_tgemv(1., problem->H, globalVelocity, 1, psi + n);

  /* compute AC function */
  fc3d_AlartCurnierFunction(m,
                            computeACFun3x3,
                            reaction,
                            velocity, problem->mu, rho,
                            psi+n+m,
                            NULL, NULL);

}


/* init memory for jacobian */
CS_INT initACPsiJacobian(
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
    /* DEBUG_PRINTF("e=%d, M->i[e]=%td, M->p[e]=%td, M->x[e]=%g\n", e, M->i[e], M->p[e], M->x[e]); */
    CHECK_RETURN(cs_zentry(J, M->i[e], M->p[e], - M->x[e]));
  }

  /* H */
  assert(M->n == H->m);
  for(int e = 0; e < H->nz; ++e)
  {
    /* DEBUG_PRINTF("e=%d, H->i[e]=%td, H->p[e] + M->n + A->n=%td, H->x[e]=%g\n", */
    /*              e, H->i[e], H->p[e] + M->n + A->n , H->x[e]); */
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
  CS_INT Astart = J->nz;

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
  CS_INT Astart)
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
    if(fabs(A->x[e]) > DBL_EPSILON)
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
  double alphamin = 1e-3;
  double alphamax = inf;

  double m1 = 0.01, m2 = 0.99;

  unsigned int n = (unsigned)NM_triplet(problem->M)->m;

  unsigned int m = problem->H->size1;

  unsigned int problem_size = n+2*m;

  // Computation of q(t) and q'(t) for t =0

  double q0 = 0.5 * cblas_ddot(problem_size, psi, 1, psi, 1);

  //  tmp <- J * direction
  cblas_dscal(problem_size, 0., tmp, 1);
  cs_gaxpy(J, direction, tmp);

  double dqdt0 = cblas_ddot(problem_size, psi, 1, tmp, 1);
  DEBUG_PRINTF("dqdt0=%e\n",dqdt0);
  DEBUG_PRINTF("q0=%e\n",q0);

  for(unsigned int iter = 0; iter < maxiter_ls; ++iter)
  {

    // tmp <- alpha*direction+solution
    cblas_dcopy(problem_size, solution, 1, tmp, 1);
    cblas_daxpy(problem_size, alpha[0], direction, 1, tmp, 1);

    ACPsi(
      problem,
      computeACFun3x3,
      tmp,  /* v */
      tmp+problem->M->size0+problem->H->size1, /* P */
      tmp+problem->M->size0, /* U */
      rho, psi);

    double q  = 0.5 * cblas_ddot(problem_size, psi, 1, psi, 1);

    assert(q >= 0);

    double slope = (q - q0) / alpha[0];

    int C1 = (slope >= m2 * dqdt0);
    int C2 = (slope <= m1 * dqdt0);

    DEBUG_PRINTF("C1=%i\t C2=%i\n",C1,C2);
    if(C1 && C2)
    {
      if(verbose > 0)
      {
        printf("---- GFC3D - NSN_AC - global line search success. Number of ls iteration = %i  alpha = %.10e, q = %.10e\n",
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
    printf("---- GFC3D - NSN_AC - global line search unsuccessful. Max number of ls iteration reached  = %i  with alpha = %.10e \n",
           maxiter_ls, alpha[0]);
  }

  return -1;
}

void gfc3d_sparseGlobalAlartCurnierInit(
  SolverOptions *SO)
{
}

/* Alart & Curnier solver for sparse global problem */
void gfc3d_nonsmooth_Newton_AlartCurnier(
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


  /* M is square */
  assert(problem->M->size0 == problem->M->size1);
  assert(problem->M->size0 == problem->H->size0);



  unsigned int iter = 0;
  unsigned int itermax = options->iparam[SICONOS_IPARAM_MAX_ITER];
  unsigned int erritermax = options->iparam[SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION];

  if (erritermax == 0)
  {
    /* output a warning here */
    erritermax = 1;
  }

  assert(itermax > 0);


  double tolerance = options->dparam[SICONOS_DPARAM_TOL];
  assert(tolerance > 0);

  if (verbose > 0)
    printf("---- GFC3D - _nonsmooth_Newton_AlartCurnier - Start with tolerance = %g\n", tolerance);



  DEBUG_PRINTF("norm of M = %e\n", NM_norm(problem->M));
  DEBUG_PRINTF("norm of H = %e\n", NM_norm(problem->H));

  /* sparse triplet storage */
  NM_triplet(problem->M);
  NM_triplet(problem->H);


  unsigned int n = (unsigned)NM_triplet(problem->M)->m;

  unsigned int m = problem->H->size1;

  unsigned int problem_size = n+2*m;




  assert((int)m == problem->numberOfContacts * problem->dimension);

  assert((int)n == problem->H->size0); /* size(velocity) ==
                                                   * Htrans*globalVelocity */


  AlartCurnierFun3x3Ptr computeACFun3x3 = NULL;

  switch (options->iparam[SICONOS_FRICTION_3D_NSN_FORMULATION])
  {
  case SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD:
  {
    computeACFun3x3 = &computeAlartCurnierSTD;
    break;
  }
  case SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD:
  {
    computeACFun3x3 = &computeAlartCurnierJeanMoreau;
    break;
  };
  case SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_GENERATED:
  {
    computeACFun3x3 = &fc3d_AlartCurnierFunctionGenerated;
    break;
  }
  case SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_GENERATED:
  {
    computeACFun3x3 = &fc3d_AlartCurnierJeanMoreauFunctionGenerated;
    break;
  }
  }

  if(options->iparam[SICONOS_FRICTION_3D_NSN_MEMORY_ALLOCATION] == 0)
  {
    /* allocate memory */
    assert(options->dWork == NULL);
    assert(options->iWork == NULL);
    options->dWork = (double *) calloc(
                       (m + /* F */
                        3 * m + /* A */
                        3 * m + /* B */
                        m + /* rho */
                        problem_size + /* psi */
                        problem_size + /* rhs */
                        problem_size + /* tmp2 */
                        problem_size + /* tmp3 */
                        problem_size   /* solution */) ,
                       sizeof(double));

    /* XXX big hack here */
    options->iWork = (int *) malloc(
                       (3 * m + /* iA */
                        3 * m + /* iB */
                        3 * m + /* pA */
                        3 * m)  /* pB */
                       * sizeof(CS_INT));

    options->iparam[SICONOS_FRICTION_3D_NSN_MEMORY_ALLOCATION] = 1;

  }

  assert(options->dWork != NULL);
  assert(options->iWork != NULL);

  double *F = options->dWork;
  double *A = F +  m;
  double *B = A +  3 * m;
  double *rho = B + 3 * m;

  double * psi = rho + m;
  double * rhs = psi + problem_size;
  double * tmp2 = rhs + problem_size;
  double * tmp3 = tmp2 + problem_size;
  double * solution = tmp3 + problem_size;

  /* XXX big hack --xhub*/
  CS_INT * iA = (CS_INT *)options->iWork;
  CS_INT * iB = iA + 3 * m;
  CS_INT * pA = iB + 3 * m;
  CS_INT * pB = pA + 3 * m;

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

  fc3d_AlartCurnierFunction(
    m,
    computeACFun3x3,
    reaction, velocity,
    problem->mu, rho,
    F, A, B);

  CS_INT Astart = initACPsiJacobian(NM_triplet(problem->M),
                                 NM_triplet(problem->H),
                                 &A_, &B_, J);

  assert(Astart > 0);

  assert(A_.m == A_.n);
  assert(B_.m == B_.n);

  assert(A_.m == problem->H->size1);

  // compute rho here
  for(unsigned int i = 0; i < m; ++i) rho[i] = options->dparam[SICONOS_FRICTION_3D_NSN_RHO];

  DEBUG_EXPR(for(unsigned int i = 0; i < m; ++i) printf("rho[%i]=%e\t",i,rho[i]); printf("\n"););

  // direction
  for(unsigned int i = 0; i < problem_size; ++i) rhs[i] = 0.;



  // quick hack to make things work
  // need to use the functions from NumericsMatrix --xhub


  NumericsSparseMatrix* SM = newNumericsSparseMatrix();
  SM->triplet = J;
  NumericsMatrix *AA = NM_create_from_data(NM_SPARSE,  (int)J->m, (int)J->n, SM);

  info[0] = 1;

  /* update local velocity from global velocity */
  /* an assertion ? */
  /* double norm_globalVelocity = cblas_dnrm2(n,globalVelocity,1); */
  /* double norm_reaction = cblas_dnrm2(m,reaction,1); */
  /* DEBUG_PRINTF("norm of globalVelocity = %e\n", norm_globalVelocity); */
  /* DEBUG_PRINTF("norm of reaction = %e\n", norm_reaction); */
  /* if ((fabs(norm_globalVelocity) < DBL_EPSILON) && (fabs(norm_reaction) < DBL_EPSILON) ) */
  /* { */
  /*   cblas_dcopy(n, problem->q, 1, globalVelocity, 1); */
  /*   int info_solver = NM_gesv(problem->M, globalVelocity, true); */
  /* } */

  cblas_dcopy(m, problem->b, 1, velocity, 1);
  NM_tgemv(1., problem->H, globalVelocity, 1, velocity);

  double norm_velocity = cblas_dnrm2(m,velocity,1);
  if (fabs(norm_velocity) < DBL_EPSILON)
  {
    for (unsigned int k =0; k < m; k++) velocity[k]=0.0;
  }

  DEBUG_PRINTF("norm of velocity = %e\n", cblas_dnrm2(m,velocity,1));
  double linear_solver_residual=0.0;

  DEBUG_EXPR(NV_display(globalVelocity,n););
  DEBUG_EXPR(NV_display(velocity,m););
  DEBUG_EXPR(NV_display(reaction,m););


  /* store the current solution */

  double* globalVelocity_k = solution;
  double* velocity_k = &(solution[n]);
  double* reaction_k = &(solution[n+m]);

  memcpy(globalVelocity_k, globalVelocity,  n * sizeof(double));
  memcpy(velocity_k, velocity,  m * sizeof(double));
  memcpy(reaction_k, reaction, m * sizeof(double));

  while(iter++ < itermax)
  {
    /* /\* store the current solution *\/ */
    /* for(unsigned int i = 0; i < n; ++i) */
    /* { */
    /*   solution[i] = globalVelocity[i]; */
    /* } */
    /* for(unsigned int i = 0; i < m; ++i) */
    /* { */
    /*   solution[i+n] = velocity[i]; */
    /*   solution[i+n+m] = reaction[i]; */
    /* } */

    /* compute psi */
    ACPsi(problem, computeACFun3x3, globalVelocity_k, reaction_k, velocity_k, rho, psi);

    /* compute A & B */
    fc3d_AlartCurnierFunction(m,
                              computeACFun3x3,
                              reaction_k, velocity_k,
                              problem->mu, rho,
                              F, A, B);
    /* update J */
    updateACPsiJacobian(NM_triplet(problem->M),
                        NM_triplet(problem->H),
                        &A_, &B_, J, Astart);

    /* rhs = -psi */
    cblas_dcopy(problem_size, psi, 1, rhs, 1);
    cblas_dscal(problem_size, -1., rhs, 1);

    /* get compress column storage for linear ops */
    CSparseMatrix* Jcsc = cs_compress(J);

    /* Solve: J X = -psi */

    /* Solve: AWpB X = -F */
    int info_solver = NM_gesv(AA, rhs, true);
    if (info_solver > 0)
    {
      fprintf(stderr, "---- GFC3D - NSN_AC - solver failed info = %d\n", info_solver);
      break;
      info[0] = 2;
      CHECK_RETURN(!cs_check_triplet(NM_triplet(AA)));
    }

    /* Check the quality of the solution */
    if (verbose > 0)
    {
      cblas_dcopy_msan(problem_size, psi, 1, tmp3, 1);
      NM_gemv(1., AA, rhs, 1., tmp3);
      linear_solver_residual = cblas_dnrm2(problem_size, tmp3, 1);

      DEBUG_PRINTF("fc3d esolve: linear equation residual = %g\n",
              cblas_dnrm2(problem_size, tmp3, 1));
      /* for the component wise scaled residual: cf mumps &
       * http://www.netlib.org/lapack/lug/node81.html */
    }

    DEBUG_EXPR(NV_display(rhs,problem_size););


    /* line search */

    double alpha = 1.0;
    int info_ls = 0;

    switch (options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH])
    {
    case SICONOS_FRICTION_3D_NSN_LINESEARCH_NO:
      /* without line search */
      info_ls = 1;
      break;

    case SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE:
      /* Goldstein Price */
      info_ls = _globalLineSearchSparseGP(problem,
                                          computeACFun3x3,
                                          solution,
                                          rhs,
                                          problem->mu, rho, F, psi, Jcsc,
                                          tmp2, &alpha, options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH_MAXITER]);
      break;
    /* case SICONOS_FRICTION_3D_NSN_LINESEARCH_ARMIJO: */
    /*   /\* FBLSA *\/ */
    /*   info_ls = frictionContactFBLSA(equation, reaction, velocity, problem->mu, rho, F, Ax, Bx, */
    /*                                  problem->M, problem->q, AWpB, tmp1, tmp2, &alpha, options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH_MAXITER]); */
    /*   break; */
    default:
    {
      numerics_error("gfc3d_nonsmooth_Newton_AlartCurnier",
                     "Unknown line search option.\n");
    }
    }

    cs_spfree(Jcsc);
    /* if(!info_ls) */
    /* { */
    /*   cblas_daxpy(problem_size, alpha, rhs, 1, solution, 1); */
    /* } */
    /* else */
    /* { */
    /*   cblas_daxpy(problem_size, 1, rhs, 1., solution, 1); */
    /* } */
    cblas_daxpy(problem_size, alpha, rhs, 1, solution, 1);


    DEBUG_EXPR(NV_display(globalVelocity_k,n););
    DEBUG_EXPR(NV_display(velocity_k,m););
    DEBUG_EXPR(NV_display(reaction_k,m););


    options->dparam[SICONOS_DPARAM_RESIDU] = INFINITY;

    if(!(iter % erritermax))
    {
      double norm_q = cblas_dnrm2(problem->M->size0 , problem->q , 1);
      gfc3d_compute_error(problem,
                          reaction_k, velocity_k, globalVelocity_k,
                          tolerance,
                          norm_q,
                          &(options->dparam[SICONOS_DPARAM_RESIDU]));
    }

    if(verbose > 0)
      printf("---- GFC3D - NSN_AC - iteration %d, residual = %g, linear solver residual = %g, tolerance = %g \n", iter, options->dparam[1],linear_solver_residual, tolerance);

    if(options->dparam[SICONOS_DPARAM_RESIDU] < tolerance)
    {
      info[0] = 0;
      break;
    }


  }

  memcpy(globalVelocity, solution,  n * sizeof(double));
  memcpy(velocity, &solution[n],  m * sizeof(double));
  memcpy(reaction, &solution[n+m], m * sizeof(double));

  if(verbose > 0)
  {
    if(!info[0])
      printf("---- GFC3D - NSN_AC - convergence after %d iterations, residual = %g\n",
             iter, options->dparam[1]);
    else
    {
      printf("---- GFC3D - NSN_AC - no convergence after %d iterations, residual = %g\n",
             iter, options->dparam[1]);
    }
  }

  options->iparam[SICONOS_IPARAM_ITER_DONE] = iter;

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

  NM_free(AA);
  free(AA);
}

int gfc3d_nonsmooth_Newton_AlartCurnier_setDefaultSolverOptions(
  SolverOptions* options)
{
  if(verbose > 0)
  {
    printf("Set the default solver options for the GLOBAL_AC Solver\n");
  }

  options->solverId = SICONOS_GLOBAL_FRICTION_3D_NSN_AC;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 20;
  options->dSize = 20;
  options->iparam = (int *) calloc(options->iSize, sizeof(int));
  options->dparam = (double *) calloc(options->dSize,  sizeof(double));
  options->dWork = NULL;
  solver_options_nullify(options);


  options->iparam[SICONOS_IPARAM_MAX_ITER] = 200;    /* input :  itermax */
  options->iparam[SICONOS_IPARAM_ITER_DONE] = 1;      /* output : #iter */

  options->iparam[SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION] = 1;      /* erritermax */

  options->iparam[SICONOS_FRICTION_3D_NSN_MPI_COM] = -1;     /* mpi com fortran */

  options->iparam[SICONOS_FRICTION_3D_NSN_MEMORY_ALLOCATION] = 0;      /* > 0 memory is allocated */

  options->iparam[SICONOS_FRICTION_3D_NSN_FORMULATION] = SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD;     /* 0 STD AlartCurnier, 1 JeanMoreau, 2 STD generated, 3 JeanMoreau generated */
  options->iparam[SICONOS_FRICTION_3D_NSN_FORMULATION] = SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD;     /* 0 STD AlartCurnier, 1 JeanMoreau, 2 STD generated, 3 JeanMoreau generated */
  options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH] = SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE;
  options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH_MAXITER]=10;



  options->dparam[SICONOS_DPARAM_TOL] = 1e-10;
  options->dparam[SICONOS_FRICTION_3D_NSN_RHO] = 1e-3;      /* default rho */


  options->internalSolvers = NULL;

  if(verbose > 0)
  {
    solver_options_print(options);
  }


  return 0;
}
