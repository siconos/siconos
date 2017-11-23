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

/* Factorisation with Newton_methods.c is needed */

#include "fc3d_nonsmooth_Newton_solvers.h"

#include "NumericsMatrix_private.h"

/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "debug.h"
#include "op3x3.h"
#include "SparseBlockMatrix.h"
#include "SparseMatrix.h"
#include "fc3d_Solvers.h"
#include "FrictionContactProblem.h"
#include "fc3d_compute_error.h"
#include "AlartCurnierGenerated.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "Friction_cst.h"
#include "SiconosLapack.h"
#include "NumericsSparseMatrix.h"
#include "sanitizer.h"
#include "numerics_verbose.h"

static void NM_dense_to_sparse_diag_t(double* A, NumericsMatrix* B, size_t block_row_size, size_t block_col_size)
{
  /* TODO  CSC, CSR version*/
  assert(A);
  if(!B->matrix2->triplet) NM_triplet_alloc(B, block_row_size*block_col_size);
  CSparseMatrix* Btriplet = B->matrix2->triplet;
  B->matrix2->origin = NS_TRIPLET;
  double* Alocal = A;
  for (size_t i = 0, j = 0; i < (size_t)B->size0; i += block_row_size, j+= block_col_size)
  {
    {
      for (size_t col_indx = j; col_indx < block_col_size+j; ++col_indx)
      {
        for (size_t row_indx = i; row_indx < block_row_size+i; ++row_indx, ++Alocal)
        {
          CHECK_RETURN(cs_zentry(Btriplet, row_indx, col_indx, *Alocal));
        }
//        Alocal = ;
      }
//      Alocal = A + block_row_size*block_col_size;
    }
  }
}

static void computeDenseAWpB(
  double *A,
  NumericsMatrix *W,
  double *B,
  NumericsMatrix *AWpB)
{
  unsigned problemSize = W->size0;

  double* result = AWpB->matrix0;
  double* Wx = W->matrix0;
  assert(result);
  assert(Wx);

  assert(problemSize >= 3);


  double Wij[9], Ai[9], Bi[9], tmp[9];

  for (unsigned int ip3 = 0, ip9 = 0; ip3 < problemSize; ip3 += 3, ip9 += 9)
  {
    assert(ip9 < 3 * problemSize - 8);

    extract3x3(3, ip9, 0, A, Ai);
    extract3x3(3, ip9, 0, B, Bi);

    for (unsigned int jp3 = 0; jp3 < problemSize; jp3 += 3)
    {
      assert(jp3 < problemSize - 2);
      assert(ip3 < problemSize - 2);

      extract3x3(problemSize, ip3, jp3, Wx, Wij);
      mm3x3(Ai, Wij, tmp);
      if (jp3 == ip3) add3x3(Bi, tmp);
      insert3x3(problemSize, ip3, jp3, result, tmp);

    }
  }
}

static void computeSparseBlockAWpB(
  double *A,
  NumericsMatrix *W,
  double *B,
  NumericsMatrix *AWpB)
{

  /* unsigned int problemSize = W->size0; */

  SparseBlockStructuredMatrix* Wb = W->matrix1;
  SparseBlockStructuredMatrix* result = AWpB->matrix1;
  assert(Wb);
  assert(result);

  /* Check for not allocated matrix  */
  if (result->nbblocks != Wb->nbblocks)
  {
    NM_copy(W, AWpB);
  }

  assert((unsigned)W->size0 >= 3);
  assert((unsigned)W->size0 / 3 >= Wb->filled1 - 1);

  double Ai[9], Bi[9], tmp[9];

  for (unsigned int row = 0, ip9 = 0, i0 = 0;
       row < Wb->filled1 - 1; ++row, ip9 += 9, i0 += 3)
  {
    assert(ip9 < 3 *  (unsigned)W->size0 - 8);

    extract3x3(3, ip9, 0, A, Ai);
    extract3x3(3, ip9, 0, B, Bi);

    for (unsigned int blockn = (unsigned int) Wb->index1_data[row];
         blockn < Wb->index1_data[row + 1]; ++blockn)
    {

      unsigned int col = (unsigned int) Wb->index2_data[blockn];

      mm3x3(Ai, Wb->block[blockn], tmp);
      if (col == row) add3x3(Bi, tmp);

      cpy3x3(tmp, result->block[blockn]);
    }
  }
  /* Invalidation of sparse storage, if any. */
  NM_clearSparseStorage(AWpB);
}

static void computeSparseAWpB(
  double *A,
  NumericsMatrix *W,
  double *B,
  NumericsMatrix *AWpB)
{
  unsigned problemSize = W->size0;
  assert(problemSize >= 3);

  assert(AWpB->matrix2);
  assert(W->matrix2);


  NM_clearSparseStorage(AWpB);
  NumericsMatrix* Amat = NM_create(NM_SPARSE, problemSize, problemSize);
  NumericsMatrix* Bmat = NM_create(NM_SPARSE, problemSize, problemSize);

  NM_dense_to_sparse_diag_t(A, Amat, 3, 3);
  NM_dense_to_sparse_diag_t(B, Bmat, 3, 3);

  /* AWpB = B */
  NM_copy(Bmat, AWpB);

  /*  AWpB += AW */
  NM_gemm(1., Amat, W, 1., AWpB);

  NM_free(Amat);
  NM_free(Bmat);
  free(Amat);
  free(Bmat);
}

void computeAWpB(
  double *A,
  NumericsMatrix *W,
  double *B,
  NumericsMatrix *AWpB)
{
  switch (W->storageType)
  {
  case NM_DENSE:
  {
    computeDenseAWpB(A, W, B, AWpB);
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    computeSparseBlockAWpB(A, W, B, AWpB);
    break;
  }
  case NM_SPARSE:
  {
    computeSparseAWpB(A, W, B, AWpB);
    break;
  }
  default:
  {
    printf("computeAWpB :: Unsupported storage type %d, exiting!\n", W->storageType);
    exit(EXIT_FAILURE);
  }
  }
}

int globalLineSearchGP(
  fc3d_nonsmooth_Newton_solvers* equation,
  double *reaction,
  double *velocity,
  double *mu,
  double *rho,
  double *F,
  double *A,
  double *B,
  NumericsMatrix *W,
  double *qfree,
  NumericsMatrix *AWpB,
  double *direction,
  double *tmp,
  double alpha[1],
  unsigned int maxiter_ls);
int globalLineSearchGP(
  fc3d_nonsmooth_Newton_solvers* equation,
  double *reaction,
  double *velocity,
  double *mu,
  double *rho,
  double *F,
  double *A,
  double *B,
  NumericsMatrix *W,
  double *qfree,
  NumericsMatrix *AWpB,
  double *direction,
  double *tmp,
  double alpha[1],
  unsigned int maxiter_ls)
{
  unsigned problemSize = W->size0;

  double inf = 1e10;
  double alphamin = 0.0;
  double alphamax = inf;

  double m1 = 0.01, m2 = 0.99;

  // Computation of q(t) and q'(t) for t =0

  double q0 = 0.5 * cblas_ddot(problemSize, F, 1, F, 1);

  if (isnan(q0) || isinf(q0))
  {
    if (verbose > 0)
    {
      fprintf(stderr, "global line search warning. q0 is not a finite number.\n");
    }
    return -1;
  }

  //  tmp <- AWpB * direction
  NM_gemv(1., AWpB, direction, 0., tmp);

  double dqdt0 = cblas_ddot(problemSize, F, 1, tmp, 1);

  for (unsigned int iter = 0; iter < maxiter_ls; ++iter)
  {

    // tmp <- alpha*direction+reaction
    cblas_dcopy(problemSize, reaction, 1, tmp, 1);
    cblas_daxpy(problemSize, alpha[0], direction, 1, tmp, 1);

    // velocity <- W*tmp + qfree
    cblas_dcopy(problemSize, qfree, 1, velocity, 1);
    NM_gemv(1., W, tmp, 1., velocity);

    equation->function(equation->data, problemSize, tmp,
                       velocity, mu, rho, F,
                       NULL, NULL);

    double q  = 0.5 * cblas_ddot(problemSize, F, 1, F, 1);

    if (isnan(q) || isinf(q))
    {
      printf("global line search warning. q is not a finite number.\n");
      return -1;
    }

    assert(q >= 0);

    double slope = (q - q0) / alpha[0];

    int C1 = (slope >= m2 * dqdt0);
    int C2 = (slope <= m1 * dqdt0);

    if (C1 && C2)
    {
      if (verbose > 0)
      {
        printf("             globalLineSearchGP. success. ls_iter = %i  alpha = %.10e, q = %.10e\n", iter, alpha[0], q);
      }

      return 0;

    }
    else if (!C1)
    {
      alphamin = alpha[0];
    }
    else
    {
      // not(C2)
      alphamax = alpha[0];
    }

    if (alpha[0] < inf)
    {
      alpha[0] = 0.5 * (alphamin + alphamax);
    }
    else
    {
      alpha[0] = alphamin;
    }

  }
  if (verbose > 0)
  {
    printf("global line search reached the  max number of iteration  = %i  with alpha = %.10e \n", maxiter_ls, alpha[0]);
  }

  return -1;
}

void fc3d_FischerBurmeisterFunctionGenerated(
  double *reaction,
  double *velocity,
  double mu,
  double *rho,
  double *f,
  double *A,
  double *B);


void fc3d_FischerBurmeisterGradFMeritGenerated(
  double rn,
  double rt1,
  double rt2,
  double un,
  double ut1,
  double ut2,
  double mu,
  double rhon,
  double rhot1,
  double rhot2,
  double *result);

void fc3d_FischerBurmeisterGradMeritFunctionGenerated(
  double *reaction,
  double *velocity,
  double mu,
  double *rho,
  double *gf);


/* cf Fachicchinei & Pang, Finite-Dimensional Variational Inequalities
 * and Complementarity Problems, Volume II, p 805. */
int frictionContactFBLSA(
  fc3d_nonsmooth_Newton_solvers* equation,
  double *reaction,
  double *velocity,
  double *mu,
  double *rho,
  double *F,
  double *A,
  double *B,
  NumericsMatrix *W,
  double *qfree,
  NumericsMatrix *blockAWpB,
  double *direction,
  double *tmp,
  double alpha[1],
  unsigned int maxiter_ls);
int frictionContactFBLSA(
  fc3d_nonsmooth_Newton_solvers* equation,
  double *reaction,
  double *velocity,
  double *mu,
  double *rho,
  double *F,
  double *A,
  double *B,
  NumericsMatrix *W,
  double *qfree,
  NumericsMatrix *blockAWpB,
  double *direction,
  double *tmp,
  double alpha[1],
  unsigned int maxiter_ls)
{
  unsigned problemSize = W->size0;

  // notes :
  // - F contains FB or grad FB merit
  // - tmp contains direction, scal*direction, reaction+scal*direction

  // cf Newton_methods.c, L59
  double p = 2.1;
  double fblsa_rho = 1e-8;
  double gamma = 1e-4;
  double scal = 1.;

  // F <- compute fb
  fc3d_FischerBurmeisterFunction(problemSize,
                                              (FischerBurmeisterFun3x3Ptr) fc3d_FischerBurmeisterFunctionGenerated,
                                              reaction,
                                              velocity,
                                              mu,
                                              rho,
                                              F,
                                              NULL,
                                              NULL);

  double thetafb0 = 0.5 * cblas_ddot(problemSize, F, 1, F, 1);

  // F <- compute gradient of fb merit function (ugly)
  fc3d_FischerBurmeisterFunction(problemSize,
                                              (FischerBurmeisterFun3x3Ptr) fc3d_FischerBurmeisterGradMeritFunctionGenerated,
                                              reaction,
                                              velocity,
                                              mu,
                                              rho,
                                              F,
                                              NULL,
                                              NULL);
  double norm_dir_exp_p = pow(cblas_dnrm2(problemSize, direction, 1), p);
  double gradmeritfb_dir = cblas_ddot(problemSize, F, 1, direction, 1);

  if (!isnan(gradmeritfb_dir) && !isinf(gradmeritfb_dir) && gradmeritfb_dir > (-fblsa_rho * norm_dir_exp_p))
  {
    if (verbose > 0)
    {
      printf("fc3d FBLSA: condition 9.1.6 unsatisfied, gradmeritfb_dir=%g, norm_r=%g\n", gradmeritfb_dir, norm_dir_exp_p);
    }

    // FIX: failure...
    if (verbose > 0)
    {
      printf("fc3d FBLSA: set d^k to - grad merit(fb)\n");
    }

    cblas_dcopy(problemSize, F, 1, direction, 1);
    cblas_dscal(problemSize, -1, direction, 1);
  }

  for (unsigned int iter = 0; iter < maxiter_ls; ++iter)
  {

    scal /= 2.;

    // tmp <- 2^(-ik)*direction+reaction
    cblas_dcopy(problemSize, reaction, 1, tmp, 1);
    cblas_daxpy(problemSize, scal, direction, 1, tmp, 1);

    // velocity <- W*tmp + qfree
    cblas_dcopy(problemSize, qfree, 1, velocity, 1);
    NM_gemv(1., W, tmp, 1., velocity);

    // compute fb
    fc3d_FischerBurmeisterFunction(problemSize,
                                                (FischerBurmeisterFun3x3Ptr) fc3d_FischerBurmeisterFunctionGenerated,
                                                tmp,
                                                velocity,
                                                mu,
                                                rho,
                                                F,
                                                NULL,
                                                NULL);

    double thetafb  = 0.5 * cblas_ddot(problemSize, F, 1, F, 1);

    // compute grad merit fb (ugly)
    fc3d_FischerBurmeisterFunction(problemSize,
                                                (FischerBurmeisterFun3x3Ptr) fc3d_FischerBurmeisterGradMeritFunctionGenerated,
                                                tmp,
                                                velocity,
                                                mu,
                                                rho,
                                                F,
                                                NULL,
                                                NULL);

    // tmp <- scal*direction
    cblas_dscal(problemSize, 0., tmp, 1);
    cblas_daxpy(problemSize, scal, direction, 1, tmp, 1);
    double grad_meritf_reaction = cblas_ddot(problemSize, F, 1, tmp, 1);

    if (!isinf(grad_meritf_reaction) && !isnan(grad_meritf_reaction) &&
        thetafb < thetafb0 + gamma * scal * grad_meritf_reaction)
    {
      if (verbose > 0)
      {
        printf("fc3d FBLSA success. iteration  = %i, thetafb=%g, thetafb0=%g, gradmeritf,reaction=%g\n", iter, thetafb, thetafb0, gamma*scal*grad_meritf_reaction);
      }
      // tmp <- reaction + tmp
      cblas_daxpy(problemSize, 1, reaction, 1, tmp, 1);

      return 0;
    }
  }

  if (verbose > 0)
  {
    printf("fc3d FBLSA reached the max number of iteration reached  = %i\n", maxiter_ls);
  }

  return -1;
}


void fc3d_nonsmooth_Newton_solvers_solve(fc3d_nonsmooth_Newton_solvers* equation,
                                      double* reaction,
                                      double* velocity,
                                      int* info,
                                      SolverOptions* options)
{


  assert(equation);
  /* verbose=1; */
  FrictionContactProblem* problem = equation->problem;

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

  assert(problem->M->matrix0 || problem->M->matrix1 || problem->M->matrix2);

  assert(!options->iparam[4]); // only host

  unsigned int problemSize = 3 * problem->numberOfContacts;

  unsigned int iter = 0;
  unsigned int itermax = options->iparam[0];
  unsigned int erritermax = options->iparam[7];

  assert(itermax > 0);


  double tolerance = options->dparam[0];
  assert(tolerance > 0);

  if (verbose > 0)
    printf("---- FC3D - _nonsmooth_Newton_solversSolve - Start with tolerance = %g\n", tolerance);

  unsigned int _3problemSize = 3 * problemSize;
  double norm_q = cblas_dnrm2(problemSize , problem->q , 1);

  void *buffer;

  if (!options->dWork)
  {
    buffer = calloc((11 * problemSize) , sizeof(double)); // F(1),
                                                          // tmp1(1),
                                                          // tmp2(1),
                                                          // tmp3(1),
                                                          // A(3),
                                                          // B(3), rho
  }
  else
  {
    buffer = options->dWork;
  }
  double *F = (double *) buffer;
  double *tmp1 = (double *) F + problemSize;
  double *tmp2 = (double *) tmp1 + problemSize;
  double *tmp3 = (double *) tmp2 + problemSize;
  double *Ax = tmp3 + problemSize;
  double *Bx = Ax + _3problemSize;
  double *rho = Bx + _3problemSize;

  NumericsMatrix *AWpB;
  if (!options->dWork)
  {
    AWpB = NM_create(problem->M->storageType,
        problem->M->size0, problem->M->size1);
  }
  else
  {
    AWpB = (NumericsMatrix*) (rho + problemSize);
  }

  /* just for allocations */
  NM_copy(problem->M, AWpB);

  if (problem->M->storageType != NM_DENSE)
  {
    switch(options->iparam[13])
    {
      case 0:
        {
          NM_linearSolverParams(AWpB)->solver = NS_CS_LUSOL;
          break;
        }
      case 1:
        {
          NM_linearSolverParams(AWpB)->solver = NS_MUMPS;

#ifdef HAVE_MPI

          assert (options->solverData);

          if ((MPI_Comm) options->solverData == MPI_COMM_NULL)
          {
            options->solverData = NM_MPI_com(MPI_COMM_NULL);
          }
          else
          {
            NM_MPI_com((MPI_Comm) options->solverData);
          }

#endif
          break;
        }
      default:
        {
          numerics_error("fc3d_nonsmooth_Newton_solvers_solve", "Unknown linear solver.\n");
        }
    }
  }

  // compute rho here
  FrictionContactProblem * localproblem =fc3d_local_problem_allocate(problem);
  assert(options->dparam[SICONOS_FRICTION_3D_NSN_RHO]>0.0);
  for (int contact = 0; contact < problem->numberOfContacts; ++contact)
  {
    if (options->iparam[SICONOS_FRICTION_3D_NSN_RHO_STRATEGY] == SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM_COND)
    {
      fc3d_local_problem_fill_M(problem, localproblem, contact);
      compute_rho_split_spectral_norm_cond(localproblem, &rho[3*contact]);
    }
    else if (options->iparam[SICONOS_FRICTION_3D_NSN_RHO_STRATEGY] == SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM)
    {
      fc3d_local_problem_fill_M(problem, localproblem, contact);
      compute_rho_split_spectral_norm(localproblem, &rho[3*contact]);
    }
    else if (options->iparam[SICONOS_FRICTION_3D_NSN_RHO_STRATEGY] == SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM)
    {
      fc3d_local_problem_fill_M(problem, localproblem, contact);
      compute_rho_spectral_norm(localproblem, &rho[3*contact]);
    }
    else if (options->iparam[SICONOS_FRICTION_3D_NSN_RHO_STRATEGY] == SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_CONSTANT)
    {
      rho[3*contact] = options->dparam[SICONOS_FRICTION_3D_NSN_RHO];
      rho[3*contact+1] = options->dparam[SICONOS_FRICTION_3D_NSN_RHO];
      rho[3*contact+2] = options->dparam[SICONOS_FRICTION_3D_NSN_RHO];
    }
    else if (options->iparam[SICONOS_FRICTION_3D_NSN_RHO_STRATEGY] == SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_ADAPTIVE)
    {
      numerics_error("fc3d_nonsmooth_Newton_solvers_solve", "Adaptive strategy for computing rho not yet implemented");
    }
    else
      numerics_error("fc3d_nonsmooth_Newton_solvers_solve", "unknown strategy for computing rho");
    numerics_printf("fc3d_AC_initialize""contact = %i, rho[0] = %4.2e, rho[1] = %4.2e, rho[2] = %4.2e", contact, rho[3*contact], rho[3*contact+1], rho[3*contact+2]);

  }
  

  // velocity <- M*reaction + qfree
  cblas_dcopy(problemSize, problem->q, 1, velocity, 1);
  NM_gemv(1., problem->M, reaction, 1., velocity);

  double linear_solver_residual=0.0;

  while (iter++ < itermax)
  {

    equation->function(equation->data,
                       problemSize,
                       reaction, velocity, equation->problem->mu,
                       rho,
                       F, Ax, Bx);
    // AW + B
    computeAWpB(Ax, problem->M, Bx, AWpB);

    cblas_dcopy_msan(problemSize, F, 1, tmp1, 1);
    cblas_dscal(problemSize, -1., tmp1, 1);

    /* Solve: AWpB X = -F */
//    NM_copy(AWpB, AWpB_backup);
    int lsi = NM_gesv(AWpB, tmp1, true);

    /* NM_copy needed here */
//    NM_copy(AWpB_backup, AWpB);

    if (lsi)
    {
      if (verbose > 0)
      {
        numerics_warning("fc3d_nonsmooth_Newton_solvers_solve -",
                         "warning! linear solver exit with code = %d\n", lsi);
      }
    }

    if (verbose > 0)
    {
      cblas_dcopy_msan(problemSize, F, 1, tmp3, 1);
      NM_gemv(1., AWpB, tmp1, 1., tmp3);
      linear_solver_residual = cblas_dnrm2(problemSize, tmp3, 1);
      /* fprintf(stderr, "fc3d esolve: linear equation residual = %g\n", */
      /*         cblas_dnrm2(problemSize, tmp3, 1)); */
      /* for the component wise scaled residual: cf mumps &
       * http://www.netlib.org/lapack/lug/node81.html */
    }
    // line search
    double alpha = 1;
    int info_ls = 0;

    cblas_dcopy_msan(problemSize, tmp1, 1, tmp3, 1);

    switch (options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH])
    {
    case SICONOS_FRICTION_3D_NSN_LINESEARCH_NO:
      /* without line search */
      info_ls = 1;
      break;

    case SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE:
      /* Goldstein Price */
      info_ls = globalLineSearchGP(equation, reaction, velocity, problem->mu, rho, F, Ax, Bx, problem->M, problem->q, AWpB, tmp1, tmp2, &alpha, options->iparam[12]);
      break;
    case SICONOS_FRICTION_3D_NSN_LINESEARCH_ARMIJO:
      /* FBLSA */
      info_ls = frictionContactFBLSA(equation, reaction, velocity, problem->mu, rho, F, Ax, Bx,
                                     problem->M, problem->q, AWpB, tmp1, tmp2, &alpha, options->iparam[12]);
      break;
    default:
    {
      numerics_error("fc3d_nonsmooth_Newton_solvers_solve",
                     "Unknown line search option.\n");
    }
    }

    if (!info_ls)
      // tmp2 should contains the reaction iterate of the line search
      //  for GP this should be the same as cblas_daxpy(problemSize, alpha, tmp1, 1, reaction, 1);
      cblas_dcopy(problemSize, tmp2, 1, reaction, 1);
    else
      cblas_daxpy(problemSize, 1., tmp3, 1., reaction, 1);

    // velocity <- M*reaction + qfree
    cblas_dcopy(problemSize, problem->q, 1, velocity, 1);
    NM_gemv(1., problem->M, reaction, 1., velocity);

    options->dparam[1] = INFINITY;

    if (!(iter % erritermax))
    {

      fc3d_compute_error(problem, reaction, velocity,
//      fc3d_FischerBurmeister_compute_error(problem, reaction, velocity,
                         tolerance, options, norm_q, &(options->dparam[1]));

      DEBUG_EXPR_WE(equation->function(equation->data, problemSize,
                                       reaction, velocity, equation->problem->mu, rho,
                                       F, NULL, NULL));


      DEBUG_EXPR_WE(assert((cblas_dnrm2(problemSize, F, 1)
                            / (1 + cblas_dnrm2(problemSize, problem->q, 1)))
                           <= (10 * options->dparam[1] + 1e-15)));

    }

    if (verbose > 0)
    {
      equation->function(equation->data, problemSize,
                         reaction, velocity, equation->problem->mu, rho,
                         F, NULL, NULL);

      printf("   ---- fc3d_nonsmooth_Newton_solvers_solve: iteration %d : , linear solver residual =%g, residual=%g, ||F||=%g\n", iter, linear_solver_residual, options->dparam[1],cblas_dnrm2(problemSize, F, 1));
    }

    if (options->callback)
    {
      options->callback->collectStatsIteration(options->callback->env, problemSize, reaction, velocity,
                                      options->dparam[1], NULL);
    }

    if (isnan(options->dparam[1]))
    {
       if (verbose > 0)
       {
         printf("            fc3d_nonsmooth_Newton_solvers_solve: iteration %d : computed residual is not a number, stop.\n", iter);
       }
       info[0] = 2;
       break;
    }

    if (options->dparam[1] < tolerance)
    {
      info[0] = 0;
      break;
    }

  }

  if (verbose > 0)
  {
    if (!info[0])
      printf("---- FC3D - NSN - convergence after %d iterations, residual : %g < %g \n",  iter, options->dparam[1],tolerance);
    else
    {
      printf("---- FC3D - NSN - no convergence after %d iterations, residual : %g  < %g \n",  iter, options->dparam[1], tolerance);
    }
  }

  options->iparam[SICONOS_IPARAM_ITER_DONE] = iter;
  if (problem->M->storageType == NM_SPARSE_BLOCK)
  {
    /* we release the pointer to avoid deallocation of the diagonal blocks of the original matrix of the problem*/
    localproblem->M->matrix0 = NULL;
  }
  freeFrictionContactProblem(localproblem);
  
  if (!options->dWork)
  {
    assert(buffer);
    free(buffer);
    options->dWork = NULL;
  }
  else
  {
    assert(buffer == options->dWork);
  }

  if (!options->dWork)
  {
    NM_free(AWpB);

    free(AWpB);
  }
  if (verbose > 0)
    printf("---- FC3D - NSN - End\n");

}
