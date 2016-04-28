/* Factorisation with Newton_Methods.c is needed */

#include "fc3d_nonsmooth_Newton_solvers.h"

#include "NumericsMatrix_private.h"

#define DEBUG_MESSAGES 1
#include "debug.h"
#include "op3x3.h"
#include "SparseBlockMatrix.h"
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

#include "sanitizer.h"

void computeDenseAWpB(
  double *A,
  NumericsMatrix *W,
  double *B,
  NumericsMatrix *AWpB);
void computeDenseAWpB(
  double *A,
  NumericsMatrix *W,
  double *B,
  NumericsMatrix *AWpB)
{
  unsigned problemSize = W->size0;

  double* result = AWpB->matrix0;

  double* Wx = W->matrix0;

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

void computeSparseAWpB(
  double *A,
  NumericsMatrix *W,
  double *B,
  NumericsMatrix *AWpB);

void computeSparseAWpB(
  double *A,
  NumericsMatrix *W,
  double *B,
  NumericsMatrix *AWpB)
{

  /* unsigned int problemSize = W->size0; */

  SparseBlockStructuredMatrix* Wb = W->matrix1;

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

      cpy3x3(tmp, AWpB->matrix1->block[blockn]);
    }
  }
  /* Invalidation of sparse storage, if any. */
  if (AWpB->matrix2)
  {
    if (AWpB->matrix2->triplet)
    {
      cs_spfree(AWpB->matrix2->triplet);
      AWpB->matrix2->triplet = NULL;
    }
    if (AWpB->matrix2->csc)
    {
      cs_spfree(AWpB->matrix2->csc);
      AWpB->matrix2->csc = NULL;
    }
    if (AWpB->matrix2->trans_csc)
    {
      cs_spfree(AWpB->matrix2->trans_csc);
      AWpB->matrix2->trans_csc = NULL;
    }
  }
}

void computeAWpB(
  double *A,
  NumericsMatrix *W,
  double *B,
  NumericsMatrix *AWpB);
void computeAWpB(
  double *A,
  NumericsMatrix *W,
  double *B,
  NumericsMatrix *AWpB)
{
  if (W->storageType == NM_DENSE)
  {
    computeDenseAWpB(A, W, B, AWpB);
  }
  else
  {
    assert (W->storageType == NM_SPARSE_BLOCK);
    computeSparseAWpB(A, W, B, AWpB);
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

  // cf Newton_Methods.c, L59
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

  int nzmax;

  if (problem->M->storageType == NM_DENSE)
  {
    nzmax = problemSize * problemSize;
  }
  else
  {
    nzmax = options->iparam[3];
  }

  assert(itermax > 0);
  assert(nzmax > 0);

  double tolerance = options->dparam[0];
  assert(tolerance > 0);
  if (verbose > 0)
    printf("------------------------ FC3D - _nonsmooth_Newton_solversSolve - Start with tolerance = %g\n", tolerance);

  unsigned int _3problemSize = 3 * problemSize;

  void *buffer;

  if (!options->dWork)
  {
    buffer = malloc((11 * problemSize) * sizeof(double)); // F(1),
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

  NumericsMatrix *AWpB, *AWpB_backup;
  if (!options->dWork)
  {
    AWpB = createNumericsMatrix(problem->M->storageType,
        problem->M->size0, problem->M->size1);

    AWpB_backup = createNumericsMatrix(problem->M->storageType,
        problem->M->size0, problem->M->size1);
  }
  else
  {
    AWpB = (NumericsMatrix*) (rho + problemSize);
    AWpB_backup = (NumericsMatrix*) (AWpB + sizeof(NumericsMatrix*));
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
          fprintf(stderr, "fc3d esolve: unknown linear solver.\n");
          return;
        }
    }
  }

  // compute rho here
  for (unsigned int i = 0; i < problemSize; ++i) rho[i] = options->dparam[3];

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
    NM_copy(AWpB, AWpB_backup);
    int lsi = NM_gesv(AWpB, tmp1);

    /* NM_copy needed here */
    NM_copy(AWpB_backup, AWpB);

    if (lsi)
    {
      if (verbose > 0)
      {
        fprintf(stderr, "fc3d esolve: warning! linear solver exit with code = %d\n", lsi);
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

    switch (options->iparam[11])
    {
    case -1:
      /* without line search */
      info_ls = 1;
      break;

    case 0:
      /* Goldstein Price */
      info_ls = globalLineSearchGP(equation, reaction, velocity, problem->mu, rho, F, Ax, Bx, problem->M, problem->q, AWpB, tmp1, tmp2, &alpha, options->iparam[12]);
      break;
    case 1:
      /* FBLSA */
      info_ls = frictionContactFBLSA(equation, reaction, velocity, problem->mu, rho, F, Ax, Bx,
                                     problem->M, problem->q, AWpB, tmp1, tmp2, &alpha, options->iparam[12]);
      break;
    default:
      {
        fprintf(stderr, "fc3d esolve: unknown line search option.\n");
        return;
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
                                      tolerance, options, &(options->dparam[1]));

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
      printf("------------------------ FC3D - NSN - convergence after %d iterations, residual : %g < %g \n",  iter, options->dparam[1],tolerance);
    else
    {
      printf("------------------------ FC3D - NSN - no convergence after %d iterations, residual : %g  < %g \n",  iter, options->dparam[1], tolerance);
    }
  }

  options->iparam[1] = iter;

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
    freeNumericsMatrix(AWpB);
    freeNumericsMatrix(AWpB_backup);

    free(AWpB);
    free(AWpB_backup);
  }
  if (verbose > 0)
    printf("------------------------ FC3D - NSN - End\n");

}
