#include "FrictionContactNonsmoothEqn.h"

#include "op3x3.h"
#include "SparseBlockMatrix.h"
#include "FrictionContact3D_Solvers.h"
#include "FrictionContactProblem.h"
#include "FrictionContact3D_compute_error.h"
#include "AlartCurnierGenerated.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "Friction_cst.h"
#include "SiconosLapack.h"
#include "FrictionContactNonsmoothEqn.h"
void computeAWpB(
  unsigned int problemSize,
  double *A,
  double *W,
  double *B,
  double *result)
{
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

      extract3x3(problemSize, ip3, jp3, W, Wij);
      mm3x3(Ai, Wij, tmp);
      if (jp3 == ip3) add3x3(Bi, tmp);
      insert3x3(problemSize, ip3, jp3, result, tmp);

#ifdef VERBOSE_DEBUG_1
      if (jp3 == ip3)
      {
        printf("Ai\n");
        print3x3(Ai);

        printf("Bi\n");
        print3x3(Bi);

        printf("Wij");
        print3x3(Wij);

        printf("result\n");
        print3x3(tmp);
      }
#endif

    }
  }
}

/* dense => merge with sparse one */
int globalLineSearchGP(
  unsigned int problemSize,
  FrictionContactNSFun3x3Ptr computeACFun3x3,
  double *reaction,
  double *velocity,
  double *mu,
  double *rho,
  double *F,
  double *A,
  double *B,
  double *W,
  double *qfree,
  double *AWpB,
  double *direction,
  double *tmp,
  double alpha[1],
  unsigned int maxiter_ls)
{
  double inf = 1e10;
  double alphamin = 0.0;
  double alphamax = inf;

  double m1 = 0.01, m2 = 0.99;

  // Computation of q(t) and q'(t) for t =0

  double q0 = 0.5 * cblas_ddot(problemSize, F, 1, F, 1);

  // useless (already computed)
  computeAWpB(problemSize, A, W, B, AWpB);

  //  tmp <- AWpB * direction
  cblas_dgemv(CblasColMajor,CblasNoTrans, problemSize, problemSize, 1., AWpB, problemSize, direction, 1, 0., tmp, 1);

  double dqdt0 = cblas_ddot(problemSize, F, 1, tmp, 1);


  for (unsigned int iter = 0; iter < maxiter_ls; ++iter)
  {

    // tmp <- alpha*direction+reaction
    cblas_dcopy(problemSize, reaction, 1, tmp, 1);
    cblas_daxpy(problemSize, alpha[0], direction, 1, tmp, 1);

    // velocity <- W*tmp + qfree
    cblas_dcopy(problemSize, qfree, 1, velocity, 1);
    cblas_dgemv(CblasColMajor,CblasNoTrans, problemSize, problemSize, 1.,
          W, problemSize, tmp, 1, 1., velocity, 1);

    frictionContact3D_FischerBurmeisterFunction(problemSize, computeACFun3x3, tmp,
        velocity, mu, rho, F,
        NULL, NULL);

    double q  = 0.5 * cblas_ddot(problemSize, F, 1, F, 1);

    assert(q >= 0);

    double slope = (q - q0) / alpha[0];

    int C1 = (slope >= m2 * dqdt0);
    int C2 = (slope <= m1 * dqdt0);

    if (C1 && C2)
    {
      if (verbose > 1)
      {
        printf("global line search success. Number of iteration = %i  alpha = %.10e, q = %.10e\n", iter, alpha[0], q);
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
    printf("global line search failed. max number of iteration reached  = %i  with alpha = %.10e \n", maxiter_ls, alpha[0]);
  }

  return -1;
}

#ifdef WITH_MUMPS
#include <mpi.h>
#include <dmumps_c.h>

#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654
#define ICNTL(I) icntl[(I)-1]
#define CNTL(I) cntl[(I)-1]


void frictionContactNonsmoothEqnInit(
  SolverOptions *options)
{


  DMUMPS_STRUC_C* mumps_id = (DMUMPS_STRUC_C*)malloc(sizeof(DMUMPS_STRUC_C));

  // options with void pointers ?
  options->dparam[7] = (long) mumps_id;

  // Initialize a MUMPS instance. Use MPI_COMM_WORLD.
  mumps_id->job = JOB_INIT;
  mumps_id->par = 1;
  mumps_id->sym = 0;
  mumps_id->comm_fortran = USE_COMM_WORLD;
  dmumps_c(mumps_id);

  options->iparam[8] = mumps_id->comm_fortran;

  if (verbose > 1)
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

  mumps_id->ICNTL(24) = 1; // Null pivot row detection see also CNTL(3) & CNTL(5)
  // ok for a cube on a plane & four contact points
  // computeAlartCurnierSTD != generated in this case...

  //mumps_id->CNTL(3) = ...;
  //mumps_id->CNTL(5) = ...;
}


void computeSparseAWpB(
  unsigned int problemSize,
  double *A,
  SparseBlockStructuredMatrix *W,
  double *B,
  SparseBlockStructuredMatrix *Wout,
  int nzmax MAYBE_UNUSED,
  int *nz,
  int *irn,
  int *jcn,
  double *AWpB)
{
  assert(problemSize >= 3);
  assert(problemSize / 3 >= W->filled1 - 1);

  double Ai[9], Bi[9], tmp[9];
  unsigned int _nz = 0;

  for (unsigned int row = 0, ip9 = 0, i0 = 0;
       row < W->filled1 - 1; ++row, ip9 += 9, i0 += 3)
  {
    assert(ip9 < 3 * problemSize - 8);

    extract3x3(3, ip9, 0, A, Ai);
    extract3x3(3, ip9, 0, B, Bi);

    for (unsigned int blockn = W->index1_data[row];
         blockn < W->index1_data[row + 1]; ++blockn)
    {

      unsigned int col = W->index2_data[blockn];
      unsigned int j0 = col * 3;

      mm3x3(Ai, W->block[blockn], tmp);
      if (col == row) add3x3(Bi, tmp);

      /* output in a BlockCSR for the line search */
      cpy3x3(tmp, Wout->block[blockn]);

      /* output in coo format for MUMPS */
      double* ptmp = tmp;
      unsigned int j = j0;
      for (unsigned int k = 0; k < 3; ++k)
      {
        unsigned int i = i0;
        OP3(if (*ptmp)
      {
        assert(_nz < nzmax);
          irn[_nz] = i++ + 1;
          jcn[_nz] = j + 1;
          AWpB[_nz] = *ptmp++;

          assert(irn[_nz] > 0);
          assert(jcn[_nz] > 0);

          assert(irn[_nz] <= problemSize);
          assert(jcn[_nz] <= problemSize);

          _nz++;
        }
        else
        {
          i++;
          ptmp++;
        });
        j++;
      }
      nz[0] = _nz;
      assert(irn[0]);
      assert(jcn[0]);
      assert(irn[_nz - 1] >= irn[0]);
      assert(jcn[_nz - 1] >= jcn[0]);
    }
  }
}


int globalLineSearchSparseGP(
  FrictionContactNonsmoothEqn* equation,
  unsigned int problemSize,
  double *reaction,
  double *velocity,
  double *mu,
  double *rho,
  double *F,
  double *A,
  double *B,
  SparseBlockStructuredMatrix *W,
  double *qfree,
  SparseBlockStructuredMatrix *blockAWpB,
  double *direction,
  double *tmp,
  double alpha[1],
  unsigned int maxiter_ls)
{
  double inf = 1e10;
  double alphamin = 0.0;
  double alphamax = inf;

  double m1 = 0.01, m2 = 0.99;

  // Computation of q(t) and q'(t) for t =0

  double q0 = 0.5 * cblas_ddot(problemSize, F, 1, F, 1);

  //  tmp <- AWpB * direction
  cblas_dscal(problemSize, 0., tmp, 1);
  prodSBM(problemSize, problemSize, 1., blockAWpB, direction, 1, tmp);

  double dqdt0 = cblas_ddot(problemSize, F, 1, tmp, 1);

  for (unsigned int iter = 0; iter < maxiter_ls; ++iter)
  {

    // tmp <- alpha*direction+reaction
    cblas_dcopy(problemSize, reaction, 1, tmp, 1);
    cblas_daxpy(problemSize, alpha[0], direction, 1, tmp, 1);

    // velocity <- W*tmp + qfree
    cblas_dcopy(problemSize, qfree, 1, velocity, 1);
    prodSBM(problemSize, problemSize, 1., W, tmp, 1., velocity);

    equation->function(equation->data, problemSize, tmp,
                       velocity, mu, rho, F,
                       NULL, NULL);

    double q  = 0.5 * cblas_ddot(problemSize, F, 1, F, 1); 

    assert(q >= 0);

    double slope = (q - q0) / alpha[0];

    int C1 = (slope >= m2 * dqdt0);
    int C2 = (slope <= m1 * dqdt0);

    if (C1 && C2)
    {
      if (verbose > 0)
      {
        printf("global line search success. Number of iteration = %i  alpha = %.10e, q = %.10e\n", iter, alpha[0], q);
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
    printf("global line search failed. max number of iteration reached  = %i  with alpha = %.10e \n", maxiter_ls, alpha[0]);
  }

  return -1;
}


void frictionContactNonsmoothEqnSolve(FrictionContactNonsmoothEqn* equation,
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

  assert(!problem->M->matrix0);
  assert(problem->M->matrix1);

  assert(!options->iparam[4]); // only host

  unsigned int problemSize = 3 * problem->numberOfContacts;

  unsigned int iter = 0;
  unsigned int itermax = options->iparam[0];
  unsigned int erritermax = options->iparam[7];
  int nzmax = options->iparam[3];
  DMUMPS_STRUC_C* mumps_id = (DMUMPS_STRUC_C*)(long) options->dparam[7];
  int mumps_com = options->iparam[8];

  if (mumps_com==-1)
  {
    int ierr, myid;
    int argc = 0;
    char **argv;
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    frictionContact3D_sparseLocalAlartCurnierInit(options);
    mumps_id = (DMUMPS_STRUC_C*)(long) options->dparam[7];
  }
  else /* we suppose mpi init has been done */
  {
    frictionContactNonsmoothEqnInit(options);
    mumps_id = (DMUMPS_STRUC_C*)(long) options->dparam[7];
    mumps_id->comm_fortran = mumps_com;
  }

  assert(mumps_id);

  assert(itermax > 0);
  assert(nzmax > 0);

  double tolerance = options->dparam[0];
  assert(tolerance > 0);

  unsigned int _3problemSize = 3 * problemSize;

  void *buffer;

  if (!options->dWork)
    buffer = malloc((10 * problemSize +          // F(1), tmp1(1), tmp2(1),
                                                 // A(3), B(3), rho (1)
                     nzmax) * sizeof(double) +   // AWpB
                    2 * nzmax * sizeof(int) +    // irn,  jcn
                    sizeof(SparseBlockStructuredMatrix)); // blockAWpB
  else
    buffer = options->dWork;

  double *F = (double *) buffer;
  double *tmp1 = (double *) F + problemSize;
  double *tmp2 = (double *) tmp1 + problemSize;
  double *A = tmp2 + problemSize;
  double *B = A + _3problemSize;
  double *rho = B + _3problemSize;
  double *AWpB = rho + problemSize;
  int *irn = (int *)(AWpB + nzmax);
  int *jcn = (int *)(irn + nzmax);

  SparseBlockStructuredMatrix *blockAWpB = (SparseBlockStructuredMatrix *)(jcn + nzmax);

  int nz[1];

  // compute rho here
  for (unsigned int i = 0; i < problemSize; ++i) rho[i] = options->dparam[3];

  mumps_id->n = problemSize;
  mumps_id->nz = nz[0];
  mumps_id->irn = irn;
  mumps_id->jcn = jcn;
  mumps_id->a = AWpB;
  mumps_id->rhs = tmp1;

  mumps_id->job = 6;

  info[0] = 1;


  // blockAWpB init
  copySBM(problem->M->matrix1, blockAWpB, 1);

  // velocity <- M*reaction + qfree
  cblas_dcopy(problemSize, problem->q, 1, velocity, 1);
  prodSBM(problemSize, problemSize, 1., problem->M->matrix1, reaction, 1., velocity);

  while (iter++ < itermax)
  {

    equation->function(equation->data,
                       problemSize,
                       reaction, velocity, equation->problem->mu,
                       rho,
                       F, A, B);

    // AW + B
    computeSparseAWpB(problemSize, A, problem->M->matrix1, B, blockAWpB, nzmax, nz, irn, jcn, AWpB);

    cblas_dcopy(problemSize, F, 1, tmp1, 1);
    cblas_dscal(problemSize, -1., tmp1, 1);

    /* Solve: AWpB X = -F */
    mumps_id->n = problemSize;
    mumps_id->nz = nz[0];
    mumps_id->irn = irn;
    mumps_id->jcn = jcn;
    mumps_id->a = AWpB;
    mumps_id->rhs = tmp1;

    dmumps_c(mumps_id);


    assert(mumps_id->info[0] >= 0);

    if (mumps_id->info[0] > 0)
      /*if (verbose>0)*/
      printf("fc3d_eSolve: MUMPS warning : info(1)=%d, info(2)=%d\n", mumps_id->info[0], mumps_id->info[1]);


    if (verbose > 0)
    {

      printf("mumps : condition number %g\n", mumps_id->rinfog[9]);
      printf("mumps : component wise scaled residual %g\n", mumps_id->rinfog[6]);
      printf("mumps : \n");

    }


    // line search
    double alpha = 1;
    int info_ls = globalLineSearchSparseGP(equation, problemSize, reaction, velocity, problem->mu, rho, F, A, B,
                                           problem->M->matrix1, problem->q, blockAWpB, tmp1, tmp2, &alpha, 100);

    if (!info_ls)
      cblas_daxpy(problemSize, alpha, tmp1, 1, reaction, 1);
    else
      cblas_daxpy(problemSize, 1, tmp1, 1., reaction, 1);


    // velocity <- M*reaction + qfree
    cblas_dcopy(problemSize, problem->q, 1, velocity, 1);
    prodSBM(problemSize, problemSize, 1., problem->M->matrix1, reaction, 1., velocity);



    options->dparam[1] = INFINITY;

    if (!(iter % erritermax))
    {
      equation->function(equation->data, problemSize,
                         reaction, velocity, equation->problem->mu, rho,
                         F, NULL, NULL);



      FrictionContact3D_compute_error(problem, reaction, velocity,
//      frictionContact3D_FischerBurmeister_compute_error(problem, reaction, velocity,
                                      tolerance, options, &(options->dparam[1]));


      assert((cblas_dnrm2(problemSize, F, 1)
              / (1 + cblas_dnrm2(problemSize, problem->q, 1)))
             <= (10 * options->dparam[1] + 1e-15));


    }

    if (verbose > 0)
      printf("fc3d_esolve: iteration %d : error=%g\n", iter, options->dparam[1]);


    if (options->callback)
    {
      options->callback->collectStatsIteration(options->callback->env, problemSize, reaction, velocity,
                                      options->dparam[1], NULL);
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
      printf("fc3d_esolve: convergence after %d iterations, error : %g\n",
             iter, options->dparam[1]);
    else
    {
      printf("fc3d_esolve: no convergence after %d iterations, error : %g\n",
             iter, options->dparam[1]);
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

}
#else
void frictionContactNonsmoothEqnInit(
  SolverOptions *options)
{
  fprintf(stderr, "frictionContactNonsmoothEqnInit needs -DWITH_MUMPS at compilation time");
}

void frictionContactNonsmoothEqnSolve(FrictionContactNonsmoothEqn* equation,
                                      double* reaction,
                                      double* velocity,
                                      int* info,
                                      SolverOptions* options)
{
  fprintf(stderr, "frictionContactNonsmoothEqnSolve needs -DWITH_MUMPS at compilation time");
}

#endif
