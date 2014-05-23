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

void frictionContact3D_AlartCurnierFunction(
  unsigned int problemSize,
  double *reaction,
  double *velocity,
  double *mu,
  double *rho,
  double *result,
  double *A,
  double *B)
{
  assert(reaction);
  assert(velocity);
  assert(rho);
  assert(mu);

  assert(problemSize / 3 > 0);
  assert(problemSize % 3 == 0);

  unsigned int i;
  for (i = 0; i < problemSize; i += 3)
  {

    computeAlartCurnierSTD
    (reaction,
     velocity,
     *mu,
     rho,
     result, A, B);


    // note: GENERATED_FUNCTION is different for a cube on plane & four
    // contact points (mu=0.8). (but test with random values ok see
    // AlartCurnierFunctions_test)
#ifdef COMPARE_WITH_GENERATED_FUNCTION
    double result_g[3];
    double A_g[9];
    double B_g[9];

    frictionContact3D_AlartCurnierFunctionGenerated(reaction,
        velocity,
        *mu,
        rho,
        result_g, A_g, B_g);
    if (result)
    {
      sub3(result, result_g);
      assert(hypot3(result_g) < 1e-7);
    }

    if (A)
    {
      sub3x3(A, A_g);
      assert(hypot9(A_g) < 1e-7);
    }

    if (B)
    {
      sub3x3(B, B_g);
      assert(hypot9(B_g) < 1e-7);
    }
#endif

    reaction += 3;
    velocity += 3;
    mu++;
    rho += 3;

    if (result)
      result += 3;

    if (A)
      A += 9;

    if (B)
      B += 9;

  }

}


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

int globalLineSearchGP(
  unsigned int problemSize,
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

    frictionContact3D_AlartCurnierFunction(problemSize, tmp,
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

void frictionContact3D_localAlartCurnier(
  FrictionContactProblem* problem,
  double *reaction,
  double *velocity,
  int *info,
  SolverOptions *options)
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

  if (!problem->M->matrix0)
  {
    frictionContact3D_sparseLocalAlartCurnier(
      problem,
      reaction,
      velocity,
      info,
      options);
    return;
  }

  assert(problem->M->matrix0);

  unsigned int problemSize = 3 * problem->numberOfContacts;

  unsigned int iter = 0;
  unsigned int itermax = options->iparam[0];
  unsigned int erritermax = options->iparam[7];

  assert(itermax > 0);

  double tolerance = options->dparam[0];
  assert(tolerance > 0);

  unsigned int problemSize2 = problemSize * problemSize;
  unsigned int _3problemSize = 3 * problemSize;

  void *buffer;

  if (!options->dWork)
  {
#ifndef NDEBUG
    buffer = malloc((14 * problemSize +
                     2 * problemSize2) * sizeof(double) +
                    problemSize * sizeof(int));
#else
    buffer = malloc((14 * problemSize +
                     problemSize2) * sizeof(double) +
                    problemSize * sizeof(int));
#endif
  }
  else
    buffer = options->dWork;

  double *F = (double *) buffer; //malloc(problemSize*sizeof(double));
  double *tmp1 = (double *) F + problemSize; //malloc(problemSize*sizeof(double));
  double *tmp2 = (double *) tmp1 + problemSize; //malloc(problemSize*sizeof(double));
  double *A = tmp2 + problemSize; //malloc(3*problemSize*sizeof(double));
  double *B = A + _3problemSize; //malloc(3*problemSize*sizeof(double));
  double *rho = B + _3problemSize; //malloc(problemSize*sizeof(double));
  double *AWpB = rho + problemSize;// malloc(problemSize*problemSize*sizeof(double));
  int *ipiv = (int *)(AWpB + problemSize2);  // malloc(problemSize*sizeof(int));
#ifndef NDEBUG
  double *AWpB_ = (double *) ipiv + problemSize;
#endif

  for (unsigned int i = 0; i < problemSize; ++i) rho[i] = 1.;

  info[0] = 1;

  // velocity <- M*reaction + qfree


  cblas_dcopy(problemSize, problem->q, 1, velocity, 1);
  cblas_dgemv(CblasColMajor,CblasNoTrans, problemSize, problemSize, 1.,
        problem->M->matrix0, problemSize, reaction, 1, 1., velocity, 1);


  while (iter++ < itermax)
  {

    frictionContact3D_AlartCurnierFunction(
      problemSize,
      reaction, velocity,
      problem->mu, rho,
      F, A, B);

    // AW + B
    computeAWpB(problemSize, A, problem->M->matrix0, B, AWpB);

    int info2 = 0;

    cblas_dcopy(problemSize, F, 1, tmp1, 1);
    cblas_dscal(problemSize, -1., tmp1, 1);

    if (options->iparam[2])
    {
      DGELS(LA_NOTRANS,problemSize, problemSize, 1, AWpB, problemSize,
            tmp1, problemSize, &info2);
    }
    else
    {

#ifndef NDEBUG
      cblas_dcopy(problemSize * problemSize, AWpB, 1, AWpB_, 1);
#endif

      // tmp1 <- sol (AWpB * tmp1 = -F)
      DGESV(problemSize, 1, AWpB, problemSize, ipiv,
            tmp1, problemSize, &info2);

#ifndef NDEBUG
      cblas_dcopy(problemSize, F, 1, tmp2, 1);
      cblas_dgemv(CblasColMajor,CblasNoTrans, problemSize, problemSize, 1., AWpB_,
            problemSize, tmp1, 1, 1., tmp2, 1);
      assert(cblas_ddot(problemSize, tmp2, 1, tmp2, 1) < 1e-10);
#endif

    }

    assert(info2 >= 0);

    if (info2 > 0)
      /*if (verbose>0)*/
      printf("LOCALAC: warning DGESV failed with U(%d,%d) == 0.\n", info2, info2);

    // line search
    double alpha = 1;
    int info_ls = globalLineSearchGP(problemSize, reaction, velocity, problem->mu, rho, F, A, B,
                                     problem->M->matrix0, problem->q, AWpB, tmp1, tmp2, &alpha, 100);

    if (!info_ls)
      cblas_daxpy(problemSize, alpha, tmp1, 1, reaction, 1);
    else
      cblas_daxpy(problemSize, 1, tmp1, 1., reaction, 1);


    // velocity <- M*reaction + qfree
    cblas_dcopy(problemSize, problem->q, 1, velocity, 1);
    cblas_dgemv(CblasColMajor,CblasNoTrans, problemSize, problemSize, 1.,
          problem->M->matrix0, problemSize, reaction, 1, 1., velocity, 1);

    options->dparam[1] = INFINITY;

    if (!(iter % erritermax))
    {
      frictionContact3D_AlartCurnierFunction(problemSize,
          reaction, velocity,
          problem->mu, rho,
          F, NULL, NULL);



      FrictionContact3D_compute_error(problem, reaction, velocity,
                                      tolerance, options, &(options->dparam[1]));


      assert((cblas_dnrm2(problemSize, F, 1)
              / (1 + cblas_dnrm2(problemSize, problem->q, 1)))
             <= (10 * options->dparam[1] + 1e-15));

    }

    if (verbose > 0)
      printf("LOCALAC: iteration %d : error=%g\n", iter, options->dparam[1]);

    if (options->dparam[1] < tolerance)
    {
      info[0] = 0;
      break;
    }


  }



  if (verbose > 0)
  {
    if (!info[0])
      printf("LOCALAC: convergence after %d iterations, error : %g\n",
             iter, options->dparam[1]);
    else
    {
      printf("LOCALAC: no convergence after %d iterations, error : %g\n",
             iter, options->dparam[1]);
    }
  }

#ifdef DUMP_PROBLEM
  if (info[0])
  {
    static int file_counter = 0;
    char filename[64];
    printf("LOCALAC: dumping problem\n");
    sprintf(filename, "LOCALAC_failure%d.dat", file_counter++);
    FILE* file = fopen(filename, "w");
    frictionContact_printInFile(problem, file);
    fclose(file);
  }
#endif

  if (!options->dWork)
  {
    assert(buffer);
    free(buffer);
  }
  else
  {
    assert(buffer == options->dWork);
  }


}

int frictionContact3D_AlartCurnier_setDefaultSolverOptions(
  SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the default solver options for the LOCALAC Solver\n");
  }

  options->solverId = SICONOS_FRICTION_3D_LOCALAC;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 9;
  options->dSize = 9;
  options->iparam = (int *) malloc(options->iSize * sizeof(int));
  options->dparam = (double *) malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  options->iWork = NULL;   options->callback = NULL; options->numericsOptions = NULL;
  for (unsigned int i = 0; i < 9; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 200;
  options->iparam[1] = 1;
  options->iparam[3] = 100000; /* nzmax*/
  options->iparam[5] = 1;  
  options->iparam[7] = 1;      /* erritermax */
  options->dparam[0] = 1e-3;
  options->dparam[3] = 1;      /* default rho */

  options->iparam[8] = -1;     /* mpi com fortran */
  options->internalSolvers = NULL;

  return 0;
}


#ifdef WITH_MUMPS
#include <mpi.h>
#include <dmumps_c.h>

#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654
#define ICNTL(I) icntl[(I)-1]
#define CNTL(I) cntl[(I)-1]

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

    frictionContact3D_AlartCurnierFunction(problemSize, tmp,
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

void frictionContact3D_sparseLocalAlartCurnierInit(
  SolverOptions *SO)
{
  DMUMPS_STRUC_C* mumps_id = malloc(sizeof(DMUMPS_STRUC_C));

  // SO with void pointers ?
  SO->dparam[7] = (long) mumps_id;

  // Initialize a MUMPS instance. Use MPI_COMM_WORLD.
  mumps_id->job = JOB_INIT;
  mumps_id->par = 1;
  mumps_id->sym = 0;
  mumps_id->comm_fortran = USE_COMM_WORLD;
  dmumps_c(mumps_id);

  SO->iparam[8] = mumps_id->comm_fortran;

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

void frictionContact3D_sparseLocalAlartCurnier(
  FrictionContactProblem* problem,
  double *reaction,
  double *velocity,
  int *info,
  SolverOptions *options)
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
    frictionContact3D_sparseLocalAlartCurnierInit(options);
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

    frictionContact3D_AlartCurnierFunction(problemSize,
        reaction, velocity,
        problem->mu, rho,
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
      printf("LOCALAC: MUMPS warning : info(1)=%d, info(2)=%d\n", mumps_id->info[0], mumps_id->info[1]);


    if (verbose > 0)
    {

      printf("mumps : condition number %g\n", mumps_id->rinfog[9]);
      printf("mumps : component wise scaled residual %g\n", mumps_id->rinfog[6]);
      printf("mumps : \n");

    }


    // line search
    double alpha = 1;
    int info_ls = globalLineSearchSparseGP(problemSize, reaction, velocity, problem->mu, rho, F, A, B,
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
      frictionContact3D_AlartCurnierFunction(problemSize,
          reaction, velocity,
          problem->mu, rho,
          F, NULL, NULL);



      FrictionContact3D_compute_error(problem, reaction, velocity,
                                      tolerance, options, &(options->dparam[1]));


      assert((cblas_dnrm2(problemSize, F, 1)
              / (1 + cblas_dnrm2(problemSize, problem->q, 1)))
             <= (10 * options->dparam[1] + 1e-15));


    }

    if (verbose > 0)
      printf("LOCALAC: iteration %d : error=%g\n", iter, options->dparam[1]);


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
      printf("LOCALAC: convergence after %d iterations, error : %g\n",
             iter, options->dparam[1]);
    else
    {
      printf("LOCALAC: no convergence after %d iterations, error : %g\n",
             iter, options->dparam[1]);
    }
  }

  options->iparam[1] = iter;

#ifdef DUMP_PROBLEM
  if (info[0])
  {
    static int file_counter = 0;
    char filename[64];
    printf("LOCALAC: dumping problem\n");
    sprintf(filename, "LOCALAC_failure%d.dat", file_counter++);
    FILE* file = fopen(filename, "w");
    frictionContact_printInFile(problem, file);
    fclose(file);
  }
#endif


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

#else /*WITH_MUMPS*/

void frictionContact3D_sparseLocalAlartCurnierInit(
  SolverOptions *SO)
{
  fprintf(stderr, "The sparse global Alart & Curnier solver needs -DWITH_MUMPS for the compilation of Siconos/Numerics\n");
}

void frictionContact3D_sparseLocalAlartCurnier(
  FrictionContactProblem* problem,
  double *reaction,
  double *velocity,
  int *info,
  SolverOptions *options)
{
  fprintf(stderr, "The sparse global Alart & Curnier solver needs -DWITH_MUMPS for the compilation of Siconos/Numerics\n");
}
#endif
