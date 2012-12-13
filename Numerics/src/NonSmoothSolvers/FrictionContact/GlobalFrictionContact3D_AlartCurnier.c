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

#ifdef WITH_MUMPS
#include "GlobalFrictionContact3D_Solvers.h"
#include "GlobalFrictionContact3D_compute_error.h"
#include "LA.h"
#include "SparseBlockMatrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "FrictionContact3D_globalAlartCurnier.h"


#include <mpi.h>
#include <dmumps_c.h>

unsigned int sizeOfACPsiJacobian(
  SparseMatrix* M,
  SparseMatrix* H)
{
  return M->n + H->n + H->m;
}

/* compute psi function */
void ACpsi(GlobalFrictionContactProblem* problem,
           double *globalVelocity,
           double *reaction,
           double *velocity,
           double *rho,
           double *psi)
{

  unsigned int localProblemSize = problem->H->size1;

  unsigned int globalProblemSize = sizeOfACPsiJacobian(problem->M->matrix2, 
                                                       problem->H->matrix2);

  /* psi <- 
     compute -problem->M * globalVelocity + problem->H * reaction + problem->q
   ... */


  /* psi + problem->M->size0 <- 
     compute -velocity + trans(problem->H) * globalVelocity + problem->b
   ... */

  /* compute AC function */
  frictionContact3D_globalAlartCurnierFunction(localProblemSize, 
                                               reaction,
                                               velocity, problem->mu, rho, 
                                               psi+problem->M->size0+
                                               problem->H->size1,
                                               NULL, NULL);

  
  

}





/* init memory for jacobian */
void initACPsiJacobian(
  SparseMatrix* M,
  SparseMatrix* H, 
  SparseMatrix *A,
  SparseMatrix *B,
  SparseMatrix *J)
{
  /* only coordinates matrix */
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

  J->n = M->n + A->m + B->m;
  J->m = M->m + A->n + B->n;
    
  J->nzmax = M->nzmax + 2*H->nzmax + 2*A->n + A->nzmax + B->nzmax;
  J->nz = J->nzmax;

  J->p = (int*) malloc((J->nz) * sizeof(int));
  J->i = (int*) malloc((J->nz) * sizeof(int));
  J->x = (double*) malloc((J->nz) * sizeof(double));
  
  unsigned int offset = 0;

  /* M */
  for (unsigned int e=0; e<M->nz; ++e)
  {
    J->x[offset+e] = M->x[e];
    J->p[offset+e] = M->p[e];
    J->i[offset+e] = M->i[e];
  }
  offset = M->nz;

  /* H */
  for (unsigned int e=0; e<H->nz; ++e)
  {
    J->x[offset+e] = H->x[e];
    J->p[offset+e] = H->p[e] + M->m + A->m;
    J->i[offset+e] = H->i[e];
  }
  offset += H->nz;

  /* Ht */
  for (unsigned int e=0; e<H->nz; ++e)
  {
    J->x[offset+e] = H->x[e];
    J->i[offset+e] = H->p[e];
    J->p[offset+e] = H->i[e] + M->n + A->n;
  }
  offset += H->nz;

  /* I */
  for (unsigned int e=0; e<A->m; ++e)
  {
    J->x[offset+e] = 1.;
    J->p[offset+e] = e + M->n;
    J->i[offset+e] = e + M->m;
  }
  offset += A->m;

  /* A */
  for (unsigned int e=0; e<A->nz; ++e)
  {
    J->x[offset+e] = A->x[e];
    J->p[offset+e] = A->p[e] + M->n;
    J->i[offset+e] = A->i[e] + M->m;
  }
  offset += A->nz;

  /* B */
  for (unsigned int e=0; e<B->nz; ++e)
  {
    J->x[offset+e] = B->x[e];
    J->p[offset+e] = B->p[e] + M->n;
    J->i[offset+e] = B->i[e] + M->m;
  }

  J->nz = offset + B->nz;

}

/* update J with new A and B */
void updateACPsiJacobian(
  SparseMatrix* M,
  SparseMatrix* H, 
  SparseMatrix *A,
  SparseMatrix *B,
  SparseMatrix *J)
{
  /* only coordinates matrix */
  assert(M->nz>=0);
  assert(H->nz>=0);
  assert(A->nz>=0);
  assert(B->nz>=0);

  /* M square */
  assert(M->m == M->n);

  /* A & B squares */
  assert(A->m == A->n);
  assert(B->m == B->n);

  assert(A->nzmax == B->nzmax);  

  assert(J->n == M->n + A->m + B->m);
  assert(J->m == M->m + A->n + B->n);
    
  assert(J->nzmax == M->nzmax + 2*H->nzmax + 2*A->n + A->nzmax + B->nzmax);
  assert(J->nz    == J->nzmax);

  assert(J->p);
  assert(J->i);
  assert(J->x);
  
  unsigned int offset = 0;

  offset = M->nz;    /* M */
  offset += H->nz;   /* H */
  offset += H->nz;   /* Ht */
  offset += A->m;    /* I */

  /* A */
  for (unsigned int e=0; e<A->nz; ++e)
  {
    J->x[offset+e] = A->x[e];
    J->p[offset+e] = A->p[e] + M->n;
    J->i[offset+e] = A->i[e] + M->m;
  }
  offset += A->nz;

  /* B */
  for (unsigned int e=0; e<B->nz; ++e)
  {
    J->x[offset+e] = B->x[e];
    J->p[offset+e] = B->p[e] + M->n;
    J->i[offset+e] = B->i[e] + M->m;
  }

}

void SparseMatrixFrom3x3DiagBlocks(int nc, double* P, SparseMatrix* R)
{
  
  R->m = 3 * nc;
  R->n = 3 * nc;
  R->x = P;

  for (unsigned int ib = 0; ib < nc; ++ib)
  {
    for (unsigned int j = 0; j < 3; ++j)
    {
      for (unsigned int i = 0; i < 3; ++i)
      {
        R->p[ib+i+3*j] = ib * 3 + i; // CHECK  
        R->i[ib+i+3*j] = ib * 3 + j;
      }
    } 
  }
}

/* the true global line search 
   (globalAlartCurnier should be renamed in localAlartCurnier)
*/
int _globalLineSearchSparseGP(
  GlobalFrictionContactProblem *problem,
  double *solution,
  double *direction,
  double *globalVelocity,
  double *reaction, 
  double *velocity,
  double *mu,
  double *rho,
  double *F,
  double *psi,
  SparseMatrix *J,
  double *tmp,
  double alpha[1],
  unsigned int maxiter_ls)
{
  double inf = 1e10;
  double alphamin = 0.0;
  double alphamax = inf;

  double m1 = 0.01, m2 = 0.99;

  unsigned int globalProblemSize = sizeOfACPsiJacobian(problem->M->matrix2, 
                                                       problem->H->matrix2);

  // Computation of q(t) and q'(t) for t =0

  double q0 = 0.5 * DDOT(globalProblemSize, psi, 1, psi, 1);

  //  tmp <- J * direction
  DSCAL(globalProblemSize, 0., tmp, 1);

  /* to be done */
  /*prodSparseMatrixCoo(globalProblemSize, globalProblemSize, 1., J->matrix2, direction, 1, tmp);*/

  double dqdt0 = DDOT(globalProblemSize, psi, 1, tmp, 1);

  for (unsigned int iter = 0; iter < maxiter_ls; ++iter)
  {

    // tmp <- alpha*direction+solution
    DCOPY(globalProblemSize, tmp, 1, tmp, 1);
    DAXPY(globalProblemSize, alpha[0], direction, 1, tmp, 1);

    /* reaction = */

    /* velocity = */ 


    ACpsi(problem, globalVelocity, reaction, velocity, rho, psi);

    double q  = 0.5 * DDOT(globalProblemSize, psi, 1, psi, 1);

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
  if (verbose > 1)
  {
    printf("global line search failed. max number of iteration reached  = %i  with alpha = %.10e \n", maxiter_ls, alpha[0]);
  }

  return -1;
}



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
  assert(problem->M->matrix1);

  assert(!options->iparam[4]); // only host
  
  unsigned int problemSize = 3 * problem->numberOfContacts;
  
  unsigned int iter = 0;
  unsigned int itermax = options->iparam[0];
  unsigned int erritermax = options->iparam[7];
  int nzmax = options->iparam[3];
  DMUMPS_STRUC_C* mumps_id = (DMUMPS_STRUC_C*)(long) options->dparam[7];

  if (!mumps_id)
  {
    /* we suppose no mpi init has been done */
    /* if this not the case you *must* call
       frictionContact3D_sparseGlobalAlartCurnierInit yourself */
    int ierr, myid;
    int argc = 0;
    char **argv;
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    frictionContact3D_sparseGlobalAlartCurnierInit(options);
    mumps_id = (DMUMPS_STRUC_C*)(long) options->dparam[7];
  }

  assert(mumps_id);

  assert(itermax > 0);
  assert(nzmax > 0);

  double tolerance = options->dparam[0];
  assert(tolerance > 0);
  
  unsigned int _3problemSize = 3 * problemSize;
  
  unsigned int globalProblemSize = sizeOfACPsiJacobian(problem->M->matrix2, 
                                                       problem->H->matrix2);

  unsigned int localProblemSize = problem->H->size1;

  void *buffer;

  if (!options->dWork)
    buffer = malloc(20 * localProblemSize * sizeof(double) +         // F(1), 
                                                                   // A(9), B(9), 
                                                                   // rho (1)
                    2 * localProblemSize * sizeof(int) +           // irn,  jcn
                    3 * sizeof(SparseMatrix *) +         // A_,B_,J
                    3 * globalProblemSize * sizeof(double) // psi, tmp1 +
      );

  else
    buffer = options->dWork;

  double *F = (double *) buffer;
  double *A = F + problemSize;
  double *B = A + _3problemSize;
  double *rho = B + _3problemSize;
  int *irn = (int *)(rho + problemSize);
  int *jcn = irn + _3problemSize;

  SparseMatrix *A_ = (SparseMatrix *)(jcn + nzmax);
  SparseMatrix *B_ = A_ + 1;
  SparseMatrix *J  = B_ + 1;

  double * psi = (double *) (J + 1);
  double * tmp1 = psi + globalProblemSize;
  double * tmp2 = tmp1 + globalProblemSize;
  double * solution = tmp2 + globalProblemSize;

  int nz[1];

  SparseMatrixFrom3x3DiagBlocks(problem->numberOfContacts, A, A_);
  SparseMatrixFrom3x3DiagBlocks(problem->numberOfContacts, B, B_);

  A_->p = irn;
  B_->p = irn;
  A_->i = jcn;
  B_->i = jcn;

  initACPsiJacobian(problem->M->matrix2, problem->H->matrix2, A_, B_, J);

  // compute rho here
  for (unsigned int i = 0; i < localProblemSize; ++i) rho[i] = 1.;

  for (unsigned int i = 0; i < globalProblemSize; ++i) tmp1[i] = 0.;

  mumps_id->n = J->m;
  mumps_id->nz = J->nz;
  mumps_id->irn = J->p;
  mumps_id->jcn = J->i;
  mumps_id->a = J;
  mumps_id->rhs = tmp1;

  mumps_id->job = 6;

  info[0] = 1;

  while (iter++ < itermax)
  {

    frictionContact3D_globalAlartCurnierFunction(problemSize,
                                                 reaction, velocity,
                                                 problem->mu, rho,
                                                 F, A, B);
    /* J */
    updateACPsiJacobian(problem->M->matrix2,
                        problem->H->matrix2,
                        A_, B_, J);

    /* rhs = -F */
    DCOPY(problemSize, F, 1, tmp1 + problem->M->size0 + problem->H->size1, 1);
    DSCAL(problemSize, -1., tmp1 + problem->M->size0 + problem->H->size1, 1);

    /* Solve: J X = -F */
    dmumps_c(mumps_id);

    assert(mumps_id->info[0] >= 0);

    if (mumps_id->info[0] > 0)
      /*if (verbose>0)*/
      printf("GLOBALAC: MUMPS warning : info(1)=%d, info(2)=%d\n", 
             mumps_id->info[0], mumps_id->info[1]);


    if (verbose > 0)
    {

      printf("mumps : condition number %g\n", mumps_id->rinfog[9]);
      printf("mumps : component wise scaled residual %g\n", 
             mumps_id->rinfog[6]);
      printf("mumps : \n");

    }


    /* line search */
    double alpha = 1;
    int info_ls = _globalLineSearchSparseGP(problem,
                                            solution,
                                            tmp1,
                                            globalVelocity,
                                            reaction, velocity, 
                                            problem->mu, rho, F, psi, J,
                                            tmp2, &alpha, 100);
    

    if (!info_ls)
    {
      DAXPY(globalProblemSize, alpha, tmp1, 1, solution, 1);
    }
    else
    {
      DAXPY(problemSize, 1, tmp1, 1., solution, 1);
    }

    options->dparam[1] = INFINITY;

    if (!(iter % erritermax))
    {
//      frictionContact3D_globalAlartCurnierFunction(problemSize,
      //                                                  reaction, velocity,
      //                                             problem->mu, rho,
//                                                   F, NULL, NULL);



      FrictionContact3D_compute_error(problem, reaction, velocity,
                                      tolerance, options, 
                                      &(options->dparam[1]));


      //     assert((DNRM2(problemSize, F, 1)
      //       / (1 + DNRM2(problemSize, problem->q, 1)))
      //        <= (10 * options->dparam[1] + 1e-15));


    }

    if (verbose > 0)
      printf("GLOBALAC: iteration %d : error=%g\n", iter, options->dparam[1]);

    if (options->dparam[1] < tolerance)
    {
      info[0] = 0;
      break;
    }


  }

  if (verbose > 0)
  {
    if (!info[0])
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
  if (info[0])
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



#endif
