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
#include "fc3d_onecontact_nonsmooth_Newton_solvers.h"
#include "fc3d_Path.h"
#include "fc3d_NCPGlockerFixedPoint.h"
#include "fc3d_projection.h"
#include "fc3d_unitary_enumerative.h"
#include "fc3d_compute_error.h"
#include "NCP_Solvers.h"
#include "SiconosBlas.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <alloca.h>
#include "op3x3.h"
#pragma GCC diagnostic ignored "-Wmissing-prototypes"


#include "fc3d_nsgs_openmp.h"


void fc3d_nsgs_openmp_ddm_naive(FrictionContactProblem* problem, double *reaction,
                               double *velocity, int* info, SolverOptions* options)
{


#if defined(USE_OPENMP) && defined(_OPENMP)

#else
  printf("fc3d_nsgs_openmp_ddm_stupid cannot be used without openmp");
  exit(1);
#endif

  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  /* Number of contacts */
  unsigned int nc = problem->numberOfContacts;
  /* Maximum number of iterations */
  int itermax = iparam[0];
  /* Tolerance */
  double tolerance = dparam[0];

  if (*info == 0)
    return;

  if (options->numberOfInternalSolvers < 1)
  {
    numericsError("fc3d_nsgs_redblack_openmp", "The NSGS method needs options for the internal solvers, options[0].numberOfInternalSolvers should be >= 1");
  }
  assert(options->internalSolvers);

  SolverOptions * localsolver_options = options->internalSolvers;

  SolverPtr local_solver = NULL;
  Update_indexPtr update_thread_problem = NULL;
  FreeSolverNSGSPtr freeSolver = NULL;
  ComputeErrorPtr computeError = NULL;

  /* Connect local solver and local problem*/
  unsigned int max_threads = omp_get_max_threads();

  FrictionContactProblem **thread_problems = alloca(max_threads*sizeof(void*));
  SolverOptions ** thread_solver_options = alloca(max_threads*sizeof(void*));
  FrictionContactProblem **local_problems = alloca(max_threads*sizeof(void*));
  SolverOptions **local_solver_options = alloca(max_threads*sizeof(void*));


  if (iparam[10] >0)
  {
    max_threads = iparam[10];
    omp_set_num_threads(max_threads);
  }

  unsigned int p = max_threads;
  unsigned int dnp = nc/(p); // domain size
  // estimation of interface size
  unsigned int size_inter = dnp/p;// dnp/4; // dnp/(p) ;

  if (p==1)
    size_inter = dnp/(p+1);

  if (size_inter%2 == 1)
  {
    size_inter ++;
  }


  printf("dnp = %i\n", dnp);
  printf("nc = %i\n", nc);
  printf("p = %i\n", p);
  printf("size_inter = %i\n", size_inter);
  unsigned int index_size =p*size_inter;
  int n_inter=p+1;
  unsigned int * index = (unsigned int*) calloc(index_size,sizeof(unsigned int));
  unsigned int * out_index = (unsigned int*) calloc((nc-index_size),sizeof(unsigned int));
  unsigned int out_index_size = 0;

  unsigned int cmp =0;
  for(unsigned int s =0; s< size_inter/2; s++)
  {
    index[cmp] = s;
    cmp++;
  }

  for (int ii =1 ; ii < n_inter-1 ; ii++)
  {
    for(int s =-size_inter/2 ; s< (int)size_inter/2; s++)
    {
      index[cmp] = (ii)*dnp+s;
      cmp++;
    }
  }
  for(int s =(int)-size_inter/2; s <0; s++)
  {
    index[cmp] = nc+s;
    cmp++;
  }

  index_size = cmp;
  printf("index_size = %i\n", index_size);

  for (unsigned int ii =0 ; ii < index_size ; ii++)
  {
    printf("index[%i] = %i\n", ii, index[ii]);
  }
  cmp=0;
  for(unsigned int i=0; i<nc; i++)
  {
    if (i == index[cmp])
    {
      cmp++;
    }
    else
    {
      out_index[out_index_size]=i;
      out_index_size++;
    }
  }

  for (unsigned int ii =0 ; ii < out_index_size ; ii++)
  {
    printf("out_index[%i] = %i\n", ii, out_index[ii]);
  }



  if (verbose > 0) printf("----------------------------------- number of threads %i\n", omp_get_max_threads()  );
  if (verbose > 0) printf("----------------------------------- number of contacts %i\n", nc );
  double * q_k = (double *) malloc(nc*3*sizeof(double));
  int thread_itermax=10, master_itermax=10,  thread_iter_total=0;


  for (unsigned int i=0; i < max_threads; i++)
  {
    printf(" initilialization of local_problem and local solver options\n");
    FrictionContactProblem *local_problem = malloc(sizeof(FrictionContactProblem));
    local_problems[i] = local_problem;
    local_problem->numberOfContacts = 1;
    local_problem->dimension = 3;
    local_problem->q = (double*)malloc(3 * sizeof(double));
    local_problem->mu = (double*)malloc(sizeof(double));

    if (problem->M->storageType == NM_DENSE || problem->M->storageType == NM_SPARSE)
    {
      local_problem->M = createNumericsMatrixFromData(NM_DENSE, 3, 3,
                                                      malloc(9 * sizeof(double)));
    }
    else /* NM_SPARSE_BLOCK */
    {
      local_problem->M = createNumericsMatrixFromData(NM_DENSE, 3, 3, NULL);
    }

    local_solver_options[i] = malloc(sizeof(SolverOptions));
    null_SolverOptions(local_solver_options[i]);
    local_solver_options[i]->dparam = NULL;
    local_solver_options[i]->iparam = NULL;
    copy_SolverOptions(localsolver_options,local_solver_options[i]);

    fc3d_nsgs_index_initialize_local_solver(&local_solver, &update_thread_problem,
                                            (FreeSolverNSGSPtr *)&freeSolver, &computeError,
                                            problem, local_problems[i],
                                            options, local_solver_options[i]);
    printf(" initilialization of thread_problem and thread solver options\n");

    FrictionContactProblem *thread_problem = malloc(sizeof(FrictionContactProblem));
    thread_problems[i] = thread_problem;
    thread_problem->M = problem->M;
    thread_problem->numberOfContacts = nc;
    thread_problem->dimension = 3;
    thread_problem->q = q_k;
    thread_problem->mu = problem->mu;;


    thread_solver_options[i] = malloc(sizeof(SolverOptions));
    null_SolverOptions(thread_solver_options[i]);
    thread_solver_options[i]->dparam = NULL;
    thread_solver_options[i]->iparam = NULL;
    copy_SolverOptions(options, thread_solver_options[i]);
    printSolverOptions(thread_solver_options[i]);
    thread_solver_options[i]->iparam[0]=thread_itermax;
    thread_solver_options[i]->iparam[1]=1; // light error

    /* fc3d_nsgs_index_initialize_local_solver(&local_solver, &update_thread_problem, */
    /*                                   (FreeSolverNSGSPtr *)&freeSolver, &computeError, */
    /*                                   problem, thread_problems[i], */
    /*                                   options, thread_solver_options[i]); */





  }


    FrictionContactProblem *master_problem = malloc(sizeof(FrictionContactProblem));
    master_problem->M = problem->M;
    master_problem->numberOfContacts = nc;
    master_problem->dimension = 3;
    master_problem->q = q_k;
    master_problem->mu = problem->mu;;

    unsigned int ** thread_index = (unsigned int **)malloc(max_threads*sizeof(unsigned int *));
    unsigned int ** thread_out_index = (unsigned int **)malloc(max_threads*sizeof(unsigned int *));

    unsigned int * thread_index_size = (unsigned int *)malloc(max_threads*sizeof(unsigned int));
    unsigned int * thread_out_index_size = (unsigned int *)malloc(max_threads*sizeof(unsigned int ));

    int  istart, istop;

    for (unsigned int i=0; i < max_threads; i++)
    {

      thread_index[i] = (unsigned int *)malloc((dnp+1)*sizeof(unsigned int));

      istart = i*dnp;
      istop = i*dnp + dnp;
      if (i == p-1)  istop =nc;
      printf("id = %i \t, istart = %i\t istop = %i \n", i, istart, istop);
      thread_index[i] = (unsigned int *)malloc((istop-istart)*sizeof(unsigned int));
      printf("thread_index building for i = %i\n", i);

      /* contruct index_local */
      int kk=0;
      printf("thread_index[%i][%i] = %i",i,kk,thread_index[i][kk]);
      for (int jj = istart; jj < istop; jj++ )
      {

        thread_index[i][kk]=jj;
        //printf("thread_index[%i][%i] = %i\n",i,kk,thread_index[i][kk]);
        kk++;
      }


      thread_index_size[i] = kk;
      printf("thread_index_size[%i] = %i",i,thread_index_size[i]);

      kk = 0;
      thread_out_index[i] = (unsigned int *)malloc(nc-(istop-istart)*sizeof(unsigned int));
      for (int jj = 0; jj < istart; jj++ )
      {
        thread_out_index[i][kk]=jj;
        kk++;
      }
      for (unsigned int jj = istop; jj < nc; jj++ )
      {
        thread_out_index[i][kk]=jj;
        kk++;
      }
      thread_out_index_size[i] = kk;

    }
    //getchar();


  /*****  NSGS Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  double error_delta_reaction=0.0;
  double error_nat=0.0;
  unsigned int *scontacts = NULL;


  double * reaction_k = (double*)malloc(nc*3*sizeof(double));
  double * velocity_k = (double*)malloc(nc*3*sizeof(double));


  while ((iter < itermax) && (hasNotConverged > 0))
  {
    printf(" \n ################### START GLOBAL LOOP ################################ \n");
    
    ++iter;
    error_delta_reaction=0.0;
    /* Loop through the contact points */
    //cblas_dcopy( n , q , incx , velocity , incy );

    /* for (unsigned int kk=0; kk < 3*nc; kk++ ) reaction_k[kk]=reaction[kk]; */
    for (unsigned int i =0 ; i < 3*nc; i++)
    {
      reaction_k[i] = reaction[i];
    }


    #pragma omp  parallel                           \
    shared(reaction_k, velocity_k)
    {
      int id = omp_get_thread_num();
      int thread_iter=0;
      int info =1;

      for (unsigned int i = 0; i < thread_index_size[id]; i++ )
      {
        int contact = thread_index[id][i];
        fc3d_nsgs_index_computeqLocal(problem, reaction, contact,
                                      thread_out_index[id], thread_out_index_size[id],
                                      &(q_k[3*contact]) );

      }

      /* call nsgs_index  with the right  fc3d_nsgs_computeqLocal*/
      fc3d_nsgs_index(thread_problems[id],
                      reaction_k, velocity_k,
                      &info, thread_solver_options[id],
                      thread_index[id], thread_index_size[id]);
      thread_iter = thread_solver_options[id]->iparam[7];
      #pragma omp critical
      thread_iter_total += thread_iter;
      error_delta_reaction +=  thread_solver_options[id]->dparam[1];;
    }

    for (unsigned int i =0 ; i < 3*nc; i++)   reaction[i] = reaction_k[i];

    printf("----------------------------------- FC3D - NSGS NAIVE - End of thread problems after %i iterations with error_delta_reaction =%e \n", thread_iter_total, error_delta_reaction);
    printf("\n \n \n");
    /* -----------------master loop ------------------------ */

    verbose=1;
    iter=0;
    double error_delta_reaction_master=0.0;


    for (unsigned int i = 0; i < out_index_size; i++ )
    {
      int  contact = out_index[i];
      fc3d_nsgs_index_computeqLocal(problem, reaction, contact,
                                    out_index, out_index_size, &(q_k[3*contact]));
    }
    
    int master_hasNotConverged =1;
    /* final loop on interface */
    while ((iter < master_itermax) && ( master_hasNotConverged > 0 ))
    {
      ++iter;
      error_delta_reaction_master=0.0;
      #pragma omp parallel for reduction(+:error_delta_reaction_master)
      for ( unsigned int i = 0 ; i < index_size ; i++)
      {
        unsigned int tid = omp_get_thread_num();
        int contact = index[i];

        if (verbose > 1) printf("----------------------------------- Contact Number %i\n", contact);
        (*update_thread_problem)(contact, master_problem, local_problems[tid],
                                 reaction, // reaction_k ??
                                 local_solver_options[tid],
                                 index, index_size);

        local_solver_options[tid]->iparam[4] = contact;

        /* version without localreaction */
        double localreaction[3];
        {
          localreaction[0] = reaction[3 * contact+0];
          localreaction[1] = reaction[3 * contact+1];
          localreaction[2] = reaction[3 * contact+2];
        };

        (*local_solver)(local_problems[tid], localreaction,
                        local_solver_options[tid]);

        /* for( int kkk=0; kkk <3 ; kkk++) printf("localreaction[%i] = %4.2e\n", kkk,localreaction[kkk] ); */
        {
          error_delta_reaction_master += pow(reaction[3 * contact] - localreaction[0], 2) +
            pow(reaction[3 * contact + 1] - localreaction[1], 2) +
            pow(reaction[3 * contact + 2] - localreaction[2], 2);

          reaction[3 * contact+0] = localreaction[0];
          reaction[3 * contact+1] = localreaction[1];
          reaction[3 * contact+2] = localreaction[2];

        }
      }

      /* if (error_delta_reaction_master < tolerance) master_hasNotConverged = 0;*/
      error_delta_reaction_master = sqrt(error_delta_reaction_master);
      double norm_r = cblas_dnrm2(nc*3 , reaction , 1);
      if (fabs(norm_r) > DBL_EPSILON)
        error_delta_reaction_master /= norm_r;
      
      if (error_delta_reaction_master < tolerance)
      {
        master_hasNotConverged = 0;
        if (verbose > 0)
          printf("----------------------------------- FC3D - NSGS NAIVE master - Iteration %i Residual = %14.7e < %7.3e\n", iter, error_delta_reaction_master, options->dparam[0]);
      }
      else
      {
        if (verbose > 0)
          printf("----------------------------------- FC3D - NSGS NAIVE naster - Iteration %i Residual = %14.7e > %7.3e\n", iter, error_delta_reaction_master, options->dparam[0]);
      }

    }
    
    error_delta_reaction+= error_delta_reaction_master;

    if (verbose > 0)  printf("----------------------------------- FC3D - NSGS NAIVE iter = %i,  error_delta_reaction = %14.7e \n", iter, error_delta_reaction);
    if (verbose > 0)  printf("----------------------------------- FC3D - NSGS NAIVE iter = %i,  error_delta_reaction_master = %14.7e \n", iter, error_delta_reaction_master);
    //double normq = cblas_dnrm2(nc*3 , problem->q , 1);
    
    /* error_delta_reaction = sqrt(error_delta_reaction); */
    /* double norm_r = cblas_dnrm2(nc*3 , reaction , 1); */

    /* if (fabs(norm_r) > DBL_EPSILON) */
    /*   error_delta_reaction /= norm_r; */

    error = error_delta_reaction;

    if (error < tolerance)
    {
      hasNotConverged = 0;
      if (verbose > 0)
      {
        printf("----------------------------------- FC3D - NSGS NAIVE - Iteration %i Residual = %14.7e < %7.3e\n", iter, error, options->dparam[0]);
        double normq = cblas_dnrm2(nc*3 , problem->q , 1);
        (*computeError)(problem, reaction , velocity, tolerance, options, normq,  &error_nat);
        printf("----------------------------------- FC3D - NSGS NAIVE - Iteration %i Full Residual = %14.7e < %7.3e\n", iter, error_nat, options->dparam[0]);
        /* test of consistency */
        double c  = 10.0;
        if   ( (error_nat/c >=  error_delta_reaction) || (error_delta_reaction >= c *error_nat))
        {
          printf("%e %e %e   \n",error_nat/c, error_delta_reaction, c *error_nat     );
          printf(" WARNING: rel error_delta_reaction is not consistent with natural map error  \n");
        }

      }
    }
    else
    {
      if (verbose > 0)
      {
        printf("----------------------------------- FC3D - NSGS DDM NAIVE- Iteration %i Residual = %14.7e > %7.3e\n", iter, error, options->dparam[0]);
      }


    }

    *info = hasNotConverged;

    if (options->callback)
    {
      options->callback->collectStatsIteration(options->callback->env, 3 * nc,
                                               reaction, velocity,
                                               error, NULL);
    }
  }

  dparam[0] = tolerance;
  dparam[1] = error;
  iparam[7] = iter;

  /***** Free memory *****/
  for (unsigned int i=0; i < max_threads; i++)
  {
    (*freeSolver)(problem,local_problems[i],local_solver_options[i]);
    if (problem->M->storageType == NM_DENSE && local_problems[i]->M->matrix0)
    {
      free(local_problems[i]->M->matrix0);
    }
    local_problems[i]->M->matrix0 = NULL;
    freeFrictionContactProblem(local_problems[i]);
    deleteSolverOptions(local_solver_options[i]);
    free(local_solver_options[i]);
  }

  if (scontacts) /* shuffle */
  {
    free(scontacts);
  }

}
