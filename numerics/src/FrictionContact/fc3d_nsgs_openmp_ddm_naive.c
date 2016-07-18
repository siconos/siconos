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
  double normq = cblas_dnrm2(nc*3 , problem->q , 1);
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

  FrictionContactProblem **interface_local_problems = alloca(max_threads*sizeof(void*));
  SolverOptions **interface_local_solver_options = alloca(max_threads*sizeof(void*));


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
  unsigned int interface_index_size =p*size_inter;
  int n_inter=p+1;
  unsigned int * interface_index = (unsigned int*) calloc(interface_index_size,sizeof(unsigned int));

  printf("out_size_inter = %i\n",(nc-interface_index_size));
  unsigned int * interface_out_index = (unsigned int*) calloc((nc-interface_index_size),sizeof(unsigned int));
  unsigned int cmp =0;
  for(unsigned int s =0; s< size_inter/2; s++)
  {
    interface_index[cmp] = s;
    cmp++;
  }


  for (int ii =1 ; ii < n_inter-1 ; ii++)
  {
    for(int s =-size_inter/2 ; s< (int)size_inter/2; s++)
    {
      interface_index[cmp] = (ii)*dnp+s;
      cmp++;
    }
  }
  for(int s =(int)-size_inter/2; s <0; s++)
  {
    interface_index[cmp] = nc+s;
    cmp++;
  }

  interface_index_size = cmp;
  printf("interface_index_size = %i\n", interface_index_size);

  for (unsigned int ii =0 ; ii < interface_index_size ; ii++)
  {
    if (ii%6 ==0)  printf("\n");
    printf("interface_index[%i] = %i\t", ii, interface_index[ii]);
  }
  printf("\n");
  unsigned int interface_out_index_size = 0;
  cmp=0;
  for(unsigned int i=0; i<nc; i++)
  {
    if (i == interface_index[cmp])
    {
      cmp++;
    }
    else
    {
     interface_out_index[interface_out_index_size]=i;
     printf("interface_out_index[%i] = %i\t", interface_out_index_size, interface_out_index[interface_out_index_size]);
     interface_out_index_size++;
    }
  }

  for (unsigned int ii =0 ; ii < interface_out_index_size ; ii++)
  {
    if (ii%6 ==0)  printf("\n");
    printf("interface_out_index[%i] = %i\t", ii, interface_out_index[ii]);
  }
  printf("\n");



  if (verbose > 0) printf("----------------------------------- number of threads %i\n", omp_get_max_threads()  );
  if (verbose > 0) printf("----------------------------------- number of contacts %i\n", nc );
  double * q_k = (double *) malloc(nc*3*sizeof(double));



  int thread_itermax=options->iparam[12], interface_itermax=options->iparam[13],  thread_iter_total=0;


  for (unsigned int i=0; i < max_threads; i++)
  {
    printf(" initilialization of interface_local_problem and local solver options\n");
    FrictionContactProblem *interface_local_problem = malloc(sizeof(FrictionContactProblem));
    interface_local_problems[i] = interface_local_problem;
    interface_local_problem->numberOfContacts = 1;
    interface_local_problem->dimension = 3;
    interface_local_problem->q = (double*)malloc(3 * sizeof(double));
    interface_local_problem->mu = (double*)malloc(sizeof(double));

    if (problem->M->storageType == NM_DENSE || problem->M->storageType == NM_SPARSE)
    {
      interface_local_problem->M = createNumericsMatrixFromData(NM_DENSE, 3, 3,
                                                      malloc(9 * sizeof(double)));
    }
    else /* NM_SPARSE_BLOCK */
    {
      interface_local_problem->M = createNumericsMatrixFromData(NM_DENSE, 3, 3, NULL);
    }

    interface_local_solver_options[i] = malloc(sizeof(SolverOptions));
    solver_options_nullify(interface_local_solver_options[i]);
    interface_local_solver_options[i]->dparam = NULL;
    interface_local_solver_options[i]->iparam = NULL;
    solver_options_copy(localsolver_options,interface_local_solver_options[i]);

    fc3d_nsgs_index_initialize_local_solver(&local_solver, &update_thread_problem,
                                            (FreeSolverNSGSPtr *)&freeSolver, &computeError,
                                            problem, interface_local_problems[i],
                                            options, interface_local_solver_options[i]);
    printf(" initilialization of thread_problem and thread solver options\n");

    FrictionContactProblem *thread_problem = malloc(sizeof(FrictionContactProblem));
    thread_problems[i] = thread_problem;
    thread_problem->M = problem->M;
    thread_problem->numberOfContacts = nc;
    thread_problem->dimension = 3;
    thread_problem->q = q_k;
    thread_problem->mu = problem->mu;;


    /* printSolverOptions(options); */
    /* getchar(); */

    thread_solver_options[i] = malloc(sizeof(SolverOptions));
    solver_options_nullify(thread_solver_options[i]);
    solver_options_copy(options, thread_solver_options[i]);
    thread_solver_options[i]->dparam[0] /= 10.0;

    thread_solver_options[i]->iparam[0]=thread_itermax;
    thread_solver_options[i]->iparam[1]=1; // light error

    /* fc3d_nsgs_index_initialize_local_solver(&local_solver, &update_thread_problem, */
    /*                                   (FreeSolverNSGSPtr *)&freeSolver, &computeError, */
    /*                                   problem, thread_problems[i], */
    /*                                   options, thread_solver_options[i]); */





  }


    FrictionContactProblem *interface_problem = malloc(sizeof(FrictionContactProblem));
    interface_problem->M = problem->M;
    interface_problem->numberOfContacts = nc;
    interface_problem->dimension = 3;
    interface_problem->q = q_k;
    interface_problem->mu = problem->mu;;

    unsigned int ** thread_index =     (unsigned int **)malloc(max_threads*sizeof(unsigned int *));
    unsigned int ** thread_out_index = (unsigned int **)malloc(max_threads*sizeof(unsigned int *));

    unsigned int * thread_index_size     = (unsigned int *)malloc(max_threads*sizeof(unsigned int));
    unsigned int * thread_out_index_size = (unsigned int *)malloc(max_threads*sizeof(unsigned int));

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
      for (int jj = istart; jj < istop; jj++ )
      {
        thread_index[i][kk]=jj;
        if (jj%6 ==0)  printf("\n");
        printf("thread_index[%i][%i] = %i\t",i,kk,thread_index[i][kk]);
        kk++;
      }

      thread_index_size[i] = kk;
      printf("\n thread_index_size[%i] = %i\n",i,thread_index_size[i]);

      kk = 0;
      printf("nc-(istop-istart) = %i \n", nc-(istop-istart));
      thread_out_index[i] = (unsigned int *)malloc( (nc-(istop-istart)) *sizeof(unsigned int));

      for (int jj = 0; jj < istart; jj++ )
      {
        thread_out_index[i][kk]=jj;
        kk++;
      }
      for (unsigned int jj = istop; jj < nc; jj++ )
      {
        thread_out_index[i][kk]=jj;
        if (jj%6 ==0)  printf("\n");
        printf("threadout__index[%i][%i] = %i\t",i,kk,thread_out_index[i][kk]);
        kk++;
      }
      thread_out_index_size[i] = kk;
      printf("\n thread_out_index_size[%i] = %i\n",i,thread_out_index_size[i]);

    }
    //getchar();


  /*****  NSGS Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  double error_delta_reaction=0.0;
  double error_nat=0.0;


  double * reaction_k = (double*)malloc(nc*3*sizeof(double));
  double * velocity_k = (double*)malloc(nc*3*sizeof(double));
  for (unsigned int i =0 ; i < 3*nc; i++)    reaction_k[i] = reaction[i];


  double normreaction_k = cblas_dnrm2(nc*3 , reaction_k , 1);
  printf("normreaction_k = %e \n", normreaction_k);

  while ((iter < itermax) && (hasNotConverged > 0))
  {
    if (verbose > 0) printf(" \n ################### START GLOBAL LOOP ################################ \n");

    ++iter;

    { /* thread block */
      error_delta_reaction=0.0;
      /* Loop through the contact points */
      //cblas_dcopy( n , q , incx , velocity , incy );
      /* for (unsigned int kk=0; kk < 3*nc; kk++ ) reaction_k[kk]=reaction[kk]; */

      #pragma omp  parallel                           \
        shared(reaction, reaction_k, velocity_k, q_k, \
               thread_problems, thread_solver_options,\
               thread_index, thread_index_size,       \
               error_delta_reaction, thread_iter_total)
      {
        int id = omp_get_thread_num();
        int thread_iter=0;
        int thread_info =1;
        for (unsigned int i = 0; i < thread_index_size[id]; i++ )
        {
          int contact = thread_index[id][i];
          fc3d_nsgs_index_computeqLocal(problem, reaction, contact,
                                        thread_out_index[id], thread_out_index_size[id],
                                        &(q_k[3*contact]) );
        }
        /* for (unsigned int i =0 ; i < 3*nc; i++) q_k[i] = problem->q[i]; */
        /* double normq_k = cblas_dnrm2(nc*3 , q_k , 1); */
        /* printf("############ normq = %e\t normq_k = %e\n", normq, normq_k); */

        /* call nsgs_index  with the right  fc3d_nsgs_computeqLocal*/
        fc3d_nsgs_index(thread_problems[id],
                        reaction_k, velocity_k,
                        &thread_info, thread_solver_options[id],
                        thread_index[id], thread_index_size[id]);

        /* call nsgs_index  with the right  fc3d_nsgs_computeqLocal*/
        /* fc3d_nsgs(thread_problems[id], */
        /*           reaction_k, velocity_k, */
        /*           &thread_info, thread_solver_options[id]); */
        thread_iter = thread_solver_options[id]->iparam[7];
        #pragma omp critical
        thread_iter_total += thread_iter;
        error_delta_reaction +=  thread_solver_options[id]->dparam[1];;
      }

      /* normreaction_k = cblas_dnrm2(nc*3 , reaction_k , 1); */
      /* printf("################### normreaction_k = %e \n", normreaction_k); */

      if (verbose > 0) printf("----------------------------------- FC3D - NSGS DDM NAIVE - End of thread problems after %i iterations with error_delta_reaction =%e \n", thread_iter_total, error_delta_reaction);
    } /* end of thread block */


    /* ------------------------------------------------------- */
    /* ----------------- interface loop ---------------------- */
    /* ------------------------------------------------------- */
    {
      double error_delta_reaction_interface=0.0;
      for (unsigned int i = 0; i < interface_index_size; i++ )
      {
        int  contact = interface_index[i];
        q_k[3*contact]=0.0;
        q_k[3*contact+1]=0.0;
        q_k[3*contact+2]=0.0;
        fc3d_nsgs_index_computeqLocal(problem, reaction_k, contact,
                                      interface_out_index, interface_out_index_size,
                                      &(q_k[3*contact]));
      }
      /* double normq_k = cblas_dnrm2(nc*3 , q_k , 1); */
      /* printf("normq = %e\t normq_k = %e\n", normq, normq_k); */

      int interface_hasNotConverged =1;
      int interface_iter=0;

      while ((interface_iter < interface_itermax) && ( interface_hasNotConverged > 0 ))
      {
        ++interface_iter;
        error_delta_reaction_interface=0.0;
        #pragma omp parallel for reduction(+:error_delta_reaction_interface)
        for ( unsigned int i = 0 ; i < interface_index_size ; i++)
        {
          unsigned int tid = omp_get_thread_num();
          int contact = interface_index[i];

          if (verbose > 1) printf("----------------------------------- Contact Number %i\n", contact);

          (*update_thread_problem)(contact, interface_problem, interface_local_problems[tid],
                                   reaction_k,
                                   interface_local_solver_options[tid],
                                   interface_index, interface_index_size);

          interface_local_solver_options[tid]->iparam[4] = contact;

          /* version without localreaction */
          double localreaction[3];
          {
            localreaction[0] = reaction_k[3 * contact+0];
            localreaction[1] = reaction_k[3 * contact+1];
            localreaction[2] = reaction_k[3 * contact+2];
          };

          (*local_solver)(interface_local_problems[tid], localreaction,
                          interface_local_solver_options[tid]);

          /* for( int kkk=0; kkk <3 ; kkk++) printf("localreaction[%i] = %4.2e\n", kkk,localreaction[kkk] ); */
          {
            error_delta_reaction_interface += pow(reaction_k[3 * contact] - localreaction[0], 2) +
              pow(reaction_k[3 * contact + 1] - localreaction[1], 2) +
              pow(reaction_k[3 * contact + 2] - localreaction[2], 2);

            reaction_k[3 * contact+0] = localreaction[0];
            reaction_k[3 * contact+1] = localreaction[1];
            reaction_k[3 * contact+2] = localreaction[2];

          }
        }

        /* if (error_delta_reaction_interface < tolerance) interface_hasNotConverged = 0;*/
        error_delta_reaction_interface = sqrt(error_delta_reaction_interface);
        double norm_r = cblas_dnrm2(nc*3 , reaction , 1);
        if (fabs(norm_r) > DBL_EPSILON)
          error_delta_reaction_interface /= norm_r;

        if (error_delta_reaction_interface < tolerance)
        {
          interface_hasNotConverged = 0;
          if (verbose > 0)
            printf("----------------------------------- FC3D - NSGS DDM NAIVE interface - Iteration %i Residual = %14.7e < %7.3e\n", interface_iter, error_delta_reaction_interface, options->dparam[0]);
        }
        else
        {
          if (verbose > 0)
            printf("----------------------------------- FC3D - NSGS DDM NAIVE interface - Iteration %i Residual = %14.7e > %7.3e\n", interface_iter, error_delta_reaction_interface, options->dparam[0]);
        }

      }

      error_delta_reaction+= error_delta_reaction_interface;

      if (verbose > 0)  printf("----------------------------------- FC3D - NSGS DDM NAIVE iter = %i,  error_delta_reaction = %14.7e \n", iter, error_delta_reaction);
      if (verbose > 0)  printf("----------------------------------- FC3D - NSGS DDM NAIVE iter = %i,  error_delta_reaction_interface = %14.7e \n", iter, error_delta_reaction_interface);
      //double normq = cblas_dnrm2(nc*3 , problem->q , 1);

      error_delta_reaction = sqrt(error_delta_reaction);
      double norm_r = cblas_dnrm2(nc*3 , reaction , 1);

      if (fabs(norm_r) > DBL_EPSILON)
        error_delta_reaction /= norm_r;


    }
    /* ----------------- end of interface loop --------------- */


    for (unsigned int i =0 ; i < 3*nc; i++)      reaction[i] = reaction_k[i];


    error = 0.0;
    double normq = cblas_dnrm2(nc*3 , problem->q , 1);
    (*computeError)(problem, reaction , velocity, tolerance, options, normq,  &error);
    if (error < tolerance)
    {
      hasNotConverged = 0;
      if (verbose > 0)
      {
        printf("----------------------------------- FC3D - NSGS DDM NAIVE - Iteration %i Residual = %14.7e < %7.3e\n", iter, error, options->dparam[0]);
        double normq = cblas_dnrm2(nc*3 , problem->q , 1);
        (*computeError)(problem, reaction , velocity, tolerance, options, normq,  &error_nat);
        printf("----------------------------------- FC3D - NSGS DDM NAIVE - Iteration %i Full Residual = %14.7e \n", iter, error_nat);
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
    (*freeSolver)(problem,interface_local_problems[i],interface_local_solver_options[i]);
    if (problem->M->storageType == NM_DENSE && interface_local_problems[i]->M->matrix0)
    {
      free(interface_local_problems[i]->M->matrix0);
    }
    interface_local_problems[i]->M->matrix0 = NULL;
    freeFrictionContactProblem(interface_local_problems[i]);
    solver_options_delete(interface_local_solver_options[i]);
    free(interface_local_solver_options[i]);
  }


}
