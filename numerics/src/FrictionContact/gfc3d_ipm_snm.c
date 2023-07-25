/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#include "CSparseMatrix_internal.h"
#include "gfc3d_Solvers.h"
#include "gfc3d_compute_error.h"
#include "SiconosLapack.h"

#include "SparseBlockMatrix.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "numerics_verbose.h"
#include "NumericsVector.h"
#include "float.h"
#include "JordanAlgebra.h"

#include "NumericsSparseMatrix.h"
#include "NumericsMatrix.h"

#include "projectionOnCone.h"

/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "siconos_debug.h"
#include "gfc3d_ipm.h"

const char* const   SICONOS_GLOBAL_FRICTION_3D_IPM_SEMISMOOTH_STR = "GFC3D IPM SNM";

/* ------------------------- Helper functions implementation ------------------------------ */

/* Returns the maximum step-length to the boundary reduced by a factor gamma. Uses long double. */
double *array_getStepLength(const double * const x, const double * const dx, const unsigned int vecSize,
                     const unsigned int varsCount, const double gamma)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  unsigned int pos;
  float_type aL, bL, cL, dL, alphaL;
  double *arr_alpha = (double*)calloc(varsCount,sizeof(double));

  double alpha = 1e20; //1.0;

  for(unsigned int i = 0; i < varsCount; ++i)
  {
    pos = i * dimension;
    aL = dnrm2l(dimension-1, dx+pos+1);
    aL = (dx[pos] - aL)*(dx[pos] + aL);
    bL = x[pos]*dx[pos];
    for (int k = 1; k < dimension; bL -= x[pos+k]*dx[pos+k], k++);
    cL = dnrm2l(dimension-1, x+pos+1);
    cL = (x[pos] - cL)*(x[pos] + cL);
    dL = bL*bL - aL*cL;
    if(aL < 0 || (bL < 0 && dL > 0 ))
      if (bL>0)
        alphaL = -(bL+sqrtl(dL))/aL;
      else
        alphaL = cL/(-bL+sqrtl(dL));
    else if((fabsl(aL) == 0.0) && (bL < 0))
      alphaL = -cL/bL/2;
    else
      alphaL = DBL_MAX;
    // printf("Cone %i: alpha = %Le\n", i, alphaL);
    // alpha = ((alphaL < alpha) ? alphaL : alpha);
    // alpha = gamma*alpha;
    // alpha = ((alpha < 1.0) ? alpha : 1.0);

    arr_alpha[i] = alphaL;
  }

  return arr_alpha;
}

/* Returns the primal constraint vector for global fricprob: out = velocity - H x globalVelocity - w - phi(s)*/
/* and the relative 2-norm of this vector: |out|/max{|velocity|, |H x globalVelocity|, |w|, |phi(s)|} */
void primalResidual_s(const double * velocity, NumericsMatrix * H, const double * globalVelocity, const double * w,
                    const double * s, double * out, double * rnorm, const double tol)
{
  size_t nd = H->size0;
  double rn;


  /* The memory for the result vectors should be allocated using calloc
   * since H is a sparse matrix. In other case the behaviour will be undefined.*/
  //  double *Hv = (double*)calloc(nd, sizeof(double));
  //double *u_minus_Hv = (double*)calloc(nd, sizeof(double));

  NM_gemv(-1.0, H, globalVelocity, 0.0, out);
  rn = cblas_dnrm2(nd, out, 1);
  cblas_daxpy(nd, 1.0, velocity, 1, out, 1);
  cblas_daxpy(nd, -1.0, w, 1, out, 1);

  for(unsigned int i=0; i<nd; i+=3) out[i] -= s[i/3];

  rn = fmax(rn, cblas_dnrm2(nd, velocity, 1));
  rn = fmax(rn, cblas_dnrm2(nd, w, 1));
  rn = fmax(rn, cblas_dnrm2(nd/3, s, 1));
  *rnorm = (rn > tol ? cblas_dnrm2(nd, out, 1) : cblas_dnrm2(nd, out, 1));

  /* *rnorm = cblas_dnrm2(nd, out, 1);  */
  // printf("rn = %e, tol = %e\n", rn, tol);
}


/* Computation of the projection error |r - proj(r-u)|/max{|r|, |u|} */
static double projectionError(const double * velocity, const double * reaction, const unsigned int nc, const double tol)
{
   double worktmp[3];
   double out = 0.0;
   double norm_u, norm_r, relative_scaling;

   for(int ic = 0 ; ic < nc ; ic++)
     {
       worktmp[0] = reaction[3*ic] -  velocity[3*ic] ;
       worktmp[1] = reaction[3*ic+1] -  velocity[3*ic+1] ;
       worktmp[2] = reaction[3*ic+2] -  velocity[3*ic+2] ;
       projectionOnCone(worktmp, 1.0);
       worktmp[0] = reaction[3*ic] -  worktmp[0];
       worktmp[1] = reaction[3*ic+1] -  worktmp[1];
       worktmp[2] = reaction[3*ic+2] -  worktmp[2];
       out +=  worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2];
     }
   out = sqrt(out);
   norm_u = cblas_dnrm2(3*nc, velocity, 1);
   norm_r = cblas_dnrm2(3*nc, reaction, 1);
   relative_scaling = fmax(norm_u, norm_r);
   if(relative_scaling > tol)
     out = out/relative_scaling;
   return out;
}


/* Writing problem data under a Matlab format in a file  */
/* The data are printed under the form of a dense format */
/* problem: min .5 v'*M*v + f'*v, s.t. H*v + w \in K (Lorentz cone)
   d = dimenion
   n = number of contact points
   m = number of degrees of freedom
   M = m x m matrix
   f = m-vector
   H = n*d x m matrix
   w = n*d-vector */
static void printDataProbMatlabFile(NumericsMatrix * M, double * f, NumericsMatrix * H, double * w, int d, int n, int m, double * mu, FILE * file)
{
  fprintf(file,"d = %3i;\n",d);
  fprintf(file,"n = %6i;\n",n);
  fprintf(file,"m = %6i;\n",m);

  fprintf(file,"M = [\n");
  CSparseMatrix_print_in_Matlab_file(NM_triplet(M), 0, file);
  fprintf(file,"];\n");
  fprintf(file,"M = sparse(int32(M(:,1)), int32(M(:,2)), M(:,3));\n");

  fprintf(file,"H = [\n");
  CSparseMatrix_print_in_Matlab_file(NM_triplet(H), 0, file);
  fprintf(file,"];\n");
  fprintf(file,"H = sparse(int32(H(:,1)), int32(H(:,2)), H(:,3));\n");

  fprintf(file,"f = [");
  for(int i = 0; i < m; i++)
  {
    fprintf(file,"%8.20e; ",f[i]);
  }
  fprintf(file,"];\n");

  fprintf(file,"w = [");
  for(int i = 0; i < n*d; i++)
  {
    fprintf(file,"%8.20e; ",w[i]);
  }
  fprintf(file,"];\n");

  fprintf(file,"mu = [");
  for(int i = 0; i < n; i++)
  {
    fprintf(file,"%8.20e; ",mu[i]);
  }
  fprintf(file,"];\n");
}

/* print iteres under a Matlab format in a file */
/* iteration = index of the iteration
   v = global velocity
   u = velocity
   r = reaction
   d = dimenion
   n = number of contact points
   m = number of degrees of freedom
*/

static void printIteresProbMatlabFile(int iteration, double * v, double * u, double * r, double * s, double * dv, double * du, double * dr, double * ds, int d, int n, int m, FILE * file)
// static void printIteresProbMatlabFile(int iteration, double pinfeas, double dinfeas, double udotr, double smub, int d, int n, int m, FILE * file)
{

  fprintf(file,"v(%3i,:) = [",iteration+1);
  for(int i = 0; i < m; i++)
  {
    fprintf(file, "%8.20e, ", v[i]);
  }
  fprintf(file,"];\n");

  fprintf(file,"u(%3i,:) = [",iteration+1);
  for(int i = 0; i < n*d; i++)
  {
    fprintf(file, "%8.20e, ", u[i]);
  }
  fprintf(file,"];\n");

  fprintf(file,"r(%3i,:) = [",iteration+1);
  for(int i = 0; i < n*d; i++)
  {
    fprintf(file, "%8.20e, ", r[i]);
  }
  fprintf(file,"];\n");

  fprintf(file,"s(%3i,:) = [",iteration+1);
  for(int i = 0; i < n; i++)
  {
    fprintf(file, "%8.20e, ", s[i]);
  }
  fprintf(file,"];\n");

  fprintf(file,"dv(%3i,:) = [",iteration+1);
  for(int i = 0; i < m; i++)
  {
    fprintf(file, "%8.20e, ", dv[i]);
  }
  fprintf(file,"];\n");

  fprintf(file,"du(%3i,:) = [",iteration+1);
  for(int i = 0; i < n*d; i++)
  {
    fprintf(file, "%8.20e, ", du[i]);
  }
  fprintf(file,"];\n");

  fprintf(file,"dr(%3i,:) = [",iteration+1);
  for(int i = 0; i < n*d; i++)
  {
    fprintf(file, "%8.20e, ", dr[i]);
  }
  fprintf(file,"];\n");

  fprintf(file,"ds(%3i,:) = [",iteration+1);
  for(int i = 0; i < n; i++)
  {
    fprintf(file, "%8.20e, ", ds[i]);
  }
  fprintf(file,"];\n");

  // fprintf(file,"pinfeas(%3i) = %20.16e;\n",iteration+1,pinfeas);
  // fprintf(file,"dinfeas(%3i) = %20.16e;\n",iteration+1,dinfeas);
  // fprintf(file,"udotr(%3i) = %20.16e;\n",iteration+1,udotr);
  // // fprintf(file,"residu(%3i) = %20.16e;\n",iteration+1,totalresidual);
  // fprintf(file,"smub(%3i) = %20.16e;\n",iteration+1,smub);



  return;
}

static void printVectorMatlabFile(int iteration, double * vec, int vecSize, FILE * file)
{
  fprintf(file,"vector(%4i,:) = [",iteration+1);
  for(int i = 0; i < vecSize; i++)
  {
    fprintf(file, "%24.16e, ", vec[i]);
  }
  fprintf(file,"];\n");
  return;
}

// cl: classify
static void printBlockVec(double * vec, int vecSize, int sizeBlock, int cl)
{
  if (cl > 0)
  {
    for(unsigned int i=0; i<vecSize; i++)
    {
      if (i%sizeBlock == 0) printf("\n%i-%i: ", i, i+sizeBlock-1);
      printf(" %.8e,", vec[i]);
      if ((i+1)%sizeBlock == 0) printf("\t\t x0 - |xb| = %.2e", vec[((i+1)/sizeBlock-1)*sizeBlock] - cblas_dnrm2(sizeBlock-1, vec+((i+1)/sizeBlock-1)*sizeBlock+1, 1));
    }
  }

  else
  {
    for(unsigned int i=0; i<vecSize; i++)
    {
      if (i%sizeBlock == 0) printf("\n%i-%i: ", i, i+sizeBlock-1);
      printf(" %.15e", vec[i]);
    }
  }
  return;
}



static CS_INT cs_print_2 (const cs *A, CS_INT brief)
{
    CS_INT p, j, m, n, nzmax, nz, *Ap, *Ai ;
    CS_ENTRY *Ax ;
    if (!A) { printf ("(null)\n") ; return (0) ; }
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    nzmax = A->nzmax ; nz = A->nz ;
    printf ("CXSparse Version %d.%d.%d, %s.  %s\n", CS_VER, CS_SUBVER,
        CS_SUBSUB, CS_DATE, CS_COPYRIGHT) ;
    if (nz < 0)
    {
        printf ("%g-by-%g, nzmax: %g nnz: %g, 1-norm: %g\n", (double) m,
            (double) n, (double) nzmax, (double) (Ap [n]), cs_norm (A)) ;
        for (j = 0 ; j < n ; j++)
        {
            printf ("    col %g : locations %g to %g\n", (double) j,
                (double) (Ap [j]), (double) (Ap [j+1]-1)) ;
            for (p = Ap [j] ; p < Ap [j+1] ; p++)
            {
                printf ("      %g : ", (double) (Ai [p])) ;
#ifdef CS_COMPLEX
                printf ("(%g, %g)\n",
                    Ax ? CS_REAL (Ax [p]) : 1, Ax ? CS_IMAG (Ax [p]) : 0) ;
#else
                printf ("%.30g\n", Ax ? Ax [p] : 1) ;
#endif
                if (brief && p > 20) { printf ("  ...\n") ; return (1) ; }
            }
        }
    }
    else
    {
        printf ("triplet: %g-by-%g, nzmax: %g nnz: %g\n", (double) m,
            (double) n, (double) nzmax, (double) nz) ;
        for (p = 0 ; p < nz ; p++)
        {

            printf ("    %g %g : ", (double) (Ai [p]), (double) (Ap [p])) ;
#ifdef CS_COMPLEX
            printf ("(%g, %g)\n",
                Ax ? CS_REAL (Ax [p]) : 1, Ax ? CS_IMAG (Ax [p]) : 0) ;
#else
            printf ("%g\n", Ax ? Ax [p] : 1) ;
#endif
            if (brief && p > 20) { printf ("  ...\n") ; return (1) ; }
        }
    }
    return (1) ;
}



// #include "NM_conversions.h"
// static CS_INT cs_lusol_2 (CS_INT order, const cs *A, CS_ENTRY *b, double tol, FILE * file)
// {
//     CS_ENTRY *x ;
//     css *S ;
//     csn *N ;
//     CS_INT n, ok ;
//     if (!CS_CSC (A) || !b) return (0) ;     /* check inputs */
//     n = A->n ;
//     S = cs_sqr (order, A, 0) ;              /* ordering and symbolic analysis */
//     N = cs_lu (A, S, tol) ;                 /* numeric LU factorization */
//     x = cs_malloc (n, sizeof (CS_ENTRY)) ;    /* get workspace */
//     ok = (S && N && x) ;
//     if (ok)
//     {
//         cs_ipvec (N->pinv, b, x, n) ;       /* x = b(p) */
//         cs_lsolve (N->L, x) ;               /* x = L\x */
//         cs_usolve (N->U, x) ;               /* x = U\x */
//         cs_ipvec (S->q, x, b, n) ;          /* b(q) = x */
//     }

//     // printf("\nP = "); for(int i=0; i<n; i++) printf("%ld, ", S->q[i]);
//     fprintf(file,"P_C = [");
//     for(int i = 0; i < n; i++)
//     {
//       fprintf(file,"%8.20ld ",N->pinv[i]);
//     }
//     fprintf(file,"];\n");

//     fprintf(file,"S_C = [");
//     for(int i = 0; i < n; i++)
//     {
//       fprintf(file,"%8.20ld ",S->q[i]);
//     }
//     fprintf(file,"];\n");


//     // printf("\nL = "); cs_print_2(N->L, 0);  printf("\n");
//     fprintf(file,"L_C = [\n");
//     CSparseMatrix_print_in_Matlab_file(NM_csc_to_triplet(N->L), 0, file);
//     fprintf(file,"];\n");
//     fprintf(file,"L_C = sparse(int32(L_C(:,1)), int32(L_C(:,2)), L_C(:,3));\n");

//     // printf("\nU = "); cs_print_2(N->U, 0); printf("\n");
//     fprintf(file,"U_C = [\n");
//     CSparseMatrix_print_in_Matlab_file(NM_csc_to_triplet(N->U), 0, file);
//     fprintf(file,"];\n");
//     fprintf(file,"U_C = sparse(int32(U_C(:,1)), int32(U_C(:,2)), U_C(:,3));\n");
//     cs_free (x) ;
//     cs_sfree (S) ;
//     cs_nfree (N) ;
//     return (ok) ;
// }








/* static int saveMatrix(NumericsMatrix* m, const char * filename) */
/* { */
/*     NumericsMatrix * md = NM_create(NM_DENSE, m->size0, m->size1); */
/*     NM_to_dense(m, md); */
/*     FILE *f; */
/*     f = fopen(filename, "wb"); */
/*     if (!f) */
/*         return 1; */
/*     for (int i = 0; i < m->size0; ++i) */
/*         for (int j = 0; j < m->size1; ++j) */
/*             fwrite(&(md->matrix0[i+j*md->size0]), sizeof(double), 1, f); */
/*     fclose(f); */
/*     NM_free(md); */
/*     return 0; */
/* } */


/* static int saveVector(double * vec, const unsigned int vecSize, const char * filename) */
/* { */
/*     FILE *f; */
/*     f = fopen(filename, "wb"); */
/*     if (!f) */
/*         return 1; */
/*     for (unsigned int i = 0; i < vecSize; ++i) */
/*         fwrite(&(vec[i]), sizeof(double), 1, f); */
/*     fclose(f); */
/*     return 0; */
/* } */


static float randomFloat(float min, float max) {
    return min + (float)rand() / ((float)RAND_MAX / (max - min));
}


/* --------------------------- Interior-point method implementation ------------------------------ */
/*
 * Implementation contains the following functions:
 *  - gfc3d_IPM_SNM_init - initialize solver (allocate memory)
 *  - gfc3d_IPM_SNM_free - deallocate memory
 *  - gfc3d_IPM_SNM_setDefaultSolverOptions - setup default solver parameters
 *  - gfc3d_IPM_SNM - optimization method
 */
void gfc3d_IPM_SNM_init(GlobalFrictionContactProblem* problem, SolverOptions* options)
{
  unsigned int m = problem->M->size0;
  unsigned int nd = problem->H->size1;
  unsigned int d = problem->dimension;

  if(!options->dWork || options->dWorkSize != (size_t)(m + nd + nd + nd/d))
  {
    options->dWork = (double*)calloc(m + nd + nd + nd/d, sizeof(double));
    options->dWorkSize = m + nd + nd + nd/d;
  }


  /* ------------- initialize starting point ------------- */
  options->solverData=(Gfc3d_IPM_init_data *)malloc(sizeof(Gfc3d_IPM_init_data));
  Gfc3d_IPM_init_data * data = (Gfc3d_IPM_init_data *)options->solverData;

  /* --------- allocate memory for tmp point ----------- */
  data->tmp_point = (IPM_tmp_point*)malloc(sizeof(IPM_tmp_point));
  data->tmp_point->t_globalVelocity = (double*)calloc(m, sizeof(double));
  data->tmp_point->t_velocity = (double*)calloc(nd, sizeof(double));
  data->tmp_point->t_reaction = (double*)calloc(nd, sizeof(double));

  /* 1. v */
  data->starting_point = (IPM_starting_point*)malloc(sizeof(IPM_starting_point));
  data->starting_point->globalVelocity = (double*)calloc(m, sizeof(double));
  for(unsigned int i = 0; i < m; ++ i)
    data->starting_point->globalVelocity[i] = 0.01;

  /* 2. u */
  data->starting_point->velocity = (double*)calloc(nd, sizeof(double));
  for(unsigned int i = 0; i < nd; ++ i)
  {
    data->starting_point->velocity[i] = 0.01;
    if(i % d == 0)
      data->starting_point->velocity[i] = 0.1;
  }

  /* 3. r */
  data->starting_point->reaction = (double*)calloc(nd, sizeof(double));
  for(unsigned int i = 0; i < nd; ++ i)
  {
    data->starting_point->reaction[i] = 0.01; //0.0351;
    if(i % d == 0)
      data->starting_point->reaction[i] = 0.1; //0.2056;
  }

  /* ------ initialize the change of variable matrix P_mu ------- */
  data->P_mu = (IPM_change_of_variable*)malloc(sizeof(IPM_change_of_variable));
  data->P_mu->mat = NM_create(NM_SPARSE, nd, nd);
  NM_triplet_alloc(data->P_mu->mat, nd);
  data->P_mu->mat->matrix2->origin = NSM_TRIPLET;
  for(unsigned int i = 0; i < nd; ++i)
    if(i % d == 0)
      /* NM_entry(data->P_mu->mat, i, i, 1. / problem->mu[(int)(i/d)]); */
      NM_entry(data->P_mu->mat, i, i, 1.);
    else
      /* NM_entry(data->P_mu->mat, i, i, 1.); */
      NM_entry(data->P_mu->mat, i, i, problem->mu[(int)(i/d)]);

  /* ------ initialize the inverse P_mu_inv of the change of variable matrix P_mu ------- */
  data->P_mu->inv_mat = NM_create(NM_SPARSE, nd, nd);
  NM_triplet_alloc(data->P_mu->inv_mat, nd);
  data->P_mu->inv_mat->matrix2->origin = NSM_TRIPLET;
  for(unsigned int i = 0; i < nd; ++i)
    if(i % d == 0)
      /* NM_entry(data->P_mu->inv_mat, i, i, problem->mu[(int)(i/d)]); */
      NM_entry(data->P_mu->inv_mat, i, i, 1.);
    else
      /* NM_entry(data->P_mu->inv_mat, i, i, 1.); */
      NM_entry(data->P_mu->inv_mat, i, i, 1.0/problem->mu[(int)(i/d)]);
  /* ------ initial parameters initialization ---------- */
  data->internal_params = (IPM_internal_params*)malloc(sizeof(IPM_internal_params));
  data->internal_params->alpha_primal = 1.0;
  data->internal_params->alpha_dual = 1.0;
  data->internal_params->sigma = 0.1;
  data->internal_params->barr_param = 1.0;


  /* ----- temporary vaults initialization ------- */
  data->tmp_vault_nd = (double**)malloc(17 * sizeof(double*));
  for(unsigned int i = 0; i < 17; ++i)
    data->tmp_vault_nd[i] = (double*)calloc(nd, sizeof(double));

  data->tmp_vault_m = (double**)malloc(2 * sizeof(double*));
  for(unsigned int i = 0; i < 2; ++i)
    data->tmp_vault_m[i] = (double*)calloc(m, sizeof(double));

}


void gfc3d_IPM_SNM_free(GlobalFrictionContactProblem* problem, SolverOptions* options)
{
  if(options->dWork)
  {
    free(options->dWork);
    options->dWork=NULL;
    options->dWorkSize = 0;
  }
  if(options->solverData)
  {
    Gfc3d_IPM_init_data * data = (Gfc3d_IPM_init_data *)options->solverData;

    free(data->starting_point->globalVelocity);
    data->starting_point->globalVelocity = NULL;

    free(data->starting_point->velocity);
    data->starting_point->velocity = NULL;

    free(data->starting_point->reaction);
    data->starting_point->reaction = NULL;

    free(data->starting_point);

    NM_clear(data->P_mu->mat);
    free(data->P_mu->mat);
    data->P_mu->mat = NULL;

    NM_clear(data->P_mu->inv_mat);
    free(data->P_mu->inv_mat);
    data->P_mu->inv_mat = NULL;

    free(data->P_mu);

    for(unsigned int i = 0; i < 17; ++i)
      free(data->tmp_vault_nd[i]);
    free(data->tmp_vault_nd);
    data->tmp_vault_nd = NULL;

    for(unsigned int i = 0; i < 2; ++i)
      free(data->tmp_vault_m[i]);
    free(data->tmp_vault_m);
    data->tmp_vault_m = NULL;

    free(data->tmp_point->t_globalVelocity);
    data->tmp_point->t_globalVelocity = NULL;

    free(data->tmp_point->t_velocity);
    data->tmp_point->t_velocity = NULL;

    free(data->tmp_point->t_reaction);
    data->tmp_point->t_reaction = NULL;

    free(data->tmp_point);

    free(data->internal_params);
  }

}

void gfc3d_IPM_SNM(GlobalFrictionContactProblem* restrict problem, double* restrict reaction,
               double* restrict velocity, double* restrict globalVelocity,
               int* restrict info, SolverOptions* restrict options, const char* problem_name)
{
  // verbose = 3;
  // printf("DBL_EPSILON %25.15e\n",DBL_EPSILON);
  // the size of the problem detection
  unsigned int m = problem->M->size0;
  unsigned int nd = problem->H->size1;
  unsigned int d = problem->dimension;
  unsigned int n = problem->numberOfContacts;

  NumericsMatrix* M = NULL;
  NumericsMatrix* H_tilde = NULL;

  /* globalFrictionContact_display(problem); */

  /* symmetrization of the matrix M */
  if(!(NM_is_symmetric(problem->M)))
  {
    printf("#################### SYMMETRIZATION ####################\n");
    NumericsMatrix *MT = NM_transpose(problem->M);
    NumericsMatrix * MSym = NM_add(1 / 2., problem->M, 1 / 2., MT);
    NM_free(problem->M);
    problem->M = MSym;
    NM_free(MT);
  }



  //for(int i = 0; i < n ; i++) printf("mu[%d] = %g\n", i, problem->mu[i]);
  //for(int i = 0; i < n ; i++) problem->mu[i]=0.75;

  /* if SICONOS_FRICTION_3D_IPM_FORCED_SPARSE_STORAGE = SICONOS_FRICTION_3D_IPM_FORCED_SPARSE_STORAGE,
     we force the copy into a NM_SPARSE storageType */
  DEBUG_PRINTF("problem->M->storageType : %i\n",problem->M->storageType);
  if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_SPARSE_STORAGE] == SICONOS_FRICTION_3D_IPM_FORCED_SPARSE_STORAGE
     && problem->M->storageType == NM_SPARSE_BLOCK)
  {
    DEBUG_PRINT("Force a copy to sparse storage type\n");
    printf("\n\n\n######################### FORCE SPARSE STORAGE #########################\n\n\n");
    M = NM_create(NM_SPARSE,  problem->M->size0,  problem->M->size1);
    NM_copy_to_sparse(problem->M, M, DBL_EPSILON);
  }
  else
  {
    M = problem->M;
  }

  NumericsMatrix * Minv = NULL;
  if ( options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_QPH )
  {
    int block_number_of_M = M->size0/3;
    unsigned int * blocksizes_of_M = (unsigned int*)malloc(block_number_of_M * sizeof(unsigned int));
    for (int i = 0; i < block_number_of_M; i++)
      *(blocksizes_of_M + i) = 3;
    Minv = NM_inverse_diagonal_block_matrix(M, block_number_of_M, blocksizes_of_M);
    free(blocksizes_of_M);
  }



  DEBUG_PRINTF("problem->M->storageType : %i\n",problem->H->storageType);
  if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_SPARSE_STORAGE] == SICONOS_FRICTION_3D_IPM_FORCED_SPARSE_STORAGE
     && problem->H->storageType == NM_SPARSE_BLOCK)
  {
    DEBUG_PRINT("Force a copy to sparse storage type\n");
    H_tilde = NM_create(NM_SPARSE,  problem->H->size1,  problem->H->size0);
    NM_copy_to_sparse(NM_transpose(problem->H), H_tilde, DBL_EPSILON);
  }
  else
  {
    H_tilde = NM_transpose(problem->H);
  }

  printf("nnz H = %8zu density = %9.4f\n",NM_nnz(problem->H), NM_nnz(problem->H)/1.0/nd/m);

  // initialize solver if it is not set
  int internal_allocation=0;
  if(!options->dWork || (options->dWorkSize != (size_t)(m + nd + nd)))
  {
    gfc3d_IPM_SNM_init(problem, options);
    internal_allocation = 1;
  }

  Gfc3d_IPM_init_data * data = (Gfc3d_IPM_init_data *)options->solverData;
  NumericsMatrix *P_mu = data->P_mu->mat;
  NumericsMatrix *P_mu_inv = data->P_mu->inv_mat;

  double *w_tilde = problem->b;
  double *w = data->tmp_vault_nd[0];
  double *f = problem->q;

  double *iden;

  // change of variable to eliminate the friction coefficients: H_tilde --> H and w_tilde --> w
  NumericsMatrix *H = NM_multiply(P_mu, H_tilde);
  NM_gemv(1.0, P_mu, w_tilde, 0.0, w);

  double *w_ori = (double*)calloc(nd, sizeof(double));
  NV_copy(w, nd, w_ori);

  // compute -H
  NumericsMatrix *minus_H = NM_create(H->storageType, H->size0, H->size1);
  NM_copy(H, minus_H);
  NM_scal(-1.0, minus_H);
  NumericsMatrix * minus_Ht = NM_transpose(minus_H);

  // compute H'
  NumericsMatrix * Ht = NM_transpose(H);


  double alpha_primal = data->internal_params->alpha_primal;
  double alpha_dual = data->internal_params->alpha_dual;
  double barr_param = data->internal_params->barr_param;
  double sigma = data->internal_params->sigma;
  sigma = 0.6;

  cblas_dcopy(nd, data->starting_point->reaction, 1, reaction, 1);
  cblas_dcopy(nd, data->starting_point->velocity, 1, velocity, 1);
  cblas_dcopy(m, data->starting_point->globalVelocity, 1, globalVelocity, 1);

  NumericsMatrix *minus_M = NM_create(
              M->storageType, M->size0,
              M->size1);  // store the matrix -M to build the matrix of the Newton linear system
  /* Create the matrix -M to build the matrix of the reduced linear system */
  NM_copy(M, minus_M);
  NM_scal(-1.0, minus_M);


  /* COMPUTATION OF A NEW STARTING POINT */

  // set the reaction vector to an arbitrary value in the interior of the cone
  for (unsigned int  i = 0; i<nd; i++)
      if (i % d == 0) reaction[i] = 0.1;
      else reaction[i] = 0.01;

  // computation of the global velocity vector: v = M\(H'*r+f)
  for (unsigned int  i = 0; i<m; i++) globalVelocity[i] = f[i];
  NM_tgemv(1.0, H, reaction, 1.0, globalVelocity);
  NM_Cholesky_solve(NM_preserve(M), globalVelocity, 1);

  // computation of the velocity u = proj(H*v + w), then move u in the interior of the cone
  /* for (unsigned int  i = 0; i<nd; i++) */
  /*   velocity[i] = w[i]; */
  /* NM_gemv(1.0, H, globalVelocity, 1.0, velocity); */
  /* for (unsigned int i = 0; i < n; i++) */
  /*   { */
  /*     projectionOnCone(velocity+3*i, 0.1); */
  /*     if (velocity[3*i] == 0) */
  /*     { */
  /*       velocity[3*i] = 0.1; */
  /*       velocity[3*i+1] = 0.01; */
  /*       velocity[3*i+2] = 0.01; */
  /*     } */
  /*   } */

  for (unsigned int  i = 0; i<nd; i++)
      if (i % d == 0) velocity[i] = 0.1;
      else velocity[i] = 0.01;


  double tol = options->dparam[SICONOS_DPARAM_TOL];
  unsigned int max_iter = options->iparam[SICONOS_IPARAM_MAX_ITER];

  double sgmp1 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1];
  double sgmp2 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2];
  double sgmp3 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3];
  double gmmp1 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1];
  double gmmp2 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2];

  int hasNotConverged = 1;
  unsigned int iteration = 0;
  double pinfeas = 1e300; double pinfeas_new = 1e300;
  double dinfeas = 1e300;
  double complem = 1e300;
  double complem_p = 1e300;
  double dualgap = 1e300;
  double udotr = 1e300, udotr_mu = 1e300; double *udotr_mu_vec = (double*)calloc(nd,sizeof(double));
  double projerr = 1e300;
  double error[6];
  double totalresidual = 1e300, totalresidual_mu = 1e300;

  double diff_fixp = 1e300;

  double *primalConstraint = data->tmp_vault_nd[1];
  double *dualConstraint = data->tmp_vault_m[0];
  double *complemConstraint = data->tmp_vault_nd[2];
  double *fixpConstraint = (double*)calloc(n,sizeof(double));

  double gmm = gmmp1+gmmp2;
  double barr_param_a, e;
  double norm_f = cblas_dnrm2(m, f, 1);
  double norm_w = cblas_dnrm2(nd, w, 1);

  double *velocity_inv = (double*)calloc(nd,sizeof(double));
  double *d_globalVelocity = (double*)calloc(m,sizeof(double));
  double *d_velocity = (double*)calloc(nd,sizeof(double));
  double *d_reaction = (double*)calloc(nd,sizeof(double));
  double *dw = (double*)calloc(nd,sizeof(double));
  double *d_s = (double*)calloc(n,sizeof(double));
  double *s = (double*)calloc(n,sizeof(double));
  for (unsigned int  i = 0; i<n; i++) s[i] = 0.014;

  double *rhs = options->dWork;
  double *rhs_2 = (double*)calloc(m+2*nd+n, sizeof(double));
  double *sol = (double*)calloc(m+2*nd+n, sizeof(double));
  double *gv_plus_dgv = data->tmp_vault_m[1];
  double *vr_jprod = data->tmp_vault_nd[3];
  double *v_plus_dv = data->tmp_vault_nd[4];
  double *r_plus_dr = data->tmp_vault_nd[5];
  double *vr_prod_sub_iden = data->tmp_vault_nd[6];
  double *dvdr_jprod = data->tmp_vault_nd[7];

  //double * r_p = (double*)calloc(nd,sizeof(double));                          // scaling vector p
  NumericsMatrix* r_Qp = NULL;                                                // matrix Qp
  double * r_rhs = (double*)calloc(m+nd, sizeof(double));
  double * r_rhs_2 = (double*)calloc(m+nd, sizeof(double));
  double *sr_rhs = (double*)calloc(nd,sizeof(double));
  double *sr_rhs_2 = (double*)calloc(nd,sizeof(double));
  double * r_dv = (double*)calloc(m,sizeof(double));
  double * r_dr = (double*)calloc(nd,sizeof(double));
  double * r_du = (double*)calloc(nd,sizeof(double));
  double * r_dv_a = (double*)calloc(m,sizeof(double));
  double * r_dr_a = (double*)calloc(nd,sizeof(double));
  double * r_du_a = (double*)calloc(nd,sizeof(double));
  double * r_adu = (double*)calloc(nd, sizeof(double));
  double * r_adr = (double*)calloc(nd, sizeof(double));
  double r_alpha_p, r_alpha_d; /* primal and dual steplengths */
  double *arr_alpha_primal = NULL, *arr_alpha_dual = NULL;
  double r_mu, r_mu_a; /* duality gap, affine duality gap */
  NumericsMatrix *JR = NULL; /* Reduced Jacobian with NT scaling */
  long JR_nzmax;
  double * Hvw = (double*)calloc(nd, sizeof(double));
  double err = 1e300;
  char fws = ' '; /* finish without scaling */

  /* list of active constraints : = 0 if x_0 <= epsilon, = 1 if lambda_2 <= epsilon , = 3 either */
  short * a_velo = (short*)calloc(n, sizeof(short));
  short * a_reac = (short*)calloc(n, sizeof(short));

  /* norm of the residuals of teh second linear system */
  double LS_norm_p = 0.; // primal feasibility
  double LS_norm_d = 0.; // dual feaqsibility
  double LS_norm_c = 0.; // complementarity
  double LS_norm_f = 0.; // fixed point

  NumericsMatrix *J = NULL, *J_dense = NULL;

  long J_nzmax;
  size_t H_nzmax = NM_nnz(H);
  size_t M_nzmax = NM_nnz(M);

  if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO] == SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_YES)
  {
    numerics_printf_verbose(1,"---- GFC3D - IPM - Problem information");
    numerics_printf_verbose(1,"---- GFC3D - IPM - 1-norm of M = %g norm of f = %g ", NM_norm_1(M), norm_f);
    numerics_printf_verbose(1,"---- GFC3D - IPM - inf-norm of M = %g ", NM_norm_inf(M));

    numerics_printf_verbose(1,"---- GFC3D - IPM - 1-norm of H = %g norm of w = %g ", NM_norm_1(problem->H), norm_w);
    numerics_printf_verbose(1,"---- GFC3D - IPM - inf-norm of H = %g ", NM_norm_inf(problem->H));
    numerics_printf_verbose(1,"---- GFC3D - IPM - M is symmetric = %i ", NM_is_symmetric(M));

    numerics_printf_verbose(1,"---- GFC3D - IPM - M size = (%i, %i) ", M->size0, M->size1);
    numerics_printf_verbose(1,"---- GFC3D - IPM - H size = (%i, %i) ", problem->H->size0, problem->H->size1);
  }

  /* ---- IPM iterations ---- */
  numerics_printf_verbose(-1, "problem dimensions n, nd x m: %1i, %6i x %-6i",n, nd, m);
      switch ( options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] )
    {
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL:
    {
      numerics_printf_verbose(-1,"LS solution: 3x3 no scaling\n");
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QP2:
    {
      numerics_printf_verbose(-1,"LS solution: 3x3 NT scaling with Qp2\n");
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QPH:
    {
      numerics_printf_verbose(-1,"LS solution: 3x3 NT scaling with QpH\n");
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QP2:
    {
      numerics_printf_verbose(-1,"LS solution: 2x2 NT scaling with Qp2\n");
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QPH:
    {
      numerics_printf_verbose(-1,"LS solution: 2x2 NT scaling with QpH\n");
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_QPH:
    {
      numerics_printf_verbose(-1,"LS solution: 1x1 NT scaling with QpH\n");
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_NOSCAL:
    {
      numerics_printf_verbose(-1,"LS solution: 4x4 no scaling\n");
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_QP2:
    {
      numerics_printf_verbose(-1,"LS solution: 4x4 NT scaling with Qp2\n");
      break;
    }
    default:
    {
      printf("ERROR\n");
    }
    }

  numerics_printf_verbose(-1, "| it  | pinfeas | dinfeas |  |s-ub| | |uor-mu||2max|uor-mu||   4*mu  |  u'r/n  | prj err | barpram |  alpha  |  |dv|   |  |du|   |  |dr|   |  |ds|   | ls prim | ls dual | ls comp | ls fixP |");
  numerics_printf_verbose(-1, "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");

  double * p = data->tmp_vault_nd[8];
  double * p2 = data->tmp_vault_nd[9];
  double * pinv = data->tmp_vault_nd[10];
  NumericsMatrix* Qp = NULL;
  NumericsMatrix* Qpinv = NULL;
  NumericsMatrix* F = NULL;
  NumericsMatrix* Finv = NULL;
  NumericsMatrix* tmpmat = NULL;
  NumericsMatrix* Qp_F = NULL;
  NumericsMatrix* Qp2 = NULL;
  NumericsMatrix* F2 = NULL;
  NumericsMatrix * eye_nd = NM_eye(nd);
  NumericsMatrix * eye_n = NM_eye(n);
  NumericsMatrix * subdiff_u = NULL;

  NumericsMatrix * minus_e = NULL;
  minus_e = NM_create(NM_SPARSE, nd, n);
  size_t minus_e_nzmax = n;
  NM_triplet_alloc(minus_e, minus_e_nzmax);
  NM_fill(minus_e, NM_SPARSE, nd, n, minus_e->matrix2);
  for(size_t i = 0; i < n; ++i)
  {
    NM_entry(minus_e, i*d, i, -1.);
  }

  double * velocity_t = data->tmp_vault_nd[11];
  double * d_velocity_t = data->tmp_vault_nd[12];
  double * d_reaction_t = data->tmp_vault_nd[13];
  double * velocity_t_inv = data->tmp_vault_nd[14];
  double * Qp_velocity_t_inv = data->tmp_vault_nd[15];
  //double * tmp1 = data->tmp_vault_nd[16];
  FILE * iterates;
  FILE * matrixH;

  char *str = (char *) malloc(200);
  strcpy( str, problem_name );
  const char * separators = "/";
  char *strToken = strtok( str, separators );
  for(int i=0; i<5; i++)
  {
    if(strToken != NULL) strToken = strtok ( NULL, separators );
  }

  strToken = strtok ( strToken, "." );
  // for(int i=0; i<strlen(strToken); i++)
  // {
  //   if(strToken[i] == '-') strToken[i] = '_';
  // }

  char matlab_name[100], probName[100];

  // int count=0; for (int i = 13; problem_name[i] != '.'; i++) {probName[count] = problem_name[i]; count++;} probName[count] = '\0';
  // sprintf(matlab_name, "%s.m",probName);

  // sprintf(matlab_name, "%s.m",strToken);
  sprintf(matlab_name, "iterates.m");

  iterates = fopen(matlab_name, "r+");
 
  if (iterates == NULL)
    {
    /* file doesn't exist */
      iterates = fopen(matlab_name, "w");
      fprintf(iterates,"data = [];\n");
    }
  fclose(iterates);

  /* writing data in a Matlab file */
  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
  {
    // iterates = fopen("iterates.m", "w");
    // remove(matlab_name);
    iterates = fopen(matlab_name, "a+");
    //fprintf(iterates,"%% data = struct;\n");
    fprintf(iterates,"data(end+1).name = \"%s\";\n", strToken);
    fprintf(iterates,"data(end).val = [\n");
    // printDataProbMatlabFile(M, f, H, w, d, n, m, problem->mu, iterates);
  }

  // ComputeErrorGlobalPtr computeError = NULL;

  // if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] == 1)
  // {
  //   computeError = (ComputeErrorGlobalPtr)&gfc3d_compute_error;
  // }
  // else
  // {
  //   computeError = (ComputeErrorGlobalPtr)&gfc3d_compute_error_convex;
  // }
  /* check the full criterion */
  double norm_q = cblas_dnrm2(m, problem->q, 1);
  double norm_b = cblas_dnrm2(nd, problem->b, 1);

  double *phiu = (double*)calloc(n,sizeof(double));
  double phiu_val = 0.0, alpha_diff_phi = 0.0;
  double tol_int = 1.;
  double kappa_eps = 0.1, kappa_mu = 1.;
  double tmp_barr_param = 0.;
  double max_uor_2mu = 0., tmp_uor_2mu = 0.;
  int findKappa = 1;
  double old_diff_fixp = 0.;


  ComputeErrorGlobalPtr computeError = NULL;
  computeError = (ComputeErrorGlobalPtr)&gfc3d_compute_error;


while(hasNotConverged != 0 && findKappa)
{
  findKappa = 0;
  // Reset vars
  // r
  for (unsigned int  i = 0; i<nd; i++)
      if (i % d == 0) reaction[i] = 0.1;
      else reaction[i] = 0.01;
  // v
  for (unsigned int  i = 0; i<m; i++) globalVelocity[i] = f[i];
  NM_tgemv(1.0, H, reaction, 1.0, globalVelocity);
  NM_Cholesky_solve(M, globalVelocity, 1);
  // u
  for (unsigned int  i = 0; i<nd; i++)
      if (i % d == 0) velocity[i] = 0.1;
      else velocity[i] = 0.01;

  // Reset params
  iteration = 0;
  hasNotConverged = 1;
  pinfeas = dinfeas = complem = udotr = projerr = diff_fixp = totalresidual = old_diff_fixp = 1e300;
  alpha_primal = alpha_dual = 1.;
  barr_param = 1.0;
  // barr_param = 1e-5;
  // kappa_mu = randomFloat(0.1, 0.99);
  // kappa_eps = randomFloat(0.7, 3.0)*n;
  // kappa_mu = 0.99;
  kappa_mu = 0.5;
  kappa_eps = n;
  // kappa_eps = 10;



  primalResidual_s(velocity, H, globalVelocity, w, s, primalConstraint, &pinfeas, tol);
  dualResidual(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol);
  diff_fixp = 0.;
  for (unsigned int i = 0; i<nd; i+=d)
  {
    diff_fixp += (s[i/d] - cblas_dnrm2(2, velocity+i+1, 1))*(s[i/d] - cblas_dnrm2(2, velocity+i+1, 1));
  }
  diff_fixp = sqrt(diff_fixp);
  old_diff_fixp = diff_fixp;
  JA_prod(velocity, reaction, nd, n, complemConstraint);
  for (int k = 0; k < nd; complemConstraint[k] -= 2*barr_param, k+=d);
  // complem = complemResidualNorm(velocity, reaction, nd, n);
  // udotr = cblas_ddot(nd, velocity, 1, reaction, 1)/n;
  // projerr = projectionError(velocity, reaction, n, tol);


  while(iteration < max_iter)
  {
    // if (totalresidual >= 1e-2 && totalresidual < 1e-1)
    // {
    //   kappa_mu = 0.6;
    // }

    // if (totalresidual >= 1e-3 && totalresidual < 1e-2)
    // {
    //   kappa_mu = 0.4;
    // }

    // if (barr_param < 1e-3)
    // {
    //   kappa_mu = 0.99;
    // }



    // /* Computation of the values of
    //  - primal residual: u - H*v - w - phi(s)
    //  - dual residual: M*v - f - H'*r
    //  - duality gap: u'*r
    //  - true duality gap: (value_of_primal - value_of_dual)
    //  - complementarity: u o r
    //  - projection error: r - proj(r-u)
    // */
    // primalResidual(velocity, H, globalVelocity, w, s, primalConstraint, &pinfeas, tol);
    // dualResidual(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol);
    // complem = complemResidualNorm(velocity, reaction, nd, n);
    // udotr = cblas_ddot(nd, velocity, 1, reaction, 1)/n;


    // JA_prod(velocity, reaction, nd, n, udotr_mu_vec);
    // for (int k = 0; k < nd; udotr_mu_vec[k] -= 2*barr_param, k+=d);
    // udotr_mu = cblas_dnrm2(2, udotr_mu_vec, 1);

    // diff_fixp = 0.;
    // for (unsigned int i = 0; i<nd; i+=d)
    // {
    //   diff_fixp += (s[i/d] - cblas_dnrm2(2, velocity+i+1, 1))*(s[i/d] - cblas_dnrm2(2, velocity+i+1, 1));
    // }
    // diff_fixp = sqrt(diff_fixp);


    // // for (unsigned int i = 0; i<nd; i+=d)
    // // {
    // //   printf("s[%i] = %8.20e,\t |ub[%i]| = %8.20e\n", i/d, s[i/d], i/d, cblas_dnrm2(2, velocity+i+1, 1));
    // // }
    // // printf("\n");



    // projerr = projectionError(velocity, reaction, n, tol);

    // totalresidual = fmax(fmax(fmax(pinfeas, dinfeas),udotr),diff_fixp);
    // totalresidual_mu = fmax(fmax(fmax(pinfeas, dinfeas),udotr_mu),diff_fixp);

    // // if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
    // //   printIteresProbMatlabFile(iteration, globalVelocity, velocity, reaction, s, d_globalVelocity, d_velocity, d_reaction, d_s, d, n, m, iterates);


    // // kappa_mu = 0.852;
    // // kappa_eps = 2*n;
    // // scale = 100;
    // // if (totalresidual_mu <= kappa_eps*barr_param)
    // if (totalresidual_mu <= 1e-8)
    // {
    //   barr_param *= kappa_mu;
    //   // tmp_barr_param = barr_param;
    //   // barr_param = fmax(scale_tol*tmp_barr_param, pow(tmp_barr_param,1.3));
    //   // barr_param = udotr / n / 3;
    //   // barr_param = fmax(diff_fixp, udotr) / n / 3;
    //   // tol_int = scale*barr_param;
    //   printf("\n");
    //   numerics_printf_verbose(-1, "| it  | pinfeas | dinfeas |  |s-ub| | |uor-mu||  u'r/n  | complem | prj err | barpram |  alpha  |  |dv|   |  |du|   |  |dr|   |  |ds|   | ls prim | ls dual | ls comp | ls fixP |");
    // }

    // // if (totalresidual <=)


    // if ( totalresidual <= tol )
    // {
    //   numerics_printf_verbose(-1, "| %3i%c| %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |",
    //                           iteration, fws, pinfeas, dinfeas, diff_fixp, udotr_mu, udotr, complem, projerr, barr_param);

    //   if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
    //     printIteresProbMatlabFile(iteration, globalVelocity, velocity, reaction, s, d_globalVelocity, d_velocity, d_reaction, d_s, d, n, m, iterates);
    //     // printIteresProbMatlabFile(iteration, pinfeas, dinfeas, udotr, diff_fixp, d, n, m, iterates);

    //   double unitur;
    //   for (int i = 0; i < n; i++)
    //   {
    //      unitur = cblas_ddot(3, velocity+3*i, 1, reaction+3*i, 1);
    //      if (unitur<0)
    //        printf("UR NEGATIF %9.2e\n", unitur);
    //   }

    //   hasNotConverged = 0;
    //   break;
    // }
    // if (alpha_primal < 1e-8)
    // {
    //   printf("\nfailure\n\n");

    //   printf("u o r - 2 * mu * e = \n");
    //   printBlockVec(udotr_mu_vec, nd, 3, 1);
    //   printf("\n\n");
    //   break;
    // }
    /*  Build the Jacobian matrix without reducing the linear system.
     *
     *
     *  In the case where the NT scaling is not used, the matrix to factorize is the following:
     *
     *         m        nd        nd      n
     *      |  M        0        -H^T     0 | m
     *      |                               |
     *  J = |  0      Arw(r)    Arw(u)    0 | nd
     *      |                               |
     *      | -H        I          0     -E | nd
     *      |                               |
     *      |  0  [0 -ub'/|ub|]    0      I | n
     *
     *
     * Case of a non-reduced system. Building the right-hand side related to the first system

       without NT scaling
       rhs = -
       [      M*v - H'*r - f      ]  m         dualConstraint
       [      u o r - 2 mu e      ]  nd        complemConstraint
       [     u - H*v - w - Es     ]  nd        primalConstraint
       [         s - |ub|         ]  n
    */


    int jacobian_is_nan = 0;
    switch ( options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] )
    {
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_NOSCAL:
    {
      // First linear linear system
      J = NM_create(NM_SPARSE, m + 2*nd + n, m + 2*nd + n);
      J_nzmax = M_nzmax + H_nzmax + 2*(d*3-2)*n + H_nzmax + nd + n + 2*n + n;
      NM_triplet_alloc(J, J_nzmax);
      J->matrix2->origin = NSM_TRIPLET;

      NumericsMatrix * arrow_r = Arrow_repr(reaction, nd, n);
      NumericsMatrix * arrow_u = Arrow_repr(velocity, nd, n) ;

      // Create subdiff_u
      subdiff_u = NM_create(NM_SPARSE, n, nd);
      size_t subdiff_u_nzmax = 2*n;
      NM_triplet_alloc(subdiff_u, subdiff_u_nzmax);
      NM_fill(subdiff_u, NM_SPARSE, n, nd, subdiff_u->matrix2);

      /* Matrix filling */
      size_t pos; double ub;
      for(size_t i = 0; i < n; ++i)
      {
        pos = i * d;
        ub = sqrt(velocity[pos+1]*velocity[pos+1]+velocity[pos+2]*velocity[pos+2]);

        // if (ub > sqrt(DBL_EPSILON))    // subdiff_u = ub/|ub|
        // {
        //   NM_entry(subdiff_u, i, pos+1, -velocity[pos+1]/ub);
        //   NM_entry(subdiff_u, i, pos+2, -velocity[pos+2]/ub);
        // }

        // else      // // subdiff_u = arbitrary {x, |x| <= 1}
        // {
        //   NM_entry(subdiff_u, i, pos+1, -1/sqrt(2));
        //   NM_entry(subdiff_u, i, pos+2, -1/sqrt(2));
        //   // printf("Cone %zu: ub = 0\n",i);
        // }

        // if (ub < DBL_EPSILON)
        // {
        //   printf("Cone %zu: ub = %.1e, %.1e,\tub/|ub| = %.1e, %.1e\n",i, velocity[pos+1], velocity[pos+2], -velocity[pos+1]/ub, -velocity[pos+2]/ub);
        // }

        NM_entry(subdiff_u, i, pos+1, -velocity[pos+1]/ub);
        NM_entry(subdiff_u, i, pos+2, -velocity[pos+2]/ub);

        // fixpConstraint[i] = s[i] - ub;  // fixpConstraint = s - |u_bar|
        fixpConstraint[i] = s[i] - ub;  // fixpConstraint = s - |u_bar|
      }

      NM_insert(J, M, 0, 0);
      NM_insert(J, minus_H, m + nd, 0);

      NM_insert(J, arrow_r, m, m);
      NM_insert(J, eye_nd, m + nd, m);
      NM_insert(J, subdiff_u, m + 2*nd, m);

      NM_insert(J, minus_Ht, 0, m + nd);
      NM_insert(J, arrow_u, m, m + nd);

      NM_insert(J, minus_e, m + nd, m + 2*nd);
      NM_insert(J, eye_n, m + 2*nd, m + 2*nd);

      /* regularization */
      /* NM_insert(J, NM_scalar(nd, -1e-7), m + nd, m + nd); */

      // if (iteration == 0)
      // {
      //   printf("\n\n-E = \n"); NM_display(minus_e);
      //   printf("\n\nsubdiff_u = \n"); NM_display(subdiff_u);
      // }


      if(arrow_r) { NM_free(arrow_r); arrow_r = NULL; }
      if(arrow_u) { NM_free(arrow_u); arrow_u = NULL; }
      if(subdiff_u) { NM_free(subdiff_u); subdiff_u = NULL; }



      jacobian_is_nan = NM_isnan(J);
      if (jacobian_is_nan)
      {
        numerics_printf_verbose(0, "The Jacobian matrix J contains NaN");
        break;
      }



      cblas_dcopy(m, dualConstraint, 1, rhs, 1);

      JA_prod(velocity, reaction, nd, n, complemConstraint);
      for (int k = 0; k < nd; complemConstraint[k] -= 2*barr_param, k+=d);
      cblas_dcopy(nd, complemConstraint, 1, rhs+m, 1);


      cblas_dcopy(nd, primalConstraint, 1, rhs+m+nd, 1);
      cblas_dcopy(n, fixpConstraint, 1, rhs+m+2*nd, 1);



      cblas_dscal(m + 2*nd + n, -1.0, rhs, 1);
      cblas_dcopy(m + 2*nd + n, rhs, 1, rhs_2, 1);    // rhs_2 = old rhs


      // int target = -5;
      // if(iteration==target)
      // {
      //   fprintf(iterates,"Jac_C = [\n");
      //   CSparseMatrix_print_in_Matlab_file(NM_triplet(J), 0, iterates);
      //   fprintf(iterates,"];\n");
      //   fprintf(iterates,"Jac_C = sparse(int32(Jac_C(:,1)), int32(Jac_C(:,2)), Jac_C(:,3));\n");


      //   fprintf(iterates,"rhs_C = [");
      //   for(int i = 0; i < m + 2*nd + n; i++)
      //   {
      //     fprintf(iterates,"%8.20e; ",rhs[i]);
      //   }
      //   fprintf(iterates,"];\n");


      //   // Compute L, U by LAPACK using dense matrix
      //   J_dense = NM_create(NM_SPARSE, J->size0, J->size1);
      //   NM_to_dense(J, J_dense);
      //   // printf("\nJ_dense"); NM_display(J_dense); printf("\n");
      //   // printf("\nJ_dense->matrix0\n");
      //   // for (int i = 0; i < J_dense->size0; i++)
      //   // {
      //   //   for (int j = 0; j < J_dense->size1; j++)
      //   //   {
      //   //     printf("%.2f ", J_dense->matrix0[i + j * J_dense->size0]);
      //   //   }
      //   //   printf("\n");
      //   // }

      //   lapack_int* ipiv = (lapack_int*)NM_iWork(J_dense, J_dense->size0, sizeof(lapack_int));
      //   lapack_int info_J = 0;

      //   DGETRF(J_dense->size0, J_dense->size1, J_dense->matrix0, J_dense->size0, ipiv, &info_J);

      //   // Display the L and U factors
      //   fprintf(iterates,"P_LAP= [\n");
      //   //printf("P_LAP = ");
      //   for(int i = 0; i < J_dense->size0; i++)
      //   {
      //     fprintf(iterates,"%d ",ipiv[i]);
      //     // printf("%d ", ipiv[i]);
      //   }
      //   fprintf(iterates,"];\n");

      //   fprintf(iterates,"L_LAP = [\n");
      //   // printf("\n\nL factor:\n");
      //   for (int i = 0; i < J_dense->size0; i++)
      //   {
      //     for (int j = 0; j < J_dense->size1; j++)
      //     {
      //         if (i > j)
      //         {
      //           fprintf(iterates,"%8.20e ",J_dense->matrix0[i + j * J_dense->size0]);
      //           // printf("%.2f ", J_dense->matrix0[i + j * J_dense->size0]);
      //         }
      //         else if (i == j)
      //         {
      //           if (J_dense->matrix0[i + j * J_dense->size0] != 0)
      //           {
      //             fprintf(iterates,"%e ", 1.);
      //             // printf("1.00 ");
      //           }
      //           else
      //           {
      //             fprintf(iterates,"%e ", 0.);
      //             // printf("0.00 ");
      //           }
      //         }
      //         else
      //         {
      //           fprintf(iterates,"%e ", 0.);
      //           // printf("0.00 ");
      //         }
      //     }
      //     fprintf(iterates,";\n");
      //     // printf("\n");
      //   }
      //   fprintf(iterates,"];\n");


      //   fprintf(iterates,"U_LAP = [\n");
      //   // printf("\nU factor:\n");
      //   for (int i = 0; i < J_dense->size0; i++)
      //   {
      //     for (int j = 0; j < J_dense->size1; j++)
      //     {
      //         if (i <= j)
      //         {
      //           fprintf(iterates,"%8.20e ",J_dense->matrix0[i + j * J_dense->size0]);
      //           // printf("%.2f ", J_dense->matrix0[i + j * J_dense->size0]);
      //         }

      //         else
      //         {
      //           fprintf(iterates,"%e ", 0.);
      //           // printf("0.00 ");
      //         }

      //     }
      //     fprintf(iterates,";\n");
      //     // printf("\n");
      //   }
      //   fprintf(iterates,"];\n");
      // }



      // SOLVE
      NM_LU_solve(J, rhs, 1);
      // cs_lusol_2(1, NM_csc(J), rhs, DBL_EPSILON);

      // if(iteration==target) cs_lusol_2(1, NM_csc(J), rhs, DBL_EPSILON, iterates);
      // else  NM_LU_solve(J, rhs, 1);

      // J_dense = NM_create(NM_SPARSE, J->size0, J->size1);
      // NM_to_dense(J, J_dense);
      // NM_LU_solve(J_dense, rhs, 1);


      // if(iteration==target)
      // {
      //   fprintf(iterates,"sol_C = [");
      //   for(int i = 0; i < m + 2*nd + n; i++)
      //   {
      //     fprintf(iterates,"%8.20e; ",rhs[i]);
      //   }
      //   fprintf(iterates,"];\n");
      // }


      cblas_dcopy(m + 2*nd + n, rhs, 1, sol, 1);
      // cblas_dcopy(m + 2*nd + n, rhs_2, 1, rhs, 1);
      NM_gemv(1.0, J, sol, -1.0, rhs_2);

      LS_norm_d = cblas_dnrm2(m, rhs_2, 1);
      LS_norm_c = cblas_dnrm2(nd, rhs_2+m, 1);
      LS_norm_p = cblas_dnrm2(nd, rhs_2+m+nd, 1);
      LS_norm_f = cblas_dnrm2(n, rhs_2+m+2*nd, 1);

      // printf("\n\ndiff = "); printBlockVec(rhs_2, m + 2*nd + n, 5);
      // printf("\n\n");

      // double LS_norm = cblas_dnrm2(m+2*nd+n, rhs_2, 1);
      // printf("\t\t\t\t\t\t\t\t\t\t\tLS_norm = %.1e\n", LS_norm);

      cblas_dcopy(m, sol, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, sol+m, 1, d_velocity, 1);
      cblas_dcopy(nd, sol+m+nd, 1, d_reaction, 1);
      cblas_dcopy(n, sol+m+2*nd, 1, d_s, 1);
      break;
    }




     /*  Build the Jacobian matrix without reducing the linear system.
     *
     *  NT scaling is used, the matrix to factorize is the following:
     *
     *         m        nd        nd      s
     *      |  M        0        -H^T     0 | m
     *      |                               |
     *  J = |  0       Qp2         I      0 | nd
     *      |                               |
     *      | -H        I          0     -e | nd
     *      |                               |
     *      |  0  [0 -ub'/|ub|]    0      I | n
     *
     *
     * Case of a non-reduced system. Building the right-hand side related to the first system

       with NT scaling
       rhs = -
       [      M*v - H'*r - f      ]  m         dualConstraint
       [            r             ]  nd        complemConstraint    r - 2 mu u^-1
       [   u - H*v - w - phi(s)   ]  nd        primalConstraint
       [         s - |ub|         ]  n
    */
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QP2:
    {
      // First linear linear system
      J = NM_create(NM_SPARSE, m + 2*nd + n, m + 2*nd + n);
      J_nzmax = M_nzmax + H_nzmax + 2*(d*3-2)*n + H_nzmax + nd + n + 2*n + n;
      NM_triplet_alloc(J, J_nzmax);
      J->matrix2->origin = NSM_TRIPLET;

      // NumericsMatrix * arrow_r = Arrow_repr(reaction, nd, n);
      // NumericsMatrix * arrow_u = Arrow_repr(velocity, nd, n) ;

      // Create subdiff_u
      subdiff_u = NM_create(NM_SPARSE, n, nd);
      size_t subdiff_u_nzmax = 2*n;
      NM_triplet_alloc(subdiff_u, subdiff_u_nzmax);
      NM_fill(subdiff_u, NM_SPARSE, n, nd, subdiff_u->matrix2);

      /* Matrix filling */
      size_t pos; double ub;
      for(size_t i = 0; i < n; ++i)
      {
        pos = i * d;
        ub = sqrt(velocity[pos+1]*velocity[pos+1]+velocity[pos+2]*velocity[pos+2]);

        NM_entry(subdiff_u, i, pos+1, -velocity[pos+1]/ub);
        NM_entry(subdiff_u, i, pos+2, -velocity[pos+2]/ub);

        fixpConstraint[i] = s[i] - ub;  // fixpConstraint = s - |u_bar|
      }

      Nesterov_Todd_vector(2, velocity, reaction, nd, n, p2);
      Qp2 = QRmat(p2, nd, n);


      NM_insert(J, M, 0, 0);
      NM_insert(J, minus_H, m + nd, 0);

      NM_insert(J, Qp2, m, m);
      NM_insert(J, eye_nd, m + nd, m);
      NM_insert(J, subdiff_u, m + 2*nd, m);

      NM_insert(J, minus_Ht, 0, m + nd);
      NM_insert(J, eye_nd, m, m + nd);

      NM_insert(J, minus_e, m + nd, m + 2*nd);
      NM_insert(J, eye_n, m + 2*nd, m + 2*nd);



      // if(arrow_r) { NM_free(arrow_r); arrow_r = NULL; }
      // if(arrow_u) { NM_free(arrow_u); arrow_u = NULL; }
      if(subdiff_u) { NM_free(subdiff_u); subdiff_u = NULL; }
      if(Qp2) {NM_free(Qp2); Qp2 = NULL;}


      jacobian_is_nan = NM_isnan(J);
      if (jacobian_is_nan)
      {
        numerics_printf_verbose(0, "The Jacobian matrix J contains NaN");
        break;
      }


      cblas_dcopy(m, dualConstraint, 1, rhs, 1);

      cblas_dcopy(nd, reaction, 1, rhs+m, 1);
      JA_inv(velocity, nd, n, velocity_inv);
      for (int k = 0; k < nd; rhs[m+k] -= 2*barr_param*velocity_inv[k], k++);

      cblas_dcopy(nd, primalConstraint, 1, rhs+m+nd, 1);
      cblas_dcopy(n, fixpConstraint, 1, rhs+m+2*nd, 1);



      cblas_dscal(m + 2*nd + n, -1.0, rhs, 1);
      cblas_dcopy(m + 2*nd + n, rhs, 1, rhs_2, 1);    // rhs_2 = old rhs




      // SOLVE
      NM_LU_solve(J, rhs, 1);




      cblas_dcopy(m + 2*nd + n, rhs, 1, sol, 1);
      // cblas_dcopy(m + 2*nd + n, rhs_2, 1, rhs, 1);
      NM_gemv(1.0, J, sol, -1.0, rhs_2);

      LS_norm_d = cblas_dnrm2(m, rhs_2, 1);
      LS_norm_c = cblas_dnrm2(nd, rhs_2+m, 1);
      LS_norm_p = cblas_dnrm2(nd, rhs_2+m+nd, 1);
      LS_norm_f = cblas_dnrm2(n, rhs_2+m+2*nd, 1);

      cblas_dcopy(m, sol, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, sol+m, 1, d_velocity, 1);
      cblas_dcopy(nd, sol+m+nd, 1, d_reaction, 1);
      cblas_dcopy(n, sol+m+2*nd, 1, d_s, 1);
      break;
    }



    /*  Build the Jacobian matrix without reducing the linear system.
     *
     *
     *  In the case where the NT scaling is not used, the matrix to factorize is the following:
     *
     *         m         nd         nd
     *      |  M         0         -H^T    | m
     *      |                              |
     *  J = |  0       Arw(r)     Arw(u)   | nd
     *      |                              |
     *      |      [1 -ub'/|ub|]           |
     *      | -H   [           ]    0      | nd
     *      |      [      I    ]           |
     *
     *
     * Case of a non-reduced system. Building the right-hand side related to the first system

       without NT scaling
       rhs = -
       [      M*v - H'*r - f      ]  m         dualConstraint
       [      u o r - 2 mu e      ]  nd        complemConstraint
       [    u - H*v - w - E|ub|   ]  nd        primalConstraint
    */
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL:
    {
      // First linear linear system
      J = NM_create(NM_SPARSE, m + 2*nd, m + 2*nd);
      J_nzmax = M_nzmax + H_nzmax + 2*(d*3-2)*n + H_nzmax + nd + 2*n;
      NM_triplet_alloc(J, J_nzmax);
      J->matrix2->origin = NSM_TRIPLET;

      NumericsMatrix * arrow_r = Arrow_repr(reaction, nd, n);
      NumericsMatrix * arrow_u = Arrow_repr(velocity, nd, n) ;

      // Create subdiff_u
      subdiff_u = NM_create(NM_SPARSE, nd, nd);
      size_t subdiff_u_nzmax = 2*n + nd;
      NM_triplet_alloc(subdiff_u, subdiff_u_nzmax);
      NM_fill(subdiff_u, NM_SPARSE, nd, nd, subdiff_u->matrix2);

      /* Matrix filling */
      size_t pos; double ub;
      for(size_t i = 0; i < n; ++i)
      {
        pos = i * d;
        ub = sqrt(velocity[pos+1]*velocity[pos+1]+velocity[pos+2]*velocity[pos+2]);

        NM_entry(subdiff_u, pos, pos, 1.);
        NM_entry(subdiff_u, pos+1, pos+1, 1.);
        NM_entry(subdiff_u, pos+2, pos+2, 1.);

        NM_entry(subdiff_u, pos, pos+1, -velocity[pos+1]/ub);
        NM_entry(subdiff_u, pos, pos+2, -velocity[pos+2]/ub);
      }

      NM_insert(J, M, 0, 0);
      NM_insert(J, minus_H, m + nd, 0);

      NM_insert(J, arrow_r, m, m);
      NM_insert(J, subdiff_u, m + nd, m);

      NM_insert(J, minus_Ht, 0, m + nd);
      NM_insert(J, arrow_u, m, m + nd);

      /* regularization */
      /* NM_insert(J, NM_scalar(nd, -barr_param), m + nd, m + nd); */

      // if (iteration == 0)
      // {
      //   printf("\n\n-E = \n"); NM_display(minus_e);
      //   printf("\n\nsubdiff_u = \n"); NM_display(subdiff_u);
      // }


      if(arrow_r) { NM_free(arrow_r); arrow_r = NULL; }
      if(arrow_u) { NM_free(arrow_u); arrow_u = NULL; }
      if(subdiff_u) { NM_free(subdiff_u); subdiff_u = NULL; }



      jacobian_is_nan = NM_isnan(J);
      if (jacobian_is_nan)
      {
        numerics_printf_verbose(0, "The Jacobian matrix J contains NaN");
        break;
      }

      // if(rhs)
      // {
      //   free(rhs); rhs = NULL;
      //   free(rhs_2); rhs_2 = NULL;
      //   free(sol); sol = NULL;
      // }
      // rhs = (double*)calloc(m + 2*nd, sizeof(double));
      // rhs_2 = (double*)calloc(m + 2*nd, sizeof(double));
      // sol = (double*)calloc(m + 2*nd, sizeof(double));

      cblas_dcopy(m, dualConstraint, 1, rhs, 1);

      JA_prod(velocity, reaction, nd, n, complemConstraint);
      for (int k = 0; k < nd; complemConstraint[k] -= 2*barr_param, k+=d);
      cblas_dcopy(nd, complemConstraint, 1, rhs+m, 1);

      cblas_dcopy(nd, primalConstraint, 1, rhs+m+nd, 1);



      cblas_dscal(m + 2*nd, -1.0, rhs, 1);
      cblas_dcopy(m + 2*nd, rhs, 1, rhs_2, 1);    // rhs_2 = old rhs


      // SOLVE
      NM_LU_solve(J, rhs, 1);

      cblas_dcopy(m + 2*nd, rhs, 1, sol, 1);
      NM_gemv(1.0, J, sol, -1.0, rhs_2);

      LS_norm_d = cblas_dnrm2(m, rhs_2, 1);
      LS_norm_c = cblas_dnrm2(nd, rhs_2+m, 1);
      LS_norm_p = cblas_dnrm2(nd, rhs_2+m+nd, 1);

      // printf("\n\ndiff = "); printBlockVec(rhs_2, m + 2*nd + n, 5);
      // printf("\n\n");

      // double LS_norm = cblas_dnrm2(m+2*nd+n, rhs_2, 1);
      // printf("\t\t\t\t\t\t\t\t\t\t\tLS_norm = %.1e\n", LS_norm);

      cblas_dcopy(m, sol, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, sol+m, 1, d_velocity, 1);
      cblas_dcopy(nd, sol+m+nd, 1, d_reaction, 1);

      break;
    }




    default:
    {
      printf("ERROR\n");
    }
    }


    if (jacobian_is_nan)
    {
     hasNotConverged = 2;
     if (J) { J = NM_free(J); J = NULL; }
     break;
    }



    /* computing the affine step-length */
    // alpha_primal = getStepLength(velocity, d_velocity, nd, n, gmm);
    // alpha_dual = getStepLength(reaction, d_reaction, nd, n, gmm);

    // printf("alpha_primal:\n");
    alpha_primal = getStepLength(velocity, d_velocity, nd, n, 0.9);

    // printf("\nalpha_dual:\n");
    alpha_dual = getStepLength(reaction, d_reaction, nd, n, 0.9);

    // arr_alpha_primal = array_getStepLength(velocity, d_velocity, nd, n, 0.999);
    // arr_alpha_dual = array_getStepLength(reaction, d_reaction, nd, n, 0.999);
    // for(unsigned int i=0; i<n; i++)
    // {
    //   if (arr_alpha_primal[i] < 1.)
    //     printf("Cone %i: alpha_primal = %.1e\n", i, arr_alpha_primal[i]);

    //   if (arr_alpha_dual[i] < 1.)
    //     printf("Cone %i:                         , alpha_dual = %.1e\n", i, arr_alpha_dual[i]);

    // }
    // printf("\n");

    // double scale = 0.5;
    // for(unsigned int i=0; i<n; i++)
    // {
    //   // u
    //   if (arr_alpha_primal[i] <= 1e-2 && arr_alpha_dual[i] > 1.)
    //   {
    //     // d_velocity[i*d] = d_velocity[i*d+1] = d_velocity[i*d+2] = 0.;
    //     // velocity[i*d] *= scale; velocity[i*d+1] *= scale; velocity[i*d+2] *= scale;
    //   }
    //   else if (arr_alpha_primal[i] <= 1e-2 && arr_alpha_dual[i] < 1.)
    //   {
    //     // d_velocity[i*d] = d_velocity[i*d+1] = d_velocity[i*d+2] = 0.;
    //     // d_velocity[i*d] *= 1.1;
    //     // d_velocity[i*d+1] = d_velocity[i*d+2] = 0.;
    //   }

    //   // r
    //   if (arr_alpha_dual[i] <= 1e-2 && arr_alpha_primal[i] > 1.)
    //   {
    //     // d_reaction[i*d] = d_reaction[i*d+1] = d_reaction[i*d+2] = 0.;
    //     // reaction[i*d] *= scale; reaction[i*d+1] *= scale; reaction[i*d+2] *= scale;
    //   }

    //   else if (arr_alpha_dual[i] <= 1e-2 && arr_alpha_primal[i] < 1.)
    //   {
    //     // d_reaction[i*d] = d_reaction[i*d+1] = d_reaction[i*d+2] = 0.;
    //     // d_reaction[i*d] *= 1.1;
    //     // d_reaction[i*d+1] = d_reaction[i*d+2] = 0.;
    //   }
    // }

    // cblas_dscal(n, 0.999, arr_alpha_primal, 1);
    // cblas_dscal(n, 0.999, arr_alpha_dual, 1);

    // alpha_primal = alpha_dual = 1e10;
    // for(unsigned int i=0; i<n; i++)
    // {
    //   if (arr_alpha_primal[i] <= 1e-2 || arr_alpha_dual[i] <= 1e-2)
    //     printf("Cone %i: alpha_p = %e,\t\talpha_d = %e\n", i, arr_alpha_primal[i], arr_alpha_dual[i]);
    //   if (arr_alpha_primal[i] < alpha_primal) alpha_primal = arr_alpha_primal[i];
    //   alpha_primal = ((alpha_primal < 1.0) ? alpha_primal : 1.0);

    //   if (arr_alpha_dual[i] < alpha_dual) alpha_dual = arr_alpha_dual[i];
    //   alpha_dual = ((alpha_dual < 1.0) ? alpha_dual : 1.0);
    // }


    // if(arr_alpha_primal) {free(arr_alpha_primal); arr_alpha_primal = NULL;}
    // if(arr_alpha_dual) {free(arr_alpha_dual); arr_alpha_dual = NULL;}

    // printf("alpha_primal = %8.20e,\t alpha_dual = %8.20e\n", alpha_primal, alpha_dual);
    // printf("iteration = %d\n", iteration);
    // if (iteration == 0)
    // {
    //   printf("d_velocity = \n"); NV_display(d_velocity, nd);
    //   printf("\nd_reaction = \n"); NV_display(d_reaction, nd);
    //   printf("\nd_s = \n"); NV_display(d_s, n);
    //   printf("\n");
    // }



    if (alpha_primal < alpha_dual)
      alpha_dual = alpha_primal;
    else
      alpha_primal = alpha_dual;









    // /* updating the gamma parameter used to compute the step-length */
    // gmm = gmmp1 + gmmp2 * fmin(alpha_primal, alpha_dual);

    // /* ----- Predictor step of Mehrotra ----- */
    // cblas_dcopy(nd, velocity, 1, v_plus_dv, 1);
    // cblas_dcopy(nd, reaction, 1, r_plus_dr, 1);
    // cblas_daxpy(nd, alpha_primal, d_velocity, 1, v_plus_dv, 1);
    // cblas_daxpy(nd, alpha_dual, d_reaction, 1, r_plus_dr, 1);

    // /* affine barrier parameter */
    // barr_param_a = cblas_ddot(nd, v_plus_dv, 1, r_plus_dr, 1) / n / 3;

    // /* computing the centralization parameter */
    // e = barr_param > sgmp1 ? fmax(1.0, sgmp2 * pow(fmin(alpha_primal, alpha_dual),2)) : sgmp3;
    // sigma = fmin(1.0, pow(barr_param_a / barr_param, e));

    // /* Computing the corrector step of Mehrotra algorithm */

    // /*      without NT scaling */
    // /*      rhs = - */
    // /*      [         M*v - H'*r - f             ]  m         dualConstraint */
    // /*      [ u o r + da^u o da^r - 2*sigma*mu*e ]  nd        complemConstraint */
    // /*      [         u - H*v - w                ]  nd        primalConstraint */


    // switch ( options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] )
    // {
    // case SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_NOSCAL:
    // {
    //   // Second linear linear system
    //   JA_prod(d_velocity, d_reaction, nd, n, dvdr_jprod);
    //   cblas_daxpy(nd, -1.0, dvdr_jprod, 1, rhs+m, 1);
    //   for (int k = 0; k < nd; rhs[m+k] += 2*sigma*barr_param, k+=d);

    //   double * rhs_tmp = (double*)calloc(m+2*nd+n,sizeof(double));
    //   cblas_dcopy(m+2*nd+n, rhs, 1, rhs_2, 1);

    //   // double * sol = (double*)calloc(m+2*nd+n,sizeof(double));

    //   // for (int itr = 0; itr < 1; itr++)
    //   // {
    //   //     NM_LU_solve(J, rhs_tmp, 1);
    //   //     cblas_daxpy(m+2*nd, 1.0, rhs_tmp, 1, sol, 1);
    //   //     cblas_dcopy(m+2*nd, rhs_2, 1, rhs_tmp, 1);
    //   //     NM_gemv(-1.0, J, sol, 1.0, rhs_tmp);
    //   //     //printf("refinement iterations = %i %8.2e\n",itr, cblas_dnrm2(m+2*nd, rhs_tmp, 1));
    //   //     // if (cblas_dnrm2(m+2*nd, rhs_tmp, 1) <= tol)
    //   //     // {
    //   //     //   break;
    //   //     // }
    //   // }

    //   NM_LU_solve(J, rhs, 1);

    //   cblas_dcopy(m + 2*nd + n, rhs, 1, sol, 1);
    //   NM_gemv(1.0, J, sol, -1.0, rhs_2);

    //   LS_norm_d = cblas_dnrm2(m, rhs_2, 1);
    //   LS_norm_c = cblas_dnrm2(nd, rhs_2+m, 1);
    //   LS_norm_p = cblas_dnrm2(nd, rhs_2+m+nd, 1);
    //   LS_norm_f = cblas_dnrm2(n, rhs_2+m+2*nd, 1);

    //   cblas_dcopy(m, sol, 1, d_globalVelocity, 1);
    //   cblas_dcopy(nd, sol+m, 1, d_velocity, 1);
    //   cblas_dcopy(nd, sol+m+nd, 1, d_reaction, 1);
    //   cblas_dcopy(n, sol+m+2*nd, 1, d_s, 1);




    //   if (rhs_tmp) { free(rhs_tmp); rhs_tmp = NULL; }
    //   // if (sol) { free(sol); sol = NULL; }

    //   break;
    // }

    // default:
    // {
    //   printf("ERROR\n");
    // }
    // }



    // alpha_primal = getStepLength(velocity, d_velocity, nd, n, gmm);
    // alpha_dual = getStepLength(reaction, d_reaction, nd, n, gmm);

    // if (alpha_primal < alpha_dual)
    //   alpha_dual = alpha_primal;
    // else
    //   alpha_primal = alpha_dual;

    // /* updating the gamma parameter used to compute the step-length at the next iteration */

    // gmm = gmmp1 + gmmp2 * fmin(alpha_primal, alpha_dual);



    /* ----- Update variables ----- */

    if (NV_isnan(d_globalVelocity, m) | NV_isnan(d_velocity, nd) | NV_isnan(d_reaction, nd) | NV_isnan(d_s, n))
    {
      hasNotConverged = 2;
      break;
    }


    // if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
    //   printIteresProbMatlabFile(iteration, globalVelocity, velocity, reaction, s, d_globalVelocity, d_velocity, d_reaction, d_s, d, n, m, iterates);
      // printIteresProbMatlabFile(iteration, pinfeas, dinfeas, udotr, diff_fixp, d, n, m, iterates);




    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL)
    {
      // s: not used -> s = |ub|, ds = 0
      for (int k = 0; k < nd; k+=d)
      {
        double ub = cblas_dnrm2(2, velocity+k+1, 1);
        s[k/d] = ub + alpha_primal*cblas_ddot(2, velocity+k+1, 1, d_velocity+k+1, 1) / ub;
        d_s[k/d] = 0.;
      }
    }


    /* double check_sub = 0.; */
    /* for (int k = 0; k < nd; k+=d) */
    /* { */
    /*   double nub = cblas_dnrm2(2, velocity+k+1, 1); */
    /*   double ndub = cblas_dnrm2(2, d_velocity+k+1, 1); */
    /*   check_sub += fabs(fabs(cblas_ddot(2, velocity+k+1, 1, d_velocity+k+1, 1))-nub*ndub); */
    /* } */
    // check_sub = sqrt(check_sub);
    // printf("\ndiff = %e\n", check_sub);


    cblas_daxpy(m, alpha_primal, d_globalVelocity, 1, globalVelocity, 1);
    cblas_daxpy(nd, alpha_primal, d_velocity, 1, velocity, 1);
    cblas_daxpy(nd, alpha_dual, d_reaction, 1, reaction, 1);
    cblas_daxpy(n, alpha_primal, d_s, 1, s, 1);
    // cblas_daxpy(n, 1., d_s, 1, s, 1);




    if (NV_isnan(globalVelocity, m) | NV_isnan(velocity, nd) | NV_isnan(reaction, nd) | NV_isnan(s, n))
    {
      hasNotConverged = 2;
      break;
    }

    // if (iteration > 20)
    // {
      // printf("\n");
      for(unsigned int i=0; i<nd; i+=d)
      {
        // if (s[i/d] < 0) printf("s[%i] = %8.20e,\t|ub[%i]| = %8.20e\n", i/d, s[i/d], i/d, cblas_dnrm2(2, velocity+i+1,1));
        // printf("s[%i] = %e,\t|ub[%i]| = %e, \ts-ub = %e\n", i/d, s[i/d], i/d, cblas_dnrm2(2, velocity+i+1,1), s[i/d]-cblas_dnrm2(2, velocity+i+1,1));
        // if (s[i/d] < 0) s[i/d] = cblas_dnrm2(2, velocity+i+1,1);
        // s[i/d] = cblas_dnrm2(2, velocity+i+1,1);
      }
      // printf("\n");
    // }






    /* Computation of the values of
     - primal residual: u - H*v - w - phi(s)
     - dual residual: M*v - f - H'*r
     - duality gap: u'*r
     - true duality gap: (value_of_primal - value_of_dual)
     - complementarity: u o r
     - projection error: r - proj(r-u)
    */
    primalResidual_s(velocity, H, globalVelocity, w, s, primalConstraint, &pinfeas, tol);
    dualResidual(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol);
    complem = complemResidualNorm(velocity, reaction, nd, n);
    udotr = cblas_ddot(nd, velocity, 1, reaction, 1)/n;

    // use complemConstraint instead of udotr_mu_vec !!!
    // JA_prod(velocity, reaction, nd, n, udotr_mu_vec);
    JA_prod(velocity, reaction, nd, n, complemConstraint);
    max_uor_2mu = 0.0;
    for (int k = 0; k < nd; k+=d)
    {
      complemConstraint[k] -= 2*barr_param;
      tmp_uor_2mu = cblas_dnrm2(3, complemConstraint+k, 1);

      // udotr_mu_vec[k] -= 2*barr_param;
      // tmp_uor_2mu = cblas_dnrm2(3, udotr_mu_vec+k, 1);
      if (tmp_uor_2mu > max_uor_2mu) max_uor_2mu = tmp_uor_2mu;
    }

    // udotr_mu = cblas_dnrm2(nd, udotr_mu_vec, 1);
    udotr_mu = cblas_dnrm2(nd, complemConstraint, 1);


    diff_fixp = 0.;
    for (unsigned int i = 0; i<nd; i+=d)
    {
      diff_fixp += (s[i/d] - cblas_dnrm2(2, velocity+i+1, 1))*(s[i/d] - cblas_dnrm2(2, velocity+i+1, 1));
    }
    diff_fixp = sqrt(diff_fixp);

    // for (unsigned int i = 0; i<nd; i+=d)
    // {
    //   printf("s[%i] = %8.20e,\t |ub[%i]| = %8.20e\n", i/d, s[i/d], i/d, cblas_dnrm2(2, velocity+i+1, 1));
    // }
    // printf("\n");



    // projerr = projectionError(velocity, reaction, n, tol);
    NM_gemv(1.0, P_mu_inv, velocity, 0.0, data->tmp_point->t_velocity);
    NM_gemv(1.0, P_mu, reaction, 0.0, data->tmp_point->t_reaction);
    (*computeError)(problem,
                    data->tmp_point->t_reaction, data->tmp_point->t_velocity, globalVelocity,
                    tol, options, norm_q, norm_b, &projerr);




    // totalresidual = fmax(fmax(fmax(pinfeas, dinfeas),udotr),diff_fixp);
    totalresidual = fmax(fmax(fmax(fmax(pinfeas, dinfeas),diff_fixp),2.*max_uor_2mu),4.*barr_param);
    totalresidual_mu = fmax(fmax(fmax(pinfeas, dinfeas),udotr_mu),diff_fixp);

    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
      // printIteresProbMatlabFile(iteration, globalVelocity, velocity, reaction, s, d_globalVelocity, d_velocity, d_reaction, d_s, d, n, m, iterates);
      fprintf(iterates,"%d %.1e %.1e %.1e %.1e %.1e %.1e %.1e %.1e %.1e %.1e %.1e %.1e %.1e %.1e %.1e %.1e %.1e %.1e;\n",
            iteration, pinfeas, dinfeas, diff_fixp, udotr_mu, 2.*max_uor_2mu, 4.*barr_param, udotr, projerr, barr_param, alpha_primal,
            fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
            fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
            fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]),
            fabs(d_s[cblas_idamax(n, d_s, 1)]),
            LS_norm_p, LS_norm_d, LS_norm_c, LS_norm_f);


    numerics_printf_verbose(-1, "| %3i%c| %.1e | %.1e | %.1e | %.1e |   %.1e  | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |",
                            iteration, fws, pinfeas, dinfeas, diff_fixp, udotr_mu, 2.*max_uor_2mu, 4.*barr_param, udotr, projerr, barr_param, alpha_primal,
                            fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
                            fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
                            fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]),
                            fabs(d_s[cblas_idamax(n, d_s, 1)]),
             LS_norm_p, LS_norm_d, LS_norm_c, LS_norm_f);




    if ( totalresidual <= tol )
    {
      // numerics_printf_verbose(-1, "| %3i%c| %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |",
      //                         iteration, fws, pinfeas, dinfeas, diff_fixp, udotr_mu, udotr, complem, projerr, barr_param);

      // if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
      //   printIteresProbMatlabFile(iteration, globalVelocity, velocity, reaction, s, d_globalVelocity, d_velocity, d_reaction, d_s, d, n, m, iterates);
        // printIteresProbMatlabFile(iteration, pinfeas, dinfeas, udotr, diff_fixp, d, n, m, iterates);

      double unitur;
      for (int i = 0; i < n; i++)
      {
         unitur = cblas_ddot(3, velocity+3*i, 1, reaction+3*i, 1);
         if (unitur<0)
           printf("UR NEGATIF %9.2e\n", unitur);
      }

      hasNotConverged = 0;
      break;
    }

    if (alpha_primal < 1e-8)
    {
      printf("\nfailure\n\n");

      // printf("u o r - 2 * mu * e = \n");
      // printBlockVec(udotr_mu_vec, nd, 3, 1);
      // printf("\n\n");
      break;
    }


    // kappa_mu = 0.852;
    // kappa_eps = 2*n;
    // scale = 100;
    // if (totalresidual_mu <= kappa_eps*barr_param)
    if (totalresidual_mu <= 0.1*barr_param) // && diff_fixp < old_diff_fixp)
    {
      barr_param *= kappa_mu;
      // printf("abs(ub'd_ub - |ub|*|d_ub|) = %e\n", check_sub);
      printf("\n");
      numerics_printf_verbose(-1, "| it  | pinfeas | dinfeas |  |s-ub| | |uor-mu||2max|uor-mu||   4*mu  |  u'r/n  | prj err | barpram |  alpha  |  |dv|   |  |du|   |  |dr|   |  |ds|   | ls prim | ls dual | ls comp | ls fixP |");
    }











    if(arr_alpha_primal) {free(arr_alpha_primal); arr_alpha_primal = NULL;}
    if(arr_alpha_dual) {free(arr_alpha_dual); arr_alpha_dual = NULL;}

    if(J) { J = NM_free(J); J = NULL;}
    if(J_dense) { J_dense = NM_free(J_dense); J_dense = NULL;}

    iteration++;
    // printf("\nw = "); NV_display(w, nd); printf("\n\n");

  } // while loop


}
// printf("\n\nFINAL \nkappa_mu = %e,\t kappa_eps = %e\n\n", kappa_mu,kappa_eps);

  /* Checking strict complementarity */
  /* For each cone i from 1 to n, one checks if u+r is in the interior of the Lorentz cone */
  /* One first computes the 3 dimensional vector somme = (u+r)/norm(u+r) */
  /* Then one checks if somme[0] > sqrt(somme[1]^2 + somme[2]^2) + ceps */

  double somme[3];
  double dur[3];
  double cesp = sqrt(DBL_EPSILON);
  double ns, ndur;
  int nsc = 0;
  int nN = 0;
  int nB = 0;
  int nR = 0;
  for (int i = 0; i < n; i++)
  {
    somme[0] = velocity[3*i] + reaction[3*i];
    somme[1] = velocity[3*i+1] + reaction[3*i+1];
    somme[2] = velocity[3*i+2] + reaction[3*i+2];
    dur[0] = velocity[3*i] - reaction[3*i];
    dur[1] = velocity[3*i+1] - reaction[3*i+1];
    dur[2] = velocity[3*i+2] - reaction[3*i+2];

    ns = somme[0] - cblas_dnrm2(2,somme+1,1);
    ndur = dur[0] - cblas_dnrm2(2,dur+1,1);
    if (ns > cesp*cblas_dnrm2(3, somme, 1))
    {
      nsc +=1;
      if (dur[0] >= cblas_dnrm2(2,dur+1,1))       nN +=1;
      else if (-dur[0] >= cblas_dnrm2(2,dur+1,1)) nB +=1;
      else                                        nR +=1;
    }
    else
      printf("cone %i %9.2e %9.2e\n", i, somme[0], cblas_dnrm2(2,somme+1, 1));
  }
  if (nsc < n)
    printf("Ratio of Strict complementarity solutions: %4i / %4i = %4.2f\n", nsc, n, (double)nsc/n);
  else
    printf("Strict complementarity satisfied: %4i / %4i  %9.2e %9.2e  %4i %4i %4i\n", nsc, n, ns, ns/cblas_dnrm2(3, somme, 1), nB, nN, nR);


  // printf("\nu_tilde = "); printBlockVec(velocity, nd, 3, 1);
  // printf("\n\nr = "); printBlockVec(reaction, nd, 3, 1);
  // printf("\n\nv = "); printBlockVec(globalVelocity, m, 3, 0);
  // printf("\n\n");

  /* printing complementarity products */

    /* double * veloprea = (double*)calloc(nd, sizeof(double)); */
    /* cblas_dcopy(nd, reaction, 1, veloprea, 1); */
    /* cblas_daxpy(nd, 1.0, velocity, 1, veloprea, 1); */
    /* for (int i = 0; i < n; i++) */
    /* { */
    /*   if ((veloprea[i*d]-cblas_dnrm2(d-1, veloprea+i*d+1, 1)) <= 1e-6) */
    /*   { */
    /*   printf("SC failure "); */
    /*   printf("%3i u: %20.14e %20.14e r: %20.14e %20.14e v+r: %20.14e\n", i, */
    /*          velocity[i*d]-cblas_dnrm2(d-1, velocity+i*d+1, 1), */
    /*          velocity[i*d]+cblas_dnrm2(d-1, velocity+i*d+1, 1), */
    /*          reaction[i*d]-cblas_dnrm2(d-1, reaction+i*d+1, 1), */
    /*          reaction[i*d]+cblas_dnrm2(d-1, reaction+i*d+1, 1), */
    /*          veloprea[i*d]-cblas_dnrm2(d-1, veloprea+i*d+1, 1)); */
    /*   getchar(); */
    /*   } */
    /* } */
    /* free(veloprea); */

  /* determining active constraints */

    /* int j; */
    /* for (int i = 0; i < n; i++) */
    /* { */
    /* j = i*d; */
    /* if (velocity[j] <= 1e-7) */
    /* a_velo[i] = 0; */
    /* else if (velocity[j] - dnrm2l(d-1, velocity+j+1) <= 1e-7) */
    /* a_velo[i] = 1; */
    /* else */
    /* a_velo[i] = 2; */

    /* if (reaction[j] <= 1e-7) */
    /* a_reac[i] = 0; */
    /* else if (reaction[j] - dnrm2l(d-1, reaction+j+1) <= 1e-7) */
    /* a_reac[i] = 1; */
    /* else */
    /* a_reac[i] = 2; */
    /* } */

    /* for (int i = 0; i < n; i++) */
    /*   printf("%d %d %d\n", i, a_velo[i], a_reac[i]); */


  /* ----- return to original variables ------ */
  NM_gemv(1.0, P_mu_inv, velocity, 0.0, data->tmp_point->t_velocity);
  cblas_dcopy(nd, data->tmp_point->t_velocity, 1, velocity, 1);

  NM_gemv(1.0, P_mu, reaction, 0.0, data->tmp_point->t_reaction);
  cblas_dcopy(nd, data->tmp_point->t_reaction, 1, reaction, 1);

  options->dparam[SICONOS_DPARAM_RESIDU] = totalresidual; //NV_max(error, 4);
  options->iparam[SICONOS_IPARAM_ITER_DONE] = iteration;

  if(internal_allocation)
  {
    gfc3d_IPM_SNM_free(problem,options);
  }

  if(H_tilde) H_tilde = NM_free(H_tilde);
  if(minus_H) minus_H = NM_free(minus_H);
  if(H) H = NM_free(H);
  if(Minv) Minv = NM_free(Minv);
  if(minus_M) minus_M = NM_free(minus_M);
  if(minus_Ht) minus_Ht = NM_free(minus_Ht);
  if(eye_nd) eye_nd = NM_free(eye_nd);
  if(eye_n) eye_n = NM_free(eye_n);
  if(Ht) Ht = NM_free(Ht);
  if(subdiff_u) subdiff_u = NM_free(subdiff_u);
  if(minus_e) { minus_e = NM_free(minus_e); minus_e = NULL; }


  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
  {
    fprintf(iterates, "];\n\n");
    fclose(iterates);
  }

  //  fclose(dfile);
  if (r_rhs) free(r_rhs);
  if (r_rhs_2) free(r_rhs_2);
  if (r_dv) free(r_dv);
  if (r_du) free(r_du);
  if (r_dr) free(r_dr);
  if (r_adu) free(r_adu);
  if (r_adr) free(r_adr);
  if (r_dr_a) free(r_dr_a);
  if (r_du_a) free(r_du_a);
  if (r_dv_a) free(r_dv_a);
  if (Hvw) free(Hvw);
  if (a_velo) free(a_velo);
  if (a_reac) free(a_reac);
  if (rhs_2) free(rhs_2);
  if (sol) free(sol);
  if (sr_rhs) free(sr_rhs);
  if (sr_rhs_2) free(sr_rhs_2);

  if (phiu) free(phiu);
  if (d_globalVelocity) free(d_globalVelocity);
  if (d_velocity) free(d_velocity);
  if (d_reaction) free(d_reaction);
  if (dw) free(dw);
  if (s) free(s);
  if (d_s) free(d_s);
  if (w_ori) free(w_ori);
  if (fixpConstraint) free(fixpConstraint);
  if (velocity_inv) free(velocity_inv);

  *info = hasNotConverged;
}

void gfc3d_ipm_snm_set_default(SolverOptions* options)
{
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO] = SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_NO;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] = 0;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] = SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_NOSCAL;

  //options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING] = SICONOS_FRICTION_3D_RESCALING_BALANCING_MHHT;
  options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING] = SICONOS_FRICTION_3D_RESCALING_NO;

  //options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] = SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE] = 1;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] = SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_NO;

  options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  options->dparam[SICONOS_DPARAM_TOL] = 1e-10;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1] = 1e-10;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2] = 3.;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3] = 1.;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1] = 0.9;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2] = 0.09; //0.095

}
