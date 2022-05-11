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

#define MIN_RELATIVE_SCALING sqrt(DBL_EPSILON)

const char* const   SICONOS_GLOBAL_FRICTION_3D_IPM_STR = "GFC3D IPM";


typedef struct
{
  double * globalVelocity;
  double * reaction;
  double * velocity;
}
  IPM_starting_point;

typedef struct
{
  double * t_globalVelocity;
  double * t_reaction;
  double * t_velocity;
}
  IPM_tmp_point;

typedef struct
{
  NumericsMatrix* mat;
  NumericsMatrix* inv_mat;
}
  IPM_change_of_variable;

typedef struct
{
  double alpha_primal; // primal step length
  double alpha_dual;   // dual step length
  double sigma;        // centering parameter
  double barr_param;   // barrier parameter
}
  IPM_internal_params;


typedef struct
{

  /* initial interior points */
  IPM_tmp_point* tmp_point;

  /* initial interior points */
  IPM_starting_point* starting_point;

  /* change of variable matrix */
  IPM_change_of_variable* P_mu;

  /* initial internal solver parameters */
  IPM_internal_params* internal_params;

  double **tmp_vault_nd;
  double **tmp_vault_m;
}
  Gfc3d_IPM_init_data;

/* typedef struct */
/* { */
/*   NumericsMatrix * mat; */
/*   NumericsMatrix * inv; */
/*   NumericsMatrix * sqr; */
/* } */
/*   NTmatrix; */


typedef long double float_type;
/* typedef double float_type; */

#include "gfc3d_ipm.h"

/* ------------------------- Helper functions implementation ------------------------------ */



/* Returns the maximum step-length to the boundary reduced by a factor gamma. Uses long double. */
double getStepLength(const double * const x, const double * const dx, const unsigned int vecSize,
                     const unsigned int varsCount, const double gamma)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  unsigned int pos;
  float_type aL, bL, cL, dL, alphaL;
  double min_alpha;

  min_alpha = 1e20; //1.0;

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
    min_alpha = ((alphaL < min_alpha) ? alphaL : min_alpha);
  }
  min_alpha = gamma*min_alpha;
  min_alpha = ((min_alpha < 1.0) ? min_alpha : 1.0);
  return min_alpha;
}

/* Returns the primal constraint vector for global fricprob: out = velocity - H x globalVelocity - w */
/* and the relative 2-norm of this vector: |out|/max{|velocity|, |H x globalVelocity|, |w|} */
void primalResidual(const double * velocity, NumericsMatrix * H, const double * globalVelocity, const double * w,
                    double * out, double * rnorm, const double tol)
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
  rn = fmax(rn, cblas_dnrm2(nd, velocity, 1));
  rn = fmax(rn, cblas_dnrm2(nd, w, 1));
  *rnorm = (rn > tol ? cblas_dnrm2(nd, out, 1)/rn : cblas_dnrm2(nd, out, 1));
  /* *rnorm = cblas_dnrm2(nd, out, 1);  */
}

/* Returns the dual constraint vector for global fricprob ( M*globalVelocity - f - H'*reaction ) */
void dualResidual(NumericsMatrix * M, const double * globalVelocity, NumericsMatrix * H, const double * reaction, const double * f,
                  double * out, double * rnorm, const double tol )
{
  double m = H->size1;
  double *HTr = (double*)calloc(m, sizeof(double));
  double rn;

  NM_gemv(1.0, M, globalVelocity, 0.0, out);
  rn = cblas_dnrm2(m, out, 1);
  cblas_daxpy(m, -1.0, f, 1, out, 1);
  NM_tgemv(1.0, H, reaction, 0.0, HTr);
  cblas_daxpy(m, -1.0, HTr, 1, out, 1);
  rn = fmax(rn, cblas_dnrm2(m, f, 1));
  rn = fmax(rn, cblas_dnrm2(m, HTr, 1));
  *rnorm = (rn > tol ? cblas_dnrm2(m, out, 1)/rn : cblas_dnrm2(m, out, 1));
  /* *rnorm = cblas_dnrm2(m, out, 1); */
  free(HTr);
}

/* Returns the 2-norm of the complementarity residual vector = 2-norm of the Jordan product velocity o reaction  */
double complemResidualNorm(const double * const velocity, const double * const reaction,
                           const unsigned int vecSize, const unsigned int varsCount)
{
  double * resid = (double*)calloc(vecSize, sizeof(double));
  JA_prod(velocity, reaction, vecSize, varsCount, resid);
  double norm2 = cblas_dnrm2(vecSize, resid, 1);
  free(resid);
  return norm2;
}

/* Returns the 2-norm of the complementarity residual vector = 2-norm of the Jordan product (Qp*velocity) o (Qp_inv * reaction)  */
double complemResidualNorm_p(const double * const velocity, const double * const reaction,
                             const unsigned int vecSize, const unsigned int varsCount)
{
  double * resid = (double*)calloc(vecSize, sizeof(double));
  double * u_p = (double*)calloc(vecSize, sizeof(double));
  double * r_p = (double*)calloc(vecSize, sizeof(double));
  double * a = (double*)calloc(vecSize, sizeof(double));
  double * b = (double*)calloc(vecSize, sizeof(double));

  Qx05y(velocity, reaction, vecSize, varsCount,a);
  Jsqrtinv(a, vecSize, varsCount, b);
  Qx05y(velocity, b, vecSize, varsCount, a);

  Qx50y(a, velocity, vecSize, varsCount, u_p);
  Qx05y(a, reaction, vecSize, varsCount, r_p);
  JA_prod(u_p, r_p, vecSize, varsCount, resid);

  double norm2 = cblas_dnrm2(vecSize, resid, 1);

  free(resid);
  free(u_p);
  free(r_p);
  free(a);
  free(b);

  return norm2;
}

/* Returns the 2-norm of the complementarity residual vector = 2-norm of the Jordan product (Qp*velocity) o (Qp_inv * reaction)  */
/* This computation is done with the formula "F" */
double complemResidualNorm_p_F(NumericsMatrix * Qp, NumericsMatrix * Qpinv,
                               const double * const velocity, const double * const reaction,
                               const unsigned int vecSize, const unsigned int varsCount)
{
  double * resid = (double*)calloc(vecSize, sizeof(double));
  double * u_p = (double*)calloc(vecSize, sizeof(double));
  double * r_p = (double*)calloc(vecSize, sizeof(double));

  NM_gemv(1.0, Qp, velocity, 0.0, u_p);
  NM_gemv(1.0, Qpinv, reaction, 0.0, r_p);
  JA_prod(u_p, r_p, vecSize, varsCount, resid);
  double norm2 = cblas_dnrm2(vecSize, resid, 1);
  free(resid);
  free(u_p);
  free(r_p);
  return norm2;
}

/* computation of the duality gap  */
/* dualgap = gapVal / (1 + (abs(primal value) + abs(dual value))/2) */
double dualGap(NumericsMatrix * M, const double * f, const double * w, const double * globalVelocity, const double * reaction, const unsigned int nd, const unsigned int m)
{
  double * Mv = (double*)calloc(m, sizeof(double));
  double vMv, pval, dval;

  NM_gemv(0.5, M, globalVelocity, 0.0, Mv);
  vMv = cblas_ddot(m, globalVelocity, 1, Mv, 1);
  free(Mv);
  pval = vMv - cblas_ddot(m, f, 1, globalVelocity, 1);
  dval = -vMv - cblas_ddot(nd, w, 1, reaction, 1);
  return (pval - dval)/ (1 + (fabs(pval) + fabs(dval))/2);
}

/* Rel gap = gapVal / (1 + abs(primal value) + abs(dual value)) */
double relGap(NumericsMatrix * M, const double * f, const double * w, const double * globalVelocity, const double * reaction, const unsigned int nd, const unsigned int m, const double gapVal)
{
  double * Mv = (double*)calloc(m, sizeof(double));
  double vMv, pval, dval;

  NM_gemv(0.5, M, globalVelocity, 0.0, Mv);
  vMv = cblas_ddot(m, globalVelocity, 1, Mv, 1);
  free(Mv);
  pval = vMv - cblas_ddot(m, f, 1, globalVelocity, 1);
  dval = -vMv - cblas_ddot(nd, w, 1, reaction, 1);
  return gapVal / (1 + fabs(pval) + fabs(dval));
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


void setErrorArray(double * error, const double pinfeas, const double dinfeas, const double udotr, const double dualgap, const double complem, const double projerr)
{
  error[0] = pinfeas;
  error[1] = dinfeas;
  error[2] = udotr;
  error[3] = dualgap;
  error[4] = complem;
  error[5] = projerr;
}

/* Return the 2-norm of the difference between two vectors */
double norm2VecDiff (const double * vec1, const double * vec2, const unsigned int vecSize)
{
  double *vecDiff;
  double nvd;
  vecDiff = (double*)calloc(vecSize,sizeof(double));
  cblas_dcopy(vecSize, vec1, 1, vecDiff, 1);
  cblas_daxpy(vecSize, -1.0, vec2, 1, vecDiff, 1);
  nvd = cblas_dnrm2(vecSize, vecDiff, 1);
  free(vecDiff);
  return nvd;
}

/* Returns the product Q_{p}*H where p is the NT vector related to the pair (x,y) and H is the matrix of the linear constraint. */
static  NumericsMatrix *  QNTpH(const double * const x, const double * const y, NumericsMatrix* H, const unsigned int vecSize, const size_t varsCount)
{
  double * a = (double*)calloc(vecSize, sizeof(double));
  double * b = (double*)calloc(vecSize, sizeof(double));

  Qx05y(x, y, vecSize, varsCount,a);
  Jsqrtinv(a, vecSize, varsCount, b);
  Qx05y(x, b, vecSize, varsCount, a);

  NumericsMatrix * QpH = NM_new();


  //Qx50y(a, z, vecSize, varsCount, out);
  NM_types storage = H->storageType;
  switch(storage)
  {
    /* case NM_DENSE: */
    /*   cblas_dgemv(CblasColMajor, CblasNoTrans, sizeY, sizeX, alpha, A->matrix0, sizeY, x, 1, beta, y, 1); */
    /*   break; */
    /* /\* SparseBlock storage *\/ */
    /* case NM_SPARSE_BLOCK: */
    /*   SBM_gemv_3x3(sizeX, sizeY, A->matrix1, x, y); */
    /*   break; */
    /* coordinate */
  case NM_SPARSE:
  {
    CSparseMatrix* H_csc = NM_csc(H);
    CSparseMatrix* QpH_csc  = cs_spalloc (H->size0, H->size1, H_csc->nzmax , 1, 0) ;        /* allocate result */

    CS_INT i, p, *Hp, *Hi ;
    Hp = H_csc->p ; Hi = H_csc->i ;

    CS_INT  *QpHp, *QpHi ;
    QpHp = QpH_csc->p ;

    CS_ENTRY *Hx, *QpHx ;
    Hx = H_csc->x ;


    unsigned int dimension = (int)(vecSize / varsCount); // should be used to work properly on different dimension of cones.

    double * z_beta = b;
    for (CS_INT k=0 ; k <  H->size0; k++)
      z_beta[k]=0;

    int * beta = (int *)malloc(H->size0/3 *sizeof(int));
    for (CS_INT k=0 ; k <  H->size0/3; k++)
      beta[k] =-1;


    CS_INT nz = 0;
    for (CS_INT k = 0 ; k < H->size1 ; k++)
    {
      /* search for beta and z_beta that are non null in the column k of H */
      CS_INT n_beta=-1, beta_old=-1;
      for (p = Hp [k] ; p < Hp [k+1] ; p++)
      {
        i = Hi[p];
        z_beta[i] =  Hx[p];
        CS_INT beta_current = i/3;
        if (beta_old != beta_current )
        {
          n_beta++;
          beta_old=beta_current;
        }
        beta[n_beta] = beta_current;
      }
      n_beta++;

      QpHp[k] = nz ;                   /* column k of QpH starts here */

      /* reallocate if needed */
      if ((nz + H->size0)> QpH_csc->nzmax && !cs_sprealloc (QpH_csc, 2*(QpH_csc->nzmax)+H->size0))
      {
        return NULL;             /* out of memory */
      }
      QpHi = QpH_csc->i ; QpHx = QpH_csc->x ;         /* C->i and C->x may be reallocated */

      /* multiplication and storage */
      for (int b =0; b< n_beta;b++)
      {
        CS_INT alpha=beta[b];

        Qx50y(&a[alpha*3], &z_beta[alpha*3], 3, 1, QpHx+nz);

        /* store out in QpH */
        QpHi[nz++]  = alpha*3;
        QpHi[nz++]  = alpha*3+1;
        QpHi[nz++]  = alpha*3+2;

        z_beta[alpha*3] =0.0;
        z_beta[alpha*3+1] =0.0;
        z_beta[alpha*3+2] =0.0;
        beta[b]=-1;
      }
    } // end loop k

    QpHp[H->size1] = nz ;
    cs_sprealloc (QpH_csc, 0) ;

    QpH->storageType=H->storageType;
    numericsSparseMatrix(QpH)->csc = QpH_csc;
    QpH->size0 = (int)QpH->matrix2->csc->m;
    QpH->size1 = (int)QpH->matrix2->csc->n;
    numericsSparseMatrix(QpH)->origin = NSM_CSC;
    free(beta);
    break;
  }
  break;
  default:
    fprintf(stderr, "Numerics, GFC3D IPM, QNTpH failed, unknown storage type for H.\n");
    exit(EXIT_FAILURE);
  }
  free(a);
  free(b);
  return QpH;
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

  fprintf(file,"H = [");
  CSparseMatrix_print_in_Matlab_file(NM_triplet(H), 0, file);
  fprintf(file,"];\n");
  fprintf(file,"H = sparse(int32(H(:,1)), int32(H(:,2)), H(:,3));\n");

  fprintf(file,"f = [");
  for(int i = 0; i < m; i++)
  {
    fprintf(file,"%22.14e; ",-f[i]);
  }
  fprintf(file,"];\n");

  fprintf(file,"w = [");
  for(int i = 0; i < n*d; i++)
  {
    fprintf(file,"%22.14e; ",w[i]);
  }
  fprintf(file,"];\n");

  fprintf(file,"mu = [");
  for(int i = 0; i < n; i++)
  {
    fprintf(file,"%22.14e; ",mu[i]);
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
static void printIteresProbMatlabFile(int iteration, double * v, double * u, double * r, int d, int n, int m, FILE * file)
{
  fprintf(file,"v(%3i,:) = [",iteration+1);
  for(int i = 0; i < m; i++)
  {
    fprintf(file, "%20.16e, ", v[i]);
  }
  fprintf(file,"];\n");

  fprintf(file,"u(%3i,:) = [",iteration+1);
  for(int i = 0; i < n*d; i++)
  {
    fprintf(file, "%20.16e, ", u[i]);
  }
  fprintf(file,"];\n");

  fprintf(file,"r(%3i,:) = [",iteration+1);
  for(int i = 0; i < n*d; i++)
  {
    fprintf(file, "%20.16e, ", r[i]);
  }
  fprintf(file,"];\n");

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



/* --------------------------- Interior-point method implementation ------------------------------ */
/*
 * Implementation contains the following functions:
 *  - gfc3d_IPM_init - initialize solver (allocate memory)
 *  - gfc3d_IPM_free - deallocate memory
 *  - gfc3d_IPM_setDefaultSolverOptions - setup default solver parameters
 *  - gfc3d_IPM - optimization method
 */
void gfc3d_IPM_init(GlobalFrictionContactProblem* problem, SolverOptions* options)
{
  unsigned int m = problem->M->size0;
  unsigned int nd = problem->H->size1;
  unsigned int d = problem->dimension;

  if(!options->dWork || options->dWorkSize != (size_t)(m + nd + nd))
  {
    options->dWork = (double*)calloc(m + nd + nd, sizeof(double));
    options->dWorkSize = m + nd + nd;
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

void gfc3d_IPM_free(GlobalFrictionContactProblem* problem, SolverOptions* options)
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

void gfc3d_IPM(GlobalFrictionContactProblem* restrict problem, double* restrict reaction,
               double* restrict velocity, double* restrict globalVelocity,
               int* restrict info, SolverOptions* restrict options)
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

  //globalFrictionContact_display(problem);

  /* symmetrization of the matrix M */
  if(!(NM_is_symmetric(problem->M)))
  {
    printf("#################### SYMMETRIZATION ####################\n");
    NumericsMatrix *MT = NM_transpose(problem->M);
    problem->M = NM_add(1/2., problem->M, 1/2., MT );
    //problem->M = Msym;
    NM_clear(MT);
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

  /* NM_display_storageType(problem->H); */
  /* NM_display(problem->H); */
  /* NV_display(problem->b, nd); */
  /* for (int k = 0; k < nd;  k++) */
  /*   problem->b[k]=0.0; */
  /* for(int  i = 0; i<nd*m; i++) */
  /*   problem->H->matrix0[i] *= 0.1; */


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

  // initialize solver if it is not set
  int internal_allocation=0;
  if(!options->dWork || (options->dWorkSize != (size_t)(m + nd + nd)))
  {
    gfc3d_IPM_init(problem, options);
    internal_allocation = 1;
  }

  Gfc3d_IPM_init_data * data = (Gfc3d_IPM_init_data *)options->solverData;
  NumericsMatrix *P_mu = data->P_mu->mat;
  NumericsMatrix *P_mu_inv = data->P_mu->inv_mat;

  double *w_tilde = problem->b;
  double *w = data->tmp_vault_nd[0];
  double *f = problem->q;

  /* TO TEST IF THE PROBLEM TO SOLVE IS FEASIBLE OR NOT */
  /* problem->M = NM_eye(m); */
  /* for(int  i = 0; i<m; i++) */
  /*   { */
  /*     f[i] = 0.0; */
  /*   } */

  /* for(int  i = 0; i<m; i++) */
  /*   { */
  /*     printf("%3.0i %10.6f\n", i, f[i]); */
  /*   } */

  /* f[1] = 1; */

  double *iden;

  // change of variable
  // H_tilde --> H
  NumericsMatrix *H = NM_multiply(P_mu, H_tilde);

  // w_tilde --> w
  NM_gemv(1.0, P_mu, w_tilde, 0.0, w);
  //cs_print(NM_triplet(H),0);
  //printf("#################### NORM(w) = %g    NORM(f) = %g\n", NV_norm_2(w, nd), NV_norm_2(f,m));

  // compute -H
  NumericsMatrix *minus_H = NM_create(H->storageType, H->size0, H->size1);
  NM_copy(H, minus_H);
  NM_scal(-1.0, minus_H);

  double alpha_primal = data->internal_params->alpha_primal;
  double alpha_dual = data->internal_params->alpha_dual;
  double barr_param = data->internal_params->barr_param;
  double sigma = data->internal_params->sigma;

  cblas_dcopy(nd, data->starting_point->reaction, 1, reaction, 1);
  cblas_dcopy(nd, data->starting_point->velocity, 1, velocity, 1);
  cblas_dcopy(m, data->starting_point->globalVelocity, 1, globalVelocity, 1);

  /* COMPUTATION OF A NEW STARTING POINT */

  // set the reaction vector to an arbitrary value in the interior of the cone
  for (unsigned int  i = 0; i<nd; i++)
      if (i % d == 0)
	reaction[i] = 0.1;
      else
	reaction[i] = 0.01;

  // computation of the global velocity vector: v = M\(H'*r+f)
  for (unsigned int  i = 0; i<m; i++) globalVelocity[i] = f[i];
  NM_tgemv(1.0, H, reaction, 1.0, globalVelocity);
  NM_Cholesky_solve(M, globalVelocity, 1);

  // computation of the velocity u = proj(H*v + w), then move u in the interior of the cone
  /* for (unsigned int  i = 0; i<nd; i++) */
  /*   velocity[i] = w[i]; */
  /* NM_gemv(1.0, H, globalVelocity, 1.0, velocity); */
  /* for (unsigned int i = 0; i < n; i++) */
  /*   { */
  /*     projectionOnCone(velocity+3*i, 0.1); */
  /*     if (velocity[3*i] == 0) */
  /* 	{ */
  /* 	  velocity[3*i] = 0.1; */
  /* 	  velocity[3*i+1] = 0.01; */
  /* 	  velocity[3*i+2] = 0.01; */
  /* 	} */
  /*   } */

  for (unsigned int  i = 0; i<nd; i++)
      if (i % d == 0)
	velocity[i] = 0.1;
      else
	velocity[i] = 0.01;


  double tol = options->dparam[SICONOS_DPARAM_TOL];
  unsigned int max_iter = options->iparam[SICONOS_IPARAM_MAX_ITER];

  double sgmp1 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1];
  double sgmp2 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2];
  double sgmp3 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3];
  double gmmp0 = 0.98;
  double gmmp1 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1];
  double gmmp2 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2];

  int hasNotConverged = 1;
  unsigned int iteration = 0;
  double pinfeas = 1e300;
  double dinfeas = 1e300;
  double complem = 1e300;
  double complem_p = 1e300;
  double dualgap = 1e300;
  double udotr = 1e300;
  double projerr = 1e300;
  double error[6];
  double totalresidual = 1e300;

  double *primalConstraint = data->tmp_vault_nd[1];
  double *dualConstraint = data->tmp_vault_m[0];
  double *complemConstraint = data->tmp_vault_nd[2];

  primalResidual(velocity, H, globalVelocity, w, primalConstraint, &pinfeas, tol);
  dualResidual(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol);
  udotr = cblas_ddot(nd, velocity, 1, reaction, 1);
  dualgap = dualGap(M, f, w, globalVelocity, reaction, nd, m);
  complem = complemResidualNorm(velocity, reaction, nd, n);
  projerr = projectionError(velocity, reaction, n, tol);

  setErrorArray(error, pinfeas, dinfeas, udotr, dualgap, complem, projerr);

  double gmm = gmmp0;
  double barr_param_a, e;
  double norm_f = cblas_dnrm2(m, f, 1);
  double norm_w = cblas_dnrm2(nd, w, 1);


  double *d_globalVelocity = (double*)calloc(m,sizeof(double));
  double *d_velocity = (double*)calloc(nd,sizeof(double));
  double *d_reaction = (double*)calloc(nd,sizeof(double));

  double *tmpsol = (double*)calloc(m+2*nd,sizeof(double)); // temporary solution

  double *rhs = options->dWork;
  double *gv_plus_dgv = data->tmp_vault_m[1];
  double *vr_jprod = data->tmp_vault_nd[3];
  double *v_plus_dv = data->tmp_vault_nd[4];
  double *r_plus_dr = data->tmp_vault_nd[5];
  double *vr_prod_sub_iden = data->tmp_vault_nd[6];
  double *dvdr_jprod = data->tmp_vault_nd[7];

  //double * r_p = (double*)calloc(nd,sizeof(double));                          // scaling vector p
  NumericsMatrix* r_Qp = NULL;                                                // matrix Qp
  NumericsMatrix *minus_M = NM_create(M->storageType, M->size0, M->size1);    // store the matrix -M to build the matrix of the Newton linear system
  //NumericsMatrix *QpH = NM_create(H->storageType, H->size0, H->size1);        // store the matrix Qp*H
  double * r_rhs = (double*)calloc(m+nd, sizeof(double));
  double * r_rhs_2 = (double*)calloc(m+nd, sizeof(double));
  double * r_dv = (double*)calloc(m,sizeof(double));
  double * r_dr = (double*)calloc(nd,sizeof(double));
  double * r_du = (double*)calloc(nd,sizeof(double));
  double * r_dv_a = (double*)calloc(m,sizeof(double));
  double * r_dr_a = (double*)calloc(nd,sizeof(double));
  double * r_du_a = (double*)calloc(nd,sizeof(double));
  double * r_adu = (double*)calloc(nd, sizeof(double));
  double * r_adr = (double*)calloc(nd, sizeof(double));
  double r_alpha_p, r_alpha_d; /* primal and dual steplengths */
  double r_alpha_primal, r_alpha_dual;
  double r_mu, r_mu_a; /* duality gap, affine duality gap */
  NumericsMatrix *JR; /* Reduced Jacobian with NT scaling */
  long JR_nzmax;
  double * Hvw = (double*)calloc(nd, sizeof(double));
  double err = 1e300;
  char fws = ' '; /* finish without scaling */

  /* list of active constraints : = 0 if x_0 <= epsilon, = 1 if lambda_2 <= epsilon , = 3 either */
  short * a_velo = (short*)calloc(n, sizeof(short));
  short * a_reac = (short*)calloc(n, sizeof(short));

  if ((options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING] == 1) && (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM] == 1))
  {
    /* Create the matrix -M to build the matrix of the reduced linear system */
    NM_copy(M, minus_M);
    NM_scal(-1.0, minus_M);
  }

  NumericsMatrix *J;

  long J_nzmax;

  size_t H_nzmax = NM_nnz(H);

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
  numerics_printf_verbose(-1, "problem dimensions d, n, m: %1i, %6i, %6i\n",d, n, m);
  numerics_printf_verbose(-1, "| it  |  dualgap  | pinfeas  | dinfeas  |  <u, r>  | complem  | proj err | barparam | alpha_p  | alpha_d  |  sigma   | |dv|/|v| | |du|/|u| | |dr|/|r| |");
  numerics_printf_verbose(-1, "--------------------------------------------------------------------------------------------------------------------------------------------");

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
  double * velocity_t = data->tmp_vault_nd[11];
  double * d_velocity_t = data->tmp_vault_nd[12];
  double * d_reaction_t = data->tmp_vault_nd[13];
  double * velocity_t_inv = data->tmp_vault_nd[14];
  double * Qp_velocity_t_inv = data->tmp_vault_nd[15];
  //double * tmp1 = data->tmp_vault_nd[16];
  FILE * iterates;
  FILE * matrixH;

  /* write matrix H in file */
  /* matrixH = fopen("matrixH.m", "w"); */
  /* CSparseMatrix_print_in_Matlab_file(NM_triplet(H), 0, matrixH); */
  /* /\* /\\* /\\\* NM_write_in_file(H, matrixH); *\\\/ *\\/ *\/ */
  /* fclose(matrixH); */

  //  FILE * dfile;

  //  dfile = fopen("dfile.m", "w");

  /* writing data in a Matlab file */
  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
  {
    iterates = fopen("iterates.m", "w");
    printDataProbMatlabFile(M, f, H, w, d, n, m, problem->mu, iterates);
  }

  ComputeErrorGlobalPtr computeError = NULL;

  if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] == 1)
  {
    computeError = (ComputeErrorGlobalPtr)&gfc3d_compute_error;
  }
  else
  {
    computeError = (ComputeErrorGlobalPtr)&gfc3d_compute_error_convex;
  }
  /* check the full criterion */
  double norm_q = cblas_dnrm2(m, problem->q, 1);
  double norm_b = cblas_dnrm2(nd, problem->b, 1);

  while(iteration < max_iter)
  {
    if ((options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING] == 1) && (totalresidual <= 1e-9) && (fws==' '))
    {
      // To solve the problem very accurately, the algorithm switches to a direct solution of the linear system without scaling and without reduction //
      options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING] = 0;
      options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM] = 0;
      fws = '*';
      // copy of the current solution into a temporary vector to evaluate the distance of this solution to the final one
      cblas_dcopy(m, globalVelocity, 1, tmpsol, 1);
      cblas_dcopy(nd, velocity, 1, tmpsol+m, 1);
      cblas_dcopy(nd, reaction, 1, tmpsol+m+nd, 1);
    }

    if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING])
    {
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD]==SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F)
      {
        Qp = NTmat(velocity, reaction, nd, n);
        Qpinv = NTmatinv(velocity, reaction, nd, n);
      }
    }

    primalResidual(velocity, H, globalVelocity, w, primalConstraint, &pinfeas, tol); // velocity - H * globalVelocity - w
    dualResidual(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol);  // M*globalVelocity - H'*reaction - f
    dualgap = dualGap(M, f, w, globalVelocity, reaction, nd, m);
    barr_param = cblas_ddot(nd, reaction, 1, velocity, 1) / n;
    //barr_param = complemResidualNorm(velocity, reaction, nd, n);
    //barr_param = (fws=='*' ? complemResidualNorm(velocity, reaction, nd, n)/n : complemResidualNorm_p(velocity, reaction, nd, n)/n) ;
    //barr_param = fabs(barr_param);

    /*
      double * rr = (double*)calloc(nd,sizeof(double));
      cblas_dcopy(nd, reaction, 1, rr, 1);
      cblas_daxpy(nd, -1.0, velocity, 1, rr, 1);
      for (int k = 0; k < n; projectionOnCone(rr+k*d, 1), k++);
      cblas_daxpy(nd, -1.0, reaction, 1, rr, 1);
      double nproj = cblas_dnrm2(nd, rr, 1);
      free(rr);
    */


    if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING])
    {
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD]==SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP)
      {
        //Nesterov_Todd_vector(0, velocity, reaction, nd, n, r_p);       // Nesterov and Todd scaling p-vector
        complem_p = complemResidualNorm_p(velocity, reaction, nd, n);  // Norm of the Jordan product of the scaled vectors velocity and reaction
      }
      else
      {
        complem_p = complemResidualNorm_p_F(Qp, Qpinv, velocity, reaction, nd, n);
      }
    }
    complem = complemResidualNorm(velocity, reaction, nd, n);
    udotr = cblas_ddot(nd, velocity, 1, reaction, 1);
    projerr = projectionError(velocity, reaction, n, tol);

    setErrorArray(error, pinfeas, dinfeas, udotr, dualgap, complem, projerr);

    /* ----- return to original variables ------ */
    //    NM_gemv(1.0, P_mu_inv, velocity, 0.0, data->tmp_point->t_velocity);
    /* cblas_dcopy(nd, data->tmp_point->t_velocity, 1, velocity, 1); */
    //    NM_gemv(1.0, P_mu, reaction, 0.0, data->tmp_point->t_reaction);
    /* cblas_dcopy(nd, data->tmp_point->t_reaction, 1, reaction, 1); */

    /* cblas_daxpy(nd, -1.0, data->tmp_point->t_velocity, 1, data->tmp_point->t_reaction, 1); */
    /* for( int k = 0; k < n; projectionOnCone(data->tmp_point->t_reaction+k*d, problem->mu[k]), k++); */

    /* (*computeError)(problem, */
    /*                 data->tmp_point->t_reaction, data->tmp_point->t_velocity, globalVelocity, */
    /*                 tol, options, */
    /*                 norm_q, norm_b,  &err); */

    // check exit condition
    totalresidual = fmax(fmax(error[0], error[1]),fmin(error[2], error[5]));
    if ( totalresidual <= tol)
    {
      numerics_printf_verbose(-1, "| %3i%c| %9.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e |",
                              iteration, fws, dualgap, pinfeas, dinfeas, udotr, complem, projerr, barr_param);

      /* if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE]) */
      /*   printIteresProbMatlabFile(iteration, globalVelocity, velocity, reaction, d, n, m, iterates); */

      hasNotConverged = 0;
      if (Qp)
	Qp = NM_free(Qp);
      if (Qpinv)
	Qpinv = NM_free(Qpinv);
      // numerics_printf_verbose(-1, "%9.2e %9.2e %9.2e\n", norm2VecDiff(tmpsol, globalVelocity, m), norm2VecDiff(tmpsol+m, velocity, nd), norm2VecDiff(tmpsol+m+nd, reaction, nd));
      break;
    }

    /*  Build the Jacobian matrix without reducing the linear system.
     *
     *
     *  In the case where the NT scaling is not used, the matrix to factorize is the following:
     *
     *         m     nd       nd
     *      |  M     0      -H^T  | m
     *      |                     |
     *  J = |  0   Arw(r)  Arw(u) | nd
     *      |                     |
     *      | -H     I        0   | nd
     *
     *
     *  In the case where the NT scaling is used, the matrix to factorize is the following:
     *
     *         m     nd       nd
     *      |  M     0      -H^T  | m
     *      |                     |
     *  J = |  0     Qp2      I   | nd
     *      |                     |
     *      | -H     I        0   | nd
     *
     *  where Qp2 is the square of the quadratic representation of the Nesterov-Todd scaling vector p
     *
     */

    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM] == 0)
    {
      J = NM_create(NM_SPARSE, m + nd + nd, m + nd + nd);
      J_nzmax = (d * d) * (m / d) + H_nzmax + 2 * (d * 3 - 2) * n + H_nzmax + nd;
      NM_triplet_alloc(J, J_nzmax);
      J->matrix2->origin = NSM_TRIPLET;

      NM_insert(J, M, 0, 0);
      NumericsMatrix * minus_Ht = NM_transpose(minus_H); // should be done only once
      NM_insert(J, minus_Ht, 0, m + nd);
      NM_free(minus_Ht);

      if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING] == 0)
      {
        NumericsMatrix * arrow_r = Arrow_repr(reaction, nd, n);
	NM_insert(J, arrow_r, m, m);
	NM_free(arrow_r);
	NumericsMatrix * arrow_v = Arrow_repr(velocity, nd, n) ;
	NM_insert(J, arrow_v, m, m + nd);
	NM_free(arrow_v);
      }
      else
      {
        Nesterov_Todd_vector(2, velocity, reaction, nd, n, p2);
        Qp2 = QRmat(p2, nd, n);
        NM_insert(J, Qp2, m, m);
        NM_insert(J, eye_nd, m, m + nd);
      }
      NM_insert(J, minus_H, m + nd, 0);
      NM_insert(J, eye_nd, m + nd, m);
      //NM_insert(J, NM_scalar(nd, -1e-8), m + nd, m + nd);
    }
    else
    {
      /* Build the Jacobian matrix corresponding to a reduced linear system. The variable du (velocity) is eliminated.
       *
       * In this case, NT scaling is always performed, but two posible ways can be done to solve the linear system.
       *
       * In a first approach, the reduced matrix is of the following form:
       *
       *         m       nd
       *      |  -M      H^T   | m
       * JR = |                |
       *      |   H     Qp^2   | nd
       *
       * where Qp^2 is the square of the quadratic representation of the vector p.
       *
       * A second possibility, which deals with a better contioned matrix, is to factorize the following reduced matrix:
       *
       *         m       nd
       *      |  -M     QpH^T  | m
       * JR = |                |
       *      |  QpH     I     | nd
       *
       *  where QpH = Qp * H.
       *
       *  The matrix QpH is computed by means of the function QNTpH.
       */
      JR = NM_create(NM_SPARSE, m + nd, m + nd);
      JR_nzmax = (d * d) * (m / d) + H_nzmax + H_nzmax + nd;
      NM_triplet_alloc(JR, JR_nzmax);
      JR->matrix2->origin = NSM_TRIPLET;
      NM_insert(JR, minus_M, 0, 0);

      /* without regularization */
      NM_insert(JR, eye_nd, m, m);

      /* add a regularization term  */
      /* NM_insert(JR, NM_add(1.0, NM_eye(nd), delta, Qp2),m, m); */

      //r_Qp = Quad_repr(r_p, nd, n);
      //NM_copy(H, QpH); /* useful ? */
      //NM_gemm(1.0, Qp, H, 0.0, QpH);

      //QpH = NM_multiply(Qp,H); // This product should be replaced by a function returning the product Qp * vector

      NumericsMatrix * QpH = QNTpH(velocity, reaction, H, nd, n);


      //NM_compare(QpH, QpH_new,1e-10);
      //printf("compare  %i", NM_compare(QpH, QpH_new,1e-10));

      /* if (!NM_compare(QpH, QpH_new,1e-10)) */
      /* { */
      /*   exit(1); */
      /* } */
      /* //getchar(); */

      NumericsMatrix * QpHt = NM_transpose(QpH);
      NM_insert(JR,QpHt, 0, m);      // Should be useless when unsing a symmetric factorization procedure
      NM_insert(JR, QpH, m, 0);
      NM_free(QpH);
      NM_free(QpHt);
    }

    /* for(unsigned int i = 0; i < nd; ++ i) */
    /* 	w[i] = 0.0; */

    /** Correction of w to take into account the dependence
        on the tangential velocity */
    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] == 1)
    {
      printf("update w\n");
      for(unsigned int i = 0; i < nd; ++ i)
      {
        if(i % d == 0)
          {/* w[i] = w_tilde[i]/(problem->mu[(int)(i/d)]) */
	    w[i] = w_tilde[i] + sqrt(velocity[i+1]*velocity[i+1]+velocity[i+2]*velocity[i+2]);
	  }
      }
      /* for(unsigned int i = 0; i < nd; ++ i) */
      /* 	printf("%20.14e\n",w[i]); */
    }
    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
      printIteresProbMatlabFile(iteration, globalVelocity, velocity, reaction, d, n, m, iterates);


    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM] == 0)
    {
      /* Case of a non-reduced system. Building the right-hand side related to the first system

         without NT scaling
         rhs = -
         [ M*v - H'*r - f ]  m         dualConstraint
         [     u o r      ]  nd        complemConstraint
         [  u - H*v - w   ]  nd        primalConstraint

         with NT scaling
         rhs = -
         [ M*v - H'*r - f ]  m         dualConstraint
         [       r        ]  nd        complemConstraint
         [  u - H*v - w   ]  nd        primalConstraint
      */

      cblas_dcopy(m, dualConstraint, 1, rhs, 1);
      cblas_dcopy(nd, primalConstraint, 1, rhs+m+nd, 1);
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING] == 0)
      {
        JA_prod(velocity, reaction, nd, n, complemConstraint);
        cblas_dcopy(nd, complemConstraint, 1, rhs+m, 1);
      }
      else
      {
        cblas_dcopy(nd, reaction, 1, rhs+m, 1);
      }
      cblas_dscal(m + nd + nd, -1.0, rhs, 1);
    }
    else
    {
      /* building the reduced right-hand side                              */
      /*         [ M*v - f - H'*r ]  m         dualConstraint              */
      /* r_rhs = [                ]                                        */
      /*         [ -Qp*(H*v + w)  ]  nd        primalConstraint - velocity */

      NV_insert(r_rhs, m + nd, dualConstraint, m, 0);
      NM_gemv(1.0, H, globalVelocity, 0.0, Hvw);
      cblas_daxpy(nd, 1.0, w, 1, Hvw, 1);              // Hvw <- H*v + w, kept in memory for the computation of r_du


      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD]==SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP)
      {
        QNTpz(velocity, reaction, Hvw, nd, n, r_rhs+m);  // Could also be done with Qxy(p, Hvw, nd, n, r_rhs+m);
      }
      else
      {
        NM_gemv(1.0, Qp, Hvw, 0.0, r_rhs+m);
      }
      cblas_dscal(nd, -1.0, r_rhs+m, 1);
      cblas_dcopy(m+nd, r_rhs, 1, r_rhs_2, 1);         // kept a copy of r_rhs that will be used for the solution of the second linear system
    }

    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM] == 0)
    {
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING] == 0)
      {
        /* Solving non-symmetric Newton system without NT scaling via LU factorization */
        //NM_gesv_expert(J, rhs, NM_KEEP_FACTORS);
        NM_LU_solve(J, rhs, 1);
      }
      else
      {
        /* Solving full symmetric Newton system with NT scaling via LDLT factorization */
        NSM_linearSolverParams(J)->solver = NSM_HSL;
        NM_LDLT_solve(J, rhs, 1);
      }
      cblas_dcopy(m, rhs, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, rhs+m, 1, d_velocity, 1);
      cblas_dcopy(nd, rhs+m+nd, 1, d_reaction, 1);
      /* d_globalVelocity = rhs; */
      /* d_velocity = rhs + m; */
      /* d_reaction = rhs + m + nd; */
    }
    else
    {
      /* solution of the first reduced Newton system with NT scaling via LDLT factorization */
      NSM_linearSolverParams(JR)->solver = NSM_HSL;
      NM_LDLT_solve(JR, r_rhs, 1);

      /* retrieve the solutions */
      cblas_dcopy(m, r_rhs, 1, d_globalVelocity, 1);
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD]==SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP)
      {
        QNTpz(velocity, reaction, r_rhs+m, nd, n, d_reaction);
      }
      else
      {
        NM_gemv(1.0, Qp, r_rhs+m, 0.0, d_reaction);
      }

      NM_gemv(1.0, H, d_globalVelocity, 0.0, d_velocity);  // d_velocity <- H*dv
      cblas_daxpy(nd, 1.0, Hvw, 1, d_velocity, 1);         // d_velocity <- H*v + w + H*dv
      cblas_daxpy(nd, -1.0, velocity, 1, d_velocity, 1);   // d_velocity <- H*(v + dv) + w  - u
    }

    /* computing the affine step-length */
    alpha_primal = getStepLength(velocity, d_velocity, nd, n, gmm);
    alpha_dual = getStepLength(reaction, d_reaction, nd, n, gmm);

    if (alpha_primal < alpha_dual)
      alpha_dual = alpha_primal;
    else
      alpha_primal = alpha_dual;

    /* updating the gamma parameter used to compute the step-length */
    gmm = gmmp1 + gmmp2 * fmin(alpha_primal, alpha_dual);

    /* ----- Predictor step of Mehrotra ----- */
    cblas_dcopy(nd, velocity, 1, v_plus_dv, 1);
    cblas_dcopy(nd, reaction, 1, r_plus_dr, 1);
    cblas_daxpy(nd, alpha_primal, d_velocity, 1, v_plus_dv, 1);
    cblas_daxpy(nd, alpha_dual, d_reaction, 1, r_plus_dr, 1);

    /* affine barrier parameter */
    barr_param_a = cblas_ddot(nd, v_plus_dv, 1, r_plus_dr, 1) / n;
    //barr_param_a = complemResidualNorm(v_plus_dv, r_plus_dr, nd, n);
    // barr_param_a = (fws=='*' ? complemResidualNorm(v_plus_dv, r_plus_dr, nd, n)/n : complemResidualNorm_p(v_plus_dv, r_plus_dr, nd, n)/n);

    /* computing the centralization parameter */
    e = barr_param > sgmp1 ? fmax(1.0, sgmp2 * pow(fmin(alpha_primal, alpha_dual),2)) : sgmp3;
    sigma = fmin(1.0, pow(barr_param_a / barr_param, e));

    /* Computing the corrector step of Mehrotra algorithm */
    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM] == 0)
    {
    /* Case of a non-reduced system. Building the right-hand side related to the second system.

         without NT scaling
         rhs = -
         [         M*v - H'*r - f             ]  m         dualConstraint
         [ u o r + da^u o da^r - 2*sigma*mu*e ]  nd        complemConstraint
         [         u - H*v - w                ]  nd        primalConstraint

         with NT scaling
         rhs = -
         [ M*v - H'*r - f ]  m         dualConstraint
         [       r        ]  nd        complemConstraint
         [  u - H*v - w   ]  nd        primalConstraint
      */

      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING] == 0)
      {
        /* Right-hand side for non-symmetric Newton system without NT scaling */
        iden = JA_iden(nd, n);
        cblas_dscal(nd, 2 * barr_param * sigma, iden, 1);
        JA_prod(velocity, reaction, nd, n, vr_jprod);
	NV_sub(vr_jprod, iden, nd, complemConstraint);
	JA_prod(d_velocity, d_reaction, nd, n, dvdr_jprod);
	NV_add(complemConstraint, dvdr_jprod, nd, complemConstraint);
	free(iden);
      }
      else
      {
        /* Right-hand side for symmetric Newton system with NT scaling */
	QNTpz(velocity, reaction, d_velocity, nd, n, d_velocity_t);
        QNTpinvz(velocity, reaction, d_reaction, nd, n, d_reaction_t);
	JA_prod(d_velocity_t, d_reaction_t, nd, n, dvdr_jprod);
	QNTpz(velocity, reaction, velocity, nd, n, velocity_t);
	JA_inv(velocity_t, nd, n, velocity_t_inv);
	JA_prod(velocity_t_inv, dvdr_jprod, nd, n, d_velocity_t);
	QNTpz(velocity, reaction, d_velocity_t, nd, n, complemConstraint);

	NV_add(reaction, complemConstraint, nd, complemConstraint);

	JA_inv(velocity, nd, n, d_velocity_t);
	cblas_dscal(nd, 2 * barr_param * sigma, d_velocity_t, 1);
	NV_sub(complemConstraint, d_velocity_t, nd, complemConstraint);
      }
      cblas_dcopy(m, dualConstraint, 1, rhs, 1);
      cblas_dcopy(nd, complemConstraint, 1, rhs+m, 1);
      cblas_dcopy(nd, primalConstraint, 1, rhs+m+nd, 1);
      cblas_dscal(m + nd + nd, -1.0, rhs, 1);

      /* Newton system solving */
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING] == 0)
      {
        //NM_gesv_expert(J, rhs, NM_KEEP_FACTORS);
        NM_LU_solve(J, rhs, 1);
      }
      else
      {
        NSM_linearSolverParams(J)->solver = NSM_HSL;
        NM_LDLT_solve(J, rhs, 1);
      }

      if (J)
	J = NM_free(J);

      /* NM_clear(J); */
      /* free(J); */

      //d_globalVelocity = rhs;
      //d_velocity = rhs + m;
      //d_reaction = rhs + m + nd;
      cblas_dcopy(m, rhs, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, rhs+m, 1, d_velocity, 1);
      cblas_dcopy(nd, rhs+m+nd, 1, d_reaction, 1);
    }
    else
    {
      /* right-hand side for the second reduced Newton system with NT scaling */
      double *r_Qp_u = (double*)calloc(nd,sizeof(double));
      double *r_Qp_u_inv = (double*)calloc(nd,sizeof(double));
      double * r_Qp_du = (double*)calloc(nd,sizeof(double));
      double * r_dudr = (double*)calloc(nd,sizeof(double));
      double * r_ududr = (double*)calloc(nd,sizeof(double));

      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD]==SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP)
      {
	QNTpz(velocity, reaction, velocity, nd, n, r_Qp_u);  // r_Qp_u <- u_hat  (Qp formula)
      }
      else
      {
        NM_gemv(1.0, Qp, velocity, 0.0, r_Qp_u); // r_Qp_u <- u_hat   (F formula)
      }

      cblas_dcopy(nd, r_Qp_u, 1, r_Qp_du, 1);                            // r_Qp_du <- u_hat
      cblas_daxpy(nd, 1.0, r_rhs+m, 1, r_Qp_du, 1);                      // r_Qp_du <- dr_check + u_hat
      cblas_dscal(nd, -1.0, r_Qp_du, 1);                                 // r_Qp_du <- -(dr_check + u_hat) = du_hat.
                                                                         // Here we use the following formula satisfied by the affine step
                                                                         //    u_hat o dr_check + r_check o du_hat = - u_hat o r_check
                                                                         // Then we use the fact that for the NT scaling we have: u_hat = r_check.
                                                                         // Therefore du_hat = -dr_check - r_check = -dr_check - u_hat

      JA_prod(r_Qp_du, r_rhs+m, nd, n, r_dudr);                         // r_dudr <- du_hat o dr_check

      for (int k = 0; k < nd; r_dudr[k] -= 2*sigma*barr_param, k+=d);   // r_dudr <- du_hat o dr_check - 2*sigma*mu*e
      Jxinvprody(r_Qp_u, r_dudr, nd, n, r_ududr);                       // r_ududr <- u_hat^{-1} o  (du_hat o dr_check - 2*sigma*mu*e)

      cblas_daxpy(nd, -1.0, r_ududr, 1, r_rhs_2+m, 1);                  // r_rhs_2+m <- -Qp*(Hv+w) - u_hat^{-1} o (du_hat o dr_check) + 2*sigma*mu*u_hat^{-1}

      free(r_Qp_u);
      free(r_Qp_u_inv);
      free(r_Qp_du);
      free(r_dudr);
      free(r_ududr);

      /* solution of the second linear system */
      NSM_linearSolverParams(JR)->solver = NSM_HSL;
      NM_LDLT_solve(JR, r_rhs_2, 1);
      NM_clear(JR);
      free(JR);

      /* retrieve the solutions */
      NV_copy(r_rhs_2, m, r_dv);
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD]==SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP)
      {
	//       Qxy(p, r_rhs_2+m, nd, n, r_dr);
	QNTpz(velocity, reaction, r_rhs_2+m, nd, n, r_dr);
      }
      else
      {
        NM_gemv(1.0, Qp, r_rhs_2+m, 0.0, r_dr);
      }
      NM_gemv(1.0, H, r_dv, 0.0, r_du);
      NV_add(r_du, Hvw, nd, r_du);
      NV_sub(r_du, velocity, nd, r_du);

      //d_globalVelocity = r_dv;
      //d_velocity = r_du;
      //d_reaction = r_dr;
      cblas_dcopy(m, r_dv, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, r_du, 1, d_velocity, 1);
      cblas_dcopy(nd, r_dr, 1, d_reaction, 1);
    }
    if (Qp)
      Qp = NM_free(Qp);
    if (Qp2)
      Qp2 = NM_free(Qp2);
    if (Qpinv)
      Qpinv = NM_free(Qpinv);

    alpha_primal = getStepLength(velocity, d_velocity, nd, n, gmm);
    alpha_dual = getStepLength(reaction, d_reaction, nd, n, gmm);

    if (alpha_primal < alpha_dual)
      alpha_dual = alpha_primal;
    else
      alpha_primal = alpha_dual;

    /* updating the gamma parameter used to compute the step-length */
    gmm = gmmp1 + gmmp2 * fmin(alpha_primal, alpha_dual);

    numerics_printf_verbose(-1, "| %3i%c| %9.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e |",
                            iteration, fws, dualgap, pinfeas, dinfeas, udotr, complem, projerr, barr_param, alpha_primal, alpha_dual, sigma,
                            cblas_dnrm2(m, d_globalVelocity, 1)/cblas_dnrm2(m, globalVelocity, 1),
                            cblas_dnrm2(nd, d_velocity, 1)/cblas_dnrm2(nd, velocity, 1),
                            cblas_dnrm2(nd, d_reaction, 1)/cblas_dnrm2(nd, reaction, 1));

    /* ----- Update variables ----- */

    if (NV_isnan(d_globalVelocity, m) | NV_isnan(d_velocity, nd) | NV_isnan(d_reaction, nd))
    {
      hasNotConverged = 2;
      break;
    }

    int j;

    cblas_daxpy(m, alpha_primal, d_globalVelocity, 1, globalVelocity, 1);

    /* for (int k = 0; k < nd; k++) */
    /* { */
    /*   j = k*d; */
    /*   if (cblas_dnrm2(d, velocity+j, 1) < DBL_EPSILON) */
    /*   { */
    /* 	// printf("%i %e\n", k, cblas_dnrm2(d-1, velocity+j+1, 1)); */
    /* 	for (int i = 0; i < d; d_velocity[j+i] = 0.0, i++); */
    /*   } */
    /* } */

    cblas_daxpy(nd, alpha_primal, d_velocity, 1, velocity, 1);
    cblas_daxpy(nd, alpha_dual, d_reaction, 1, reaction, 1);

    /* printf("Norm velocity:  %12.8e\n", cblas_dnrm2(nd, velocity, 1)); */
    /* printf("Norm reaction:  %12.8e\n", cblas_dnrm2(nd, reaction, 1)); */
    /* printf("Norm GlobalVe:  %12.8e\n", cblas_dnrm2(m, globalVelocity, 1)); */

    /* for (int k = 0; k < n; k++) */
    /* { */
    /*   j = k*d; */
    /*   if (velocity[j] <= cblas_dnrm2(d-1, velocity+j+1, 1)) */
    /*   { */
    /* 	printf("u[%i] %e\n", j, velocity[j]-cblas_dnrm2(d-1, velocity+j+1, 1)); */
    /* 	//velocity[j] = ((1+DBL_EPSILON)*cblas_dnrm2(d-1, velocity+j+1, 1)); */
    /* 	//getchar(); */
    /*   } */
    /*   if (reaction[j] <= cblas_dnrm2(d-1, reaction+j+1, 1)) */
    /*   { */
    /* 	printf("r[%i] %e\n", j, reaction[j]-cblas_dnrm2(d-1, reaction+j+1, 1)); */
    /* 	//reaction[j] = ((1+DBL_EPSILON)*cblas_dnrm2(d-1, reaction+j+1, 1)); */
    /* 	//getchar(); */
    /*   } */
    /* } */

    if (NV_isnan(globalVelocity, m) | NV_isnan(velocity, nd) | NV_isnan(reaction, nd))
    {
      hasNotConverged = 2;
      break;
    }

    iteration++;

  } // while loop

  /* printing complementarity products */

    /* double * veloprea = (double*)calloc(nd, sizeof(double)); */
    /* cblas_dcopy(nd, reaction, 1, veloprea, 1); */
    /* cblas_daxpy(nd, 1.0, velocity, 1, veloprea, 1); */
    /* for (int i = 0; i < n; i++) */
    /* { */
    /*   if ((veloprea[i*d]-cblas_dnrm2(d-1, veloprea+i*d+1, 1)) <= 1e-6) */
    /*   { */
    /* 	printf("SC failure "); */
    /* 	printf("%3i u: %20.14e %20.14e r: %20.14e %20.14e v+r: %20.14e\n", i, */
    /* 	       velocity[i*d]-cblas_dnrm2(d-1, velocity+i*d+1, 1), */
    /* 	       velocity[i*d]+cblas_dnrm2(d-1, velocity+i*d+1, 1), */
    /* 	       reaction[i*d]-cblas_dnrm2(d-1, reaction+i*d+1, 1), */
    /* 	       reaction[i*d]+cblas_dnrm2(d-1, reaction+i*d+1, 1), */
    /* 	       veloprea[i*d]-cblas_dnrm2(d-1, veloprea+i*d+1, 1)); */
    /* 	getchar(); */
    /*   } */
    /* } */
    /* free(veloprea); */

  /* determining active constraints */
  /*
    int j;
    for (int i = 0; i < n; i++)
    {
    j = i*d;
    if (velocity[j] <= 1e-8)
    a_velo[i] = 0;
    else if (velocity[j] - dnrm2l(d-1, velocity+j+1) <= 1e-8)
    a_velo[i] = 1;
    else
    a_velo[i] = 2;

    if (reaction[j] <= 1e-8)
    a_reac[i] = 0;
    else if (reaction[j] - dnrm2l(d-1, reaction+j+1) <= 1e-8)
    a_reac[i] = 1;
    else
    a_reac[i] = 2;
    }
  */

  /* ----- return to original variables ------ */
  NM_gemv(1.0, P_mu_inv, velocity, 0.0, data->tmp_point->t_velocity);
  cblas_dcopy(nd, data->tmp_point->t_velocity, 1, velocity, 1);

  NM_gemv(1.0, P_mu, reaction, 0.0, data->tmp_point->t_reaction);
  cblas_dcopy(nd, data->tmp_point->t_reaction, 1, reaction, 1);

  options->dparam[SICONOS_DPARAM_RESIDU] = totalresidual; //NV_max(error, 4);
  options->iparam[SICONOS_IPARAM_ITER_DONE] = iteration;

  if(internal_allocation)
  {
    gfc3d_IPM_free(problem,options);
  }

  NM_clear(H_tilde);
  free(H_tilde);
  NM_clear(minus_H);
  free(minus_H);
  NM_clear(H);
  free(H);
  NM_clear(minus_M);
  free(minus_M);

  eye_nd = NM_free(eye_nd);

  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
    fclose(iterates);

  //  fclose(dfile);
  free(r_rhs);
  free(r_rhs_2);
  free(r_dv);
  free(r_du);
  free(r_dr);
  free(r_adu);
  free(r_adr);
  free(r_dr_a);
  free(r_du_a);
  free(r_dv_a);
  free(Hvw);
  free(a_velo);
  free(a_reac);

  free(tmpsol);

  free(d_globalVelocity);
  free(d_velocity);
  free(d_reaction);


  *info = hasNotConverged;
}

void gfc3d_ipm_set_default(SolverOptions* options)
{

  options->iparam[SICONOS_IPARAM_MAX_ITER] = 200;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO] = SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_NO;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] = 0;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING] = 1;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE] = 0;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM] = 0;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING] = 0;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] = SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP;

  options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING] = SICONOS_FRICTION_3D_RESCALING_BALANCING_MHHT;

  options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1] = 1e-10;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2] = 3.;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3] = 1.;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1] = 0.9;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2] = 0.095; //0.095

}
