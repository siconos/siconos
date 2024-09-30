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
#include "fc3d_compute_error.h"
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
#include "NumericsMatrix_internal.h"

#include "projectionOnCone.h"

/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "siconos_debug.h"
#include "gfc3d_ipm.h"
#include <stdarg.h>         // for va_list, va_start, va_end
#include "fc3d_Solvers.h"
#include "FrictionContactProblem.h"              // for FrictionContactProblem
#include "GlobalFrictionContactProblem.h"        // for GlobalFrictionContac...

#include "io_tools.h"
#include <time.h>

#if defined(WITH_FCLIB)
#include <hdf5.h>
#include <hdf5_hl.h>
#include <fclib.h>
#endif

#define NUM_BLOCKS 0
#define CONTACT_INDEX 1
#define BODY_INDEX 2
#define RANK_Hc 3

const char* const   SICONOS_GLOBAL_FRICTION_3D_IPM_SNM_STR = "GFC3D IPM SNM";

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

void primalResidual_s_type(const double * velocity, NumericsMatrix * H, const double * globalVelocity, const double * w,
                    const double * s, double * out, double * rnorm, const double tol, const int type)
{
  size_t nd = H->size0;
  double rn;


  /* The memory for the result vectors should be allocated using calloc
   * since H is a sparse matrix. In other case the behaviour will be undefined.*/
  //  double *Hv = (double*)calloc(nd, sizeof(double));
  //double *u_minus_Hv = (double*)calloc(nd, sizeof(double));

  NM_gemv(-1.0, H, globalVelocity, 0.0, out);
  rn = NV_norm_type(nd, out, type);
  cblas_daxpy(nd, 1.0, velocity, 1, out, 1);
  cblas_daxpy(nd, -1.0, w, 1, out, 1);

  for(unsigned int i=0; i<nd; i+=3) out[i] -= s[i/3];

  rn = fmax(rn, NV_norm_type(nd, velocity, type));
  rn = fmax(rn, NV_norm_type(nd, w, type));
  rn = fmax(rn, NV_norm_type(nd/3, s, type));
  *rnorm = (rn > tol ? NV_norm_type(nd, out, type) : NV_norm_type(nd, out, type));
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

  // fprintf(file,"v(%3i,:) = [",iteration+1);
  // for(int i = 0; i < m; i++)
  // {
  //   fprintf(file, "%8.20e, ", v[i]);
  // }
  // fprintf(file,"];\n");



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

  // fprintf(file,"dv(%3i,:) = [",iteration+1);
  // for(int i = 0; i < m; i++)
  // {
  //   fprintf(file, "%8.20e, ", dv[i]);
  // }
  // fprintf(file,"];\n");

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

  // fprintf(file,"ds(%3i,:) = [",iteration+1);
  // for(int i = 0; i < n; i++)
  // {
  //   fprintf(file, "%8.20e, ", ds[i]);
  // }
  // fprintf(file,"];\n");

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


static float randomFloat(float min, float max) {
    return min + (float)rand() / ((float)RAND_MAX / (max - min));
}


// [OLD FUNCTION] This function is used to read blocks info stored each FCLIB problem test
// Need the FCLIB dataset with blocks info stored inside
int *read_fricprob_block(const char* path, int type, int blk_index)
{
  int *out = NULL;

  if (type < NUM_BLOCKS || type > RANK_Hc)
  {
    fprintf (stderr, "ERROR: out of \"type\"\n");
    return NULL;
  }


  int is_hdf5 = check_hdf5_file(path);
  if(is_hdf5)
  {
#if defined(WITH_FCLIB)
    hid_t  file_id, group_id, subgroup_id, dataset_id, dataspace_id;
    hssize_t num_elements;
    if ((file_id = H5Fopen (path, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
    {
      fprintf (stderr, "ERROR: opening file failed\n");
      return NULL;
    }

    // Get number of blocks
    int numBlk = 0;
    group_id = H5Gopen (file_id, "/fclib_global/blocks", H5P_DEFAULT);
    H5LTread_dataset_int (file_id, "/fclib_global/blocks/N", &numBlk);


    if (type == NUM_BLOCKS)
    {
      out = (int *)calloc(1, sizeof(int));
      *out = numBlk;
    }

    else
    {
      if (blk_index < 0)
      {
        fprintf (stderr, "ERROR: block index must not be negative\n");
        return NULL;
      }
      else if (blk_index >= numBlk)
      {
        fprintf (stderr, "ERROR: block index must be in range [0, %d)\n", numBlk);
        return NULL;
      }

      char block_X[20] = "block_";
      char dest_path[50] = "";
      sprintf(block_X + strlen(block_X), "%d", blk_index);
      sprintf(dest_path + strlen(dest_path), "%s%s", "/fclib_global/blocks/", block_X);
      subgroup_id = H5Gopen (file_id, dest_path, H5P_DEFAULT);

      if (type == CONTACT_INDEX)
      {
        sprintf(dest_path + strlen(dest_path), "%s", "/contact");
        dataset_id = H5Dopen1(file_id, dest_path);
        dataspace_id = H5Dget_space(dataset_id);
        num_elements = H5Sget_simple_extent_npoints(dataspace_id); // Get the number of elements of contact vector

        out = (int *)calloc(num_elements+1, sizeof(int));
        out[0] = (int)num_elements;                                // 1st element of returned vector is the number of data elements that the vector contains
        H5LTread_dataset_int (file_id, dest_path, out+1);          // Get data
      }

      else if (type == BODY_INDEX)
      {
        sprintf(dest_path + strlen(dest_path), "%s", "/body");
        dataset_id = H5Dopen1(file_id, dest_path);
        dataspace_id = H5Dget_space(dataset_id);
        num_elements = H5Sget_simple_extent_npoints(dataspace_id);

        out = (int *)calloc(num_elements+1, sizeof(int));
        out[0] = (int)num_elements;
        H5LTread_dataset_int (file_id, dest_path, out+1);
      }

      else if (type == RANK_Hc)
      {
        sprintf(dest_path + strlen(dest_path), "%s", "/rank_H_blk");
        out = (int *)calloc(1, sizeof(int));
        H5LTread_dataset_int (file_id, dest_path, out);
      }

      H5Gclose (subgroup_id);
    }

    H5Gclose (group_id);
    H5Fclose (file_id);


#else
    numerics_error("gfc3d_IPM_SNM",
                   "Try to read an hdf5 file, while fclib interface is not active. Recompile Siconos with fclib.",
                   path);
#endif
  }
  else
    numerics_error("gfc3d_IPM_SNM", "Not a hdf5 file ", path);

  return out;
}


static void NM_insert_Arrow_to_Triplet(CSparseMatrix *triplet, const unsigned int start_i, const unsigned int start_j,
                                      const double* const vec, const unsigned int vecSize, const size_t varsCount)
{
  size_t dimension = (size_t)(vecSize / varsCount);
  size_t pos;
  size_t total_element = (dimension * 3 - 2) * varsCount;

  if (triplet->nzmax < (total_element + triplet->nz))
  {
    fprintf(stderr,
            "NM_insert_Arrow_to_Triplet: Size of allocated triplet memory is not sufficient.\n");
    exit(EXIT_FAILURE);
  }

  for (size_t i = 0; i < varsCount; ++i)
  {
    pos = i * dimension;

    triplet->x [triplet->nz] = vec[pos];
    triplet->i [triplet->nz] = start_i + pos;
    triplet->p [triplet->nz++] = start_j + pos;

    for (size_t j = 1; j < dimension; ++j)
    {
      triplet->x [triplet->nz] = vec[pos + j];
      triplet->i [triplet->nz] = start_i + pos;
      triplet->p [triplet->nz++] = start_j + pos + j;

      triplet->x [triplet->nz] = vec[pos + j];
      triplet->i [triplet->nz] = start_i + pos + j;
      triplet->p [triplet->nz++] = start_j + pos;

      triplet->x [triplet->nz] = vec[pos];
      triplet->i [triplet->nz] = start_i + pos + j;
      triplet->p [triplet->nz++] = start_j + pos + j;
    }
  }
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
               int* restrict info, SolverOptions* restrict options)
{
  // verbose = 3;

  // int type = NORM_2;
  int type = NORM_INF;

  char *blk_num_name = NULL;
  if (options->solverId == SICONOS_GLOBAL_FRICTION_3D_IPM_SNM_SEP)
  {
    blk_num_name = (char *)malloc(10*sizeof(char));
    char *blk_num_ptr = options->solverData;
    strcpy(blk_num_name, blk_num_ptr);
    free(options->solverData); options->solverData = NULL;
  }



  // the size of the problem detection
  unsigned int m = problem->M->size0;
  unsigned int nd = problem->H->size1;
  unsigned int d = problem->dimension;
  unsigned int n = problem->numberOfContacts;



  unsigned int mp2nd = m + 2*nd;
  size_t no_n = 0, no_m = 0, no_nd = 0;

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
  // for(int i = 0; i < n ; i++) problem->mu[i]=0.3;


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









  // printf("nnz H = %8zu density = %9.4f\n",NM_nnz(problem->H), NM_nnz(problem->H)/1.0/nd/m);
  printf("nnz H = %8zu density = %9.4f, \t mu = %.2e\n",NM_nnz(problem->H), NM_nnz(problem->H)/1.0/nd/m, problem->mu[0]);

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
  double *w = data->tmp_vault_nd[no_nd++];
  double *f = problem->q;

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
  sigma = 0.3;
  double alpha_complem = 0., alpha_diffixP = 0.;
  double min_alpha_primal = 1.;

  // cblas_dcopy(nd, data->starting_point->reaction, 1, reaction, 1);
  // cblas_dcopy(nd, data->starting_point->velocity, 1, velocity, 1);
  // cblas_dcopy(m, data->starting_point->globalVelocity, 1, globalVelocity, 1);

  NumericsMatrix *minus_M = NM_create(
              M->storageType, M->size0,
              M->size1);  // store the matrix -M to build the matrix of the Newton linear system
  /* Create the matrix -M to build the matrix of the reduced linear system */
  NM_copy(M, minus_M);
  NM_scal(-1.0, minus_M);


  // /* COMPUTATION OF A NEW STARTING POINT */
  // // set the reaction vector to an arbitrary value in the interior of the cone
  // for (unsigned int  i = 0; i<nd; i++)
  //     if (i % d == 0) reaction[i] = 0.1;
  //     else reaction[i] = 0.01;

  // // computation of the global velocity vector: v = M\(H'*r+f)
  // for (unsigned int  i = 0; i<m; i++) globalVelocity[i] = f[i];
  // NM_tgemv(1.0, H, reaction, 1.0, globalVelocity);
  // NM_Cholesky_solve(NM_preserve(M), globalVelocity, 1);

  // for (unsigned int  i = 0; i<nd; i++)
  //     if (i % d == 0) velocity[i] = 0.1;
  //     else velocity[i] = 0.01;


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
  double udotr = 1e300, uor_mu = 1e300;
  double projerr = 1e300, projerr_u = 1e300;
  double error[6];
  double totalresidual = 1e300, totalresidual_mu = 1e300, totalresidual_Jor = 1e300;

  double diff_fixp = 1e300, nub = 1e300;
  double *diff_fixp_vec = (double*)calloc(n,sizeof(double));

  double *primalConstraint = data->tmp_vault_nd[no_nd++];
  double *dualConstraint = data->tmp_vault_m[no_m++];
  double *complemConstraint = data->tmp_vault_nd[no_nd++];
  double *complemConstraint_mu = data->tmp_vault_nd[no_nd++];
  double *fixpConstraint = (double*)calloc(n,sizeof(double));
  double *arr_norm_theta = (double*)calloc(max_iter,sizeof(double));
  double norm_theta = 1e300;

  double gmm = gmmp1+gmmp2;
  double barr_param_a, e;
  double norm_f = cblas_dnrm2(m, f, 1);
  double norm_w = cblas_dnrm2(nd, w, 1);

  double *velocity_inv = (double*)calloc(nd,sizeof(double));
  double *d_globalVelocity = (double*)calloc(m,sizeof(double));
  double *d_velocity = (double*)calloc(nd,sizeof(double));
  double *d_reaction = (double*)calloc(nd,sizeof(double));
  double *d_s = (double*)calloc(n,sizeof(double));
  double *s = (double*)calloc(n,sizeof(double));


  double *u_plus_du = data->tmp_vault_nd[no_nd++];  // for Mehrotra
  double *r_plus_dr = data->tmp_vault_nd[no_nd++];  // for Mehrotra
  double *dudr_jprod = data->tmp_vault_nd[no_nd++];  // for Mehrotra
  double *s_plus_ds = data->tmp_vault_nd[no_n++];  // for Mehrotra
  double diff_fixp_plus = 1e300;
  double *rhs_tmp = NULL;

  double *rhs = options->dWork;
  double *rhs_2 = (double*)calloc(m+2*nd+n, sizeof(double));
  // double *rhs_tmp = (double*)calloc(m+2*nd+n,sizeof(double));
  double *sol = (double*)calloc(m+2*nd+n, sizeof(double));

  double *nub_vec = (double*)calloc(n, sizeof(double));
  int pos = 0;

  char fws = ' '; /* finish without scaling */

  /* norm of the residuals of teh second linear system */
  double LS_norm_p = 0.; // primal feasibility
  double LS_norm_d = 0.; // dual feaqsibility
  double LS_norm_c = 0.; // complementarity
  double LS_norm_f = 0.; // fixed point

  NumericsMatrix *J = NULL, *J_dense = NULL;
  NumericsMatrix *arrow_r = NULL, *arrow_u = NULL;

  long J_nzmax;
  size_t H_nzmax = NM_nnz(H);
  size_t M_nzmax = NM_nnz(M);
  size_t arrow_nzmax = 0, subdiff_u_nzmax = 0;

  /* For CLASSIFICATION BNRT */
  int nB, nN, nR, nT;
  nB = nN = nR = nT = 0;
  int *setR = (int*)calloc(n, sizeof(int));

  size_t J_nz_captured = 0, J_nz_final = 0;

  // Timer
  long clk_tck = CLOCKS_PER_SEC;
  clock_t t1, t2;
  double total_time = 0., max_time = 10;
  // ###############################


  // For QP2
  double * Qp_u = data->tmp_vault_nd[no_nd++];
  double * Qp_du = data->tmp_vault_nd[no_nd++];
  double * Qpinv_dr = data->tmp_vault_nd[no_nd++];
  double * ududr = data->tmp_vault_nd[no_nd++];
  double * Qp_ududr = data->tmp_vault_nd[no_nd++];
  int max_refine = 1;



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
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_NOSCAL:
    {
      numerics_printf_verbose(-1,"Global friction contact problem - LS solution: 4x4 no scaling\n");
      arrow_nzmax = (d*3-2)*n;
      subdiff_u_nzmax = 2*n;
      J_nzmax = M_nzmax + H_nzmax + 2*(d*3-2)*n + H_nzmax + nd + n + 2*n + n;
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_QP2:
    {
      numerics_printf_verbose(-1,"Global friction contact problem - LS solution: 4x4 NT scaling with Qp2 and Mehrotra\n");
      break;
    }

    default:
    {
      printf("ERROR\n");
    }
    }


  numerics_printf_verbose(-1, "| it  | pinfeas | dinfeas |  |s-ub| |Ma|ui'ri||  |uor|  | prj err | barpram |  sigma  ||  alpha  |  |dv|   |  |du|   |  |dr|   |  |ds|   | ls prim | ls dual | ls comp | ls fixP |");
  numerics_printf_verbose(-1, "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");


  double * p2 = data->tmp_vault_nd[no_nd++];
  NumericsMatrix* Qp2 = NULL;
  NumericsMatrix * eye_nd = NM_eye(nd);
  NumericsMatrix * eye_n = NM_eye(n);
  NumericsMatrix * subdiff_u = NULL, * mat_ub = NULL, * mat_S = NULL;


  NumericsMatrix * minus_e = NULL;
  minus_e = NM_create(NM_SPARSE, nd, n);
  size_t minus_e_nzmax = n;
  NM_triplet_alloc(minus_e, minus_e_nzmax);
  NM_fill(minus_e, NM_SPARSE, nd, n, minus_e->matrix2);
  for(size_t i = 0; i < n; ++i)
  {
    NM_entry(minus_e, i*d, i, -1.);
  }

  FILE * iterates;
  FILE * iterates_2;
  FILE *sol_file;


  char *strToken = NULL;
  char *str = (char *) malloc(200);
  strcpy( str, problem->name );
  const char * separators = "/";
  strToken = strtok( str, separators );
  for(int i=0; i<5; i++)
  {
    if(strToken != NULL) strToken = strtok ( NULL, separators );
  }



  strToken = strtok ( strToken, "." );
  // for(int i=0; i<strlen(strToken); i++)
  // {
  //   if(strToken[i] == '-') strToken[i] = '_';
  // }

  // Append the block number into the test name
  if (options->solverId == SICONOS_GLOBAL_FRICTION_3D_IPM_SNM_SEP)
  {
    strcat(strToken, blk_num_name);
  }


  char matlab_name[100], probName[100];
  sprintf(matlab_name, "iterates_Spheres_no_s.m");

  /* writing data in a Matlab file */
  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
  {
    iterates = fopen(matlab_name, "w");
    printDataProbMatlabFile(M, f, H, w, d, n, m, problem->mu, iterates);
  }

  /* check the full criterion */
  double norm_q = cblas_dnrm2(m, problem->q, 1);
  double norm_b = cblas_dnrm2(nd, problem->b, 1);

  double tmp_barr_param = 0.;
  double max_uor_2mu = 0., tmp_uor_2mu = 0.;
  double scale_sub_diff = 0.91; //1.05;
  double max_val = 0, max_val_old = 0, scale_sigma = 0, old_sigma = 0;

  ComputeErrorGlobalPtr computeError = NULL;
  computeError = (ComputeErrorGlobalPtr)&gfc3d_compute_error;

  int load_starting_point = 0, save_sol_point = 0, pertu_point = 0;

  // Reset vars
  // r
  for (unsigned int  i = 0; i<nd; i++)
    if (i % d == 0) reaction[i] = 1.;
    else reaction[i] = 0.1;

  // v
  cblas_dcopy(m, f, 1, globalVelocity, 1);
  NM_tgemv(1.0, H, reaction, 1.0, globalVelocity);
  NM_Cholesky_solve(M, globalVelocity, 1);

  // u
  for (unsigned int  i = 0; i<nd; i++)
      if (i % d == 0) velocity[i] = 1.;
      else velocity[i] = 0.1;
  // s
  for (unsigned int  i = 0; i<n; i++)
  {
    nub = cblas_dnrm2(2, velocity+i*d+1, 1);
    s[i] = nub;
  }


  if (load_starting_point)
  {
    sol_file = fopen("sol_data_copy.res", "r");
    if (!sol_file) printf("\n\ngfc3d_ipm_snm: Solution data file is not available!!! \n\n");
    else
    {
      char read_prob_name[200];
      int load_v = 0, load_u = 0, load_r = 0, c = 0, len_prob_name = 0, newlineCount = 0;

      // Traverse the problem names in data file for a match
      for (int i=0; i<1091; i++)
      {
        if (fgets(read_prob_name, sizeof(read_prob_name), sol_file) != NULL)
        {
          len_prob_name = strlen(read_prob_name);
          if (len_prob_name > 0 && read_prob_name[len_prob_name - 1] == '\n')
          {
              read_prob_name[len_prob_name - 1] = '\0';  // Replace the newline character with null terminator
          }
          if (strcmp(read_prob_name, strToken) == 0) // Problem names are matched
          {
            load_v = load_u = load_r = 1;
            break;
          }
        }

        // Go to the next problem name
        newlineCount = 0;
        while ((c = fgetc(sol_file)) != EOF)
        {
          if (c == '\n')
          {
            newlineCount++;
            if (newlineCount == 3) break; // Stop reading after 2 lines
          }
        }
      }

      // load v
      if (load_v)
      {
        for (int i=0; i < m; i++)
        {
          fscanf(sol_file, "%lf ", globalVelocity+i);
        }
        fscanf(sol_file, "\n");
      }

      // load u
      if (load_u)
      {
        for (int i=0; i < nd; i++)
        {
          fscanf(sol_file, "%lf ", velocity+i);
        }
        fscanf(sol_file, "\n");
      }

      // load r
      if (load_r)
      {
        for (int i=0; i < nd; i++)
        {
          fscanf(sol_file, "%lf ", reaction+i);
        }
        fscanf(sol_file, "\n");
      }

      // // load s
      // for (int i=0; i < n; i++)
      // {
      //   fscanf(sol_file, "%lf ", s+i);
      // }
      // fscanf(sol_file, "\n");
      printf("\ngfc3d_ipm_snm: Starting point is successfully loaded.\n\n");
    }
    fclose(sol_file);

    // Sol perturbation
    if (pertu_point)
    {
      for (int i=0; i < m; i++)
      {
        globalVelocity[i] *= 1.1;
      }

      for (int i=0; i < nd; i++)
      {
        if (i%d == 0)
        {
          velocity[i] *= 1.2;
          reaction[i] *= 1.1;
        }
        else
        {
          velocity[i] *= 0.9;
          reaction[i] *= 0.8;
        }
      }

      for (int i=0; i < n; i++)
      {
        s[i] *= 1.05;
      }
      printf("\nThe point is successfully perturbed.\n\n");
    }
  }

  // Reset params
  iteration = 0;
  hasNotConverged = 1;
  pinfeas = dinfeas = complem = udotr = projerr = diff_fixp = totalresidual = 1e300;
  alpha_primal = alpha_dual = 1.;
  barr_param = 1.;
  sigma = 0.1;


  while(iteration < max_iter)
  {
    int jacobian_is_nan = 0;
    switch ( options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] )
    {
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
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_NOSCAL:
    {
      t1 = clock();
      // Stopping test using norm type
      primalResidual_s_type(velocity, H, globalVelocity, w, s, primalConstraint, &pinfeas, tol, type);
      dualResidual_type(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol, type);
      // udotr = xdoty_type(n, nd, velocity, reaction, type);
      udotr = xdoty_type(n, nd, velocity, reaction, NORM_2_INF);
      // complem = complemResidualNorm(velocity, reaction, nd, n);
      complem = complemResidualNorm_type(velocity, reaction, nd, n, NORM_2_INF);
      barr_param = cblas_ddot(nd, reaction, 1, velocity, 1) / n;

      diff_fixp = 0.;
      for (unsigned int i = 0; i<n; i++)
      {
        pos = i*d;

        nub_vec[i] = cblas_dnrm2(2, velocity+pos+1, 1);
        diff_fixp_vec[i] = s[i]-nub_vec[i];
      }

      // Stopping test using norm type
      if (type == NORM_2)
        diff_fixp = cblas_dnrm2(n, diff_fixp_vec, 1);
      else if (type == NORM_INF)
      {
        diff_fixp = fabs(diff_fixp_vec[0]);
        for (int i=1; i<n; i++)
        {
          diff_fixp = fmax(diff_fixp, fabs(diff_fixp_vec[i]));
        }
      }
      else
      {
        fprintf(stderr, "type = %d is undefined.\n", type);
        exit(EXIT_FAILURE);
      }


      t2 = clock();
      total_time += (double)(t2-t1)/(double)clk_tck;




      // totalresidual = fmax(fmax(fmax(pinfeas, dinfeas),diff_fixp),complem);
      totalresidual = fmax(fmax(fmax(pinfeas, dinfeas),diff_fixp),udotr);
      totalresidual_Jor = fmax(fmax(fmax(pinfeas, dinfeas),diff_fixp),complem);


      // // compute Projection Error

      // projerr = projectionError(velocity, reaction, n, tol);
      projerr = projectionError_based_reaction_norm_infinity_conic(H, M, f, w, reaction, n, NOT_ON_DUAL_CONE);

      if (total_time > max_time)
      {
        hasNotConverged = 5;
        numerics_printf_verbose(-1, "| %3i%c| %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e ||",
                              iteration, fws, pinfeas, dinfeas, diff_fixp, udotr, complem, projerr, barr_param, sigma);
        printf("\nREQUEST TIMEOUT\n");
        break;
      }
      else if (alpha_primal <= 1e-13)
      {
        hasNotConverged = 6;
        numerics_printf_verbose(-1, "| %3i%c| %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e ||",
                              iteration, fws, pinfeas, dinfeas, diff_fixp, udotr, complem, projerr, barr_param, sigma);
        printf("\nTiny step-length\n\n");
        break;
      }


      // if ( projerr <= tol && totalresidual_Jor <= tol )
      // if ( projerr <= tol )
      if ( totalresidual_Jor <= tol )
      {
        double unitur;
        for (int i = 0; i < n; i++)
        {
           unitur = cblas_ddot(3, velocity+3*i, 1, reaction+3*i, 1);
           if (unitur<0)
             printf("UR NEGATIF %9.2e\n", unitur);
        }

        numerics_printf_verbose(-1, "| %3i%c| %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e ||",
                              iteration, fws, pinfeas, dinfeas, diff_fixp, udotr, complem, projerr, barr_param, sigma);

        hasNotConverged = 0;
        break;
      }

      t1 = clock();
      // scale_sub_diff = fmax(0.9, alpha_primal);
      scale_sub_diff = 1.;
      if (iteration == 0)
      {
        // First linear linear system
        // Create Jac
        J = NM_create(NM_SPARSE, m + 2*nd + n, m + 2*nd + n);
        NM_triplet_alloc(J, J_nzmax);
        J->matrix2->origin = NSM_TRIPLET;

        // Insert fixed data
        NM_insert(J, M, 0, 0);
        NM_insert(J, minus_H, m + nd, 0);
        NM_insert(J, minus_Ht, 0, m + nd);
        NM_insert(J, eye_nd, m + nd, m);
        NM_insert(J, minus_e, m + nd, m + 2*nd);
        NM_insert(J, eye_n, m + 2*nd, m + 2*nd);

        // /* regularization */
        // // double regul = -1.*fmin(1e-6, 1.*fmax(diff_fixp, udotr));
        // double regul = -1.*1e-8;
        // NM_insert(J, NM_scalar(nd, regul), m + nd, m + nd);


        J_nz_captured = J->matrix2->triplet->nz;  // Save nz at this moment for resetting not fixed data after

        // Compute NOT fixed data
        arrow_r = Arrow_repr(reaction, nd, n);    // Allocation is done in Arrow_repr
        arrow_u = Arrow_repr(velocity, nd, n);

        subdiff_u = NM_create(NM_SPARSE, n, nd);  // Allocation of subdiff_u
        NM_triplet_alloc(subdiff_u, subdiff_u_nzmax);
        NM_fill(subdiff_u, NM_SPARSE, n, nd, subdiff_u->matrix2);

        for(size_t i = 0; i < n; ++i)
        {
          pos = i * d;
          nub = nub_vec[i];

          NM_entry(subdiff_u, i, pos+1, -1.*scale_sub_diff*velocity[pos+1]/nub);
          NM_entry(subdiff_u, i, pos+2, -1.*scale_sub_diff*velocity[pos+2]/nub);

          fixpConstraint[i] = diff_fixp_vec[i];
        }


        // Insert NOT fixed data to Jac
        NM_insert(J, arrow_r, m, m);
        NM_insert(J, subdiff_u, m + 2*nd, m);
        NM_insert(J, arrow_u, m, m + nd);

        J_nz_final = J->matrix2->triplet->nz;
      } // end if iter == 0

      else
      {
        // Clear memory no longer used
        if(arrow_r) { NM_free(arrow_r); arrow_r = NULL; }
        if(arrow_u) { NM_free(arrow_u); arrow_u = NULL; }
        if(subdiff_u) { NM_free(subdiff_u); subdiff_u = NULL; }

        // Clear internal data
        NM_clearCSC(J);
        NM_clearDense(J);
        NM_clearCSCTranspose(J);
        NM_clearCSR(J);
        if(J->matrix2->linearSolverParams)
        {
          NSM_linearSolverParams_free(J->matrix2->linearSolverParams);
          J->matrix2->linearSolverParams = NULL;
        }
        NM_internalData_free(J);
        if (!NM_destructible(J))
        {
          NM_clear(J->destructible);
          J->destructible = J;
        }
        J->matrix2->triplet->nz = J_nz_captured;    // nz goes back the end point of fixed data

        // Replace the matrix content, without reallocation
        NM_insert_Arrow_to_Triplet(J->matrix2->triplet, m, m, reaction, nd, n);
        NM_insert_Arrow_to_Triplet(J->matrix2->triplet, m, m+nd, velocity, nd, n);

        CSparseMatrix *J_cs = J->matrix2->triplet;

        for (size_t i = 0; i < n; i++)
        {
          pos = i * d;
          // nub = cblas_dnrm2(2, velocity+pos+1, 1);
          nub = nub_vec[i];

          for (size_t j = 1; j < d; j++)
          {
            J_cs->x [J_cs->nz] = -1.*scale_sub_diff*velocity[pos+j]/nub;
            J_cs->i [J_cs->nz] = mp2nd + i;
            J_cs->p [J_cs->nz++] = m + pos + j;
          }

          // fixpConstraint[i] = s[i] - nub;  // fixpConstraint = s - |u_bar|
          fixpConstraint[i] = diff_fixp_vec[i];
        }

      }


      jacobian_is_nan = NM_isnan(J);
      if (jacobian_is_nan)
      {
        numerics_printf_verbose(0, "The Jacobian matrix J contains NaN");
        break;
      }

      // Compute rhs
      cblas_dcopy(m, dualConstraint, 1, rhs, 1);
      JA_prod(velocity, reaction, nd, n, complemConstraint);
      cblas_dcopy(nd, complemConstraint, 1, rhs+m, 1);
      cblas_dcopy(nd, primalConstraint, 1, rhs+m+nd, 1);
      cblas_dcopy(n, fixpConstraint, 1, rhs+m+2*nd, 1);

      cblas_dscal(m + 2*nd + n, -1.0, rhs, 1);
      cblas_dcopy(m + 2*nd + n, rhs, 1, rhs_2, 1);    // rhs_2 = old rhs (already had minus sign)


      if (cblas_dnrm2(n, fixpConstraint, 1) > complemResidualNorm(velocity, reaction, nd, n))
      {
        sigma = 0.499;
      }
      else
      {
        // Solve
        rhs_tmp = (double*)calloc(m+2*nd+n,sizeof(double));
        cblas_dcopy(m+2*nd+n, rhs_2, 1, rhs_tmp, 1);
        for (int k=0; k<m+2*nd+n; sol[k] = 0., k++);  // reset sol

        max_refine = 1;
        if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] == SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES) max_refine = 10;

        for (int itr = 0; itr < max_refine; itr++)
        {
          NM_LU_solve(J, rhs_tmp, 1);
          cblas_daxpy(m+2*nd+n, 1.0, rhs_tmp, 1, sol, 1);
          cblas_dcopy(m+2*nd+n, rhs_2, 1, rhs_tmp, 1);
          NM_gemv(-1.0, J, sol, 1.0, rhs_tmp);
          //printf("refinement iterations = %i %8.2e\n",itr, cblas_dnrm2(m+2*nd, rhs_tmp, 1));
          if (cblas_dnrm2(m+2*nd+n, rhs_tmp, 1) <= 1e-14)
          {
            // printf("\nrefinement iterations = %d %8.2e\n",itr+1, cblas_dnrm2(m+2*nd+n, rhs_tmp, 1));
            free(rhs_tmp);
            break;
          }
        }

        // Compute direction
        cblas_dcopy(m, sol, 1, d_globalVelocity, 1);
        cblas_dcopy(nd, sol+m, 1, d_velocity, 1);
        cblas_dcopy(nd, sol+m+nd, 1, d_reaction, 1);
        cblas_dcopy(n, sol+m+2*nd, 1, d_s, 1);

        // Compute step-length
        alpha_primal = getStepLength(velocity, d_velocity, nd, n, 0.999);
        alpha_dual = getStepLength(reaction, d_reaction, nd, n, 0.999);
        alpha_primal = alpha_dual = fmin(alpha_primal, alpha_dual);

        // gmm = gmmp1 + gmmp2 * alpha_primal;

        /* ----- Corrector step of Mehrotra ----- */
        cblas_dcopy(nd, velocity, 1, u_plus_du, 1);
        cblas_dcopy(nd, reaction, 1, r_plus_dr, 1);
        cblas_daxpy(nd, alpha_primal, d_velocity, 1, u_plus_du, 1);
        cblas_daxpy(nd, alpha_dual, d_reaction, 1, r_plus_dr, 1);

        barr_param_a = cblas_ddot(nd, u_plus_du, 1, r_plus_dr, 1) / n ;
        e = barr_param > sgmp1 ? fmax(1.0, sgmp2 * pow(alpha_primal,2)) : sgmp3;
        // sigma = 0.49*fmin(1.0, pow(barr_param_a / barr_param, e));
        sigma = 0.499*fmin(1., pow(barr_param_a / barr_param, e));

        // Compute rhs for 2nd linear system
        cblas_dcopy(m+2*nd+n, rhs_2, 1, rhs, 1);    // Get back value for rhs


        if (sigma < 0.25)
        {
          JA_prod(d_velocity, d_reaction, nd, n, dudr_jprod);
          cblas_daxpy(nd, -1., dudr_jprod, 1, rhs + m, 1);
        }
      }


      for (int k = 0; k < nd; rhs[m+k] += 2*sigma*barr_param, k+=d);
      cblas_dcopy(m + 2*nd + n, rhs, 1, rhs_2, 1);    // rhs_2 = old rhs

      // SOLVE
      // NM_LU_solve(J, rhs, 1);


      rhs_tmp = (double*)calloc(m+2*nd+n,sizeof(double));
      cblas_dcopy(m+2*nd+n, rhs_2, 1, rhs_tmp, 1);
      for (int k=0; k<m+2*nd+n; sol[k] = 0., k++);  // reset sol

      max_refine = 1;
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] == SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES) max_refine = 10;

      for (int itr = 0; itr < max_refine; itr++)
      {
        NM_LU_solve(J, rhs_tmp, 1);
        cblas_daxpy(m+2*nd+n, 1.0, rhs_tmp, 1, sol, 1);
        cblas_dcopy(m+2*nd+n, rhs_2, 1, rhs_tmp, 1);
        NM_gemv(-1.0, J, sol, 1.0, rhs_tmp);
        //printf("refinement iterations = %i %8.2e\n",itr, cblas_dnrm2(m+2*nd, rhs_tmp, 1));
        if (cblas_dnrm2(m+2*nd+n, rhs_tmp, 1) <= 1e-14)
        {
          // printf("\nrefinement iterations = %d %8.2e\n",itr+1, cblas_dnrm2(m+2*nd+n, rhs_tmp, 1));
          // free(rhs_tmp);
          break;
        }
      }

      NM_gemv(1.0, J, sol, -1.0, rhs_2);

      LS_norm_d = cblas_dnrm2(m, rhs_2, 1);
      LS_norm_c = cblas_dnrm2(nd, rhs_2+m, 1);
      LS_norm_p = cblas_dnrm2(nd, rhs_2+m+nd, 1);
      LS_norm_f = cblas_dnrm2(n, rhs_2+m+2*nd, 1);

      cblas_dcopy(m, sol, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, sol+m, 1, d_velocity, 1);
      cblas_dcopy(nd, sol+m+nd, 1, d_reaction, 1);
      cblas_dcopy(n, sol+m+2*nd, 1, d_s, 1);

      /* computing the step-length */
      alpha_primal = getStepLength(velocity, d_velocity, nd, n, gmm);
      alpha_dual = getStepLength(reaction, d_reaction, nd, n, gmm);
      alpha_primal = alpha_dual = fmin(alpha_primal, alpha_dual);
      gmm = gmmp1 + gmmp2 * alpha_primal;

      cblas_daxpy(m, alpha_primal, d_globalVelocity, 1, globalVelocity, 1);
      cblas_daxpy(nd, alpha_primal, d_velocity, 1, velocity, 1);
      cblas_daxpy(nd, alpha_primal, d_reaction, 1, reaction, 1);
      cblas_daxpy(n, alpha_primal, d_s, 1, s, 1);

      t2 = clock();
      total_time += (double)(t2-t1)/(double)clk_tck;

      if (NV_isnan(globalVelocity, m) | NV_isnan(velocity, nd) | NV_isnan(reaction, nd) | NV_isnan(s, n))
      {
        hasNotConverged = 2;
        break;
      }

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
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_QP2:
    {
      // Stopping test using norm type
      primalResidual_s_type(velocity, H, globalVelocity, w, s, primalConstraint, &pinfeas, tol, type);
      dualResidual_type(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol, type);
      // udotr = xdoty_type(n, nd, velocity, reaction, type);
      udotr = xdoty_type(n, nd, velocity, reaction, NORM_2_INF);
      complem = complemResidualNorm(velocity, reaction, nd, n);
      barr_param = cblas_ddot(nd, reaction, 1, velocity, 1) / n;

      diff_fixp = 0.;
      for (unsigned int i = 0; i<nd; i+=d)
      {
        nub = cblas_dnrm2(2, velocity+i+1, 1);
        diff_fixp_vec[i/d] = fabs(s[i/d]-nub);
      }

      // Stopping test using norm type
      if (type == NORM_2)
        diff_fixp = cblas_dnrm2(n, diff_fixp_vec, 1);
      else if (type == NORM_INF)
      {
        int maxIndex = cblas_idamax(n, diff_fixp_vec, 1);
        diff_fixp = fabs(diff_fixp_vec[maxIndex]);
      }
      else
      {
        fprintf(stderr, "type = %d is undefined.\n", type);
        exit(EXIT_FAILURE);
      }

      // totalresidual = fmax(fmax(fmax(pinfeas, dinfeas),diff_fixp),complem);
      totalresidual = fmax(fmax(fmax(pinfeas, dinfeas),diff_fixp),udotr);

      // compute Projection Error
      // projerr = projectionError(velocity, reaction, n, tol);
      projerr = projectionError_based_reaction_norm_infinity_conic(H, M, f, w, reaction, n, NOT_ON_DUAL_CONE);



      if ( totalresidual <= tol )
      {
        double unitur;
        for (int i = 0; i < n; i++)
        {
           unitur = cblas_ddot(3, velocity+3*i, 1, reaction+3*i, 1);
           if (unitur<0)
             printf("UR NEGATIF %9.2e\n", unitur);
        }

        if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
        {
          fprintf(iterates,"%d %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e;\n",
                iteration, pinfeas, dinfeas, diff_fixp, uor_mu, udotr, projerr, barr_param, alpha_primal,
                fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
                fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
                fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]),
                fabs(d_s[cblas_idamax(n, d_s, 1)]),
                LS_norm_p, LS_norm_d, LS_norm_c, LS_norm_f);
        }

        numerics_printf_verbose(-1, "| %3i%c| %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e ||",
                              iteration, fws, pinfeas, dinfeas, diff_fixp, udotr, complem, projerr, barr_param, sigma);

        hasNotConverged = 0;
        break;
      }


      // First linear linear system
      J = NM_create(NM_SPARSE, m + 2*nd + n, m + 2*nd + n);
      J_nzmax = M_nzmax + H_nzmax + 2*(d*3-2)*n + H_nzmax + nd + n + 2*n + n;
      NM_triplet_alloc(J, J_nzmax);
      J->matrix2->origin = NSM_TRIPLET;

      // arrow_r = Arrow_repr(reaction, nd, n);
      // arrow_u = Arrow_repr(velocity, nd, n) ;

      // Create subdiff_u
      subdiff_u = NM_create(NM_SPARSE, n, nd);
      size_t subdiff_u_nzmax = 2*n;
      NM_triplet_alloc(subdiff_u, subdiff_u_nzmax);
      NM_fill(subdiff_u, NM_SPARSE, n, nd, subdiff_u->matrix2);

      /* Matrix filling */
      scale_sub_diff = fmax(0.9, alpha_primal);
      for(size_t i = 0; i < n; ++i)
      {
        pos = i * d;
        nub = cblas_dnrm2(2, velocity+pos+1, 1);

        NM_entry(subdiff_u, i, pos+1, -1.*scale_sub_diff*velocity[pos+1]/nub);
        NM_entry(subdiff_u, i, pos+2, -1.*scale_sub_diff*velocity[pos+2]/nub);

        fixpConstraint[i] = s[i] - nub;  // fixpConstraint = s - |u_bar|
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

      // /* regularization */
      // // double regul = -1.*fmin(1e-6, 1.*fmax(diff_fixp, udotr));
      // double regul = -1.*1e-6;
      // NM_insert(J, NM_scalar(nd, regul), m + nd, m + nd);

      if(subdiff_u) { NM_free(subdiff_u); subdiff_u = NULL; }
      if(Qp2) {NM_free(Qp2); Qp2 = NULL;}


      jacobian_is_nan = NM_isnan(J);
      if (jacobian_is_nan)
      {
        numerics_printf_verbose(0, "The Jacobian matrix J contains NaN");
        break;
      }


      // Compute rhs
      cblas_dcopy(m, dualConstraint, 1, rhs, 1);
      cblas_dcopy(nd, reaction, 1, rhs+m, 1);
      cblas_dcopy(nd, primalConstraint, 1, rhs+m+nd, 1);
      cblas_dcopy(n, fixpConstraint, 1, rhs+m+2*nd, 1);

      cblas_dscal(m + 2*nd + n, -1.0, rhs, 1);
      cblas_dcopy(m + 2*nd + n, rhs, 1, rhs_2, 1);    // rhs_2 = old rhs


      // Solve
      rhs_tmp = (double*)calloc(m+2*nd+n,sizeof(double));
      cblas_dcopy(m+2*nd+n, rhs_2, 1, rhs_tmp, 1);
      for (int k=0; k<m+2*nd+n; sol[k] = 0., k++);  // reset sol

      max_refine = 1;
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] == SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES) max_refine = 10;

      for (int itr = 0; itr < max_refine; itr++)
      {
        NM_LU_solve(J, rhs_tmp, 1);
        cblas_daxpy(m+2*nd+n, 1.0, rhs_tmp, 1, sol, 1);
        cblas_dcopy(m+2*nd+n, rhs_2, 1, rhs_tmp, 1);
        NM_gemv(-1.0, J, sol, 1.0, rhs_tmp);
        //printf("refinement iterations = %i %8.2e\n",itr, cblas_dnrm2(m+2*nd, rhs_tmp, 1));
        if (cblas_dnrm2(m+2*nd+n, rhs_tmp, 1) <= 1e-14)
        {
          // printf("\nrefinement iterations = %d %8.2e\n",itr+1, cblas_dnrm2(m+2*nd+n, rhs_tmp, 1));
          free(rhs_tmp);
          break;
        }
      }

      // Compute direction
      cblas_dcopy(m, sol, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, sol+m, 1, d_velocity, 1);
      cblas_dcopy(nd, sol+m+nd, 1, d_reaction, 1);
      cblas_dcopy(n, sol+m+2*nd, 1, d_s, 1);

      // Compute step-length
      alpha_primal = getStepLength(velocity, d_velocity, nd, n, 0.999);
      alpha_dual = getStepLength(reaction, d_reaction, nd, n, 0.999);
      alpha_primal = alpha_dual = fmin(alpha_primal, alpha_dual);

      gmm = gmmp1 + gmmp2 * alpha_primal;



      /* ----- Corrector step of Mehrotra ----- */
      cblas_dcopy(nd, velocity, 1, u_plus_du, 1);
      cblas_dcopy(nd, reaction, 1, r_plus_dr, 1);
      cblas_daxpy(nd, alpha_primal, d_velocity, 1, u_plus_du, 1);
      cblas_daxpy(nd, alpha_dual, d_reaction, 1, r_plus_dr, 1);

      barr_param_a = cblas_ddot(nd, u_plus_du, 1, r_plus_dr, 1) / n ;
      e = barr_param > sgmp1 ? fmax(1.0, sgmp2 * pow(alpha_primal,2)) : sgmp3;
      sigma = 0.5*fmin(1.0, pow(barr_param_a / barr_param, e));

      // if (diff_fixp > 10.*udotr)
      if (cblas_dnrm2(n, diff_fixp_vec, 1) > 1.*complem)
      {
        sigma = 0.5;
      }



      // Compute rhs for 2nd linear system
      cblas_dcopy(m+2*nd+n, rhs_2, 1, rhs, 1);    // Get back value for rhs (already had minus sign inside)


      QNTpz(velocity, reaction, d_velocity, nd, n, Qp_du);          // Qp_du = du_affine_hat
      QNTpinvz(velocity, reaction, d_reaction, nd, n, Qpinv_dr);    // Qpinv_dr = dr_affine_check
      JA_prod(Qp_du, Qpinv_dr, nd, n, dudr_jprod);                    // dudr_jprod = Qp_du o Qpinv_dr
      QNTpz(velocity, reaction, velocity, nd, n, Qp_u);             // Qp_u = u_hat
      Jxinvprody(Qp_u, dudr_jprod, nd, n, ududr);                     // ududr = u_hat_inv o (Qp_du o Qpinv_dr)
      QNTpz(velocity, reaction, ududr, nd, n, Qp_ududr);             // Qp_ududr = Qp * [ u_hat_inv o (Qp_du o Qpinv_dr) ]

      JA_inv(velocity, nd, n, velocity_inv);                    // velocity_inv = u^-1
      cblas_dscal(nd, 2.*sigma*barr_param, velocity_inv, 1);    // velocity_inv = 2 * sigma * mu * u^-1

      cblas_daxpy(nd, -1.0, Qp_ududr, 1, rhs+m, 1);
      cblas_daxpy(nd, 1.0, velocity_inv, 1, rhs+m, 1);

      cblas_dcopy(m + 2*nd + n, rhs, 1, rhs_2, 1);    // rhs_2 = old rhs


      // SOLVE
      // NM_LU_solve(J, rhs, 1);


      rhs_tmp = (double*)calloc(m+2*nd+n,sizeof(double));
      cblas_dcopy(m+2*nd+n, rhs_2, 1, rhs_tmp, 1);
      for (int k=0; k<m+2*nd+n; sol[k] = 0., k++);  // reset sol

      max_refine = 1;
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] == SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES) max_refine = 10;

      for (int itr = 0; itr < max_refine; itr++)
      {
        NM_LU_solve(J, rhs_tmp, 1);
        cblas_daxpy(m+2*nd+n, 1.0, rhs_tmp, 1, sol, 1);
        cblas_dcopy(m+2*nd+n, rhs_2, 1, rhs_tmp, 1);
        NM_gemv(-1.0, J, sol, 1.0, rhs_tmp);
        //printf("refinement iterations = %i %8.2e\n",itr, cblas_dnrm2(m+2*nd, rhs_tmp, 1));
        if (cblas_dnrm2(m+2*nd+n, rhs_tmp, 1) <= 1e-14)
        {
          // printf("\nrefinement iterations = %d %8.2e\n",itr+1, cblas_dnrm2(m+2*nd+n, rhs_tmp, 1));
          free(rhs_tmp);
          break;
        }
      }

      NM_gemv(1.0, J, sol, -1.0, rhs_2);

      LS_norm_d = cblas_dnrm2(m, rhs_2, 1);
      LS_norm_c = cblas_dnrm2(nd, rhs_2+m, 1);
      LS_norm_p = cblas_dnrm2(nd, rhs_2+m+nd, 1);
      LS_norm_f = cblas_dnrm2(n, rhs_2+m+2*nd, 1);

      cblas_dcopy(m, sol, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, sol+m, 1, d_velocity, 1);
      cblas_dcopy(nd, sol+m+nd, 1, d_reaction, 1);
      cblas_dcopy(n, sol+m+2*nd, 1, d_s, 1);

      /* computing the step-length */
      alpha_primal = getStepLength(velocity, d_velocity, nd, n, gmm);
      alpha_dual = getStepLength(reaction, d_reaction, nd, n, gmm);
      alpha_primal = alpha_dual = fmin(alpha_primal, alpha_dual);

      cblas_daxpy(m, alpha_primal, d_globalVelocity, 1, globalVelocity, 1);
      cblas_daxpy(nd, alpha_primal, d_velocity, 1, velocity, 1);
      cblas_daxpy(nd, alpha_primal, d_reaction, 1, reaction, 1);
      cblas_daxpy(n, alpha_primal, d_s, 1, s, 1);

      if (NV_isnan(globalVelocity, m) | NV_isnan(velocity, nd) | NV_isnan(reaction, nd) | NV_isnan(s, n))
      {
        hasNotConverged = 2;
        break;
      }

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


    if (NV_isnan(globalVelocity, m) | NV_isnan(velocity, nd) | NV_isnan(reaction, nd) | NV_isnan(s, n))
    {
      hasNotConverged = 2;
      break;
    }

    if (hasNotConverged == 0 || hasNotConverged >= 2)
      break;

    numerics_printf_verbose(-1, "| %3i%c| %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e || %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |",
                          iteration, fws, pinfeas, dinfeas, diff_fixp, udotr, complem, projerr, barr_param, sigma, alpha_primal,
                          fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
                          fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
                          fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]),
                          fabs(d_s[cblas_idamax(n, d_s, 1)]),
           LS_norm_p, LS_norm_d, LS_norm_c, LS_norm_f);


    classify_BNRT(velocity, reaction, nd, n, &nB, &nN, &nR, &nT);
    // printf("B/N/R/T = %d %d %d %d\n", nB, nN, nR, nT);


    if(J_dense) { J_dense = NM_free(J_dense); J_dense = NULL;}

    iteration++;

  } // while loop
  t2 = clock();




  /* Checking strict complementarity */
  classify_BNRT(velocity, reaction, nd, n, &nB, &nN, &nR, &nT);
  if (nT > 0)
    printf("Ratio of Strict complementarity solutions: %4i / %4i = %4.2f \t %4i %4i %4i\n", n-nT, n, (double)(n-nT)/(double)n, nB, nN, nR);
  else
    printf("Strict complementarity satisfied: %4i %4i %4i\n", nB, nN, nR);




  // Store solution into file
  if (save_sol_point)
  {
    sol_file = fopen("sol_data.res", "w");
    // store v
    for (int i=0; i < m; i++)
    {
      fprintf(sol_file, "%8.20e ", globalVelocity[i]);
    }
    fprintf(sol_file, "\n");

    // store u
    for (int i=0; i < nd; i++)
    {
      fprintf(sol_file, "%8.20e ", velocity[i]);
    }
    fprintf(sol_file, "\n");

    // store r
    for (int i=0; i < nd; i++)
    {
      fprintf(sol_file, "%8.20e ", reaction[i]);
    }
    fprintf(sol_file, "\n");

    // store s
    for (int i=0; i < n; i++)
    {
      fprintf(sol_file, "%8.20e ", s[i]);
    }
    fprintf(sol_file, "\n");

    fclose(sol_file);
  }



  /* ----- return to original variables ------ */
  NM_gemv(1.0, P_mu_inv, velocity, 0.0, data->tmp_point->t_velocity);
  cblas_dcopy(nd, data->tmp_point->t_velocity, 1, velocity, 1);

  NM_gemv(1.0, P_mu, reaction, 0.0, data->tmp_point->t_reaction);
  cblas_dcopy(nd, data->tmp_point->t_reaction, 1, reaction, 1);

  options->dparam[SICONOS_DPARAM_RESIDU] = totalresidual; //NV_max(error, 4);
  // options->iparam[SICONOS_IPARAM_ITER_DONE] = iteration;


  if(internal_allocation)
  {
    gfc3d_IPM_SNM_free(problem,options);
  }

  options->solverData = (double *)malloc(sizeof(double));
  double *projerr_ptr = (double *)options->solverData;
  *projerr_ptr = projerr;


  if(H_tilde) H_tilde = NM_free(H_tilde);
  if(minus_H) minus_H = NM_free(minus_H);
  if(H) H = NM_free(H);
  if(minus_M) minus_M = NM_free(minus_M);
  if(minus_Ht) minus_Ht = NM_free(minus_Ht);
  if(eye_nd) eye_nd = NM_free(eye_nd);
  if(eye_n) eye_n = NM_free(eye_n);
  if(Ht) Ht = NM_free(Ht);
  if(subdiff_u) subdiff_u = NM_free(subdiff_u);
  if(minus_e) { minus_e = NM_free(minus_e); minus_e = NULL; }
  if(arrow_r) { NM_free(arrow_r); arrow_r = NULL; }
  if(arrow_u) { NM_free(arrow_u); arrow_u = NULL; }
  if(subdiff_u) { NM_free(subdiff_u); subdiff_u = NULL; }
  if(J) { J = NM_free(J); J = NULL; }
  if (blk_num_name) {free(blk_num_name); blk_num_name = NULL;}

  if (rhs_tmp) {free(rhs_tmp); rhs_tmp = NULL;}
  if (nub_vec) {free(nub_vec); nub_vec = NULL;}
  if (setR) {free(setR); setR = NULL;}

  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
  {
    // fprintf(iterates, "];\n\n");
    fclose(iterates);
    if (iterates_2) fclose(iterates_2);
  }

  if (rhs_2) free(rhs_2);
  if (sol) free(sol);

  if (d_globalVelocity) free(d_globalVelocity);
  if (d_velocity) free(d_velocity);
  if (d_reaction) free(d_reaction);
  if (s) free(s);
  if (d_s) free(d_s);
  if (w_ori) free(w_ori);
  if (fixpConstraint) free(fixpConstraint);
  if (velocity_inv) free(velocity_inv);
  if (arr_norm_theta) free(arr_norm_theta);
  if (str) {free(str); str = NULL;}

  *info = hasNotConverged;



}

void gfc3d_ipm_snm_set_default(SolverOptions* options)
{
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO] = SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_NO;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] = SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_NOSCAL;
  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] = SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_QP2;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE] = 0;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] = SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES;

  options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000000;
  options->dparam[SICONOS_DPARAM_TOL] = 1e-14;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1] = 1e-10;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2] = 3.;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3] = 1.;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1] = 0.9;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2] = 0.09; //0.095

}
