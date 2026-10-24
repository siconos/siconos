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

#include <assert.h>
#include <float.h>
#include <limits.h>  // For LLONG_MAX
#include <math.h>
#include <stdarg.h>  // for va_list, va_start, va_end
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "FrictionContact_options.h"
#include "JordanAlgebra.h"
#include "NumericsMatrix.h"
#include "NumericsMatrix_internal.h"
#include "NumericsSparseMatrix.h"
#include "NumericsVector.h"
#include "SiconosBlas.h"
#include "fc3d_short_names.h"
#include "float.h"
#include "gfc3d_Solvers.h"
#include "gfc3d_ipm.h"
#include "gfc3d_ipm_internal.h"
#include "io_tools.h"
#include "numerics_errors.h"
#include "numerics_verbose.h"
#include "safe_casts.h"
#include "solver_registry.h"
/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "siconos_debug.h"

#if defined(WITH_FCLIB) && defined(WITH_HDF5)
#include <fclib.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

#define NUM_BLOCKS 0
#define CONTACT_INDEX 1
#define BODY_INDEX 2
#define RANK_Hc 3

/* ------------------------- Helper functions implementation ------------------------------ */
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
static void printDataProbMatlabFile(NumericsMatrix* M, double* f, NumericsMatrix* H, double* w,
                                    int d, int n, int m, double* mu, FILE* file) {
  fprintf(file, "d = %3i;\n", d);
  fprintf(file, "n = %6i;\n", n);
  fprintf(file, "m = %6i;\n", m);

  fprintf(file, "M = [\n");
  CSparseMatrix_print_in_Matlab_file(NM_triplet(M), 0, file);
  fprintf(file, "];\n");
  fprintf(file, "M = sparse(int32(M(:,1)), int32(M(:,2)), M(:,3));\n");

  fprintf(file, "H = [\n");
  CSparseMatrix_print_in_Matlab_file(NM_triplet(H), 0, file);
  fprintf(file, "];\n");
  fprintf(file, "H = sparse(int32(H(:,1)), int32(H(:,2)), H(:,3));\n");

  fprintf(file, "f = [");
  for (int i = 0; i < m; i++) {
    fprintf(file, "%8.20e; ", f[i]);
  }
  fprintf(file, "];\n");

  fprintf(file, "w = [");
  for (int i = 0; i < n * d; i++) {
    fprintf(file, "%8.20e; ", w[i]);
  }
  fprintf(file, "];\n");

  fprintf(file, "mu = [");
  for (int i = 0; i < n; i++) {
    fprintf(file, "%8.20e; ", mu[i]);
  }
  fprintf(file, "];\n");
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

// static void printIteresProbMatlabFile(int iteration, double* v, double* u, double* r,
//                                       double* s, double* dv, double* du, double* dr,
//                                       double* ds, int d, int n, int m, FILE* file)
// // static void printIteresProbMatlabFile(int iteration, double pinfeas, double dinfeas,
// double
// // udotr, double smub, int d, int n, int m, FILE * file)
// {
//   // fprintf(file,"v(%3i,:) = [",iteration+1);
//   // for(int i = 0; i < m; i++)
//   // {
//   //   fprintf(file, "%8.20e, ", v[i]);
//   // }
//   // fprintf(file,"];\n");

//   fprintf(file, "u(%3i,:) = [", iteration + 1);
//   for (int i = 0; i < n * d; i++) {
//     fprintf(file, "%8.20e, ", u[i]);
//   }
//   fprintf(file, "];\n");

//   fprintf(file, "r(%3i,:) = [", iteration + 1);
//   for (int i = 0; i < n * d; i++) {
//     fprintf(file, "%8.20e, ", r[i]);
//   }
//   fprintf(file, "];\n");

//   fprintf(file, "s(%3i,:) = [", iteration + 1);
//   for (int i = 0; i < n; i++) {
//     fprintf(file, "%8.20e, ", s[i]);
//   }
//   fprintf(file, "];\n");

//   // fprintf(file,"dv(%3i,:) = [",iteration+1);
//   // for(int i = 0; i < m; i++)
//   // {
//   //   fprintf(file, "%8.20e, ", dv[i]);
//   // }
//   // fprintf(file,"];\n");

//   fprintf(file, "du(%3i,:) = [", iteration + 1);
//   for (int i = 0; i < n * d; i++) {
//     fprintf(file, "%8.20e, ", du[i]);
//   }
//   fprintf(file, "];\n");

//   fprintf(file, "dr(%3i,:) = [", iteration + 1);
//   for (int i = 0; i < n * d; i++) {
//     fprintf(file, "%8.20e, ", dr[i]);
//   }
//   fprintf(file, "];\n");

//   // fprintf(file,"ds(%3i,:) = [",iteration+1);
//   // for(int i = 0; i < n; i++)
//   // {
//   //   fprintf(file, "%8.20e, ", ds[i]);
//   // }
//   // fprintf(file,"];\n");

//   // fprintf(file,"pinfeas(%3i) = %20.16e;\n",iteration+1,pinfeas);
//   // fprintf(file,"dinfeas(%3i) = %20.16e;\n",iteration+1,dinfeas);
//   // fprintf(file,"udotr(%3i) = %20.16e;\n",iteration+1,udotr);
//   // // fprintf(file,"residu(%3i) = %20.16e;\n",iteration+1,totalresidual);
//   // fprintf(file,"smub(%3i) = %20.16e;\n",iteration+1,smub);

//   return;
// }

// static void printVectorMatlabFile(int iteration, double* vec, int vecSize, FILE* file) {
//   fprintf(file, "vector(%4i,:) = [", iteration + 1);
//   for (int i = 0; i < vecSize; i++) {
//     fprintf(file, "%24.16e, ", vec[i]);
//   }
//   fprintf(file, "];\n");
//   return;
// }

// static float randomFloat(float min, float max) {
//   return min + (float)rand() / ((float)RAND_MAX / (max - min));
// }

// [OLD FUNCTION] This function is used to read blocks info stored each FCLIB problem test
// Need the FCLIB dataset with blocks info stored inside
int* read_fricprob_block(const char* path, int type, int blk_index) {
  int* out = NULL;

  if (type < NUM_BLOCKS || type > RANK_Hc) {
    fprintf(stderr, "ERROR: out of \"type\"\n");
    return NULL;
  }

  int is_hdf5 = check_hdf5_file(path);
  if (is_hdf5) {
#if defined(WITH_FCLIB) && defined(WITH_HDF5)
    hid_t file_id, group_id, subgroup_id, dataset_id, dataspace_id;
    hssize_t num_elements;
    if ((file_id = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0) {
      fprintf(stderr, "ERROR: opening file failed\n");
      return NULL;
    }

    // Get number of blocks
    int numBlk = 0;
    group_id = H5Gopen(file_id, "/fclib_global/blocks", H5P_DEFAULT);
    H5LTread_dataset_int(file_id, "/fclib_global/blocks/N", &numBlk);

    if (type == NUM_BLOCKS) {
      out = (int*)calloc(1, sizeof(int));
      *out = numBlk;
    }

    else {
      if (blk_index < 0) {
        fprintf(stderr, "ERROR: block index must not be negative\n");
        return NULL;
      } else if (blk_index >= numBlk) {
        fprintf(stderr, "ERROR: block index must be in range [0, %d)\n", numBlk);
        return NULL;
      }

      char block_X[20] = "block_";
      char dest_path[50] = "";
      sprintf(block_X + strlen(block_X), "%d", blk_index);
      sprintf(dest_path + strlen(dest_path), "%s%s", "/fclib_global/blocks/", block_X);
      subgroup_id = H5Gopen(file_id, dest_path, H5P_DEFAULT);

      if (type == CONTACT_INDEX) {
        sprintf(dest_path + strlen(dest_path), "%s", "/contact");
        dataset_id = H5Dopen1(file_id, dest_path);
        dataspace_id = H5Dget_space(dataset_id);
        num_elements = H5Sget_simple_extent_npoints(
            dataspace_id);  // Get the number of elements of contact vector

        out = (int*)calloc(num_elements + 1, sizeof(int));
        out[0] = (int)num_elements;  // 1st element of returned vector is the number of data
                                     // elements that the vector contains
        H5LTread_dataset_int(file_id, dest_path, out + 1);  // Get data
      }

      else if (type == BODY_INDEX) {
        sprintf(dest_path + strlen(dest_path), "%s", "/body");
        dataset_id = H5Dopen1(file_id, dest_path);
        dataspace_id = H5Dget_space(dataset_id);
        num_elements = H5Sget_simple_extent_npoints(dataspace_id);

        out = (int*)calloc(num_elements + 1, sizeof(int));
        out[0] = (int)num_elements;
        H5LTread_dataset_int(file_id, dest_path, out + 1);
      }

      else if (type == RANK_Hc) {
        sprintf(dest_path + strlen(dest_path), "%s", "/rank_H_blk");
        out = (int*)calloc(1, sizeof(int));
        H5LTread_dataset_int(file_id, dest_path, out);
      }

      H5Gclose(subgroup_id);
    }

    H5Gclose(group_id);
    H5Fclose(file_id);

#else
    numerics_error_log(
        "gfc3d_IPM_SNM",
        "Try to read an hdf5 file, while fclib interface is not active. Recompile "
        "Siconos with fclib.",
        path);
    return NULL;
#endif
  } else {
    numerics_error_log("gfc3d_IPM_SNM", "Not a hdf5 file ", path);
    return NULL;
  }
  return out;
}

static int NM_insert_Arrow_to_Triplet(CSparseMatrix* triplet, const size_t start_i,
                                      const size_t start_j, const double* const vec,
                                      const size_t vecSize, const size_t varsCount) {
  size_t dimension = (size_t)(vecSize / varsCount);
  size_t pos;
  size_t total_element = (dimension * 3 - 2) * varsCount;
  if (!(total_element <= LLONG_MAX))  // To detect bad conversion with nzmax type
    return numerics_error("NM_insert_Arrow_to_Triplet", "value too large for an int64_t");

  if (triplet->nzmax < ((int64_t)total_element + triplet->nz)) {
    return numerics_error("NM_insert_Arrow_to_Triplet",
                          " Size of allocated triplet memory is not sufficient.");
  }

  for (size_t i = 0; i < varsCount; ++i) {
    pos = i * dimension;

    triplet->x[triplet->nz] = vec[pos];
    triplet->i[triplet->nz] = to_csint(start_i + pos);
    triplet->p[triplet->nz++] = to_csint(start_j + pos);

    for (size_t j = 1; j < dimension; ++j) {
      triplet->x[triplet->nz] = vec[pos + j];
      triplet->i[triplet->nz] = to_csint(start_i + pos);
      triplet->p[triplet->nz++] = to_csint(start_j + pos + j);

      triplet->x[triplet->nz] = vec[pos + j];
      triplet->i[triplet->nz] = to_csint(start_i + pos + j);
      triplet->p[triplet->nz++] = to_csint(start_j + pos);

      triplet->x[triplet->nz] = vec[pos];
      triplet->i[triplet->nz] = to_csint(start_i + pos + j);
      triplet->p[triplet->nz++] = to_csint(start_j + pos + j);
    }
  }
  return 0;
}

/* --------------------------- Interior-point method implementation
 * ------------------------------ */

/* Prepare the problem matrices: symmetrize M and, if requested, convert M and H to sparse
 * storage. Outputs M (may point to problem->M or a fresh copy), the original
 * problem->M if it was symmetrized (NULL otherwise), and H_tilde (a transpose or a
 * sparse copy of the transpose). */
static void gfc3d_IPM_SNM_prepare_problem(GlobalFrictionContactProblem* problem,
                                          SolverOptions* options, size_t m, size_t nd,
                                          NumericsMatrix** M_out,
                                          NumericsMatrix** M_original_out,
                                          NumericsMatrix** H_tilde_out) {
  /* symmetrization of the matrix M */
  *M_original_out = NULL;
  if (!(NM_is_symmetric(problem->M))) {
    printf("#################### SYMMETRIZATION ####################\n");
    NumericsMatrix* MT = NM_transpose(problem->M);
    NumericsMatrix* MSym = NM_add(1 / 2., problem->M, 1 / 2., MT);
    *M_original_out = problem->M; /* keep the user-provided matrix for restoration */
    problem->M = MSym;
    MT = NM_free(MT);
  }

  /* if SICONOS_FRICTION_3D_FORCED_SPARSE_STORAGE is set, force a copy into NM_SPARSE
   * storageType */
  DEBUG_PRINTF("problem->M->storageType : %i\n", problem->M->storageType);
  if (options->iparam[SICONOS_FRICTION_3D_IPARAM_SPARSE_STORAGE] ==
          SICONOS_FRICTION_3D_FORCED_SPARSE_STORAGE &&
      problem->M->storageType == NM_SPARSE_BLOCK) {
    DEBUG_PRINT("Force a copy to sparse storage type\n");
    printf(
        "\n\n\n######################### FORCE SPARSE STORAGE "
        "#########################\n\n\n");
    *M_out = NM_copy_to_sparse(problem->M, DBL_EPSILON);
  } else {
    *M_out = problem->M;
  }

  DEBUG_PRINTF("problem->M->storageType : %i\n", problem->H->storageType);
  if (options->iparam[SICONOS_FRICTION_3D_IPARAM_SPARSE_STORAGE] ==
          SICONOS_FRICTION_3D_FORCED_SPARSE_STORAGE &&
      problem->H->storageType == NM_SPARSE_BLOCK) {
    DEBUG_PRINT("Force a copy to sparse storage type\n");
    *H_tilde_out = NM_copy_to_sparse(NM_transpose(problem->H), DBL_EPSILON);
  } else {
    *H_tilde_out = NM_transpose(problem->H);
  }

  numerics_printf_verbose(1, "nnz M = %8zu density = %9.4f, \t size = (%i, %i)",
                          NM_nnz(problem->M),
                          NM_nnz(problem->M) / 1.0 / problem->M->size0 / problem->M->size1,
                          problem->M->size0, problem->M->size1);
  numerics_printf_verbose(1, "nnz H = %8zu density = %9.4f, \t mu = %.2e", NM_nnz(problem->H),
                          NM_nnz(problem->H) / 1.0 / nd / m, problem->mu[0]);
}

/* Apply the friction-coefficient scaling (P_mu) to H_tilde and w, then build the derived
 * matrices needed by the IPM iterations (-H, -H^T, H^T, -M). */
static void gfc3d_IPM_SNM_scale_problem(NumericsMatrix* M, NumericsMatrix* H_tilde,
                                        GlobalFrictionContactProblem* problem,
                                        Gfc3d_IPM_init_data* data, size_t nd,
                                        NumericsMatrix** H_out, double** w_out,
                                        double** w_ori_out, NumericsMatrix** minus_H_out,
                                        NumericsMatrix** minus_Ht_out, NumericsMatrix** Ht_out,
                                        NumericsMatrix** minus_M_out) {
  NumericsMatrix* P_mu = data->P_mu->mat;

  double* w_tilde = problem->b;
  double* w = data->tmp_vault_nd[0];

  NumericsMatrix* H = NM_multiply(P_mu, H_tilde);
  NM_gemv(1.0, P_mu, w_tilde, 0.0, w);

  double* w_ori = (double*)calloc(nd, sizeof(double));
  NV_copy(w, nd, w_ori);

  NumericsMatrix* minus_H = NM_copy(H);
  NM_scal(-1.0, minus_H);
  NumericsMatrix* minus_Ht = NM_transpose(minus_H);
  NumericsMatrix* Ht = NM_transpose(H);

  NumericsMatrix* minus_M = NM_copy(M);
  NM_scal(-1.0, minus_M);

  *H_out = H;
  *w_out = w;
  *w_ori_out = w_ori;
  *minus_H_out = minus_H;
  *minus_Ht_out = minus_Ht;
  *Ht_out = Ht;
  *minus_M_out = minus_M;
}

/* Initialize the IPM iterate (reaction, velocity, globalVelocity, s) with a default interior
 * point and, if load_starting_point is enabled, overwrite it from a solution file. */
static void gfc3d_IPM_SNM_compute_starting_point(size_t m, size_t nd, size_t d, size_t n,
                                                 NumericsMatrix* M, NumericsMatrix* H,
                                                 const double* f, const char* strToken,
                                                 double* reaction, double* velocity,
                                                 double* globalVelocity, double* s,
                                                 Gfc3d_IPM_init_data* data) {
  int load_starting_point = 0, pertu_point = 0;

  // Reset vars
  for (size_t i = 0; i < nd; i++)
    if (i % d == 0)
      reaction[i] = 1.;
    else
      reaction[i] = 0.1;

  cblas_dcopy(m, f, 1, globalVelocity, 1);
  NM_tgemv(1.0, H, reaction, 1.0, globalVelocity);
  int info_cholesky = 1;
  if (!data->cholesky_failed) {
    info_cholesky = NM_Cholesky_solve(M, globalVelocity, 1);
    if (info_cholesky) {
      data->cholesky_failed = 1;
      numerics_printf_verbose(1,
                              "Cholesky failed. We continue as if M is the identity matrix");
    }
  } else {
    numerics_printf_verbose(
        2, "gfc3d_IPM_SNM_compute_starting_point: using LU after previous Cholesky failure");
    int info_lu = NM_LU_solve(M, globalVelocity, 1);
    if (info_lu) {
      numerics_printf_verbose(1, "LU failed. We continue as if M is the identity matrix");
    }
  }

  for (size_t i = 0; i < nd; i++)
    if (i % d == 0)
      velocity[i] = 1.;
    else
      velocity[i] = 0.1;

  for (size_t i = 0; i < n; i++) {
    double nub = cblas_dnrm2(2, velocity + i * d + 1, 1);
    s[i] = nub;
  }

  if (load_starting_point) {
    FILE* sol_file = fopen("sol_data_copy.res", "r");
    if (!sol_file)
      printf("\n\ngfc3d_ipm_snm: Solution data file is not available!!! \n\n");
    else {
      char read_prob_name[200];
      int load_v = 0, load_u = 0, load_r = 0, c = 0, len_prob_name = 0, newlineCount = 0;

      // Traverse the problem names in data file for a match
      for (int i = 0; i < 1091; i++) {
        if (fgets(read_prob_name, sizeof(read_prob_name), sol_file) != NULL) {
          len_prob_name = strlen(read_prob_name);
          if (len_prob_name > 0 && read_prob_name[len_prob_name - 1] == '\n') {
            read_prob_name[len_prob_name - 1] =
                '\0';  // Replace the newline character with null terminator
          }
          if (strcmp(read_prob_name, strToken) == 0)  // Problem names are matched
          {
            load_v = load_u = load_r = 1;
            break;
          }
        }

        // Go to the next problem name
        newlineCount = 0;
        while ((c = fgetc(sol_file)) != EOF) {
          if (c == '\n') {
            newlineCount++;
            if (newlineCount == 3) break;  // Stop reading after 2 lines
          }
        }
      }

      // load v
      if (load_v) {
        for (size_t i = 0; i < m; i++) {
          fscanf(sol_file, "%lf ", globalVelocity + i);
        }
        fscanf(sol_file, "\n");
      }

      // load u
      if (load_u) {
        for (size_t i = 0; i < nd; i++) {
          fscanf(sol_file, "%lf ", velocity + i);
        }
        fscanf(sol_file, "\n");
      }

      // load r
      if (load_r) {
        for (size_t i = 0; i < nd; i++) {
          fscanf(sol_file, "%lf ", reaction + i);
        }
        fscanf(sol_file, "\n");
      }

      printf("\ngfc3d_ipm_snm: Starting point is successfully loaded.\n\n");
    }
    fclose(sol_file);

    // Sol perturbation
    if (pertu_point) {
      for (size_t i = 0; i < m; i++) {
        globalVelocity[i] *= 1.1;
      }

      for (size_t i = 0; i < nd; i++) {
        if (i % d == 0) {
          velocity[i] *= 1.2;
          reaction[i] *= 1.1;
        } else {
          velocity[i] *= 0.9;
          reaction[i] *= 0.8;
        }
      }

      for (size_t i = 0; i < n; i++) {
        s[i] *= 1.05;
      }
      printf("\nThe point is successfully perturbed.\n\n");
    }
  }
}

/* Bundle of resources that must be released at the end of gfc3d_IPM_SNM.  Keeping them in one
 * place guarantees a single cleanup path and makes future error propagation easier. */
typedef struct {
  GlobalFrictionContactProblem* problem;
  SolverOptions* options;
  Gfc3d_IPM_init_data* data;
  int internal_allocation;
  int hasNotConverged;
  size_t nd;
  size_t n;
  double* velocity;
  double* reaction;
  double totalresidual;
  double projerr;

  NumericsMatrix* M;
  NumericsMatrix* M_original;
  NumericsMatrix* H_tilde;
  NumericsMatrix* H;
  NumericsMatrix* minus_H;
  NumericsMatrix* minus_Ht;
  NumericsMatrix* Ht;
  NumericsMatrix* minus_M;
  NumericsMatrix* eye_nd;
  NumericsMatrix* eye_n;
  NumericsMatrix* subdiff_u;
  NumericsMatrix* minus_e;
  NumericsMatrix* arrow_r;
  NumericsMatrix* arrow_u;
  NumericsMatrix* J;
  NumericsMatrix* J_dense;

  double* w_ori;
  double* diff_fixp_vec;
  double* fixpConstraint;
  double* arr_norm_theta;
  double* d_globalVelocity;
  double* d_velocity;
  double* d_reaction;
  double* d_s;
  double* s;
  double* velocity_inv;
  double* rhs_2;
  double* rhs_tmp;
  double* sol;
  double* nub_vec;
  int* setR;
  char* str;
  char* blk_num_name;
  FILE* iterates;
  FILE* iterates_2;
} gfc3d_IPM_SNM_cleanup_data;

/* Single cleanup routine for gfc3d_IPM_SNM.  Restores the original scaling, frees memory, and
 * writes back convergence info. */
static void gfc3d_IPM_SNM_cleanup(gfc3d_IPM_SNM_cleanup_data* cd, int* info) {
  /* ----- return to original variables ------ */
  NM_gemv(1.0, cd->data->P_mu->inv_mat, cd->velocity, 0.0, cd->data->tmp_point->t_velocity);
  cblas_dcopy(cd->nd, cd->data->tmp_point->t_velocity, 1, cd->velocity, 1);

  NM_gemv(1.0, cd->data->P_mu->mat, cd->reaction, 0.0, cd->data->tmp_point->t_reaction);
  cblas_dcopy(cd->nd, cd->data->tmp_point->t_reaction, 1, cd->reaction, 1);

  cd->options->dparam[SICONOS_DPARAM_RESIDU] = cd->totalresidual;

  if (cd->internal_allocation) {
    gfc3d_IPM_free(cd->problem, cd->options);
  }

  cd->options->solverData = (double*)malloc(sizeof(double));
  double* projerr_ptr = (double*)cd->options->solverData;
  *projerr_ptr = cd->projerr;

  /* Free any working copy of M that is not the problem's own matrix (e.g. a
   * sparse-storage copy), then restore the original user-provided M if
   * symmetrization replaced it. */
  if (cd->M && cd->M != cd->problem->M) {
    cd->M = NM_free(cd->M);
  }
  if (cd->M_original) {
    cd->problem->M = NM_free(cd->problem->M);
    cd->problem->M = cd->M_original;
  }

  if (cd->H_tilde) cd->H_tilde = NM_free(cd->H_tilde);
  if (cd->minus_H) cd->minus_H = NM_free(cd->minus_H);
  if (cd->H) cd->H = NM_free(cd->H);
  if (cd->minus_M) cd->minus_M = NM_free(cd->minus_M);
  if (cd->minus_Ht) cd->minus_Ht = NM_free(cd->minus_Ht);
  if (cd->eye_nd) cd->eye_nd = NM_free(cd->eye_nd);
  if (cd->eye_n) cd->eye_n = NM_free(cd->eye_n);
  if (cd->Ht) cd->Ht = NM_free(cd->Ht);
  if (cd->subdiff_u) cd->subdiff_u = NM_free(cd->subdiff_u);
  if (cd->minus_e) cd->minus_e = NM_free(cd->minus_e);
  if (cd->arrow_r) cd->arrow_r = NM_free(cd->arrow_r);
  if (cd->arrow_u) cd->arrow_u = NM_free(cd->arrow_u);
  if (cd->J) cd->J = NM_free(cd->J);
  if (cd->blk_num_name) {
    free(cd->blk_num_name);
    cd->blk_num_name = NULL;
  }

  if (cd->iterates) {
    fclose(cd->iterates);
    if (cd->iterates_2) fclose(cd->iterates_2);
  }

  if (cd->rhs_tmp) free(cd->rhs_tmp);
  if (cd->rhs_2) free(cd->rhs_2);
  if (cd->sol) free(cd->sol);

  if (cd->d_globalVelocity) free(cd->d_globalVelocity);
  if (cd->d_velocity) free(cd->d_velocity);
  if (cd->d_reaction) free(cd->d_reaction);
  if (cd->s) free(cd->s);
  if (cd->d_s) free(cd->d_s);
  if (cd->w_ori) free(cd->w_ori);
  if (cd->fixpConstraint) free(cd->fixpConstraint);
  if (cd->velocity_inv) free(cd->velocity_inv);
  if (cd->arr_norm_theta) free(cd->arr_norm_theta);
  if (cd->nub_vec) free(cd->nub_vec);
  if (cd->setR) free(cd->setR);
  if (cd->diff_fixp_vec) free(cd->diff_fixp_vec);
  if (cd->str) {
    free(cd->str);
    cd->str = NULL;
  }

  *info = cd->hasNotConverged;
}

/* Workspace used by gfc3d_IPM_SNM.  Grouping the state in one struct lets us split the
 * per-linear-system-formulation switch bodies into separate functions without passing dozens
 * of arguments. */
typedef struct {
  /* problem dimensions */
  size_t m, nd, d, n, mp2nd;
  int type;
  double tol;
  int max_iter;

  /* algorithm parameters */
  double sgmp1, sgmp2, sgmp3;
  double gmmp1, gmmp2;
  double gmm;
  double scale_sub_diff;

  /* handles */
  GlobalFrictionContactProblem* problem;
  SolverOptions* options;
  Gfc3d_IPM_init_data* data;

  /* matrices */
  NumericsMatrix *M, *H_tilde, *H, *minus_H, *minus_Ht, *Ht, *minus_M;
  NumericsMatrix *eye_nd, *eye_n, *minus_e;
  NumericsMatrix *J, *J_dense, *arrow_r, *arrow_u, *subdiff_u, *Qp2;

  /* vectors */
  double *w, *w_ori, *f;
  double *reaction, *velocity, *globalVelocity, *s;
  double *primalConstraint, *dualConstraint, *complemConstraint, *fixpConstraint;
  double *nub_vec, *diff_fixp_vec;
  double *d_globalVelocity, *d_velocity, *d_reaction, *d_s;
  double *u_plus_du, *r_plus_dr, *dudr_jprod;
  double *p2, *Qp_u, *Qp_du, *Qpinv_dr, *ududr, *Qp_ududr, *velocity_inv;
  double *rhs, *rhs_2, *rhs_tmp, *sol;
  double* arr_norm_theta;
  int* setR;

  /* iteration state */
  int internal_allocation;
  int hasNotConverged;
  int iteration;
  int jacobian_is_nan;
  double pinfeas, dinfeas, complem, udotr, uor_mu, projerr, diff_fixp;
  double totalresidual, totalresidual_Jor;
  double barr_param, sigma, alpha_primal, alpha_dual;
  double LS_norm_p, LS_norm_d, LS_norm_c, LS_norm_f;
  long J_nzmax, J_nz_captured;
  size_t subdiff_u_nzmax, H_nzmax, M_nzmax;

  /* file / string handles */
  FILE* iterates;
  FILE* iterates_2;
  char* strToken;
  char* str;
  char fws;
} gfc3d_IPM_SNM_workspace;

/* ---- per-linear-system-formulation step helpers ---- */
static void gfc3d_IPM_SNM_step_NOSCAL(gfc3d_IPM_SNM_workspace* ws) {
  /* locals not carried in the workspace struct */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
  int jacobian_is_nan = 0;
  size_t pos = 0;
  double nub = 0.0;
  int max_refine = 1;
  double barr_param_a = 0.0, e = 0.0;

  /* problem dimensions and parameters */
  const size_t m = ws->m;
  const size_t nd = ws->nd;
  const size_t d = ws->d;
  const size_t n = ws->n;
  const size_t mp2nd = ws->mp2nd;
  const int type = ws->type;
  const double tol = ws->tol;
  const double sgmp1 = ws->sgmp1;
  const double sgmp2 = ws->sgmp2;
  const double sgmp3 = ws->sgmp3;
  const double gmmp1 = ws->gmmp1;
  const double gmmp2 = ws->gmmp2;
  SolverOptions* options = ws->options;

  /* handles and matrices */
  NumericsMatrix* M = ws->M;
  NumericsMatrix* H = ws->H;
  NumericsMatrix* minus_H = ws->minus_H;
  NumericsMatrix* minus_Ht = ws->minus_Ht;
  NumericsMatrix* minus_M = ws->minus_M;
  NumericsMatrix* eye_nd = ws->eye_nd;
  NumericsMatrix* eye_n = ws->eye_n;
  NumericsMatrix* minus_e = ws->minus_e;
  NumericsMatrix* J = ws->J;
  NumericsMatrix* J_dense = ws->J_dense;
  NumericsMatrix* arrow_r = ws->arrow_r;
  NumericsMatrix* arrow_u = ws->arrow_u;
  NumericsMatrix* subdiff_u = ws->subdiff_u;
  NumericsMatrix* Qp2 = ws->Qp2;

  /* vectors */
  double* w = ws->w;
  double* f = ws->f;
  double* reaction = ws->reaction;
  double* velocity = ws->velocity;
  double* globalVelocity = ws->globalVelocity;
  double* s = ws->s;
  double* primalConstraint = ws->primalConstraint;
  double* dualConstraint = ws->dualConstraint;
  double* complemConstraint = ws->complemConstraint;
  double* fixpConstraint = ws->fixpConstraint;
  double* nub_vec = ws->nub_vec;
  double* diff_fixp_vec = ws->diff_fixp_vec;
  double* d_globalVelocity = ws->d_globalVelocity;
  double* d_velocity = ws->d_velocity;
  double* d_reaction = ws->d_reaction;
  double* d_s = ws->d_s;
  double* u_plus_du = ws->u_plus_du;
  double* r_plus_dr = ws->r_plus_dr;
  double* dudr_jprod = ws->dudr_jprod;
  double* rhs = ws->rhs;
  double* rhs_2 = ws->rhs_2;
  double* rhs_tmp = ws->rhs_tmp;
  double* sol = ws->sol;

  /* iteration state scalars - local copies, copied back at the end */
  int iteration = ws->iteration;
  int hasNotConverged = ws->hasNotConverged;
  double pinfeas = ws->pinfeas;
  double dinfeas = ws->dinfeas;
  double complem = ws->complem;
  double udotr = ws->udotr;
  double uor_mu = ws->uor_mu;
  double projerr = ws->projerr;
  double diff_fixp = ws->diff_fixp;
  double totalresidual = ws->totalresidual;
  double totalresidual_Jor = ws->totalresidual_Jor;
  double barr_param = ws->barr_param;
  double sigma = ws->sigma;
  double alpha_primal = ws->alpha_primal;
  double alpha_dual = ws->alpha_dual;
  double gmm = ws->gmm;
  double scale_sub_diff = ws->scale_sub_diff;
  double LS_norm_p = ws->LS_norm_p;
  double LS_norm_d = ws->LS_norm_d;
  double LS_norm_c = ws->LS_norm_c;
  double LS_norm_f = ws->LS_norm_f;

  /* integer/long bookkeeping */
  long J_nzmax = ws->J_nzmax;
  long J_nz_captured = ws->J_nz_captured;
  size_t subdiff_u_nzmax = ws->subdiff_u_nzmax;

  /* misc */
  FILE* iterates = ws->iterates;
  char fws = ws->fws;

#pragma GCC diagnostic pop

  switch (0) {
    case 0: {
      // Stopping test using norm type
      primalResidual_s_type(velocity, H, globalVelocity, w, s, primalConstraint, &pinfeas, tol,
                            type);
      dualResidual_type(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol,
                        type);
      // udotr = xdoty_type(n, nd, velocity, reaction, type);
      udotr = xdoty_type(n, nd, velocity, reaction, NORM_2_INF);
      // complem = complemResidualNorm(velocity, reaction, nd, n);
      complem = complemResidualNorm_type(velocity, reaction, nd, n, NORM_2_INF);
      barr_param = cblas_ddot(nd, reaction, 1, velocity, 1) / n;

      diff_fixp = 0.;
      for (size_t i = 0; i < n; i++) {
        pos = i * d;

        nub_vec[i] = cblas_dnrm2(2, velocity + pos + 1, 1);
        diff_fixp_vec[i] = s[i] - nub_vec[i];
      }

      // Stopping test using norm type
      if (type == NORM_2)
        diff_fixp = cblas_dnrm2(n, diff_fixp_vec, 1);
      else if (type == NORM_INF) {
        diff_fixp = fabs(diff_fixp_vec[0]);
        for (size_t i = 1; i < n; i++) {
          diff_fixp = fmax(diff_fixp, fabs(diff_fixp_vec[i]));
        }
      } else {
        assert(0);
      }

      // totalresidual = fmax(fmax(fmax(pinfeas, dinfeas),diff_fixp),complem);
      totalresidual = fmax(fmax(fmax(pinfeas, dinfeas), diff_fixp), udotr);
      totalresidual_Jor = fmax(fmax(fmax(pinfeas, dinfeas), diff_fixp), complem);

      // // compute Projection Error

      // projerr = projectionError(velocity, reaction, n, tol);
      projerr = projectionError_based_reaction_norm_infinity_conic(H, M, f, w, reaction, n,
                                                                   NOT_ON_DUAL_CONE, ws->data);

      // if (total_time > max_time)
      if (0 == 1) {
        // FP: add comment to be able to compile when -Werror,-Wunreachable-code
        //   hasNotConverged = 5;
        //   numerics_printf_verbose(
        //       -1, "| %3i%c| %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e ||",
        //       iteration, fws, pinfeas, dinfeas, diff_fixp, udotr, complem, projerr,
        //       barr_param, sigma);
        //   printf("\nREQUEST TIMEOUT\n");
        //   break;
        // } else if (alpha_primal <= 1e-13) {
        //   hasNotConverged = 6;
        //   numerics_printf_verbose(
        //       -1, "| %3i%c| %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e ||",
        //       iteration, fws, pinfeas, dinfeas, diff_fixp, udotr, complem, projerr,
        //       barr_param, sigma);
        //   printf("\nTiny step-length\n\n");
        //   break;
      }

      // if ( projerr <= tol && totalresidual_Jor <= tol )
      // if ( projerr <= tol )
      if (totalresidual_Jor <= tol) {
        double unitur;
        for (size_t i = 0; i < n; i++) {
          unitur = cblas_ddot(3, velocity + 3 * i, 1, reaction + 3 * i, 1);
          if (unitur < 0) printf("UR NEGATIF %9.2e\n", unitur);
        }

        numerics_printf_verbose(
            1, "| %3i%c| %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e ||", iteration,
            fws, pinfeas, dinfeas, diff_fixp, udotr, complem, projerr, barr_param, sigma);

        hasNotConverged = 0;
        break;
      }

      // scale_sub_diff = fmax(0.9, alpha_primal);
      scale_sub_diff = 1.;
      if (iteration == 0) {
        // First linear linear system
        // Create Jac
        J = NM_create(NM_SPARSE, m + 2 * nd + n, m + 2 * nd + n);
        NM_triplet_alloc(J, J_nzmax);
        J->sparse->origin = NSM_TRIPLET;

        // Insert fixed data
        NM_insert(J, M, 0, 0);
        NM_insert(J, minus_H, m + nd, 0);
        NM_insert(J, minus_Ht, 0, m + nd);
        NM_insert(J, eye_nd, m + nd, m);
        NM_insert(J, minus_e, m + nd, m + 2 * nd);
        NM_insert(J, eye_n, m + 2 * nd, m + 2 * nd);

        // /* regularization */
        // // double regul = -1.*fmin(1e-6, 1.*fmax(diff_fixp, udotr));
        // double regul = -1.*1e-8;
        // NM_insert(J, NM_scalar(nd, regul), m + nd, m + nd);

        J_nz_captured = J->sparse->triplet
                            ->nz;  // Save nz at this moment for resetting not fixed data after

        // Compute NOT fixed data
        arrow_r = Arrow_repr(reaction, nd, n);  // Allocation is done in Arrow_repr
        arrow_u = Arrow_repr(velocity, nd, n);

        subdiff_u = NM_create(NM_SPARSE, n, nd);  // Allocation of subdiff_u
        NM_triplet_alloc(subdiff_u, subdiff_u_nzmax);

        for (size_t i = 0; i < n; ++i) {
          pos = i * d;
          nub = nub_vec[i];

          NM_entry(subdiff_u, i, pos + 1, -1. * scale_sub_diff * velocity[pos + 1] / nub);
          NM_entry(subdiff_u, i, pos + 2, -1. * scale_sub_diff * velocity[pos + 2] / nub);

          fixpConstraint[i] = diff_fixp_vec[i];
        }

        // Insert NOT fixed data to Jac
        NM_insert(J, arrow_r, m, m);
        NM_insert(J, subdiff_u, m + 2 * nd, m);
        NM_insert(J, arrow_u, m, m + nd);

        // J_nz_final = J->sparse->triplet->nz;
      }  // end if iter == 0

      else {
        // Clear memory no longer used
        if (arrow_r) {
          arrow_r = NM_free(arrow_r);
        }
        if (arrow_u) {
          arrow_u = NM_free(arrow_u);
        }
        if (subdiff_u) {
          subdiff_u = NM_free(subdiff_u);
        }

        // Clear internal data
        NM_clearCSC(J);
        NM_clearDense(J);
        NM_clearCSCTranspose(J);
        NM_clearCSR(J);
        if (J->sparse->linearSolverParams) {
          NSM_linearSolverParams_free(J->sparse->linearSolverParams);
          J->sparse->linearSolverParams = NULL;
        }
        NM_internalData_free(J);
        if (!NM_destructible(J)) {
          NM_clear(J->destructible);
          J->destructible = J;
        }
        J->sparse->triplet->nz = J_nz_captured;  // nz goes back the end point of fixed data

        // Replace the matrix content, without reallocation
        NM_insert_Arrow_to_Triplet(J->sparse->triplet, m, m, reaction, nd, n);
        NM_insert_Arrow_to_Triplet(J->sparse->triplet, m, m + nd, velocity, nd, n);

        CSparseMatrix* J_cs = J->sparse->triplet;

        for (size_t i = 0; i < n; i++) {
          pos = i * d;
          // nub = cblas_dnrm2(2, velocity+pos+1, 1);
          nub = nub_vec[i];

          for (size_t j = 1; j < d; j++) {
            J_cs->x[J_cs->nz] = -1. * scale_sub_diff * velocity[pos + j] / nub;
            J_cs->i[J_cs->nz] = mp2nd + i;
            J_cs->p[J_cs->nz++] = m + pos + j;
          }

          // fixpConstraint[i] = s[i] - nub;  // fixpConstraint = s - |u_bar|
          fixpConstraint[i] = diff_fixp_vec[i];
        }
      }

      jacobian_is_nan = NM_isnan(J);
      if (jacobian_is_nan) {
        numerics_printf_verbose(0, "The Jacobian matrix J contains NaN");
        break;
      }

      // Compute rhs
      cblas_dcopy(m, dualConstraint, 1, rhs, 1);
      JA_prod(velocity, reaction, nd, n, complemConstraint);
      cblas_dcopy(nd, complemConstraint, 1, rhs + m, 1);
      cblas_dcopy(nd, primalConstraint, 1, rhs + m + nd, 1);
      cblas_dcopy(n, fixpConstraint, 1, rhs + m + 2 * nd, 1);

      cblas_dscal(m + 2 * nd + n, -1.0, rhs, 1);
      cblas_dcopy(m + 2 * nd + n, rhs, 1, rhs_2,
                  1);  // rhs_2 = old rhs (already had minus sign)

      if (cblas_dnrm2(n, fixpConstraint, 1) > complemResidualNorm(velocity, reaction, nd, n)) {
        sigma = 0.499;
      } else {
        // Solve
        cblas_dcopy(m + 2 * nd + n, rhs_2, 1, rhs_tmp, 1);
        for (size_t k = 0; k < m + 2 * nd + n; sol[k] = 0., k++);  // reset sol

        max_refine = 1;
        if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] ==
            SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES)
          max_refine = 10;

        for (int itr = 0; itr < max_refine; itr++) {
          NM_LU_solve(J, rhs_tmp, 1);
          cblas_daxpy(m + 2 * nd + n, 1.0, rhs_tmp, 1, sol, 1);
          cblas_dcopy(m + 2 * nd + n, rhs_2, 1, rhs_tmp, 1);
          NM_gemv(-1.0, J, sol, 1.0, rhs_tmp);
          // printf("refinement iterations = %i %8.2e\n",itr, cblas_dnrm2(m+2*nd, rhs_tmp,
          // 1));
          if (cblas_dnrm2(m + 2 * nd + n, rhs_tmp, 1) <= 1e-14) {
            // printf("\nrefinement iterations = %d %8.2e\n",itr+1, cblas_dnrm2(m+2*nd+n,
            // rhs_tmp, 1));
            //              free(rhs_tmp);
            break;
          }
        }

        // Compute direction
        cblas_dcopy(m, sol, 1, d_globalVelocity, 1);
        cblas_dcopy(nd, sol + m, 1, d_velocity, 1);
        cblas_dcopy(nd, sol + m + nd, 1, d_reaction, 1);
        cblas_dcopy(n, sol + m + 2 * nd, 1, d_s, 1);

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

        barr_param_a = cblas_ddot(nd, u_plus_du, 1, r_plus_dr, 1) / n;
        e = barr_param > sgmp1 ? fmax(1.0, sgmp2 * pow(alpha_primal, 2)) : sgmp3;
        // sigma = 0.49*fmin(1.0, pow(barr_param_a / barr_param, e));
        sigma = 0.499 * fmin(1., pow(barr_param_a / barr_param, e));

        // Compute rhs for 2nd linear system
        cblas_dcopy(m + 2 * nd + n, rhs_2, 1, rhs, 1);  // Get back value for rhs

        if (sigma < 0.25) {
          JA_prod(d_velocity, d_reaction, nd, n, dudr_jprod);
          cblas_daxpy(nd, -1., dudr_jprod, 1, rhs + m, 1);
        }
      }

      for (size_t k = 0; k < nd; rhs[m + k] += 2 * sigma * barr_param, k += d);
      cblas_dcopy(m + 2 * nd + n, rhs, 1, rhs_2, 1);  // rhs_2 = old rhs

      // SOLVE
      // NM_LU_solve(J, rhs, 1);

      cblas_dcopy(m + 2 * nd + n, rhs_2, 1, rhs_tmp, 1);
      for (size_t k = 0; k < m + 2 * nd + n; sol[k] = 0., k++);  // reset sol

      max_refine = 1;
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] ==
          SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES)
        max_refine = 10;

      for (int itr = 0; itr < max_refine; itr++) {
        NM_LU_solve(J, rhs_tmp, 1);
        cblas_daxpy(m + 2 * nd + n, 1.0, rhs_tmp, 1, sol, 1);
        cblas_dcopy(m + 2 * nd + n, rhs_2, 1, rhs_tmp, 1);
        NM_gemv(-1.0, J, sol, 1.0, rhs_tmp);
        // printf("refinement iterations = %i %8.2e\n",itr, cblas_dnrm2(m+2*nd, rhs_tmp, 1));
        if (cblas_dnrm2(m + 2 * nd + n, rhs_tmp, 1) <= 1e-14) {
          // printf("\nrefinement iterations = %d %8.2e\n",itr+1, cblas_dnrm2(m+2*nd+n,
          // rhs_tmp, 1)); free(rhs_tmp);
          break;
        }
      }

      NM_gemv(1.0, J, sol, -1.0, rhs_2);

      LS_norm_d = cblas_dnrm2(m, rhs_2, 1);
      LS_norm_c = cblas_dnrm2(nd, rhs_2 + m, 1);
      LS_norm_p = cblas_dnrm2(nd, rhs_2 + m + nd, 1);
      LS_norm_f = cblas_dnrm2(n, rhs_2 + m + 2 * nd, 1);

      cblas_dcopy(m, sol, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, sol + m, 1, d_velocity, 1);
      cblas_dcopy(nd, sol + m + nd, 1, d_reaction, 1);
      cblas_dcopy(n, sol + m + 2 * nd, 1, d_s, 1);

      /* computing the step-length */
      alpha_primal = getStepLength(velocity, d_velocity, nd, n, gmm);
      alpha_dual = getStepLength(reaction, d_reaction, nd, n, gmm);
      alpha_primal = alpha_dual = fmin(alpha_primal, alpha_dual);
      gmm = gmmp1 + gmmp2 * alpha_primal;

      cblas_daxpy(m, alpha_primal, d_globalVelocity, 1, globalVelocity, 1);
      cblas_daxpy(nd, alpha_primal, d_velocity, 1, velocity, 1);
      cblas_daxpy(nd, alpha_primal, d_reaction, 1, reaction, 1);
      cblas_daxpy(n, alpha_primal, d_s, 1, s, 1);

      // total_time += (double)(t2 - t1) / (double)clk_tck;

      if (NV_isnan(globalVelocity, m) | NV_isnan(velocity, nd) | NV_isnan(reaction, nd) |
          NV_isnan(s, n)) {
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
      * Case of a non-reduced system. Building the right-hand side related to the first
      system

        with NT scaling
        rhs = -
        [      M*v - H'*r - f      ]  m         dualConstraint
        [            r             ]  nd        complemConstraint    r - 2 mu u^-1
        [   u - H*v - w - phi(s)   ]  nd        primalConstraint
        [         s - |ub|         ]  n
     */
  }
  /* ---- copy back modified state ---- */
  ws->hasNotConverged = hasNotConverged;
  ws->jacobian_is_nan = jacobian_is_nan;
  ws->pinfeas = pinfeas;
  ws->dinfeas = dinfeas;
  ws->complem = complem;
  ws->udotr = udotr;
  ws->uor_mu = uor_mu;
  ws->projerr = projerr;
  ws->diff_fixp = diff_fixp;
  ws->totalresidual = totalresidual;
  ws->totalresidual_Jor = totalresidual_Jor;
  ws->barr_param = barr_param;
  ws->sigma = sigma;
  ws->alpha_primal = alpha_primal;
  ws->alpha_dual = alpha_dual;
  ws->gmm = gmm;
  ws->scale_sub_diff = scale_sub_diff;
  ws->LS_norm_p = LS_norm_p;
  ws->LS_norm_d = LS_norm_d;
  ws->LS_norm_c = LS_norm_c;
  ws->LS_norm_f = LS_norm_f;
  ws->J_nzmax = J_nzmax;
  ws->J_nz_captured = J_nz_captured;

  ws->J = J;
  ws->J_dense = J_dense;
  ws->arrow_r = arrow_r;
  ws->arrow_u = arrow_u;
  ws->subdiff_u = subdiff_u;
  ws->Qp2 = Qp2;
}

static void gfc3d_IPM_SNM_step_QP2(gfc3d_IPM_SNM_workspace* ws) {
  /* locals not carried in the workspace struct */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
  int jacobian_is_nan = 0;
  size_t pos = 0;
  double nub = 0.0;
  int max_refine = 1;
  double barr_param_a = 0.0, e = 0.0;

  /* problem dimensions and parameters */
  const size_t m = ws->m;
  const size_t nd = ws->nd;
  const size_t d = ws->d;
  const size_t n = ws->n;
  const size_t mp2nd = ws->mp2nd;
  const int type = ws->type;
  const double tol = ws->tol;
  const double sgmp1 = ws->sgmp1;
  const double sgmp2 = ws->sgmp2;
  const double sgmp3 = ws->sgmp3;
  const double gmmp1 = ws->gmmp1;
  const double gmmp2 = ws->gmmp2;
  SolverOptions* options = ws->options;

  /* handles and matrices */
  NumericsMatrix* M = ws->M;
  NumericsMatrix* H = ws->H;
  NumericsMatrix* minus_H = ws->minus_H;
  NumericsMatrix* minus_Ht = ws->minus_Ht;
  NumericsMatrix* minus_M = ws->minus_M;
  NumericsMatrix* eye_nd = ws->eye_nd;
  NumericsMatrix* eye_n = ws->eye_n;
  NumericsMatrix* minus_e = ws->minus_e;
  NumericsMatrix* J = ws->J;
  NumericsMatrix* J_dense = ws->J_dense;
  NumericsMatrix* arrow_r = ws->arrow_r;
  NumericsMatrix* arrow_u = ws->arrow_u;
  NumericsMatrix* subdiff_u = ws->subdiff_u;

  /* vectors */
  double* w = ws->w;
  double* f = ws->f;
  double* reaction = ws->reaction;
  double* velocity = ws->velocity;
  double* globalVelocity = ws->globalVelocity;
  double* s = ws->s;
  double* primalConstraint = ws->primalConstraint;
  double* dualConstraint = ws->dualConstraint;
  double* complemConstraint = ws->complemConstraint;
  double* fixpConstraint = ws->fixpConstraint;
  double* nub_vec = ws->nub_vec;
  double* diff_fixp_vec = ws->diff_fixp_vec;
  double* d_globalVelocity = ws->d_globalVelocity;
  double* d_velocity = ws->d_velocity;
  double* d_reaction = ws->d_reaction;
  double* d_s = ws->d_s;
  double* u_plus_du = ws->u_plus_du;
  double* r_plus_dr = ws->r_plus_dr;
  double* dudr_jprod = ws->dudr_jprod;
  double* rhs = ws->rhs;
  double* rhs_2 = ws->rhs_2;
  double* rhs_tmp = ws->rhs_tmp;
  double* sol = ws->sol;

  /* iteration state scalars - local copies, copied back at the end */
  int iteration = ws->iteration;
  int hasNotConverged = ws->hasNotConverged;
  double pinfeas = ws->pinfeas;
  double dinfeas = ws->dinfeas;
  double complem = ws->complem;
  double udotr = ws->udotr;
  double uor_mu = ws->uor_mu;
  double projerr = ws->projerr;
  double diff_fixp = ws->diff_fixp;
  double totalresidual = ws->totalresidual;
  double totalresidual_Jor = ws->totalresidual_Jor;
  double barr_param = ws->barr_param;
  double sigma = ws->sigma;
  double alpha_primal = ws->alpha_primal;
  double alpha_dual = ws->alpha_dual;
  double gmm = ws->gmm;
  double scale_sub_diff = ws->scale_sub_diff;
  double LS_norm_p = ws->LS_norm_p;
  double LS_norm_d = ws->LS_norm_d;
  double LS_norm_c = ws->LS_norm_c;
  double LS_norm_f = ws->LS_norm_f;

  /* integer/long bookkeeping */
  long J_nzmax = ws->J_nzmax;
  long J_nz_captured = ws->J_nz_captured;

  /* misc */
  FILE* iterates = ws->iterates;
  char fws = ws->fws;

  NumericsMatrix* Qp2 = ws->Qp2;
  double* p2 = ws->p2;
  double* Qp_u = ws->Qp_u;
  double* Qp_du = ws->Qp_du;
  double* Qpinv_dr = ws->Qpinv_dr;
  double* ududr = ws->ududr;
  double* Qp_ududr = ws->Qp_ududr;
  double* velocity_inv = ws->velocity_inv;
  size_t H_nzmax = ws->H_nzmax;
  size_t M_nzmax = ws->M_nzmax;

#pragma GCC diagnostic pop

  switch (0) {
    case 0: {
      // Stopping test using norm type
      primalResidual_s_type(velocity, H, globalVelocity, w, s, primalConstraint, &pinfeas, tol,
                            type);
      dualResidual_type(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol,
                        type);
      // udotr = xdoty_type(n, nd, velocity, reaction, type);
      udotr = xdoty_type(n, nd, velocity, reaction, NORM_2_INF);
      complem = complemResidualNorm(velocity, reaction, nd, n);
      barr_param = cblas_ddot(nd, reaction, 1, velocity, 1) / n;

      diff_fixp = 0.;
      for (size_t i = 0; i < nd; i += d) {
        nub = cblas_dnrm2(2, velocity + i + 1, 1);
        diff_fixp_vec[i / d] = fabs(s[i / d] - nub);
      }

      // Stopping test using norm type
      if (type == NORM_2)
        diff_fixp = cblas_dnrm2(n, diff_fixp_vec, 1);
      else if (type == NORM_INF) {
        int maxIndex = cblas_idamax(n, diff_fixp_vec, 1);
        diff_fixp = fabs(diff_fixp_vec[maxIndex]);
      } else {
        assert(0);
      }

      // totalresidual = fmax(fmax(fmax(pinfeas, dinfeas),diff_fixp),complem);
      totalresidual = fmax(fmax(fmax(pinfeas, dinfeas), diff_fixp), udotr);

      // compute Projection Error
      // projerr = projectionError(velocity, reaction, n, tol);
      projerr = projectionError_based_reaction_norm_infinity_conic(H, M, f, w, reaction, n,
                                                                   NOT_ON_DUAL_CONE, ws->data);

      if (totalresidual <= tol) {
        double unitur;
        for (size_t i = 0; i < n; i++) {
          unitur = cblas_ddot(3, velocity + 3 * i, 1, reaction + 3 * i, 1);
          if (unitur < 0) printf("UR NEGATIF %9.2e\n", unitur);
        }

        if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE]) {
          fprintf(iterates,
                  "%d %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e "
                  "%.10e %.10e %.10e %.10e %.10e;\n",
                  iteration, pinfeas, dinfeas, diff_fixp, uor_mu, udotr, projerr, barr_param,
                  alpha_primal, fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
                  fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
                  fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]),
                  fabs(d_s[cblas_idamax(n, d_s, 1)]), LS_norm_p, LS_norm_d, LS_norm_c,
                  LS_norm_f);
        }

        numerics_printf_verbose(
            1, "| %3i%c| %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e ||", iteration,
            fws, pinfeas, dinfeas, diff_fixp, udotr, complem, projerr, barr_param, sigma);

        hasNotConverged = 0;
        break;
      }

      // First linear linear system
      J = NM_create(NM_SPARSE, m + 2 * nd + n, m + 2 * nd + n);
      J_nzmax = M_nzmax + H_nzmax + 2 * (d * 3 - 2) * n + H_nzmax + nd + n + 2 * n + n;
      NM_triplet_alloc(J, J_nzmax);
      J->sparse->origin = NSM_TRIPLET;

      // arrow_r = Arrow_repr(reaction, nd, n);
      // arrow_u = Arrow_repr(velocity, nd, n) ;

      // Create subdiff_u
      subdiff_u = NM_create(NM_SPARSE, n, nd);
      size_t subdiff_u_nzmax = 2 * n;
      NM_triplet_alloc(subdiff_u, subdiff_u_nzmax);

      /* Matrix filling */
      scale_sub_diff = fmax(0.9, alpha_primal);
      for (size_t i = 0; i < n; ++i) {
        pos = i * d;
        nub = cblas_dnrm2(2, velocity + pos + 1, 1);

        NM_entry(subdiff_u, i, pos + 1, -1. * scale_sub_diff * velocity[pos + 1] / nub);
        NM_entry(subdiff_u, i, pos + 2, -1. * scale_sub_diff * velocity[pos + 2] / nub);

        fixpConstraint[i] = s[i] - nub;  // fixpConstraint = s - |u_bar|
      }

      Nesterov_Todd_vector(2, velocity, reaction, nd, n, p2);
      Qp2 = QRmat(p2, nd, n);

      NM_insert(J, M, 0, 0);
      NM_insert(J, minus_H, m + nd, 0);

      NM_insert(J, Qp2, m, m);
      NM_insert(J, eye_nd, m + nd, m);
      NM_insert(J, subdiff_u, m + 2 * nd, m);

      NM_insert(J, minus_Ht, 0, m + nd);
      NM_insert(J, eye_nd, m, m + nd);

      NM_insert(J, minus_e, m + nd, m + 2 * nd);
      NM_insert(J, eye_n, m + 2 * nd, m + 2 * nd);

      // /* regularization */
      // // double regul = -1.*fmin(1e-6, 1.*fmax(diff_fixp, udotr));
      // double regul = -1.*1e-6;
      // NM_insert(J, NM_scalar(nd, regul), m + nd, m + nd);

      if (subdiff_u) {
        subdiff_u = NM_free(subdiff_u);
      }
      if (Qp2) {
        Qp2 = NM_free(Qp2);
      }

      jacobian_is_nan = NM_isnan(J);
      if (jacobian_is_nan) {
        numerics_printf_verbose(0, "The Jacobian matrix J contains NaN");
        break;
      }

      // Compute rhs
      cblas_dcopy(m, dualConstraint, 1, rhs, 1);
      cblas_dcopy(nd, reaction, 1, rhs + m, 1);
      cblas_dcopy(nd, primalConstraint, 1, rhs + m + nd, 1);
      cblas_dcopy(n, fixpConstraint, 1, rhs + m + 2 * nd, 1);

      cblas_dscal(m + 2 * nd + n, -1.0, rhs, 1);
      cblas_dcopy(m + 2 * nd + n, rhs, 1, rhs_2, 1);  // rhs_2 = old rhs

      // Solve
      // rhs_tmp = (double*)calloc(m + 2 * nd + n, sizeof(double));
      cblas_dcopy(m + 2 * nd + n, rhs_2, 1, rhs_tmp, 1);
      for (size_t k = 0; k < m + 2 * nd + n; sol[k] = 0., k++);  // reset sol

      max_refine = 1;
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] ==
          SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES)
        max_refine = 10;

      for (int itr = 0; itr < max_refine; itr++) {
        NM_LU_solve(J, rhs_tmp, 1);
        cblas_daxpy(m + 2 * nd + n, 1.0, rhs_tmp, 1, sol, 1);
        cblas_dcopy(m + 2 * nd + n, rhs_2, 1, rhs_tmp, 1);
        NM_gemv(-1.0, J, sol, 1.0, rhs_tmp);
        // printf("refinement iterations = %i %8.2e\n",itr, cblas_dnrm2(m+2*nd, rhs_tmp, 1));
        if (cblas_dnrm2(m + 2 * nd + n, rhs_tmp, 1) <= 1e-14) {
          // printf("\nrefinement iterations = %d %8.2e\n",itr+1, cblas_dnrm2(m+2*nd+n,
          // rhs_tmp, 1));
          break;
        }
      }

      // Compute direction
      cblas_dcopy(m, sol, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, sol + m, 1, d_velocity, 1);
      cblas_dcopy(nd, sol + m + nd, 1, d_reaction, 1);
      cblas_dcopy(n, sol + m + 2 * nd, 1, d_s, 1);

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

      barr_param_a = cblas_ddot(nd, u_plus_du, 1, r_plus_dr, 1) / n;
      e = barr_param > sgmp1 ? fmax(1.0, sgmp2 * pow(alpha_primal, 2)) : sgmp3;
      sigma = 0.5 * fmin(1.0, pow(barr_param_a / barr_param, e));

      // if (diff_fixp > 10.*udotr)
      if (cblas_dnrm2(n, diff_fixp_vec, 1) > 1. * complem) {
        sigma = 0.5;
      }

      // Compute rhs for 2nd linear system
      cblas_dcopy(m + 2 * nd + n, rhs_2, 1, rhs,
                  1);  // Get back value for rhs (already had minus sign inside)

      QNTpz(velocity, reaction, d_velocity, nd, n, Qp_du);  // Qp_du = du_affine_hat
      QNTpinvz(velocity, reaction, d_reaction, nd, n,
               Qpinv_dr);                                // Qpinv_dr = dr_affine_check
      JA_prod(Qp_du, Qpinv_dr, nd, n, dudr_jprod);       // dudr_jprod = Qp_du o Qpinv_dr
      QNTpz(velocity, reaction, velocity, nd, n, Qp_u);  // Qp_u = u_hat
      Jxinvprody(Qp_u, dudr_jprod, nd, n, ududr);  // ududr = u_hat_inv o (Qp_du o Qpinv_dr)
      QNTpz(velocity, reaction, ududr, nd, n,
            Qp_ududr);  // Qp_ududr = Qp * [ u_hat_inv o (Qp_du o Qpinv_dr) ]

      JA_inv(velocity, nd, n, velocity_inv);  // velocity_inv = u^-1
      cblas_dscal(nd, 2. * sigma * barr_param, velocity_inv,
                  1);  // velocity_inv = 2 * sigma * mu * u^-1

      cblas_daxpy(nd, -1.0, Qp_ududr, 1, rhs + m, 1);
      cblas_daxpy(nd, 1.0, velocity_inv, 1, rhs + m, 1);

      cblas_dcopy(m + 2 * nd + n, rhs, 1, rhs_2, 1);  // rhs_2 = old rhs

      // SOLVE
      // NM_LU_solve(J, rhs, 1);

      cblas_dcopy(m + 2 * nd + n, rhs_2, 1, rhs_tmp, 1);
      for (size_t k = 0; k < m + 2 * nd + n; sol[k] = 0., k++);  // reset sol

      max_refine = 1;
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] ==
          SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES)
        max_refine = 10;

      for (int itr = 0; itr < max_refine; itr++) {
        NM_LU_solve(J, rhs_tmp, 1);
        cblas_daxpy(m + 2 * nd + n, 1.0, rhs_tmp, 1, sol, 1);
        cblas_dcopy(m + 2 * nd + n, rhs_2, 1, rhs_tmp, 1);
        NM_gemv(-1.0, J, sol, 1.0, rhs_tmp);
        // printf("refinement iterations = %i %8.2e\n",itr, cblas_dnrm2(m+2*nd, rhs_tmp, 1));
        if (cblas_dnrm2(m + 2 * nd + n, rhs_tmp, 1) <= 1e-14) {
          // printf("\nrefinement iterations = %d %8.2e\n",itr+1, cblas_dnrm2(m+2*nd+n,
          // rhs_tmp, 1));
          break;
        }
      }

      NM_gemv(1.0, J, sol, -1.0, rhs_2);

      LS_norm_d = cblas_dnrm2(m, rhs_2, 1);
      LS_norm_c = cblas_dnrm2(nd, rhs_2 + m, 1);
      LS_norm_p = cblas_dnrm2(nd, rhs_2 + m + nd, 1);
      LS_norm_f = cblas_dnrm2(n, rhs_2 + m + 2 * nd, 1);

      cblas_dcopy(m, sol, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, sol + m, 1, d_velocity, 1);
      cblas_dcopy(nd, sol + m + nd, 1, d_reaction, 1);
      cblas_dcopy(n, sol + m + 2 * nd, 1, d_s, 1);

      /* computing the step-length */
      alpha_primal = getStepLength(velocity, d_velocity, nd, n, gmm);
      alpha_dual = getStepLength(reaction, d_reaction, nd, n, gmm);
      alpha_primal = alpha_dual = fmin(alpha_primal, alpha_dual);

      cblas_daxpy(m, alpha_primal, d_globalVelocity, 1, globalVelocity, 1);
      cblas_daxpy(nd, alpha_primal, d_velocity, 1, velocity, 1);
      cblas_daxpy(nd, alpha_primal, d_reaction, 1, reaction, 1);
      cblas_daxpy(n, alpha_primal, d_s, 1, s, 1);

      if (NV_isnan(globalVelocity, m) | NV_isnan(velocity, nd) | NV_isnan(reaction, nd) |
          NV_isnan(s, n)) {
        hasNotConverged = 2;
        break;
      }

      break;
    }
  }
  /* ---- copy back modified state ---- */
  ws->hasNotConverged = hasNotConverged;
  ws->jacobian_is_nan = jacobian_is_nan;
  ws->pinfeas = pinfeas;
  ws->dinfeas = dinfeas;
  ws->complem = complem;
  ws->udotr = udotr;
  ws->uor_mu = uor_mu;
  ws->projerr = projerr;
  ws->diff_fixp = diff_fixp;
  ws->totalresidual = totalresidual;
  ws->totalresidual_Jor = totalresidual_Jor;
  ws->barr_param = barr_param;
  ws->sigma = sigma;
  ws->alpha_primal = alpha_primal;
  ws->alpha_dual = alpha_dual;
  ws->gmm = gmm;
  ws->scale_sub_diff = scale_sub_diff;
  ws->LS_norm_p = LS_norm_p;
  ws->LS_norm_d = LS_norm_d;
  ws->LS_norm_c = LS_norm_c;
  ws->LS_norm_f = LS_norm_f;
  ws->J_nzmax = J_nzmax;
  ws->J_nz_captured = J_nz_captured;

  ws->J = J;
  ws->J_dense = J_dense;
  ws->arrow_r = arrow_r;
  ws->arrow_u = arrow_u;
  ws->subdiff_u = subdiff_u;
  ws->Qp2 = Qp2;
}

/*
 * Implementation contains the following functions:
 *  - init ---> same routine as for gfc3d_IPM
 *  - free ---> deallocate memory, same routine as for gfc3d_ipm
 *  - gfc3d_IPM_SNM_setDefaultSolverOptions - setup default solver parameters
 *  - gfc3d_IPM_SNM - optimization method
 */
void gfc3d_IPM_SNM(GlobalFrictionContactProblem* restrict problem, double* restrict reaction,
                   double* restrict velocity, double* restrict globalVelocity,
                   int* restrict info, SolverOptions* restrict options) {
  // verbose = 3;

  // int type = NORM_2;
  int type = NORM_INF;

  char* blk_num_name = NULL;

  // the size of the problem detection
  size_t m = problem->M->size0;
  size_t nd = problem->H->size1;
  size_t d = problem->dimension;
  size_t n = problem->numberOfContacts;

  size_t mp2nd = m + 2 * nd;
  size_t no_m = 0, no_nd = 0;

  NumericsMatrix* M = NULL;
  NumericsMatrix* M_original = NULL;
  NumericsMatrix* H_tilde = NULL;

  /* globalFrictionContactProblem_display(problem); */

  gfc3d_IPM_SNM_prepare_problem(problem, options, m, nd, &M, &M_original, &H_tilde);

  // initialize solver if it is not set
  int internal_allocation = 0;
  size_t work_size = m + nd + nd + nd / d;
  if (!options->dWork || (options->dWorkSize < work_size)) {
    gfc3d_IPM_init(problem, options, work_size);
    internal_allocation = 1;
  }

  Gfc3d_IPM_init_data* data = (Gfc3d_IPM_init_data*)options->solverData;
  data->cholesky_failed = 0;

  double* f = problem->q;

  NumericsMatrix* H = NULL;
  double* w = NULL;
  double* w_ori = NULL;
  NumericsMatrix* minus_H = NULL;
  NumericsMatrix* minus_Ht = NULL;
  NumericsMatrix* Ht = NULL;
  NumericsMatrix* minus_M = NULL;
  gfc3d_IPM_SNM_scale_problem(M, H_tilde, problem, data, nd, &H, &w, &w_ori, &minus_H,
                              &minus_Ht, &Ht, &minus_M);
  no_nd++;

  double alpha_primal = data->internal_params->alpha_primal;
  double alpha_dual = data->internal_params->alpha_dual;
  double barr_param = data->internal_params->barr_param;
  double sigma = data->internal_params->sigma;
  sigma = 0.3;
  // double alpha_complem = 0., alpha_diffixP = 0.;
  // double min_alpha_primal = 1.;

  // /* COMPUTATION OF A NEW STARTING POINT */
  // // set the reaction vector to an arbitrary value in the interior of the cone
  // for (size_t  i = 0; i<nd; i++)
  //     if (i % d == 0) reaction[i] = 0.1;
  //     else reaction[i] = 0.01;

  // // computation of the global velocity vector: v = M\(H'*r+f)
  // for (size_t  i = 0; i<m; i++) globalVelocity[i] = f[i];
  // NM_tgemv(1.0, H, reaction, 1.0, globalVelocity);
  // NM_Cholesky_solve(NM_preserve(M), globalVelocity, 1);

  // for (size_t  i = 0; i<nd; i++)
  //     if (i % d == 0) velocity[i] = 0.1;
  //     else velocity[i] = 0.01;

  double tol = options->dparam[SICONOS_DPARAM_TOL];
  int max_iter = options->iparam[SICONOS_IPARAM_MAX_ITER];

  double sgmp1 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1];
  double sgmp2 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2];
  double sgmp3 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3];
  double gmmp1 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1];
  double gmmp2 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2];

  int hasNotConverged = 1;
  int iteration = 0;
  double pinfeas = 1e300;
  //  double pinfeas_new = 1e300;
  double dinfeas = 1e300;
  double complem = 1e300;
  // double complem_p = 1e300;
  // double dualgap = 1e300;
  double udotr = 1e300, uor_mu = 1e300;
  double projerr = 1e300;  //, projerr_u = 1e300;
  // double error[6];
  double totalresidual = 1e300, /*totalresidual_mu = 1e300,*/ totalresidual_Jor = 1e300;

  double diff_fixp = 1e300;
  double* diff_fixp_vec = (double*)calloc(n, sizeof(double));

  double* primalConstraint = data->tmp_vault_nd[no_nd++];
  double* dualConstraint = data->tmp_vault_m[no_m++];
  double* complemConstraint = data->tmp_vault_nd[no_nd++];
  // double* complemConstraint_mu = data->tmp_vault_nd[no_nd++];
  double* fixpConstraint = (double*)calloc(n, sizeof(double));
  double* arr_norm_theta = (double*)calloc(max_iter, sizeof(double));
  //  double norm_theta = 1e300;

  double gmm = gmmp1 + gmmp2;
  double norm_f = cblas_dnrm2(m, f, 1);
  double norm_w = cblas_dnrm2(nd, w, 1);

  double* velocity_inv = (double*)calloc(nd, sizeof(double));
  double* d_globalVelocity = (double*)calloc(m, sizeof(double));
  double* d_velocity = (double*)calloc(nd, sizeof(double));
  double* d_reaction = (double*)calloc(nd, sizeof(double));
  double* d_s = (double*)calloc(n, sizeof(double));
  double* s = (double*)calloc(n, sizeof(double));

  double* u_plus_du = data->tmp_vault_nd[no_nd++];   // for Mehrotra
  double* r_plus_dr = data->tmp_vault_nd[no_nd++];   // for Mehrotra
  double* dudr_jprod = data->tmp_vault_nd[no_nd++];  // for Mehrotra
  // double* s_plus_ds = data->tmp_vault_nd[no_n++];    // for Mehrotra
  // double diff_fixp_plus = 1e300;

  double* rhs = options->dWork;
  double* rhs_2 = (double*)calloc(m + 2 * nd + n, sizeof(double));
  double* rhs_tmp = (double*)calloc(m + 2 * nd + n, sizeof(double));
  double* sol = (double*)calloc(m + 2 * nd + n, sizeof(double));

  double* nub_vec = (double*)calloc(n, sizeof(double));

  char fws = ' '; /* finish without scaling */

  /* norm of the residuals of teh second linear system */
  double LS_norm_p = 0.;  // primal feasibility
  double LS_norm_d = 0.;  // dual feaqsibility
  double LS_norm_c = 0.;  // complementarity
  double LS_norm_f = 0.;  // fixed point

  NumericsMatrix *J = NULL, *J_dense = NULL;
  NumericsMatrix *arrow_r = NULL, *arrow_u = NULL;

  long J_nzmax;
  size_t H_nzmax = NM_nnz(H);
  size_t M_nzmax = NM_nnz(M);
  //  size_t arrow_nzmax = 0;
  size_t subdiff_u_nzmax = 0;
  /* For CLASSIFICATION BNRT */
  int nB, nN, nR, nT;
  nB = nN = nR = nT = 0;
  int* setR = (int*)calloc(n, sizeof(int));

  size_t J_nz_captured = 0;

  // For QP2
  double* Qp_u = data->tmp_vault_nd[no_nd++];
  double* Qp_du = data->tmp_vault_nd[no_nd++];
  double* Qpinv_dr = data->tmp_vault_nd[no_nd++];
  double* ududr = data->tmp_vault_nd[no_nd++];
  double* Qp_ududr = data->tmp_vault_nd[no_nd++];

  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO] ==
      SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_YES) {
    numerics_printf_verbose(1, "---- GFC3D - IPM - Problem information");
    numerics_printf_verbose(1, "---- GFC3D - IPM - 1-norm of M = %g norm of f = %g ",
                            NM_norm_1(M), norm_f);
    numerics_printf_verbose(1, "---- GFC3D - IPM - inf-norm of M = %g ", NM_norm_inf(M));

    numerics_printf_verbose(1, "---- GFC3D - IPM - 1-norm of H = %g norm of w = %g ",
                            NM_norm_1(problem->H), norm_w);
    numerics_printf_verbose(1, "---- GFC3D - IPM - inf-norm of H = %g ",
                            NM_norm_inf(problem->H));
    numerics_printf_verbose(1, "---- GFC3D - IPM - M is symmetric = %i ", NM_is_symmetric(M));

    numerics_printf_verbose(1, "---- GFC3D - IPM - M size = (%i, %i) ", M->size0, M->size1);
    numerics_printf_verbose(1, "---- GFC3D - IPM - H size = (%i, %i) ", problem->H->size0,
                            problem->H->size1);
  }

  /* ---- IPM iterations ---- */
  numerics_printf_verbose(1, "problem dimensions n, nd x m: %1i, %6i x %-6i", n, nd, m);
  switch (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM]) {
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_NOSCAL: {
      numerics_printf_verbose(
          -1, "Global friction contact problem - LS solution: 4x4 no scaling\n");
      // arrow_nzmax = (d * 3 - 2) * n;
      subdiff_u_nzmax = 2 * n;
      J_nzmax = M_nzmax + H_nzmax + 2 * (d * 3 - 2) * n + H_nzmax + nd + n + 2 * n + n;
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_QP2: {
      numerics_printf_verbose(1,
                              "Global friction contact problem - LS solution: 4x4 NT scaling "
                              "with Qp2 and Mehrotra\n");
      break;
    }

    default: {
      printf("ERROR\n");
    }
  }

  numerics_printf_verbose(1,
                          "| it  | pinfeas | dinfeas |  |s-ub| |Ma|ui'ri||  |uor|  | prj err "
                          "| barpram |  sigma  ||  alpha  |  |dv|   |  |du|   |  |dr|   |  "
                          "|ds|   | ls prim | ls dual | ls comp | ls fixP |");
  numerics_printf_verbose(1,
                          "-------------------------------------------------------------------"
                          "-------------------------------------------------------------------"
                          "--------------------------------------------");

  double* p2 = data->tmp_vault_nd[no_nd++];
  NumericsMatrix* Qp2 = NULL;
  NumericsMatrix* eye_nd = NM_eye(nd);
  NumericsMatrix* eye_n = NM_eye(n);
  NumericsMatrix* subdiff_u = NULL;

  NumericsMatrix* minus_e = NULL;
  minus_e = NM_create(NM_SPARSE, nd, n);
  size_t minus_e_nzmax = n;
  NM_triplet_alloc(minus_e, minus_e_nzmax);
  for (size_t i = 0; i < n; ++i) {
    NM_entry(minus_e, i * d, i, -1.);
  }

  FILE* iterates = NULL;
  FILE* iterates_2 = NULL;
  FILE* sol_file = NULL;

  char* strToken = NULL;
  char* str = (char*)malloc(200);
  if (problem->name) {
    strcpy(str, problem->name);
  } else {
    strcpy(str, "foo_");
  }
  const char* separators = "/";
  strToken = strtok(str, separators);
  for (int i = 0; i < 5; i++) {
    if (strToken != NULL) strToken = strtok(NULL, separators);
  }

  strToken = strtok(strToken, ".");
  // for(int i=0; i<strlen(strToken); i++)
  // {
  //   if(strToken[i] == '-') strToken[i] = '_';
  // }

  char matlab_name[100];
  snprintf(matlab_name, sizeof(matlab_name), "iterates_Spheres_no_s.m");

  /* writing data in a Matlab file */
  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE]) {
    iterates = fopen(matlab_name, "w");
    printDataProbMatlabFile(M, f, H, w, d, n, m, problem->mu, iterates);
  }

  /* check the full criterion */
  // double norm_q = cblas_dnrm2(m, problem->q, 1);
  // double norm_b = cblas_dnrm2(nd, problem->b, 1);

  // double tmp_barr_param = 0.;
  // double max_uor_2mu = 0.;       //, tmp_uor_2mu = 0.;
  double scale_sub_diff = 0.91;  // 1.05;
  // double max_val_old = 0, scale_sigma = 0, old_sigma = 0;

  //  ComputeErrorGlobalPtr computeError = NULL;
  // computeError = (ComputeErrorGlobalPtr)&gfc3d_compute_error;

  gfc3d_IPM_SNM_compute_starting_point(m, nd, d, n, M, H, f, strToken, reaction, velocity,
                                       globalVelocity, s, data);

  // Reset params
  iteration = 0;
  hasNotConverged = 1;
  pinfeas = dinfeas = complem = udotr = projerr = diff_fixp = totalresidual = 1e300;
  alpha_primal = alpha_dual = 1.;
  barr_param = 1.;
  sigma = 0.1;

  /* ---- Pack all workspace into a single struct for the per-LS-form step functions ---- */
  gfc3d_IPM_SNM_workspace ws = {0};
  ws.m = m;
  ws.nd = nd;
  ws.d = d;
  ws.n = n;
  ws.mp2nd = mp2nd;
  ws.type = type;
  ws.tol = tol;
  ws.max_iter = max_iter;
  ws.sgmp1 = sgmp1;
  ws.sgmp2 = sgmp2;
  ws.sgmp3 = sgmp3;
  ws.gmmp1 = gmmp1;
  ws.gmmp2 = gmmp2;
  ws.gmm = gmm;
  ws.scale_sub_diff = scale_sub_diff;
  ws.problem = problem;
  ws.options = options;
  ws.data = data;
  ws.M = M;
  ws.H_tilde = H_tilde;
  ws.H = H;
  ws.minus_H = minus_H;
  ws.minus_Ht = minus_Ht;
  ws.Ht = Ht;
  ws.minus_M = minus_M;
  ws.eye_nd = eye_nd;
  ws.eye_n = eye_n;
  ws.minus_e = minus_e;
  ws.J = J;
  ws.J_dense = J_dense;
  ws.arrow_r = arrow_r;
  ws.arrow_u = arrow_u;
  ws.subdiff_u = subdiff_u;
  ws.Qp2 = Qp2;
  ws.w = w;
  ws.w_ori = w_ori;
  ws.f = f;
  ws.reaction = reaction;
  ws.velocity = velocity;
  ws.globalVelocity = globalVelocity;
  ws.s = s;
  ws.primalConstraint = primalConstraint;
  ws.dualConstraint = dualConstraint;
  ws.complemConstraint = complemConstraint;
  ws.fixpConstraint = fixpConstraint;
  ws.nub_vec = nub_vec;
  ws.diff_fixp_vec = diff_fixp_vec;
  ws.d_globalVelocity = d_globalVelocity;
  ws.d_velocity = d_velocity;
  ws.d_reaction = d_reaction;
  ws.d_s = d_s;
  ws.u_plus_du = u_plus_du;
  ws.r_plus_dr = r_plus_dr;
  ws.dudr_jprod = dudr_jprod;
  ws.p2 = p2;
  ws.Qp_u = Qp_u;
  ws.Qp_du = Qp_du;
  ws.Qpinv_dr = Qpinv_dr;
  ws.ududr = ududr;
  ws.Qp_ududr = Qp_ududr;
  ws.velocity_inv = velocity_inv;
  ws.rhs = rhs;
  ws.rhs_2 = rhs_2;
  ws.rhs_tmp = rhs_tmp;
  ws.sol = sol;
  ws.arr_norm_theta = arr_norm_theta;
  ws.setR = setR;
  ws.internal_allocation = internal_allocation;
  ws.hasNotConverged = hasNotConverged;
  ws.iteration = iteration;
  ws.pinfeas = pinfeas;
  ws.dinfeas = dinfeas;
  ws.complem = complem;
  ws.udotr = udotr;
  ws.uor_mu = uor_mu;
  ws.projerr = projerr;
  ws.diff_fixp = diff_fixp;
  ws.totalresidual = totalresidual;
  ws.totalresidual_Jor = totalresidual_Jor;
  ws.barr_param = barr_param;
  ws.sigma = sigma;
  ws.alpha_primal = alpha_primal;
  ws.alpha_dual = alpha_dual;
  ws.LS_norm_p = LS_norm_p;
  ws.LS_norm_d = LS_norm_d;
  ws.LS_norm_c = LS_norm_c;
  ws.LS_norm_f = LS_norm_f;
  ws.J_nzmax = J_nzmax;
  ws.J_nz_captured = J_nz_captured;
  ws.subdiff_u_nzmax = subdiff_u_nzmax;
  ws.H_nzmax = H_nzmax;
  ws.M_nzmax = M_nzmax;
  ws.iterates = iterates;
  ws.iterates_2 = iterates_2;
  ws.strToken = strToken;
  ws.str = str;
  ws.fws = fws;

  while (iteration < max_iter) {
    int jacobian_is_nan = 0;
    switch (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM]) {
      /*  Build the Jacobian matrix without reducing the linear system.
       *
       *
       *  In the case where the NT scaling is not used, the matrix to factorize is the
       following:
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
        gfc3d_IPM_SNM_step_NOSCAL(&ws);
        break;

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
        * Case of a non-reduced system. Building the right-hand side related to the first
        system

          with NT scaling
          rhs = -
          [      M*v - H'*r - f      ]  m         dualConstraint
          [            r             ]  nd        complemConstraint    r - 2 mu u^-1
          [   u - H*v - w - phi(s)   ]  nd        primalConstraint
          [         s - |ub|         ]  n
       */
      case SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_QP2:
        gfc3d_IPM_SNM_step_QP2(&ws);
        break;

      default: {
        printf("ERROR\n");
      }
    }

    /* ---- sync main-loop locals from workspace ---- */
    jacobian_is_nan = ws.jacobian_is_nan;
    hasNotConverged = ws.hasNotConverged;
    J = ws.J;
    J_dense = ws.J_dense;
    arrow_r = ws.arrow_r;
    arrow_u = ws.arrow_u;
    subdiff_u = ws.subdiff_u;
    Qp2 = ws.Qp2;
    pinfeas = ws.pinfeas;
    dinfeas = ws.dinfeas;
    complem = ws.complem;
    udotr = ws.udotr;
    uor_mu = ws.uor_mu;
    projerr = ws.projerr;
    diff_fixp = ws.diff_fixp;
    totalresidual = ws.totalresidual;
    totalresidual_Jor = ws.totalresidual_Jor;
    barr_param = ws.barr_param;
    sigma = ws.sigma;
    alpha_primal = ws.alpha_primal;
    alpha_dual = ws.alpha_dual;
    gmm = ws.gmm;
    scale_sub_diff = ws.scale_sub_diff;
    LS_norm_p = ws.LS_norm_p;
    LS_norm_d = ws.LS_norm_d;
    LS_norm_c = ws.LS_norm_c;
    LS_norm_f = ws.LS_norm_f;
    J_nzmax = ws.J_nzmax;
    J_nz_captured = ws.J_nz_captured;

    if (jacobian_is_nan) {
      hasNotConverged = 2;
      if (J) {
        J = NM_free(J);
      }
      break;
    }

    if (NV_isnan(globalVelocity, m) | NV_isnan(velocity, nd) | NV_isnan(reaction, nd) |
        NV_isnan(s, n)) {
      hasNotConverged = 2;
      break;
    }

    if (hasNotConverged == 0 || hasNotConverged >= 2) break;

    numerics_printf_verbose(
        -1,
        "| %3i%c| %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e || %.1e | %.1e | %.1e "
        "| %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |",
        iteration, fws, pinfeas, dinfeas, diff_fixp, udotr, complem, projerr, barr_param,
        sigma, alpha_primal, fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
        fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
        fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]), fabs(d_s[cblas_idamax(n, d_s, 1)]),
        LS_norm_p, LS_norm_d, LS_norm_c, LS_norm_f);

    classify_BNRT(velocity, reaction, nd, n, &nB, &nN, &nR, &nT);
    // printf("B/N/R/T = %d %d %d %d\n", nB, nN, nR, nT);

    if (J_dense) {
      J_dense = NM_free(J_dense);
    }

    iteration++;

  }  // while loop

  /* Checking strict complementarity */
  classify_BNRT(velocity, reaction, nd, n, &nB, &nN, &nR, &nT);
  if (nT > 0)
    printf("Ratio of Strict complementarity solutions: %4zu / %4zu = %4.2f \t %4i %4i %4i\n",
           n - nT, n, (double)(n - nT) / (double)n, nB, nN, nR);
  else
    printf("Strict complementarity satisfied: %4i %4i %4i\n", nB, nN, nR);

  // Store solution into file
  int save_sol_point = 0;
  if (save_sol_point) {
    sol_file = fopen("sol_data.res", "w");
    // store v
    for (size_t i = 0; i < m; i++) {
      fprintf(sol_file, "%8.20e ", globalVelocity[i]);
    }
    fprintf(sol_file, "\n");

    // store u
    for (size_t i = 0; i < nd; i++) {
      fprintf(sol_file, "%8.20e ", velocity[i]);
    }
    fprintf(sol_file, "\n");

    // store r
    for (size_t i = 0; i < nd; i++) {
      fprintf(sol_file, "%8.20e ", reaction[i]);
    }
    fprintf(sol_file, "\n");

    // store s
    for (size_t i = 0; i < n; i++) {
      fprintf(sol_file, "%8.20e ", s[i]);
    }
    fprintf(sol_file, "\n");

    fclose(sol_file);
  }

  /* ----- cleanup: restore scaling and free memory ------ */
  gfc3d_IPM_SNM_cleanup_data cd = {0};
  cd.problem = problem;
  cd.options = options;
  cd.data = data;
  cd.internal_allocation = internal_allocation;
  cd.hasNotConverged = hasNotConverged;
  cd.nd = nd;
  cd.n = n;
  cd.velocity = velocity;
  cd.reaction = reaction;
  cd.totalresidual = totalresidual;
  cd.projerr = projerr;
  cd.M = M;
  cd.M_original = M_original;
  cd.H_tilde = H_tilde;
  cd.H = H;
  cd.minus_H = minus_H;
  cd.minus_Ht = minus_Ht;
  cd.Ht = Ht;
  cd.minus_M = minus_M;
  cd.eye_nd = eye_nd;
  cd.eye_n = eye_n;
  cd.subdiff_u = subdiff_u;
  cd.minus_e = minus_e;
  cd.arrow_r = arrow_r;
  cd.arrow_u = arrow_u;
  cd.J = J;
  cd.J_dense = J_dense;
  cd.w_ori = w_ori;
  cd.diff_fixp_vec = diff_fixp_vec;
  cd.fixpConstraint = fixpConstraint;
  cd.arr_norm_theta = arr_norm_theta;
  cd.d_globalVelocity = d_globalVelocity;
  cd.d_velocity = d_velocity;
  cd.d_reaction = d_reaction;
  cd.d_s = d_s;
  cd.s = s;
  cd.velocity_inv = velocity_inv;
  cd.rhs_2 = rhs_2;
  cd.rhs_tmp = rhs_tmp;
  cd.sol = sol;
  cd.nub_vec = nub_vec;
  cd.setR = setR;
  cd.str = str;
  cd.blk_num_name = blk_num_name;
  iterates = NULL;
  iterates_2 = NULL;
  cd.iterates = iterates;
  cd.iterates_2 = iterates_2;

  gfc3d_IPM_SNM_cleanup(&cd, info);
}

void gfc3d_IPM_SNM_set_default(SolverOptions* options) {
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO] =
      SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_NO;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] =
      SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_NOSCAL;
  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] =
  // SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_QP2;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE] = 0;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] =
      SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES;

  options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000000;
  options->dparam[SICONOS_DPARAM_TOL] = 1e-14;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1] = 1e-10;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2] = 3.;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3] = 1.;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1] = 0.9;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2] = 0.09;  // 0.095
}

static int gfc3d_ipm_snm_init_wrap(void* problem, SolverOptions* options) {
  /* set_default is now called by the registry before init_wrap */
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int gfc3d_ipm_snm_solve_wrap(void* problem, double* reaction, double* velocity,
                                    double* globalVelocity, SolverOptions* options) {
  int info = NUMERICS_OK;
  gfc3d_IPM((GlobalFrictionContactProblem*)problem, reaction, velocity, globalVelocity, &info,
            options);
  return info;
}

static void gfc3d_ipm_snm_free_wrap(void* problem, SolverOptions* options) {
  /* Cleanup if needed */
  (void)problem;
  (void)options;
}
REGISTER_SOLVER_3VAR(
    GFC3D_IPM_SNM, "GFC3D_IPM_SNM",
    "Interior Point Method with Smoothing and Newton for Global 3D Friction Contact",
    gfc3d_ipm_snm_init_wrap, gfc3d_ipm_snm_solve_wrap, gfc3d_ipm_snm_free_wrap, NULL,
    gfc3d_IPM_SNM_set_default, /* set_default */
    500,                       /* default_max_iter - from set_default */
    1e-10,                     /* default_tol - from set_default */
    0 /* is_local_solver */)
