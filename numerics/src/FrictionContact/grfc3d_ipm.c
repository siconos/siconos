/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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
#include "grfc3d_Solvers.h"             // for GRFCProb, SolverOpt, Friction_cst, grfc3d_...
#include "grfc3d_compute_error.h"       // for grfc3d_compute_error
#include "SiconosLapack.h"

#include "SparseBlockMatrix.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#include "numerics_verbose.h"
#include "NumericsVector.h"
#include "float.h"
#include "JordanAlgebra.h"              // for JA functions

#include "NumericsSparseMatrix.h"
#include "NumericsMatrix.h"

#include <time.h>                       // for clock()
#include "gfc3d_ipm.h"                  // for primalResidual, dualResidual, ...
#include "grfc3d_ipm.h"                 // for dnrm2sqrl
#include "cs.h"

/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "siconos_debug.h"


typedef struct
{
  double * globalVelocity;  // v
  double * velocity;        // u
  double * reaction;        // r
}
  IPM_point;

typedef struct
{
  double * velocity_1; // velocity_1 = (t, u_bar)
  double * velocity_2; // velocity_2 = (t_prime, u_tilde)
  double * reaction_1; // reaction_1 = (r0, r_bar)
  double * reaction_2; // reaction_2 = (r0, r_tilde)
  double * t;
  double * t_prime;
}
  IPM_grfc3d_point;

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
  IPM_point* starting_point;        // initial point
  IPM_point* original_point;        // original point which is not changed by the matrix P_mu
  IPM_grfc3d_point* grfc3d_point;

  /* change of variable matrix */
  IPM_change_of_variable* P_mu;

  /* initial internal solver parameters */
  IPM_internal_params* internal_params;

  /* problem parameters related to their size
   * tmp_vault_m[0] = dualConstraint
   * tmp_vault_m[1] = gv_plus_dgv
   * tmp_vault_nd[0] = w
   * tmp_vault_nd[1] = primalConstraint
   * tmp_vault_nd[2] = vr_jprod
   * tmp_vault_nd[3] = v_plus_dv = u + alpha_p * du
   * tmp_vault_nd[4] = r_plus_dr = r + alpha_d * dr
   * tmp_vault_nd[0] = velocity_1t = F_1 * u_1
   * tmp_vault_nd[0] = velocity_2t = F_2 * u_2
   * tmp_vault_nd[0] = d_velocity_1t = F_1 * du_1
   * tmp_vault_nd[0] = d_velocity_2t = F_2 * du_2
   * tmp_vault_nd[0] = d_reaction_1t = Finv_1 * dr_1
   * tmp_vault_nd[0] = d_reaction_2t = Finv_2 * dr_2
   * tmp_vault_n_dminus2[0] = complemConstraint_1
   * tmp_vault_n_dminus2[1] = complemConstraint_2
   * tmp_vault_n_dminus2[2] = dvdr_jprod_1
   * tmp_vault_n_dminus2[3] = dvdr_jprod_2
   * tmp_vault_n_dminus2[4] = p_bar
   * tmp_vault_n_dminus2[5] = p_tilde


   * tmp_vault_n[0] = d_t
   * tmp_vault_n[1] = d_t_prime
   */
  double **tmp_vault_m;
  double **tmp_vault_nd;
  double **tmp_vault_n_dminus2;
  double **tmp_vault_n;
}
  Grfc3d_IPM_data;







/* ------------------------- Helper functions ------------------------------ */
/** Return a speacial sub-vector such that
 * the 1st element is always taken and
 * so do from i-th to j-th elements,
 * starting index is 1
 */
static void extract_vector(const double * const vec, const unsigned int vecSize, const int varsCount, const unsigned int i, const unsigned int j, double * out)
{
  assert(vec);
  assert(i >= 1);
  assert(i <= j);
  assert(out);

  unsigned int vec_dim = (int)(vecSize / varsCount);
  assert(j <= vec_dim);

  unsigned int out_dim = vec_dim - 2;
  assert(out_dim > 0);

  size_t posX = 0;
  size_t posY = 0;

  for(size_t k = 0; k < varsCount; k++)
  {
    out[posX++] = vec[posY];
    for(size_t l = i-1; l < j; l++)
    {
      out[posX++] = vec[posY + l];
    }
    posY += vec_dim;
  }
}


/** Compute a block matrix J of form
 *      |  1     0     0                 1     0     0                    ...  |
 *      |                                                                      |
 *      |  0     1     0                 0     0     0                    ...  |
 *      |                                                                      |
 *      |  0     0     1                 0     0     0                    ...  |
 *      |                                                                      |
 *      |  0     0     0                 0     1     0                    ...  |
 *      |                                                                      |
 *      |  0     0     0                 0     0     1                    ...  |
 *  J = |                                                                      |
 *      |                 .     .     .                  .     .     .    ...  |
 *      |                                                                      |
 *      |                 .     .     .                  .     .     .    ...  |
 *      |                                                                      |
 *      |                 .     .     .                  .     .     .    ...  |
 *      |                                                                      |
 *      |                 .     .     .                  .     .     .    ...  |
 *      |                                                                      |
 *      |                 .     .     .                  .     .     .    ...  |
 *      |                                                                      |
 *      | ...   ...   ...   ...   ...   ...   ...   ...   ...   ...   ... ...  |
 */
static NumericsMatrix* compute_J_matrix(const unsigned int varsCount)
{
  assert(varsCount > 0);

  NumericsMatrix * J = NM_create(NM_SPARSE, 5*varsCount, 3*varsCount*2);
  NumericsMatrix * J_1 = NM_create(NM_SPARSE, 5, 3);
  NumericsMatrix * J_2 = NM_create(NM_SPARSE, 5, 3);

  long J_nzmax = 3*2*varsCount;
  long J_1_nzmax = 3;
  long J_2_nzmax = 3;

  NM_triplet_alloc(J, J_nzmax);
  NM_triplet_alloc(J_1, J_1_nzmax);
  NM_triplet_alloc(J_2, J_1_nzmax);

  J->matrix2->origin = NSM_TRIPLET;
  J_1->matrix2->origin = NSM_TRIPLET;
  J_2->matrix2->origin = NSM_TRIPLET;

  NM_insert(J_1, NM_eye(3), 0, 0);
  NM_insert(J_2, NM_eye(1), 0, 0);
  NM_insert(J_2, NM_eye(2), 3, 1);

  for (unsigned i = 0; i < varsCount; i++)
  {
    NM_insert(J, J_1, i*5, i*3);
    NM_insert(J, J_2, i*5, varsCount*3 + i*3);
  }

  return J;
}


/*
 * Return members of NT matrix F family (Tutuncu, Toh and Todd, Math.Prog 2003, pp. 195-196)
 * f[in,out]: element of NT matrix matrix F
 * wf[in,out]: a coefficient related to f
 * F2[in,out] = F^2        Finv[in,out] = F^-1      Finv2[in,out] = F^-2
 * A member computed must be allocated before.
 * If an input is NULL, then this member will not be computed.
 */
static void family_of_F(const double * const x, const double * const z, const unsigned int vecSize, const size_t varsCount,
                      double * f, float_type * wf, NumericsMatrix * F, NumericsMatrix * Finv, NumericsMatrix * F2, NumericsMatrix * Finv2)
{
  size_t dimension = (size_t)(vecSize / varsCount);
  int f_NULL, wf_NULL;
  f_NULL = 0; wf_NULL = 0;
  NumericsMatrix *F_block, *Finv_block, *F2_block, *Finv2_block;

  if (!f)   // f is always allocated to use
  {
    f_NULL = 1;
    f = (double*)calloc(vecSize, sizeof(double));
  }
  if (!wf)
  {
    wf_NULL = 1;
    wf = (float_type*)calloc(varsCount, sizeof(float_type));
  }
  if (F)
  {
    F_block = NM_create(NM_DENSE, dimension, dimension);
  }
  if (Finv)
  {
    Finv_block = NM_create(NM_DENSE, dimension, dimension);
  }
  if (F2)
  {
    F2_block = NM_create(NM_DENSE, dimension, dimension);
  }
  if (Finv2)
  {
    Finv2_block = NM_create(NM_DENSE, dimension, dimension);
  }

  float_type gamx, gamz, gamf, w, w2;
  for(unsigned int i = 0; i < vecSize; i += dimension)
  {
    gamx = gammal(x+i, dimension);
    gamz = gammal(z+i, dimension);
    w = sqrtl(gamz/gamx);

    wf[(int)(i/dimension)] = w;
    if (F2 || Finv2)
      w2 = gamz/gamx;

    f[i] = z[i]/w + w*x[i];
    for(unsigned int j = 1; j < dimension; ++j)
    {
      f[i+j] = z[i+j]/w - w*x[i+j];
    }

    gamf = gammal(f+i, dimension);
    for(size_t j = 0; j < dimension; ++j)
    {
      f[i+j] = f[i+j]/gamf;
    }

    // Construct matrix F
    if (F)
    {
      F_block->matrix0[0] = f[i]*w;
      for(size_t j = 1; j < dimension; ++j)
      {
        F_block->matrix0[j] = f[i+j]*w;
        F_block->matrix0[j*dimension] = F_block->matrix0[j];
      }
      for(size_t j = 1; j < dimension; ++j)
      {
        for(size_t k = 1; k < dimension; ++k)
          F_block->matrix0[j+dimension*k] = (j==k) ? (1+f[i+j]*f[i+k]/(1+f[i]))*w : f[i+j]*f[i+k]/(1+f[i])*w;
      }
      NM_insert(F, F_block, i, i);
    }

    // Construct matrix F^-1
    if (Finv)
    {
      Finv_block->matrix0[0] = f[i]/w;
      for(size_t j = 1; j < dimension; ++j)
      {
        Finv_block->matrix0[j] = -f[i+j]/w;
        Finv_block->matrix0[j*dimension] = Finv_block->matrix0[j];
      }
      for(size_t j = 1; j < dimension; ++j)
      {
        for(size_t k = 1; k < dimension; ++k)
          Finv_block->matrix0[j+dimension*k] = (j==k) ? (1+f[i+j]*f[i+k]/(1+f[i]))/w : f[i+j]*f[i+k]/(1+f[i])/w;
      }
      NM_insert(Finv, Finv_block, i, i);
    }

    // Construct matrix F^2
    if (F2)
    {
      F2_block->matrix0[0] = w2*cblas_ddot(dimension, f+i, 1, f+i, 1);
      for(size_t j = 1; j < dimension; ++j)
      {
        F2_block->matrix0[j] = 2*f[i]*f[i+j]*w2;
        F2_block->matrix0[j*dimension] = F2_block->matrix0[j];
      }
      for(size_t j = 1; j < dimension; ++j)
      {
        for(size_t k = 1; k < dimension; ++k)
          F2_block->matrix0[j+dimension*k] = (j==k) ? (1+2*f[i+j]*f[i+k])*w2 : 2*f[i+j]*f[i+k]*w2;
      }
      NM_insert(F2, F2_block, i, i);
    }

    // Construct matrix F^-2
    if (Finv2)
    {
      Finv2_block->matrix0[0] = cblas_ddot(dimension, f+i, 1, f+i, 1)/w2;
      for(size_t j = 1; j < dimension; ++j)
      {
        Finv2_block->matrix0[j] = -2*f[i]*f[i+j]/w2;
        Finv2_block->matrix0[j*dimension] = Finv2_block->matrix0[j];
      }
      for(size_t j = 1; j < dimension; ++j)
      {
        for(size_t k = 1; k < dimension; ++k)
          Finv2_block->matrix0[j+dimension*k] = (j==k) ? (1+2*f[i+j]*f[i+k])/w2 : 2*f[i+j]*f[i+k]/w2;
      }
      NM_insert(Finv2, Finv2_block, i, i);
    }
  }

  if (f_NULL > 0)
  {
    free(f); f = NULL;
  }
  if (wf_NULL > 0)
  {
    free(wf); wf = NULL;
  }
  if (F)
  {
    NM_clear(F_block); free(F_block);
  }
  if (Finv)
  {
    NM_clear(Finv_block); free(Finv_block);
  }
  if (F2)
  {
    NM_clear(F2_block); free(F2_block);
  }
  if (Finv2)
  {
    NM_clear(Finv2_block); free(Finv2_block);
  }
}




/* Return the matrix P^-1 where P is the matrix satisfying Jac = P*P' */
static  NumericsMatrix *  Pinv(const double * const f, const double * const g, const float_type * const wf, const float_type * const wg, const unsigned int vecSize, const size_t varsCount)
{
  size_t dim = (size_t)(vecSize / varsCount); // dim must be 3
  size_t d5 = dim+2;  // d5 must be 5

  NumericsMatrix * out = NM_create(NM_SPARSE, 5*varsCount, 5*varsCount);
  NM_triplet_alloc(out, (5+2*(2*2))*varsCount);

  NumericsMatrix * out15 = NM_create(NM_DENSE, 1, 5);
  NumericsMatrix * out22 = NM_create(NM_DENSE, 2, 2);
  double * othor = (double*)calloc(2, sizeof(double));

  float_type coef, coef_tmp;
  coef = 1.0; coef_tmp = 1.0;

  for(size_t i = 0; i < varsCount; i++)
  {
    // For matrix 1x5 at (0,0)
    out15->matrix0[0] = 1./sqrtl( dnrm2sqrl(dim,f+i*dim)/(wf[i]*wf[i]) + dnrm2sqrl(dim,g+i*dim)/(wg[i]*wg[i])
                                - 4*f[i*dim]*f[i*dim]*dnrm2sqrl(dim-1,f+i*dim+1)/(wf[i]*wf[i]*dnrm2sqrl(dim,f+i*dim))
                                - 4*g[i*dim]*g[i*dim]*dnrm2sqrl(dim-1,g+i*dim+1)/(wg[i]*wg[i]*dnrm2sqrl(dim,g+i*dim)) );

    coef = 2*f[i*dim]*out15->matrix0[0]/dnrm2sqrl(dim,f+i*dim);
    for(size_t k = 1; k < dim; k++)
    {
      out15->matrix0[k] = coef*f[i*dim+k];
    }

    coef = 2*g[i*dim]*out15->matrix0[0]/dnrm2sqrl(dim,g+i*dim);
    for(size_t k = 1; k < dim; k++)
    {
      out15->matrix0[k+2] = coef*g[i*dim+k];
    }
    NM_insert(out, out15, i*d5, i*d5);

    // For matrix 2x2 at (1,1)
    coef = wf[i]/dnrm2sqrl(dim-1,f+i*dim+1);
    coef_tmp = 1./dnrm2l(dim,f+i*dim);
    othor[0] = -f[i*dim+2];
    othor[1] = f[i*dim+1];
    size_t pos=0;
    for(size_t k = 1; k < dim; k++)
    {
      for(size_t l = 1; l < dim; l++)
      {
        out22->matrix0[pos++] = coef*(coef_tmp*f[i*dim+k]*f[i*dim+l]+othor[k-1]*othor[l-1]);
      }
    }
    NM_insert(out, out22, i*d5+1, i*d5+1);

    // For matrix 2x2 at (3,3)
    coef = wg[i]/dnrm2sqrl(dim-1,g+i*dim+1);
    coef_tmp = 1./dnrm2l(dim,g+i*dim);
    othor[0] = -g[i*dim+2];
    othor[1] = g[i*dim+1];
    pos=0;
    for(size_t k = 1; k < dim; k++)
    {
      for(size_t l = 1; l < dim; l++)
      {
        out22->matrix0[pos++] = coef*(coef_tmp*g[i*dim+k]*g[i*dim+l]+othor[k-1]*othor[l-1]);
      }
    }
    NM_insert(out, out22, i*d5+3, i*d5+3);
  }

  NM_clear(out15); free(out15);
  NM_clear(out22); free(out22);
  free(othor);
  return out;
}







/* Return the matrix P^-1*H where P is the matrix satisfying Jac = P*P' */ // TO DO HERE
static  NumericsMatrix *  multiply_PinvH(const double * const f, const double * const g, const float_type * const wf, const float_type * const wg, const unsigned int vecSize, const size_t varsCount, NumericsMatrix *H)
{
  size_t dim = (size_t)(vecSize / varsCount); // dim must be 3
  size_t d5 = dim+2;  // d5 must be 5

  NumericsMatrix * out = NM_new();

  if(H->storageType != NM_SPARSE)
  {
    fprintf(stderr, "Numerics, GRFC3D IPM, PinvH failed, only accept for NM_SPARSE of H.\n");
    exit(EXIT_FAILURE);
  }


  CSparseMatrix* H_csc = NM_csc(H);
  CSparseMatrix* out_csc  = cs_spalloc (H->size0, H->size1, H_csc->nzmax , 1, 0) ;        /* allocate result */

  CS_INT *Hp, *Hi ;
  Hp = H_csc->p ; Hi = H_csc->i ;

  CS_INT  *outp, *outi ;
  outp = out_csc->p ;

  CS_ENTRY *Hx, *outx ;
  Hx = H_csc->x ;

  CS_INT nz = 0;

  float_type P0=1.0, coef=1.0, coef_tmp=1.0;
  double * othor = (double*)calloc(2, sizeof(double));

  bool multiplied = 0;

  for (CS_INT k = 0 ; k < H->size1 ; k++) // traverse all cols of H
  {
    outp[k] = nz ;                   /* column k of out starts here */

    /* reallocate if needed */
    if (nz + H->size1 > out_csc->nzmax && !cs_sprealloc (out_csc, 2*(out_csc->nzmax)+H->size1))
    {
      return NULL;             /* out of memory */
    }
    outi = out_csc->i ; outx = out_csc->x ;   /* out->i and out->x may be reallocated */

    for(size_t i = 0; i < varsCount; i++) // traverse all blocks[5x5] of Pinv
    {
      for (size_t j = 0; j < d5; j++)  // traverse all rows-block of Pinv
      {
        /* multiplication and storage */
        if (j==0) // 1st row of P^-1
        {
          outx[nz] = 0.;
          for (CS_INT p = Hp [k] ; p < Hp [k+1] ; p++)  // traverse all existential rows of H
          {
            if(i*d5<=Hi[p] && Hi[p]<(i+1)*d5) // rows of H such that they belongs to each block of P^-1
            {
              multiplied = 1;
              P0 = 1./sqrtl( dnrm2sqrl(dim,f+i*dim)/(wf[i]*wf[i]) + dnrm2sqrl(dim,g+i*dim)/(wg[i]*wg[i])
                           - 4*f[i*dim]*f[i*dim]*dnrm2sqrl(dim-1,f+i*dim+1)/(wf[i]*wf[i]*dnrm2sqrl(dim,f+i*dim))
                           - 4*g[i*dim]*g[i*dim]*dnrm2sqrl(dim-1,g+i*dim+1)/(wg[i]*wg[i]*dnrm2sqrl(dim,g+i*dim)) );

              if (Hi[p] == i*d5) {outx[nz] += P0*Hx[p];}

              if (i*d5+1<=Hi[p] && Hi[p]<=i*d5+2)
              {
                coef = 2.*f[i*dim]*P0/dnrm2sqrl(dim,f+i*dim);
                if (Hi[p] == i*d5+1) {outx[nz] += coef*f[i*dim+1]*Hx[p];}
                if (Hi[p] == i*d5+2) {outx[nz] += coef*f[i*dim+2]*Hx[p];}
              }

              if (i*d5+3<=Hi[p] && Hi[p]<=i*d5+4)
              {
                coef = 2.*g[i*dim]*P0/dnrm2sqrl(dim,g+i*dim);
                if (Hi[p] == i*d5+3) {outx[nz] += coef*g[i*dim+1]*Hx[p];}
                if (Hi[p] == i*d5+4) {outx[nz] += coef*g[i*dim+2]*Hx[p];}
              }
              // printf("\n\ni = %zu, k = %ld, j= %zu, p = %ld; \nP0 = %3.20Lf, Hx[p] = %3.20e; \nnz = %ld, outi[nz] = %lu, outp[k] = %lu, outx[nz] = %3.20e\n",
                      // i, k, j, p, P0, Hx[p], nz, j+i*d5, outp[k], outx[nz]);
            }
          } // end rows of H
          if (multiplied)
          {
            outi[nz++] = j+i*d5; multiplied = 0;
// printf("\n\ni = %zu, k = %ld, j= %zu; \nHp [k] = %li, Hp [k+1] = %li; \nnz = %ld, outi[nz] = %lu, outp[k] = %lu\n", i, k, j, Hp [k], Hp [k+1], nz-1, j+i*d5, outp[k]);
          }
        } // end 1st row of P^-1


        else if (j==1 || j==2) // A^-1/2 of P^-1
        {
          outx[nz] = 0.;
          for (CS_INT p = Hp [k] ; p < Hp [k+1] ; p++)  // traverse all existential rows of H
          {
            if (i*d5+1<=Hi[p] && Hi[p]<=i*d5+2) // get the rows of H that belongs to A^-1/2 of each block
            {
              multiplied = 1;
              coef = wf[i]/dnrm2sqrl(dim-1,f+i*dim+1);
              coef_tmp = 1./dnrm2l(dim,f+i*dim);
              othor[0] = -f[i*dim+2];
              othor[1] = f[i*dim+1];

              if(Hi[p]==i*d5+1)
              {
                if (j==1) outx[nz] += coef*(coef_tmp*f[i*dim+1]*f[i*dim+1]+othor[0]*othor[0])*Hx[p]; // 1st row of A^-1/2
                if (j==2) outx[nz] += coef*(coef_tmp*f[i*dim+2]*f[i*dim+1]+othor[1]*othor[0])*Hx[p]; // 2nd row
              }

              if(Hi[p]==i*d5+2)
              {
                if (j==1) outx[nz] += coef*(coef_tmp*f[i*dim+1]*f[i*dim+2]+othor[0]*othor[1])*Hx[p]; // 1st row of A^-1/2
                if (j==2) outx[nz] += coef*(coef_tmp*f[i*dim+2]*f[i*dim+2]+othor[1]*othor[1])*Hx[p]; // 2nd row
              }
              // printf("\n\ni = %zu, k = %ld, j= %zu, p = %ld; \ncoef = %3.20Lf, coef_tmp = %3.20Lf, Hx[p] = %3.20e; \nnz = %ld, outi[nz] = %lu, outp[k] = %lu, outx[nz] = %3.20e\n",
              //         i, k, j, p, coef, coef_tmp, Hx[p], nz, j+i*d5, outp[k], outx[nz]);
            }
          } // end rows of H
          if (multiplied) { outi[nz++] = j+i*d5; multiplied = 0;}
        } // end A^-1/2


        else if (j==3 || j==4) // C^-1/2 of P^-1
        {
          outx[nz] = 0.;
          for (CS_INT p = Hp [k] ; p < Hp [k+1] ; p++)  // traverse all existential rows of H
          {
            if (i*d5+3<=Hi[p] && Hi[p]<=i*d5+4) // get the rows of H that belongs to C^-1/2 of each block
            {
              multiplied = 1;
              coef = wg[i]/dnrm2sqrl(dim-1,g+i*dim+1);
              coef_tmp = 1./dnrm2l(dim,g+i*dim);
              othor[0] = -g[i*dim+2];
              othor[1] = g[i*dim+1];

              if(Hi[p]==i*d5+3)
              {
                if (j==3) outx[nz] += coef*(coef_tmp*g[i*dim+1]*g[i*dim+1]+othor[0]*othor[0])*Hx[p]; // 1st row of C^-1/2
                if (j==4) outx[nz] += coef*(coef_tmp*g[i*dim+2]*g[i*dim+1]+othor[1]*othor[0])*Hx[p]; // 2nd row
              }

              if(Hi[p]==i*d5+4)
              {
                if (j==3) outx[nz] += coef*(coef_tmp*g[i*dim+1]*g[i*dim+2]+othor[0]*othor[1])*Hx[p]; // 1st row of C^-1/2
                if (j==4) outx[nz] += coef*(coef_tmp*g[i*dim+2]*g[i*dim+2]+othor[1]*othor[1])*Hx[p]; // 2nd row
              }
            }
          } // end rows of H
          if (multiplied) { outi[nz++] = j+i*d5; multiplied = 0;}
        } // end C^-1/2
      } // end traverse all rows-block[5x5] of Pinv
    } // end traverse all blocks[5x5] of Pinv
  } // end traverse all cols of H

  outp[H->size1] = nz ;
  cs_sprealloc (out_csc, 0) ;

  out->storageType = H->storageType;
  numericsSparseMatrix(out)->csc = out_csc;
  out->size0 = (int)out->matrix2->csc->m;
  out->size1 = (int)out->matrix2->csc->n;
  numericsSparseMatrix(out)->origin = NSM_CSC;

  return out;
}









/* Return the matrix P^-1'*x where P is the matrix satisfying Jac = P*P' */
static  void  PinvTx(const double * const f, const double * const g, const float_type * const wf, const float_type * const wg, const unsigned int vecSize, const size_t varsCount, const double * const x, double * out)
{
  size_t dim = (size_t)(vecSize / varsCount);
  assert(dim == 3);  // dim must be 3
  size_t d5 = dim+2;  // d5 must be 5
  assert(x);

  float_type P0=1., coef=1., coef_tmp=1., ddot_1=1., ddot_2=1.;
  size_t pos=0;
  double * othor = (double*)calloc(2, sizeof(double));

  for(size_t i = 0; i < varsCount; i++)
  {
    // For x_0
    P0 = 1./sqrtl( dnrm2sqrl(dim,f+i*dim)/(wf[i]*wf[i]) + dnrm2sqrl(dim,g+i*dim)/(wg[i]*wg[i])
                                - 4*f[i*dim]*f[i*dim]*dnrm2sqrl(dim-1,f+i*dim+1)/(wf[i]*wf[i]*dnrm2sqrl(dim,f+i*dim))
                                - 4*g[i*dim]*g[i*dim]*dnrm2sqrl(dim-1,g+i*dim+1)/(wg[i]*wg[i]*dnrm2sqrl(dim,g+i*dim)) );
    out[i*d5] = x[i*d5]*P0;


    // For x_bar
    coef = 2.*f[i*dim]*P0/dnrm2sqrl(dim,f+i*dim);
    for(size_t k = 1; k < dim; k++)
    {
      out[i*d5+k] = coef*f[i*dim+k]*x[i*d5];
    }

    coef = wf[i]/dnrm2sqrl(dim-1,f+i*dim+1);
    coef_tmp = 1./dnrm2l(dim,f+i*dim);
    othor[0] = -f[i*dim+2];
    othor[1] = f[i*dim+1];

    for(size_t k = 1; k < dim; k++)
    {
      ddot_1 = 0.; ddot_2 = 0.;
      for(size_t l = 1; l < dim; l++)
      {
        ddot_1 += f[i*dim+l]*x[i*d5+l];        // Compute f_bar'*x_bar
        ddot_2 += othor[l-1]*x[i*d5+l];  // Compute othor'*x_bar
      }

      out[i*d5+k] += coef*(coef_tmp*ddot_1*f[i*dim+k]+ddot_2*othor[k-1]);
    }


    // For x_tilde
    coef = 2.*g[i*dim]*P0/dnrm2sqrl(dim,g+i*dim);
    for(size_t k = 1; k < dim; k++)
    {
      out[i*d5+k+2] = coef*g[i*dim+k]*x[i*d5];
    }

    coef = wg[i]/dnrm2sqrl(dim-1,g+i*dim+1);
    coef_tmp = 1./dnrm2l(dim,g+i*dim);
    othor[0] = -g[i*dim+2];
    othor[1] = g[i*dim+1];

    for(size_t k = 1; k < dim; k++)
    {
      ddot_1 = 0.; ddot_2 = 0.;
      for(size_t l = 1; l < dim; l++)
      {
        ddot_1 += g[i*dim+l]*x[i*d5+l+2];          // Compute g_bar'*x_bar
        ddot_2 += othor[l-1]*x[i*d5+l+2];  // Compute othor'*x_bar
      }
      out[i*d5+k+2] += coef*(coef_tmp*ddot_1*g[i*dim+k]+ddot_2*othor[k-1]);
    }
  }

  free(othor);
}







/* Return the matrix (J*Q^-2)'*x */
static  void JQinv2Tx(const double * const f, const double * const g, const float_type * const wf, const float_type * const wg, const unsigned int vecSize, const size_t varsCount, const double * const x, double * out)
{
  size_t dim = (size_t)(vecSize / varsCount); // dim must be 3
  assert(dim == 3);
  size_t d5 = dim+2;  // d5 must be 5
  size_t n_d3 = varsCount*dim;  // n_d3 = x*dim


  float_type coef = 1., w2f = 1., w2g = 1., ddot = 1.;

  for(size_t i = 0; i < varsCount; i++)
  {
    w2f = wf[i]*wf[i];
    // For a*x0 + b'*x_bar
    out[i*dim] = dnrm2sqrl(dim,f+i*dim)*x[i*d5]/w2f; // a*x0

    coef = -2.*f[i*dim]/w2f;
    for(size_t k = 1; k < dim; k++)
    {
      out[i*dim] += coef*f[i*dim+k]*x[i*d5+k]; // + b'*x_bar
    }

    // For x0*b + A*x_bar
    ddot = 0.;
    for(size_t l = 1; l < dim; l++)
    {
      ddot += f[i*dim+l]*x[i*d5+l]; // Compute f_bar'*x_bar
    }
    for(size_t k = 1; k < dim; k++)
    {
      out[i*dim+k] = (x[i*d5+k]+2.*(ddot-f[i*dim]*x[i*d5])*f[i*dim+k])/w2f;
    }
  }

  for(size_t i = 0; i < varsCount; i++)
  {
    w2g = wg[i]*wg[i];
    // For c*x0 + d'*x_tilde
    out[i*dim+n_d3] = dnrm2sqrl(dim,g+i*dim)*x[i*d5]/w2g; // c*x0

    coef = -2.*g[i*dim]/w2g;
    for(size_t k = 1; k < dim; k++)
    {
      out[i*dim+n_d3] += coef*g[i*dim+k]*x[i*d5+k+2]; // + d'*x_tilde
    }


    // For x0*d + C*x_bar
    ddot = 0.;
    for(size_t l = 1; l < dim; l++)
    {
      ddot += g[i*dim+l]*x[i*d5+l+2]; // Compute g_bar'*x_tilde
    }
    for(size_t k = 1; k < dim; k++)
    {
      out[i*dim+k+n_d3] = (x[i*d5+k+2]+2.*(ddot-g[i*dim]*x[i*d5])*g[i*dim+k])/w2g;
    }
  }
}






/* Return the matrix Q^-2*x */
static  void Qinv2x(const double * const f, const double * const g, const float_type * const wf, const float_type * const wg, const unsigned int vecSize, const size_t varsCount, const double * const x, double * out)
{
  size_t dim = (size_t)(vecSize / varsCount); // dim must be 3
  assert(dim == 3);
  size_t d5 = dim+2;  // d5 must be 5
  size_t d6 = dim+3;  // d6 must be 6
  size_t n_d3 = varsCount*dim;  // n_d3 = n*dim


  float_type coef = 1., w2f = 1., w2g = 1., ddot = 1.;

  for(size_t i = 0; i < varsCount; i++)
  {
    w2f = wf[i]*wf[i];
    // For a*x0 + b'*x_bar
    out[i*dim] = dnrm2sqrl(dim,f+i*dim)*x[i*dim]/w2f; // a*x0

    coef = -2.*f[i*dim]/w2f;
    for(size_t k = 1; k < dim; k++)
    {
      out[i*dim] += coef*f[i*dim+k]*x[i*dim+k]; // + b'*x_bar
    }

    // For x0*b + A*x_bar
    ddot = 0.;
    for(size_t l = 1; l < dim; l++)
    {
      ddot += f[i*dim+l]*x[i*dim+l]; // Compute f_bar'*x_bar
    }
    for(size_t k = 1; k < dim; k++)
    {
      out[i*dim+k] = (x[i*dim+k]+2.*(ddot-f[i*dim]*x[i*dim])*f[i*dim+k])/w2f;
    }
  }

  for(size_t i = 0; i < varsCount; i++)
  {
    w2g = wg[i]*wg[i];
    // For c*x0' + d'*x_tilde
    out[i*dim+n_d3] = dnrm2sqrl(dim,g+i*dim)*x[i*dim+n_d3]/w2g; // c*x0'

    coef = -2.*g[i*dim]/w2g;
    for(size_t k = 1; k < dim; k++)
    {
      out[i*dim+n_d3] += coef*g[i*dim+k]*x[i*dim+n_d3+k]; // + d'*x_tilde
    }


    // For x0'*d + C*x_bar
    ddot = 0.;
    for(size_t l = 1; l < dim; l++)
    {
      ddot += g[i*dim+l]*x[i*dim+n_d3+l]; // Compute g_bar'*x_tilde
    }
    for(size_t k = 1; k < dim; k++)
    {
      out[i*dim+k+n_d3] = (x[i*dim+n_d3+k]+2.*(ddot-g[i*dim]*x[i*dim+n_d3])*g[i*dim+k])/w2g;
    }
  }
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
static void printDataProbMatlabFile(NumericsMatrix * M, double * f, NumericsMatrix * H, double * w, int d, int n, int m, double * mu, double * mu_r, FILE * file)
{
  fprintf(file,"d = %3i;\n",d);
  fprintf(file,"n = %6i;\n",n);
  fprintf(file,"m = %6i;\n",m);

  fprintf(file,"M = [\n");
  CSparseMatrix_print_in_Matlab_file(NM_triplet(M), 0, file);
  fprintf(file,"];\n");
  // fprintf(file,"M = sparse(int32(M(:,1)), int32(M(:,2)), M(:,3));\n");
  fprintf(file,"M = sparse(M(:,1), M(:,2), M(:,3));\n");

  fprintf(file,"H = [");
  CSparseMatrix_print_in_Matlab_file(NM_triplet(H), 0, file);
  fprintf(file,"];\n");
  // fprintf(file,"H = sparse(int32(H(:,1)), int32(H(:,2)), H(:,3));\n");
  fprintf(file,"H = sparse(H(:,1), H(:,2), H(:,3));\n");

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

  fprintf(file,"mu_r = [");
  for(int i = 0; i < n; i++)
  {
    fprintf(file,"%22.14e; ",mu_r[i]);
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
static void printInteresProbMatlabFile(int iteration, double * v, double * u_1, double * u_2, double * r_1, double * r_2, int d, int n, int m, FILE * file)
{
  fprintf(file,"v(%3i,:) = [",iteration+1);
  for(int i = 0; i < m; i++)
  {
    fprintf(file, "%20.16e, ", v[i]);
  }
  fprintf(file,"];\n");

  fprintf(file,"u1(%3i,:) = [",iteration+1);
  for(int i = 0; i < n*(d-2); i++)
  {
    fprintf(file, "%20.16e, ", u_1[i]);
  }
  fprintf(file,"];\n");

  fprintf(file,"u2(%3i,:) = [",iteration+1);
  for(int i = 0; i < n*(d-2); i++)
  {
    fprintf(file, "%20.16e, ", u_2[i]);
  }
  fprintf(file,"];\n");

  fprintf(file,"r1(%3i,:) = [",iteration+1);
  for(int i = 0; i < n*(d-2); i++)
  {
    fprintf(file, "%20.16e, ", r_1[i]);
  }
  fprintf(file,"];\n");

  fprintf(file,"r2(%3i,:) = [",iteration+1);
  for(int i = 0; i < n*(d-2); i++)
  {
    fprintf(file, "%20.16e, ", r_2[i]);
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





















/* --------------------------- Interior-point method implementation ------------------------------ */
/* General rolling friction contact problem */
/* Convex case: */
/* problem: min .5 v'*M*v + f'*v, s.t. H*v + w \in F (rolling friction cone)
   d = dimenion
   n = number of contact points
   m = number of degrees of freedom
   M = m x m matrix
   f = m-vector
   H = n*d x m matrix
   w = n*d-vector */
void grfc3d_IPM(GlobalRollingFrictionContactProblem* restrict problem, double* restrict reaction,
               double* restrict velocity, double* restrict globalVelocity,
               int* restrict info, SolverOptions* restrict options)
{
  clock_t t1 = clock();
  printf("\n\n#################### grfc3d_ipm.c 001 OK ####################\n\n");

 // verbose = 3;

  // the size of the problem detection
  unsigned int m = problem->M->size0;
  unsigned int nd = problem->H->size1;
  unsigned int d = problem->dimension;
  unsigned int n = problem->numberOfContacts;
  unsigned int m_plus_nd = m+nd;
  unsigned int d_minus_2 = d-2;
  unsigned int d_plus_1 = d+1;
  unsigned int n_dminus2 = n*d_minus_2;
  unsigned int n_dplus1 = n*d_plus_1;




  size_t posX = 0; // Used as the index in the source array
  size_t posY = 0; // Used as the index in the destination array

  NumericsMatrix* M = NULL, *minus_M = NULL;
  NumericsMatrix* H_origin = NULL, *minus_H = NULL;

  // globalRollingFrictionContact_display(problem);

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
    H_origin = NM_create(NM_SPARSE,  problem->H->size1,  problem->H->size0);
    NM_copy_to_sparse(NM_transpose(problem->H), H_origin, DBL_EPSILON);
  }
  else
  {
    H_origin = NM_transpose(problem->H);  // H <== H' because some different storages in fclib and paper
  }

  // initialize solver if it is not set
  int internal_allocation=0;
  if(!options->dWork || (options->dWorkSize != (size_t)(m + nd + n_dplus1)))
  {
    grfc3d_IPM_init(problem, options);
    internal_allocation = 1;
  }

  Grfc3d_IPM_data * data = (Grfc3d_IPM_data *)options->solverData;
  NumericsMatrix *P_mu = data->P_mu->mat;
  NumericsMatrix *P_mu_inv = data->P_mu->inv_mat;

  double *w_origin = problem->b;
  double *w = data->tmp_vault_nd[0];

  // compute -f
  double *f = problem->q;
  // cblas_dscal(m, -1.0, f, 1); // f <== -f because some different storages in fclib and paper
  // double *minus_f = (double*)calloc(m, sizeof(double));
  // cblas_dcopy(m, f, 1, minus_f, 1);
  // cblas_dscal(m, -1.0, minus_f, 1);

  double *iden;

  // change of variable
  // H_origin --> H
  NumericsMatrix *H = NM_multiply(P_mu, H_origin);

  // compute -H
  minus_H = NM_create(H->storageType, H->size0, H->size1);
  NM_copy(H, minus_H);
  NM_gemm(-1.0, H, NM_eye(H->size1), 0.0, minus_H);

  // w_origin --> w
  NM_gemv(1.0, P_mu, w_origin, 0.0, w);


  double alpha_primal_1 = 0;
  double alpha_primal_2 = 0;
  double alpha_dual_1 = 0;
  double alpha_dual_2 = 0;
  double alpha_primal = data->internal_params->alpha_primal;
  double alpha_dual = data->internal_params->alpha_dual;
  double barr_param = data->internal_params->barr_param;
  double sigma = data->internal_params->sigma;



  double alpha_primal_1_CHECK = 0;
  double alpha_primal_2_CHECK = 0;
  double alpha_dual_1_CHECK = 0;
  double alpha_dual_2_CHECK = 0;
  double alpha_primal_CHECK = 0;
  double alpha_dual_CHECK = 0;
  double barr_param_CHECK = 0;
  double sigma_CHECK = 0;



  cblas_dcopy(nd, data->starting_point->reaction, 1, reaction, 1);
  cblas_dcopy(nd, data->starting_point->velocity, 1, velocity, 1);
  cblas_dcopy(m, data->starting_point->globalVelocity, 1, globalVelocity, 1);


  /* 1. t, t_prime */
  double * t = (double*)calloc(n, sizeof(double));
  double * t_prime = (double*)calloc(n, sizeof(double));
  for(unsigned int i = 0; i < n; ++ i)
  {
    t[i] = 2.0;
    t_prime[i] = 1.0;
  }

  double * velocity_1 = (double*)calloc(n_dminus2, sizeof(double));   // = (t, u_bar)
  double * velocity_2 = (double*)calloc(n_dminus2, sizeof(double));   // = (t', u_tilde)
  double * reaction_1 = (double*)calloc(n_dminus2, sizeof(double));   // = (r0, r_bar)
  double * reaction_2 = (double*)calloc(n_dminus2, sizeof(double));   // = (r0, r_tilde)










  double * globalVelocity_CHECK = (double*)calloc(m, sizeof(double));   // = (t', u_tilde)
  double * velocity_CHECK = (double*)calloc(nd, sizeof(double));   // = (t, u_bar)
  double * reaction_CHECK = (double*)calloc(nd, sizeof(double));   // = (r0, r_bar)

  double * velocity_1_CHECK = (double*)calloc(n_dminus2, sizeof(double));   // = (t, u_bar)
  double * velocity_2_CHECK = (double*)calloc(n_dminus2, sizeof(double));   // = (t', u_tilde)
  double * reaction_1_CHECK = (double*)calloc(n_dminus2, sizeof(double));   // = (r0, r_bar)
  double * reaction_2_CHECK = (double*)calloc(n_dminus2, sizeof(double));   // = (r0, r_tilde)

  double * t_CHECK = (double*)calloc(n, sizeof(double));
  double * t_prime_CHECK = (double*)calloc(n, sizeof(double));

  cblas_dcopy(nd, data->starting_point->reaction, 1, reaction_CHECK, 1);
  cblas_dcopy(nd, data->starting_point->velocity, 1, velocity_CHECK, 1);
  cblas_dcopy(m, data->starting_point->globalVelocity, 1, globalVelocity_CHECK, 1);
  cblas_dcopy(n, t, 1, t_CHECK, 1);
  cblas_dcopy(n, t_prime, 1, t_prime_CHECK, 1);






  double tol = options->dparam[SICONOS_DPARAM_TOL];
  unsigned int max_iter = options->iparam[SICONOS_IPARAM_MAX_ITER];
  double sgmp1 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1];
  double sgmp2 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2];
  double sgmp3 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3];
  double gmmp0 = 0.999;
  double gmmp1 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1];
  double gmmp2 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2];

  int hasNotConverged = 1;
  unsigned int iteration = 0;
  double pinfeas = 1e300;
  double dinfeas = 1e300;
  double complem_1 = 1e300;
  double complem_2 = 1e300;
  double gapVal = 1e300;
  double relgap = 1e300;
  double u1dotr1 = 1e300;     // u1 = velecity_1, r1 = reaction_1
  double u2dotr2 = 1e300;     // u2 = velecity_2, r2 = reaction_2
  double error_array[6];

  double pinfeas_CHECK = 1e300;
  double dinfeas_CHECK = 1e300;
  double complem_1_CHECK = 1e300;
  double complem_2_CHECK = 1e300;
  double gapVal_CHECK = 1e300;
  double relgap_CHECK = 1e300;
  double u1dotr1_CHECK = 1e300;     // u1 = velecity_1, r1 = reaction_1
  double u2dotr2_CHECK = 1e300;     // u2 = velecity_2, r2 = reaction_2


  double gmm_CHECK = gmmp0;
  double barr_param_a_CHECK, e_CHECK;


  double gmm = gmmp0;
  double barr_param_a, e;
  double norm_f = cblas_dnrm2(m, f, 1);
  double norm_w = cblas_dnrm2(nd, w, 1);




  double *primalConstraint = data->tmp_vault_nd[1];
  double *dualConstraint = data->tmp_vault_m[0];
  double *complemConstraint_1 = data->tmp_vault_n_dminus2[0];
  double *complemConstraint_2 = data->tmp_vault_n_dminus2[1];

  double *primalConstraint_CHECK = (double*)calloc(nd,sizeof(double));
  double *dualConstraint_CHECK = (double*)calloc(m,sizeof(double));
  double *complemConstraint_1_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double *complemConstraint_2_CHECK = (double*)calloc(n_dminus2,sizeof(double));

  double *d_globalVelocity = (double*)calloc(m,sizeof(double));
  double *d_velocity = (double*)calloc(nd,sizeof(double));
  double *d_velocity_1 = (double*)calloc(n_dminus2,sizeof(double));
  double *d_velocity_2 = (double*)calloc(n_dminus2,sizeof(double));
  double *d_reaction = (double*)calloc(nd,sizeof(double));
  double *d_reaction_1 = (double*)calloc(n_dminus2,sizeof(double));
  double *d_reaction_2 = (double*)calloc(n_dminus2,sizeof(double));
  double *d_t = (double*)calloc(n, sizeof(double));
  double *d_t_prime = (double*)calloc(n, sizeof(double));


  double *d_globalVelocity_CHECK = (double*)calloc(m,sizeof(double));
  double *d_velocity_CHECK = (double*)calloc(nd,sizeof(double));
  double *d_velocity_1_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double *d_velocity_2_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double *d_reaction_CHECK = (double*)calloc(nd,sizeof(double));
  double *d_reaction_1_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double *d_reaction_2_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double *d_t_CHECK = (double*)calloc(n, sizeof(double));
  double *d_t_prime_CHECK = (double*)calloc(n, sizeof(double));

  double *diff_d_globalVelocity = (double*)calloc(m, sizeof(double));
  double *diff_d_velocity = (double*)calloc(nd, sizeof(double));
  double *diff_d_reaction = (double*)calloc(nd, sizeof(double));



  double *rhs = options->dWork;
  double *rhs_tmp = NULL;
  double *rhs_CHECK = (double*)calloc(m + nd + n*(d+1), sizeof(double));
  double *rhs_CHECK_save = (double*)calloc(m + nd + n*(d+1), sizeof(double));
  double *rhs_NT_nonRe = (double*)calloc(m + nd + n*(d+1), sizeof(double));
  double *rhs_NT_nonRe_save = (double*)calloc(m + nd + n*(d+1), sizeof(double));
  double *rhs_NT_Re_save = (double*)calloc(m + nd, sizeof(double));


  // TO DO: Allocate temporarily here, need to move to another place after
  double *r1r2 =(double*)calloc(n_dplus1,sizeof(double));
  double *u1u2 = (double*)calloc(n_dplus1,sizeof(double));
  double *velocity_1_inv = (double*)calloc(n_dminus2,sizeof(double));
  double *velocity_2_inv = (double*)calloc(n_dminus2,sizeof(double));
  double *Qpvelocity_1 = (double*)calloc(n_dminus2,sizeof(double));
  double *Qpvelocity_2 = (double*)calloc(n_dminus2,sizeof(double));
  double *Qpvelocity_1_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double *Qpvelocity_2_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double *QpinvReaction_1 = (double*)calloc(n_dminus2,sizeof(double));
  double *QpinvReaction_2 = (double*)calloc(n_dminus2,sizeof(double));
  double *Qp_uQpinv_r_1 = (double*)calloc(n_dminus2,sizeof(double));
  double *Qp_uQpinv_r_2 = (double*)calloc(n_dminus2,sizeof(double));
  double *Qpvelocity_inv_dvdr_1 = (double*)calloc(n_dminus2,sizeof(double));
  double *Qpvelocity_inv_dvdr_2 = (double*)calloc(n_dminus2,sizeof(double));
  double *Qpvelocity_inv_dvdr_1_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double *Qpvelocity_inv_dvdr_2_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double *QpQpvelocity_inv_dvdr_1_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double *QpQpvelocity_inv_dvdr_2_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double *u_inv_dvdr_1 = (double*)calloc(n_dminus2,sizeof(double));
  double *u_inv_dvdr_2 = (double*)calloc(n_dminus2,sizeof(double));
  double *QpQpvelocity_inv_dvdr_1 = (double*)calloc(n_dminus2,sizeof(double));
  double *QpQpvelocity_inv_dvdr_2 = (double*)calloc(n_dminus2,sizeof(double));
  double *Qpvelocity_1_inv = (double*)calloc(n_dminus2,sizeof(double));
  double *Qpvelocity_2_inv = (double*)calloc(n_dminus2,sizeof(double));
  double *Qpvelocity_1_inv_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double *Qpvelocity_2_inv_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double *dvdr_jprod = (double*)calloc(n_dplus1,sizeof(double));
  double *tmp_n_dplus1 = (double*)calloc(n_dplus1,sizeof(double));
  double *Qinv2_x = (double*)calloc(n_dplus1,sizeof(double));
  double *Hdv = (double*)calloc(nd,sizeof(double));
  double *tmp_nd = (double*)calloc(nd,sizeof(double));
  double *tmp_nd_2 = (double*)calloc(nd,sizeof(double));
  // double *d_reaction_reduced = (double*)calloc(nd,sizeof(double));

  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING] && options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM])
  {
    minus_M = NM_create(M->storageType, M->size0, M->size1);    // store the matrix -M to build the matrix of the Newton linear system
    /* Create the matrix -M to build the matrix of the reduced linear system */
    NM_copy(M, minus_M); /* useful ? */
    NM_gemm(-1.0, M, NM_eye(M->size1), 0.0, minus_M);
  }
  double* f_NT = NULL;
  double* g_NT = NULL;
  float_type* wf_NT = NULL;
  float_type* wg_NT = NULL;



  double *gv_plus_dgv = data->tmp_vault_m[1];
  //double *vr_jprod = data->tmp_vault_nd[2];
  double *v_plus_dv = data->tmp_vault_nd[3];
  double *r_plus_dr = data->tmp_vault_nd[4];
  double *dvdr_jprod_1 = data->tmp_vault_n_dminus2[2];
  double *dvdr_jprod_2 = data->tmp_vault_n_dminus2[3];
  double *vr_jprod_1 = data->tmp_vault_n_dminus2[22];
  double *vr_jprod_2 = data->tmp_vault_n_dminus2[23];



  double *gv_plus_dgv_CHECK = (double*)calloc(m,sizeof(double));
  double *v_plus_dv_CHECK = (double*)calloc(nd,sizeof(double));
  double *r_plus_dr_CHECK = (double*)calloc(nd,sizeof(double));
  double *dvdr_jprod_1_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double *dvdr_jprod_2_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double *vr_jprod_1_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double *vr_jprod_2_CHECK = (double*)calloc(n_dminus2,sizeof(double));







  NumericsMatrix *Jac; /* Jacobian matrix */
  NumericsMatrix *Jac_CHECK, *Jac_NT_nonRe; /* Jacobian matrix */
  long Jac_nzmax,  Jac_CHECK_nzmax, Jac_NT_nonRe_nzmax;
  size_t M_nzmax = NM_nnz(M);
  size_t H_nzmax = NM_nnz(H);

  NumericsMatrix *J = compute_J_matrix(n); /* use for Jac in the NT scaling case */

  double full_error = 1e300;
  char fws = ' '; /* finish without scaling */

  // if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING] > 0)
  // {
  //   J = compute_J_matrix(n);
  // }



  /* -------------------------- Display problem info -------------------------- */
  if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO] == SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_YES)
  {
    numerics_printf_verbose(1,"---- GRFC3D - IPM - Problem information");
    numerics_printf_verbose(1,"---- GRFC3D - IPM - 1-norm of M = %g norm of f = %g ", NM_norm_1(M), norm_f);
    numerics_printf_verbose(1,"---- GRFC3D - IPM - inf-norm of M = %g ", NM_norm_inf(M));

    numerics_printf_verbose(1,"---- GRFC3D - IPM - 1-norm of H = %g norm of w = %g ", NM_norm_1(problem->H), norm_w);
    numerics_printf_verbose(1,"---- GRFC3D - IPM - inf-norm of H = %g ", NM_norm_inf(problem->H));
    numerics_printf_verbose(1,"---- GRFC3D - IPM - M is symmetric = %i ", NM_is_symmetric(M));

    numerics_printf_verbose(1,"---- GRFC3D - IPM - M size = (%i, %i) ", M->size0, M->size1);
    numerics_printf_verbose(1,"---- GRFC3D - IPM - H size = (%i, %i) ", problem->H->size0, problem->H->size1);
  }


  /* -------------------------- Display IPM iterations -------------------------- */
  numerics_printf_verbose(-1, "problem dimensions d, n, m: %1i, %6i, %6i\n",d, n, m);
  numerics_printf_verbose(-1, "| it  |  rel gap  | pinfeas  | dinfeas  | <u1, r1> | <u2, r2> | complem1 | complem2 | full err | barparam | alpha_p  | alpha_d  |  sigma   | |dv|/|v| | |du|/|u| | |dr|/|r| |");
  numerics_printf_verbose(-1, "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------");

  double * p_bar = data->tmp_vault_n_dminus2[4];
  double * p_tilde = data->tmp_vault_n_dminus2[5];
  double * p2_bar = data->tmp_vault_n_dminus2[20];
  double * p2_tilde = data->tmp_vault_n_dminus2[21];
  NumericsMatrix* Qp_bar = NULL;
  NumericsMatrix* Qp_tilde = NULL;
  NumericsMatrix* Qpinv_bar = NULL;
  NumericsMatrix* Qpinv_tilde = NULL;
  NumericsMatrix* Qpinv2_bar = NULL;
  NumericsMatrix* Qpinv2_tilde = NULL;
  NumericsMatrix* Qinv = NULL;
  NumericsMatrix* Qinv2 = NULL;
  NumericsMatrix* JQinv = NULL;
  NumericsMatrix* JQinv2 = NULL;
  NumericsMatrix* Qp2_bar = NULL;
  NumericsMatrix* Qp2_tilde = NULL;
  NumericsMatrix* P_inv = NULL;
  NumericsMatrix* PinvH = NULL;



  NumericsMatrix* Qp_bar_CHECK = NULL;
  NumericsMatrix* Qp_tilde_CHECK = NULL;
  NumericsMatrix* Qpinv_bar_CHECK = NULL;
  NumericsMatrix* Qpinv_tilde_CHECK = NULL;
  NumericsMatrix* Qp2_bar_CHECK = NULL;
  NumericsMatrix* Qp2_tilde_CHECK = NULL;





  double * velocity_1t = data->tmp_vault_n_dminus2[6];
  double * velocity_2t = data->tmp_vault_n_dminus2[7];
  double * d_velocity_1t = data->tmp_vault_n_dminus2[8];
  double * d_velocity_2t = data->tmp_vault_n_dminus2[9];
  double * d_reaction_1t = data->tmp_vault_n_dminus2[10];
  double * d_reaction_2t = data->tmp_vault_n_dminus2[11];
  double * velocity_1t_inv = data->tmp_vault_n_dminus2[12];
  double * velocity_2t_inv = data->tmp_vault_n_dminus2[13];
  double * Qp_velocity_1t_inv = data->tmp_vault_n_dminus2[14];
  double * Qp_velocity_2t_inv = data->tmp_vault_n_dminus2[15];
  double * F_velocity_1t_inv = data->tmp_vault_n_dminus2[16];
  double * F_velocity_2t_inv = data->tmp_vault_n_dminus2[17];
  double * tmp1 = data->tmp_vault_n_dminus2[18];
  double * tmp2 = data->tmp_vault_n_dminus2[19];




  double * velocity_1t_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double * velocity_2t_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double * d_velocity_1t_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double * d_velocity_2t_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double * d_reaction_1t_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double * d_reaction_2t_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double * velocity_1t_inv_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double * velocity_2t_inv_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double * Qp_velocity_1t_inv_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double * Qp_velocity_2t_inv_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double * F_velocity_1t_inv_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double * F_velocity_2t_inv_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double * tmp1_CHECK = (double*)calloc(n_dminus2,sizeof(double));
  double * tmp2_CHECK = (double*)calloc(n_dminus2,sizeof(double));









  FILE * iterates;
  FILE * matrixH;
  char matlab_name[200];
  sprintf(matlab_name, "iteratesNC%d.m",n);


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
    iterates = fopen(matlab_name, "w");
    printDataProbMatlabFile(M, f, H, w, d, n, m, problem->mu, problem->mu_r, iterates);
  }




  ComputeErrorGlobalRollingPtr computeError = NULL;
  computeError = (ComputeErrorGlobalRollingPtr)&grfc3d_compute_error;




  /* -------------------------- Check the full criterion -------------------------- */
  while(iteration < max_iter)
  {

    /* -------------------------- Extract vectors -------------------------- */
    /* 2. velocity_1 = (t, u_bar), velocity_2 = (t_prime, u_tilde) */
    extract_vector(velocity, nd, n, 2, 3, velocity_1);
    extract_vector(velocity, nd, n, 4, 5, velocity_2);
    for(unsigned int i = 0; i < n; i++)
    {
      velocity_1[i*d_minus_2] = t[i];
      velocity_2[i*d_minus_2] = t_prime[i];
    }

    /* 3. reaction_1 = (r0, r_bar), reaction_2 = (r0, r_tilde) */
    extract_vector(reaction, nd, n, 2, 3, reaction_1);
    extract_vector(reaction, nd, n, 4, 5, reaction_2);



    // // For computing NT Non-reduced at the same time
    // extract_vector(velocity_CHECK, nd, n, 2, 3, velocity_1_CHECK);
    // extract_vector(velocity_CHECK, nd, n, 4, 5, velocity_2_CHECK);
    // for(unsigned int i = 0; i < n; i++)
    // {
    //   velocity_1_CHECK[i*d_minus_2] = t_CHECK[i];
    //   velocity_2_CHECK[i*d_minus_2] = t_prime_CHECK[i];
    // }

    // /* 3. reaction_1 = (r0, r_bar), reaction_2 = (r0, r_tilde) */
    // extract_vector(reaction_CHECK, nd, n, 2, 3, reaction_1_CHECK);
    // extract_vector(reaction_CHECK, nd, n, 4, 5, reaction_2_CHECK);









    /* writing data in a Matlab file */
    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
    {
      printInteresProbMatlabFile(iteration, globalVelocity, velocity_1, velocity_2, reaction_1, reaction_2, d, n, m, iterates);
    }





    /* -------------------------- ?????? -------------------------- */
    if ((options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING] == 1) && (full_error <= 1e-6) && (fws==' '))
    {
      // To solve the problem very accurately, the algorithm switches to a direct solution of the linear system without scaling//
      options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING] = 0;
      fws = '*';
    }



    /* Primal residual = velocity - H * globalVelocity - w */
    primalResidual(velocity, H, globalVelocity, w, primalConstraint, &pinfeas);
    /* Dual residual = M*globalVelocity - H'*reaction + f */
    dualResidual(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas);



    // // For computing NT Non-reduced at the same time
    // primalResidual(velocity_CHECK, H, globalVelocity_CHECK, w, primalConstraint_CHECK, &pinfeas_CHECK);\
    // dualResidual(M, globalVelocity_CHECK, H, reaction_CHECK, f, dualConstraint_CHECK, &dinfeas_CHECK);





    // if (iteration <= 1)
    // {
    //   printf("\n\n========== PRINTING FOR DEBUG 1 ==========\n");
    //   printf("iteration = %i\n", iteration);
    //   // printf("Vector v:\n");
    //   // NM_vector_display(globalVelocity,m);
    //   // printf("\n\nVector u:\n");
    //   // NM_vector_display(velocity,nd);
    //   // printf("\n\nVector r:\n");
    //   // NM_vector_display(reaction,nd);

    //   printf("primalConstraint = %6.20e\n\n", pinfeas);
    //   printf("dualConstraint = %6.20e\n\n", dinfeas);
    //   printf("========== END PRINTING FOR DEBUG 1 ==========\n\n");
    //   // break;
    // }




    // /* Dual gap = (primal value - dual value)/ (1 + abs(primal value) + abs(dual value)) */
    // // Note: primal objectif func = 1/2 * v' * M *v + f' * v
    // dualgap = dualGap(M, f, w, globalVelocity, reaction, nd, m);


    /* Gap value = u'.v */
    gapVal = cblas_ddot(nd, reaction, 1, velocity, 1);

    // Note: primal objectif func = 1/2 * v' * M *v + f' * v
    relgap = relGap(M, f, w, globalVelocity, reaction, nd, m, gapVal);

    // barr_param = (gapVal / nd)*sigma;
    // barr_param = gapVal / nd;
    barr_param = gapVal / (n);
    //barr_param = complemResidualNorm(velocity, reaction, nd, n);
    //barr_param = (fws=='*' ? complemResidualNorm(velocity, reaction, nd, n)/n : complemResidualNorm_p(velocity, reaction, nd, n)/n) ;
    //barr_param = fabs(barr_param);


    complem_1 = complemResidualNorm(velocity_1, reaction_1, n_dminus2, n); // (t, u_bar) o (r0, r_bar)
    complem_2 = complemResidualNorm(velocity_2, reaction_2, n_dminus2, n); // (t', u_tilde) o (r0, r_tilde)





    // // For computing NT Non-reduced at the same time
    // gapVal_CHECK = cblas_ddot(nd, reaction_CHECK, 1, velocity_CHECK, 1);
    // relgap_CHECK = relGap(M, f, w, globalVelocity_CHECK, reaction_CHECK, nd, m, gapVal_CHECK);
    // barr_param_CHECK = gapVal_CHECK / (n);
    // complem_1_CHECK = complemResidualNorm(velocity_1_CHECK, reaction_1_CHECK, n_dminus2, n); // (t, u_bar) o (r0, r_bar)
    // complem_2_CHECK = complemResidualNorm(velocity_2_CHECK, reaction_2_CHECK, n_dminus2, n); // (t', u_tilde) o (r0, r_tilde)






    /* ----- return to original variables ------ */
    NM_gemv(1.0, P_mu_inv, velocity, 0.0, data->original_point->velocity);
    NM_gemv(1.0, P_mu, reaction, 0.0, data->original_point->reaction);








    if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] == 1)         // non-smooth case
    {
      (*computeError)(problem,
                    data->original_point->reaction, data->original_point->velocity, globalVelocity,
                    tol, &full_error, 1);
    }
    else if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] == 0)    // convex case
    {
      (*computeError)(problem,
                    data->original_point->reaction, data->original_point->velocity, globalVelocity,
                    tol, &full_error, 0);
    }
    // error = fmax(barr_param, fmax(complem_1, fmax(complem_2, fmax(pinfeas, dinfeas))));





    u1dotr1 = cblas_ddot(n_dminus2, velocity_1, 1, reaction_1, 1);
    u2dotr2 = cblas_ddot(n_dminus2, velocity_2, 1, reaction_2, 1);



    // // For computing NT Non-reduced at the same time
    // u1dotr1_CHECK = cblas_ddot(n_dminus2, velocity_1_CHECK, 1, reaction_1_CHECK, 1);
    // u2dotr2_CHECK = cblas_ddot(n_dminus2, velocity_2_CHECK, 1, reaction_2_CHECK, 1);





    error_array[0] = pinfeas;
    error_array[1] = dinfeas;
    error_array[2] = u1dotr1;
    error_array[3] = u2dotr2;
    error_array[4] = complem_1;
    error_array[5] = complem_2;







    // check exit condition
    //if (error <= tol) //((NV_max(error, 4) <= tol) || (err <= tol))
    if (NV_max(error_array, 4) <= tol)
    {
      numerics_printf_verbose(-1, "| %3i%c| %9.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e |",
                              iteration, fws, relgap, pinfeas, dinfeas, u1dotr1, u2dotr2, complem_1, complem_2, full_error, barr_param);

      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
        printInteresProbMatlabFile(iteration, globalVelocity, velocity_1, velocity_2, reaction_1, reaction_2, d, n, m, iterates);

      hasNotConverged = 0;
      break;
    }















    /* -------------------------- Compute Jacobian -------------------------- */
    /* with NT scaling: /scaling/ */
    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING])
    {
      /* with reduced system: /scaling/reduced/ */
      if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM])
      {
        if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] == SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F)
        {
          /*   1. Build the reduced Jacobian matrix
           *
           *            m        nd
           *        |  -M     H'P^-1' | m
           *  Jac = |                 |
           *        | P^-1H      I    | nd
           */
          Qp_bar = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
          Qpinv_bar = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
          Qpinv2_bar = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
          Qp_tilde = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
          Qpinv_tilde = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
          Qpinv2_tilde = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
          NM_triplet_alloc(Qp_bar, d_minus_2 * d_minus_2 * n);
          NM_triplet_alloc(Qpinv_bar, d_minus_2 * d_minus_2 * n);
          NM_triplet_alloc(Qpinv2_bar, d_minus_2 * d_minus_2 * n);
          NM_triplet_alloc(Qp_tilde, d_minus_2 * d_minus_2 * n);
          NM_triplet_alloc(Qpinv_tilde, d_minus_2 * d_minus_2 * n);
          NM_triplet_alloc(Qpinv2_tilde, d_minus_2 * d_minus_2 * n);

          f_NT = (double*)calloc(n_dminus2, sizeof(double));
          g_NT = (double*)calloc(n_dminus2, sizeof(double));
          wf_NT = (float_type*)calloc(n, sizeof(float_type));
          wg_NT = (float_type*)calloc(n, sizeof(float_type));

          family_of_F(velocity_1, reaction_1, n_dminus2, n, f_NT, wf_NT, Qp_bar, Qpinv_bar, NULL, Qpinv2_bar);
          family_of_F(velocity_2, reaction_2, n_dminus2, n, g_NT, wg_NT, Qp_tilde, Qpinv_tilde, NULL, Qpinv2_tilde);

          // Qp2_bar = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
          // Qp2_tilde = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
          // NM_triplet_alloc(Qp2_bar, d_minus_2 * d_minus_2 * n);
          // NM_triplet_alloc(Qp2_tilde, d_minus_2 * d_minus_2 * n);
          // family_of_F(velocity_1, reaction_1, n_dminus2, n, f_NT, wf_NT, Qp_bar, Qpinv_bar, Qp2_bar, Qpinv2_bar);
          // family_of_F(velocity_2, reaction_2, n_dminus2, n, g_NT, wg_NT, Qp_tilde, Qpinv_tilde, Qp2_tilde, Qpinv2_tilde);


          P_inv = Pinv(f_NT, g_NT, wf_NT, wg_NT, n_dminus2, n);
          PinvH = NM_multiply(P_inv,H);
          // PinvH = multiply_PinvH(f_NT, g_NT, wf_NT, wg_NT, n_dminus2, n, H);


          Jac = NM_create(NM_SPARSE, m + nd, m + nd);
          Jac_nzmax = M_nzmax + 2*H_nzmax + nd;
          NM_triplet_alloc(Jac, Jac_nzmax);
          Jac->matrix2->origin = NSM_TRIPLET;
          NM_insert(Jac, minus_M, 0, 0);
          NM_insert(Jac, NM_transpose(PinvH), 0, m);
          NM_insert(Jac, PinvH, m, 0);
          NM_insert(Jac, NM_eye(nd), m, m);
          // if(iteration==0)
          // {
          //   NM_display(Jac);
          //   break;
          // }
          //free(f_NT); free(g_NT); free(wf_NT); free(wg_NT);






          // //For computing NT Non-reduced at the same time
          // Qp_bar_CHECK = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
          // Qpinv_bar_CHECK = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
          // Qp2_bar_CHECK = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
          // Qp_tilde_CHECK = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
          // Qpinv_tilde_CHECK = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
          // Qp2_tilde_CHECK = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
          // NM_triplet_alloc(Qp_bar_CHECK, d_minus_2 * d_minus_2 * n);
          // NM_triplet_alloc(Qpinv_bar_CHECK, d_minus_2 * d_minus_2 * n);
          // NM_triplet_alloc(Qp2_bar_CHECK, d_minus_2 * d_minus_2 * n);
          // NM_triplet_alloc(Qp_tilde_CHECK, d_minus_2 * d_minus_2 * n);
          // NM_triplet_alloc(Qpinv_tilde_CHECK, d_minus_2 * d_minus_2 * n);
          // NM_triplet_alloc(Qp2_tilde_CHECK, d_minus_2 * d_minus_2 * n);

          // family_of_F(velocity_1_CHECK, reaction_1_CHECK, n_dminus2, n, NULL, NULL, Qp_bar_CHECK, Qpinv_bar_CHECK, Qp2_bar_CHECK, NULL);
          // family_of_F(velocity_2_CHECK, reaction_2_CHECK, n_dminus2, n, NULL, NULL, Qp_tilde_CHECK, Qpinv_tilde_CHECK, Qp2_tilde_CHECK, NULL);

          // Jac_CHECK = NM_create(NM_SPARSE, m + nd + n_dplus1, m + nd + n_dplus1);
          // Jac_CHECK_nzmax = M_nzmax + 2*H_nzmax + 2*9*n*n + 6*n;
          // NM_triplet_alloc(Jac_CHECK, Jac_CHECK_nzmax);
          // Jac_CHECK->matrix2->origin = NSM_TRIPLET;
          // NM_insert(Jac_CHECK, M, 0, 0);
          // NM_insert(Jac_CHECK, NM_transpose(minus_H), 0, m);
          // NM_insert(Jac_CHECK, minus_H, m, 0);
          // NM_insert(Jac_CHECK, J, m, m_plus_nd);
          // NM_insert(Jac_CHECK, NM_transpose(J), m_plus_nd, m);
          // // NM_insert(Jac_CHECK, Qp2_bar_CHECK, m_plus_nd, m_plus_nd);
          // // NM_insert(Jac_CHECK, Qp2_tilde_CHECK, m_plus_nd+n_dminus2, m_plus_nd+n_dminus2);
          // NM_insert(Jac_CHECK, Qp2_bar, m_plus_nd, m_plus_nd);
          // NM_insert(Jac_CHECK, Qp2_tilde, m_plus_nd+n_dminus2, m_plus_nd+n_dminus2);











        }
        else if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] == SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP)
        {
          // TO DO: not yet
          printf("\n\nTO DO: not yet\n\n");
          break;
        }
      } // END OF with reduced system

      else /* with NOT-reduced system:  /scaling/NOT-reduced/ */
      {
        /* -------------------------- Compute NT directions -------------------------- */
        if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] == SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F)
        {
          // Qp_bar = NTmat(velocity_1, reaction_1, n_dminus2, n);
          // Qp2_bar = NTmatsqr(velocity_1, reaction_1, n_dminus2, n);

          // Qp_tilde = NTmat(velocity_2, reaction_2, n_dminus2, n);
          // Qp2_tilde = NTmatsqr(velocity_2, reaction_2, n_dminus2, n);

          Qp_bar = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
          Qpinv_bar = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
          Qp2_bar = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
          Qp_tilde = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
          Qpinv_tilde = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
          Qp2_tilde = NM_create(NM_SPARSE, n_dminus2, n_dminus2);
          NM_triplet_alloc(Qp_bar, d_minus_2 * d_minus_2 * n);
          NM_triplet_alloc(Qpinv_bar, d_minus_2 * d_minus_2 * n);
          NM_triplet_alloc(Qp2_bar, d_minus_2 * d_minus_2 * n);
          NM_triplet_alloc(Qp_tilde, d_minus_2 * d_minus_2 * n);
          NM_triplet_alloc(Qpinv_tilde, d_minus_2 * d_minus_2 * n);
          NM_triplet_alloc(Qp2_tilde, d_minus_2 * d_minus_2 * n);

          family_of_F(velocity_1, reaction_1, n_dminus2, n, NULL, NULL, Qp_bar, Qpinv_bar, Qp2_bar, NULL);
          family_of_F(velocity_2, reaction_2, n_dminus2, n, NULL, NULL, Qp_tilde, Qpinv_tilde, Qp2_tilde, NULL);


          // if (iteration == 0)
          // {
          //   double* f_NT = (double*)calloc(n_dminus2, sizeof(double));
          //   double* g_NT = (double*)calloc(n_dminus2, sizeof(double));
          //   float_type* wf_NT = (float_type*)calloc(n, sizeof(float_type));
          //   float_type* wg_NT = (float_type*)calloc(n, sizeof(float_type));

          //   family_of_F(velocity_1, reaction_1, n_dminus2, n, f_NT, wf_NT, NULL, NULL, NULL, NULL);
          //   family_of_F(velocity_2, reaction_2, n_dminus2, n, g_NT, wg_NT, NULL, NULL, NULL, NULL);

          //   printf("\n\nVector u_1:\n");
          //   NM_vector_display(velocity_1,n_dminus2);
          //   printf("\n\nVector u_2:\n");
          //   NM_vector_display(velocity_2,n_dminus2);
          //   printf("\n\nVector r_1:\n");
          //   NM_vector_display(reaction_1,n_dminus2);
          //   printf("\n\nVector r_2:\n");
          //   NM_vector_display(reaction_2,n_dminus2);

          //   printf("\n\nf_NT = \n");
          //   NV_display(f_NT, n_dminus2);
          //   printf("\n\ng_NT = \n");
          //   NV_display(g_NT, n_dminus2);

          //   printf("\n\nwf_NT = \n");
          //   for(int i=0; i<n; i++){printf("%Lf\n",wf_NT[i]);}
          //   printf("\n\nwg_NT = \n");
          //   for(int i=0; i<n; i++){printf("%Lf\n",wg_NT[i]);}

          //   NumericsMatrix * P_inv = Pinv(f_NT, g_NT, wf_NT, wg_NT, n_dminus2, n);
          //   printf("\n\nP_inv = \n");
          //   NM_display(P_inv);
          //   break;
          // }
        }
        else if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] == SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP)
        {
          Nesterov_Todd_vector(0, velocity_1, reaction_1, n_dminus2, n, p_bar);
          Qp_bar = QRmat(p_bar, n_dminus2, n);
          Nesterov_Todd_vector(2, velocity_1, reaction_1, n_dminus2, n, p2_bar);
          Qp2_bar = QRmat(p2_bar, n_dminus2, n);

          Nesterov_Todd_vector(0, velocity_2, reaction_2, n_dminus2, n, p_tilde);
          Qp_tilde = QRmat(p_tilde, n_dminus2, n);
          Nesterov_Todd_vector(2, velocity_2, reaction_2, n_dminus2, n, p2_tilde);
          Qp2_tilde = QRmat(p2_tilde, n_dminus2, n);
        }
        /* -------------------------- Compute Jacobian -------------------------- */
        /*   1. Build the Jacobian matrix
         *
         *           m    nd   n(d+1)
         *        |  M    -H'    0   | m
         *        |                  |
         *  Jac = | -H     0     J   | nd
         *        |                  |
         *        |  0     J'    Q^2 | n(d+1)
         *
         *          n(d-2)    n(d-2)
         *        | Qp_bar      0    | n(d-2)
         *   Q  = |                  |
         *        |  0      Qp_tilde | n(d-2)
         */

        Jac = NM_create(NM_SPARSE, m + nd + n_dplus1, m + nd + n_dplus1);
        Jac_nzmax = M_nzmax + 2*H_nzmax + 2*9*n*n + 6*n;
        NM_triplet_alloc(Jac, Jac_nzmax);
        Jac->matrix2->origin = NSM_TRIPLET;
        NM_insert(Jac, M, 0, 0);
        NM_insert(Jac, NM_transpose(minus_H), 0, m);
        NM_insert(Jac, minus_H, m, 0);
        NM_insert(Jac, J, m, m_plus_nd);
        NM_insert(Jac, NM_transpose(J), m_plus_nd, m);
        NM_insert(Jac, Qp2_bar, m_plus_nd, m_plus_nd);
        NM_insert(Jac, Qp2_tilde, m_plus_nd+n_dminus2, m_plus_nd+n_dminus2);

      } // END OF /scaling/NOT-reduced/
    }   // END OF /scaling/

    else /* without NT scaling: /NOT-scaling/ */
    {
      /* with reduced system: /NOT-scaling/reduced/ */
      if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM])
      {
        printf("\n\nNo work 1\n\n");
        break;
      } // END OF /NOT-scaling/reduced/

      /* with NOT-reduced system: /NOT-scaling/NOT-reduced/ */
      else
      {
        /*   1. Build the Jacobian matrix
         *
         *               m     nd      n(d-2)    n(d-2)
         *            |  M     -H'        0         0    | m
         *            |                                  |
         *    Jac   = | -H      0      [...... J ......] | nd
         *            |                                  |
         *            |  0   block_1   arw(r1)      0    | n(d-2)
         *            |                                  |
         *            |  0   block_2      0      arw(r2) | n(d-2)
         *
         *  where
         *
         *                1       2        2
         *            |   t     u_bar'     0  | 1
         *  block_1 = |                       |         x n blocks
         *            | u_bar     tI       0  | 2
         *
         *  and
         *
         *                1       2        2
         *            |   t'      0     u_tilde' | 1
         *  block_2 = |                          |      x n blocks
         *            | u_tilde   0       t'I    | 2
         *
         */
        Jac = NM_create(NM_SPARSE, m + nd + n_dplus1, m + nd + n_dplus1);
        Jac_nzmax = M_nzmax + 2*H_nzmax + 2*9*n*n + 6*n;
        NM_triplet_alloc(Jac, Jac_nzmax);
        Jac->matrix2->origin = NSM_TRIPLET;
        NM_insert(Jac, M, 0, 0);
        NM_insert(Jac, NM_transpose(minus_H), 0, m);
        NM_insert(Jac, minus_H, m, 0);
        NM_insert(Jac, J, m, m_plus_nd);

        /* Create matrices block_1, block_2 */
        int n3, n5;
        n3 = 3*n;
        n5 = 5*n;
        NumericsMatrix * block_1 = NM_create(NM_SPARSE, n5, n3);
        NumericsMatrix * block_2 = NM_create(NM_SPARSE, n5, n3);
        long blocks_nzmax = 3*2*n;
        NM_triplet_alloc(block_1, blocks_nzmax);
        NM_triplet_alloc(block_2, blocks_nzmax);
        NM_fill(block_1, NM_SPARSE, n5, n3, block_1->matrix2);
        NM_fill(block_2, NM_SPARSE, n5, n3, block_2->matrix2);
        block_1->matrix2->origin = NSM_TRIPLET;
        block_2->matrix2->origin = NSM_TRIPLET;

        for(size_t i = 0; i < n; ++i)
        {
          posX = i * d_minus_2;   // row = 3*i
          posY = i * d;           // col = 5*i
          NM_entry(block_1, posX, posY, velocity_1[posX]);
          NM_entry(block_2, posX, posY, velocity_2[posX]);

          for(size_t j = 1; j < d_minus_2; ++j)
          {
            NM_entry(block_1, posX, posY + j, velocity_1[posX + j]);
            NM_entry(block_1, posX + j, posY, velocity_1[posX + j]);
            NM_entry(block_1, posX + j, posY + j, velocity_1[posX]);

            NM_entry(block_2, posX, posY + j + 2, velocity_2[posX + j]);
            NM_entry(block_2, posX + j, posY, velocity_2[posX + j]);
            NM_entry(block_2, posX + j, posY + j + 2, velocity_2[posX]);
          }
        }

        /* Continue adding into Jac matrix */
        NM_insert(Jac, block_1, m_plus_nd, m);
        NM_insert(Jac, Arrow_repr(reaction_1, n_dminus2, n), m_plus_nd, m_plus_nd);
        NM_insert(Jac, block_2, m_plus_nd+n_dminus2, m);
        NM_insert(Jac, Arrow_repr(reaction_2, n_dminus2, n), m_plus_nd+n_dminus2, m_plus_nd+n_dminus2);

        NM_clear(block_1);
        free(block_1);
        NM_clear(block_2);
        free(block_2);
      } // END OF /NOT-scaling/NOT-reduced/
    }   // END OF /NOT-scaling/









    /** Correction of w to take into account the dependence
        on the tangential velocity */
    if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] == 1)
    {
      for(unsigned int i = 0; i < nd; ++ i)
      {
        if(i % d == 0)
          /* w[i] = w_origin[i]/(problem->mu[(int)(i/d)]) */
          w[i] = w_origin[i] + sqrt(velocity[i+1]*velocity[i+1]+velocity[i+2]*velocity[i+2]) + sqrt(velocity[i+3]*velocity[i+3]+velocity[i+4]*velocity[i+4]);
      }

    }









    /* -------------- Build the right-hand side related to the first system (sigma = 0) -------------- */
    /* NT scaling: /scaling/ */
    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING])
    {
      /* Reduced system: /scaling/reduced/ */
      if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM])
      {






        // cblas_dcopy(m, dualConstraint, 1, rhs_CHECK_save, 1);
        // cblas_dcopy(nd, primalConstraint, 1, rhs_CHECK_save+m, 1);

        // cblas_dcopy(n_dminus2, reaction_1, 1, rhs_CHECK_save+m_plus_nd, 1);
        // cblas_dcopy(n_dminus2, reaction_2, 1, rhs_CHECK_save+m_plus_nd+n_dminus2, 1);

        // cblas_dscal(m + nd + n_dplus1, -1.0, rhs_CHECK_save, 1);





        // // For computing NT Non-reduced at the same time
        // cblas_dcopy(m, dualConstraint_CHECK, 1, rhs_NT_nonRe, 1);
        // cblas_dcopy(nd, primalConstraint_CHECK, 1, rhs_NT_nonRe+m, 1);
        // cblas_dcopy(n_dminus2, reaction_1_CHECK, 1, rhs_NT_nonRe+m_plus_nd, 1);
        // cblas_dcopy(n_dminus2, reaction_2_CHECK, 1, rhs_NT_nonRe+m_plus_nd+n_dminus2, 1);
        // cblas_dscal(m + nd + n_dplus1, -1.0, rhs_NT_nonRe, 1);

        // cblas_dcopy(m + nd + n_dplus1, rhs_NT_nonRe, 1, rhs_NT_nonRe_save, 1);










        /*  rhs = +     <==== positive sign
         *        [             M*v - H'*r + f            ]  m
         *        [ P^-1*( [u-H*v-w] - J*Q^-2*[r_1; r_2]) ]  nd
         */
        cblas_dcopy(m, dualConstraint, 1, rhs, 1);

        // Regroup r_1, r_2 into only one vector [r_1; r_2]
        cblas_dcopy(n_dminus2, reaction_1, 1, r1r2, 1);
        cblas_dcopy(n_dminus2, reaction_2, 1, r1r2+n_dminus2, 1);             // r1r2 = [r_1; r_2]

        Qinv2 = NM_create(NM_SPARSE, n_dplus1, n_dplus1);
        NM_triplet_alloc(Qinv2, 2 * d_minus_2 * d_minus_2 * n);
        NM_insert(Qinv2, Qpinv2_bar, 0, 0);
        NM_insert(Qinv2, Qpinv2_tilde, n_dminus2, n_dminus2);

        JQinv2 = NM_multiply(J, Qinv2);
        NM_gemv(-1.0, JQinv2, r1r2, 0.0, tmp_nd);           // tmp_nd = - J*Q^-2*[r_1; r_2]
        cblas_daxpy(nd, 1.0, primalConstraint, 1, tmp_nd, 1);          // tmp_nd = [u-H*v-w] - J*Q^-2*[r_1; r_2]
        NM_gemv(1.0, P_inv, tmp_nd, 0.0, tmp_nd_2);  // tmp_nd_2 = P^-1*( [u-H*v-w] - J*Q^-2*[r_1; r_2])

        cblas_dcopy(nd, tmp_nd_2, 1, rhs+m, 1);

        cblas_dcopy(m+nd, rhs, 1, rhs_NT_Re_save, 1);
      } // End of with reduced system

      else /* NOT-reduced system:  /scaling/NOT-reduced/ */
      {
        /*  rhs = -     <==== ATTENTION to negative sign
         *        [ M*v - H'*r + f ]  m         dualConstraint        = a (No have negative sign at the beginning)
         *        [  u - H*v - w   ]  nd        primalConstraint      = b
         *        [      r_1       ]  n(d-2)    complemConstraint 1 |
         *        [      r_2       ]  n(d-2)    complemConstraint 2 | = c
         */
        cblas_dcopy(m, dualConstraint, 1, rhs, 1);
        cblas_dcopy(nd, primalConstraint, 1, rhs+m, 1);

        cblas_dcopy(n_dminus2, reaction_1, 1, rhs+m_plus_nd, 1);
        assert((m_plus_nd+2*n_dminus2) == (m_plus_nd+n_dplus1)); // m + nd + n(d+1): size of rhs
        cblas_dcopy(n_dminus2, reaction_2, 1, rhs+m_plus_nd+n_dminus2, 1);

        cblas_dscal(m + nd + n_dplus1, -1.0, rhs, 1);
      } // End of /scaling/NOT-reduced/
    }   // End of /scaling/

    else /* NOT NT scaling: /NOT-scaling/ */
    {
      /* Reduced system: /NOT-scaling/reduced/ */
      if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM])
      {
        printf("\n\nNo work 2\n\n");
        break;
      } // End of /NOT-scaling/reduced/

      /* NOT-reduced system: /NOT-scaling/NOT-reduced/ */
      else
      {
        /*
         *  rhs = -     <==== ATTENTION to negative sign
         *        [ M*v - H'*r + f ]  m         dualConstraint
         *        [  u - H*v - w   ]  nd        primalConstraint
         *        [   u_1 o r_1    ]  n(d-2)    complemConstraint 1
         *        [   u_2 o r_2    ]  n(d-2)    complemConstraint 2
         */
        cblas_dcopy(m, dualConstraint, 1, rhs, 1);
        cblas_dcopy(nd, primalConstraint, 1, rhs+m, 1);

        JA_prod(velocity_1, reaction_1, n_dminus2, n, complemConstraint_1);
        JA_prod(velocity_2, reaction_2, n_dminus2, n, complemConstraint_2);
        cblas_dcopy(n_dminus2, complemConstraint_1, 1, rhs+m_plus_nd, 1);
        assert((m_plus_nd+2*n_dminus2) == (m_plus_nd+n_dplus1)); // m + nd + n(d+1): size of rhs
        cblas_dcopy(n_dminus2, complemConstraint_2, 1, rhs+m_plus_nd+n_dminus2, 1);

        cblas_dscal(m + nd + n_dplus1, -1.0, rhs, 1);
      } // End of /NOT-scaling/NOT-reduced/
    }   // End of /NOT-scaling/







    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING])
    {
      /* Solving full symmetric Newton system with NT scaling via LDLT factorization */
      NSM_linearSolverParams(Jac)->solver = NSM_HSL;
      NM_LDLT_solve(Jac, rhs, 1);


      // // For computing NT Non-reduced at the same time
      // if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM])
      // {
      //   NSM_linearSolverParams(Jac_CHECK)->solver = NSM_HSL;
      //   NM_LDLT_solve(Jac_CHECK, rhs_NT_nonRe, 1);
      // }
    }
    else
    {
      /* Solving non-symmetric Newton system without NT scaling via LU factorization */
      NM_LU_solve(Jac, rhs, 1);
    }



    /* ---------------------------- Retrieve the solutions for predictor step ---------------------------- */
    /* Reduced system */
    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM])
    {
      /*
       * Note that the results are stored in rhs, then
       *
       * Results rhs  =  [ d_globalVelocity ]  m
       *                 [ d_reaction       ]  nd
       *
       * d_reaction   = (d_reaction_0, d_reaction_bar, d_reaction_tilde)
       */
      cblas_dcopy(m, rhs, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, rhs+m, 1, tmp_nd, 1);               // tmp_nd <-- d_reaction_reduced = P'*d_reaction

      // Recover d_reaction
      NM_tgemv(1.0, P_inv, tmp_nd, 0.0, d_reaction); // d_reaction = P^-1'*d_reaction_reduced
      // PinvTx(f_NT, g_NT, wf_NT, wg_NT, n_dminus2, n, tmp_nd, d_reaction);

      extract_vector(d_reaction, nd, n, 2, 3, d_reaction_1);
      extract_vector(d_reaction, nd, n, 4, 5, d_reaction_2);


      // Recover d_velocity = (dt + dt', d_u_bar, d_u_tilde)
      NM_tgemv(-1.0, JQinv2, d_reaction, 0.0, u1u2);    // u1u2 = -(J*Q^-2)'*d_reaction
      // JQinv2Tx(f_NT, g_NT, wf_NT, wg_NT, n_dminus2, n, d_reaction, u1u2);
      // cblas_dscal(n_dplus1, -1.0, u1u2, 1);


      NM_gemv(-1.0, Qinv2, r1r2, 1.0, u1u2);            // u1u2 <-- du1du2 = -(J*Q^-2)'*d_reaction - Q^-2*[r1;r2]
      // Qinv2x(f_NT, g_NT, wf_NT, wg_NT, n_dminus2, n, r1r2, Qinv2_x);
      // cblas_daxpy(n_dplus1, -1.0, Qinv2_x, 1, u1u2, 1);



      // Recover d_u_bar & d_u_tilde
      NM_gemv(1.0, H, d_globalVelocity, 0.0, d_velocity);              // d_velocity = H*dv
      cblas_daxpy(nd, -1.0, primalConstraint, 1, d_velocity, 1);       // d_velocity = H*dv - (u-Hv-w)
      //NV_copy(Hdv, nd, d_velocity);
      extract_vector(d_velocity, nd, n, 2, 3, d_velocity_1);
      extract_vector(d_velocity, nd, n, 4, 5, d_velocity_2);

      // Recover d_t & d_t'
      for(size_t i = 0; i < n; i++)
      {
        d_velocity_1[i*d_minus_2] = u1u2[i*d_minus_2];
        d_velocity_2[i*d_minus_2] = u1u2[i*d_minus_2+n_dminus2];
      }








      // cblas_dcopy(m, d_globalVelocity, 1, rhs_CHECK, 1);
      // cblas_dcopy(nd, d_reaction, 1, rhs_CHECK+m, 1);
      // cblas_dcopy(n_dminus2, d_velocity_1, 1, rhs_CHECK+m+nd, 1);
      // cblas_dcopy(n_dminus2, d_velocity_2, 1, rhs_CHECK+m+nd+n_dminus2, 1);








      // // For computing NT Non-reduced at the same time
      // cblas_dcopy(m, rhs_NT_nonRe, 1, d_globalVelocity_CHECK, 1);
      // cblas_dcopy(nd, rhs_NT_nonRe+m, 1, d_reaction_CHECK, 1);
      // extract_vector(d_reaction_CHECK, nd, n, 2, 3, d_reaction_1_CHECK);
      // extract_vector(d_reaction_CHECK, nd, n, 4, 5, d_reaction_2_CHECK);

      // cblas_dcopy(n_dminus2, rhs_NT_nonRe+m_plus_nd, 1, d_velocity_1_CHECK, 1);
      // cblas_dcopy(n_dminus2, rhs_NT_nonRe+m_plus_nd+n_dminus2, 1, d_velocity_2_CHECK, 1);

      // for(size_t i = 0; i < n; i++)
      // {
      //   posX = i*d;
      //   posY = i*d_minus_2;
      //   d_velocity_CHECK[posX] = d_velocity_1_CHECK[posY] + d_velocity_2_CHECK[posY];
      //   d_velocity_CHECK[posX+1] = d_velocity_1_CHECK[posY + 1];
      //   d_velocity_CHECK[posX+2] = d_velocity_1_CHECK[posY + 2];
      //   d_velocity_CHECK[posX+3] = d_velocity_2_CHECK[posY + 1];
      //   d_velocity_CHECK[posX+4] = d_velocity_2_CHECK[posY + 2];
      // }


    } // End of Reduced system

    else /* NOT-reduced system */
    {
      /*
       *                 [ d_globalVelocity ]  m
       *                 [ d_reaction_0     ]  1 x n  |
       *                 [ d_reaction_bar   ]  2 x n  | = nd
       * Results rhs  =  [ d_reaction_tilde ]  2 x n  |
       *                 [ d_t              ]  1 x n  }
       *                 [ d_velocity_bar   ]  2 x n  } = n(d-2) = n_dminus2
       *                 [ d_t'             ]  1 x n  |
       *                 [ d_velocity_tilde ]  2 x n  | = n(d-2)
       *
       * d_reaction   = (d_reaction_0, d_reaction_bar   , d_reaction_tilde)
       * d_reaction_1 = (d_reaction_0, d_reaction_bar   )
       * d_reaction_2 = (d_reaction_0, d_reaction_tilde )
       *
       * d_velocity   = (d_t + d_t'  , d_velocity_bar   , d_velocity_tilde)
       * d_velocity_1 = (d_t         , d_velocity_bar   )
       * d_velocity_2 = (d_t'        , d_velocity_tilde )
       *
       */
      cblas_dcopy(m, rhs, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, rhs+m, 1, d_reaction, 1);
      extract_vector(d_reaction, nd, n, 2, 3, d_reaction_1);
      extract_vector(d_reaction, nd, n, 4, 5, d_reaction_2);

      cblas_dcopy(n_dminus2, rhs+m_plus_nd, 1, d_velocity_1, 1);
      assert((m_plus_nd+2*n_dminus2) == (m_plus_nd+n_dplus1)); // m + nd + n(d+1): size of rhs
      cblas_dcopy(n_dminus2, rhs+m_plus_nd+n_dminus2, 1, d_velocity_2, 1);

      for(size_t i = 0; i < n; i++)
      {
        posX = i*d;
        posY = i*d_minus_2;
        d_velocity[posX] = d_velocity_1[posY] + d_velocity_2[posY];
        d_velocity[posX+1] = d_velocity_1[posY + 1];
        d_velocity[posX+2] = d_velocity_1[posY + 2];
        d_velocity[posX+3] = d_velocity_2[posY + 1];
        d_velocity[posX+4] = d_velocity_2[posY + 2];
      }
    } // End of NOT-Reduced system











    /* computing the affine step-length */
    alpha_primal_1 = getStepLength(velocity_1, d_velocity_1, n_dminus2, n, gmm);
    alpha_primal_2 = getStepLength(velocity_2, d_velocity_2, n_dminus2, n, gmm);
    alpha_dual_1 = getStepLength(reaction_1, d_reaction_1, n_dminus2, n, gmm);
    alpha_dual_2 = getStepLength(reaction_2, d_reaction_2, n_dminus2, n, gmm);

    alpha_primal = fmin(alpha_primal_1, fmin(alpha_primal_2, fmin(alpha_dual_1, alpha_dual_2)));
    alpha_dual = alpha_primal;

    /* updating the gamma parameter used to compute the step-length */
    // gmm = gmmp1 + gmmp2 * fmin(alpha_primal, alpha_dual);
    gmm = gmmp1 + gmmp2 * alpha_primal;



    // // For computing NT Non-reduced at the same time
    // alpha_primal_1_CHECK = getStepLength(velocity_1_CHECK, d_velocity_1_CHECK, n_dminus2, n, gmm_CHECK);
    // alpha_primal_2_CHECK = getStepLength(velocity_2_CHECK, d_velocity_2_CHECK, n_dminus2, n, gmm_CHECK);
    // alpha_dual_1_CHECK = getStepLength(reaction_1_CHECK, d_reaction_1_CHECK, n_dminus2, n, gmm_CHECK);
    // alpha_dual_2_CHECK = getStepLength(reaction_2_CHECK, d_reaction_2_CHECK, n_dminus2, n, gmm_CHECK);

    // alpha_primal_CHECK = fmin(alpha_primal_1_CHECK, fmin(alpha_primal_2_CHECK, fmin(alpha_dual_1_CHECK, alpha_dual_2_CHECK)));
    // alpha_dual_CHECK = alpha_primal_CHECK;

    // gmm_CHECK = gmmp1 + gmmp2 * alpha_primal_CHECK;









    // if (iteration >= 0)
    // {
    //   printf("\n\n========== PRINTING FOR DEBUG 2 ==========\n");
    //   printf("iteration = %i\t1st RESULT\n", iteration);
    //   NV_sub(d_globalVelocity, d_globalVelocity_CHECK, m, diff_d_globalVelocity);
    //   NV_sub(d_velocity, d_velocity_CHECK, nd, diff_d_velocity);
    //   NV_sub(d_reaction, d_reaction_CHECK, nd, diff_d_reaction);
    //   printf("| dv - dv_CHECK | = %10.30Le\n", dnrm2l(m, diff_d_globalVelocity));
    //   printf("| du - du_CHECK | = %10.30Le\n", dnrm2l(nd, diff_d_velocity));
    //   printf("| dr - dr_CHECK | = %10.30Le\n", dnrm2l(nd, diff_d_reaction));
    //   // printf("RHS reduced | ");
    //   // NV_display(rhs, m+nd);
    //   // printf("\n\nRHS NonRe | ");
    //   // NV_display(rhs_NT_nonRe, m+nd+n_dplus1);

    //   // NM_gemv(1.0, Jac_CHECK, rhs_CHECK, -1.0, rhs_CHECK_save);
    //   // printf("1st RESULT: \n");
    //   // printf("Jac_NT_NonRe * rhs_NT_Re >= 1e-15 :\n");
    //   // for (size_t i=0; i<m + nd + n*(d+1); i++)
    //   // {
    //   //   if (rhs_CHECK_save[i] >= 1e-12)
    //   //   {
    //   //     printf("rhs[%zu] = %5.30e\n", i, rhs_CHECK_save[i]);
    //   //   }
    //   // }

    //   // NM_gemv(1.0, Jac, rhs, -1.0, rhs_NT_Re_save);
    //   // printf("|      Jac_NT_Re * rhs_NT_Re - rhs_NT_Re_save       | = %5.30Le\n", dnrm2l(m_plus_nd, rhs_NT_Re_save));
    //   // NM_gemv(1.0, Jac_CHECK, rhs_CHECK, -1.0, rhs_CHECK_save);
    //   // printf("| Jac_NT_NonRe * rhs_NT_Recover - rhs_NT_NonRe_save | = %5.30Le\n", dnrm2l(m_plus_nd, rhs_CHECK_save));
    //   // printf("\n\nJac_NT_Re * rhs_NT_Re >= 1e-15 :\n");
    //   // for (size_t i=0; i<m + nd; i++)
    //   // {
    //   //   if (rhs_NT_Re_save[i] >= 1e-15)
    //   //   {
    //   //     printf("rhs[%zu] = %5.30e\n", i, rhs_NT_Re_save[i]);
    //   //   }
    //   // }

    //   // cblas_daxpy(m, alpha_primal, d_globalVelocity, 1, globalVelocity, 1);
    //   // cblas_daxpy(nd, alpha_primal, d_velocity, 1, velocity, 1);
    //   // cblas_daxpy(nd, alpha_dual, d_reaction, 1, reaction, 1);

    //   // primalResidual(velocity, H, globalVelocity, w, primalConstraint_CHECK, &pinfeas_CHECK);
    //   // dualResidual(M, globalVelocity, H, reaction, f, dualConstraint_CHECK, &dinfeas_CHECK);

    //   // printf("\nprimalConstraint = %6.20e\n", pinfeas_CHECK);
    //   // printf("dualConstraint = %6.20e\n", dinfeas_CHECK);


    //   // NM_gemv(1.0, Jac_CHECK, rhs_NT_nonRe, -1.0, rhs_NT_nonRe_save);
    //   // printf("\n\nJac_NT_NonRe * rhs_NT_NonRe >= 1e-15 :\n");
    //   // for (size_t i=0; i<m + nd + n*(d+1); i++)
    //   // {
    //   //   if (rhs_NT_nonRe_save[i] >= 1e-15)
    //   //   {
    //   //     printf("rhs[%zu] = %5.30e\n", i, rhs_NT_nonRe_save[i]);
    //   //   }
    //   // }

    //   // printf("\n\nabs(rhs_CHECK - rhs_NT_nonRe >= 1e-15 :\n");
    //   // for (size_t i=0; i<m + nd + n*(d+1); i++)
    //   // {
    //   //   if (fabsl(rhs_CHECK[i] - rhs_NT_nonRe[i]) >= 1e-10)
    //   //   {
    //   //     printf("rhs_CHECK[%zu]    = %5.30e\n", i, rhs_CHECK[i]);
    //   //     printf("rhs_NT_nonRe[%zu] = %5.30e\n\n", i, rhs_NT_nonRe[i]);
    //   //   }
    //   // }


    //   // if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM])
    //   // {
    //   //   printf("\n\nRHS reduced =\n");
    //   //   NV_display(rhs_CHECK, m + nd + n*(d+1));
    //   // }
    //   // else
    //   // {
    //   //   printf("RHS non reduced =\n");
    //   //   NV_display(rhs, m + nd + n*(d+1));
    //   // }


    //   // cblas_daxpy(m, -alpha_primal, d_globalVelocity, 1, globalVelocity, 1);
    //   // cblas_daxpy(nd, -alpha_primal, d_velocity, 1, velocity, 1);
    //   // cblas_daxpy(nd, -alpha_dual, d_reaction, 1, reaction, 1);
    //   printf("========== END PRINTING FOR DEBUG 2 ==========\n\n");
    //   // break;
    // }










    /* -------------------------- Predictor step of Mehrotra -------------------------- */
    cblas_dcopy(nd, velocity, 1, v_plus_dv, 1);
    cblas_dcopy(nd, reaction, 1, r_plus_dr, 1);
    cblas_daxpy(nd, alpha_primal, d_velocity, 1, v_plus_dv, 1);
    cblas_daxpy(nd, alpha_dual, d_reaction, 1, r_plus_dr, 1);

    /* affine barrier parameter */
    // barr_param_a = (cblas_ddot(nd, v_plus_dv, 1, r_plus_dr, 1) / nd)*sigma;
    // barr_param_a = cblas_ddot(nd, v_plus_dv, 1, r_plus_dr, 1) / nd;
    barr_param_a = cblas_ddot(nd, v_plus_dv, 1, r_plus_dr, 1) / (n);
    //barr_param_a = complemResidualNorm(v_plus_dv, r_plus_dr, nd, n);
    // barr_param_a = (fws=='*' ? complemResidualNorm(v_plus_dv, r_plus_dr, nd, n)/n : complemResidualNorm_p(v_plus_dv, r_plus_dr, nd, n)/n);

    /* computing the centralization parameter */
    e = barr_param > sgmp1 ? fmax(1.0, sgmp2 * pow(fmin(alpha_primal, alpha_dual),2)) : sgmp3;
    sigma = fmin(1.0, pow(barr_param_a / barr_param, e))/d;




    // // For computing NT Non-reduced at the same time
    // cblas_dcopy(nd, velocity_CHECK, 1, v_plus_dv_CHECK, 1);
    // cblas_dcopy(nd, reaction_CHECK, 1, r_plus_dr_CHECK, 1);
    // cblas_daxpy(nd, alpha_primal_CHECK, d_velocity_CHECK, 1, v_plus_dv_CHECK, 1);
    // cblas_daxpy(nd, alpha_dual_CHECK, d_reaction_CHECK, 1, r_plus_dr_CHECK, 1);

    // barr_param_a_CHECK = cblas_ddot(nd, v_plus_dv_CHECK, 1, r_plus_dr_CHECK, 1) / (n);
    // e_CHECK = barr_param_CHECK > sgmp1 ? fmax(1.0, sgmp2 * pow(fmin(alpha_primal_CHECK, alpha_dual_CHECK),2)) : sgmp3;
    // sigma_CHECK = fmin(1.0, pow(barr_param_a_CHECK / barr_param_CHECK, e_CHECK))/d;













    /* -------------------------- Corrector step of Mehrotra -------------------------- */
    cblas_dcopy(m, dualConstraint, 1, rhs, 1);

    /* NT scaling: /scaling/ */
    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING] > 0)
    {
      /* Reduced system: /scaling/reduced/ */
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM])
      {
        /* Right-hand side for symmetric Newton system with NT scaling */

        /* 1st terms: [r_1; r_2] */
        // already in var r1r2 = [r_1; r_2]


        /* 2nd terms: 2 * barr_param * sigma * [(u_1)^-1; (u_2)^-1]  */
        JA_inv(velocity_1, n_dminus2, n, velocity_1_inv);                         // velocity_1_inv    = (u_1)^-1
        JA_inv(velocity_2, n_dminus2, n, velocity_2_inv);                         // velocity_2_inv    = (u_2)^-1
        cblas_dcopy(n_dminus2, velocity_1_inv, 1, u1u2, 1);
        cblas_dcopy(n_dminus2, velocity_2_inv, 1, u1u2+n_dminus2, 1);             // u1u2 = [(u_1)^-1; (u_2)^-1]
        cblas_dscal(n_dplus1, 2 * barr_param * sigma, u1u2, 1);                   // u1u2 = 2nd term


        /* 3rd terms: (Qp_bar * u_1)^-1 o [(Qp_bar * du_1) o (Qpinv_bar * dr_1)]  */
        NM_gemv(1.0, Qp_bar, velocity_1, 0.0, Qpvelocity_1);                      // Qpvelocity_1 = Qp_bar * u_1
        NM_gemv(1.0, Qp_tilde, velocity_2, 0.0, Qpvelocity_2);                    // Qpvelocity_2 = Qp_tilde * u_2
        JA_inv(Qpvelocity_1, n_dminus2, n, Qpvelocity_1_inv);                     // Qpvelocity_1_inv = (Qp_bar * u_1)^-1
        JA_inv(Qpvelocity_2, n_dminus2, n, Qpvelocity_2_inv);                     // Qpvelocity_2_inv = (Qp_tilde * u_2)^-1

        NM_gemv(1.0, Qp_bar, d_velocity_1, 0.0, d_velocity_1t);                      // d_velocity_1t     = Qp_bar * du_1
        NM_gemv(1.0, Qp_tilde, d_velocity_2, 0.0, d_velocity_2t);                    // d_velocity_2t     = Qp_tilde * du_2
        NM_gemv(1.0, Qpinv_bar, d_reaction_1, 0.0, d_reaction_1t);                   // d_reaction_1t     = Qpinv_bar * dr_1
        NM_gemv(1.0, Qpinv_tilde, d_reaction_2, 0.0, d_reaction_2t);                 // d_reaction_2t     = Qpinv_tilde * dr_2

        JA_prod(d_velocity_1t, d_reaction_1t, n_dminus2, n, dvdr_jprod_1);           // dvdr_jprod_1      = (Qp_bar * du_1) o (Qpinv_bar * dr_1)
        JA_prod(d_velocity_2t, d_reaction_2t, n_dminus2, n, dvdr_jprod_2);           // dvdr_jprod_2      = (Qp_tilde * du_2) o (Qpinv_tilde * dr_2)

        JA_prod(Qpvelocity_1_inv, dvdr_jprod_1, n_dminus2, n, Qpvelocity_inv_dvdr_1);   // Qpvelocity_inv_dvdr_1 = 3rd term_1
        JA_prod(Qpvelocity_2_inv, dvdr_jprod_2, n_dminus2, n, Qpvelocity_inv_dvdr_2);   // Qpvelocity_inv_dvdr_2 = 3rd term_2

        cblas_dcopy(n_dminus2, Qpvelocity_inv_dvdr_1, 1, dvdr_jprod, 1);
        cblas_dcopy(n_dminus2, Qpvelocity_inv_dvdr_2, 1, dvdr_jprod+n_dminus2, 1);            // dvdr_jprod = 3rd term


        /* updated rhs = P^-1*([u-H*v-w] - J*Q^-2*{[r_1; r_2] - 2nd terms} - J*Q^-1*3rd terms ) */
        NV_sub(r1r2, u1u2, n_dplus1, tmp_n_dplus1);                  // tmp_n_dplus1 = 1st terms - 2nd terms
        NM_gemv(-1.0, JQinv2, tmp_n_dplus1, 0.0, tmp_nd);            // tmp_nd = - J*Q^-2*{1st terms - 2nd terms}


        Qinv = NM_create(NM_SPARSE, n_dplus1, n_dplus1);
        NM_triplet_alloc(Qinv, 2 * d_minus_2 * d_minus_2 * n);
        NM_insert(Qinv, Qpinv_bar, 0, 0);
        NM_insert(Qinv, Qpinv_tilde, n_dminus2, n_dminus2);
        JQinv = NM_multiply(J, Qinv);
        NM_gemv(-1.0, JQinv, dvdr_jprod, 0.0, tmp_nd_2);            // tmp_nd_2 = - J*Q^-1*3rd terms


        cblas_daxpy(nd, 1.0, primalConstraint, 1, tmp_nd, 1);   // tmp_nd   =        [u-H*v-w] - J*Q^-2*{1st terms - 2nd terms}
        cblas_daxpy(nd, 1.0, tmp_nd, 1, tmp_nd_2, 1);           // tmp_nd_2 =        [u-H*v-w] - J*Q^-2*{1st terms - 2nd terms} - J*Q^-1*3rd terms

        NM_gemv(1.0, P_inv, tmp_nd_2, 0.0, tmp_nd);             // tmp_nd   = P^-1*( [u-H*v-w] - J*Q^-2*{1st terms - 2nd terms} - J*Q^-1*3rd terms )

        cblas_dcopy(nd, tmp_nd, 1, rhs+m, 1);


        cblas_dcopy(m+nd, rhs, 1, rhs_NT_Re_save, 1);





        // cblas_dcopy(m, dualConstraint, 1, rhs_CHECK_save, 1);
        // cblas_dcopy(nd, primalConstraint, 1, rhs_CHECK_save+m, 1);

        // double *iden_1 = JA_iden(n_dminus2, n);
        // double *iden_2 = JA_iden(n_dminus2, n);
        // cblas_dscal(n_dminus2, -2 * barr_param * sigma, iden_1, 1);
        // cblas_dcopy(n_dminus2, iden_1, 1, iden_2, 1);

        // NM_gemv(1.0, Qpinv_bar, reaction_1, 0.0, QpinvReaction_1);
        // NM_gemv(1.0, Qpinv_tilde, reaction_2, 0.0, QpinvReaction_2);
        // JA_prod(Qpvelocity_1, QpinvReaction_1, n_dminus2, n, Qp_uQpinv_r_1);
        // JA_prod(Qpvelocity_2, QpinvReaction_2, n_dminus2, n, Qp_uQpinv_r_2);

        // cblas_daxpy(n_dminus2, 1.0, dvdr_jprod_1, 1, iden_1, 1);
        // cblas_daxpy(n_dminus2, 1.0, dvdr_jprod_2, 1, iden_2, 1);
        // cblas_daxpy(n_dminus2, 1.0, Qp_uQpinv_r_1, 1, iden_1, 1);
        // cblas_daxpy(n_dminus2, 1.0, Qp_uQpinv_r_2, 1, iden_2, 1);

        // cblas_dcopy(n_dminus2, iden_1, 1, rhs_CHECK_save+m_plus_nd, 1);
        // cblas_dcopy(n_dminus2, iden_2, 1, rhs_CHECK_save+m_plus_nd+n_dminus2, 1);
        // cblas_dscal(m + nd + n_dplus1, -1.0, rhs_CHECK_save, 1);

        // free(iden_1); free(iden_2);




        // // For computing NT Non-reduced at the same time
        // cblas_dcopy(m, dualConstraint_CHECK, 1, rhs_NT_nonRe, 1);
        // cblas_dcopy(nd, primalConstraint_CHECK, 1, rhs_NT_nonRe+m, 1);
        // JA_inv(velocity_1_CHECK, n_dminus2, n, velocity_1t_inv_CHECK);                         // velocity_1t_inv    = u_1^-1
        // JA_inv(velocity_2_CHECK, n_dminus2, n, velocity_2t_inv_CHECK);                         // velocity_2t_inv    = u_2^-1
        // cblas_dscal(n_dminus2, 2 * barr_param_CHECK * sigma_CHECK, velocity_1t_inv_CHECK, 1);      // velocity_1t_inv = 2nd term
        // cblas_dscal(n_dminus2, 2 * barr_param_CHECK * sigma_CHECK, velocity_2t_inv_CHECK, 1);      // velocity_2t_inv = 2nd term

        // NM_gemv(1.0, Qp_bar_CHECK, velocity_1_CHECK, 0.0, Qpvelocity_1_CHECK);                      // Qpvelocity_1 = Qp_bar * u_1
        // NM_gemv(1.0, Qp_tilde_CHECK, velocity_2_CHECK, 0.0, Qpvelocity_2_CHECK);                    // Qpvelocity_2 = Qp_tilde * u_2
        // JA_inv(Qpvelocity_1_CHECK, n_dminus2, n, Qpvelocity_1_inv_CHECK);                     // Qpvelocity_1_inv = (Qp_bar * u_1)^-1
        // JA_inv(Qpvelocity_2_CHECK, n_dminus2, n, Qpvelocity_2_inv_CHECK);                     // Qpvelocity_2_inv = (Qp_tilde * u_2)^-1

        // NM_gemv(1.0, Qp_bar_CHECK, d_velocity_1_CHECK, 0.0, d_velocity_1t_CHECK);                      // d_velocity_1t     = Qp_bar * du_1
        // NM_gemv(1.0, Qp_tilde_CHECK, d_velocity_2_CHECK, 0.0, d_velocity_2t_CHECK);                    // d_velocity_2t     = Qp_tilde * du_2
        // NM_gemv(1.0, Qpinv_bar_CHECK, d_reaction_1_CHECK, 0.0, d_reaction_1t_CHECK);                   // d_reaction_1t     = Qpinv_1 * dr_1
        // NM_gemv(1.0, Qpinv_tilde_CHECK, d_reaction_2_CHECK, 0.0, d_reaction_2t_CHECK);                 // d_reaction_2t     = Qpinv_2 * dr_2

        // JA_prod(d_velocity_1t_CHECK, d_reaction_1t_CHECK, n_dminus2, n, dvdr_jprod_1_CHECK);           // dvdr_jprod_1      = (Qp_bar * du_1) o (Qpinv_bar * dr_1)
        // JA_prod(d_velocity_2t_CHECK, d_reaction_2t_CHECK, n_dminus2, n, dvdr_jprod_2_CHECK);           // dvdr_jprod_2      = (Qp_tilde * du_2) o (Qpinv_tilde * dr_2)

        // JA_prod(Qpvelocity_1_inv_CHECK, dvdr_jprod_1_CHECK, n_dminus2, n, Qpvelocity_inv_dvdr_1_CHECK);
        // JA_prod(Qpvelocity_2_inv_CHECK, dvdr_jprod_2_CHECK, n_dminus2, n, Qpvelocity_inv_dvdr_2_CHECK);
        // NM_gemv(1.0, Qp_bar_CHECK, Qpvelocity_inv_dvdr_1_CHECK, 0.0, QpQpvelocity_inv_dvdr_1_CHECK);   // QpQpvelocity_inv_dvdr_1_CHECK = 3rd term_1
        // NM_gemv(1.0, Qp_tilde_CHECK, Qpvelocity_inv_dvdr_2_CHECK, 0.0, QpQpvelocity_inv_dvdr_2_CHECK);// QpQpvelocity_inv_dvdr_2_CHECK = 3rd term_2

        // NV_sub(reaction_1_CHECK, velocity_1t_inv_CHECK, n_dminus2, tmp1_CHECK);                   // tmp1 = r_1 - 2nd term
        // NV_sub(reaction_2_CHECK, velocity_2t_inv_CHECK, n_dminus2, tmp2_CHECK);                   // tmp2 = r_2 - 2nd term
        // NV_add(tmp1_CHECK, QpQpvelocity_inv_dvdr_1_CHECK, n_dminus2, complemConstraint_1_CHECK);
        // NV_add(tmp2_CHECK, QpQpvelocity_inv_dvdr_2_CHECK, n_dminus2, complemConstraint_2_CHECK);

        // cblas_dcopy(n_dminus2, complemConstraint_1_CHECK, 1, rhs_NT_nonRe+m_plus_nd, 1);
        // cblas_dcopy(n_dminus2, complemConstraint_2_CHECK, 1, rhs_NT_nonRe+m_plus_nd+n_dminus2, 1);
        // cblas_dscal(m + nd + n_dplus1, -1.0, rhs_NT_nonRe, 1);

      } // End of Reduced system

      else /* NOT-reduced system: /scaling/NOT-reduced/ */
      {
        cblas_dcopy(nd, primalConstraint, 1, rhs+m, 1);

        /* Right-hand side for symmetric Newton system with NT scaling */
        /* 2nd terms: 2 * barr_param * sigma * (u_1)^-1  */
        JA_inv(velocity_1, n_dminus2, n, velocity_1t_inv);                         // velocity_1t_inv    = u_1^-1
        JA_inv(velocity_2, n_dminus2, n, velocity_2t_inv);                         // velocity_2t_inv    = u_2^-1
        cblas_dscal(n_dminus2, 2 * barr_param * sigma, velocity_1t_inv, 1);      // velocity_1t_inv = 2nd term
        cblas_dscal(n_dminus2, 2 * barr_param * sigma, velocity_2t_inv, 1);      // velocity_2t_inv = 2nd term

        /* 3rd terms: Qp_bar * [ (Qp_bar * u_1)^-1 o (Qp_bar * du_1) o (Qpinv_1 * dr_1) ] */
        NM_gemv(1.0, Qp_bar, velocity_1, 0.0, Qpvelocity_1);                      // Qpvelocity_1 = Qp_bar * u_1
        NM_gemv(1.0, Qp_tilde, velocity_2, 0.0, Qpvelocity_2);                    // Qpvelocity_2 = Qp_tilde * u_2
        JA_inv(Qpvelocity_1, n_dminus2, n, Qpvelocity_1_inv);                     // Qpvelocity_1_inv = (Qp_bar * u_1)^-1
        JA_inv(Qpvelocity_2, n_dminus2, n, Qpvelocity_2_inv);                     // Qpvelocity_2_inv = (Qp_tilde * u_2)^-1

        NM_gemv(1.0, Qp_bar, d_velocity_1, 0.0, d_velocity_1t);                      // d_velocity_1t     = Qp_bar * du_1
        NM_gemv(1.0, Qp_tilde, d_velocity_2, 0.0, d_velocity_2t);                    // d_velocity_2t     = Qp_tilde * du_2
        if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] == SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP)
        {
          QNTpinvz(velocity_1, reaction_1, d_reaction_1, n_dminus2, n, d_reaction_1t); // d_reaction_1t     = Qpinv_bar * dr_1
          QNTpinvz(velocity_2, reaction_2, d_reaction_2, n_dminus2, n, d_reaction_2t); // d_reaction_2t     = Qpinv_tilde * dr_2
        }
        else if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] == SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F)
        {
          NM_gemv(1.0, Qpinv_bar, d_reaction_1, 0.0, d_reaction_1t);                   // d_reaction_1t     = Qpinv_bar * dr_1
          NM_gemv(1.0, Qpinv_tilde, d_reaction_2, 0.0, d_reaction_2t);                 // d_reaction_2t     = Qpinv_tilde * dr_2
        }

        JA_prod(d_velocity_1t, d_reaction_1t, n_dminus2, n, dvdr_jprod_1);           // dvdr_jprod_1      = (Qp_bar * du_1) o (Qpinv_bar * dr_1)
        JA_prod(d_velocity_2t, d_reaction_2t, n_dminus2, n, dvdr_jprod_2);           // dvdr_jprod_2      = (Qp_tilde * du_1) o (Qpinv_tilde * dr_1)

        JA_prod(Qpvelocity_1_inv, dvdr_jprod_1, n_dminus2, n, Qpvelocity_inv_dvdr_1); // Qpvelocity_inv_dvdr_1 = (Qp_bar * u_1)^-1 o (Qp_bar * du_1) o (Qpinv_bar * dr_1)
        JA_prod(Qpvelocity_2_inv, dvdr_jprod_2, n_dminus2, n, Qpvelocity_inv_dvdr_2); // Qpvelocity_inv_dvdr_2 = (Qp_tilde * u_2)^-1 o (Qp_tilde * du_2) o (Qpinv_tilde * dr_1)

        NM_gemv(1.0, Qp_bar, Qpvelocity_inv_dvdr_1, 0.0, QpQpvelocity_inv_dvdr_1);   // QpQpvelocity_inv_dvdr_1 = 3rd term_1
        NM_gemv(1.0, Qp_tilde, Qpvelocity_inv_dvdr_2, 0.0, QpQpvelocity_inv_dvdr_2);// QpQpvelocity_inv_dvdr_2 = 3rd term_2




        // NM_gemv(1.0, NM_multiply(Arrow_repr(d_velocity_1t, n_dminus2, n), Arrow_repr(d_reaction_1t, n_dminus2, n)), velocity_1t_inv, 0.0, dvdr_jprod_1);
        // NM_gemv(1.0, NM_multiply(Arrow_repr(d_velocity_2t, n_dminus2, n), Arrow_repr(d_reaction_2t, n_dminus2, n)), velocity_2t_inv, 0.0, dvdr_jprod_2);



        /* updated rhs = r_1 - 2nd term + 3rd term */
        NV_sub(reaction_1, velocity_1t_inv, n_dminus2, tmp1);                   // tmp1 = r_1 - 2nd term
        NV_sub(reaction_2, velocity_2t_inv, n_dminus2, tmp2);                   // tmp2 = r_2 - 2nd term
        // NV_add(tmp1, dvdr_jprod_1, n_dminus2, complemConstraint_1);
        // NV_add(tmp2, dvdr_jprod_2, n_dminus2, complemConstraint_2);
        NV_add(tmp1, QpQpvelocity_inv_dvdr_1, n_dminus2, complemConstraint_1);
        NV_add(tmp2, QpQpvelocity_inv_dvdr_2, n_dminus2, complemConstraint_2);

        /* Update only 2 complementarities of rhs  */
        cblas_dcopy(n_dminus2, complemConstraint_1, 1, rhs+m_plus_nd, 1);
        cblas_dcopy(n_dminus2, complemConstraint_2, 1, rhs+m_plus_nd+n_dminus2, 1);
        cblas_dscal(m + nd + n_dplus1, -1.0, rhs, 1);
      } // End of NOT-Reduced system
    } // End of NT scaling

    else/* NOT scaling: /NOT-scaling/ */
    {
      /* Reduced system: /NOT-scaling/reduced/ */
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM])
      {
        printf("\n\nNo work 3\n\n");
        break;
      } // End of Reduced system

      else /* NOT-reduced system: /NOT-scaling/NOT-reduced/ */
      {
        cblas_dcopy(nd, primalConstraint, 1, rhs+m, 1);

        /* Right-hand side for non-symmetric Newton system without NT scaling */
        /* 1st terms: u_1 o r_1  */
        JA_prod(velocity_1, reaction_1, n_dminus2, n, vr_jprod_1);
        JA_prod(velocity_2, reaction_2, n_dminus2, n, vr_jprod_2);

        /* 2nd terms: 2 * barr_param * sigma * e  */
        iden = JA_iden(n_dminus2, n);
        cblas_dscal(n_dminus2, 2 * barr_param * sigma, iden, 1);

        /* 3rd terms: du_1 o dr_1  */
        JA_prod(d_velocity_1, d_reaction_1, n_dminus2, n, dvdr_jprod_1);
        JA_prod(d_velocity_2, d_reaction_2, n_dminus2, n, dvdr_jprod_2);

        /* updated rhs = 1st term - 2nd term + 3rd term */
        NV_sub(vr_jprod_1, iden, n_dminus2, tmp1);
        NV_sub(vr_jprod_2, iden, n_dminus2, tmp2);
        NV_add(tmp1, dvdr_jprod_1, n_dminus2, complemConstraint_1);
        NV_add(tmp2, dvdr_jprod_2, n_dminus2, complemConstraint_2);
        free(iden);

        /* Update only 2 complementarities of rhs  */
        cblas_dcopy(n_dminus2, complemConstraint_1, 1, rhs+m_plus_nd, 1);
        cblas_dcopy(n_dminus2, complemConstraint_2, 1, rhs+m_plus_nd+n_dminus2, 1);
        cblas_dscal(m + nd + n_dplus1, -1.0, rhs, 1);
      } // End of NOT-Reduced system
    } // End of NOT scaling














    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING] > 0)
    {
      /* Solving full symmetric Newton system with NT scaling via LDLT factorization */
      NM_LDLT_solve(Jac, rhs, 1);


      // // For computing NT Non-reduced at the same time
      // if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM])
      // {
      //   NM_LDLT_solve(Jac_CHECK, rhs_NT_nonRe, 1);
      // }


    }
    else
    {
      /* Solving non-symmetric Newton system without NT scaling via LU factorization */
      NM_LU_solve(Jac, rhs, 1);
    }




    // NM_clear(Jac);
    // free(Jac);





    /* ---------------------------- Retrieve the solutions for corrector step ---------------------------- */
    /* Reduced system */
    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM])
    {
      cblas_dcopy(m, rhs, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, rhs+m, 1, tmp_nd, 1);         // tmp_nd <-- d_reaction_reduced = P'*d_reaction

      // Recover d_reaction
      NM_tgemv(1.0, P_inv, tmp_nd, 0.0, d_reaction); // d_reaction = P^-1'*d_reaction_reduced
      // PinvTx(f_NT, g_NT, wf_NT, wg_NT, n_dminus2, n, tmp_nd, d_reaction);
      extract_vector(d_reaction, nd, n, 2, 3, d_reaction_1);
      extract_vector(d_reaction, nd, n, 4, 5, d_reaction_2);

      // Recover d_velocity = (dt + dt', d_u_bar, d_u_tilde)
      NM_tgemv(-1.0, JQinv2, d_reaction, 0.0, u1u2);    // u1u2 = -(J*Q^-2)'*d_reaction
      // JQinv2Tx(f_NT, g_NT, wf_NT, wg_NT, n_dminus2, n, d_reaction, u1u2);
      // cblas_dscal(n_dplus1, -1.0, u1u2, 1);


      // Recover d_u_bar & d_u_tilde
      NM_gemv(1.0, H, d_globalVelocity, 0.0, d_velocity);              // d_velocity = H*dv
      cblas_daxpy(nd, -1.0, primalConstraint, 1, d_velocity, 1);       // d_velocity = H*dv - (u-Hv-w)
      extract_vector(d_velocity, nd, n, 2, 3, d_velocity_1);
      extract_vector(d_velocity, nd, n, 4, 5, d_velocity_2);

      // Recover d_t & d_t'
      NM_gemv(-1.0, Qinv2, tmp_n_dplus1, 1.0, u1u2);            // u1u2 = -(J*Q^-2)'*d_reaction - Q^-2*(1st terms - 2nd terms)  // tmp_n_dplus1 = 1st terms - 2nd terms
      // Qinv2x(f_NT, g_NT, wf_NT, wg_NT, n_dminus2, n, tmp_n_dplus1, Qinv2_x);
      // cblas_daxpy(n_dplus1, -1.0, Qinv2_x, 1, u1u2, 1);
      NM_gemv(-1.0, Qinv, dvdr_jprod, 1.0, u1u2);           // u1u2 = -(J*Q^-2)'*d_reaction - Q^-2*(1st terms - 2nd terms) - Q^-1*3rd term // dvdr_jprod = 3rd term

      for(size_t i = 0; i < n; i++)
      {
        d_velocity_1[i*d_minus_2] = u1u2[i*d_minus_2];
        d_t[i] = u1u2[i*d_minus_2];

        d_velocity_2[i*d_minus_2] = u1u2[i*d_minus_2+n_dminus2];
        d_t_prime[i] = u1u2[i*d_minus_2+n_dminus2];
      }






      // cblas_dcopy(m, d_globalVelocity, 1, rhs_CHECK, 1);
      // cblas_dcopy(nd, d_reaction, 1, rhs_CHECK+m, 1);
      // cblas_dcopy(n_dminus2, d_velocity_1, 1, rhs_CHECK+m+nd, 1);
      // cblas_dcopy(n_dminus2, d_velocity_2, 1, rhs_CHECK+m+nd+n_dminus2, 1);





      // // For computing NT Non-reduced at the same time
      // cblas_dcopy(m, rhs_NT_nonRe, 1, d_globalVelocity_CHECK, 1);
      // cblas_dcopy(nd, rhs_NT_nonRe+m, 1, d_reaction_CHECK, 1);
      // extract_vector(d_reaction_CHECK, nd, n, 2, 3, d_reaction_1_CHECK);
      // extract_vector(d_reaction_CHECK, nd, n, 4, 5, d_reaction_2_CHECK);

      // cblas_dcopy(n_dminus2, rhs_NT_nonRe+m_plus_nd, 1, d_velocity_1_CHECK, 1);
      // cblas_dcopy(n_dminus2, rhs_NT_nonRe+m_plus_nd+n_dminus2, 1, d_velocity_2_CHECK, 1);

      // for(size_t i = 0; i < n; i++)
      // {
      //   posX = i*d;
      //   posY = i*d_minus_2;
      //   d_velocity_CHECK[posX] = d_velocity_1_CHECK[posY] + d_velocity_2_CHECK[posY];
      //   d_velocity_CHECK[posX+1] = d_velocity_1_CHECK[posY + 1];
      //   d_velocity_CHECK[posX+2] = d_velocity_1_CHECK[posY + 2];
      //   d_velocity_CHECK[posX+3] = d_velocity_2_CHECK[posY + 1];
      //   d_velocity_CHECK[posX+4] = d_velocity_2_CHECK[posY + 2];
      //   d_t_CHECK[i] = d_velocity_1_CHECK[posY];
      //   d_t_prime_CHECK[i] = d_velocity_2_CHECK[posY];
      // }







    } // End of Reduced system

    else /* NOT-reduced system */
    {


      cblas_dcopy(m, rhs, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, rhs+m, 1, d_reaction, 1);
      extract_vector(d_reaction, nd, n, 2, 3, d_reaction_1);
      extract_vector(d_reaction, nd, n, 4, 5, d_reaction_2);

      cblas_dcopy(n_dminus2, rhs+m_plus_nd, 1, d_velocity_1, 1);
      assert((m_plus_nd+2*n_dminus2) == (m_plus_nd+n_dplus1)); // m + nd + n(d+1): size of rhs
      cblas_dcopy(n_dminus2, rhs+m_plus_nd+n_dminus2, 1, d_velocity_2, 1);

      for(size_t i = 0; i < n; i++)
      {
        posX = i*d;
        posY = i*d_minus_2;
        d_velocity[posX] = d_velocity_1[posY] + d_velocity_2[posY];
        d_velocity[posX+1] = d_velocity_1[posY + 1];
        d_velocity[posX+2] = d_velocity_1[posY + 2];
        d_velocity[posX+3] = d_velocity_2[posY + 1];
        d_velocity[posX+4] = d_velocity_2[posY + 2];
        d_t[i] = d_velocity_1[posY];
        d_t_prime[i] = d_velocity_2[posY];
      }
    } // End of NOT-Reduced system












    // if (iteration == 0)
    // {
    //   printf("\n\n========== PRINTING FOR DEBUG ==========\n");
    //   printf("\n\niteration: %d\n", iteration);
    //   // printf("Vector v:\n");
    //   // NM_vector_display(globalVelocity,m);
    //   // printf("\n\nVector u:\n");
    //   // NM_vector_display(velocity,nd);
    //   // printf("\n\nVector r:\n");
    //   // NM_vector_display(reaction,nd);
    //   // printf("\n\nVector t:\n");
    //   // NM_vector_display(t,n);
    //   // printf("\n\nVector t_prime:\n");
    //   // NM_vector_display(t_prime,n);
    //   // printf("\n\nVector u_1:\n");
    //   // NM_vector_display(velocity_1,n_dminus2);
    //   // printf("\n\nVector u_2:\n");
    //   // NM_vector_display(velocity_2,n_dminus2);
    //   // printf("\n\nVector r_1:\n");
    //   // NM_vector_display(reaction_1,n_dminus2);
    //   // printf("\n\nVector r_2:\n");
    //   // NM_vector_display(reaction_2,n_dminus2);

    //   printf("Vector dv:\n");
    //   NM_vector_display(d_globalVelocity,m);
    //   printf("\n\nVector du:\n");
    //   NM_vector_display(d_velocity,nd);
    //   printf("\n\nVector dr:\n");
    //   NM_vector_display(d_reaction,nd);
    //   printf("\n\nVector dt:\n");
    //   NM_vector_display(d_t,n);
    //   printf("\n\nVector dt_prime:\n");
    //   NM_vector_display(d_t_prime,n);
    //   printf("\n\nVector du_1:\n");
    //   NM_vector_display(d_velocity_1,n_dminus2);
    //   printf("\n\nVector du_2:\n");
    //   NM_vector_display(d_velocity_2,n_dminus2);
    //   printf("\n\nVector dr_1:\n");
    //   NM_vector_display(d_reaction_1,n_dminus2);
    //   printf("\n\nVector dr_2:\n");
    //   NM_vector_display(d_reaction_2,n_dminus2);

    //   printf("\n\nJac:\n");
    //   NM_display(Jac);
    //   printf("\n\nrhs:\n");
    //   NM_vector_display(rhs, m_plus_nd+n_dplus1);

    //   // printf("\n\nalpha_primal_1 = %f\n", alpha_primal_1);
    //   // printf("\n\nalpha_primal_2 = %f\n", alpha_primal_2);
    //   // printf("\n\nalpha_primal = %f\n", alpha_primal);
    //   // printf("\n\nalpha_dual_1 = %f\n", alpha_dual_1);
    //   // printf("\n\nalpha_dual_2 = %f\n", alpha_dual_2);
    //   // printf("\n\nalpha_dual = %f\n", alpha_dual);

    //   printf("\n\ncomplem_1 = %9.16f\n", complem_1);
    //   printf("\n\ncomplem_2 = %9.16f\n", complem_2);
    //   // printf("\n\nerror = %f\n", error);
    //   printf("========== END PRINTING FOR DEBUG ==========\n\n");
    //   break;
    // }





















    /* computing the affine step-length */
    alpha_primal_1 = getStepLength(velocity_1, d_velocity_1, n_dminus2, n, gmm);
    alpha_primal_2 = getStepLength(velocity_2, d_velocity_2, n_dminus2, n, gmm);
    alpha_dual_1   = getStepLength(reaction_1, d_reaction_1, n_dminus2, n, gmm);
    alpha_dual_2   = getStepLength(reaction_2, d_reaction_2, n_dminus2, n, gmm);

    alpha_primal = fmin(alpha_primal_1, fmin(alpha_primal_2, fmin(alpha_dual_1, alpha_dual_2)));
    alpha_dual = alpha_primal;

    /* updating the gamma parameter used to compute the step-length */
    gmm = gmmp1 + gmmp2 * alpha_primal;




    // // For computing NT Non-reduced at the same time
    // alpha_primal_1_CHECK = getStepLength(velocity_1_CHECK, d_velocity_1_CHECK, n_dminus2, n, gmm_CHECK);
    // alpha_primal_2_CHECK = getStepLength(velocity_2_CHECK, d_velocity_2_CHECK, n_dminus2, n, gmm_CHECK);
    // alpha_dual_1_CHECK = getStepLength(reaction_1_CHECK, d_reaction_1_CHECK, n_dminus2, n, gmm_CHECK);
    // alpha_dual_2_CHECK = getStepLength(reaction_2_CHECK, d_reaction_2_CHECK, n_dminus2, n, gmm_CHECK);

    // alpha_primal_CHECK = fmin(alpha_primal_1_CHECK, fmin(alpha_primal_2_CHECK, fmin(alpha_dual_1_CHECK, alpha_dual_2_CHECK)));
    // alpha_dual_CHECK = alpha_primal_CHECK;

    // gmm_CHECK = gmmp1 + gmmp2 * alpha_primal_CHECK;
















    /* print out all useful parameters */
    numerics_printf_verbose(-1, "| %3i%c| %9.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e |",
                            iteration, fws, relgap, pinfeas, dinfeas, u1dotr1, u2dotr2, complem_1, complem_2, full_error, barr_param, alpha_primal, alpha_dual, sigma,
                            cblas_dnrm2(m, d_globalVelocity, 1)/cblas_dnrm2(m, globalVelocity, 1),
                            cblas_dnrm2(nd, d_velocity, 1)/cblas_dnrm2(nd, velocity, 1),
                            cblas_dnrm2(nd, d_reaction, 1)/cblas_dnrm2(nd, reaction, 1));


















    /* ----- Update variables ----- */
    cblas_daxpy(m, alpha_primal, d_globalVelocity, 1, globalVelocity, 1);
    cblas_daxpy(nd, alpha_primal, d_velocity, 1, velocity, 1);
    cblas_daxpy(nd, alpha_dual, d_reaction, 1, reaction, 1);
    cblas_daxpy(n, alpha_dual, d_t, 1, t, 1);
    cblas_daxpy(n, alpha_dual, d_t_prime, 1, t_prime, 1);




    // // For computing NT Non-reduced at the same time
    // cblas_daxpy(m, alpha_primal_CHECK, d_globalVelocity_CHECK, 1, globalVelocity_CHECK, 1);
    // cblas_daxpy(nd, alpha_primal_CHECK, d_velocity_CHECK, 1, velocity_CHECK, 1);
    // cblas_daxpy(nd, alpha_dual_CHECK, d_reaction_CHECK, 1, reaction_CHECK, 1);
    // cblas_daxpy(n, alpha_dual_CHECK, d_t_CHECK, 1, t_CHECK, 1);
    // cblas_daxpy(n, alpha_dual_CHECK, d_t_prime_CHECK, 1, t_prime_CHECK, 1);










    if (NV_isnan(globalVelocity, m) | NV_isnan(velocity, nd) | NV_isnan(reaction, nd))
    {
      hasNotConverged = 2;
      break;
    }




    // if (iteration >= 0)
    // {

    //   printf("\n\n========== PRINTING FOR DEBUG 3 ==========\n");
    //   printf("iteration = %i\t2nd RESULT\n", iteration);
    //   NV_sub(d_globalVelocity, d_globalVelocity_CHECK, m, diff_d_globalVelocity);
    //   NV_sub(d_velocity, d_velocity_CHECK, nd, diff_d_velocity);
    //   NV_sub(d_reaction, d_reaction_CHECK, nd, diff_d_reaction);
    //   printf("| dv - dv_CHECK | = %10.30Le\n", dnrm2l(m, diff_d_globalVelocity));
    //   printf("| du - du_CHECK | = %10.30Le\n", dnrm2l(nd, diff_d_velocity));
    //   printf("| dr - dr_CHECK | = %10.30Le\n", dnrm2l(nd, diff_d_reaction));
    //   // printf("RHS reduced | ");
    //   // NV_display(rhs, m+nd);
    //   // printf("\n\nRHS NonRe | ");
    //   // NV_display(rhs_NT_nonRe, m+nd+n_dplus1);
    //   numerics_printf_verbose(-1, "| %3i%c| %9.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e |",
    //                         iteration, fws, relgap_CHECK, pinfeas_CHECK, dinfeas_CHECK, u1dotr1_CHECK, u2dotr2_CHECK, complem_1_CHECK, complem_2_CHECK, full_error, barr_param_CHECK, alpha_primal_CHECK, alpha_dual_CHECK, sigma_CHECK);


    //   // NM_gemv(1.0, Jac_CHECK, rhs_CHECK, -1.0, rhs_CHECK_save);
    //   // printf("2nd RESULT: \n");
    //   // printf("Jac_NT_NonRe * rhs_NT_Re >= 1e-15 :\n");
    //   // for (size_t i=0; i<m + nd + n*(d+1); i++)
    //   // {
    //   //   if (rhs_CHECK_save[i] >= 1e-12)
    //   //   {
    //   //     printf("rhs[%zu] = %5.30e\n", i, rhs_CHECK_save[i]);
    //   //   }
    //   // }

    //   // NM_gemv(1.0, Jac, rhs, -1.0, rhs_NT_Re_save);
    //   // printf("|       Jac_NT_Re * rhs_NT_Re - rhs_NT_Re_save      | = %5.30Le\n", dnrm2l(m_plus_nd, rhs_NT_Re_save));
    //   // NM_gemv(1.0, Jac_CHECK, rhs_CHECK, -1.0, rhs_CHECK_save);
    //   // printf("| Jac_NT_NonRe * rhs_NT_Recover - rhs_NT_NonRe_save | = %5.30Le\n", dnrm2l(m_plus_nd, rhs_CHECK_save));
    //   // printf("\n\nJac_NT_Re * rhs_NT_Re >= 1e-15 :\n");
    //   // for (size_t i=0; i<m + nd; i++)
    //   // {
    //   //   if (rhs_NT_Re_save[i] >= 1e-15)
    //   //   {
    //   //     printf("rhs[%zu] = %5.30e\n", i, rhs_NT_Re_save[i]);
    //   //   }
    //   // }

    //   // primalResidual(velocity, H, globalVelocity, w, primalConstraint_CHECK, &pinfeas_CHECK);
    //   // dualResidual(M, globalVelocity, H, reaction, f, dualConstraint_CHECK, &dinfeas_CHECK);

    //   // printf("\nprimalConstraint = %6.20e\n", pinfeas_CHECK);
    //   // printf("dualConstraint = %6.20e\n", dinfeas_CHECK);
    //   // NV_display(rhs_CHECK_save, m + nd + n*(d+1));

    //   // NM_gemv(1.0, Jac_CHECK, rhs_NT_nonRe, -1.0, rhs_NT_nonRe_save);
    //   // printf("\n\nJac_NT_NonRe * rhs_NT_NonRe >= 1e-15 :\n");
    //   // for (size_t i=0; i<m + nd + n*(d+1); i++)
    //   // {
    //   //   if (rhs_NT_nonRe_save[i] >= 1e-15)
    //   //   {
    //   //     printf("rhs[%zu] = %5.30e\n", i, rhs_NT_nonRe_save[i]);
    //   //   }
    //   // }

    //   // printf("\n\nabs(rhs_CHECK - rhs_NT_nonRe >= 1e-15 :\n");
    //   // for (size_t i=0; i<m + nd + n*(d+1); i++)
    //   // {
    //   //   if (fabsl(rhs_CHECK[i] - rhs_NT_nonRe[i]) >= 1e-10)
    //   //   {
    //   //     printf("rhs_CHECK[%zu]    = %5.30e\n", i, rhs_CHECK[i]);
    //   //     printf("rhs_NT_nonRe[%zu] = %5.30e\n\n", i, rhs_NT_nonRe[i]);
    //   //   }
    //   // }

    //   // if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM])
    //   // {
    //   //   printf("\n\nRHS reduced =\n");
    //   //   NV_display(rhs_CHECK, m + nd + n*(d+1));
    //   // }
    //   // else
    //   // {
    //   //   printf("RHS non reduced =\n");
    //   //   NV_display(rhs, m + nd + n*(d+1));
    //   // }

    //   printf("========== END PRINTING FOR DEBUG 3 ==========\n\n");
    //   if (iteration == 10)break;
    // }




    NM_clear(Jac);
    free(Jac);
    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM])
    {
      // NM_clear(Jac_CHECK);
      // free(Jac_CHECK);
    }

    iteration++;
  } // end of while loop


  /* -------------------------- Return to original variables -------------------------- */
  // TO DO: If we need this point, plz uncomment
  // NM_gemv(1.0, P_mu_inv, velocity, 0.0, data->original_point->velocity);
  // NM_gemv(1.0, P_mu, reaction, 0.0, data->original_point->reaction);
  // cblas_dcopy(m, globalVelocity, 1, data->original_point->globalVelocity, 1);

  options->dparam[SICONOS_DPARAM_RESIDU] = full_error; //NV_max(error, 4);
  options->iparam[SICONOS_IPARAM_ITER_DONE] = iteration;




  clock_t t2 = clock();
  long clk_tck = CLOCKS_PER_SEC;

  /* writing data in a Matlab file */
  // if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
  // {
  //   char matlab_file_name[256];
  //   sprintf(matlab_file_name,"sigma_nc-%d-.m", problem->numberOfContacts);
  //   matlab_file = fopen(matlab_file_name, "w");
  //   printInteresProbMatlabFile(iteration, globalVelocity, velocity, reaction, d, n, m, (double)(t2-t1)/(double)clk_tck, matlab_file);
  //   fclose(matlab_file);
  // }





  if(internal_allocation)
  {
    grfc3d_IPM_free(problem,options);
  }

  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING])
  {
    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM])
    {
      NM_clear(minus_M);
      free(minus_M);
    }
    else
    {

    }
    NM_clear(Qp_bar);
    free(Qp_bar);
    // NM_clear(Qp2_bar);
    // free(Qp2_bar);
    NM_clear(Qp_tilde);
    free(Qp_tilde);
    // NM_clear(Qp2_tilde);
    // free(Qp2_tilde);
  }







  NM_clear(H_origin);
  free(H_origin);
  NM_clear(minus_H);
  free(minus_H);
  NM_clear(H);
  free(H);

  NM_clear(J);
  free(J);

  // NM_clear(Jac); // already
  // free(Jac);

  free(t);
  free(t_prime);
  free(velocity_1);
  free(velocity_2);
  free(reaction_1);
  free(reaction_2);
  free(d_globalVelocity);
  free(d_velocity);
  free(d_velocity_1);
  free(d_velocity_2);
  free(d_reaction);
  free(d_reaction_1);
  free(d_reaction_2);
  free(d_t);
  free(d_t_prime);




















  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
    fclose(iterates);

  //  fclose(dfile);

  *info = hasNotConverged;








} // end of grfc3d_IPM











/* initialize solver (allocate memory) */
void grfc3d_IPM_init(GlobalRollingFrictionContactProblem* problem, SolverOptions* options)
{
  unsigned int m = problem->M->size0;
  unsigned int nd = problem->H->size1;
  unsigned int d = problem->dimension;

  unsigned int n = (int)(nd / d);
  unsigned int n_dminus2 = n*(d-2);


  // TO DO: need to review this. Ex. if options->dWork existed and options->dWorkSize is unexpected, then what will we do ?
  if(!options->dWork || options->dWorkSize != (size_t)(m + nd + n*(d+1)))
  {
    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM])
    {
      options->dWork = (double*)calloc(m + nd, sizeof(double));
      options->dWorkSize = m + nd;
    }
    else
    {
      options->dWork = (double*)calloc(m + nd + n*(d+1), sizeof(double));
      options->dWorkSize = m + nd + n*(d+1);
    }

  }


  /* ------------- initialize starting point ------------- */
  options->solverData=(Grfc3d_IPM_data *)malloc(sizeof(Grfc3d_IPM_data));
  Grfc3d_IPM_data * data = (Grfc3d_IPM_data *)options->solverData;


  /* --------- allocate memory for IPM point ----------- */
  data->starting_point = (IPM_point*)malloc(sizeof(IPM_point));

  /* 1. v */
  data->starting_point->globalVelocity = (double*)calloc(m, sizeof(double));
  for(unsigned int i = 0; i < m; ++ i)
    data->starting_point->globalVelocity[i] = 0.01;

  /* 2. u */
  data->starting_point->velocity = (double*)calloc(nd, sizeof(double));
  for(unsigned int i = 0; i < nd; ++ i)
  {
    data->starting_point->velocity[i] = 0.001;
    if(i % d == 0)
      data->starting_point->velocity[i] = 3.0;
  }

  /* 3. r */
  data->starting_point->reaction = (double*)calloc(nd, sizeof(double));
  for(unsigned int i = 0; i < nd; ++ i)
  {
    data->starting_point->reaction[i] = 0.04;
    if(i % d == 0)
      data->starting_point->reaction[i] = 0.5;
  }


  /* original point which is not changed by the matrix P_mu */
  data->original_point = (IPM_point*)malloc(sizeof(IPM_point));
  data->original_point->globalVelocity = (double*)calloc(m, sizeof(double));
  data->original_point->velocity = (double*)calloc(nd, sizeof(double));
  data->original_point->reaction = (double*)calloc(nd, sizeof(double));



  /* --------- allocate memory for IPM grfc3d point ----------- */
  // data->grfc3d_point = (IPM_grfc3d_point*)malloc(sizeof(IPM_grfc3d_point));

  /* 1. t, t_prime */
  // data->grfc3d_point->t = (double*)calloc(n, sizeof(double));
  // data->grfc3d_point->t_prime = (double*)calloc(n, sizeof(double));
  // for(unsigned int i = 0; i < n; ++ i)
  // {
  //   data->grfc3d_point->t[i] = 1.0;
  //   data->grfc3d_point->t_prime[i] = 1.0;
  // }

  /*
   * 2. velocity_1 = (t, u_bar), velocity_2 = (t_prime, u_tilde)
   * We will assign these variables later, in the while loop
   */
  // data->grfc3d_point->velocity_1 = (double*)calloc(n_dminus2, sizeof(double));
  // data->grfc3d_point->velocity_2 = (double*)calloc(n_dminus2, sizeof(double));

  /*
   * 3. reaction_1 = (r0, r_bar), reaction_2 = (r0, r_tilde)
   * We will assign these variables later, in the while loop
   */
  // data->grfc3d_point->reaction_1 = (double*)calloc(n_dminus2, sizeof(double));
  // data->grfc3d_point->reaction_2 = (double*)calloc(n_dminus2, sizeof(double));


  /* ------ initialize the change of variable matrix P_mu ------- */
  data->P_mu = (IPM_change_of_variable*)malloc(sizeof(IPM_change_of_variable));
  data->P_mu->mat = NM_create(NM_SPARSE, nd, nd);
  NM_triplet_alloc(data->P_mu->mat, nd);
  data->P_mu->mat->matrix2->origin = NSM_TRIPLET;
  for(unsigned int i = 0; i < nd; ++i)
  {
    if(i % d == 0)
      /* NM_entry(data->P_mu->mat, i, i, 1. / problem->mu[(int)(i/d)]); */
      NM_entry(data->P_mu->mat, i, i, 1.);
    else if(i % d == 1 || i % d == 2)
      /* NM_entry(data->P_mu->mat, i, i, 1.); */
      NM_entry(data->P_mu->mat, i, i, problem->mu[(int)(i/d)]);
    else
      /* NM_entry(data->P_mu->mat, i, i, 1.); */
      NM_entry(data->P_mu->mat, i, i, problem->mu_r[(int)(i/d)]);
  }


  /* ------ initialize the inverse P_mu_inv of the change of variable matrix P_mu ------- */
  data->P_mu->inv_mat = NM_create(NM_SPARSE, nd, nd);
  NM_triplet_alloc(data->P_mu->inv_mat, nd);
  data->P_mu->inv_mat->matrix2->origin = NSM_TRIPLET;
  for(unsigned int i = 0; i < nd; ++i)
  {
    if(i % d == 0)
      /* NM_entry(data->P_mu->inv_mat, i, i, problem->mu[(int)(i/d)]); */
      NM_entry(data->P_mu->inv_mat, i, i, 1.);
    else if(i % d == 1 || i % d == 2)
      /* NM_entry(data->P_mu->inv_mat, i, i, 1.); */
      NM_entry(data->P_mu->inv_mat, i, i, 1.0/problem->mu[(int)(i/d)]);
    else
      /* NM_entry(data->P_mu->inv_mat, i, i, 1.); */
      NM_entry(data->P_mu->inv_mat, i, i, 1.0/problem->mu_r[(int)(i/d)]);
  }


  /* ------ initial parameters initialization ---------- */
  data->internal_params = (IPM_internal_params*)malloc(sizeof(IPM_internal_params));
  data->internal_params->alpha_primal = 1.0;
  data->internal_params->alpha_dual = 1.0;
  data->internal_params->sigma = 0.1;
  data->internal_params->barr_param = 1.0;


  /* ----- temporary vaults initialization ------- */
  data->tmp_vault_m = (double**)malloc(2 * sizeof(double*));
  for(unsigned int i = 0; i < 2; ++i)
    data->tmp_vault_m[i] = (double*)calloc(m, sizeof(double));

  data->tmp_vault_nd = (double**)malloc(10 * sizeof(double*));
  for(unsigned int i = 0; i < 10; ++i)
    data->tmp_vault_nd[i] = (double*)calloc(nd, sizeof(double));

  data->tmp_vault_n_dminus2 = (double**)malloc(25 * sizeof(double*));
  for(unsigned int i = 0; i < 25; ++i)
    data->tmp_vault_n_dminus2[i] = (double*)calloc(n_dminus2, sizeof(double));

  data->tmp_vault_n = (double**)malloc(2 * sizeof(double*));
  for(unsigned int i = 0; i < 2; ++i)
    data->tmp_vault_n[i] = (double*)calloc(n, sizeof(double));
} // end of grfc3d_IPM_init



/* deallocate memory */
void grfc3d_IPM_free(GlobalRollingFrictionContactProblem* problem, SolverOptions* options)
{
  if(options->dWork)
  {
    free(options->dWork);
    options->dWork = NULL;
    options->dWorkSize = 0;
  }

  if(options->solverData)
  {
    Grfc3d_IPM_data * data = (Grfc3d_IPM_data *)options->solverData;

    free(data->starting_point->globalVelocity);
    data->starting_point->globalVelocity = NULL;

    free(data->starting_point->velocity);
    data->starting_point->velocity = NULL;

    free(data->starting_point->reaction);
    data->starting_point->reaction = NULL;

    free(data->starting_point);
    data->starting_point = NULL;

    free(data->original_point->globalVelocity);
    data->original_point->globalVelocity = NULL;

    free(data->original_point->velocity);
    data->original_point->velocity = NULL;

    free(data->original_point->reaction);
    data->original_point->reaction = NULL;

    free(data->original_point);
    data->original_point = NULL;

    NM_clear(data->P_mu->mat);
    free(data->P_mu->mat);
    data->P_mu->mat = NULL;

    NM_clear(data->P_mu->inv_mat);
    free(data->P_mu->inv_mat);
    data->P_mu->inv_mat = NULL;

    free(data->P_mu);
    data->P_mu = NULL;

    for(unsigned int i = 0; i < 2; ++i)
      free(data->tmp_vault_m[i]);
    free(data->tmp_vault_m);
    data->tmp_vault_m = NULL;

    for(unsigned int i = 0; i < 10; ++i)
      free(data->tmp_vault_nd[i]);
    free(data->tmp_vault_nd);
    data->tmp_vault_nd = NULL;

    for(unsigned int i = 0; i < 25; ++i)
      free(data->tmp_vault_n_dminus2[i]);
    free(data->tmp_vault_n_dminus2);
    data->tmp_vault_n_dminus2 = NULL;

    for(unsigned int i = 0; i < 2; ++i)
      free(data->tmp_vault_n[i]);
    free(data->tmp_vault_n);
    data->tmp_vault_n = NULL;

    // free(data->grfc3d_point->velocity_1);
    // data->grfc3d_point->velocity_1 = NULL;

    // free(data->grfc3d_point->velocity_2);
    // data->grfc3d_point->velocity_2 = NULL;

    // free(data->grfc3d_point->reaction_1);
    // data->grfc3d_point->reaction_1 = NULL;

    // free(data->grfc3d_point->reaction_2);
    // data->grfc3d_point->reaction_2 = NULL;

    // free(data->grfc3d_point->t);
    // data->grfc3d_point->t = NULL;

    // free(data->grfc3d_point->t_prime);
    // data->grfc3d_point->t_prime = NULL;

    // free(data->grfc3d_point);

    free(data->internal_params);
    data->internal_params = NULL;
  }

  free(options->solverData);
  options->solverData = NULL;

} // end of grfc3d_IPM_free


/* setup default solver parameters */
void grfc3d_IPM_set_default(SolverOptions* options)
{

  options->iparam[SICONOS_IPARAM_MAX_ITER] = 300;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO] = SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_NO;

  /* 0: convex case;  1: non-smooth case */
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] = 0;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING] = 1;
  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] = SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP;
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] = SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_F;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE] = 0;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REDUCED_SYSTEM] = 1;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING] = 0;

  options->dparam[SICONOS_DPARAM_TOL] = 1e-10;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1] = 1e-5;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2] = 3.;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3] = 1.;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1] = 0.9;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2] = 0.09; //0.09

} // end of grfc3d_IPM_set_default

