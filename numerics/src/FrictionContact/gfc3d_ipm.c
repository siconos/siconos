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
#include <time.h>


const char* const   SICONOS_GLOBAL_FRICTION_3D_IPM_STR = "GFC3D IPM";


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

// if (alphaL < 0.1) printf("Cone %d: alpha = %Le\n",i,alphaL);

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
  // printf("rn = %e, tol = %e\n", rn, tol);
}

void primalResidual_type(const double * velocity, NumericsMatrix * H, const double * globalVelocity, const double * w,
                    double * out, double * rnorm, const double tol, const int type)
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
  rn = fmax(rn, NV_norm_type(nd, velocity, type));
  rn = fmax(rn, NV_norm_type(nd, w, type));
  *rnorm = (rn > tol ? NV_norm_type(nd, out, type) : NV_norm_type(nd, out, type));
}

static void primalResidual2(const double * velocity, NumericsMatrix * H, const double * globalVelocity, const double * w,
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
  *rnorm = (rn > tol ? cblas_dnrm2(nd, out, 1)/rn : cblas_dnrm2(nd, out, 1));

  /* *rnorm = cblas_dnrm2(nd, out, 1);  */
  // printf("rn = %e, tol = %e\n", rn, tol);
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
void dualResidual_type(NumericsMatrix * M, const double * globalVelocity, NumericsMatrix * H, const double * reaction, const double * f,
                  double * out, double * rnorm, const double tol, const int type)
{
  double m = H->size1;
  double *HTr = (double*)calloc(m, sizeof(double));
  double rn;

  NM_gemv(1.0, M, globalVelocity, 0.0, out);
  rn = NV_norm_type(m, out, type);
  cblas_daxpy(m, -1.0, f, 1, out, 1);
  NM_tgemv(1.0, H, reaction, 0.0, HTr);
  cblas_daxpy(m, -1.0, HTr, 1, out, 1);
  rn = fmax(rn, NV_norm_type(m, f, type));
  rn = fmax(rn, NV_norm_type(m, HTr, type));
  *rnorm = (rn > tol ? NV_norm_type(m, out, type) : NV_norm_type(m, out, type));
  free(HTr);
}

double NV_norm_type(const unsigned int vecSize, const double * const vec, const int type)
{
  double norm = -1;

  if (type == NORM_2) // L-2 norm
  {
    norm = cblas_dnrm2(vecSize, vec, 1);
  }

  else if (type == NORM_INF) // L-inf norm
  {
    int maxIndex = cblas_idamax(vecSize, vec, 1);
    norm = fabs(vec[maxIndex]);
  }

  else
  {
    fprintf(stderr, "NV_norm_type: type = %d is undefined.\n", type);
    exit(EXIT_FAILURE);
  }

  return norm;
}

double xdoty_type(const unsigned int varsCount, const unsigned int vecSize, const double * x, const double * y, const int type)
{
  double xdoty = -1;

  if (type == NORM_2)
    xdoty = cblas_ddot(vecSize, x, 1, y, 1);

  else if (type == NORM_INF)
    for (int i = 0; i<varsCount; i++)
    {
      xdoty = fmax(xdoty, fabs(cblas_ddot(3, x+i*3, 1, y+i*3, 1)));
    }

  else if (type == NORM_2_INF)
  {
    double *xioyi = (double*)calloc(3, sizeof(double));
    for (int i = 0; i<varsCount; i++)
    {
      JA_prod(x+i*3, y+i*3, 3, 1, xioyi);
      xdoty = fmax(xdoty, cblas_dnrm2(3, xioyi, 1));
    }
    free(xioyi);
  }

  else
  {
    fprintf(stderr, "xdoty_type: type = %d is undefined.\n", type);
    exit(EXIT_FAILURE);
  }

  return xdoty;
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

/* Returns the type-norm of the complementarity residual vector = type-norm of the Jordan product velocity o reaction  */
double complemResidualNorm_type(const double * const velocity, const double * const reaction,
                           const unsigned int vecSize, const unsigned int varsCount, const int type)
{
  double * resid = NULL;
  double norm = -1;
  unsigned int d = vecSize/varsCount;

  if (type == NORM_2)
  {
    resid = (double*)calloc(vecSize, sizeof(double));
    JA_prod(velocity, reaction, vecSize, varsCount, resid);
    norm = cblas_dnrm2(vecSize, resid, 1);
  }

  else if (type == NORM_2_INF)
  {
    resid = (double*)calloc(d, sizeof(double));

    for (int i=0; i<vecSize; i+=d)
    {
      JA_prod(velocity+i, reaction+i, d, 1, resid);
      norm = fmax(norm, cblas_dnrm2(d, resid, 1));
    }
  }

  else
  {
    numerics_error("complemResidualNorm_type", "unknown norm type");
  }

  if (resid) free(resid);
  return norm;
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
double projectionError(const double * velocity, const double * reaction, const unsigned int nc, const double tol)
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
   relative_scaling = fmax(norm_u, norm_r); relative_scaling = 1.;
   // if(relative_scaling > tol)
   if(relative_scaling > DBL_EPSILON)
     out = out/relative_scaling;
   return out;
}


/* Computation of the projection error |proj(r-u)|/max{|r|, |u|} */
static double projectionError2(const double * velocity, const double * reaction, const unsigned int nc, const double tol)
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
       // worktmp[0] = reaction[3*ic] -  worktmp[0];
       // worktmp[1] = reaction[3*ic+1] -  worktmp[1];
       // worktmp[2] = reaction[3*ic+2] -  worktmp[2];
       out +=  worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2];
     }
   out = sqrt(out);
   norm_u = cblas_dnrm2(3*nc, velocity, 1);
   norm_r = cblas_dnrm2(3*nc, reaction, 1);
   relative_scaling = fmax(norm_u, norm_r); relative_scaling = 1.;
   // if(relative_scaling > tol)
   if(relative_scaling > DBL_EPSILON)
     out = out/relative_scaling;
   return out;
}

/* Computation of the projection error |r - proj(r-u)|_inf */
double projectionError_norm_infinity_conic(const double * velocity, const double * reaction, const unsigned int nc)
{
  double worktmp[3];
  double out = -1.0;

  for(int ic = 0 ; ic < nc ; ic++)
  {
    worktmp[0] = reaction[3*ic]   -  velocity[3*ic];
    worktmp[1] = reaction[3*ic+1] -  velocity[3*ic+1];
    worktmp[2] = reaction[3*ic+2] -  velocity[3*ic+2];
    projectionOnCone(worktmp, 1.0);
    worktmp[0] = reaction[3*ic]   -  worktmp[0];
    worktmp[1] = reaction[3*ic+1] -  worktmp[1];
    worktmp[2] = reaction[3*ic+2] -  worktmp[2];
    out =  fmax(out, worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2]);
  }

  return sqrt(out);
}


/* Computation of the projection error max { |ri - proj(ri-ui)| / max { |ri|, |ui| } } */
double projectionError_relative_norm_infinity_conic(const double * velocity, const double * reaction, const unsigned int nc)
{
  double worktmp[3];
  double out = -1.0;
  double nr = -1.0, nu = -1.0, error;
  int idx = 0;

  for(int ic = 0 ; ic < nc ; ic++)
  {
    idx = 3*ic;

    nr = cblas_dnrm2(3, reaction + idx, 1);
    nu = cblas_dnrm2(3, velocity + idx, 1);

    worktmp[0] = reaction[idx]   -  velocity[idx];
    worktmp[1] = reaction[idx+1] -  velocity[idx+1];
    worktmp[2] = reaction[idx+2] -  velocity[idx+2];
    projectionOnCone(worktmp, 1.0);
    worktmp[0] = reaction[idx]   -  worktmp[0];
    worktmp[1] = reaction[idx+1] -  worktmp[1];
    worktmp[2] = reaction[idx+2] -  worktmp[2];

    error = cblas_dnrm2(3, worktmp, 1) / fmax(nr, nu);

    out =  fmax(out, error);
  }

  return out;
}


/* Computation of the projection error on dual cone |u - proj(u-r)|_inf */
double projectionError_dual_norm_infinity_conic(const double * velocity, const double * reaction, const unsigned int nc)
{
  double worktmp[3];
  double out = 0.0;

  for(int ic = 0 ; ic < nc ; ic++)
  {
    worktmp[0] = velocity[3*ic]   - reaction[3*ic];
    worktmp[1] = velocity[3*ic+1] - reaction[3*ic+1];
    worktmp[2] = velocity[3*ic+2] - reaction[3*ic+2];
    projectionOnDualCone(worktmp, 1.0);
    worktmp[0] = velocity[3*ic]   -  worktmp[0];
    worktmp[1] = velocity[3*ic+1] -  worktmp[1];
    worktmp[2] = velocity[3*ic+2] -  worktmp[2];
    out =  fmax(out, worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2]);
  }

  return sqrt(out);
}


/* Computation of the projection error with recomputation of u and v */
double projectionError_based_reaction_norm_infinity_conic(NumericsMatrix *H, NumericsMatrix *M,
                    const double * f, const double * w, const double * reaction, const unsigned int varsCount,
                    const int on_dual_cone)
{
  size_t m = M->size0;
  size_t d = 3;
  size_t nd = varsCount*d;

  double *globalVelocity = (double*)calloc(m, sizeof(double));
  double *velocity = (double*)calloc(nd, sizeof(double));

  // Re-compute v = M\(H'*r+f)
  cblas_dcopy(m, f, 1, globalVelocity, 1);
  NM_tgemv(1.0, H, reaction, 1.0, globalVelocity);
  NM_Cholesky_solve(M, globalVelocity, 1);

  // Re-compute u = Hv + w + phi(Hv + w)
  cblas_dcopy(nd, w, 1, velocity, 1);
  NM_gemv(1, H, globalVelocity, 1.0, velocity);
  for (int i = 0; i < nd; i += d)
  {
    velocity[i] += cblas_dnrm2(2, velocity + i + 1, 1);
  }

  double out = -1;

  if (on_dual_cone)
  {
    out = projectionError_dual_norm_infinity_conic(velocity, reaction, varsCount);
  }
  else
  {
    out = projectionError_norm_infinity_conic(velocity, reaction, varsCount);
  }

  free(globalVelocity);
  free(velocity);

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

  fprintf(file,"H = [\n");
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
static void printIteresProbMatlabFile(int iteration, double * v, double * u, double * r, double * dv, double * du, double * dr, int d, int n, int m, FILE * file)
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

  // fprintf(file,"pinfeas(%3i) = %20.16e;\n",iteration+1,pinfeas);
  // fprintf(file,"dinfeas(%3i) = %20.16e;\n",iteration+1,dinfeas);
  // fprintf(file,"udotr(%3i) = %20.16e;\n",iteration+1,udotr);
  // fprintf(file,"residu(%3i) = %20.16e;\n",iteration+1,totalresidual);
  // fprintf(file,"phiu(%3i) = %20.16e;\n",iteration+1,norm_phiu);
  // // fprintf(file,"alpha_diff_phi(%3i) = %20.16e;\n",iteration+1,alpha_diff_phi);
  // fprintf(file,"alpha_u(%3i) = %20.16e;\n",iteration+1,alpha_diff_phi);



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


void classify_BNRT(const double * velocity, const double * reaction, const unsigned int vecSize, const unsigned int varsCount,
                   int *nB, int *nN, int *nR, int *nT)
{
  size_t d = (size_t)(vecSize / varsCount);
  if (d != 3) {
    fprintf(stderr,
            "classify_BNRT: This function ONLY supports for 3-dimensional model.\n");
    exit(EXIT_FAILURE);
  }

  double somme[3];
  double dur[3];
  double cesp = sqrt(DBL_EPSILON);
  double ns, ndur;
  *nB = *nN = *nR = *nT = 0;
  for (int i = 0; i < varsCount; i++)
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
      if (dur[0] >= cblas_dnrm2(2,dur+1,1))       *nN +=1;
      else if (-dur[0] >= cblas_dnrm2(2,dur+1,1)) *nB +=1;
      else                                        *nR +=1;
    }
  }

  *nT = varsCount - *nB - *nN - *nR;
}


/* Return the classification BNRT of the input: orignal u (uN;uT) and r (rN; rT)
 * We need first a change of velocity : u_tilde = (uN + mu*|uT|; mu*uT)
 * these u_tilde and r belong to friction cones
 */
void classify_BNRT_velocity_original(const double *mu, const double * velocity, const double * reaction, const unsigned int vecSize, const unsigned int varsCount,
                   int *nB, int *nN, int *nR, int *nT)
{
  size_t d = (size_t)(vecSize / varsCount);
  if (d != 3) {
    fprintf(stderr,
            "classify_BNRT_velocity_original: This function ONLY supports for 3-dimensional model.\n");
    exit(EXIT_FAILURE);
  }

  double *velocity_no_mu = (double*)calloc(vecSize, sizeof(double)); // = (uN + mu*|uT|; mu*uT)
  double *reaction_no_mu = (double*)calloc(vecSize, sizeof(double)); // = (rN          ; rT/mu)


  for (unsigned int i = 0; i<vecSize; i+=d)
  {
    velocity_no_mu[i] = velocity[i] + mu[i/d]*cblas_dnrm2(2, velocity+i+1, 1);
    velocity_no_mu[i+1] = mu[i/d]*velocity[i+1];
    velocity_no_mu[i+2] = mu[i/d]*velocity[i+2];

    reaction_no_mu[i] = reaction[i];
    reaction_no_mu[i+1] = reaction[i+1]/mu[i/d];
    reaction_no_mu[i+2] = reaction[i+2]/mu[i/d];
  }

  classify_BNRT(velocity_no_mu, reaction_no_mu, vecSize, varsCount, nB, nN, nR, nT);

  free(velocity_no_mu);
  free(reaction_no_mu);
}


/* Return the classification BNRT of the input: modified u (uN + mu*|uT|; uT) and r (rN; rT)
 * We need first a change of velocity : u_tilde = (uN + mu*|uT|; mu*uT)
 * these u_tilde and r belong to friction cones
 */
void classify_BNRT_velocity_modified(const double *mu, const double * velocity, const double * reaction, const unsigned int vecSize, const unsigned int varsCount,
                   int *nB, int *nN, int *nR, int *nT)
{
  size_t d = (size_t)(vecSize / varsCount);
  if (d != 3) {
    fprintf(stderr,
            "classify_BNRT_velocity_modified: This function ONLY supports for 3-dimensional model.\n");
    exit(EXIT_FAILURE);
  }

  double *velocity_no_mu = (double*)calloc(vecSize, sizeof(double)); // = (uN + mu*|uT|; mu*uT)
  double *reaction_no_mu = (double*)calloc(vecSize, sizeof(double)); // = (rN          ; rT/mu)


  for (unsigned int i = 0; i<vecSize; i+=d)
  {
    velocity_no_mu[i] = velocity[i];
    velocity_no_mu[i+1] = mu[i/d]*velocity[i+1];
    velocity_no_mu[i+2] = mu[i/d]*velocity[i+2];

    reaction_no_mu[i] = reaction[i];
    reaction_no_mu[i+1] = reaction[i+1]/mu[i/d];
    reaction_no_mu[i+2] = reaction[i+2]/mu[i/d];
  }

  classify_BNRT(velocity_no_mu, reaction_no_mu, vecSize, varsCount, nB, nN, nR, nT);

  free(velocity_no_mu);
  free(reaction_no_mu);
}


// By defaut, array setR is set to 0. If a cone is classified in R, the corresponding index in setR is set to 1.
void classify_indices_R(const double * velocity, const double * reaction, const unsigned int vecSize, const unsigned int varsCount,
                   int *nR, int *setR)
{
  // Note that setR must be allocated varsCount memory space
  // Reset setR
  for (int i = 0; i < varsCount; i++)
  {
    setR[i] = 0;
  }
  *nR = 0;

  size_t d = (size_t)(vecSize / varsCount);
  if (d != 3) {
    fprintf(stderr,
            "classify_indices_R: This function ONLY supports for 3-dimensional model.\n");
    exit(EXIT_FAILURE);
  }

  double somme[3];
  double dur[3];
  double cesp = sqrt(DBL_EPSILON);
  double ns, ndur;

  for (int i = 0; i < varsCount; i++)
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
      if (dur[0] >= cblas_dnrm2(2,dur+1,1))
        continue;
      else if (-dur[0] >= cblas_dnrm2(2,dur+1,1))
        continue;
      else
      {
        setR[i] = 1;
        (*nR)++;
      }
    }
  }

  // for (int i = 0; i < varsCount; i++)
  // {
  //   printf("classify_indices_R: setR[%d] = %d\n", i, setR[i]);
  // }
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

// cl: classify
void printBlockVec(double * vec, int vecSize, int sizeBlock, int cl)
{
  if (cl > 0)
  {
    for(unsigned int i=0; i<vecSize; i++)
    {
      if (i%sizeBlock == 0) printf("\n%4i-%4i: ", i, i+sizeBlock-1);
      printf(" %.8e,", vec[i]);
      if ((i+1)%sizeBlock == 0) printf("\t\t x0 - |xb| = %.2e", vec[((i+1)/sizeBlock-1)*sizeBlock] - cblas_dnrm2(sizeBlock-1, vec+((i+1)/sizeBlock-1)*sizeBlock+1, 1));
    }
  }

  else
  {
    for(unsigned int i=0; i<vecSize; i++)
    {
      if (i%sizeBlock == 0) printf("\n%4i-%4i: ", i, i+sizeBlock-1);
      printf(" %.15e", vec[i]);
    }
  }
  return;
}



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

/* check the solution of the linear system A*x = b  */
/* return |A*x[p:p+n-1]-b|/|b[p:p+n-1]|            */
static double relative_error_linear_system_solution(NumericsMatrix* const A, const double * x, const double * b, int b_size, int p, int n)
{
  double out;
  double norm_b = cblas_dnrm2(n, b+p, 1);
  double * rhs = (double*)calloc(b_size, sizeof(double));
  cblas_dcopy(b_size, b, 1, rhs, 1);
  NM_gemv(1.0, A, x, -1.0, rhs);
  out = cblas_dnrm2(n, rhs+p, 1);
  printf("out=%9.2e\n", out);
  out = out/norm_b;
  free(rhs);
  return out;
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

  // int type = NORM_2;
  int type = NORM_INF;




  // the size of the problem detection
  unsigned int m = problem->M->size0;
  unsigned int nd = problem->H->size1;
  unsigned int d = problem->dimension;
  unsigned int n = problem->numberOfContacts;



  // int internal_allocation_tmp=0;
  // if(!options->dWork || (options->dWorkSize != (size_t)(m + nd + nd)))
  // {
  //   gfc3d_IPM_init(problem, options);
  //   internal_allocation_tmp = 1;
  // }


  // for (double mu=1e-2; mu>1e-10; mu/=10.)
  // {
  //   // gfc3d_IPM_fixed(problem, reaction, velocity, globalVelocity, info, options, 1e-3, problem_name);
  //   // printf("mu = %e\n", mu);
  //   gfc3d_IPM_fixed(problem, reaction, velocity, globalVelocity, info, options, mu, problem_name);
  //   if (*info == 0) printf("\nSUCCESS\n\n");
  //   else printf("\nFAILURE\n\n");
  // }

  // if (internal_allocation_tmp)
  // {
  //   gfc3d_IPM_free(problem,options);
  // }
  // return;







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
  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] > 0)
  {
    for(int i = 0; i < n ; i++) problem->mu[i] = 0.3;
  }

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

  // printf("nnz H = %8zu density = %9.4f\n",NM_nnz(problem->H), NM_nnz(problem->H)/1.0/nd/m);
  printf("nnz H = %8zu density = %9.4f, \t mu = %.2e\n",NM_nnz(problem->H), NM_nnz(problem->H)/1.0/nd/m, problem->mu[0]);

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

  // change of variable to eliminate the friction coefficients: H_tilde --> H and w_tilde --> w
  NumericsMatrix *H = NM_multiply(P_mu, H_tilde);
  NM_gemv(1.0, P_mu, w_tilde, 0.0, w);

  //cs_print(NM_triplet(H),0);
  //printf("#################### NORM(w) = %g    NORM(f) = %g\n", NV_norm_2(w, nd), NV_norm_2(f,m));

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
  double alpha_tmp;

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
      if (i % d == 0)
	reaction[i] = 0.1;
      else
	reaction[i] = 0.01;

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
  double gmmp1 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1];
  double gmmp2 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2];

  int hasNotConverged = 1;
  unsigned int iteration = 0;
  double pinfeas = 1e300; double pinfeas_new = 1e300;
  double dinfeas = 1e300;
  double complem = 1e300;
  double complem_p = 1e300;
  double dualgap = 1e300;
  double udotr = 1e300;
  double projerr = 1e300;
  double error[6];
  double totalresidual = 1e300, totalresidual_Jor = 1e300;
  double diff_fixp = 1e300, nub = 1e300;
  double *diff_fixp_vec = (double*)calloc(n,sizeof(double));

  double *primalConstraint = data->tmp_vault_nd[1];
  double *dualConstraint = data->tmp_vault_m[0];
  double *complemConstraint = data->tmp_vault_nd[2];

  double gmm = gmmp1+gmmp2;
  double barr_param_a, e, barr_param_fixed;
  double norm_f = cblas_dnrm2(m, f, 1);
  double norm_w = cblas_dnrm2(nd, w, 1);


  double *d_globalVelocity = (double*)calloc(m,sizeof(double));
  double *d_velocity = (double*)calloc(nd,sizeof(double));
  double *d_reaction = (double*)calloc(nd,sizeof(double));
  double *dw = (double*)calloc(nd,sizeof(double));
  double *u_plus_phiu = (double*)calloc(nd,sizeof(double));
  double *d_s = (double*)calloc(n,sizeof(double));
  double *s = (double*)calloc(n,sizeof(double));
  for (unsigned int  i = 0; i<n; i++) s[i] = cblas_dnrm2(2, velocity+i*d+1, 1);

  /* double *tmpsol = (double*)calloc(m+2*nd,sizeof(double)); // temporary solution */

  double *rhs = options->dWork;
  double * rhs_2 = (double*)calloc(m+2*nd, sizeof(double));
  double *gv_plus_dgv = data->tmp_vault_m[1];
  double *vr_jprod = data->tmp_vault_nd[3];
  double *v_plus_dv = data->tmp_vault_nd[4];
  double *r_plus_dr = data->tmp_vault_nd[5];
  double *vr_prod_sub_iden = data->tmp_vault_nd[6];
  double *dvdr_jprod = data->tmp_vault_nd[7];


  //double * r_p = (double*)calloc(nd,sizeof(double));                          // scaling vector p
  NumericsMatrix* r_Qp = NULL;                                                // matrix Qp
  //NumericsMatrix *QpH = NM_create(H->storageType, H->size0, H->size1);        // store the matrix Qp*H
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
  double r_alpha_primal, r_alpha_dual;
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
  double LS_norm_p = 0; // primal feasibility
  double LS_norm_d = 0; // dual feaqsibility
  double LS_norm_c = 0; // complementarity

  // Timer
  long clk_tck = CLOCKS_PER_SEC;
  clock_t t1, t2;
  double total_time = 0, max_time = 1.;

  // For TEST7
  double tols_TEST7_scalarP = 1e-3, tols_TEST7_prj = 1e-3, tols_TEST7_JordanP = 1e-8;  // Var used for checking stopping test satisfying different tols
  double projerr_TEST7 = 0;
  double norm_q_TEST7 = 0;
  int num_TEST7_scalarP = 30, num_TEST7_prj = 30, num_TEST7_JordanP = 8;
  size_t J_nz_captured = 0, J_nz_final = 0;

  FILE *file_TEST7_scalarP = NULL, *all_TEST7_scalarP = NULL, *file_TEST7_prj = NULL, *all_TEST7_prj = NULL, *file_TEST7_JordanP = NULL, *all_TEST7_JordanP = NULL;
  char name_file_TEST7_scalarP[30], name_file_TEST7_prj[30], name_file_TEST7_JordanP[300];


  int nB, nN, nR, nT;
  nB = nN = nR = nT = 0;


  NumericsMatrix *J = NULL;

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
      // numerics_printf_verbose(-1,"LS solution: 2x2 NT scaling with QpH\n");
      if (type == NORM_2)
        numerics_printf_verbose(-1,"LS solution: 2x2 NT scaling with QpH,\t stopping test: norm 2\n");
      else if (type == NORM_INF)
        numerics_printf_verbose(-1,"LS solution: 2x2 NT scaling with QpH,\t stopping test: norm inf\n");
      else
        numerics_printf_verbose(-1,"LS solution: 2x2 NT scaling with QpH\n");
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_QPH:
    {
      numerics_printf_verbose(-1,"LS solution: 1x1 NT scaling with QpH\n");
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_TEST6:
    {
      numerics_printf_verbose(-1,"TESTING: 3x3 no scaling, mu fixed, s fixed\n");
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QPH_TEST:
    {
      numerics_printf_verbose(-1,"TESTING: 2x2 NT scaling with QpH, mu fixed, s fixed\n");
      break;
    }
    default:
    {
      printf("ERROR\n");
    }
  }

  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_TEST6 ||
      options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QPH_TEST)
  {
    numerics_printf_verbose(-1, "| it  |  dualgap | pinfeas | dinfeas |  <u, r> |  |s-ub| | complem | prj err | barpram |  2*n*mu |  alpha  |  |dv|   |  |du|   |  |dr|   | ls prim | ls dual | ls comp |");
    numerics_printf_verbose(-1, "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------");
  }
  else
  {
    numerics_printf_verbose(-1, "| it  |  dualgap | pinfeas | dinfeas |  <u, r> | complem | prj err | barpram | alpha_p | alpha_d |  sigma  |  |dv|   |  |du|   |  |dr|   | ls prim | ls dual | ls comp |");
    numerics_printf_verbose(-1, "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------");
  }

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
  NumericsMatrix * eye_phiu = NULL;
  double * velocity_t = data->tmp_vault_nd[11];
  double * d_velocity_t = data->tmp_vault_nd[12];
  double * d_reaction_t = data->tmp_vault_nd[13];
  double * velocity_t_inv = data->tmp_vault_nd[14];
  double * Qp_velocity_t_inv = data->tmp_vault_nd[15];
  //double * tmp1 = data->tmp_vault_nd[16];
  FILE * iterates;
  FILE * matrixH;

  char *str = (char *) malloc(200);
  if (problem->name)
    {
      strcpy( str, problem->name );
    }
  else
    {
      strcpy( str, "foo_" );
    }
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

  char matlab_name[100];
  // sprintf(matlab_name, "iterates.m",strToken);
  // sprintf(matlab_name, "iterates__%.10s.m",strToken);


  /* writing data in a Matlab file */
  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
  {
    iterates = fopen("./mat/NOSCALEXT.m", "w");
    // iterates = fopen("./mat/QPHEXT.m", "w");
    // iterates = fopen(matlab_name, "w");
    printDataProbMatlabFile(M, f, H, w, d, n, m, problem->mu, iterates);
    // sprintf(matlab_name, "iterates_muF_sF_1.m");
    // iterates = fopen(matlab_name, "a+");
    // fprintf(iterates,"data(end+1).name = \"%s\";\n", strToken);
    // fprintf(iterates,"data(end).val = [\n");
  }

  ComputeErrorGlobalPtr computeError = NULL;
  computeError = (ComputeErrorGlobalPtr)&gfc3d_compute_error_convex;
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
  int load_starting_point = 0, save_sol_point = 0;
  FILE *sol_file = NULL;
  // barr_param = cblas_ddot(nd, reaction, 1, velocity, 1) / n ;
  while(iteration < max_iter)
  {
    // if ( options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING] == SICONOS_FRICTION_3D_IPM_IPARAM_FINISH_WITHOUT_SCALING_YES )
    //   if ( (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] != SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL) && (totalresidual <= 1e-10) && (fws==' ') )
    //   {
    //   	// To solve the problem very accurately, the algorithm switches to a direct solution of the linear system without scaling and without reduction //
    //   	options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] = SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL;
    //   	fws = '*';
    //   	// copy of the current solution into a temporary vector to evaluate the distance of this solution to the final one
    //   	/* cblas_dcopy(m, globalVelocity, 1, tmpsol, 1); */
    //   	/* cblas_dcopy(nd, velocity, 1, tmpsol+m, 1); */
    //   	/* cblas_dcopy(nd, reaction, 1, tmpsol+m+nd, 1); */
    //   }

    t1 = clock();
    /** Correction of w to take into account the dependence on the tangential velocity */
    // if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] == 1)
    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] == SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S_AT_EACH_ITE)
    {
      for(unsigned int i = 0; i < nd; ++ i)
      {
        if(i % d == 0)
        {
          w[i] = w_tilde[i] + cblas_dnrm2(2, velocity+i+1, 1);
        }
      }
    }
    else if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] == SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S_AFTER_SOLVING_LS)
    {
      primalResidual_type(velocity, H, globalVelocity, w, primalConstraint, &pinfeas, tol, type);
      dualResidual_type(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol, type);
      udotr = xdoty_type(n, nd, velocity, reaction, type);
      complem = complemResidualNorm_type(velocity, reaction, nd, n, NORM_2_INF);

      if (fmax(pinfeas, fmax(dinfeas, complem)) < tol)
      {
        numerics_printf_verbose(-1, "| %3i%c| %8.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |",
                              iteration, fws, dualgap, pinfeas, dinfeas, udotr, complem, projerr, barr_param);
        // int nB, nN, nR, nT;
        // classify_BNRT(velocity, reaction, nd, n, &nB, &nN, &nR, &nT);
        // printf("B/N/R/T = %d %d %d %d\n", nB, nN, nR, nT);

        if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
        {
          double worktmp[3];
          printIteresProbMatlabFile(iteration, globalVelocity, velocity, reaction, d_globalVelocity, d_velocity, d_reaction, d, n, m, iterates);
          fprintf(iterates, "pinfeas(%3i) = %.20e;\n", iteration+1, pinfeas);
          fprintf(iterates, "dinfeas(%3i) = %.20e;\n", iteration+1, dinfeas);
          fprintf(iterates, "udotr(%3i) = %.20e;\n", iteration+1, udotr);
          fprintf(iterates, "complem(%3i) = %.20e;\n", iteration+1, complem);
          fprintf(iterates, "LSp(%3i) = %.20e;\n", iteration+1, LS_norm_p);
          fprintf(iterates, "LSd(%3i) = %.20e;\n", iteration+1, LS_norm_d);
          fprintf(iterates, "LSc(%3i) = %.20e;\n", iteration+1, LS_norm_c);
          fprintf(iterates, "alpha(%3i) = %.20e;\n", iteration+1, alpha_primal);
          fprintf(iterates,"projerr(%3i,:) = [",iteration+1);
          for(int i = 0; i < n; i++)
          {
            fprintf(iterates, "%.20e, ", projectionError(velocity+i*d, reaction+i*d, 1, tol));
          }
          fprintf(iterates,"];\n");
          fprintf(iterates,"proj(%3i,:) = [",iteration+1);
          for(int i = 0; i < n; i++)
          {
            worktmp[0] = reaction[3*i] -  velocity[3*i] ;
            worktmp[1] = reaction[3*i+1] -  velocity[3*i+1] ;
            worktmp[2] = reaction[3*i+2] -  velocity[3*i+2] ;
            projectionOnCone(worktmp, 1.0);
            fprintf(iterates, "%.20e, %.20e, %.20e, ", worktmp[0], worktmp[1], worktmp[2]);
          }
          fprintf(iterates,"];\n");
          fprintf(iterates,"nrmU(%3i,:) = [",iteration+1);
          for(int i = 0; i < n; i++)
          {
            fprintf(iterates, "%.20e, ", cblas_dnrm2(d, velocity+i*d, 1));
          }
          fprintf(iterates,"];\n");
          fprintf(iterates,"nrmR(%3i,:) = [",iteration+1);
          for(int i = 0; i < n; i++)
          {
            fprintf(iterates, "%.20e, ", cblas_dnrm2(d, reaction+i*d, 1));
          }
          fprintf(iterates,"];\n");
          fprintf(iterates,"nrmDU(%3i,:) = [",iteration+1);
          for(int i = 0; i < n; i++)
          {
            fprintf(iterates, "%.20e, ", cblas_dnrm2(d, d_velocity+i*d, 1));
          }
          fprintf(iterates,"];\n");
          fprintf(iterates,"nrmDR(%3i,:) = [",iteration+1);
          for(int i = 0; i < n; i++)
          {
            fprintf(iterates, "%.20e, ", cblas_dnrm2(d, d_reaction+i*d, 1));
          }
          fprintf(iterates,"];\n");
          fprintf(iterates,"eig1U(%3i,:) = [",iteration+1);
          for(int i = 0; i < n; i++)
          {
            fprintf(iterates, "%.20e, ", velocity[i*d] + cblas_dnrm2(2, velocity+i*d+1, 1));
          }
          fprintf(iterates,"];\n");
          fprintf(iterates,"eig1R(%3i,:) = [",iteration+1);
          for(int i = 0; i < n; i++)
          {
            fprintf(iterates, "%.20e, ", reaction[i*d] + cblas_dnrm2(2, reaction+i*d+1, 1));
          }
          fprintf(iterates,"];\n");
          fprintf(iterates,"eig2U(%3i,:) = [",iteration+1);
          for(int i = 0; i < n; i++)
          {
            fprintf(iterates, "%.20e, ", velocity[i*d] - cblas_dnrm2(2, velocity+i*d+1, 1));
          }
          fprintf(iterates,"];\n");
          fprintf(iterates,"eig2R(%3i,:) = [",iteration+1);
          for(int i = 0; i < n; i++)
          {
            fprintf(iterates, "%.20e, ", reaction[i*d] - cblas_dnrm2(2, reaction+i*d+1, 1));
          }
          fprintf(iterates,"];\n");
        }
        iteration++;



        for(unsigned int i = 0; i < nd; ++ i)
        {
          if(i % d == 0)
          {
            w[i] = w_tilde[i] + cblas_dnrm2(2, velocity+i+1, 1);
          }
        }
        printf("\n");
        numerics_printf_verbose(-1, "| it  |  dualgap | pinfeas | dinfeas |  <u, r> | complem | prj err | barpram | alpha_p | alpha_d |  sigma  |  |dv|   |  |du|   |  |dr|   | ls prim | ls dual | ls comp |");
        numerics_printf_verbose(-1, "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------");

      }
    }
    // else if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] == SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S_AFTER_SOLVING_LS_AND_CONVERGE_TO_MU)
    // {
    //   if (iteration == 0)
    //   {
    //     for(unsigned int i = 0; i < nd; ++ i)
    //     {
    //       if(i % d == 0)
    //       {
    //         w[i] = w_tilde[i] + cblas_dnrm2(2, velocity+i+1, 1);
    //       }
    //     }
    //     barr_param_fixed = 1.;
    //   }
    // }


    // if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] == 1)
    // {
      // if (iteration == 0)
      // {
      //   if (load_starting_point)
      //   {
      //     sol_file = fopen("sol_data.res", "r");
      //     if (!sol_file) printf("\n\ngfc3d_ipm: Solution data file is not available!!! \n\n");
      //     else
      //     {
      //       // load v
      //       for (int i=0; i < m; i++)
      //       {
      //         fscanf(sol_file, "%lf ", globalVelocity+i);
      //       }
      //       fscanf(sol_file, "\n");

      //       // load u
      //       for (int i=0; i < nd; i++)
      //       {
      //         fscanf(sol_file, "%lf ", velocity+i);
      //       }
      //       fscanf(sol_file, "\n");

      //       // load r
      //       for (int i=0; i < nd; i++)
      //       {
      //         fscanf(sol_file, "%lf ", reaction+i);
      //       }
      //       fscanf(sol_file, "\n");

      //       printf("\ngfc3d_ipm: Starting point is successfully loaded.\n");
      //       for(unsigned int i=0; i<nd; i+=d)
      //       {
      //         printf("\n%4i-%4i: u =", i, i+d-1);
      //         for(unsigned int j=0; j<d; j++) printf(" %.2e,", velocity[i+j]);
      //         printf("\t u0 - |ub| = %.2e \t r =", velocity[i] - cblas_dnrm2(2, velocity+i+1, 1));
      //         for(unsigned int j=0; j<d; j++) printf(" %.2e,", reaction[i+j]);
      //         printf("\t r0 - |rb| = %.2e", reaction[i] - cblas_dnrm2(2, reaction+i+1, 1));
      //       }
      //       printf("\n\n");
      //       load_starting_point = 0;
      //     }
      //     fclose(sol_file);
      //   }


      //   // barr_param = cblas_ddot(nd, reaction, 1, velocity, 1) / n ;
      //   barr_param = 1e-5;
      //   for(unsigned int i = 0; i < nd; ++ i)
      //   {
      //     if(i % d == 0)
      //     {
      //       phiu_val = cblas_dnrm2(2, velocity+i+1, 1);
      //       w[i] = w_tilde[i] + phiu_val;
      //       phiu[i/d] = phiu_val;
      //     }
      //   }
      // }
      // for(unsigned int i = 0; i < nd; ++ i)
      // {
      //   if(i % d == 0)
      //   {
      //     phiu_val = cblas_dnrm2(2, velocity+i+1, 1);
      //     w[i] = w_tilde[i] + phiu_val;
      //     s[i/d] = phiu_val;
      //     phiu[i/d] = phiu_val;
      //   }
      // }
      // primalResidual(velocity, H, globalVelocity, w, primalConstraint, &pinfeas, tol); // primalConstraint = u - H*v -w
      // dualResidual(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol);  // dualConstraint =  M*v - H'*r - f
      // dualgap = dualGap(M, f, w, globalVelocity, reaction, nd, m);
      // JA_prod(velocity, reaction, nd, n, complemConstraint);
      // for (unsigned int k = 0; k < nd; complemConstraint[k] -= 2 * barr_param, k += d);
      // complem = cblas_dnrm2(nd, complemConstraint, 1);
      // diff_fixp = 0.;
      // for (unsigned int i = 0; i<nd; i+=d)
      // {
      //   nub = cblas_dnrm2(2, velocity+i+1, 1);
      //   diff_fixp_vec[i/d] = s[i/d]-nub;
      // }
      // diff_fixp = cblas_dnrm2(n, diff_fixp_vec, 1);
      // totalresidual_mu = fmax(fmax(pinfeas, dinfeas),complem);


      // if (totalresidual_mu <= tol)
      // {
      //   numerics_printf_verbose(-1, "| %3i%c| %8.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |",
      //                         iteration, fws, dualgap, pinfeas, dinfeas, udotr, diff_fixp, complem, projerr, barr_param);
      //   // for(unsigned int i=0; i<nd; i+=d)
      //   // {
      //   //   printf("\n%4i-%4i: u =", i, i+d-1);
      //   //   for(unsigned int j=0; j<d; j++) printf(" %.2e,", velocity[i+j]);
      //   //   printf("\t u0 - |ub| = %.2e \t r =", velocity[i] - cblas_dnrm2(2, velocity+i+1, 1));
      //   //   for(unsigned int j=0; j<d; j++) printf(" %.2e,", reaction[i+j]);
      //   //   printf("\t r0 - |rb| = %.2e", reaction[i] - cblas_dnrm2(2, reaction+i+1, 1));
      //   // }

      //   if (diff_fixp <= tol)
      //   {
      //     if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
      //     {
      //       fprintf(iterates,"%d %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e;\n",
      //               iteration, pinfeas, dinfeas, udotr, complem, projerr, barr_param, alpha_primal,
      //               fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
      //               fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
      //               fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]),
      //               fabs(d_s[cblas_idamax(n, d_s, 1)]),
      //               LS_norm_p, LS_norm_d, LS_norm_c);
      //     }
      //     hasNotConverged = 0;
      //     break;
      //   }

      //   // Store solution into file
      //   if (save_sol_point)
      //   {
      //     sol_file = fopen("sol_data.res", "w");
      //     // store v
      //     for (int i=0; i < m; i++)
      //     {
      //       fprintf(sol_file, "%8.20e ", globalVelocity[i]);
      //     }
      //     fprintf(sol_file, "\n");

      //     // store u
      //     for (int i=0; i < nd; i++)
      //     {
      //       fprintf(sol_file, "%8.20e ", velocity[i]);
      //     }
      //     fprintf(sol_file, "\n");

      //     // store r
      //     for (int i=0; i < nd; i++)
      //     {
      //       fprintf(sol_file, "%8.20e ", reaction[i]);
      //     }
      //     fprintf(sol_file, "\n");
      //     fclose(sol_file);
      //     printf("\n\n Saved successfully solution points. \n\n");
      //     save_sol_point = 0;
      //   }


      //   if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
      //   {
      //     fprintf(iterates,"%d %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e;\n",
      //             iteration, pinfeas, dinfeas, udotr, complem, projerr, barr_param, alpha_primal,
      //             fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
      //             fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
      //             fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]),
      //             fabs(d_s[cblas_idamax(n, d_s, 1)]),
      //             LS_norm_p, LS_norm_d, LS_norm_c);
      //   }

      //   double unitur;
      //   for (int i = 0; i < n; i++)
      //   {
      //     unitur = cblas_ddot(3, velocity+3*i, 1, reaction+3*i, 1);
      //     // unitur = cblas_ddot(3, u_plus_phiu+3*i, 1, reaction+3*i, 1);
      //     if (unitur<0)
      //       printf("UR NEGATIF %9.2e\n", unitur);
      //   }
      //   // hasNotConverged = 0;
      //   // break;
      //   printf("\n");
      //   numerics_printf_verbose(-1, "| it  |  dualgap | pinfeas | dinfeas |  <u, r> |  |s-ub| | complem | prj err | barpram |  2*n*mu |  alpha  |  |dv|   |  |du|   |  |dr|   | ls prim | ls dual | ls comp |");

        // cblas_dcopy(nd, w, 1, dw, 1); /* dw is used to compute an extrapolation step for the velocity after an update of w */
        // for(unsigned int i = 0; i < nd; ++ i)
        // {
        //   if(i % d == 0)
        //   {
        //     phiu_val = cblas_dnrm2(2, velocity+i+1, 1);
        //     w[i] = w_tilde[i] + phiu_val;
        //     // s[i/d] = phiu_val;
        //     // phiu[i/d] = phiu_val;
        //   }
        // }
      // }






      // /* computation of the extrapolation step for u */
      // cblas_daxpy(nd, -1.0, w, 1, dw, 1);
      // cblas_dscal(nd, -1.0, dw, 1);                            /* dw = w_new - w_old */
      // alpha_primal = getStepLength(velocity, dw, nd, n, 0.9);  /* u_new = u_old + alpha*dw */
      // cblas_daxpy(nd, alpha_primal, dw, 1, velocity, 1);

      // double diff_pinfeas = 0.0;
      // for(unsigned int i = 0; i < nd; ++ i)
      // {
      //   if(i % d == 0)
      //   {
      //     diff_pinfeas += (1+alpha_primal)*phiu[i] - alpha_primal*sqrt(velocity[i+1]*velocity[i+1]+velocity[i+2]*velocity[i+2]);
      //   }
      // }

      //printf("alpha_u = %e,\t", alpha_primal); //printf("diff_pinfeas = %e\n\n", diff_pinfeas);
    // }





    /* Computation of the values of
     - primal residual: u - H*v - w
     - dual residual: M*v - f - H'*r
     - duality gap: u'*r
     - true duality gap: (value_of_primal - value_of_dual)
     - complementarity: u o r
     - projection error: r - proj(r-u)
    */







    // Stopping test using norm type
    primalResidual_type(velocity, H, globalVelocity, w, primalConstraint, &pinfeas, tol, type);
    dualResidual_type(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol, type);
    udotr = xdoty_type(n, nd, velocity, reaction, type);
    dualgap = dualGap(M, f, w, globalVelocity, reaction, nd, m);
    // complem = complemResidualNorm(velocity, reaction, nd, n);
    complem = complemResidualNorm_type(velocity, reaction, nd, n, NORM_2_INF);



    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_TEST6 ||
        options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QPH_TEST)
    {
      JA_prod(velocity, reaction, nd, n, complemConstraint);
      for (unsigned int k = 0; k < nd; complemConstraint[k] -= 2 * barr_param, k += d);
      complem = cblas_dnrm2(nd, complemConstraint, 1);

      diff_fixp = 0.;
      for (unsigned int i = 0; i<nd; i+=d)
      {
        nub = cblas_dnrm2(2, velocity+i+1, 1);
        diff_fixp_vec[i/d] = s[i/d]-nub;
      }
      diff_fixp = cblas_dnrm2(n, diff_fixp_vec, 1);
    }
    else
    {
      // barr_param = cblas_ddot(nd, reaction, 1, velocity, 1) / n / 3.;
      // barr_param = cblas_ddot(nd, reaction, 1, velocity, 1) / n;
      // barr_param = cblas_ddot(nd, reaction, 1, velocity, 1) / (3*n);
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] > 0)
      {
        barr_param = complemResidualNorm(velocity, reaction, nd, n)/n;
        // barr_param = cblas_ddot(nd, reaction, 1, velocity, 1) / n;
      }
      else
        barr_param = cblas_ddot(nd, reaction, 1, velocity, 1) / n / 3.;
    }


    t2 = clock();
    total_time += (double)(t2-t1)/(double)clk_tck;


    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] > 0)
    {
      projerr = projectionError_norm_infinity_conic(velocity, reaction, n);
      projerr_TEST7 = projerr;
    }
    else
      projerr = projectionError(velocity, reaction, n, tol);



    // projerr = projectionError(u_plus_phiu, reaction, n, tol);
    // projerr = projectionError(velocity, reaction, n, tol);
    // projerr = projectionError_norm_infinity_conic(velocity, reaction, n);
    // projerr = projectionError_relative_norm_infinity_conic(velocity, reaction, n);


    setErrorArray(error, pinfeas, dinfeas, udotr, dualgap, complem, projerr);

    // check exit condition
    totalresidual = fmax(fmax(error[0], error[1]),error[2]);
    totalresidual_Jor = fmax(fmax(pinfeas, dinfeas),complem);

    // if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
    //   printIteresProbMatlabFile(iteration, globalVelocity, velocity, reaction, d_globalVelocity, d_velocity, d_reaction, d, n, m, iterates);



    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] > 0)
    {

      // Criteria complementarity
      int repeat = 0;
      if (total_time <= max_time && alpha_primal > 1e-13)
      {
        // // Check Stopping test for different tols, vary from 1e-3 to 1e-10
        // while ( totalresidual <= tols_TEST7_scalarP && num_TEST7_scalarP < 11 )
        // {
        //   repeat++;
        //   if (repeat == 1)
        //   {
        //     classify_BNRT(velocity, reaction, nd, n, &nB, &nN, &nR, &nT);
        //   }

        //   sprintf(name_file_TEST7_scalarP, "ipm-scalarP-%02d.pp", num_TEST7_scalarP);
        //   file_TEST7_scalarP = fopen(name_file_TEST7_scalarP, "a+");
        //   fprintf(file_TEST7_scalarP, "sumry: 0  %.2e   %5i %5i %5i %4i %4i %4i %4i  %.6f   %s\n",
        //           totalresidual, iteration, n, m, nB, nN, nR, n-nB-nN-nR,
        //           total_time, strToken);
        //   fclose(file_TEST7_scalarP);

        //   all_TEST7_scalarP = fopen("ipm-scalarP.pp", "a+");
        //   if (num_TEST7_scalarP == 10)
        //   {
        //     fprintf(all_TEST7_scalarP, "sumry: 0  %.2e  %.2e  %.2e   %5i %5i %5i %4i %4i %4i %4i  %.6f   %s\n",
        //           totalresidual, complem, projerr, iteration, n, m, nB, nN, nR, n-nB-nN-nR,
        //           total_time, strToken);
        //   }
        //   else
        //   {
        //     fprintf(all_TEST7_scalarP, "sumry: 0  %.2e  %.2e  %.2e   %5i %5i %5i %4i %4i %4i %4i  %.6f\n",
        //           totalresidual, complem, projerr, iteration, n, m, nB, nN, nR, n-nB-nN-nR,
        //           total_time);
        //   }
        //   fclose(all_TEST7_scalarP);

        //   tols_TEST7_scalarP /= 10.;
        //   num_TEST7_scalarP++;
        // }

        // // Criteria projection
        // repeat = 0;
        // while ( projerr_TEST7 <= tols_TEST7_prj && num_TEST7_prj < 11 )
        // {
        //   repeat++;
        //   if (repeat == 1)
        //   {
        //     classify_BNRT(velocity, reaction, nd, n, &nB, &nN, &nR, &nT);
        //   }

        //   sprintf(name_file_TEST7_prj, "ipm-prjerr-%02d.pp", num_TEST7_prj);
        //   file_TEST7_prj = fopen(name_file_TEST7_prj, "a+");
        //   fprintf(file_TEST7_prj, "sumry: 0  %.2e   %5i %5i %5i %4i %4i %4i %4i  %.6f   %s\n",
        //           projerr_TEST7, iteration, n, m, nB, nN, nR, n-nB-nN-nR,
        //           total_time, strToken);
        //   fclose(file_TEST7_prj);

        //   all_TEST7_prj = fopen("ipm-prjerr.pp", "a+");
        //   if (num_TEST7_prj == 10)
        //   {
        //     fprintf(all_TEST7_prj, "sumry: 0  %.2e  %.2e  %.2e   %5i %5i %5i %4i %4i %4i %4i  %.6f   %s\n",
        //           totalresidual, complem, projerr, iteration, n, m, nB, nN, nR, n-nB-nN-nR,
        //           total_time, strToken);
        //   }
        //   else
        //   {
        //     fprintf(all_TEST7_prj, "sumry: 0  %.2e  %.2e  %.2e   %5i %5i %5i %4i %4i %4i %4i  %.6f\n",
        //           totalresidual, complem, projerr, iteration, n, m, nB, nN, nR, n-nB-nN-nR,
        //           total_time);
        //   }
        //   fclose(all_TEST7_prj);

        //   tols_TEST7_prj /= 10.;
        //   num_TEST7_prj++;
        // }

        // totalresidual with Jordan product
        repeat = 0;
        while ( totalresidual_Jor <= tols_TEST7_JordanP && num_TEST7_JordanP < 9 )
        {
          repeat++;
          if (repeat == 1)
          {
            classify_BNRT(velocity, reaction, nd, n, &nB, &nN, &nR, &nT);
          }

          sprintf(name_file_TEST7_JordanP, "./res_20240706/update/external/NOSCAL/ipm-JordanP-%02d.pp", num_TEST7_JordanP);
          file_TEST7_JordanP = fopen(name_file_TEST7_JordanP, "a+");
          fprintf(file_TEST7_JordanP, "sumry: 0  %.2e   %5i %5i %5i %4i %4i %4i %4i  %.6f   %s\n",
                  totalresidual_Jor, iteration, n, m, nB, nN, nR, n-nB-nN-nR,
                  total_time, strToken);
          fclose(file_TEST7_JordanP);

          // all_TEST7_JordanP = fopen("ipm-JordanP.pp", "a+");
          // if (num_TEST7_JordanP == 10)
          // {
          //   fprintf(all_TEST7_JordanP, "sumry: 0  %.2e  %.2e  %.2e   %5i %5i %5i %4i %4i %4i %4i  %.6f   %s\n",
          //         totalresidual, complem, projerr, iteration, n, m, nB, nN, nR, n-nB-nN-nR,
          //         total_time, strToken);
          // }
          // else
          // {
          //   fprintf(all_TEST7_JordanP, "sumry: 0  %.2e  %.2e  %.2e   %5i %5i %5i %4i %4i %4i %4i  %.6f\n",
          //         totalresidual, complem, projerr, iteration, n, m, nB, nN, nR, n-nB-nN-nR,
          //         total_time);
          // }
          // fclose(all_TEST7_JordanP);

          tols_TEST7_JordanP /= 10.;
          num_TEST7_JordanP++;
        }
      } // end if (total_time <= max_time && alpha_primal > 1e-13)
    } // end if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] > 0)


    if (total_time > max_time)
    {
      hasNotConverged = 5;
      numerics_printf_verbose(-1, "| %3i%c| %8.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |",
                              iteration, fws, dualgap, pinfeas, dinfeas, udotr, complem, projerr, barr_param);
      printf("\nREQUEST TIMEOUT\n");
      double unitur;
      for (int i = 0; i < n; i++)
      {
        unitur = cblas_ddot(3, velocity+3*i, 1, reaction+3*i, 1);
        // unitur = cblas_ddot(3, u_plus_phiu+3*i, 1, reaction+3*i, 1);
        if (unitur<0)
          printf("UR NEGATIF %9.2e\n", unitur);
      }
      break;
    }

    // if (alpha_primal < 1e-10)
    // {
    //   printf("\nSMALL Step-length\n");
    //   break;
    // }



    double worktmp[3];
    // if ( totalresidual <= tol )
    if ( totalresidual_Jor <= tol )
    // if ( projerr <= tol )
    // if ( fmax(totalresidual,projerr) <= tol )
    // if ( fmax(fmax(totalresidual,projerr),complem) <= tol )
    {
      numerics_printf_verbose(-1, "| %3i%c| %8.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |",
                              iteration, fws, dualgap, pinfeas, dinfeas, udotr, complem, projerr, barr_param);

      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
      {
        printIteresProbMatlabFile(iteration, globalVelocity, velocity, reaction, d_globalVelocity, d_velocity, d_reaction, d, n, m, iterates);
        fprintf(iterates, "pinfeas(%3i) = %.20e;\n", iteration+1, pinfeas);
        fprintf(iterates, "dinfeas(%3i) = %.20e;\n", iteration+1, dinfeas);
        fprintf(iterates, "udotr(%3i) = %.20e;\n", iteration+1, udotr);
        fprintf(iterates, "complem(%3i) = %.20e;\n", iteration+1, complem);
        fprintf(iterates, "LSp(%3i) = %.20e;\n", iteration+1, LS_norm_p);
        fprintf(iterates, "LSd(%3i) = %.20e;\n", iteration+1, LS_norm_d);
        fprintf(iterates, "LSc(%3i) = %.20e;\n", iteration+1, LS_norm_c);
        fprintf(iterates, "alpha(%3i) = %.20e;\n", iteration+1, alpha_primal);
        fprintf(iterates,"projerr(%3i,:) = [",iteration+1);
        for(int i = 0; i < n; i++)
        {
          fprintf(iterates, "%.20e, ", projectionError(velocity+i*d, reaction+i*d, 1, tol));
        }
        fprintf(iterates,"];\n");
        fprintf(iterates,"proj(%3i,:) = [",iteration+1);
        for(int i = 0; i < n; i++)
        {
          worktmp[0] = reaction[3*i] -  velocity[3*i] ;
          worktmp[1] = reaction[3*i+1] -  velocity[3*i+1] ;
          worktmp[2] = reaction[3*i+2] -  velocity[3*i+2] ;
          projectionOnCone(worktmp, 1.0);
          fprintf(iterates, "%.20e, %.20e, %.20e, ", worktmp[0], worktmp[1], worktmp[2]);
        }
        fprintf(iterates,"];\n");
        fprintf(iterates,"nrmU(%3i,:) = [",iteration+1);
        for(int i = 0; i < n; i++)
        {
          fprintf(iterates, "%.20e, ", cblas_dnrm2(d, velocity+i*d, 1));
        }
        fprintf(iterates,"];\n");
        fprintf(iterates,"nrmR(%3i,:) = [",iteration+1);
        for(int i = 0; i < n; i++)
        {
          fprintf(iterates, "%.20e, ", cblas_dnrm2(d, reaction+i*d, 1));
        }
        fprintf(iterates,"];\n");
        fprintf(iterates,"nrmDU(%3i,:) = [",iteration+1);
        for(int i = 0; i < n; i++)
        {
          fprintf(iterates, "%.20e, ", cblas_dnrm2(d, d_velocity+i*d, 1));
        }
        fprintf(iterates,"];\n");
        fprintf(iterates,"nrmDR(%3i,:) = [",iteration+1);
        for(int i = 0; i < n; i++)
        {
          fprintf(iterates, "%.20e, ", cblas_dnrm2(d, d_reaction+i*d, 1));
        }
        fprintf(iterates,"];\n");
        fprintf(iterates,"eig1U(%3i,:) = [",iteration+1);
        for(int i = 0; i < n; i++)
        {
          fprintf(iterates, "%.20e, ", velocity[i*d] + cblas_dnrm2(2, velocity+i*d+1, 1));
        }
        fprintf(iterates,"];\n");
        fprintf(iterates,"eig1R(%3i,:) = [",iteration+1);
        for(int i = 0; i < n; i++)
        {
          fprintf(iterates, "%.20e, ", reaction[i*d] +  cblas_dnrm2(2, reaction+i*d+1, 1));
        }
        fprintf(iterates,"];\n");
        fprintf(iterates,"eig2U(%3i,:) = [",iteration+1);
        for(int i = 0; i < n; i++)
        {
          fprintf(iterates, "%.20e, ", velocity[i*d] - cblas_dnrm2(2, velocity+i*d+1, 1));
        }
        fprintf(iterates,"];\n");
        fprintf(iterates,"eig2R(%3i,:) = [",iteration+1);
        for(int i = 0; i < n; i++)
        {
          fprintf(iterates, "%.20e, ", reaction[i*d] - cblas_dnrm2(2, reaction+i*d+1, 1));
        }
        fprintf(iterates,"];\n");
      }

      double unitur;
      for (int i = 0; i < n; i++)
      {
      	unitur = cblas_ddot(3, velocity+3*i, 1, reaction+3*i, 1);
        // unitur = cblas_ddot(3, u_plus_phiu+3*i, 1, reaction+3*i, 1);
      	if (unitur<0)
      	  printf("UR NEGATIF %9.2e\n", unitur);
      }

      hasNotConverged = 0;
      if (Qp) Qp = NM_free(Qp);
      if (Qpinv) Qpinv = NM_free(Qpinv);
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
    /* Build the Jacobian matrix corresponding to a reduced linear system. The variable du (velocity) is eliminated.
     *
     * In this case, NT scaling is always performed, but two posible ways can be done to solve the linear system.
     *
     * In a first approach, the reduced matrix is of the following form:
     *
     *         m       nd
     *      |  -M      H^T    | m
     * JR = |                 |
     *      |   H     Qp^{-2} | nd
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
     *
     * Case of a non-reduced system. Building the right-hand side related to the first system

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

    /* Building the rhs for the first linear system in case of a reduced linear system  */
    /* In the case where Qp2 is inside the matrix, the rhs is the of the form           */
    /*         [ M*v - f - H'*r                ]  m         dualConstraint              */
    /* r_rhs = [                               ]                                        */
    /*         [ u - H*v - w - Qp_inv*r_check  ]  nd                                    */
    /*                                                                                  */
    /* In the case where QpH is inside the matrix, the rhs is the of the form           */
    /*        [ M*v - f - H'*r ]  m         dualConstraint                             */
    /* r_rhs = [                ]                                                       */
    /*         [ -Qp*(H*v + w)  ]  nd                                                   */

    int jacobian_is_nan = 0;
    switch ( options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] )
    {
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL:
    {
      t1 = clock();
      // First linear linear system
      J = NM_create(NM_SPARSE, m + 2*nd, m + 2*nd);
      J_nzmax = (d * d) * (m / d) + H_nzmax + (3*d-2)*n + 3*d*n + H_nzmax + nd;
      NM_triplet_alloc(J, J_nzmax);
      J->matrix2->origin = NSM_TRIPLET;

      NumericsMatrix * arrow_r = Arrow_repr(reaction, nd, n);
      NumericsMatrix * arrow_u = Arrow_repr(velocity, nd, n) ;

      NM_insert(J, M, 0, 0);
      NM_insert(J, minus_H, m + nd, 0);
      NM_insert(J, arrow_r, m, m);
      NM_insert(J, eye_nd, m + nd, m);
      NM_insert(J, minus_Ht, 0, m + nd);
      NM_insert(J, arrow_u, m, m + nd);

      /* regularization */
       // NM_insert(J, NM_scalar(nd, -fmin(1e-7,barr_param)), m + nd, m + nd);

      NM_free(arrow_r);
      NM_free(arrow_u);


      cblas_dcopy(m + 2 * nd, rhs, 1, rhs_2, 1);
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] == SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES) {
        	double *rhs_save = (double *)calloc(m + 2 * nd, sizeof(double));
        	cblas_dcopy(m + 2 * nd, rhs, 1, rhs_save, 1);
        	//NM_LU_refine(J, rhs, rhs_save, 1, tol, 10, 0);
        	double residu;
        	NM_LU_refine(J, rhs,  1e-14, 10,  &residu);
        	free(rhs_save);
      }
      else
        NM_LU_solve(J, rhs, 1);

      // // Create Arw(u) + |ub|I
      // NumericsMatrix * arrow_u = NM_create(NM_SPARSE, nd, nd);
      // size_t arrow_u_nzmax = (d * 3 - 2) * n;
      // NM_triplet_alloc(arrow_u, arrow_u_nzmax);
      // NM_fill(arrow_u, NM_SPARSE, nd, nd, arrow_u->matrix2);

      // // Create Arw(r) + r*(0 ub/|ub|)
      // NumericsMatrix * arrow_r = NM_create(NM_SPARSE, nd, nd);
      // size_t arrow_r_nzmax = d * 3 * n;
      // NM_triplet_alloc(arrow_r, arrow_r_nzmax);
      // NM_fill(arrow_r, NM_SPARSE, nd, nd, arrow_r->matrix2);

      // /* Arrow matrix filling */
      // size_t pos; double ub;
      // for(size_t i = 0; i < n; ++i)
      // {
      //   pos = i * d;
      //   ub = sqrt(velocity[pos+1]*velocity[pos+1]+velocity[pos+2]*velocity[pos+2]);

      //   NM_entry(arrow_u, pos, pos, velocity[pos]+ub);  // Arw(u) + |ub|I
      //   NM_entry(arrow_r, pos, pos, reaction[pos]);     // Arw(r) + r*(0 ub/|ub|)

      //   for(size_t j = 1; j < d; ++j)
      //   {
      //     NM_entry(arrow_u, pos, pos + j, velocity[pos + j]);
      //     NM_entry(arrow_u, pos + j, pos, velocity[pos + j]);
      //     NM_entry(arrow_u, pos + j, pos + j, velocity[pos]+ub);

      //     NM_entry(arrow_r, pos, pos + j, reaction[pos + j] + reaction[pos]*velocity[pos + j]/ub);
      //     NM_entry(arrow_r, pos + j, pos, reaction[pos + j]);
      //     NM_entry(arrow_r, pos + j, pos + j, reaction[pos] + reaction[pos + j]*velocity[pos + j]/ub);
      //   }
      //   NM_entry(arrow_r, pos + 1, pos + 2, reaction[pos + 1]*velocity[pos + 2]/ub);
      //   NM_entry(arrow_r, pos + 2, pos + 1, reaction[pos + 2]*velocity[pos + 1]/ub);
      // }


      // NM_insert(J, M, 0, 0);
      // NM_insert(J, minus_H, m + nd, 0);
      // NM_insert(J, arrow_r, m, m);
      // NM_insert(J, eye_nd, m + nd, m);
      // NM_insert(J, minus_Ht, 0, m + nd);
      // NM_insert(J, arrow_u, m, m + nd);

      // /* regularization */
      // /* NM_insert(J, NM_scalar(nd, -barr_param), m + nd, m + nd); */

      // NM_free(arrow_r);
      // NM_free(arrow_u);




      jacobian_is_nan = NM_isnan(J);
      if (jacobian_is_nan)
      {
        numerics_printf_verbose(0, "The Jacobian matrix J contains NaN");
        break;
      }

      JA_prod(velocity, reaction, nd, n, complemConstraint);
      // JA_prod(u_plus_phiu, reaction, nd, n, complemConstraint);

      // if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] == SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S_AFTER_SOLVING_LS_AND_CONVERGE_TO_MU)
      // {
      //   for (unsigned int i = 0; i<nd; i+=d)
      //   {
      //     complemConstraint[i] -= barr_param_fixed;
      //   }
      // }

      cblas_dcopy(m, dualConstraint, 1, rhs, 1);
      cblas_dcopy(nd, complemConstraint, 1, rhs+m, 1);
      cblas_dcopy(nd, primalConstraint, 1, rhs+m+nd, 1);

      cblas_dscal(m + nd + nd, -1.0, rhs, 1);

      cblas_dcopy(m+2*nd, rhs, 1, rhs_2, 1);

      NM_LU_solve(J, rhs, 1);

      cblas_dcopy(m, rhs, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, rhs+m, 1, d_velocity, 1);
      cblas_dcopy(nd, rhs+m+nd, 1, d_reaction, 1);

      t2 = clock();
      total_time += (double)(t2-t1)/(double)clk_tck;
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QP2:
    {
      // First linear linear system
      J = NM_create(NM_SPARSE, m + 2*nd, m + 2*nd);
      J_nzmax = (d * d) * (m / d) + H_nzmax + 2 * (d * 3 - 2) * n + H_nzmax + nd;
      NM_triplet_alloc(J, J_nzmax);
      J->matrix2->origin = NSM_TRIPLET;

      Nesterov_Todd_vector(2, velocity, reaction, nd, n, p2);
      Qp2 = QRmat(p2, nd, n);

      NM_insert(J, M, 0, 0);
      NM_insert(J, minus_H, m + nd, 0);
      NM_insert(J, Qp2, m, m);
      NM_insert(J, eye_nd, m + nd, m);
      NM_insert(J, minus_Ht, 0, m + nd);
      NM_insert(J, eye_nd, m, m + nd);

      if(Qp2) Qp2 = NM_free(Qp2);

      jacobian_is_nan = NM_isnan(J);
      if (jacobian_is_nan)
  {
    numerics_printf_verbose(0, "The Jacobian matrix J contains NaN");
    break;
  }

      cblas_dcopy(m, dualConstraint, 1, rhs, 1);
      cblas_dcopy(nd, reaction, 1, rhs+m, 1);
      cblas_dcopy(nd, primalConstraint, 1, rhs+m+nd, 1);

      cblas_dscal(m + nd + nd, -1.0, rhs, 1);

      cblas_dcopy(m+2*nd, rhs, 1, rhs_2, 1);

      NM_LDLT_solve(J, rhs, 1);

      cblas_dcopy(m, rhs, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, rhs+m, 1, d_velocity, 1);
      cblas_dcopy(nd, rhs+m+nd, 1, d_reaction, 1);

      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QPH:
    {
      // First linear linear system
      J = NM_create(NM_SPARSE, m + 2*nd, m + 2*nd);
      J_nzmax = (d * d) * (m / d) + H_nzmax + 2 * (d * 3 - 2) * n + H_nzmax + nd;
      NM_triplet_alloc(J, J_nzmax);
      J->matrix2->origin = NSM_TRIPLET;

      NumericsMatrix * minusQpH = QNTpH(velocity, reaction, H, nd, n);
      NM_scal(-1.0, minusQpH);
      NumericsMatrix * minusQpHt = NM_transpose(minusQpH);

      NM_insert(J, M, 0, 0);
      NM_insert(J, minusQpH, m + nd, 0);
      NM_insert(J, eye_nd, m, m);
      NM_insert(J, eye_nd, m + nd, m);
      NM_insert(J, minusQpHt, 0, m + nd);
      NM_insert(J, eye_nd, m, m + nd);


      NM_free(minusQpH);
      NM_free(minusQpHt);

      double *r_check = (double*)calloc(nd,sizeof(double));
      QNTpinvz(velocity, reaction, reaction, nd, n, r_check);

      double *Qp_primalConstraint = (double*)calloc(nd,sizeof(double));
      QNTpz(velocity, reaction, primalConstraint, nd, n, Qp_primalConstraint);
      cblas_dcopy(nd, Qp_primalConstraint, 1, primalConstraint, 1);

      cblas_dcopy(m, dualConstraint, 1, rhs, 1);
      cblas_dcopy(nd, r_check, 1, rhs+m, 1);
      cblas_dcopy(nd, primalConstraint, 1, rhs+m+nd, 1);
      cblas_dscal(m + nd + nd, -1.0, rhs, 1);

      cblas_dcopy(m+2*nd, rhs, 1, rhs_2, 1);

      free(r_check);
      free(Qp_primalConstraint);
      jacobian_is_nan = NM_isnan(J);
      if (jacobian_is_nan)
    	{
    	  numerics_printf_verbose(0, "The Jacobian matrix J contains NaN");
    	  break;
    	}

      if ( options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] == SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES )
      {
      	double *rhs_save = (double*)calloc(m+2*nd,sizeof(double));
      	cblas_dcopy(m+2*nd, rhs, 1, rhs_save, 1);
      	NM_LDLT_refine(J, rhs, rhs_save, 1, 1e-14, 10, 0);
      	free(rhs_save);
      }
      else
      {
        NM_LDLT_solve(J, rhs, 1);
      }

      cblas_dcopy(m, rhs, 1, d_globalVelocity, 1);
      QNTpinvz(velocity, reaction, rhs+m, nd, n, d_velocity);
      QNTpz(velocity, reaction, rhs+m+nd, nd, n, d_reaction);

      /* NM_gemv(1.0, H, globalVelocity, 0.0, Hvw); */
      /* cblas_daxpy(nd, 1.0, w, 1, Hvw, 1); */
      /* NM_gemv(1.0, H, d_globalVelocity, 0.0, d_velocity);  // d_velocity <- H*dv */
      /* cblas_daxpy(nd, 1.0, Hvw, 1, d_velocity, 1);         // d_velocity <- H*v + w + H*dv */
      /* cblas_daxpy(nd, -1.0, velocity, 1, d_velocity, 1);   // d_velocity <- H*(v + dv) + w  - u */


      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QP2:
    {
      // First linear linear system
      JR = NM_create(NM_SPARSE, m + nd, m + nd);
      JR_nzmax = (d * d) * (m / d) + H_nzmax + H_nzmax + nd;
      NM_triplet_alloc(JR, JR_nzmax);
      JR->matrix2->origin = NSM_TRIPLET;

      Nesterov_Todd_vector(3, velocity, reaction, nd, n, p2);
      Qp2 = QRmat(p2, nd, n);

      NM_insert(JR, minus_M, 0, 0);
      NM_insert(JR, Ht, 0, m);
      NM_insert(JR, H, m, 0);
      NM_insert(JR, Qp2, m, m);

      if(Qp2) Qp2 = NM_free(Qp2);

      jacobian_is_nan = NM_isnan(JR);
      if (jacobian_is_nan)
  {
    numerics_printf_verbose(0, "The Jacobian matrix J contains NaN");
    break;
  }

      NV_insert(r_rhs, m + nd, dualConstraint, m, 0);
      NM_gemv(1.0, H, globalVelocity, 0.0, Hvw);       // Hvw <- H*v
      cblas_daxpy(nd, 1.0, w, 1, Hvw, 1);              // Hvw <- H*v + w. Hvw will be used for the computation of r_du
      NV_insert(r_rhs, m + nd, Hvw, nd, m);
      cblas_dscal(nd, -1.0, r_rhs+m, 1);

      cblas_dcopy(m+nd, r_rhs, 1, r_rhs_2, 1);

      NSM_linearSolverParams(JR)->solver = NSM_HSL;
      NM_LDLT_solve(JR, r_rhs, 1);

      cblas_dcopy(m, r_rhs, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, r_rhs+m, 1, d_reaction, 1);
      NM_gemv(1.0, H, d_globalVelocity, 0.0, d_velocity);  // d_velocity <- H*dv
      cblas_daxpy(nd, 1.0, Hvw, 1, d_velocity, 1);         // d_velocity <- H*v + w + H*dv
      cblas_daxpy(nd, -1.0, velocity, 1, d_velocity, 1);   // d_velocity <- H*(v + dv) + w  - u

      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QPH:
    {
      t1 = clock();
      // First linear linear system
      JR = NM_create(NM_SPARSE, m + nd, m + nd);
      JR_nzmax = (d * d) * (m / d) + H_nzmax + H_nzmax + nd;
      NM_triplet_alloc(JR, JR_nzmax);
      JR->matrix2->origin = NSM_TRIPLET;

      NumericsMatrix * QpH = QNTpH(velocity, reaction, H, nd, n);
      NumericsMatrix * QpHt = NM_transpose(QpH);

      NM_insert(JR, minus_M, 0, 0);
      NM_insert(JR, QpH, m, 0);
      NM_insert(JR,QpHt, 0, m);
      NM_insert(JR, eye_nd, m, m);

      //      if (iteration == 0) printf("NNZ2X2_QPH %zu\n",NM_nnz(JR));

      if(QpH) QpH = NM_free(QpH);
      if(QpHt) QpHt = NM_free(QpHt);

      jacobian_is_nan = NM_isnan(JR);
      if (jacobian_is_nan)
      {
        numerics_printf_verbose(0, "The Jacobian matrix J contains NaN");
        break;
      }

      NV_insert(r_rhs, m + nd, dualConstraint, m, 0);
      NM_gemv(1.0, H, globalVelocity, 0.0, Hvw);
      cblas_daxpy(nd, 1.0, w, 1, Hvw, 1);              // Hvw will be used for the computation of r_du
      QNTpz(velocity, reaction, Hvw, nd, n, r_rhs+m);
      cblas_dscal(nd, -1.0, r_rhs+m, 1);

      cblas_dcopy(m+nd, r_rhs, 1, r_rhs_2, 1);


      if ( options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] == SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES )
      {
        double *rhs_save = (double*)calloc(m+nd,sizeof(double));
        cblas_dcopy(m+nd, r_rhs, 1, rhs_save, 1);
        NM_LDLT_refine(JR, r_rhs, rhs_save, 1, 1e-14, 10, 0);
        free(rhs_save);
      }
      else
      {
        NM_LDLT_solve(JR, r_rhs, 1);
      }


      // NSM_linearSolverParams(JR)->solver = NSM_HSL;
      // NM_LDLT_solve(JR, r_rhs, 1);




      cblas_dcopy(m, r_rhs, 1, d_globalVelocity, 1);
      QNTpz(velocity, reaction, r_rhs+m, nd, n, d_reaction);
      NM_gemv(1.0, H, d_globalVelocity, 0.0, d_velocity);  // d_velocity <- H*dv
      cblas_daxpy(nd, 1.0, Hvw, 1, d_velocity, 1);         // d_velocity <- H*v + w + H*dv
      cblas_daxpy(nd, -1.0, velocity, 1, d_velocity, 1);   // d_velocity <- H*(v + dv) + w  - u

      t2 = clock();
      total_time += (double)(t2-t1)/(double)clk_tck;

      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_QPH:
    {
      // First linear linear system
      JR = NM_create(NM_SPARSE, nd, nd);
      //JR_nzmax = 2*H_nzmax + nd;
      //NM_triplet_alloc(JR, JR_nzmax);
      //JR->matrix2->origin = NSM_TRIPLET;

      NumericsMatrix * JR_a = NULL;
      NumericsMatrix * JR_b = NULL;

      NumericsMatrix * QpH = QNTpH(velocity, reaction, H, nd, n);
      NumericsMatrix * QpHt = NM_transpose(QpH);

      JR_a = NM_multiply(Minv, QpHt);
      JR_b = NM_multiply(QpH, JR_a);
      JR = NM_add(1.0, JR_b, 1.0, eye_nd);
      JR_a = NM_free(JR_a);
      JR_b = NM_free(JR_b);

      jacobian_is_nan = NM_isnan(JR);
      if (jacobian_is_nan)
      {
        numerics_printf_verbose(0, "The Jacobian matrix J contains NaN");
        break;
      }

      //      if (iteration == 0) printf("NNZ1X1_QPH %zu\n",NM_nnz(JR));

      double * fHr = (double*)calloc(m,sizeof(double));
      double * MfHr = (double*)calloc(m,sizeof(double));
      double * vdv = (double*)calloc(m,sizeof(double));
      double * rdr = (double*)calloc(nd,sizeof(double));
      double * w_h = (double*)calloc(nd,sizeof(double));

      cblas_dcopy(m, f, 1, fHr, 1);
      NM_gemv(1.0, Ht, reaction, 1.0, fHr);         // fHr <- H'*r + f
      NM_gemv(1.0, Minv, fHr, 0.0, MfHr);           // MfHr <- M^{-1}* (H'*r + f)
      NM_gemv(-1.0, QpH, MfHr, 0.0, sr_rhs);
      QNTpz(velocity, reaction, w, nd, n, w_h);
      cblas_daxpy(nd, -1.0, w_h, 1, sr_rhs, 1);

      cblas_dcopy(nd, sr_rhs, 1, sr_rhs_2, 1);

      /* NM_Cholesky_solve(JR, sr_rhs, 1); */

      NSM_linearSolverParams(JR)->solver = NSM_HSL;
      NM_LDLT_solve(JR, sr_rhs, 1);

      QNTpz(velocity, reaction, sr_rhs, nd, n, d_reaction);

      NV_add(reaction, d_reaction, nd, rdr);
      cblas_dcopy(m, f, 1, MfHr, 1);
      NM_gemv(1.0, Ht, rdr, 1.0, MfHr);
      NM_gemv(1.0, Minv, MfHr, 0.0, d_globalVelocity);
      cblas_daxpy(m, -1.0, globalVelocity, 1, d_globalVelocity, 1);

      NV_add(globalVelocity, d_globalVelocity, m, vdv);
      cblas_dcopy(nd, w, 1, d_velocity, 1);
      NM_gemv(1.0, H, vdv, 1.0, d_velocity);
      cblas_daxpy(nd, -1.0, velocity, 1, d_velocity, 1);

      free(fHr);
      free(MfHr);
      free(vdv);
      free(rdr);
      free(w_h);
      QpH = NM_free(QpH);
      QpHt = NM_free(QpHt);

      break;
    }

    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_TEST6:
    {
      // First linear linear system
      J = NM_create(NM_SPARSE, m + 2*nd, m + 2*nd);
      J_nzmax = (d * d) * (m / d) + H_nzmax + (3*d-2)*n + 3*d*n + H_nzmax + nd;
      NM_triplet_alloc(J, J_nzmax);
      J->matrix2->origin = NSM_TRIPLET;

      NumericsMatrix * arrow_r = Arrow_repr(reaction, nd, n);
      NumericsMatrix * arrow_u = Arrow_repr(velocity, nd, n) ;

      NM_insert(J, M, 0, 0);
      NM_insert(J, minus_H, m + nd, 0);
      NM_insert(J, arrow_r, m, m);
      NM_insert(J, eye_nd, m + nd, m);
      NM_insert(J, minus_Ht, 0, m + nd);
      NM_insert(J, arrow_u, m, m + nd);

      NM_free(arrow_r);
      NM_free(arrow_u);


      jacobian_is_nan = NM_isnan(J);
      if (jacobian_is_nan)
      {
        numerics_printf_verbose(0, "The Jacobian matrix J contains NaN");
        break;
      }

      cblas_dcopy(m, dualConstraint, 1, rhs, 1);
      cblas_dcopy(nd, complemConstraint, 1, rhs+m, 1);
      cblas_dcopy(nd, primalConstraint, 1, rhs+m+nd, 1);

      cblas_dscal(m + nd + nd, -1.0, rhs, 1);

      cblas_dcopy(m+2*nd, rhs, 1, rhs_2, 1);

      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] == SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES)
      {
        double *rhs_save = (double *)calloc(m + 2 * nd, sizeof(double));
        double residu;
        cblas_dcopy(m + 2 * nd, rhs, 1, rhs_save, 1);
        NM_LU_refine(J, rhs,  1e-14, 10,  &residu);
        free(rhs_save);
      }
      else
        NM_LU_solve(J, rhs, 1);

      cblas_dcopy(m, rhs, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, rhs+m, 1, d_velocity, 1);
      cblas_dcopy(nd, rhs+m+nd, 1, d_reaction, 1);

      double * rhs_tmp = (double*)calloc(m+2*nd,sizeof(double));
      cblas_dcopy(m+2*nd, rhs_2, 1, rhs_tmp, 1);

      NM_gemv(1.0, J, rhs, -1.0, rhs_tmp);
      LS_norm_d = cblas_dnrm2(m, rhs_tmp, 1);
      LS_norm_c = cblas_dnrm2(nd, rhs_tmp+m, 1);
      LS_norm_p = cblas_dnrm2(nd, rhs_tmp+m+nd, 1);

      free(rhs_tmp);

      break;
    }

    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QPH_TEST:
    {
      // First linear linear system
      JR = NM_create(NM_SPARSE, m + nd, m + nd);
      JR_nzmax = (d * d) * (m / d) + H_nzmax + H_nzmax + nd;
      NM_triplet_alloc(JR, JR_nzmax);
      JR->matrix2->origin = NSM_TRIPLET;

      NumericsMatrix * QpH = QNTpH(velocity, reaction, H, nd, n);
      NumericsMatrix * QpHt = NM_transpose(QpH);

      NM_insert(JR, minus_M, 0, 0);
      NM_insert(JR, QpH, m, 0);
      NM_insert(JR,QpHt, 0, m);
      NM_insert(JR, eye_nd, m, m);

      if(QpH) QpH = NM_free(QpH);
      if(QpHt) QpHt = NM_free(QpHt);

      jacobian_is_nan = NM_isnan(JR);
      if (jacobian_is_nan)
      {
        numerics_printf_verbose(0, "The Jacobian matrix J contains NaN");
        break;
      }

      double * r_Qp_u = (double*)calloc(nd,sizeof(double));
      double * r_Qp_u_inv = (double*)calloc(nd,sizeof(double));
      QNTpz(velocity, reaction, velocity, nd, n, r_Qp_u);  // r_Qp_u <- u_hat
      JA_inv(r_Qp_u, nd, n, r_Qp_u_inv);                   // r_Qp_u_inv <- u_hat_inv
      cblas_dscal(nd, -2.*barr_param, r_Qp_u_inv, 1);      // r_Qp_u_inv <= -2 * mu * u_hat_inv

      cblas_dcopy(m, dualConstraint, 1, r_rhs, 1);
      // NV_insert(r_rhs, m + nd, dualConstraint, m, 0);
      NM_gemv(1.0, H, globalVelocity, 0.0, Hvw);
      cblas_daxpy(nd, 1.0, w, 1, Hvw, 1);              // Hvw will be used for the computation of r_du
      QNTpz(velocity, reaction, Hvw, nd, n, r_rhs+m);
      cblas_daxpy(nd, 1.0, r_Qp_u_inv, 1, r_rhs+m, 1);
      cblas_dscal(nd, -1.0, r_rhs+m, 1);

      cblas_dcopy(m+nd, r_rhs, 1, r_rhs_2, 1);

      NSM_linearSolverParams(JR)->solver = NSM_HSL;
      NM_LDLT_solve(JR, r_rhs, 1);

      cblas_dcopy(m, r_rhs, 1, d_globalVelocity, 1);
      QNTpz(velocity, reaction, r_rhs+m, nd, n, d_reaction);
      NM_gemv(1.0, H, d_globalVelocity, 0.0, d_velocity);  // d_velocity <- H*dv
      cblas_daxpy(nd, 1.0, Hvw, 1, d_velocity, 1);         // d_velocity <- H*v + w + H*dv
      cblas_daxpy(nd, -1.0, velocity, 1, d_velocity, 1);   // d_velocity <- H*(v + dv) + w  - u

      double * rhs_tmp = (double*)calloc(m+nd,sizeof(double));
      cblas_dcopy(m+nd, r_rhs_2, 1, rhs_tmp, 1);
      NM_gemv(1.0, JR, r_rhs, -1.0, rhs_tmp);
      LS_norm_d = cblas_dnrm2(m, rhs_tmp, 1);
      LS_norm_c = cblas_dnrm2(nd, rhs_tmp+m, 1);

      free(rhs_tmp);
      free(r_Qp_u);
      free(r_Qp_u_inv);
      if (JR) JR = NM_free(JR);

      break;
    }

    default:
    {
      printf("ERROR 1\n");
    }
    }
    if (jacobian_is_nan)
    {
    	hasNotConverged = 2;
    	if (J) J = NM_free(J);
      if (JR) JR = NM_free(JR);
    	break;
    }

    /* computing the affine step-length */
    alpha_primal = getStepLength(velocity, d_velocity, nd, n, gmm);
    // alpha_primal = getStepLength(u_plus_phiu, d_velocity, nd, n, 1.0); //gmm);
    alpha_dual = getStepLength(reaction, d_reaction, nd, n, gmm);


    if (alpha_primal < alpha_dual)
      alpha_dual = alpha_primal;
    else
      alpha_primal = alpha_dual;

    /* updating the gamma parameter used to compute the step-length */
    gmm = gmmp1 + gmmp2 * fmin(alpha_primal, alpha_dual);

    /* ----- Predictor step of Mehrotra ----- */
    // cblas_dcopy(nd, u_plus_phiu, 1, v_plus_dv, 1);
    cblas_dcopy(nd, velocity, 1, v_plus_dv, 1);
    cblas_dcopy(nd, reaction, 1, r_plus_dr, 1);
    cblas_daxpy(nd, alpha_primal, d_velocity, 1, v_plus_dv, 1);
    cblas_daxpy(nd, alpha_dual, d_reaction, 1, r_plus_dr, 1);

    /* affine barrier parameter */
    // barr_param_a = cblas_ddot(nd, v_plus_dv, 1, r_plus_dr, 1) / n / 3.;
    // barr_param_a = cblas_ddot(nd, v_plus_dv, 1, r_plus_dr, 1) / n ;

    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] > 0)
    {
      barr_param_a = complemResidualNorm(v_plus_dv, r_plus_dr, nd, n)/n;
      // barr_param_a = cblas_ddot(nd, v_plus_dv, 1, r_plus_dr, 1) / n;
    }
    else
      barr_param_a = cblas_ddot(nd, v_plus_dv, 1, r_plus_dr, 1) / n / 3.;



    /* computing the centralization parameter */
    e = barr_param > sgmp1 ? fmax(1.0, sgmp2 * pow(fmin(alpha_primal, alpha_dual),2)) : sgmp3;
    sigma = fmin(1.0, pow(barr_param_a / barr_param, e));

    /* Computing the corrector step of Mehrotra algorithm */

    /*      without NT scaling */
    /*      rhs = - */
    /*      [         M*v - H'*r - f             ]  m         dualConstraint */
    /*      [ u o r + da^u o da^r - 2*sigma*mu*e ]  nd        complemConstraint */
    /*      [         u - H*v - w                ]  nd        primalConstraint */

    /*      with NT scaling */
    /*      rhs = - */
    /*      [ M*v - H'*r - f ]  m         dualConstraint */
    /*      [       r        ]  nd        complemConstraint */
    /*      [  u - H*v - w   ]  nd        primalConstraint */

    switch ( options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] )
    {
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL:
      {
        t1 = clock();
      	// Second linear linear system

      	// // Create eye_phiu = I + phi(u)
      	// eye_phiu = NM_create(NM_SPARSE, nd, nd);
      	// size_t eye_phiu_nzmax = nd + 2*n;
      	// NM_triplet_alloc(eye_phiu, eye_phiu_nzmax);
      	// NM_fill(eye_phiu, NM_SPARSE, nd, nd, eye_phiu->matrix2);

      	// /* Arrow matrix filling */
      	// size_t pos; double ub;
      	// for(size_t i = 0; i < n; ++i)
      	// {
      	//   pos = i * d;
      	//   ub = sqrt(velocity[pos+1]*velocity[pos+1]+velocity[pos+2]*velocity[pos+2]);

      	//   NM_entry(eye_phiu, pos, pos, 1.);
      	//   NM_entry(eye_phiu, pos, pos+1, velocity[pos+1]/ub);
      	//   NM_entry(eye_phiu, pos, pos+2, velocity[pos+2]/ub);
      	//   NM_entry(eye_phiu, pos+1, pos+1, 1.);
      	//   NM_entry(eye_phiu, pos+2, pos+2, 1.);
      	// }

      	// // d_velocity = (I + phi(u))d^u
      	// NM_gemv(1.0, eye_phiu, rhs+m, 0.0, d_velocity);

      	JA_prod(d_velocity, d_reaction, nd, n, dvdr_jprod);
        cblas_daxpy(nd, -1.0, dvdr_jprod, 1, rhs_2 + m, 1);
        for (unsigned int k = 0; k < nd; rhs_2[m + k] += 2 * sigma * barr_param, k += d);


        double *rhs_tmp = (double *)calloc(m + 2 * nd, sizeof(double));
        cblas_dcopy(m + 2 * nd, rhs_2, 1, rhs_tmp, 1);

        /* double *sol = (double *)calloc(m + 2 * nd, sizeof(double)); */

        /* for (int itr = 0; itr < 10; itr++) { */
        /*   NM_LU_solve(J, rhs_tmp, 1); */

        /*   cblas_daxpy(m + 2 * nd, 1.0, rhs_tmp, 1, sol, 1); */
        /*   cblas_dcopy(m + 2 * nd, rhs_2, 1, rhs_tmp, 1); */
        /*   NM_gemv(-1.0, J, sol, 1.0, rhs_tmp); */
        /*   printf("refinement iterations = %i %8.2e < %8.2e\n",itr, cblas_dnrm2(m+2*nd, rhs_tmp, 1), tol); */
	/*   printf("J size %i %i\n", J->size0, J->size1); */
	/*   NM_write_in_filename(J, "J_IPM_LU.txt"); */
        /*   if (cblas_dnrm2(m + 2 * nd, rhs_tmp, 1) <= tol) { */
        /*     break; */
        /*   } */
	/*   getchar(); */
        /* } */
        /* cblas_dcopy(m, sol, 1, d_globalVelocity, 1); */
        /* cblas_dcopy(nd, sol + m, 1, d_velocity, 1); */
        /* cblas_dcopy(nd, sol + m + nd, 1, d_reaction, 1); */

        /* NM_gemv(1.0, J, sol, -1.0, rhs_2); */

        if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] ==
            SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES)
        {
	       double *rhs_save = (double *)calloc(m + 2 * nd, sizeof(double));
          cblas_dcopy(m + 2 * nd, rhs_tmp, 1, rhs_save, 1);
          //NM_LU_refine(J, rhs, rhs_save, 1, tol, 10, 0);
      	  double residu;
      	  NM_LU_refine(J, rhs_tmp,  1e-14, 10,  &residu);
          free(rhs_save);
      	}
      	else
      	  NM_LU_solve(J, rhs_tmp, 1);

        cblas_dcopy(m, rhs_tmp, 1, d_globalVelocity, 1);
        cblas_dcopy(nd, rhs_tmp + m, 1, d_velocity, 1);
        cblas_dcopy(nd, rhs_tmp + m + nd, 1, d_reaction, 1);
        NM_gemv(1.0, J, rhs_tmp, -1.0, rhs_2);

        LS_norm_d = cblas_dnrm2(m, rhs_2, 1);
        LS_norm_c = cblas_dnrm2(nd, rhs_2 + m, 1);
        LS_norm_p = cblas_dnrm2(nd, rhs_2 + m + nd, 1);

        t2 = clock();
        total_time += (double)(t2-t1)/(double)clk_tck;

        free(rhs_tmp);
        //free(sol);

        J = NM_free(J);

        break;
      }

    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QP2:
    {
      // Second linear linear system
      QNTpz(velocity, reaction, d_velocity, nd, n, d_velocity_t);
      QNTpinvz(velocity, reaction, d_reaction, nd, n, d_reaction_t);
      JA_prod(d_velocity_t, d_reaction_t, nd, n, dvdr_jprod);
      QNTpz(velocity, reaction, velocity, nd, n, velocity_t);
      JA_inv(velocity_t, nd, n, velocity_t_inv);
      JA_prod(velocity_t_inv, dvdr_jprod, nd, n, d_velocity_t);
      QNTpz(velocity, reaction, d_velocity_t, nd, n, dvdr_jprod);
      cblas_daxpy(nd, -1.0, dvdr_jprod, 1, rhs_2 + m, 1);

      JA_inv(velocity, nd, n, d_velocity_t);
      for (unsigned int k = 0; k < nd;
           rhs_2[m + k] += 2 * sigma * barr_param * d_velocity_t[k], k++);

      /* double *rhs_save = (double*)calloc(m+2*nd,sizeof(double)); */
      /* cblas_dcopy(m+2*nd, rhs_2, 1, rhs_save, 1); */

      double *rhs_tmp = (double *)calloc(m + 2 * nd, sizeof(double));
      cblas_dcopy(m + 2 * nd, rhs_2, 1, rhs_tmp, 1);


      if ( options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] == SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES )
      {
        double *rhs_save = (double*)calloc(m+2*nd,sizeof(double));
        cblas_dcopy(m+2*nd, rhs_2, 1, rhs_save, 1);
        NM_LDLT_refine(J, rhs_2, rhs_save, 1, 1e-14, 10, 0);
        free(rhs_save);
      }
      else
      {
        NM_LDLT_solve(J, rhs_2, 1);
      }

      cblas_dcopy(m, rhs_2, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, rhs_2 + m, 1, d_velocity, 1);
      cblas_dcopy(nd, rhs_2 + m + nd, 1, d_reaction, 1);

      NM_gemv(1.0, J, rhs_2, -1.0, rhs_tmp);

      LS_norm_d = cblas_dnrm2(m, rhs_tmp, 1);
      LS_norm_c = cblas_dnrm2(nd, rhs_tmp + m, 1);
      LS_norm_p = cblas_dnrm2(nd, rhs_tmp + m + nd, 1);

      free(rhs_tmp);

      J = NM_free(J);

      break;
    }


    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QPH:
    {
      // Second linear linear system
      JA_prod(rhs+m, rhs+m+nd, nd, n, dvdr_jprod);
      for (int k = 0; k < nd; dvdr_jprod[k] -=  2*sigma*barr_param, k+=d);

      QNTpz(velocity, reaction, velocity, nd, n, velocity_t);
      JA_inv(velocity_t, nd, n, velocity_t_inv);
      JA_prod(velocity_t_inv, dvdr_jprod, nd, n, d_velocity_t);
      cblas_daxpy(nd, -1.0, d_velocity_t, 1, rhs_2+m, 1);

      double * rhs_tmp = (double*)calloc(m+2*nd,sizeof(double));
      cblas_dcopy(m+2*nd, rhs_2, 1, rhs_tmp, 1);

      if ( options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] == SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES )
      {
      	double *rhs_save = (double*)calloc(m+2*nd,sizeof(double));
      	cblas_dcopy(m+2*nd, rhs_2, 1, rhs_save, 1);
      	NM_LDLT_refine(J, rhs_2, rhs_save, 1, 1e-14, 10, 0);
      	free(rhs_save);
      }
      else
      {
        NM_LDLT_solve(J, rhs_2, 1);
      }

      NM_gemv(1.0, J, rhs_2, -1.0, rhs_tmp);

      LS_norm_d = cblas_dnrm2(m, rhs_tmp, 1);
      LS_norm_c = cblas_dnrm2(nd, rhs_tmp+m, 1);
      LS_norm_p = cblas_dnrm2(nd, rhs_tmp+m+nd, 1);

      free(rhs_tmp);

      /* cblas_dcopy(m+2*nd, rhs_tmp_2, 1, rhs_tmp, 1); */
      /* NM_gemv(1.0, J, rhs_2, -1.0, rhs_tmp); */
      /* cblas_dscal(m+2*nd, -1.0, rhs_tmp, 1); */
      /* NM_LDLT_solve(J, rhs_tmp, 1); */
      /* cblas_daxpy(m+2*nd, 1.0, rhs_tmp, 1, rhs_2, 1); */

      //NM_LU_solve(J, rhs_2, 1);


      cblas_dcopy(m, rhs_2, 1, d_globalVelocity, 1);
      QNTpinvz(velocity, reaction, rhs_2+m, nd, n, d_velocity);
      QNTpz(velocity, reaction, rhs_2+m+nd, nd, n, d_reaction);

      /* NM_gemv(1.0, H, globalVelocity, 0.0, Hvw); */
      /* cblas_daxpy(nd, 1.0, w, 1, Hvw, 1); */
      /* NM_gemv(1.0, H, d_globalVelocity, 0.0, d_velocity);  // d_velocity <- H*dv */
      /* cblas_daxpy(nd, 1.0, Hvw, 1, d_velocity, 1);         // d_velocity <- H*v + w + H*dv */
      /* cblas_daxpy(nd, -1.0, velocity, 1, d_velocity, 1);   // d_velocity <- H*(v + dv) + w  - u */

      /* double * du = (double*)calloc(nd,sizeof(double)); */
      /* QNTpz(velocity, reaction, d_velocity, nd, n, du); */
      /* //printf("|d_velocity-du| = %9.2e\n",norm2VecDiff(rhs_2+m, du, nd)); */
      /* free(du); */

      //double toto = cblas_dnrm2(m+2*nd, rhs_tmp, 1);
      //printf("norm(rhs_tmp) = %9.2e\n", cblas_dnrm2(m+2*nd, rhs_tmp, 1));
      //NM_gemv(1.0, J, rhs_2, -1.0, rhs_tmp);
      /* printf("rhs_tmp = %9.2e\n",relative_error_linear_system_solution(J, rhs_2, rhs_tmp_2, m+2*nd, 0, m+2*nd)); */
      /* printf("error 1 = %9.2e\n", relative_error_linear_system_solution(J, rhs_2, rhs_tmp_2, m+2*nd, 0, m)); */
      /* printf("error 2 = %9.2e\n", relative_error_linear_system_solution(J, rhs_2, rhs_tmp_2, m+2*nd, m, nd)); */
      /* printf("error 3 = %9.2e\n", relative_error_linear_system_solution(J, rhs_2, rhs_tmp_2, m+2*nd, m+nd, nd)); */

      /* double * du = (double*)calloc(nd,sizeof(double)); */
      /* NM_gemv(1.0, H, globalVelocity, 0.0, Hvw); */
      /* cblas_daxpy(nd, 1.0, w, 1, Hvw, 1); */
      /* NM_gemv(1.0, H, d_globalVelocity, 0.0, du);  // d_velocity <- H*dv */
      /* cblas_daxpy(nd, 1.0, Hvw, 1, du, 1);         // d_velocity <- H*v + w + H*dv */
      /* cblas_daxpy(nd, -1.0, velocity, 1, du, 1);   // d_velocity <- H*(v + dv) + w  - u */
      /* double * pinv = (double*)calloc(nd,sizeof(double)); */
      /* Nesterov_Todd_vector(1, velocity, reaction, nd, n, pinv); */

      /* double l_pinv = 1e300; */
      /* double lp; */
      /* int j;  */
      /* for (int k = 0; k < n; k++) */
      /* { */
      /* 	j = k*d; */
      /* 	lp = pinv[j] - cblas_dnrm2(d-1, pinv+j+1, 1); */
      /* 	l_pinv = (l_pinv>lp) ? lp : l_pinv;  */
      /* } */
      /* printf("|d_velocity-du| = %9.2e %9.2e\n",norm2VecDiff(d_velocity,du,nd), l_pinv); */

      /* free(pinv); */
      /* free(du); */

      /* double * vdv = (double*)calloc(m, sizeof(double)); */
      /* NV_add(globalVelocity, d_globalVelocity, m, vdv); */
      /* NM_gemv(1.0, H, vdv, 0.0, d_velocity); */
      /* cblas_daxpy(nd, 1.0, w, 1, d_velocity, 1); */
      /* cblas_daxpy(nd, -1.0, velocity, 1, d_velocity, 1); */
      /* free(vdv); */

      if(J) J = NM_free(J);

      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QP2:
    {
      // Second linear linear system
      double * r_Qp_u = (double*)calloc(nd,sizeof(double));
      double * r_Qp_du = (double*)calloc(nd,sizeof(double));
      double * r_dudr = (double*)calloc(nd,sizeof(double));
      double * r_ududr = (double*)calloc(nd,sizeof(double));
      double * r_Qpinv_dr = (double*)calloc(nd,sizeof(double));

      QNTpz(velocity, reaction, d_velocity, nd, n, r_Qp_du);          // du_affine_hat
      QNTpinvz(velocity, reaction, d_reaction, nd, n, r_Qpinv_dr);     // dr_affine_check
      JA_prod(r_Qp_du, r_Qpinv_dr, nd, n, r_dudr);                    // Jordan product du_affine__hat o dr_affine_check
      for (int k = 0; k < nd; r_dudr[k] -= 2*sigma*barr_param, k+=d); // du_affine_hat o dr_affine_check - 2*sigma*mu*e
      QNTpz(velocity, reaction, velocity, nd, n, r_Qp_u);             // u_hat
      Jxinvprody(r_Qp_u, r_dudr, nd, n, r_ududr);                     // Jordan product u_hat_inv o (du_affine__hat o dr_affine_check)
      QNTpinvz(velocity, reaction, r_ududr, nd, n, r_ududr);           // Qp_inv * u_hat_inv o du_affine_hat o dr_affine_check - 2*sigma*mu*u_hatinv)
      cblas_daxpy(nd, -1.0, r_ududr, 1, r_rhs_2+m, 1);

      double * rhs_tmp = (double*)calloc(m+nd,sizeof(double));
      cblas_dcopy(m+nd, r_rhs_2, 1, rhs_tmp, 1);

      NSM_linearSolverParams(JR)->solver = NSM_HSL;
      NM_LDLT_solve(JR, r_rhs_2, 1);

      cblas_dcopy(m, r_rhs_2, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, r_rhs_2+m, 1, d_reaction, 1);
      NM_gemv(1.0, H, d_globalVelocity, 0.0, d_velocity);
      cblas_daxpy(nd, 1.0, Hvw, 1, d_velocity, 1);
      cblas_daxpy(nd, -1.0, velocity, 1, d_velocity, 1);

      NM_gemv(1.0, JR, r_rhs_2, -1.0, rhs_tmp);
      LS_norm_d = cblas_dnrm2(m, rhs_tmp, 1);
      LS_norm_c = cblas_dnrm2(nd, rhs_tmp+m, 1);

      free(rhs_tmp);

      if(JR) JR = NM_free(JR);
      free(r_Qp_u);
      free(r_Qp_du);
      free(r_dudr);
      free(r_ududr);
      free(r_Qpinv_dr);

      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QPH:
    {
      t1 = clock();

      // Second linear linear system
      double * r_Qp_u = (double*)calloc(nd,sizeof(double));
      double * r_Qp_du = (double*)calloc(nd,sizeof(double));
      double * r_dudr = (double*)calloc(nd,sizeof(double));
      double * r_ududr = (double*)calloc(nd,sizeof(double));


      QNTpz(velocity, reaction, velocity, nd, n, r_Qp_u);  // r_Qp_u <- u_hat  (Qp formula)
      cblas_dcopy(nd, r_Qp_u, 1, r_Qp_du, 1);         // r_Qp_du <- u_hat
      cblas_daxpy(nd, 1.0, r_rhs+m, 1, r_Qp_du, 1);   // r_Qp_du <- dr_check + u_hat
      cblas_dscal(nd, -1.0, r_Qp_du, 1);              // r_Qp_du <- -(dr_check + u_hat) = du_hat.
	                                                  // Here we use the following formula satisfied by the affine step
	                                                  //    u_hat o dr_check + r_check o du_hat = - u_hat o r_check
	                                                  // Then we use the fact that for the NT scaling we have: u_hat = r_check.
	                                                  // Therefore du_hat = -dr_check - r_check = -dr_check - u_hat

      JA_prod(r_Qp_du, r_rhs+m, nd, n, r_dudr);                         // r_dudr <- du_hat o dr_check

      //cblas_dscal(nd, 2, r_dudr, 1);

      for (int k = 0; k < nd; r_dudr[k] -= 2*sigma*barr_param, k+=d);   // r_dudr <- du_hat o dr_check - 2*sigma*mu*e
      Jxinvprody(r_Qp_u, r_dudr, nd, n, r_ududr);                       // r_ududr <- u_hat^{-1} o  (du_hat o dr_check - 2*sigma*mu*e)

      cblas_daxpy(nd, -1.0, r_ududr, 1, r_rhs_2+m, 1);                  // r_rhs_2+m <- -Qp*(Hv+w) - u_hat^{-1} o (du_hat o dr_check) + 2*sigma*mu*u_hat^{-1}

      double * rhs_tmp = (double*)calloc(m+nd,sizeof(double));
      cblas_dcopy(m+nd, r_rhs_2, 1, rhs_tmp, 1);

      if ( options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] == SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES )
      {
      	double *rhs_save = (double*)calloc(m+nd,sizeof(double));
      	cblas_dcopy(m+nd, r_rhs_2, 1, rhs_save, 1);
      	NM_LDLT_refine(JR, r_rhs_2, rhs_save, 1, 1e-14, 10, 0);
      	free(rhs_save);
      }
      else
      {
        NM_LDLT_solve(JR, r_rhs_2, 1);
      }

      /* NSM_linearSolverParams(JR)->solver = NSM_HSL; */
      /* NM_LDLT_solve(JR, r_rhs_2, 1); */

      cblas_dcopy(m, r_rhs_2, 1, d_globalVelocity, 1);
      QNTpz(velocity, reaction, r_rhs_2+m, nd, n, d_reaction);
      NM_gemv(1.0, H, d_globalVelocity, 0.0, d_velocity);
      cblas_daxpy(nd, 1.0, Hvw, 1, d_velocity, 1);
      cblas_daxpy(nd, -1.0, velocity, 1, d_velocity, 1);

      NM_gemv(1.0, JR, r_rhs_2, -1.0, rhs_tmp);
      LS_norm_d = cblas_dnrm2(m, rhs_tmp, 1);
      LS_norm_c = cblas_dnrm2(nd, rhs_tmp+m, 1);
      free(rhs_tmp);


      t2 = clock();
      total_time += (double)(t2-t1)/(double)clk_tck;



      if (JR) JR = NM_free(JR);
      free(r_Qp_u);
      free(r_Qp_du);
      free(r_dudr);
      free(r_ududr);

      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_1X1_QPH:
    {
      // Second linear linear system
      double * vdv = (double*)calloc(m,sizeof(double));
      double * rdr = (double*)calloc(nd,sizeof(double));
      double * rhs_cor = (double*)calloc(nd,sizeof(double));
      double * MfHr = (double*)calloc(m,sizeof(double));

      /* NV_display(d_reaction, nd); */
      /* NV_display(d_globalVelocity, m); */
      /* NV_display(d_velocity, nd); */

      QNTpz(velocity, reaction, d_velocity, nd, n, d_velocity_t);
      JA_prod(d_velocity_t, sr_rhs, nd, n, dvdr_jprod);
      for (int k = 0; k < nd; dvdr_jprod[k] -= 2*sigma*barr_param, k+=d);
      QNTpz(velocity, reaction, velocity, nd, n, velocity_t);
      JA_inv(velocity_t, nd, n, velocity_t_inv);
      JA_prod(velocity_t_inv, dvdr_jprod, nd, n, rhs_cor);
      cblas_daxpy(nd, -1.0, rhs_cor, 1, sr_rhs_2, 1);

      double * rhs_tmp = (double*)calloc(nd,sizeof(double));
      cblas_dcopy(nd, sr_rhs_2, 1, rhs_tmp, 1);

      /* NM_Cholesky_solve(JR, sr_rhs_2, 1); */

      if ( options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] == SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES )
      {
      	double *rhs_save = (double*)calloc(nd,sizeof(double));
      	cblas_dcopy(nd, sr_rhs_2, 1, rhs_save, 1);
      	NM_LDLT_refine(JR, sr_rhs_2, rhs_save, 1, 1e-10, 10, 0);
      	free(rhs_save);
      }
      else
      {
        NM_LDLT_solve(JR, sr_rhs_2, 1);
      }


      /* NSM_linearSolverParams(JR)->solver = NSM_HSL; */
      /* NM_LDLT_solve(JR, sr_rhs_2, 1); */

      QNTpz(velocity, reaction, sr_rhs_2, nd, n, d_reaction);

      NV_add(reaction, d_reaction, nd, rdr);
      cblas_dcopy(m, f, 1, MfHr, 1);
      NM_gemv(1.0, Ht, rdr, 1.0, MfHr);
      NM_gemv(1.0, Minv, MfHr, 0.0, d_globalVelocity);
      cblas_daxpy(m, -1.0, globalVelocity, 1, d_globalVelocity, 1);

      NV_add(globalVelocity, d_globalVelocity, m, vdv);
      cblas_dcopy(nd, w, 1, d_velocity, 1);
      NM_gemv(1.0, H, vdv, 1.0, d_velocity);
      cblas_daxpy(nd, -1.0, velocity, 1, d_velocity, 1);

      NM_gemv(1.0, JR, sr_rhs_2, -1.0, rhs_tmp);
      LS_norm_c = cblas_dnrm2(nd, rhs_tmp, 1);
      free(rhs_tmp);

      JR = NM_free(JR);
      free(vdv);
      free(rdr);
      free(rhs_cor);
      free(MfHr);

      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_TEST6:
    {
      gmm = 0.99;
      J = NM_free(J);
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QPH_TEST:
    {
      gmm = 0.99;
      if (JR) JR = NM_free(JR);
      break;
    }
    default:
    {
      printf("ERROR 2\n");
    }
    }

    if (Qp)
      Qp = NM_free(Qp);
    if (Qpinv)
      Qpinv = NM_free(Qpinv);

    t1 = clock();


    alpha_primal = getStepLength(velocity, d_velocity, nd, n, gmm);
    // alpha_primal = getStepLength(u_plus_phiu, d_velocity, nd, n, gmm);
    alpha_dual = getStepLength(reaction, d_reaction, nd, n, gmm);
    alpha_tmp = alpha_dual;

    if (alpha_primal < alpha_dual)
      alpha_dual = alpha_primal;
    else
      alpha_primal = alpha_dual;

    /* updating the gamma parameter used to compute the step-length at the next iteration */

    gmm = gmmp1 + gmmp2 * fmin(alpha_primal, alpha_dual);

    t2 = clock();
    total_time += (double)(t2-t1)/(double)clk_tck;


    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_TEST6 ||
        options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QPH_TEST)
    {
      numerics_printf_verbose(-1, "| %3i%c| %8.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |",
                            iteration, fws, dualgap, pinfeas, dinfeas, udotr, diff_fixp, complem, projerr, barr_param, 2*n*barr_param, alpha_primal,
                            fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
                            fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
                            fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]),
          LS_norm_p, LS_norm_d, LS_norm_c);
    }
    else
      numerics_printf_verbose(-1, "| %3i%c| %8.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |",
                            iteration, fws, dualgap, pinfeas, dinfeas, udotr, complem, projerr, barr_param, alpha_primal, alpha_tmp, sigma,
                            fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
                            fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
                            fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]),
			    LS_norm_p, LS_norm_d, LS_norm_c);
      // numerics_printf_verbose(-1, "| %3i%c| %8.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |",
      //                       iteration, fws, dualgap, pinfeas, dinfeas, udotr, complem, projerr, barr_param, alpha_primal, alpha_tmp, sigma,
      //                       cblas_dnrm2(m, d_globalVelocity, 1),
      //                       cblas_dnrm2(nd, d_velocity, 1),
      //                       cblas_dnrm2(nd, d_reaction, 1),
      //     LS_norm_p, LS_norm_d, LS_norm_c);


    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
    {
      // fprintf(iterates,"%d %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e;\n",
      //         iteration, pinfeas, dinfeas, udotr, complem, projerr, barr_param, alpha_primal,
      //         fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
      //         fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
      //         fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]),
      //         fabs(d_s[cblas_idamax(n, d_s, 1)]),
      //         LS_norm_p, LS_norm_d, LS_norm_c);
      printIteresProbMatlabFile(iteration, globalVelocity, velocity, reaction, d_globalVelocity, d_velocity, d_reaction, d, n, m, iterates);
      fprintf(iterates, "pinfeas(%3i) = %.20e;\n", iteration+1, pinfeas);
      fprintf(iterates, "dinfeas(%3i) = %.20e;\n", iteration+1, dinfeas);
      fprintf(iterates, "udotr(%3i) = %.20e;\n", iteration+1, udotr);
      fprintf(iterates, "complem(%3i) = %.20e;\n", iteration+1, complem);
      fprintf(iterates, "LSp(%3i) = %.20e;\n", iteration+1, LS_norm_p);
      fprintf(iterates, "LSd(%3i) = %.20e;\n", iteration+1, LS_norm_d);
      fprintf(iterates, "LSc(%3i) = %.20e;\n", iteration+1, LS_norm_c);
      fprintf(iterates, "alpha(%3i) = %.20e;\n", iteration+1, alpha_primal);
      fprintf(iterates,"projerr(%3i,:) = [",iteration+1);
      for(int i = 0; i < n; i++)
      {
        fprintf(iterates, "%.20e, ", projectionError(velocity+i*d, reaction+i*d, 1, tol));
      }
      fprintf(iterates,"];\n");
      fprintf(iterates,"proj(%3i,:) = [",iteration+1);
      for(int i = 0; i < n; i++)
      {
        worktmp[0] = reaction[3*i] -  velocity[3*i] ;
        worktmp[1] = reaction[3*i+1] -  velocity[3*i+1] ;
        worktmp[2] = reaction[3*i+2] -  velocity[3*i+2] ;
        projectionOnCone(worktmp, 1.0);
        fprintf(iterates, "%.20e, %.20e, %.20e, ", worktmp[0], worktmp[1], worktmp[2]);
      }
      fprintf(iterates,"];\n");
      fprintf(iterates,"nrmU(%3i,:) = [",iteration+1);
      for(int i = 0; i < n; i++)
      {
        fprintf(iterates, "%.20e, ", cblas_dnrm2(d, velocity+i*d, 1));
      }
      fprintf(iterates,"];\n");
      fprintf(iterates,"nrmR(%3i,:) = [",iteration+1);
      for(int i = 0; i < n; i++)
      {
        fprintf(iterates, "%.20e, ", cblas_dnrm2(d, reaction+i*d, 1));
      }
      fprintf(iterates,"];\n");
      fprintf(iterates,"nrmDU(%3i,:) = [",iteration+1);
      for(int i = 0; i < n; i++)
      {
        fprintf(iterates, "%.20e, ", cblas_dnrm2(d, d_velocity+i*d, 1));
      }
      fprintf(iterates,"];\n");
      fprintf(iterates,"nrmDR(%3i,:) = [",iteration+1);
      for(int i = 0; i < n; i++)
      {
        fprintf(iterates, "%.20e, ", cblas_dnrm2(d, d_reaction+i*d, 1));
      }
      fprintf(iterates,"];\n");
      fprintf(iterates,"eig1U(%3i,:) = [",iteration+1);
      for(int i = 0; i < n; i++)
      {
        fprintf(iterates, "%.20e, ", velocity[i*d] + cblas_dnrm2(2, velocity+i*d+1, 1));
      }
      fprintf(iterates,"];\n");
      fprintf(iterates,"eig1R(%3i,:) = [",iteration+1);
      for(int i = 0; i < n; i++)
      {
        fprintf(iterates, "%.20e, ", reaction[i*d] + cblas_dnrm2(2, reaction+i*d+1, 1));
      }
      fprintf(iterates,"];\n");
      fprintf(iterates,"eig2U(%3i,:) = [",iteration+1);
      for(int i = 0; i < n; i++)
      {
        fprintf(iterates, "%.20e, ", velocity[i*d] - cblas_dnrm2(2, velocity+i*d+1, 1));
      }
      fprintf(iterates,"];\n");
      fprintf(iterates,"eig2R(%3i,:) = [",iteration+1);
      for(int i = 0; i < n; i++)
      {
        fprintf(iterates, "%.20e, ", reaction[i*d] - cblas_dnrm2(2, reaction+i*d+1, 1));
      }
      fprintf(iterates,"];\n");
    }


    /* ----- Update variables ----- */

    if (NV_isnan(d_globalVelocity, m) | NV_isnan(d_velocity, nd) | NV_isnan(d_reaction, nd))
    {
      hasNotConverged = 2;
      break;
    }

    // if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
    //   printIteresProbMatlabFile(iteration, globalVelocity, velocity, reaction, d_globalVelocity, d_velocity, d_reaction, d, n, m, iterates);


    for (unsigned int i=0; i<n; i++)
    {
      d_s[i] = cblas_ddot(2, velocity+i*d+1, 1, d_velocity+i*d+1, 1)/cblas_dnrm2(2, velocity+i*d+1, 1) - s[i] +  cblas_dnrm2(2, velocity+i*d+1, 1);
    }



    cblas_daxpy(m, alpha_primal, d_globalVelocity, 1, globalVelocity, 1);
    cblas_daxpy(nd, alpha_primal, d_velocity, 1, velocity, 1);
    cblas_daxpy(nd, alpha_dual, d_reaction, 1, reaction, 1);

    // cblas_daxpy(n, alpha_dual, d_s, 1, s, 1);

    /* for (unsigned int  i = 0; i<nd; i+=d) velocity[i] = (1+barr_param)*velocity[i]; */
    /* for (unsigned int  i = 0; i<nd; i+=d) reaction[i] = (1+barr_param)*reaction[i]; */

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
// NM_display(M);
  /* Checking strict complementarity */
  /* For each cone i from 1 to n, one checks if u+r is in the interior of the Lorentz cone */
  /* One first computes the 3 dimensional vector somme = (u+r)/norm(u+r) */
  /* Then one checks if somme[0] > sqrt(somme[1]^2 + somme[2]^2) + ceps */

  double somme[3];
  double dur[3];
  double cesp = sqrt(DBL_EPSILON);
  double ns, ndur;
  int nsc = 0;
  // int nN = 0;
  // int nB = 0;
  // int nR = 0;
  nN = nB = nR = 0;
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
    printf("Strict complementarity satisfied: %4i / %4i  %9.2e %9.2e  %4i %4i %4i\n", nsc, n, ns, ns/cblas_dnrm2(3, s, 1), nB, nN, nR);

  if (hasNotConverged == 0)
    printf("test IPM: success, ");
  else
    printf("test IPM: failure, ");

  printf("%.1e  %.1e,     iter=%3d, n=%4d, \t%s \n", totalresidual, projerr, iteration, n, strToken);


  printf("\n\n");



  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] > 0)
  {
    // if (totalresidual > tols_TEST7_scalarP || total_time > max_time || alpha_primal <= 1e-13)
    // {
    //   while ( num_TEST7_scalarP < 11)
    //   {
    //     sprintf(name_file_TEST7_scalarP, "ipm-scalarP-%02d.pp", num_TEST7_scalarP);
    //     file_TEST7_scalarP = fopen(name_file_TEST7_scalarP, "a+");
    //     fprintf(file_TEST7_scalarP, "sumry: 1  %.2e   %5i %5i %5i %4i %4i %4i %4i  %.6f   %s\n",
    //             totalresidual, iteration, n, m, nB, nN, nR, n-nB-nN-nR,
    //             total_time, strToken);
    //     fclose(file_TEST7_scalarP);

    //     all_TEST7_scalarP = fopen("ipm-scalarP.pp", "a+");
    //     if (num_TEST7_scalarP == 10)
    //     {
    //       fprintf(all_TEST7_scalarP, "sumry: 1  %.2e  %.2e  %.2e   %5i %5i %5i %4i %4i %4i %4i  %.6f   %s\n",
    //             totalresidual, complem, projerr, iteration, n, m, nB, nN, nR, n-nB-nN-nR,
    //             total_time, strToken);
    //     }
    //     else
    //     {
    //       fprintf(all_TEST7_scalarP, "sumry: 1  %.2e  %.2e  %.2e   %5i %5i %5i %4i %4i %4i %4i  %.6f\n",
    //             totalresidual, complem, projerr, iteration, n, m, nB, nN, nR, n-nB-nN-nR,
    //             total_time);
    //     }
    //     fclose(all_TEST7_scalarP);

    //     num_TEST7_scalarP++;
    //   }
    // }

    // if (projerr_TEST7 > tols_TEST7_prj || total_time > max_time || alpha_primal <= 1e-13)
    // {
    //   while ( num_TEST7_prj < 11 )
    //   {
    //     sprintf(name_file_TEST7_prj, "ipm-prjerr-%02d.pp", num_TEST7_prj);
    //     file_TEST7_prj = fopen(name_file_TEST7_prj, "a+");
    //     fprintf(file_TEST7_prj, "sumry: 1  %.2e   %5i %5i %5i %4i %4i %4i %4i  %.6f   %s\n",
    //             projerr_TEST7, iteration, n, m, nB, nN, nR, n-nB-nN-nR,
    //             total_time, strToken);
    //     fclose(file_TEST7_prj);

    //     all_TEST7_prj = fopen("ipm-prjerr.pp", "a+");
    //     if (num_TEST7_prj == 10)
    //     {
    //       fprintf(all_TEST7_prj, "sumry: 1  %.2e  %.2e  %.2e   %5i %5i %5i %4i %4i %4i %4i  %.6f   %s\n",
    //             totalresidual, complem, projerr, iteration, n, m, nB, nN, nR, n-nB-nN-nR,
    //             total_time, strToken);
    //     }
    //     else
    //     {
    //       fprintf(all_TEST7_prj, "sumry: 1  %.2e  %.2e  %.2e   %5i %5i %5i %4i %4i %4i %4i  %.6f\n",
    //             totalresidual, complem, projerr, iteration, n, m, nB, nN, nR, n-nB-nN-nR,
    //             total_time);
    //     }
    //     fclose(all_TEST7_prj);

    //     num_TEST7_prj++;
    //   }
    // }

    if (totalresidual_Jor > tols_TEST7_JordanP || total_time > max_time || alpha_primal <= 1e-13)
    {
      while ( num_TEST7_JordanP < 9)
      {
        sprintf(name_file_TEST7_JordanP, "./res_20240706/update/external/NOSCAL/ipm-JordanP-%02d.pp", num_TEST7_JordanP);
        file_TEST7_JordanP = fopen(name_file_TEST7_JordanP, "a+");
        fprintf(file_TEST7_JordanP, "sumry: 1  %.2e   %5i %5i %5i %4i %4i %4i %4i  %.6f   %s\n",
                totalresidual_Jor, iteration, n, m, nB, nN, nR, n-nB-nN-nR,
                total_time, strToken);
        fclose(file_TEST7_JordanP);

        // all_TEST7_JordanP = fopen("ipm-JordanP.pp", "a+");
        // if (num_TEST7_JordanP == 10)
        // {
        //   fprintf(all_TEST7_JordanP, "sumry: 1  %.2e  %.2e  %.2e   %5i %5i %5i %4i %4i %4i %4i  %.6f   %s\n",
        //         totalresidual_Jor, complem, projerr, iteration, n, m, nB, nN, nR, n-nB-nN-nR,
        //         total_time, strToken);
        // }
        // else
        // {
        //   fprintf(all_TEST7_JordanP, "sumry: 1  %.2e  %.2e  %.2e   %5i %5i %5i %4i %4i %4i %4i  %.6f\n",
        //         totalresidual_Jor, complem, projerr, iteration, n, m, nB, nN, nR, n-nB-nN-nR,
        //         total_time);
        // }
        // fclose(all_TEST7_JordanP);

        num_TEST7_JordanP++;
      }
    }
  }





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

  options->solverData = (double *)malloc(sizeof(double));
  double *projerr_ptr = (double *)options->solverData;
  *projerr_ptr = projerr;


  if(H_tilde) H_tilde = NM_free(H_tilde);
  if(minus_H) minus_H = NM_free(minus_H);
  if(H) H = NM_free(H);
  if (Minv) Minv = NM_free(Minv);
  if(minus_M) minus_M = NM_free(minus_M);
  if(minus_Ht) minus_Ht = NM_free(minus_Ht);
  if(eye_nd) eye_nd = NM_free(eye_nd);
  if(Ht) Ht = NM_free(Ht);

  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
  {
    // fprintf(iterates, "];\n\n");
    fclose(iterates);
  }

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
  free(rhs_2);
  free(sr_rhs);
  free(sr_rhs_2);

  if (phiu) free(phiu);

  //  free(tmpsol);

  free(d_globalVelocity);
  free(d_velocity);
  free(d_reaction);
  free(diff_fixp_vec);


  *info = hasNotConverged;
}




void gfc3d_IPM_fixed(GlobalFrictionContactProblem* restrict problem, double* restrict reaction,
               double* restrict velocity, double* restrict globalVelocity,
               int* restrict info, SolverOptions* restrict options, double barr_param, const char* problem_name)
{
  unsigned int m = problem->M->size0;
  unsigned int nd = problem->H->size1;
  unsigned int d = problem->dimension;
  unsigned int n = problem->numberOfContacts;

  size_t no_m = 0, no_nd = 0;

  NumericsMatrix * M = NULL;
  NumericsMatrix * H_tilde = NULL;
  NumericsMatrix * eye_nd = NM_eye(nd);

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

  M = problem->M;
  H_tilde = NM_transpose(problem->H);

  NumericsMatrix *minus_M = NM_create( M->storageType, M->size0, M->size1);
  NM_copy(M, minus_M);
  NM_scal(-1.0, minus_M);

  // initialize solver if it is not set
  int internal_allocation=0;
  if(!options->dWork || (options->dWorkSize != (size_t)(m + nd + nd)))
  // if(!options->dWork)
  {
    gfc3d_IPM_init(problem, options);
    internal_allocation = 1;
    printf("\n\n gfc3d_IPM_fixed: allocation !!! \n\n");
  }

  Gfc3d_IPM_init_data * data = (Gfc3d_IPM_init_data *)options->solverData;
  NumericsMatrix *P_mu = data->P_mu->mat;
  NumericsMatrix *P_mu_inv = data->P_mu->inv_mat;

  double *w_tilde = problem->b;
  double *w = data->tmp_vault_nd[no_nd++];
  double *f = problem->q;

  NumericsMatrix *H = NM_multiply(P_mu, H_tilde);
  NM_gemv(1.0, P_mu, w_tilde, 0.0, w);

  double alpha_primal = data->internal_params->alpha_primal;
  double alpha_dual = data->internal_params->alpha_dual;

  // Copy the initial data into executive vars
  cblas_dcopy(nd, data->starting_point->reaction, 1, reaction, 1);
  cblas_dcopy(nd, data->starting_point->velocity, 1, velocity, 1);
  cblas_dcopy(m, data->starting_point->globalVelocity, 1, globalVelocity, 1);

  // // computation of the global velocity vector: v = M\(H'*r+f)
  // for (unsigned int  i = 0; i<m; i++) globalVelocity[i] = f[i];
  // NM_tgemv(1.0, H, reaction, 1.0, globalVelocity);
  // NM_Cholesky_solve(NM_preserve(M), globalVelocity, 1);

  double tol = options->dparam[SICONOS_DPARAM_TOL];
  unsigned int max_iter = options->iparam[SICONOS_IPARAM_MAX_ITER];

  double gmmp1 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1];
  double gmmp2 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2];

  int hasNotConverged = 1;
  unsigned int iteration = 0;
  double pinfeas = 1e300;
  double dinfeas = 1e300;
  double complem = 1e300;
  double dualgap = 1e300;
  double udotr = 1e300;
  double projerr = 1e300;
  double unitur = 1e300;
  double totalresidual = 1e300, totalresidual_mu = 1e300;
  double diff_fixp = 1e300, nub = 1e300;
  double *diff_fixp_vec = (double*)calloc(n,sizeof(double));

  double *primalConstraint = data->tmp_vault_nd[no_nd++];
  double *dualConstraint = data->tmp_vault_m[no_m++];
  double *complemConstraint = data->tmp_vault_nd[no_nd++];

  double gmm = gmmp1+gmmp2;
  double barr_param_a, e;
  double norm_f = cblas_dnrm2(m, f, 1);
  double norm_w = cblas_dnrm2(nd, w, 1);

  double *d_globalVelocity = data->tmp_vault_m[no_m++];
  double *d_velocity = data->tmp_vault_nd[no_nd++];
  double *d_reaction = data->tmp_vault_nd[no_nd++];
  double *s = (double*)calloc(n,sizeof(double));

  double *rhs = options->dWork;
  double *rhs_2 = (double*)calloc(m+2*nd, sizeof(double));

  double * Hvw = data->tmp_vault_nd[no_nd++];
  double * r_Qp_u = data->tmp_vault_nd[no_nd++];
  double * r_Qp_u_inv = data->tmp_vault_nd[no_nd++];

  double norm_q = cblas_dnrm2(m, problem->q, 1);
  double norm_b = cblas_dnrm2(nd, problem->b, 1);

  /* norm of the residuals of teh second linear system */
  double LS_norm_d = 0; // dual feaqsibility
  double LS_norm_c = 0; // complementarity

  NumericsMatrix *J = NULL;
  long J_nzmax;
  size_t H_nzmax = NM_nnz(H);

  numerics_printf_verbose(-1,"LS solution: 2x2 NT scaling with QpH, mu fixed, s fixed\n");

  numerics_printf_verbose(-1, "| it  |  dualgap | pinfeas | dinfeas |  <u, r> |  |s-ub| | complem | prj err | barpram |  2*n*mu |  alpha  |  |dv|   |  |du|   |  |dr|   | ls dual | ls comp |");
  numerics_printf_verbose(-1, "--------------------------------------------------------------------------------------------------------------------------------------------------------------");

  FILE * iterates;

  char *str = (char *) malloc(200);
  strcpy( str, problem_name );
  const char * separators = "/";
  char *strToken = strtok( str, separators );
  for(int i=0; i<5; i++)
  {
    if(strToken != NULL) strToken = strtok ( NULL, separators );
  }
  strToken = strtok ( strToken, "." );

  char matlab_name[100];

  /* writing data in a Matlab file */
  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
  {
    sprintf(matlab_name, "iterates_muF_sF.m");
    iterates = fopen(matlab_name, "a+");
    fprintf(iterates,"data(end+1).name = \"%s\";\n", strToken);
    fprintf(iterates,"data(end).val = [\n");
  }


  while(iteration < max_iter)
  {
    switch(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S])
    {
      case SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S_AFTER_SOLVING_LS:
      {
        if (iteration == 0)
        {
          for(unsigned int i = 0; i < nd; i+=d)
          {
            nub = cblas_dnrm2(2, velocity+i+1, 1);
            w[i] = w_tilde[i] + nub;
            s[i/d] = nub;
          }
        }

        primalResidual(velocity, H, globalVelocity, w, primalConstraint, &pinfeas, tol); // primalConstraint = u - H*v -w
        dualResidual(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol);  // dualConstraint =  M*v - H'*r - f
        dualgap = dualGap(M, f, w, globalVelocity, reaction, nd, m);
        JA_prod(velocity, reaction, nd, n, complemConstraint);
        for (unsigned int k = 0; k < nd; complemConstraint[k] -= 2.*barr_param, k += d);
        complem = cblas_dnrm2(nd, complemConstraint, 1);
        udotr = cblas_ddot(nd, velocity, 1, reaction, 1);
        projerr = projectionError(velocity, reaction, n, tol);
        for (unsigned int i = 0; i<nd; i+=d)
        {
          nub = cblas_dnrm2(2, velocity+i+1, 1);
          diff_fixp_vec[i/d] = s[i/d]-nub;
        }
        diff_fixp = cblas_dnrm2(n, diff_fixp_vec, 1);

        totalresidual_mu = fmax(fmax(pinfeas, dinfeas),complem);
        totalresidual = fmax(totalresidual_mu, diff_fixp);

        if (totalresidual <= tol)   // SUCCESS => need to quit the main loop
        {
          numerics_printf_verbose(-1, "| %3i | %8.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |",
                                iteration, dualgap, pinfeas, dinfeas, udotr, diff_fixp, complem, projerr, barr_param);

          for (int i = 0; i < n; i++)
          {
            unitur = cblas_ddot(3, velocity+3*i, 1, reaction+3*i, 1);
            if (unitur<0) printf("Cone %d: UR NEGATIF ui'ri = %9.2e\n", i, unitur);
          }

          hasNotConverged = 0;
        }

        else if (totalresidual_mu <= tol) // Finish solving ls (1 time) => update s after
        {
          numerics_printf_verbose(-1, "| %3i | %8.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |",
                                iteration, dualgap, pinfeas, dinfeas, udotr, diff_fixp, complem, projerr, barr_param);

          printf("\n");
          numerics_printf_verbose(-1, "| it  |  dualgap | pinfeas | dinfeas |  <u, r> |  |s-ub| | complem | prj err | barpram |  2*n*mu |  alpha  |  |dv|   |  |du|   |  |dr|   | ls prim | ls dual | ls comp |");

          // update s
          for(unsigned int i = 0; i < nd; i+=d)
          {
            nub = cblas_dnrm2(2, velocity+i+1, 1);
            w[i] = w_tilde[i] + nub;
            s[i/d] = nub;
          }
        }

        break;
      }

      case SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S_AT_EACH_ITE:
      {
        // update s
        for(unsigned int i = 0; i < nd; i+=d)
        {
          nub = cblas_dnrm2(2, velocity+i+1, 1);
          w[i] = w_tilde[i] + nub;
          s[i/d] = nub;
        }
        diff_fixp = 0.;

        primalResidual(velocity, H, globalVelocity, w, primalConstraint, &pinfeas, tol); // primalConstraint = u - H*v -w
        dualResidual(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol);  // dualConstraint =  M*v - H'*r - f
        dualgap = dualGap(M, f, w, globalVelocity, reaction, nd, m);
        JA_prod(velocity, reaction, nd, n, complemConstraint);
        for (unsigned int k = 0; k < nd; complemConstraint[k] -= 2.*barr_param, k += d);
        complem = cblas_dnrm2(nd, complemConstraint, 1);
        udotr = cblas_ddot(nd, velocity, 1, reaction, 1);
        projerr = projectionError(velocity, reaction, n, tol);

        totalresidual = fmax(fmax(pinfeas, dinfeas),complem);

        if (totalresidual <= tol)   // SUCCESS => need to quit the main loop
        {
          numerics_printf_verbose(-1, "| %3i | %8.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |",
                                iteration, dualgap, pinfeas, dinfeas, udotr, diff_fixp, complem, projerr, barr_param);

          for (int i = 0; i < n; i++)
          {
            unitur = cblas_ddot(3, velocity+3*i, 1, reaction+3*i, 1);
            if (unitur<0) printf("Cone %d: UR NEGATIF ui'ri = %9.2e\n", i, unitur);
          }

          hasNotConverged = 0;
        }
        break;
      }

      default:
      {
        printf("gfc3d_IPM_fixed: update_s strategy is unknown. \n\n");
      }

    }


    if (hasNotConverged == 0)
    {
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
      {
        fprintf(iterates,"%d %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e;\n",
                iteration, pinfeas, dinfeas, udotr, complem, projerr, barr_param, alpha_primal,
                fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
                fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
                fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]),
                LS_norm_d, LS_norm_c);
      }
      break;
    }


    J = NM_create(NM_SPARSE, m + nd, m + nd);
    J_nzmax = (d * d) * (m / d) + H_nzmax + H_nzmax + nd;
    NM_triplet_alloc(J, J_nzmax);
    J->matrix2->origin = NSM_TRIPLET;

    NumericsMatrix * QpH = QNTpH(velocity, reaction, H, nd, n);
    NumericsMatrix * QpHt = NM_transpose(QpH);

    NM_insert(J, minus_M, 0, 0);
    NM_insert(J, QpH, m, 0);
    NM_insert(J,QpHt, 0, m);
    NM_insert(J, eye_nd, m, m);

    if(QpH) QpH = NM_free(QpH);
    if(QpHt) QpHt = NM_free(QpHt);

    if (NM_isnan(J))
    {
      numerics_printf_verbose(0, "The Jacobian matrix J contains NaN");
      hasNotConverged = 2;
      if (J) J = NM_free(J);
      break;
    }

    QNTpz(velocity, reaction, velocity, nd, n, r_Qp_u);  // r_Qp_u <- u_hat
    JA_inv(r_Qp_u, nd, n, r_Qp_u_inv);                   // r_Qp_u_inv <- u_hat_inv
    cblas_dscal(nd, -2.*barr_param, r_Qp_u_inv, 1);      // r_Qp_u_inv <= -2 * mu * u_hat_inv

    cblas_dcopy(m, dualConstraint, 1, rhs, 1);
    NM_gemv(1.0, H, globalVelocity, 0.0, Hvw);
    cblas_daxpy(nd, 1.0, w, 1, Hvw, 1);                  // Hvw will be used for the computation of r_du
    QNTpz(velocity, reaction, Hvw, nd, n, rhs+m);
    cblas_daxpy(nd, 1.0, r_Qp_u_inv, 1, rhs+m, 1);
    cblas_dscal(nd, -1.0, rhs+m, 1);

    cblas_dcopy(m+nd, rhs, 1, rhs_2, 1);

    NSM_linearSolverParams(J)->solver = NSM_HSL;
    NM_LDLT_solve(J, rhs, 1);

    cblas_dcopy(m, rhs, 1, d_globalVelocity, 1);
    QNTpz(velocity, reaction, rhs+m, nd, n, d_reaction);
    NM_gemv(1.0, H, d_globalVelocity, 0.0, d_velocity);  // d_velocity <- H*dv
    cblas_daxpy(nd, 1.0, Hvw, 1, d_velocity, 1);         // d_velocity <- H*v + w + H*dv
    cblas_daxpy(nd, -1.0, velocity, 1, d_velocity, 1);   // d_velocity <- H*(v + dv) + w  - u

    double * rhs_tmp = (double*)calloc(m+nd,sizeof(double));
    cblas_dcopy(m+nd, rhs_2, 1, rhs_tmp, 1);
    NM_gemv(1.0, J, rhs, -1.0, rhs_tmp);
    LS_norm_d = cblas_dnrm2(m, rhs_tmp, 1);
    LS_norm_c = cblas_dnrm2(nd, rhs_tmp+m, 1);

    free(rhs_tmp);

    if (J) J = NM_free(J);

    /* computing the affine step-length */
    alpha_primal = getStepLength(velocity, d_velocity, nd, n, gmm);
    alpha_dual = getStepLength(reaction, d_reaction, nd, n, gmm);

    if (alpha_primal < alpha_dual) alpha_dual = alpha_primal;
    else  alpha_primal = alpha_dual;

    /* updating the gamma parameter used to compute the step-length */
    gmm = gmmp1 + gmmp2 * fmin(alpha_primal, alpha_dual);

    numerics_printf_verbose(-1, "| %3i | %8.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |",
                            iteration, dualgap, pinfeas, dinfeas, udotr, diff_fixp, complem, projerr, barr_param, 2*n*barr_param, alpha_primal,
                            fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
                            fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
                            fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]),
                            LS_norm_d, LS_norm_c);

    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
    {
      fprintf(iterates,"%d %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e;\n",
              iteration, pinfeas, dinfeas, udotr, complem, projerr, barr_param, alpha_primal,
              fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
              fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
              fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]),
              LS_norm_d, LS_norm_c);
    }

    /* ----- Update variables ----- */

    if (NV_isnan(d_globalVelocity, m) | NV_isnan(d_velocity, nd) | NV_isnan(d_reaction, nd))
    {
      hasNotConverged = 2;
      break;
    }
    cblas_daxpy(m, alpha_primal, d_globalVelocity, 1, globalVelocity, 1);
    cblas_daxpy(nd, alpha_primal, d_velocity, 1, velocity, 1);
    cblas_daxpy(nd, alpha_dual, d_reaction, 1, reaction, 1);

    iteration++;
  } // while loop

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
    printf("Strict complementarity satisfied: %4i / %4i  %9.2e %9.2e  %4i %4i %4i\n", nsc, n, ns, ns/cblas_dnrm2(3, s, 1), nB, nN, nR);


  options->dparam[SICONOS_DPARAM_RESIDU] = totalresidual;
  options->iparam[SICONOS_IPARAM_ITER_DONE] = iteration;
  cblas_dcopy(nd, velocity, 1, data->starting_point->velocity, 1);
  cblas_dcopy(nd, reaction, 1, data->starting_point->reaction, 1);
  cblas_dcopy(m, globalVelocity, 1, data->starting_point->globalVelocity, 1);


  if(internal_allocation)
  {
    gfc3d_IPM_free(problem,options);
  }

  options->solverData = (double *)malloc(sizeof(double));
  double *projerr_ptr = (double *)options->solverData;
  *projerr_ptr = projerr;


  if(H_tilde) H_tilde = NM_free(H_tilde);
  if(H) H = NM_free(H);
  if(minus_M) minus_M = NM_free(minus_M);
  if(eye_nd) eye_nd = NM_free(eye_nd);
  free(diff_fixp_vec);

  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
  {
    fprintf(iterates, "];\n\n");
    fclose(iterates);
  }
  *info = hasNotConverged;
}


void gfc3d_ipm_set_default(SolverOptions* options)
{
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO] = SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_NO;

  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] = 0;
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] = SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S_AFTER_SOLVING_LS;
  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] = SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S_AT_EACH_ITE;
  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] = SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S_AFTER_SOLVING_LS_AND_CONVERGE_TO_MU;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] = SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL;
  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] = SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_QP2;
  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] = SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QPH;
  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] = SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_TEST6;
  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] = SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QPH_TEST;


  //options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING] = SICONOS_FRICTION_3D_RESCALING_BALANCING_MHHT;
  options->iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING] = SICONOS_FRICTION_3D_RESCALING_NO;

  //options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING_METHOD] = SICONOS_FRICTION_3D_IPM_NESTEROV_TODD_SCALING_WITH_QP;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE] = 0;

  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] = SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_NO;
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] = SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES;

  options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000000;
  options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1] = 1e-10;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2] = 3.;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3] = 1.;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1] = 0.9;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2] = 0.09; //0.095

}
