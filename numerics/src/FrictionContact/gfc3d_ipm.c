/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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
#include "CSparseMatrix_internal.h"
#include "NumericsSparseMatrix.h"
#include "NumericsMatrix.h"


/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "debug.h"


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


/* ------------------------- Helper functions declaration ------------------------------ */

/* Returns the step length for variables update in IPM [1, p. 29] */
/* \param x is the initial point to update. */
/* \param dx is the Newton step. */
/* \param vecSize is the size of the vectors x and dx. */
/* \param varsCount is the count of variables concatenated into vector x. */
/* \param gamma is the safety parameter. */
/* \return scalar, the step length */

/* \cite 1. K.C. Toh, R.H. Tutuncu, M.J. Todd, */
/*          On the implementation and usage of SDPT3 - a Matlab software package */
/*          for semidefinite-quadratic-linear programming, version 4.0 */
/*          Draft, 17 July 2006 */
static double getNewtonStepLength(const double * const x, const double * const dx,
				  const unsigned int vecSize, const unsigned int varsCount, const double gamma);

/* Returns the maximum step-length to the boundary reduced by a factor gamma. Uses long double. */
static double getStepLength(const double * const x, const double * const dx, const unsigned int vecSize,
                                  const unsigned int varsCount, const double gamma);

/**
 * Returns the primal constraint vector for global fricprob ( velocity - H @ globalVelocity - w )
 * \param velocity is the vector of relative velocities.
 * \param H is the constraint matrix.
 * \param globalVelocity is the vector of generalized velocities.
 * \param w is the constraint vector.
 * \param out is the result velocity - H @ globalVelocity - w vector.
 */
static void primalResidualVector(const double * velocity, NumericsMatrix * H,
                                 const double * globalVelocity, const double * w, double * out);


/**
 * Returns the dual constraint vector for global fricprob ( M @ globalVelocity + f - H @ reaction )
 * \param M is the mass matrix.
 * \param globalVelocity is the vector of generalized velocities.
 * \param H is the constraint matrix.
 * \param reaction is the vector of reaction forces at each contact point.
 * \param f is the constraint vector (vector of internal and external forces).
 * \param out os the result M @ globalVelocity + f - H @ reaction vector.
 */
static void dualResidualVector(NumericsMatrix * M, const double * globalVelocity,
                               NumericsMatrix * H, const double * reaction, const double * f,
                               double * out);

/**
 * Returns the Inf-norm of primal residual vecor ( velocity - H @ globalVelocity - w )
 * \param velocity is the vector of relative velocities.
 * \param H is the constraint matrix.
 * \param globalVelocity is the vector of generalized velocities.
 * \param w is the constraint vector.
 */
static double primalResidualNorm(const double * velocity, NumericsMatrix * H,
                                 const double * globalVelocity, const double * w);

/**
 * Returns the Inf-norm of the dual residual vector ( M @ globalVelocity + f - H @ reaction )
 * \param M is the mass matrix.
 * \param globalVelocity is the vector of generalized velocities.
 * \param H is the constraint matrix.
 * \param reaction is the vector of reaction forces at each contact point.
 * \param f is the constraint vector (vector of internal and external forces).
 */
static double dualResidualNorm(NumericsMatrix * M, const double * globalVelocity,
                               NumericsMatrix * H, const double * reaction, const double * f);

/**
 * Returns the Inf-norm of the cemplementarity residual vector ( velocity o reaction )
 * \param velocity is the vector of relative velocities at each contact point.
 * \param reaction is the vector of reaction forces at each contact point.
 * \param vecSize is the size of reaction and velocity vector.
 * \param varsCount is a count of variables concatenated into vectors reaction and velocity.
 */
static double complemResidualNorm(const double * const velocity, const double * const reaction,
                                  const unsigned int vecSize, const unsigned int varsCount);

/* Returns the 2-norm of the complementarity residual vector = 2-norm of the Jordan product (Qp*velocity) o (Qp_inv * reaction)  */
static double complemResidualNorm_p(const double * const velocity, const double * const reaction,
                                    const unsigned int vecSize, const unsigned int varsCount);

/* computation of the duality gap  */
static double dualGap(NumericsMatrix * M, const double * f, const double * w, const double * globalVelocity, const double * reaction, const unsigned int nd, const unsigned int m);

/* Return the 2-norm of the difference between two vectors */ 
double norm2VecDiff (const double * vec1, const double * vec2, const unsigned int vecSize);

/* PA: Return the product Q_sqrt(x)*y */
void Qx05y(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * out);

/* PA: Return the product Q_inv_sqrt(x)*y */
void Qx50y(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * out);

/* PA: Return J_sqrt(x) */
void Jsqrt(const double * const x, const unsigned int vecSize, const size_t varsCount, double * out);

/* PA: Return J_sqrtinv(x) */
void Jsqrtinv(const double * const x, const unsigned int vecSize, const size_t varsCount, double * out);

/* PA: Jordan algebra, returns inv(x) */
void Jinv(const double * const x, const unsigned int vecSize, const size_t varsCount, double * out);

/* PA: Return the Nesterov-Todd vector */
void Nesterov_Todd_vector(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * p);

/* Computation of Qx*y by means of the formula 2*(x'*y)*x - det(x)*R*y */
void Qxy(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * z);

/* Returns the product Q_{p}*z where p is the NT vector related to the pair (x,y) */ 
void QNTpz(const double * const x, const double * const y,const double * const z, const unsigned int vecSize, const size_t varsCount, double * out);

/* Returns the product Q_{p^{-1}}*z where p is the NT vector related to the pair (x,y) */ 
void QNTpinvz(const double * const x, const double * const y,const double * const z, const unsigned int vecSize, const size_t varsCount, double * out);

/* returns the Jordan product x^{-1} o y by using the formula x^{-1} = R*x/det(x), where R is the reflection matrix */
void Jxinvprody(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * out);

/* Returns the 2-norm of a vector - uses long double - based on blas_dnrm2 */
static long double dnrm2l(const unsigned int n, const double * x);

/* ------------------------- Helper functions implementation ------------------------------ */

/* Returns the 2-norm of a vector - uses long double - based on blas_dnrm2 */
static long double dnrm2l(const unsigned int n, const double * x)
{
  long double norm, scale, ssq, absxi, quo;

  if (n < 1)
    norm = 0.0;
  else if (n == 1)
    norm = fabs(x[0]);
  else
    {
      scale = 0.0;
      ssq = 1.0;
      for (size_t i = 0; i < n; i++)
	{
	  if (x[i] != 0)
		{
		  absxi = fabsl(x[i]);
		  if (scale < absxi)
		    {
		      quo = scale/absxi;
		      ssq = 1.0 + ssq * (quo * quo);
		      scale = absxi;
		    }
		  else
		    {
		      quo = absxi/scale;
		      ssq = ssq + (quo * quo);
		    }		  
		}
	  norm = scale * sqrtl(ssq);
	}
    }
  return norm;
}

/* Returns the step length for variables update in IPM */
static double getNewtonStepLength(const double * const x, const double * const dx, const unsigned int vecSize,
                                  const unsigned int varsCount, const double gamma)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  double * alpha_list = (double*)calloc(varsCount, sizeof(double));

  unsigned int pos;
  double ai, bi, ci, di, alpha, min_alpha;
  double  *xi2, *dxi2, *xi_dxi;

  const double *dxi, *xi;

  dxi2 = (double*)calloc(dimension, sizeof(double));
  xi2 = (double*)calloc(dimension, sizeof(double));
  xi_dxi = (double*)calloc(dimension, sizeof(double));

  for(unsigned int i = 0; i < varsCount; ++i)
  {
    pos = i * dimension;
    xi = x + pos;
    dxi = dx + pos;

    NV_power2(dxi, dimension, dxi2);
    ai = dxi2[0] - NV_reduce((dxi2 + 1), dimension - 1);

    NV_prod(xi, dxi, dimension, xi_dxi);
    bi = xi_dxi[0] - NV_reduce((xi_dxi + 1), dimension - 1);
    //   bi = gamma*bi;

    NV_power2(xi, dimension, xi2);
    ci = xi2[0] - NV_reduce((xi2 + 1), dimension - 1);
    //    ci = gamma*gamma*ci;
    
    di = bi * bi - ai * ci;

    if(ai < 0 || (bi < 0 && ai < (bi * bi) / ci))
      alpha = ((-bi - sqrt(di)) / ai);
    else if((fabs(ai) < DBL_EPSILON) && (bi < 0))
      alpha = (-ci / (2 * bi));
    else
      alpha = DBL_MAX;
    //NV_display(xi2, dimension);
    //printf("**************** %3i ai = %9.2e b = %9.2e ci = %9.2e alpha = %9.2e\n",i, ai, bi, ci, alpha);

    if(fabs(alpha) < DBL_EPSILON)
      alpha = 0.0;

    alpha_list[i] = alpha;
  }

  min_alpha = NV_min(alpha_list, varsCount);

  free(xi2);
  free(dxi2);
  free(xi_dxi);
  free(alpha_list);

  return fmin(1.0, gamma * min_alpha);
  //return fmin(1.0, min_alpha);
}

/* Returns the maximum step-length to the boundary reduced by a factor gamma. Uses long double. */
static double getStepLength(const double * const x, const double * const dx, const unsigned int vecSize,
                                  const unsigned int varsCount, const double gamma)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  unsigned int pos;
  long double aL, bL, cL, dL, alphaL;
  double min_alpha;

  min_alpha = 1.0;
  
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
  return gamma*min_alpha;
}

/* Returns the primal constraint vector for global fricprob ( velocity - H @ globalVelocity - w ) */
static void primalResidualVector(const double * velocity, NumericsMatrix * H,
                                 const double * globalVelocity, const double * w, double * out)
{
  double nd = H->size0;

  /* The memory for the result vectors should be allocated using calloc
   * since H is a sparse matrix. In other case the behaviour will be undefined.*/
  //  double *Hv = (double*)calloc(nd, sizeof(double));
  //double *u_minus_Hv = (double*)calloc(nd, sizeof(double));

  // Hv = H @ globalVelocity
  NM_gemv(-1.0, H, globalVelocity, 0.0, out);
  cblas_daxpy(nd, 1.0, velocity, 1, out, 1);
  cblas_daxpy(nd, -1.0, w, 1, out, 1);
  // u_minus_Hv = velocity - H @ globalVelocity
  //NV_sub(velocity, Hv, nd, u_minus_Hv);

  // out = velocity - H @ globalVelocity - w
  // NV_sub(u_minus_Hv, w, nd, out);

  // free allocated memory
  //free(Hv);
  //free(u_minus_Hv);
}

/* Returns the dual constraint vector for global fricprob ( M*globalVelocity - f - H'*reaction ) */
static void dualResidualVector(NumericsMatrix * M, const double * globalVelocity,
                               NumericsMatrix * H, const double * reaction, const double * f,
                               double * out)
{
  double m = H->size1;
  double *HTr = (double*)calloc(m, sizeof(double));

  NM_gemv(1.0, M, globalVelocity, 0.0, out);
  cblas_daxpy(m, -1.0, f, 1, out, 1);
  NM_tgemv(1.0, H, reaction, 0.0, HTr);
  cblas_daxpy(m, -1.0, HTr, 1, out, 1);
  free(HTr);
}

/* Returns the 2-norm of primal residual vector = | H * globalVelocity + w - velocity |_2 / (1 + |w|_inf) */
static double primalResidualNorm(const double * velocity, NumericsMatrix * H,
                                 const double * globalVelocity, const double * w)
{
  double * resid = (double*)calloc(H->size0, sizeof(double));
  primalResidualVector(velocity, H, globalVelocity, w, resid);
  double norm_2 = cblas_dnrm2(H->size0, resid, 1);
  free(resid);
  return norm_2 / (1 + NV_norm_inf(w, H->size0));
}

/* Returns the 2-norm of the dual residual vector  = | M * globalVelocity - H * reaction + f |_2 / (1 + |f|_inf  */
static double dualResidualNorm(NumericsMatrix * M, const double * globalVelocity,
                               NumericsMatrix * H, const double * reaction, const double * f)
{
  double * resid = (double*)calloc(H->size1, sizeof(double));
  dualResidualVector(M, globalVelocity, H, reaction, f, resid);
  double norm_2 = cblas_dnrm2(H->size1, resid, 1);
  free(resid);
  return norm_2 / (1 + NV_norm_inf(f, H->size1));
}

/* Returns the 2-norm of the complementarity residual vector = 2-norm of the Jordan product velocity o reaction  */
static double complemResidualNorm(const double * const velocity, const double * const reaction,
                                  const unsigned int vecSize, const unsigned int varsCount)
{
  double * resid = (double*)calloc(vecSize, sizeof(double));
  JA_prod(velocity, reaction, vecSize, varsCount, resid);
  double norm2 = cblas_dnrm2(vecSize, resid, 1);
  free(resid);
  return norm2;
}

/* Returns the 2-norm of the complementarity residual vector = 2-norm of the Jordan product (Qp*velocity) o (Qp_inv * reaction)  */
static double complemResidualNorm_p(const double * const velocity, const double * const reaction,
                                    const unsigned int vecSize, const unsigned int varsCount)
{
  double * resid = (double*)calloc(vecSize, sizeof(double));
  double * u_p = (double*)calloc(vecSize, sizeof(double));
  double * r_p = (double*)calloc(vecSize, sizeof(double));
  //double * p_inv = (double*)calloc(vecSize, sizeof(double));
  double * a = (double*)calloc(vecSize, sizeof(double));
  double * b = (double*)calloc(vecSize, sizeof(double));

  Qx05y(velocity, reaction, vecSize, varsCount,a);
  Jsqrtinv(a, vecSize, varsCount, b);
  Qx05y(velocity, b, vecSize, varsCount, a);
  
  Qx50y(a, velocity, vecSize, varsCount, u_p);
  Qx05y(a, reaction, vecSize, varsCount, r_p);
  JA_prod(u_p, r_p, vecSize, varsCount, resid);

  double norm2 = cblas_dnrm2(vecSize, resid, 1);
  
  /* Qxy(p, velocity, vecSize, varsCount, u_p); */
  /* JA_inv(p, vecSize, varsCount, p_inv); */
  /* Qxy(p_inv, reaction, vecSize, varsCount, r_p); */
  /* JA_prod(u_p, r_p, vecSize, varsCount, resid); */
  
  /* norm2 = cblas_dnrm2(vecSize, resid, 1); */
  /* printf("complemnt-2 = %.15e\n",norm2); */
  //norm2 = norm2/ cblas_dnrm2(vecSize, velocity, 1);
  //norm2 = norm2 / cblas_dnrm2(vecSize, reaction, 1);
  free(resid);
  free(u_p);
  free(r_p);
  //free(p_inv);
  free(a);
  free(b);
  
  return norm2;
}

/* computation of the duality gap  */
static double dualGap(NumericsMatrix * M, const double * f, const double * w, const double * globalVelocity, const double * reaction, const unsigned int nd, const unsigned int m)
{
  double * Mv = (double*)calloc(m, sizeof(double));
  double vMv, pval, dval; 

  NM_gemv(0.5, M, globalVelocity, 0.0, Mv);
  vMv = cblas_ddot(m, globalVelocity, 1, Mv, 1);
  free(Mv);
  pval = vMv - cblas_ddot(m, f, 1, globalVelocity, 1);
  dval = -vMv - cblas_ddot(nd, w, 1, reaction, 1);
  return (pval - dval)/ (1 + fabs(pval) + fabs(dval));
}

static void setErrorArray(double * error, const double pinfeas, const double dinfeas,
                          const double dualgap, const double complem)
{
  error[0] = pinfeas;
  error[1] = dinfeas;
  error[2] = dualgap;
  error[3] = complem;
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

/* PA: Return the product Q_sqrt(x)*y */
void Qx05y(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * out)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  long double l1, l2, c1y, c2y, nxb, fx1, fx2, dx; 
  size_t j;
  long double *xb = (long double*)calloc(dimension-1, sizeof(long double));

  for (int i = 0; i < dimension - 1; xb[i] = 1/sqrtl(dimension-1), i++);
  
  for(size_t i = 0; i < varsCount; i++)
    {
      j = i*dimension;
      nxb = dnrm2l(dimension-1, x+j+1);
      if (nxb > 0)
	for (int k = 0; k < dimension-1; xb[k] = x[j+1+k]/nxb, k++);
      l1 = x[j]+nxb;
      l2 = x[j]-nxb;
      dx = sqrtl(l1*l2);
      c1y = y[j];
      for (int k = 0; k < dimension-1; c1y += xb[k]*y[j+1+k], k++);
      c2y = 2*y[j] - c1y;
      fx1 = (l1*c1y + dx*c2y)/2;
      fx2 = (dx*c1y + l2*c2y)/2;
      out[j] = fx1 + fx2 - dx*y[j];
      for (int k = 0; k < dimension-1; out[j+k+1] = fx1*xb[k] - fx2*xb[k] + dx*y[j+k+1], k++);
    }
  free(xb);
}

/* PA: Return the product Q_inv_sqrt(x)*y */
void Qx50y(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * out)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  long double l1, l2, c1y, c2y, nxb, fx1, fx2, dx; 
  size_t j;
  long double *xb = (long double*)calloc(dimension-1, sizeof(long double));

  for (int i = 0; i < dimension - 1; xb[i] = 1/sqrtl(dimension-1), i++);
  
  for(size_t i = 0; i < varsCount; i++)
    {
      j = i*dimension;
      nxb = dnrm2l(dimension-1, x+j+1);
      if (nxb > 0)
	for (int k = 0; k < dimension-1; xb[k] = x[j+1+k]/nxb, k++);
      
      l1 = x[j]+nxb;
      l2 = x[j]-nxb;
      dx = 1/sqrtl(l1*l2);
      c1y = y[j];
      for (int k = 0; k < dimension-1; c1y += xb[k]*y[j+1+k], k++);
      c2y = 2*y[j] - c1y;
      fx1 = (c1y/l1 + dx*c2y)/2;
      fx2 = (dx*c1y + c2y/l2)/2;
      out[j] = fx1 + fx2 - dx*y[j];
      for (int k = 0; k < dimension-1; out[j+k+1] = fx1*xb[k] - fx2*xb[k] + dx*y[j+k+1], k++);
    }
  free(xb);
}
/* void Qx50y(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * out) */
/* { */
/*   unsigned int dimension = (int)(vecSize / varsCount); */
/*   double l1, l2, c1y, c2y, normx, fx1, fx2, dx;  */
/*   size_t j; */
/*   double *xn = (double*)calloc(dimension-1, sizeof(double)); */

/*   for (int i = 0; i < dimension - 1; xn[i] = 1/sqrt(dimension-1), i++); */

/*   for(size_t i = 0; i < varsCount; i++) */
/*     { */
/*       j = i*dimension; */
/*       normx = cblas_dnrm2(dimension-1, x+j+1, 1); */
/*       if (normx > 0) */
/* 	{ */
/* 	  cblas_dcopy(dimension-1, x+j+1, 1, xn, 1); */
/* 	  cblas_dscal(dimension-1, 1.0/normx, xn, 1); */
/* 	} */
/*       cblas_dcopy(dimension-1, x+j+1, 1, xn, 1); */
/*       cblas_dscal(dimension-1, 1.0/normx, xn, 1); */
/*       l1 = 1/(x[j]+normx); */
/*       l2 = 1/(x[j]-normx); */
/*       dx = sqrt(l1*l2); */
/*       c1y = y[j] + cblas_ddot(dimension-1, xn, 1, y+j+1, 1); */
/*       c2y = 2*y[j] - c1y; */
/*       fx1 = (l1*c1y + dx*c2y)/2; */
/*       fx2 = (dx*c1y + l2*c2y)/2; */
/*       out[j] = fx1 + fx2 - dx*y[j]; */
/*       for (int k = 0; k < dimension-1; out[j+k+1] = fx1*xn[k] - fx2*xn[k] + dx*y[j+k+1], k++); */
/*     } */
/*   free(xn); */
/* } */

/* PA: Jordan algebra, returns inv(x) */
void Jinv(const double * const x, const unsigned int vecSize, const size_t varsCount, double * out)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  long double l1, l2, normx; 
  size_t j;

  for(size_t i = 0; i < varsCount; i++)
    {
      j = i*dimension;
      normx = dnrm2l(dimension-1, x+j+1);
      l1 = 1/(x[j]+normx)/2;
      l2 = 1/(x[j]-normx)/2;
      out[j] = l1+l2;
      for (int k = 1; k < dimension; out[j+k] = l1*(x[j+k]/normx) - l2*(x[j+k]/normx), k++);
    }
}

/* /\* PA: Return J_sqrt(x) *\/ */
/* void Jsqrt(const double * const x, const unsigned int vecSize, const size_t varsCount, double * out) */
/* { */
/*   unsigned int dimension = (int)(vecSize / varsCount); */
/*   double l1, l2, normx; */
/*   size_t j; */

/*   for(size_t i = 0; i < varsCount; i++) */
/*     { */
/*       j = i*dimension; */
/*       normx = cblas_dnrm2(dimension-1, x+j+1, 1); */
/*       l1 = sqrt(x[j]+normx)/2; */
/*       l2 = sqrt(x[j]-normx)/2; */
/*       out[j] = l1+l2; */
/*       for (int k = 1; k < dimension; out[j+k] = l1*(x[j+k]/normx) - l2*(x[j+k]/normx), k++); */
/*     } */
/* } */
/* /\* PA: Return J_sqrtinv(x) *\/ */
/* void Jsqrtinv(const double * const x, const unsigned int vecSize, const size_t varsCount, double * out) */
/* { */
/*   unsigned int dimension = (int)(vecSize / varsCount); */
/*   double l1, l2, normx; */
/*   size_t j; */

/*   for(size_t i = 0; i < varsCount; i++) */
/*     { */
/*       j = i*dimension; */
/*       normx = cblas_dnrm2(dimension-1, x+j+1, 1); */
/*       l1 = 1/sqrt(x[j]+normx)/2; */
/*       l2 = 1/sqrt(x[j]-normx)/2; */
/*       out[j] = l1+l2; */
/*       for (int k = 1; k < dimension; out[j+k] = l1*(x[j+k]/normx) - l2*(x[j+k]/normx), k++); */
/*     } */
/* } */

/* PA: Return J_sqrt(x) */
void Jsqrt(const double * const x, const unsigned int vecSize, const size_t varsCount, double * out)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  long double * xx = (long double*)calloc(dimension, sizeof(long double));
  long double l1, l2, normx;
  size_t j;

  for(size_t i = 0; i < varsCount; i++)
    {
      j = i*dimension;
      for (int k = 0; k < dimension; xx[k] = x[j+k], k++);
      normx = 0;
      for (int k = 1; k < dimension; normx += xx[k]*xx[k], k++);
      normx = sqrtl(normx);
      l1 = sqrtl(xx[0]+normx)/2;
      l2 = sqrtl(xx[0]-normx)/2;
      out[j] = l1+l2;
      for (int k = 1; k < dimension; out[j+k] = l1*(xx[k]/normx) - l2*(xx[k]/normx), k++);
    }
  free(xx);
}
/* PA: Return J_sqrtinv(x) */
void Jsqrtinv(const double * const x, const unsigned int vecSize, const size_t varsCount, double * out)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  long double l1, l2, normx;
  size_t j;

  for(size_t i = 0; i < varsCount; i++)
    {
      j = i*dimension;
      normx = sqrtl(normx);
      normx = dnrm2l(dimension-1, x+j+1);
      l1 = 1/sqrtl(x[j]+normx)/2;
      l2 = 1/sqrtl(x[j]-normx)/2;
      out[j] = l1+l2;
      for (int k = 1; k < dimension; out[j+k] = l1*(x[j+k]/normx) - l2*(x[j+k]/normx), k++);
    }
}

/* PA: Return the Nesterov-Todd vector */
void Nesterov_Todd_vector(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * p)
{
  double * a = (double*)calloc(vecSize, sizeof(double));
  double * b = (double*)calloc(vecSize, sizeof(double));
  
  Qx05y(x, y, vecSize, varsCount,a);
  Jsqrtinv(a, vecSize, varsCount, b);
  Qx05y(x, b, vecSize, varsCount, a);
  Jsqrtinv(a, vecSize, varsCount, p);

  free(a);
  free(b);
}

/* PA: Return the Nesterov-Todd vector by means of the second formula */
void Nesterov_Todd_vector_b(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * p)
{
  double * a = (double*)calloc(vecSize, sizeof(double));
  double * b = (double*)calloc(vecSize, sizeof(double));
  
  Qx05y(y, x, vecSize, varsCount,a);
  Jsqrt(a, vecSize, varsCount, b);
  Qx50y(y, b, vecSize, varsCount, a);
  Jsqrtinv(a, vecSize, varsCount, p);

  free(a);
  free(b);
}

/* Computation of Qx*y by means of the formula 2*(x'*y)*x - det(x)*R*y */
void Qxy(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * z)

{
  unsigned int dimension = (int)(vecSize / varsCount);
  size_t j;
  double xy;
  long double dx, nxb;
  
  for(size_t i = 0; i < varsCount; i++)
    {
      j = i*dimension;
      xy = 2*cblas_ddot(dimension, x+j, 1, y+j, 1);
      /* dx = x[j]*x[j]; */
      /* for (int k = 1; k < dimension; dx -= x[j+k]*x[j+k], k++); */
      nxb = dnrm2l(dimension-1, x+j+1);
      dx = (x[j] + nxb)*(x[j] - nxb);
      z[j] = xy*x[j] - dx*y[j];
      for (int k = 1; k < dimension; z[j+k] = xy*x[j+k] + dx*y[j+k], k++);
    }
}

/* Returns the product Q_{p}*z where p is the NT vector related to the pair (x,y) */ 
void QNTpz(const double * const x, const double * const y,const double * const z, const unsigned int vecSize, const size_t varsCount, double * out)
{
  double * a = (double*)calloc(vecSize, sizeof(double));
  double * b = (double*)calloc(vecSize, sizeof(double));
  
  Qx05y(x, y, vecSize, varsCount,a);
  Jsqrtinv(a, vecSize, varsCount, b);
  Qx05y(x, b, vecSize, varsCount, a);
  Qx50y(a, z, vecSize, varsCount, out);

  free(a);
  free(b);
}

/* Returns the product Q_{p^{-1}}*z where p is the NT vector related to the pair (x,y) */ 
void QNTpinvz(const double * const x, const double * const y,const double * const z, const unsigned int vecSize, const size_t varsCount, double * out)
{
  double * a = (double*)calloc(vecSize, sizeof(double));
  double * b = (double*)calloc(vecSize, sizeof(double));
  
  Qx05y(x, y, vecSize, varsCount,a);
  Jsqrtinv(a, vecSize, varsCount, b);
  Qx05y(x, b, vecSize, varsCount, a);
  Qx05y(a, z, vecSize, varsCount, out);

  free(a);
  free(b);
}

/* returns the Jordan product x^{-1} o y by using the formula x^{-1} = R*x/det(x), where R is the reflection matrix */
void Jxinvprody(const double * const x, const double * const y, const unsigned int vecSize, const size_t varsCount, double * out)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  long double nxb, detx, tmp; 
  size_t j;

  for(size_t i = 0; i < varsCount; i++)
    {
      j = i*dimension;
      
      nxb = dnrm2l(dimension-1, x+j+1);
      detx = (x[j] + nxb) * (x[j] - nxb);
      
      tmp = x[j]*y[j];
      for (int k = 1; k < dimension; tmp -= x[j+k]*y[j+k], k++);
      out[j] = tmp/detx;

      for (int k = 1; k < dimension; out[j+k] = (x[j]*y[j+k] - y[j]*x[j+k])/detx, k++);
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
void printDataProbMatlabFile(NumericsMatrix * M, double * f, NumericsMatrix * H, double * w, int d, int n, int m, FILE * file)
{
  NumericsMatrix* M_dense = NM_create(NM_DENSE, m, m);
  NumericsMatrix* H_dense = NM_create(NM_DENSE, n*d, m);
  
  fprintf(file,"d = %3i;\n",d);
  fprintf(file,"n = %6i;\n",n);
  fprintf(file,"m = %6i;\n",m);
    
  NM_to_dense(M, M_dense);
  fprintf(file,"M = [");
  for(int i = 0; i < m; i++)
    {
      for(int j = 0; j < m; j++)
	{
	  fprintf(file,"%22.14e ", M_dense->matrix0[i+j*m]);
	}
      fprintf(file,";\n");
    }
  fprintf(file,"];\n");
  
  NM_to_dense(H, H_dense);
  fprintf(file,"H = [");
  for(int i = 0; i < n*d; i++)
    {
      for(int j = 0; j < m; j++)
	{
	  fprintf(file,"%22.14e ", H_dense->matrix0[i+j*n*d]);
	}
      fprintf(file,";\n");
    }
  fprintf(file,"];\n");
  
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

  NM_clear(M_dense);
  NM_clear(H_dense);
  free(M_dense);
  free(H_dense);
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
void printIteresProbMatlabFile(int iteration, double * v, double * u, double * r, int d, int n, int m, FILE * file)
{
  fprintf(file,"v(%3i,:) = [",iteration+1);
  for(int i = 0; i < m; i++)
    {
      fprintf(file, "%22.14e, ", v[i]);
    }
  fprintf(file,"];\n");
  
  fprintf(file,"u(%3i,:) = [",iteration+1);
  for(int i = 0; i < n*d; i++)
    {
      fprintf(file, "%22.14e, ", u[i]);
    }
  fprintf(file,"];\n");
  
  fprintf(file,"r(%3i,:) = [",iteration+1);
  for(int i = 0; i < n*d; i++)
    {
      fprintf(file, "%22.14e, ", r[i]);
    }
  fprintf(file,"];\n");
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
    data->starting_point->reaction[i] = 0.0351;
    if(i % d == 0)
      data->starting_point->reaction[i] = 0.2056;
  }

  /* ------ initialize the change of variable matrix P_mu ------- */
  data->P_mu = (IPM_change_of_variable*)malloc(sizeof(IPM_change_of_variable));
  data->P_mu->mat = NM_create(NM_SPARSE, nd, nd);
  NM_triplet_alloc(data->P_mu->mat, nd);
  data->P_mu->mat->matrix2->origin = NSM_TRIPLET;
  for(unsigned int i = 0; i < nd; ++i)
    if(i % d == 0)
      NM_entry(data->P_mu->mat, i, i, 1. / problem->mu[(int)(i/d)]);
    else
      NM_entry(data->P_mu->mat, i, i, 1.);

  /* ------ initialize the inverse P_mu_inv of the change of variable matrix P_mu ------- */
  data->P_mu->inv_mat = NM_create(NM_SPARSE, nd, nd);
  NM_triplet_alloc(data->P_mu->inv_mat, nd);
  data->P_mu->inv_mat->matrix2->origin = NSM_TRIPLET;
  for(unsigned int i = 0; i < nd; ++i)
    if(i % d == 0)
      NM_entry(data->P_mu->inv_mat, i, i, problem->mu[(int)(i/d)]);
    else
      NM_entry(data->P_mu->inv_mat, i, i, 1.);

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
  // the size of the problem detection
  unsigned int m = problem->M->size0;
  unsigned int nd = problem->H->size1;
  unsigned int d = problem->dimension;
  unsigned int n = problem->numberOfContacts;

  NumericsMatrix* M = NULL;
  NumericsMatrix* H_tilde = NULL;

  //globalFrictionContact_display(problem);
  
  /* if SICONOS_FRICTION_3D_IPM_FORCED_SPARSE_STORAGE = SICONOS_FRICTION_3D_IPM_FORCED_SPARSE_STORAGE,
     we force the copy into a NM_SPARSE storageType */
  DEBUG_PRINTF("problem->M->storageType : %i\n",problem->M->storageType);
  if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_SPARSE_STORAGE] == SICONOS_FRICTION_3D_IPM_FORCED_SPARSE_STORAGE
      && problem->M->storageType == NM_SPARSE_BLOCK)
  {
    DEBUG_PRINT("Force a copy to sparse storage type\n");
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

  double *iden;

  // change of variable
  // H_tilde --> H
  NumericsMatrix *H = NM_multiply(P_mu, H_tilde);
  // w_tilde --> w
  NM_gemv(1.0, P_mu, w_tilde, 0.0, w);

  // compute -H
  NumericsMatrix *minus_H = NM_create(H->storageType, H->size0, H->size1);
  NM_copy(H, minus_H);
  NM_gemm(-1.0, H, NM_eye(H->size1), 0.0, minus_H);

  double alpha_primal = data->internal_params->alpha_primal;
  double alpha_dual = data->internal_params->alpha_dual;
  double barr_param = data->internal_params->barr_param;
  double sigma = data->internal_params->sigma;

  cblas_dcopy(nd, data->starting_point->reaction, 1, reaction, 1);
  cblas_dcopy(nd, data->starting_point->velocity, 1, velocity, 1);
  cblas_dcopy(m, data->starting_point->globalVelocity, 1, globalVelocity, 1);

  double tol = options->dparam[SICONOS_DPARAM_TOL];
  unsigned int max_iter = options->iparam[SICONOS_IPARAM_MAX_ITER];

  double sgmp1 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1];
  double sgmp2 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2];
  double sgmp3 = options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3];
  double gmmp0 = 0.98;
  double gmmp1 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1];
  double gmmp2 = options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2];

  bool hasNotConverged = 1;
  unsigned int iteration = 0;
  double pinfeas = 1e300;
  double dinfeas = 1e300;
  double complem = 1e300;
  double dualgap = 1e300;
  double error[4];
  error[0] = pinfeas;
  error[1] = dinfeas;
  error[2] = dualgap;
  error[3] = complem;

  double gmm = gmmp0;
  double barr_param_a, e;
  double norm_f = cblas_dnrm2(m, f, 1);
  double norm_w = cblas_dnrm2(nd, w, 1);

  double *primalConstraint = data->tmp_vault_nd[1];
  double *dualConstraint = data->tmp_vault_m[0];
  double *complemConstraint = data->tmp_vault_nd[2];

  double *d_globalVelocity = (double*)calloc(m,sizeof(double));
  double *d_velocity = (double*)calloc(nd,sizeof(double));
  double *d_reaction = (double*)calloc(nd,sizeof(double));

  double *rhs = options->dWork;
  double *gv_plus_dgv = data->tmp_vault_m[1];
  double *vr_jprod = data->tmp_vault_nd[3];
  double *v_plus_dv = data->tmp_vault_nd[4];
  double *r_plus_dr = data->tmp_vault_nd[5];
  double *vr_prod_sub_iden = data->tmp_vault_nd[6];
  double *dvdr_jprod = data->tmp_vault_nd[7];

  int ReducedSystem = 1;
  int PrintMatlab = 0;
  
  double * r_p = (double*)calloc(nd,sizeof(double));                          // scaling vector p
  NumericsMatrix* r_Qp = NULL;                                                // matrix Qp 
  NumericsMatrix *minus_M = NM_create(M->storageType, M->size0, M->size1);    // store the matrix -M to build the matrix of the Newton linear system
  NumericsMatrix *QpH = NM_create(H->storageType, H->size0, H->size1);        // store the matrix Qp*H
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
  double err;
  char noNT = ' ';

  /* list of active constraints : = 0 if x_0 <= epsilon, = 1 if lambda_2 <= epsilon , = 3 either */
  short * a_velo = (short*)calloc(n, sizeof(short)); 
  short * a_reac = (short*)calloc(n, sizeof(short));
  
  if ((options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING]) && (ReducedSystem == 1))
    {
      /* Create the matrix -M to build the matrix of the reduced linear system */
      NM_copy(M, minus_M); /* useful ? */
      NM_gemm(-1.0, M, NM_eye(M->size1), 0.0, minus_M);
    }
  
  NumericsMatrix *J;
  
  long J_nzmax;

  size_t H_nzmax = NM_nnz(H);

  if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO] ==
      SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_YES)
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
  numerics_printf_verbose(-1, "(d, n, m) = (%1i, %6i, %6i)\n",d, n, m);
  numerics_printf_verbose(-1, "| it  |  dualgap  | pinfeas  | dinfeas  | complem  | full err | barparam | alpha_p  | alpha_d  |  sigma   | |dv|/|v| | |du|/|u| | |dr|/|r| |");
  numerics_printf_verbose(-1, "--------------------------------------------------------------------------------------------------------------------------------------------");

  double * p = data->tmp_vault_nd[8];
  double * p2 = data->tmp_vault_nd[9];
  double * pinv = data->tmp_vault_nd[10];
  NumericsMatrix* Qp = NULL;
  NumericsMatrix* Qp2 = NULL;
  NumericsMatrix* Qpinv = NULL;
  double * velocity_t = data->tmp_vault_nd[11];
  double * d_velocity_t = data->tmp_vault_nd[12];
  double * d_reaction_t = data->tmp_vault_nd[13];
  double * velocity_t_inv = data->tmp_vault_nd[14];
  double * Qp_velocity_t_inv = data->tmp_vault_nd[15];
  double * tmp1 = data->tmp_vault_nd[16];
  FILE * iteres; 

  /* writing data in a Matlab file */
  if (PrintMatlab == 1)
    {
      iteres = fopen("iteres.m", "w");
      printDataProbMatlabFile(M, f, H, w, d, n, m, iteres);
    }

  ComputeErrorGlobalPtr computeError = NULL;

  if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_UPDATE_S]==
      SICONOS_FRICTION_3D_ADMM_UPDATE_S_YES)
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
    if (NV_max(error, 4) <= tol)
      {
    	options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING] = 0;
    	noNT = '*';
      }
    
    ReducedSystem = (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING] == 0) ? 0 : ReducedSystem;
    
    if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING])
    {
      NesterovToddVector(velocity, reaction, nd, n, p);
      JA_power2(p, nd, n, p2);
      Qp = Quad_repr(p, nd, n);
      Qp2 = Quad_repr(p2, nd, n);
      JA_inv(p, nd, n, pinv);
      Qpinv = Quad_repr(pinv, nd, n);
     }

    primalResidualVector(velocity, H, globalVelocity, w, primalConstraint); // velocity - H * globalVelocity - w
    dualResidualVector(M, globalVelocity, H, reaction, f, dualConstraint);  // M*globalVelocity - H'*reaction - f
    /***** These two vectors should be used to calculate residual norms  *****/
    
    pinfeas = primalResidualNorm(velocity, H, globalVelocity, w);
    dinfeas = dualResidualNorm(M, globalVelocity, H, reaction, f);
    dualgap = dualGap(M, f, w, globalVelocity, reaction, nd, m);
    barr_param = cblas_ddot(nd, reaction, 1, velocity, 1) / nd;

    if(options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING])
      {
	Nesterov_Todd_vector(velocity, reaction, nd, n, r_p);           // Nesterov and Todd scaling p-vector
	complem = complemResidualNorm_p(velocity, reaction, nd, n);  // Norm of the Jordan product of the scaled vectors velocity and reaction
      }
    else
      {
	complem = complemResidualNorm(velocity, reaction, nd, n);
      }
    
    setErrorArray(error, pinfeas, dinfeas, dualgap, complem);


    /* ----- return to original variables ------ */
    NM_gemv(1.0, P_mu_inv, velocity, 0.0, data->tmp_point->t_velocity);
    /* cblas_dcopy(nd, data->tmp_point->t_velocity, 1, velocity, 1); */
    NM_gemv(1.0, P_mu, reaction, 0.0, data->tmp_point->t_reaction);
    /* cblas_dcopy(nd, data->tmp_point->t_reaction, 1, reaction, 1); */

    (*computeError)(problem,
                    data->tmp_point->t_reaction, data->tmp_point->t_velocity, globalVelocity,
                    tol, options,
                    norm_q, norm_b,  &err);

    // check exit condition
    if ((NV_max(error, 4) <= tol) && (err <= tol))
    {
      numerics_printf_verbose(-1, "| %3i%c| %9.2e | %.2e | %.2e | %.2e | %.2e | %.2e |",
			      iteration, noNT, dualgap, pinfeas, dinfeas, complem, err, barr_param);
      PrintMatlab == 1 ? printIteresProbMatlabFile(iteration, globalVelocity, velocity, reaction, d, n, m, iteres) : ' ' ;
      hasNotConverged = 0;
      break;
    }

    /*    1. Build the Jacobian matrix
     *
     *         m     nd       nd        
     *      |  M     0      -H^T  | m
     *      |                     |
     *  J = |  0   Arw(r)  Arw(u) | nd
     *      |                     |
     *      | -H     I        0   | nd
     *
     *  or
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

    if (ReducedSystem !=1) 
      {
	J = NM_create(NM_SPARSE, m + nd + nd, m + nd + nd);
	J_nzmax = (d * d) * (m / d) + H_nzmax + 2 * (d * 3 - 2) * n + H_nzmax + nd;
	NM_triplet_alloc(J, J_nzmax);
	J->matrix2->origin = NSM_TRIPLET;
	
	NM_insert(J, M, 0, 0);
	NM_insert(J, NM_transpose(minus_H), 0, m + nd);
	if(!options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING])
	  {
	    NM_insert(J, Arrow_repr(reaction, nd, n), m, m);
	    NM_insert(J, Arrow_repr(velocity, nd, n), m, m + nd);
	  }
	else
	  {
	    NM_insert(J, Qp2, m, m);
	    NM_insert(J, NM_eye(nd), m, m + nd);
	  }
	NM_insert(J, minus_H, m + nd, 0);
	NM_insert(J, NM_eye(nd), m + nd, m);
      }
    else
      {
	/* Create the 2x2 blocks reduced matrix
	 *         m       nd
	 *      |  -M     QpH^T  | m
	 * JR = |                |
	 *      |  QpH     I     | nd
	 *
	 *  where QpH = Qp * H
	 */
	JR = NM_create(NM_SPARSE, m + nd, m + nd);
	JR_nzmax = (d * d) * (m / d) + H_nzmax + H_nzmax + nd;
	NM_triplet_alloc(JR, JR_nzmax);
	JR->matrix2->origin = NSM_TRIPLET;
	NM_insert(JR, minus_M, 0, 0);
	NM_insert(JR, NM_eye(nd), m, m);
	r_Qp = Quad_repr(r_p, nd, n);                           // Should be replaced by a function returning the product Qp * vector 
	NM_copy(H, QpH); /* useful ? */
	NM_gemm(1.0, r_Qp, H, 0.0, QpH);                      
	NM_insert(JR, NM_transpose(QpH), 0, m);                // Should be useless when unsing a symmetric factorization procedure 
	NM_insert(JR, QpH, m, 0);
      }
    

    /** Correction of w to take into account the dependence
        on the tangential velocity */
    if(options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_UPDATE_S]==
        SICONOS_FRICTION_3D_ADMM_UPDATE_S_YES)
    {
      for(unsigned int i = 0; i < nd; ++ i)
      {
        if(i % d == 0)
          w[i] = w_tilde[i]/(problem->mu[(int)(i/d)])
                 + sqrt(velocity[i+1]*velocity[i+1]+velocity[i+2]*velocity[i+2]);
      }
    }

    PrintMatlab == 1 ? printIteresProbMatlabFile(iteration, globalVelocity, velocity, reaction, d, n, m, iteres) : ' ' ;

    if (ReducedSystem != 1)
      {
	/* building the right-hand side related to the first system 
           
	   without NT scaling
	    
                     [ M*v - H'*r - f ]  m         dualConstraint
             rhs = - [     u o r      ]  nd        complemConstraint
                     [  u - H*v - w   ]  nd        primalConstraint
           
	   with NT scaling
	 
                     [ M*v - H'*r - f ]  m         dualConstraint
             rhs = - [       r        ]  nd        complemConstraint
                     [  u - H*v - w   ]  nd        primalConstraint
	*/

	cblas_dcopy(m, dualConstraint, 1, rhs, 1);
	cblas_dcopy(nd, primalConstraint, 1, rhs+m+nd, 1);
	if(!options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING])
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
	/* building the reduced right-hand side

                     [ M*v - f - H'*r ]  m         dualConstraint
             r_rhs = [                ]
                     [ -Qp*(H*v + w)  ]  nd        primalConstraint - velocity

	*/
	NV_insert(r_rhs, m + nd, dualConstraint, m, 0);
	NM_gemv(1.0, H, globalVelocity, 0.0, Hvw);
	cblas_daxpy(nd, 1.0, w, 1, Hvw, 1);              // Hvw <- H*v + w, kept in memory for the computation of r_du
	QNTpz(velocity, reaction, Hvw, nd, n, r_rhs+m);  // Could also be done with Qxy(r_p, Hvw, nd, n, r_rhs+m); 
        cblas_dscal(nd, -1.0, r_rhs+m, 1);         
        cblas_dcopy(m+nd, r_rhs, 1, r_rhs_2, 1);         // kept a copy of r_rhs that will be used for the solution of the second linear system
      }

    if (ReducedSystem != 1)
      {
	if (!options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING])
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
	QNTpz(velocity, reaction, r_rhs+m, nd, n, d_reaction);
	NM_gemv(1.0, H, d_globalVelocity, 0.0, d_velocity);
	cblas_daxpy(nd, 1.0, Hvw, 1, d_velocity, 1);
	cblas_daxpy(nd, -1.0, velocity, 1, d_velocity, 1);
     }

    /* computing the affine step-length */
    alpha_primal = getStepLength(velocity, d_velocity, nd, n, gmm);
    alpha_dual = getStepLength(reaction, d_reaction, nd, n, gmm);

    /* updating the gamma parameter used to compute the step-length */
    gmm = gmmp1 + gmmp2 * fmin(alpha_primal, alpha_dual);

    /* ----- Predictor step of Mehrotra ----- */
    cblas_dcopy(nd, velocity, 1, v_plus_dv, 1);
    cblas_dcopy(nd, reaction, 1, r_plus_dr, 1);
    cblas_daxpy(nd, alpha_primal, d_velocity, 1, v_plus_dv, 1);
    cblas_daxpy(nd, alpha_dual, d_reaction, 1, r_plus_dr, 1);
    
    /* affine barrier parameter */	
    barr_param_a = cblas_ddot(nd, v_plus_dv, 1, r_plus_dr, 1) / nd;
 
    /* computing the centralization parameter */
    e = barr_param > sgmp1 ? fmax(1.0, sgmp2 * pow(fmin(alpha_primal, alpha_dual),2)) : sgmp3;
    sigma = fmin(1.0, pow(barr_param_a / barr_param, e));
    
    /* computing the corrector step of Mehrotra algorithm */
    if (ReducedSystem != 1)
      {
	if (!options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING])
	  {
	    /* Right-hand side for non-symmetric Newton system without NT scaling */
	    iden = JA_iden(nd, n);
	    cblas_dscal(nd, 2 * barr_param * sigma, iden, 1);
	    JA_prod(velocity, reaction, nd, n, vr_jprod);
	    //JA_prod(d_velocity, d_reaction, nd, n, dvdr_jprod);
	    NV_sub(vr_jprod, iden, nd, complemConstraint);
	    //NV_add(vr_prod_sub_iden, dvdr_jprod, nd, complemConstraint);
	    free(iden);
	  }
	else
	  {
	    /* Right-hand side for symmetric Newton system with NT scaling */
	    NM_gemv(1.0, Qp, velocity, 0.0, velocity_t);
	    NM_gemv(1.0, Qp, d_velocity, 0.0, d_velocity_t);
	    NM_gemv(1.0, Qpinv, d_reaction, 0.0, d_reaction_t);
	    JA_inv(velocity_t, nd, n, velocity_t_inv);
	    NM_gemv(1.0, Qp, velocity_t_inv, 0.0, Qp_velocity_t_inv);
	    JA_inv(velocity, nd, n, Qp_velocity_t_inv); 
	    cblas_dscal(nd, 2 * barr_param * sigma, Qp_velocity_t_inv, 1);
	    JA_prod(d_velocity_t, d_reaction_t, nd, n, dvdr_jprod);
	    NV_sub(reaction, Qp_velocity_t_inv, nd, tmp1);
	    /* NV_add(tmp1, dvdr_jprod, nd, complemConstraint); */
	    JA_prod(velocity_t_inv, dvdr_jprod, nd, n, d_velocity_t); // PA: left-Jordan-product by velocity_t_inv
	    NM_gemv(1.0, Qp, d_velocity_t, 0.0, complemConstraint); // PA: left-Jordan procuct by Qp
	    NV_add(tmp1, complemConstraint, nd, complemConstraint);
	  }
        cblas_dcopy(m, dualConstraint, 1, rhs, 1);
	cblas_dcopy(nd, complemConstraint, 1, rhs+m, 1);
	cblas_dcopy(nd, primalConstraint, 1, rhs+m+nd, 1);
	cblas_dscal(m + nd + nd, -1.0, rhs, 1);

	/* Newton system solving */
	if (!options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING])
	  {
	    //NM_gesv_expert(J, rhs, NM_KEEP_FACTORS);
	    NM_LU_solve(J, rhs, 1);
	  }
	else
	  {
	    NSM_linearSolverParams(J)->solver = NSM_HSL;
	    NM_LDLT_solve(J, rhs, 1);
	  }
	
	NM_clear(J);
	free(J);
	
	d_globalVelocity = rhs;
	d_velocity = rhs + m;
	d_reaction = rhs + m + nd;
      }
    else
      {
	/* right-hand side for the second reduced Newton system with NT scaling */
	double *r_Qp_u = (double*)calloc(nd,sizeof(double));
	double *r_Qp_u_inv = (double*)calloc(nd,sizeof(double));
	double * r_Qp_du = (double*)calloc(nd,sizeof(double));
	double * r_dudr = (double*)calloc(nd,sizeof(double));
	double * r_ududr = (double*)calloc(nd,sizeof(double));

	Qxy(r_p, velocity, nd, n, r_Qp_u);                                 // u_hat
	
	cblas_dcopy(nd, r_Qp_u, 1, r_Qp_du, 1);
	cblas_daxpy(nd, 1.0, r_rhs+m, 1, r_Qp_du, 1); 
	cblas_dscal(nd, -1.0, r_Qp_du, 1);                                // du_hat = -(dr_check + u_hat)
	
	JA_prod(r_Qp_du, r_rhs+m, nd, n, r_dudr);                         // du_hat o dr_check
	
        for (int k = 0; k < nd; r_dudr[k] -= 2*sigma*barr_param, k+=d);   // du_hat o dr_check - 2*sigma*mu*e
	Jxinvprody(r_Qp_u, r_dudr, nd, n, r_ududr);                       // u_hat^{-1} o  (du_hat o dr_check - 2*sigma*mu*e)

	cblas_daxpy(nd, -1.0, r_ududr, 1, r_rhs_2+m, 1);                  // -Qp*(Hv+w) - u_hat^{-1} o (du_hat o dr_check) + 2*sigma*mu*u_hat^{-1}

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
	Qxy(r_p, r_rhs_2+m, nd, n, r_dr); 
	NM_gemv(1.0, H, r_dv, 0.0, r_du);
	NV_add(r_du, Hvw, nd, r_du);
	NV_sub(r_du, velocity, nd, r_du);

	d_globalVelocity = r_dv;
	d_velocity = r_du;
	d_reaction = r_dr;	
      }


    alpha_primal = getStepLength(velocity, d_velocity, nd, n, gmm);
    alpha_dual = getStepLength(reaction, d_reaction, nd, n, gmm);

    /* updating the gamma parameter used to compute the step-length */
    gmm = gmmp1 + gmmp2 * fmin(alpha_primal, alpha_dual);

    numerics_printf_verbose(-1, "| %3i%c| %9.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e | %.2e |",
                            iteration, noNT, dualgap, pinfeas, dinfeas, complem, err, barr_param, alpha_primal, alpha_dual, sigma,
			    cblas_dnrm2(m, d_globalVelocity, 1)/cblas_dnrm2(m, globalVelocity, 1),
			    cblas_dnrm2(nd, d_velocity, 1)/cblas_dnrm2(nd, velocity, 1),
			    cblas_dnrm2(nd, d_reaction, 1)/cblas_dnrm2(nd, reaction, 1));
        
    /* ----- Update variables ----- */
    cblas_daxpy(m, alpha_primal, d_globalVelocity, 1, globalVelocity, 1);
    cblas_daxpy(nd, alpha_primal, d_velocity, 1, velocity, 1);
    cblas_daxpy(nd, alpha_dual, d_reaction, 1, reaction, 1);

    iteration++;
  
  } // while loop

  /* printing complementarity products */
  printf("\nsmallest eigenvalue | ");
  printf("velocity               reaction             velocity+reaction\n");  
  double * veloprea = (double*)calloc(nd, sizeof(double));
  cblas_dcopy(nd, reaction, 1, veloprea, 1);
  cblas_daxpy(nd, 1.0, velocity, 1, veloprea, 1);
  long double vpv, vpr;
  for (int i = 0; i < (i>10 ? 10:n); i++)
    {
      vpv = velocity[i*d+1]*velocity[i*d+1];
      for (int k = 2; k < d; vpv = vpv + velocity[i*d+k]*velocity[i*d+k], k++);
      vpv = velocity[i*d] - sqrtl(vpv);
      vpr = reaction[i*d+1]*reaction[i*d+1];
      for (int k = 2; k < d; vpr = vpr + reaction[i*d+k]*reaction[i*d+k], k++);
      vpr = reaction[i*d] - sqrtl(vpr);
      printf("%19c %3i %20.14e %20.14e %20.14e\n", ' ', i,
      	     velocity[i*d]-cblas_dnrm2(d-1, velocity+i*d+1, 1),
      	     reaction[i*d]-cblas_dnrm2(d-1, reaction+i*d+1, 1),
      	     veloprea[i*d]-cblas_dnrm2(d-1, veloprea+i*d+1, 1));
      printf("%20c    %20.14Le %20.14Le\n", 'L', vpv, vpr);
    }
  free(veloprea);

  /* determining active constraints */
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
  
  
  /* ----- return to original variables ------ */
  NM_gemv(1.0, P_mu_inv, velocity, 0.0, data->tmp_point->t_velocity);
  cblas_dcopy(nd, data->tmp_point->t_velocity, 1, velocity, 1);

  NM_gemv(1.0, P_mu, reaction, 0.0, data->tmp_point->t_reaction);
  cblas_dcopy(nd, data->tmp_point->t_reaction, 1, reaction, 1);

  options->dparam[SICONOS_DPARAM_RESIDU] = NV_max(error, 4);
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

  PrintMatlab == 1 ? fclose(iteres) : ' ';
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

  *info = hasNotConverged;
}

void gfc3d_ipm_set_default(SolverOptions* options)
{

  options->iparam[SICONOS_IPARAM_MAX_ITER] = 200;
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO] =
    SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_NO;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_NESTEROV_TODD_SCALING] = 1;

  options->iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_UPDATE_S]=
    SICONOS_FRICTION_3D_ADMM_UPDATE_S_NO;


  options->dparam[SICONOS_DPARAM_TOL] = 1e-12;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1] = 1e-8;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2] = 3.;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3] = 1.;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1] = 0.9;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2] = 0.095; //0.095

}
