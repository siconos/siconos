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

#include "projectionOnCone.h"

/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "siconos_debug.h"
#include "gfc3d_ipm.h"
#include <stdarg.h>         // for va_list, va_start, va_end

const char* const   SICONOS_GLOBAL_FRICTION_3D_IPM_SNM_STR = "GFC3D IPM SNM";

/* ------------------------- Helper functions implementation ------------------------------ */
/* Declairation of variadic function pointer */
typedef int (*Check_merite_function)(const double, const double, const double, va_list);

/* BACKTRACKING LINE-SEARCH
 * Description: compute a step-length alpha such that
 *              f(x + alpha * dx) <= f(x) + alpha * omega * nabla f(x)^T dx
 * [in] _func: Merite function is used for evaluation
 * [in] lb, ub: lb <= alpha <= ub
 * [in] omega: a constant could be 0.01, 0.001, 0.0001
 * [in] numParams: number of parameters as input for _func
 * [out] alpha: step-length
 */
double backtrack_linesearch(Check_merite_function _func, const double lb, const double ub, const double omega, const double numParams, ...)
{
  assert(("Lower bound <= upper bound", lb <= ub));
  va_list args;
  va_start(args, numParams);
  double alpha = ub;

  while (_func(alpha, omega, numParams, args) == 0 && lb < alpha)
  {
    alpha /= 2.;
    // alpha -= 0.00001;
    if (alpha < lb) alpha = lb;
    va_start(args, numParams);
  }

  va_end(args);
  return alpha;
}


/* This routine is to verify the inequality :
 * | (u + alpha * du) o (r + alpha * dr) | <= (1 - alpha * omega) | u o r | + 2 * alpha * omega * sigma * (u'r)^2 / n / | u o r |
 *
 * [in] alpha: step-length
 * [in] omega: a constant could be 0.01, 0.001, 0.0001
 * [in] sigma in ]0, 0.5[: param of the algo
 * [in] numParams: number of parameters must be 7
 * [in] ATTENTION others params must be: double *u, double *du, double *r, double *dr,
                               unsigned int vecSize, unsigned int varsCount, double sigma
 * [out] 1: if satisfy, 0: it not
 */
int check_merite_function_complem(const double alpha, const double omega, const double numParams, va_list args)
{
  assert(("The number of params must be 7. WARNING: Params input must be u, du, r, dr, vecSize, varsCount, sigma", numParams  == 7));

  // Read params
  double * u = va_arg(args, double*);
  double * du = va_arg(args, double*);
  double * r = va_arg(args, double*);
  double * dr = va_arg(args, double*);
  unsigned int vecSize = va_arg(args, unsigned int);
  unsigned int varsCount = va_arg(args, unsigned int);
  double sigma = va_arg(args, double);

  double * uor = (double*)calloc(vecSize, sizeof(double));
  double * uplus = (double*)calloc(vecSize, sizeof(double));
  double * rplus = (double*)calloc(vecSize, sizeof(double));
  double rhs = 0.;

  JA_prod(u, r, vecSize, varsCount, uor);
  rhs = cblas_dnrm2(vecSize, uor, 1);             // rhs = | u o r |
  rhs = (1. - alpha * omega ) * rhs + 2 * alpha * omega * sigma * pow(cblas_ddot(vecSize, u, 1, r, 1),2) / varsCount / rhs;

  cblas_dcopy(vecSize, u, 1, uplus, 1);
  cblas_daxpy(vecSize, alpha, du, 1, uplus, 1);   // uplus = u + alpha * du
  cblas_dcopy(vecSize, r, 1, rplus, 1);
  cblas_daxpy(vecSize, alpha, dr, 1, rplus, 1);   // rplus = r + alpha * dr

  // printf("| u o r | = %e, rhs 1 = %e\n", cblas_dnrm2(vecSize, uor, 1), rhs);

  JA_prod(uplus, rplus, vecSize, varsCount, uor); // rhs   = uplus o rplus

  rhs -= cblas_dnrm2(vecSize, uor, 1);            // rhs   = (1 - alpha * omega) | u o r | + 2 * alpha * omega * sigma * (u'r)^2 / n / | u o r | - | uplus o rplus |

  // printf("|u+ o r+| = %e, rhs 2 = %e\n", cblas_dnrm2(vecSize, uor, 1), rhs);

  free(uor);
  free(uplus);
  free(rplus);

  if (rhs >= 0) return 1;
  else return 0;
}


/* This routine is to verify the inequality :
 * | s + alpha * ds - |ub + alpha * dub| | <= (1 - alpha * omega) | s - ub |
 *
 * [in] alpha: step-length
 * [in] omega: a constant could be 0.01, 0.001, 0.0001
 * [in] numParams: number of parameters must be 6
 * [in] ATTENTION others params must be: double *u, double *du, double *s, double *ds,
                               unsigned int vecSize, unsigned int varsCount
 * [out] 1: if satisfy, 0: it not
 */
int check_merite_function_diffixP(const double alpha, const double omega, const double numParams, va_list args)
{
  assert(("The number of params must be 6. WARNING: Params input must be u, du, s, ds, vecSize, varsCount", numParams  == 6));

  // Read params
  double * u = va_arg(args, double*);
  double * du = va_arg(args, double*);
  double * s = va_arg(args, double*);
  double * ds = va_arg(args, double*);
  unsigned int vecSize = va_arg(args, unsigned int);
  unsigned int varsCount = va_arg(args, unsigned int);

  double * uplus = (double*)calloc(vecSize, sizeof(double));
  double * splus = (double*)calloc(varsCount, sizeof(double));
  double rhs = 0.;

  double diff_fixp = 0., nub = 0;
  unsigned int d = vecSize / varsCount;

  for (unsigned int i = 0; i<vecSize; i+=d)
  {
    nub = cblas_dnrm2(2, u+i+1, 1);
    diff_fixp += (s[i/d] - nub)*(s[i/d] - nub);
  }
  diff_fixp = sqrt(diff_fixp);                    // diff_fixp = | s - ub |

  rhs = (1. - alpha * omega) * diff_fixp;

  // printf("| s - ub | = %e, rhs 1 = %e, alpha = %e, s = %e, nub = %e, ", diff_fixp, rhs, alpha, s[0], nub);

  cblas_dcopy(vecSize, u, 1, uplus, 1);
  cblas_daxpy(vecSize, alpha, du, 1, uplus, 1);   // uplus = u + alpha * du
  cblas_dcopy(varsCount, s, 1, splus, 1);
  cblas_daxpy(varsCount, alpha, ds, 1, splus, 1); // splus = s + alpha * ds

  diff_fixp = 0.;
  for (unsigned int i = 0; i<vecSize; i+=d)
  {
    nub = cblas_dnrm2(2, uplus+i+1, 1);
    diff_fixp += (splus[i/d] - nub)*(splus[i/d] - nub);
  }
  diff_fixp = sqrt(diff_fixp);                    // diff_fixp = | splus - uplus_bar |

  // printf("s+ = %e, nub+ = %e, |s+ - ub+| = %e\n", splus[0], nub, diff_fixp);

  rhs -= diff_fixp;                               // rhs   = (1 - alpha * omega) | s - ub | - | splus - uplus_bar |

  free(uplus);
  free(splus);

  if (rhs >= 0) return 1;
  else return 0;
}



/* This routine is to verify the inequality :
 * Theta((v,u,r,s)+alpha*(dv,du,dr,ds)) <= Theta(v,u,r,s) + alpha * omega * Theta'((v,u,r,s);(dv,du,dr,ds))
 *
 * where
 *
 *                  | [    Mv - H'r - f   ] |   dualConstraint
 * Theta(v,u,r,s) = | [        u o r      ] |   complemConstraint
 *                  | [  u - Hv - w - Es  ] |   primalConstraint
 *                  | [     s - u_tilde   ] |   fixpConstraint
 *
 * and
 *
 *   Theta'((v,u,r,s);(dv,du,dr,ds))
 * = -|Mv - H'r - f| - |u - Hv - w - Es| - |uor| + 2*sigma*(u'r)^2 / n / | u o r | - |s - u_tilde|
 *
 * [in] alpha: step-length
 * [in] omega: a constant could be 0.01, 0.001, 0.0001
 * [in] sigma in ]0, 0.5[: param of the algo
 * [in] numParams: number of parameters must be 20
 * [in] ATTENTION others params must be:
            double *v, double *dv, double *u, double *du, double *r, double *dr, double *s, double *ds,
            unsigned int nd, unsigned int n, double sigma, double pinfeas, double dinfeas, double complem, double diff_fixp,
            NumericsMatrix *M, NumericsMatrix *H, double *f, double *w, double tol
 * [out] 1: if satisfy, 0: it not
 */
int check_merite_function_Theta(const double alpha, const double omega, const double numParams, va_list args)
{
  assert(("The number of params must be 20. \nWARNING: Params input must be v, dv, u, du, r, dr, s, ds, vecSize, varsCount, \\sigma, pinfeas, dinfeas, complem, diff_fixp, \\M, H, f, w, tol", numParams  == 20));

  // Read params
  double * v = va_arg(args, double*);
  double * dv = va_arg(args, double*);
  double * u = va_arg(args, double*);
  double * du = va_arg(args, double*);
  double * r = va_arg(args, double*);
  double * dr = va_arg(args, double*);
  double * s = va_arg(args, double*);
  double * ds = va_arg(args, double*);
  unsigned int nd = va_arg(args, unsigned int);
  unsigned int n = va_arg(args, unsigned int);
  double sigma = va_arg(args, double);
  double pinfeas = va_arg(args, double);
  double dinfeas = va_arg(args, double);
  double complem = va_arg(args, double);
  double diff_fixp = va_arg(args, double);
  NumericsMatrix * M = va_arg(args, NumericsMatrix*);
  NumericsMatrix * H = va_arg(args, NumericsMatrix*);
  double * f = va_arg(args, double*);
  double * w = va_arg(args, double*);
  double tol = va_arg(args, double);

  unsigned int m = M->size0;
  unsigned int d = nd / n;

  // Initiate temporary vars
  double * vplus = (double*)calloc(m, sizeof(double));
  double * uplus = (double*)calloc(nd, sizeof(double));
  double * rplus = (double*)calloc(nd, sizeof(double));
  double * splus = (double*)calloc(n, sizeof(double));

  double * primalConstraint = (double*)calloc(nd, sizeof(double));
  double * dualConstraint = (double*)calloc(m, sizeof(double));

  double rhs = 0., nub = 0.;

  rhs = - dinfeas - pinfeas - complem + 2*sigma*pow(cblas_ddot(nd, u, 1, r, 1), 2) / n / complem - diff_fixp;

  if (rhs >= 0.)
  {
    printf("WARNING: Directional derivative Theta' = %e >= 0. alpha = %e\n", rhs, alpha);
  }


  // Compute Theta(v,u,r,s)
  rhs = alpha*omega*rhs + sqrt(pow(pinfeas,2)+pow(dinfeas,2)+pow(complem,2)+pow(diff_fixp,2));


  // Compute Theta(v+,u+,r+,s+)
  cblas_dcopy(m, v, 1, vplus, 1);
  cblas_dcopy(nd, u, 1, uplus, 1);
  cblas_dcopy(nd, r, 1, rplus, 1);
  cblas_dcopy(n, s, 1, splus, 1);

  cblas_daxpy(m, alpha, dv, 1, vplus, 1);
  cblas_daxpy(nd, alpha, du, 1, uplus, 1);
  cblas_daxpy(nd, alpha, dr, 1, rplus, 1);
  cblas_daxpy(n, alpha, ds, 1, splus, 1);

  primalResidual_s(uplus, H, vplus, w, splus, primalConstraint, &pinfeas, tol);
  dualResidual(M, vplus, H, rplus, f, dualConstraint, &dinfeas, tol);
  complem = complemResidualNorm(uplus, rplus, nd, n);
  diff_fixp = 0.;
  for (unsigned int i = 0; i<nd; i+=d)
  {
    nub = cblas_dnrm2(2, uplus+i+1, 1);
    diff_fixp += (splus[i/d] - nub)*(splus[i/d] - nub);
  }

  rhs -= sqrt(pow(pinfeas,2)+pow(dinfeas,2)+pow(complem,2)+diff_fixp);

  free(vplus); free(uplus); free(rplus); free(splus);
  free(primalConstraint); free(dualConstraint);

  if (rhs >= 0) return 1;
  else return 0;
}





/* This routine is to verify the inequality :
 * Theta(vk+,uk+,rk+,sk+) <= max_{0<=j<=p} Theta(v_{k-j},u_{k-j},r_{k-j},s_{k-j}) + alpha * omega * Theta'((vk,uk,rk,sk);(dvk,duk,drk,dsk))
 *
 * where
 *
 * (vk+,uk+,rk+,sk+) = (vk,uk,rk,sk) + alpha * (dvk,duk,drk,dsk)
 *
 * [in] alpha: step-length
 * [in] omega: a constant could be 0.01, 0.001, 0.0001
 * [in] sigma in ]0, 0.5[: param of the algo
 * [in] numParams: number of parameters must be 23
 * [in] ATTENTION others params must be:
            double *v, double *dv, double *u, double *du, double *r, double *dr, double *s, double *ds,
            unsigned int nd, unsigned int n, double sigma, double pinfeas, double dinfeas, double complem, double diff_fixp,
            unsigned int iteration, double *arr_norm_theta, unsigned int p, NumericsMatrix *M, NumericsMatrix *H, double *f, double *w, double tol
 * [out] 1: if satisfy, 0: it not
 */
int check_merite_function_Theta_nonMonotone(const double alpha, const double omega, const double numParams, va_list args)
{
  assert(("The number of params must be 23. \nWARNING: Params input must be v, dv, u, du, r, dr, s, ds, vecSize, varsCount, \\sigma, pinfeas, dinfeas, complem, diff_fixp, \\iteration, arr_norm_theta, p, M, H, f, w, tol", numParams  == 23));

  // Read params
  double * v = va_arg(args, double*);
  double * dv = va_arg(args, double*);
  double * u = va_arg(args, double*);
  double * du = va_arg(args, double*);
  double * r = va_arg(args, double*);
  double * dr = va_arg(args, double*);
  double * s = va_arg(args, double*);
  double * ds = va_arg(args, double*);
  unsigned int nd = va_arg(args, unsigned int);
  unsigned int n = va_arg(args, unsigned int);
  double sigma = va_arg(args, double);
  double pinfeas = va_arg(args, double);
  double dinfeas = va_arg(args, double);
  double complem = va_arg(args, double);
  double diff_fixp = va_arg(args, double);
  unsigned int iteration = va_arg(args, unsigned int);
  double * arr_norm_theta = va_arg(args, double*);
  unsigned int p = va_arg(args, unsigned int);
  NumericsMatrix * M = va_arg(args, NumericsMatrix*);
  NumericsMatrix * H = va_arg(args, NumericsMatrix*);
  double * f = va_arg(args, double*);
  double * w = va_arg(args, double*);
  double tol = va_arg(args, double);

  unsigned int m = M->size0;
  unsigned int d = nd / n;

  // Initiate temporary vars
  double * vplus = (double*)calloc(m, sizeof(double));
  double * uplus = (double*)calloc(nd, sizeof(double));
  double * rplus = (double*)calloc(nd, sizeof(double));
  double * splus = (double*)calloc(n, sizeof(double));

  double * primalConstraint = (double*)calloc(nd, sizeof(double));
  double * dualConstraint = (double*)calloc(m, sizeof(double));

  double rhs = 0., nub = 0., max_theta = 0.;

  // rhs = Theta'((v_k,u_k,r_k,s_k); (dv_k,du_k,dr_k,ds_k))
  rhs = - dinfeas - pinfeas - complem + 2*sigma*pow(cblas_ddot(nd, u, 1, r, 1), 2) / n / complem - diff_fixp;

  if (rhs >= 0.)
  {
    printf("WARNING: Directional derivative Theta' = %e >= 0. alpha = %e\n", rhs, alpha);
  }


  // Compute max_{0<=j<=p} Theta(v_{k-j},u_{k-j},r_{k-j},s_{k-j})
  if (iteration >= p && p >= 2)
  {
    max_theta = arr_norm_theta[iteration-1];
    for(unsigned int i=iteration-2; i>=iteration-p; i--)
    {
      if (arr_norm_theta[i] > max_theta) max_theta = arr_norm_theta[i];
    }
  }
  else
  {
    printf("\ncheck_merite_function_Theta_nonMonotone failed! Param p must be >= 2 and the current iteration must be >= p.\n");
  }

  rhs = max_theta + alpha*omega*rhs;


  // Compute Theta(vk+,uk+,rk+,sk+)
  cblas_dcopy(m, v, 1, vplus, 1);
  cblas_dcopy(nd, u, 1, uplus, 1);
  cblas_dcopy(nd, r, 1, rplus, 1);
  cblas_dcopy(n, s, 1, splus, 1);

  cblas_daxpy(m, alpha, dv, 1, vplus, 1);
  cblas_daxpy(nd, alpha, du, 1, uplus, 1);
  cblas_daxpy(nd, alpha, dr, 1, rplus, 1);
  cblas_daxpy(n, alpha, ds, 1, splus, 1);

  primalResidual_s(uplus, H, vplus, w, splus, primalConstraint, &pinfeas, tol);
  dualResidual(M, vplus, H, rplus, f, dualConstraint, &dinfeas, tol);
  complem = complemResidualNorm(uplus, rplus, nd, n);
  diff_fixp = 0.;
  for (unsigned int i = 0; i<nd; i+=d)
  {
    nub = cblas_dnrm2(2, uplus+i+1, 1);
    diff_fixp += (splus[i/d] - nub)*(splus[i/d] - nub);
  }

  rhs -= sqrt(pow(pinfeas,2)+pow(dinfeas,2)+pow(complem,2)+diff_fixp);

  free(vplus); free(uplus); free(rplus); free(splus);
  free(primalConstraint); free(dualConstraint);

  if (rhs >= 0.) return 1;
  else return 0;
}





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


// /* Computation of the projection error |r - proj(r-u)|/max{|r|, |u|} */
// static double projectionError(const double * velocity, const double * reaction, const unsigned int nc, const double tol)
// {
//    double worktmp[3];
//    double out = 0.0;
//    double norm_u, norm_r, relative_scaling;

//    for(int ic = 0 ; ic < nc ; ic++)
//      {
//        worktmp[0] = reaction[3*ic] -  velocity[3*ic] ;
//        worktmp[1] = reaction[3*ic+1] -  velocity[3*ic+1] ;
//        worktmp[2] = reaction[3*ic+2] -  velocity[3*ic+2] ;
//        projectionOnCone(worktmp, 1.0);
//        worktmp[0] = reaction[3*ic] -  worktmp[0];
//        worktmp[1] = reaction[3*ic+1] -  worktmp[1];
//        worktmp[2] = reaction[3*ic+2] -  worktmp[2];
//        out +=  worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2];
//      }
//    out = sqrt(out);
//    norm_u = cblas_dnrm2(3*nc, velocity, 1);
//    norm_r = cblas_dnrm2(3*nc, reaction, 1);
//    relative_scaling = fmax(norm_u, norm_r);
//    if(relative_scaling > tol)
//      out = out/relative_scaling;
//    return out;
// }


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

  // fprintf(file,"du(%3i,:) = [",iteration+1);
  // for(int i = 0; i < n*d; i++)
  // {
  //   fprintf(file, "%8.20e, ", du[i]);
  // }
  // fprintf(file,"];\n");

  // fprintf(file,"dr(%3i,:) = [",iteration+1);
  // for(int i = 0; i < n*d; i++)
  // {
  //   fprintf(file, "%8.20e, ", dr[i]);
  // }
  // fprintf(file,"];\n");

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

// cl: classify
static void printBlockVec(double * vec, int vecSize, int sizeBlock, int cl)
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


static float randomFloat(float min, float max) {
    return min + (float)rand() / ((float)RAND_MAX / (max - min));
}

// typedef struct Triplet {
//     int row_index;
//     int col_index;
//     double value;
// } Triplet;

// Comparison function for qsort
// int compareTriplets(const void *a, const void *b) {
//     return ((Triplet *)a)->row_index - ((Triplet *)b)->row_index;
// }

// #include <stdlib.h>
// void sortTriplets(Triplet *triplets, size_t num_triplets) {
//     qsort(triplets, num_triplets, sizeof(Triplet), compareTriplets);
// }

// // This function is to delete 3*n columns of matrix H related to an array of cones which needs to be deleted
// // The 1st cone is 0
// // Not allocation and H should be stored as in CSC type after deletion
// void NM_clear_cone_matrix_H(NumericsMatrix *H, unsigned int n_cones_to_clear, int *cones_to_clear)
// {
//   // #include "CSparseMatrix.h"
//   switch(H->storageType)
//   {
//   case NM_SPARSE:
//   {
//     assert(H->matrix2);
//     if (H->matrix2->origin == NSM_CSC)
//     {
//       NM_triplet(H);
//       H->matrix2->origin= NSM_TRIPLET;
//       NM_clearCSC(H);
//     }

//     if (H->matrix2->origin == NSM_TRIPLET)
//     {
//       CSparseMatrix *H_triplet = H->matrix2->triplet;
//       int first = 0, last = H_triplet->nz-1, delete_counter = 0, stop = 0;
//       cs_long_t *rows = H_triplet->i, *cols = H_triplet->p, target = -1;
//       double *val = H_triplet->x;
//       for (unsigned int i=0; i<n_cones_to_clear; i++)
//       {
//         target = (cs_long_t)cones_to_clear[i]*3;
//         if (target > H_triplet->m)
//         {
//           n_cones_to_clear--;
//           continue;
//         }

//         first = 0; stop = 0;
//         while (first <= last && target<H_triplet->m && !stop)
//         {
//           if(rows[last] == target || rows[last] == target+1 || rows[last] == target+2) // Last row is to be deleted
//           {
//             last--;
//             delete_counter++;
//           }
//           else
//           {
//             for (unsigned int j=first; j<last; j++)
//             {
//               if (rows[j] == target || rows[j] == target+1 || rows[j] == target+2) // Search is done
//               {
//                 // Swap this value with the last one
//                 rows[j] = rows[last];
//                 cols[j] = cols[last];
//                 val[j] = val[last];
//                 last--;
//                 delete_counter++;
//                 first = j+1;
//                 break;
//               }
//               if (j == last-1) stop = 1;
//             }
//           }
//         }
//       }
//       // Update nz of H
//       H_triplet->nz = last+1;

//       // Reduce size of H
//       H_triplet->m -= 3*n_cones_to_clear;

//       // Sort row_index
//       Triplet *triplets = malloc(H_triplet->nz*sizeof(Triplet));
//       for (int k=0; k<H_triplet->nz; k++)
//       {
//         triplets[k].row_index = rows[k];
//         triplets[k].col_index = cols[k];
//         triplets[k].value = val[k];
//       }
//       sortTriplets(triplets, H_triplet->nz);

//       for (int k=0; k<H_triplet->nz; k++)
//       {
//         rows[k] = triplets[k].row_index;
//         cols[k] = triplets[k].col_index;
//         val[k] = triplets[k].value;
//       }

//       for (unsigned int i=0; i<n_cones_to_clear; i++)
//       {
//         target = ((cs_long_t)cones_to_clear[i]-i)*3;
//         if (target > H_triplet->m) continue;

//         for (unsigned int j=0; j<=last; j++)
//         {
//           if (rows[j] >= target)
//           {
//             rows[j] -= 3;
//           }
//         }
//       }

//       NM_csc(H);
//       H->matrix2->origin= NSM_CSC;
//     }
//     else
//       assert(0 && "NM_clear_cone_matrix_H supports only NSM_TRIPLET and NSM_CSC, or unknown origin");
//   }

//   default:
//   {
//     assert(0 && "NM_clear_cone_matrix_H supports only NM_SPARSE, or unknown storageType");
//   }
//   }
// }

// /* This function is to extract some rows and columns of a matrix. The new one has a reduced size. For example:
//  *
//  *     [ 1 0 3 0 ]                             [ 1 3 ]
//  * A = [ 0 9 0 0 ]  => Ac = A(:, column 0 2) = [ 0 0 ]
//  *     [ 0 0 0 4 ]                             [ 0 0 ]
//  * An allocation is done. All matrices are stored as sparse.
//  */
// NumericsMatrix * NM_extract(NumericsMatrix *A, int n_rows, int *target_rows, int n_cols, int *target_cols)
// {
//   NumericsMatrix * Ac = NULL;
//   switch(A->storageType)
//   {
//   case NM_SPARSE:
//   {
//     assert(A->matrix2);
//     NM_triplet(A);
//     if (A->matrix2->origin == NSM_CSC) NM_clearCSC(A);
//     if (A->matrix2->origin == NSM_CSR) NM_clearCSR(A);
//     A->matrix2->origin= NSM_TRIPLET;

//     CSparseMatrix *A_triplet = A->matrix2->triplet;
//     cs_long_t *A_rows = A_triplet->i, *A_cols = A_triplet->p;
//     double *A_vals = A_triplet->x;
//     CS_INT A_nz = A_triplet->nz, Ac_nz = 0;
//     int out_of_size = 0;

//     Ac = NM_create(NM_SPARSE, n_rows, n_cols);

//     // Count the number of non-zero elements in the target rows and columns
//     for (int i=0; i<n_rows; i++)
//       for (int j=0; j<n_cols; j++)
//         for (int k=0; k<A_nz; k++)
//         {
//           if (A_rows[k] == target_rows[i] && A_cols[k] == target_cols[j])
//             Ac_nz++;
//           // if (target_rows[i] > A_rows[k] || target_cols[j] > A_cols[k])
//           //   out_of_size = 1;
//         }



//     // if (out_of_size)
//     // {
//     //   assert(0 && "NM_extract: Target rows or columns are out of size of the matrix!");
//     //   return NULL;
//     // }

//     if (Ac_nz == 0) return NULL;

//     // Create the compressed matrix
//     CSparseMatrix *Ac_triplet = cs_spalloc(n_rows, n_cols, Ac_nz, 1, 1);
//     cs_long_t *Ac_rows = Ac_triplet->i, *Ac_cols = Ac_triplet->p;
//     double *Ac_vals = Ac_triplet->x;
//     int index = 0;

//     for (int i=0; i<n_rows; i++)
//       for (int j=0; j<n_cols; j++)
//         for (int k=0; k<A_nz; k++)
//           if (A_rows[k] == target_rows[i] && A_cols[k] == target_cols[j])
//           {
//             Ac_rows[index] = i;
//             Ac_cols[index] = j;
//             Ac_vals[index++] = A_vals[k];
//             // memcpy(Ac_vals+index, A_vals+k, (size_t)sizeof(double));
//             // index++;
//           }

//     assert(index!=Ac_nz && "NM_extract: There is an error in the routine for inserting elements to the new matrix");

//     Ac->matrix2->triplet = Ac_triplet;
//     Ac->matrix2->origin = NSM_TRIPLET;
//     Ac->matrix2->triplet->nz = Ac_nz;

//     // printf("NM_extract: check Ac = \n");
//     // for (int i=0; i<index; i++)
//     //   printf("%2lld %2lld: %e\n", Ac_rows[i], Ac_cols[i], Ac_vals[i]);
//   }

//   default:
//   {
//     assert(0 && "NM_extract supports only NM_SPARSE, or unknown storageType");
//   }
//   }

//   return Ac;
// }



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

  // change of variable to eliminate the friction coefficients: H_tilde --> H and w_tilde --> w
  NumericsMatrix *H = NM_multiply(P_mu, H_tilde);
  NM_gemv(1.0, P_mu, w_tilde, 0.0, w);

  // FILE *file = NULL;
  // FILE *file = fopen("Hmat.m", "w");
  // fprintf(file,"H = [\n");
  // CSparseMatrix_print_in_Matlab_file(NM_triplet(H), 0, file);
  // fprintf(file,"];\n");
  // fprintf(file,"H = sparse(int32(H(:,1)), int32(H(:,2)), H(:,3));\n");

  // fprintf(file,"M = [\n");
  // CSparseMatrix_print_in_Matlab_file(NM_triplet(M), 0, file);
  // fprintf(file,"];\n");
  // fprintf(file,"M = sparse(int32(M(:,1)), int32(M(:,2)), M(:,3));\n");

  // int target_rows[] = {0, 1, 2, 36, 37, 38};
  // int target_cols[] = {0, 1, 2, 3, 4, 5};
  // NumericsMatrix *Hc = NM_extract(H, 6, target_rows, 6, target_cols);

  // printf("H = \n"); NM_display(H);
  // printf("\n\nM = \n"); NM_display(M);
  // printf("\n\nHc = \n"); NM_display(Hc);

  // int cones_to_clear[] = {28};
  // NM_clear_cone_matrix_H(H, 5, cones_to_clear);
  // n -= 5; nd = n*d;
  // int cones_to_clear[] = {337,352,356,384,400};
  // int n_cones_to_clear = 1;
  // NM_clear_cone_matrix_H(H, n_cones_to_clear, cones_to_clear);


  // fprintf(file,"Hc = [\n");
  // CSparseMatrix_print_in_Matlab_file(NM_triplet(Hc), 0, file);
  // fprintf(file,"];\n");
  // fprintf(file,"Hc = sparse(int32(Hc(:,1)), int32(Hc(:,2)), Hc(:,3));\n");
  // fclose(file);

  // // Delete the corresponding cones in w
  // int assign_w = 1, count_w = 0;
  // double *w_tmp = (double*)calloc(nd, sizeof(double));
  // for (int i=0; i<nd; i+=d)
  // {
  //   assign_w = 1;
  //   for(int j=0; j<n_cones_to_clear; j++)
  //   {
  //     if(i/3 == cones_to_clear[j]) {assign_w = 0; break;}
  //   }

  //   if (assign_w)
  //   {
  //     w_tmp[count_w++] = w[i];
  //     w_tmp[count_w++] = w[i+1];
  //     w_tmp[count_w++] = w[i+2];
  //   }
  // }

  // n -= n_cones_to_clear; nd = n*d;
  // NV_copy(w_tmp, nd, w);
  // free(w_tmp);


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
  double udotr = 1e300, uor_mu = 1e300;
  double projerr = 1e300;
  double error[6];
  double totalresidual = 1e300, totalresidual_mu = 1e300;

  double diff_fixp = 1e300, nub = 1e300;
  double *diff_fixp_vec = (double*)calloc(n,sizeof(double));

  double *primalConstraint = data->tmp_vault_nd[1];
  double *dualConstraint = data->tmp_vault_m[0];
  double *complemConstraint = data->tmp_vault_nd[2];
  double *complemConstraint_mu = (double*)calloc(nd,sizeof(double));
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


  double *u_plus_du = data->tmp_vault_nd[3];  // for Mehrotra
  double *r_plus_dr = data->tmp_vault_nd[4];  // for Mehrotra
  double *dudr_jprod = data->tmp_vault_nd[5];  // for Mehrotra

  double *rhs = options->dWork;
  double *rhs_2 = (double*)calloc(m+2*nd+n, sizeof(double));
  double *sol = (double*)calloc(m+2*nd+n, sizeof(double));

  char fws = ' '; /* finish without scaling */

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
  // numerics_printf_verbose(-1, "problem dimensions n, nd x m: %1i, %6i x %-6i",n, nd, m);
      switch ( options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] )
    {
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_NOSCAL:
    {
      numerics_printf_verbose(-1,"Global friction contact problem - LS solution: 4x4 no scaling\n");
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_QP2:
    {
      numerics_printf_verbose(-1,"Global friction contact problem - LS solution: 4x4 NT scaling with Qp2\n");
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST1:
    {
      numerics_printf_verbose(-1,"Global friction contact problem - TESTING: 4x4 no scaling + equation: s^2 = |ub|^2 \n");
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST2:
    {
      numerics_printf_verbose(-1,"Global friction contact problem - TESTING: 4x4 no scaling + Mehrotra \n");
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST3:
    {
      numerics_printf_verbose(-1,"Global friction contact problem - TESTING: 4x4 no scaling + backtracking line-search for |uor| and |s-ub|\n");
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST4:
    {
      numerics_printf_verbose(-1,"Global friction contact problem - TESTING: 4x4 no scaling + backtracking line-search for function Theta\n");
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST5:
    {
      numerics_printf_verbose(-1,"Global friction contact problem - TESTING: 4x4 no scaling + non-monotone line search for function Theta\n");
      break;
    }
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_TEST6:
    {
      numerics_printf_verbose(-1,"Global friction contact problem - TESTING: 3x3 no scaling (model without s) \n");
      break;
    }
    default:
    {
      printf("ERROR\n");
    }
    }

  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST3)
  {
    numerics_printf_verbose(-1, "| it  | pinfeas | dinfeas |  |s-ub| |  |uor|  |  u'r/n  | prj err | barpram | alpha 1 | alpha 2 | alpha 3 |  |dv|   |  |du|   |  |dr|   |  |ds|   | ls prim | ls dual | ls comp | ls fixP |");
    numerics_printf_verbose(-1, "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
  }
  else if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST4)
  {
    numerics_printf_verbose(-1, "| it  | pinfeas | dinfeas |  |s-ub| |  |uor|  | |Theta| |  u'r/n  | prj err | barpram | alpha 1 | alpha 2 |  |dv|   |  |du|   |  |dr|   |  |ds|   | ls prim | ls dual | ls comp | ls fixP |");
    numerics_printf_verbose(-1, "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
  }
  else if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_TEST6)
  {
    numerics_printf_verbose(-1, "| it  | pinfeas | dinfeas |  |uor|  |  u'r/n  | prj err | barpram | alpha_r | alpha_u |  |dv|   |  |du|   |  |dr|   | ls prim | ls dual | ls comp |");
    numerics_printf_verbose(-1, "---------------------------------------------------------------------------------------------------------------------------------------------------");

  }
  else
  {
    numerics_printf_verbose(-1, "| it  | pinfeas | dinfeas |  |s-ub| | |u o r| |   u'r   | prj err | barpram |  sigma  ||  alpha  |  |dv|   |  |du|   |  |dr|   |  |ds|   | ls prim | ls dual | ls comp | ls fixP |");
    numerics_printf_verbose(-1, "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
  }

  double * p2 = data->tmp_vault_nd[9];
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
  // char *str = (char *) malloc(200);
  // strcpy( str, problem_name );
  // const char * separators = "/";
  // char *strToken = strtok( str, separators );
  // for(int i=0; i<5; i++)
  // {
  //   if(strToken != NULL) strToken = strtok ( NULL, separators );
  // }

  // strToken = strtok ( strToken, "." );
  // for(int i=0; i<strlen(strToken); i++)
  // {
  //   if(strToken[i] == '-') strToken[i] = '_';
  // }

  char matlab_name[100], probName[100];

  // int count=0; for (int i = 13; problem_name[i] != '.'; i++) {probName[count] = problem_name[i]; count++;} probName[count] = '\0';
  // sprintf(matlab_name, "%s.m",probName);

  // sprintf(matlab_name, "%s.m",strToken);
  sprintf(matlab_name, "iterates_Spheres_no_s.m");

  /* writing data in a Matlab file */
  if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
  {
    // iterates_2 = fopen("box466_s_u0_ub_noSubDiff.m", "w");
    iterates = fopen(matlab_name, "a+");
    // iterates = fopen(matlab_name, "w");
    fprintf(iterates,"%% data = struct;\n");
    fprintf(iterates,"data(end+1).name = \"%s\";\n", strToken);
    fprintf(iterates,"data(end).val = [\n");
    // printDataProbMatlabFile(M, f, H, w, d, n, m, problem->mu, iterates_2);
  }

  /* check the full criterion */
  double norm_q = cblas_dnrm2(m, problem->q, 1);
  double norm_b = cblas_dnrm2(nd, problem->b, 1);

  double kappa_eps = 0.1, kappa_mu = 0.3;
  double tmp_barr_param = 0.;
  double max_uor_2mu = 0., tmp_uor_2mu = 0.;
  int findParam = 1;
  double scale_sub_diff = 0.91; //1.05;


  ComputeErrorGlobalPtr computeError = NULL;
  computeError = (ComputeErrorGlobalPtr)&gfc3d_compute_error;

  // ComputeErrorPtr = computeError = NULL;
  // computeError = (ComputeErrorPtr)&gfc3d_compute_error;


  int load_starting_point = 0, save_sol_point = 0, pertu_point = 0;

// while(hasNotConverged != 0 && findParam)
while(findParam)
{
  findParam = 0;
  // Reset vars
  // r
  for (unsigned int  i = 0; i<nd; i++)
    if (i % d == 0) reaction[i] = 0.1;
    else reaction[i] = 0.01;
  // v
  for (unsigned int  i = 0; i<m; i++) globalVelocity[i] = f[i];
  NM_tgemv(1.0, H, reaction, 1.0, globalVelocity);
  NM_Cholesky_solve(M, globalVelocity, 1);
  // for (unsigned int  i = 0; i<m; i++) globalVelocity[i] = 0.02;
  // u
  for (unsigned int  i = 0; i<nd; i++)
      if (i % d == 0) velocity[i] = 0.1;
      else velocity[i] = 0.01;
  // s
  for (unsigned int  i = 0; i<n; i++)
  {
    nub = cblas_dnrm2(2, velocity+i*d+1, 1);
    s[i] = nub;
  }

  if (findParam && options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE] && iteration > 0)
  {
    fprintf(iterates,"%% data = struct;\n");
    fprintf(iterates,"data(end+1).name = \"%s\";\n", strToken);
    fprintf(iterates,"data(end).val = [\n");
  }

  if (load_starting_point)
  {
    // sol_file = fopen("sol_data.res", "r");
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
  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] = SICONOS_FRICTION_3D_IPM_IPARAM_LS_2X2_QPH_TEST;
  // gfc3d_IPM_fixed(problem, reaction, velocity, globalVelocity, info, options, 1e-4, problem_name);
  // for(unsigned int i=0; i<n; i+=d)
  // {
  //   printf("\n%3i-%3i: r = ", i, i+d-1);
  //   for(unsigned int j=0; j<d; j++) printf(" %.2e,", reaction[i+j]);
  //   printf("\tr0-|rb| = %.2e", reaction[i]-cblas_dnrm2(2, reaction+i+1, 1));
  //   printf("\t\t u = ");
  //   for(unsigned int j=0; j<d; j++) printf(" %.2e,", velocity[i+j]);
  //   printf("\tu0-|ub| = %.2e", velocity[i]-cblas_dnrm2(2, velocity+i+1, 1));
  // }
  // printf("\n");
  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] = SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_NOSCAL;
  // numerics_printf_verbose(-1, "| it  | pinfeas | dinfeas |  |s-ub| | |u o r| |  u'r/n  | prj err | barpram |  sigma  ||  alpha  |  |dv|   |  |du|   |  |dr|   |  |ds|   | ls prim | ls dual | ls comp | ls fixP |");
  // numerics_printf_verbose(-1, "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");


  // Reset params
  iteration = 0;
  hasNotConverged = 1;
  pinfeas = dinfeas = complem = udotr = projerr = diff_fixp = totalresidual = 1e300;
  alpha_primal = alpha_dual = 1.;
  barr_param = 1.;
  sigma = 0.1;

  // if (iteration == 0)
  // {
  //   primalResidual_s(velocity, H, globalVelocity, w, s, primalConstraint, &pinfeas, tol);
  //   dualResidual(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol);
  //   JA_prod(velocity, reaction, nd, n, complemConstraint);
  //   for (int k = 0; k < nd; complemConstraint[k] -= 2*barr_param, k+=d);
  // }

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
      // if (iteration == 0)
      // {
      //   printf("H = \n"); NM_display(H);
      //   NumericsMatrix *Hs = NM_clear_zero(H, 1e-13);
      //   // NM_display_tol(H, 1e-13);
      //   printf("\nHs = \n"); NM_display(Hs);

      // }


      if (load_starting_point && iteration == 0)
      {
        double bdr = -1., bdu = -1., nr = -1, nu = -1;
        for (unsigned int i=0; i<nd; i+=d)
        {
          // if (velocity[i] < 0.) velocity[i] *= -1.; // if u0 < 0 => u0 >0
          // if (reaction[i] < 0.) reaction[i] *= -1.; // if u0 < 0 => u0 >0

          // nub = cblas_dnrm2(2, velocity+i+1, 1);
          // s[i/d] = nub;
          // velocity[i] += nub;
        }

        printf("\nBEFORE:");
        for(unsigned int i=0; i<n; i+=d)
        {
          printf("\n%3i-%3i: r = ", i, i+d-1);
          for(unsigned int j=0; j<d; j++) printf(" %.2e,", reaction[i+j]);
          printf("\tr0-|rb| = %.2e", reaction[i]-cblas_dnrm2(2, reaction+i+1, 1));
          printf("\t\t u = ");
          for(unsigned int j=0; j<d; j++) printf(" %.2e,", velocity[i+j]);
          printf("\tu0-|ub| = %.2e", velocity[i]-cblas_dnrm2(2, velocity+i+1, 1));
        }
        printf("\n");

        double scale_zero = 1e0, scale_increase = 1. + 1e-4;
        for (unsigned int i=0; i<nd; i+=d)
        {
          // pre-processing input data
          nu = cblas_dnrm2(d, velocity+i, 1);
          nr = cblas_dnrm2(d, reaction+i, 1);
          bdr = reaction[i] - cblas_dnrm2(2, reaction+i+1, 1);
          bdu = velocity[i] - cblas_dnrm2(2, velocity+i+1, 1);

          if (bdr < 0.) reaction[i] += fabs(bdr) + 1e-15;
          if (bdu < 0.) velocity[i] += fabs(bdu) + 1e-15;

          bdr = reaction[i] - cblas_dnrm2(2, reaction+i+1, 1);
          bdu = velocity[i] - cblas_dnrm2(2, velocity+i+1, 1);
          // // tol = 1e-2;
          if (bdr < tol/scale_zero/10.) reaction[i] *= scale_increase;
          if (bdu < tol/scale_zero/10.) velocity[i] *= scale_increase;

          bdr = reaction[i] - cblas_dnrm2(2, reaction+i+1, 1);
          bdu = velocity[i] - cblas_dnrm2(2, velocity+i+1, 1);

          if (nu == 0.)
          {
            velocity[i] = tol/scale_zero; velocity[i+1] = velocity[i]/100.; velocity[i+2] = velocity[i]/1000.;
            if (bdr < tol/scale_zero/10.) reaction[i] *= scale_increase;
          }

          if (nr == 0.)
          {
            reaction[i] = tol/scale_zero; reaction[i+1] = reaction[i]/100.; reaction[i+2] = reaction[i]/1000.;
            if (bdu < tol/scale_zero/10.) velocity[i] *= scale_increase;
          }

          if (nu != 0. && nu < tol/scale_zero/10. && bdr > tol*scale_zero*10.)
          {
            velocity[i] = tol/scale_zero; velocity[i+1] = velocity[i]/100.; velocity[i+2] = velocity[i]/1000.;
            reaction[i] *= scale_increase;
          }
          else if (nr != 0. && nr < tol/scale_zero/10. && bdu > tol*scale_zero*10.)
          {
            reaction[i] = tol/scale_zero; reaction[i+1] = reaction[i]/100.; reaction[i+2] = reaction[i]/1000.;
            velocity[i] *= scale_increase;
          }
          else if (nu != 0. && nr != 0. && bdr < tol && bdu < tol)
          {
            reaction[i] *= scale_increase; // for r in bd(L), increase r0 by a very small value
            velocity[i] *= scale_increase;

          }

          nub = cblas_dnrm2(2, velocity+i+1, 1);
          s[i/d] = nub*scale_increase;
          velocity[i] += nub;
        }

        printf("\nAFTER:");
        for(unsigned int i=0; i<n; i+=d)
        {
          printf("\n%3i-%3i: r = ", i, i+d-1);
          for(unsigned int j=0; j<d; j++) printf(" %.2e,", reaction[i+j]);
          printf("\tr0-|rb| = %.2e", reaction[i]-cblas_dnrm2(2, reaction+i+1, 1));
          printf("\t\t u = ");
          for(unsigned int j=0; j<d; j++) printf(" %.2e,", velocity[i+j]);
          printf("\tu0-|ub| = %.2e", velocity[i]-cblas_dnrm2(2, velocity+i+1, 1));
        }
        printf("\n");
      }
      // tol = 1e-10;
      // Check Stopping test
      /* Computation of the values of
       - primal residual: u - H*v - w - phi(s)
       - dual residual: M*v - f - H'*r
       - duality gap: u'*r
       - complementarity: u o r
       - projection error: r - proj(r-u)
      */
      // printf("u = "); for (int i = 0; i<9; i++) printf(" %.2e", velocity[i]);
      // printf("\n\nr = "); for (int i = 0; i<9; i++) printf(" %.2e", reaction[i]);
      // printf("\n");
      primalResidual_s(velocity, H, globalVelocity, w, s, primalConstraint, &pinfeas, tol);
      dualResidual(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol);
      complem = complemResidualNorm(velocity, reaction, nd, n);
      udotr = cblas_ddot(nd, velocity, 1, reaction, 1);
      // barr_param = udotr /n;
      // if (iteration == 0)
      barr_param = (udotr / n); //*kappa_mu;

      JA_prod(velocity, reaction, nd, n, complemConstraint);
      cblas_dcopy(nd, complemConstraint, 1, complemConstraint_mu, 1);
      for (int k = 0; k < nd; k+=d)
      {
        complemConstraint_mu[k] -= 2*sigma*barr_param;
      }
      uor_mu = cblas_dnrm2(nd, complemConstraint_mu, 1);


      diff_fixp = 0.;
      for (unsigned int i = 0; i<nd; i+=d)
      {
        nub = cblas_dnrm2(2, velocity+i+1, 1);
        // diff_fixp += (s[i/d] - nub)*(s[i/d] - nub);
        // diff_fixp += fabs(s[i/d] - nub);
        diff_fixp_vec[i/d] = s[i/d]-nub;

      }
      // diff_fixp = sqrt(diff_fixp);
      // diff_fixp /= n;
      diff_fixp = cblas_dnrm2(n, diff_fixp_vec, 1);
      // file = fopen("smub2.m", "a+");
      // fprintf(file,"S(%3i,:) = [",iteration+1);
      // for(int i = 0; i < n; i++)
      // {
      //   fprintf(file, "%8.20e, ", diff_fixp_vec[i]);
      // }
      // fprintf(file,"];\n");
      // fclose(file);


      // compute Projection Error
      NM_gemv(1.0, P_mu_inv, velocity, 0.0, data->tmp_point->t_velocity);
      NM_gemv(1.0, P_mu, reaction, 0.0, data->tmp_point->t_reaction);
      // (*computeError)(problem,
      //                 data->tmp_point->t_reaction, data->tmp_point->t_velocity, globalVelocity,
      //                 tol, options, norm_q, norm_b, &projerr);
      // gfc3d_compute_error_r(problem,
      //                 data->tmp_point->t_reaction, data->tmp_point->t_velocity, globalVelocity,
      //                 tol, options, norm_q, norm_b, &projerr);
      projerr = projectionError(velocity, reaction, n, tol);

      totalresidual = fmax(fmax(fmax(pinfeas, dinfeas),diff_fixp),complem);
      totalresidual_mu = fmax(fmax(fmax(pinfeas, dinfeas),diff_fixp),uor_mu);

      // if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
      // {
      //   fprintf(iterates,"%d %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e;\n",
      //         iteration, pinfeas, dinfeas, diff_fixp, uor_mu, udotr, projerr, barr_param, alpha_primal,
      //         fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
      //         fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
      //         fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]),
      //         fabs(d_s[cblas_idamax(n, d_s, 1)]),
      //         LS_norm_p, LS_norm_d, LS_norm_c, LS_norm_f);
      // }


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
                              iteration, fws, pinfeas, dinfeas, diff_fixp, complem, udotr, projerr, barr_param, sigma);

        hasNotConverged = 0;
        break;
      }



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
      size_t pos;
      scale_sub_diff = 1.;
      for(size_t i = 0; i < n; ++i)
      {
        pos = i * d;
        nub = cblas_dnrm2(2, velocity+pos+1, 1);

        // if (ub > sqrt(DBL_EPSILON))    // subdiff_u = ub/|ub|
        // {
        //   NM_entry(subdiff_u, i, pos+1, -velocity[pos+1]/ub);
        //   NM_entry(subdiff_u, i, pos+2, -velocity[pos+2]/ub);
        // }

        // else      // // subdiff_u = arbitrary {x, |x| <= 1}
        // {
        //   NM_entry(subdiff_u, i, pos+1, -1/sqrt(2));
        //   NM_entry(subdiff_u, i, pos+2, -1/sqrt(2));
        //   // printf("Cone %zu: nub = 0\n",i);
        // }

        // if (nub < DBL_EPSILON)
        // {
        //   printf("Cone %zu: nub = %.1e, %.1e,\tub/|nub| = %.1e, %.1e\n",i, velocity[pos+1], velocity[pos+2], -velocity[pos+1]/nub, -velocity[pos+2]/nub);
        // }

        // if (rbd > 1e-8)
        // {
        //   scale = 0.9;
        //   // printf("\nScale!\n");
        // }
        // else scale = 1.;

        // if (rbd < 1e-5)
        // {
          NM_entry(subdiff_u, i, pos+1, -1.*scale_sub_diff*velocity[pos+1]/nub);
          NM_entry(subdiff_u, i, pos+2, -1.*scale_sub_diff*velocity[pos+2]/nub);
          // printf("\n nub = %e \n", nub);
        // }
        fixpConstraint[i] = s[i] - nub;  // fixpConstraint = s - |u_bar|
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
      // NM_insert(J, NM_scalar(nd, -1e-7), m + nd, m + nd);
      // NumericsMatrix * delta = NM_create(NM_SPARSE, n, nd);
      // size_t delta_nzmax = n;
      // NM_triplet_alloc(delta, delta_nzmax);
      // NM_fill(delta, NM_SPARSE, n, nd, delta->matrix2);
      // for(size_t i = 0; i < n; ++i)
      // {
      //   NM_entry(delta, i, i*d, -1e-6);
      // }
      // NM_insert(J, delta, m + 2*nd, m + nd);


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
      cblas_dcopy(nd, complemConstraint, 1, rhs+m, 1);
      for (int k = 0; k < nd; rhs[m+k] -= 2*sigma*barr_param, k+=d);

      cblas_dcopy(nd, primalConstraint, 1, rhs+m+nd, 1);
      cblas_dcopy(n, fixpConstraint, 1, rhs+m+2*nd, 1);



      cblas_dscal(m + 2*nd + n, -1.0, rhs, 1);
      cblas_dcopy(m + 2*nd + n, rhs, 1, rhs_2, 1);    // rhs_2 = old rhs


      // SOLVE
      // NM_LU_solve(J, rhs, 1);


      double * rhs_tmp = (double*)calloc(m+2*nd+n,sizeof(double));
      cblas_dcopy(m+2*nd+n, rhs_2, 1, rhs_tmp, 1);
      for (int k=0; k<m+2*nd+n; sol[k] = 0., k++);  // reset sol

      int max_refine = 1;
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] == SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES) max_refine = 10;

      for (int itr = 0; itr < max_refine; itr++)
      {
        NM_LU_solve(J, rhs_tmp, 1);                       // rhs_tmp = d = solution of J*d = rhs
        cblas_daxpy(m+2*nd+n, 1.0, rhs_tmp, 1, sol, 1);   // sol = x_+ = x + d
        cblas_dcopy(m+2*nd+n, rhs_2, 1, rhs_tmp, 1);      // rhs_tmp = old rhs = b
        NM_gemv(-1.0, J, sol, 1.0, rhs_tmp);              // rhs_tmp = b - J*x_+
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

      /* computing the affine step-length */
      alpha_primal = getStepLength(velocity, d_velocity, nd, n, 0.95);
      alpha_dual = getStepLength(reaction, d_reaction, nd, n, 0.95);
      alpha_primal = alpha_dual = fmin(alpha_primal, alpha_dual);
      // if (20 < iteration && iteration < 30) {alpha_primal = alpha_dual = 1.;}
      // if (load_starting_point == 1 && (iteration == 0 || iteration == 49))
      // if (totalresidual < tol*10.)
      // {
      //   double *arr_alpha_u = array_getStepLength(velocity, d_velocity, nd, n, 0.95);
      //   double *arr_alpha_r = array_getStepLength(reaction, d_reaction, nd, n, 0.95);

      //   double somme_t[3];
      //   double dur_t[3];
      //   double ns_t, ndur_t;

      //   printf("\nFINAL:");
      //   for(unsigned int i=0; i<n; i+=d)
      //   {
      //     printf("\n%3i-%3i: r = ", i, i+d-1);
      //     for(unsigned int j=0; j<d; j++) printf(" %.2e,", reaction[i+j]);
      //     printf("\tr0-|rb| = %.2e", reaction[i]-cblas_dnrm2(2, reaction+i+1, 1));
      //     printf("\t u = ");
      //     for(unsigned int j=0; j<d; j++) printf(" %.2e,", velocity[i+j]);
      //     printf("\tu0-|ub| = %.2e", velocity[i]-cblas_dnrm2(2, velocity+i+1, 1));
      //     printf("\tal_r = %e, al_u = %e", arr_alpha_r[i/d], arr_alpha_u[i/d]);

      //     somme_t[0] = velocity[i] + reaction[i];
      //     somme_t[1] = velocity[i+1] + reaction[i+1];
      //     somme_t[2] = velocity[i+2] + reaction[i+2];
      //     dur_t[0] = velocity[i] - reaction[i];
      //     dur_t[1] = velocity[i+1] - reaction[i+1];
      //     dur_t[2] = velocity[i+2] - reaction[i+2];

      //     ns_t = somme_t[0] - cblas_dnrm2(2,somme_t+1,1);
      //     ndur_t = dur_t[0] - cblas_dnrm2(2,dur_t+1,1);
      //     if (ns_t > sqrt(DBL_EPSILON)*cblas_dnrm2(3, somme_t, 1))
      //     {
      //       if (dur_t[0] >= cblas_dnrm2(2,dur_t+1,1))  printf("\t==> N");
      //       else if (-dur_t[0] >= cblas_dnrm2(2,dur_t+1,1)) printf("\t==> B");
      //       else                                        printf("\t==> R");
      //     }
      //     else
      //       printf("\t==> NOT SC");
      //   }
      //   printf("\n");








      //   // for (int i = 0; i<n; i++)
      //   // {
      //   //   printf("Cone %d: alpha_u = %e, alpha_r = %e\n", i, arr_alpha_u[i], arr_alpha_r[i]);
      //   // }
      //   // printf("\n");
      // }
      // if (alpha_primal < 1e-8)
      // {
      //   printf("\nfailure\n\n");
      //   for(unsigned int i=0; i<n; i+=d)
      //   {
      //     printf("\n%i-%i: r = ", i, i+d-1);
      //     for(unsigned int j=0; j<d; j++) printf(" %.2e,", reaction[i+j]);
      //     printf("\tr0-|rb| = %.2e", reaction[i]-cblas_dnrm2(2, reaction+i+1, 1));
      //     printf("\t\t u = ");
      //     for(unsigned int j=0; j<d; j++) printf(" %.2e,", velocity[i+j]);
      //     printf("\tu0-|ub| = %.2e", velocity[i]-cblas_dnrm2(2, velocity+i+1, 1));
      //   }
      //   printf("\n");
      //   hasNotConverged = 2;
      //   break;
      // }


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
      NM_gemv(1.0, J, sol, -1.0, rhs_2);

      LS_norm_d = cblas_dnrm2(m, rhs_2, 1);
      LS_norm_c = cblas_dnrm2(nd, rhs_2+m, 1);
      LS_norm_p = cblas_dnrm2(nd, rhs_2+m+nd, 1);
      LS_norm_f = cblas_dnrm2(n, rhs_2+m+2*nd, 1);

      cblas_dcopy(m, sol, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, sol+m, 1, d_velocity, 1);
      cblas_dcopy(nd, sol+m+nd, 1, d_reaction, 1);
      cblas_dcopy(n, sol+m+2*nd, 1, d_s, 1);

      /* computing the affine step-length */
      alpha_primal = getStepLength(velocity, d_velocity, nd, n, 0.95);
      alpha_dual = getStepLength(reaction, d_reaction, nd, n, 0.95);
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

      primalResidual_s(velocity, H, globalVelocity, w, s, primalConstraint, &pinfeas, tol);
      dualResidual(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol);
      complem = complemResidualNorm(velocity, reaction, nd, n);
      udotr = cblas_ddot(nd, velocity, 1, reaction, 1)/n;

      JA_prod(velocity, reaction, nd, n, complemConstraint);
      max_uor_2mu = 0.0;
      for (int k = 0; k < nd; k+=d)
      {
        complemConstraint[k] -= 2*barr_param;
        tmp_uor_2mu = cblas_dnrm2(3, complemConstraint+k, 1);

        if (tmp_uor_2mu > max_uor_2mu) max_uor_2mu = tmp_uor_2mu;
      }
      uor_mu = cblas_dnrm2(nd, complemConstraint, 1);


      diff_fixp = 0.;
      for (unsigned int i = 0; i<nd; i+=d)
      {
        nub = cblas_dnrm2(2, velocity+i+1, 1);
        diff_fixp += (s[i/d] - nub)*(s[i/d] - nub);
      }
      diff_fixp = sqrt(diff_fixp);

      break;
    }


    // A COMPLETER cette option
    /*  Build the Jacobian matrix without reducing the linear system.
     *
     *
     *  In the case where the NT scaling is not used, the matrix to factorize is the following:
     *
     *         nd           nd        n
     *      |   W           -I        E | nd
     *      |                           |
     *  J = | Arw(r)      Arw(u)      0 | nd
     *      |                           |
     *      |   0     [0 -ub'/|ub|]   I | n
     *
     *
     * Case of a non-reduced system. Building the right-hand side related to the first system

       without NT scaling
       rhs = -
       [      Wr + q + Es - u     ]  nd        primalConstraint   // use var primalConstraint (size nd) for this equ to take advantage of the set of initialized vars
       [      u o r - 2 mu e      ]  nd        complemConstraint
       [         s - |ub|         ]  n
    */
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_NOSCAL_REDUCED:
    {
      // First linear linear system
      J = NM_create(NM_SPARSE, 2*nd + n, 2*nd + n);
      J_nzmax = nd*nd + nd + n + 2*(d*3-2)*n + 2*n + n;
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
      scale_sub_diff = 0.95;
      for(size_t i = 0; i < n; ++i)
      {
        pos = i * d;
        ub = sqrt(velocity[pos+1]*velocity[pos+1]+velocity[pos+2]*velocity[pos+2]);

        NM_entry(subdiff_u, i, pos+1, -1.*scale_sub_diff*velocity[pos+1]/ub);
        NM_entry(subdiff_u, i, pos+2, -1.*scale_sub_diff*velocity[pos+2]/ub);

        fixpConstraint[i] = s[i] - ub;  // fixpConstraint = s - |u_bar|
      }

      // W = H*M^-1*H'
      NumericsMatrix *MmHT = NM_transpose(H);

      NM_insert(J, M, 0, 0);
      NM_insert(J, minus_H, m + nd, 0);

      NM_insert(J, arrow_r, m, m);
      NM_insert(J, eye_nd, m + nd, m);
      NM_insert(J, subdiff_u, m + 2*nd, m);

      NM_insert(J, minus_Ht, 0, m + nd);
      NM_insert(J, arrow_u, m, m + nd);

      NM_insert(J, minus_e, m + nd, m + 2*nd);
      NM_insert(J, eye_n, m + 2*nd, m + 2*nd);

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


      // SOLVE
      // NM_LU_solve(J, rhs, 1);


      double * rhs_tmp = (double*)calloc(m+2*nd+n,sizeof(double));
      cblas_dcopy(m+2*nd+n, rhs_2, 1, rhs_tmp, 1);
      for (int k=0; k<m+2*nd+n; sol[k] = 0., k++);  // reset sol

      int max_refine = 1;
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] == SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES) max_refine = 10;

      for (int itr = 0; itr < max_refine; itr++)
      {
        NM_LU_solve(J, rhs_tmp, 1);                       // rhs_tmp = d = solution of J*d = rhs
        cblas_daxpy(m+2*nd+n, 1.0, rhs_tmp, 1, sol, 1);   // sol = x_+ = x + d
        cblas_dcopy(m+2*nd+n, rhs_2, 1, rhs_tmp, 1);      // rhs_tmp = old rhs = b
        NM_gemv(-1.0, J, sol, 1.0, rhs_tmp);              // rhs_tmp = b - J*x_+
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

      /* computing the affine step-length */
      alpha_primal = getStepLength(velocity, d_velocity, nd, n, 0.95);
      alpha_dual = getStepLength(reaction, d_reaction, nd, n, 0.95);
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

      primalResidual_s(velocity, H, globalVelocity, w, s, primalConstraint, &pinfeas, tol);
      dualResidual(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol);
      complem = complemResidualNorm(velocity, reaction, nd, n);
      udotr = cblas_ddot(nd, velocity, 1, reaction, 1)/n;

      JA_prod(velocity, reaction, nd, n, complemConstraint);
      max_uor_2mu = 0.0;
      for (int k = 0; k < nd; k+=d)
      {
        complemConstraint[k] -= 2*barr_param;
        tmp_uor_2mu = cblas_dnrm2(3, complemConstraint+k, 1);

        if (tmp_uor_2mu > max_uor_2mu) max_uor_2mu = tmp_uor_2mu;
      }
      uor_mu = cblas_dnrm2(nd, complemConstraint, 1);


      diff_fixp = 0.;
      for (unsigned int i = 0; i<nd; i+=d)
      {
        nub = cblas_dnrm2(2, velocity+i+1, 1);
        diff_fixp += (s[i/d] - nub)*(s[i/d] - nub);
      }
      diff_fixp = sqrt(diff_fixp);

      break;
    }




    /*  This attemps is for TESTING.
     *  Its description is detailed in the solver printing routine.
     *  To reach this, need to search SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST1 until the place containing numerics_printf_verbose(-1,"TESTING: ...
     *
     *  Build the Jacobian matrix without reducing the linear system.
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
     *      |  0     [0 -ub']      0      S | n
     *
     *
     * Case of a non-reduced system. Building the right-hand side related to the first system

       without NT scaling
       rhs = -
       [      M*v - H'*r - f      ]  m         dualConstraint
       [      u o r - 2 mu e      ]  nd        complemConstraint
       [     u - H*v - w - Es     ]  nd        primalConstraint
       [    0.5 (s^2 - |ub|^2)    ]  n         fixpConstraint
    */
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST1:
    {
      // First linear linear system
      J = NM_create(NM_SPARSE, m + 2*nd + n, m + 2*nd + n);
      J_nzmax = M_nzmax + H_nzmax + 2*(d*3-2)*n + H_nzmax + nd + n + 2*n + n;
      NM_triplet_alloc(J, J_nzmax);
      J->matrix2->origin = NSM_TRIPLET;

      NumericsMatrix * arrow_r = Arrow_repr(reaction, nd, n);
      NumericsMatrix * arrow_u = Arrow_repr(velocity, nd, n) ;

      // Create matrix
      //          [ 0 -ub'_1               ]
      // mat_ub = [           0 -ub'_2     ]
      //          [                     ...]
      mat_ub = NM_create(NM_SPARSE, n, nd);
      size_t mat_ub_nzmax = 2*n;
      NM_triplet_alloc(mat_ub, mat_ub_nzmax);
      NM_fill(mat_ub, NM_SPARSE, n, nd, mat_ub->matrix2);

      /* Matrix filling */
      size_t pos;
      for(size_t i = 0; i < n; ++i)
      {
        pos = i * d;

        NM_entry(mat_ub, i, pos+1, -1.*velocity[pos+1]);
        NM_entry(mat_ub, i, pos+2, -1.*velocity[pos+2]);

        // fixpConstraint = 0.5 * (s^2 - |u_bar|^2)
        fixpConstraint[i] = 0.5*(s[i]*s[i] - velocity[pos+1]*velocity[pos+1] - velocity[pos+2]*velocity[pos+2]);
      }


      // Create matrix mat_S = diag(s_i)
      mat_S = NM_create(NM_SPARSE, n, n);
      size_t mat_S_nzmax = n;
      NM_triplet_alloc(mat_S, mat_S_nzmax);
      NM_fill(mat_S, NM_SPARSE, n, n, mat_S->matrix2);
      for(size_t i = 0; i < n; ++i)
      {
        NM_entry(mat_S, i, i, s[i]);
      }

      NM_insert(J, M, 0, 0);
      NM_insert(J, minus_H, m + nd, 0);

      NM_insert(J, arrow_r, m, m);
      NM_insert(J, eye_nd, m + nd, m);
      NM_insert(J, mat_ub, m + 2*nd, m);

      NM_insert(J, minus_Ht, 0, m + nd);
      NM_insert(J, arrow_u, m, m + nd);

      NM_insert(J, minus_e, m + nd, m + 2*nd);
      NM_insert(J, mat_S, m + 2*nd, m + 2*nd);


      if(arrow_r) { NM_free(arrow_r); arrow_r = NULL; }
      if(arrow_u) { NM_free(arrow_u); arrow_u = NULL; }
      if(mat_ub)  { NM_free(mat_ub); mat_ub = NULL; }
      if(mat_S)   { NM_free(mat_S); mat_S = NULL; }



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


      // SOLVE
      // NM_LU_solve(J, rhs, 1);


      double * rhs_tmp = (double*)calloc(m+2*nd+n,sizeof(double));
      cblas_dcopy(m+2*nd+n, rhs_2, 1, rhs_tmp, 1);
      for (int k=0; k<m+2*nd+n; sol[k] = 0., k++);  // reset sol

      int max_refine = 1;
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

      /* computing the affine step-length */
      alpha_primal = getStepLength(velocity, d_velocity, nd, n, 0.95);
      alpha_dual = getStepLength(reaction, d_reaction, nd, n, 0.95);
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

      primalResidual_s(velocity, H, globalVelocity, w, s, primalConstraint, &pinfeas, tol);
      dualResidual(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol);
      complem = complemResidualNorm(velocity, reaction, nd, n);
      udotr = cblas_ddot(nd, velocity, 1, reaction, 1)/n;

      JA_prod(velocity, reaction, nd, n, complemConstraint);
      max_uor_2mu = 0.0;
      for (int k = 0; k < nd; k+=d)
      {
        complemConstraint[k] -= 2*barr_param;
        tmp_uor_2mu = cblas_dnrm2(3, complemConstraint+k, 1);

        if (tmp_uor_2mu > max_uor_2mu) max_uor_2mu = tmp_uor_2mu;
      }
      uor_mu = cblas_dnrm2(nd, complemConstraint, 1);


      diff_fixp = 0.;
      for (unsigned int i = 0; i<nd; i+=d)
      {
        diff_fixp += 0.25 * pow(s[i/d]*s[i/d] - velocity[i+1]*velocity[i+1] - velocity[i+2]*velocity[i+2], 2.);
      }
      diff_fixp = sqrt(diff_fixp);

      break;
    }

    /*  This attemps is for TESTING.
     *  Its description is detailed in the solver printing routine.
     *  To reach this, need to search SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST2 until the place containing numerics_printf_verbose(-1,"TESTING: ...
     *
     *  Build the Jacobian matrix without reducing the linear system.
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
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST2:
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
      scale_sub_diff = 0.95;

      for(size_t i = 0; i < n; ++i)
      {
        pos = i * d;
        ub = sqrt(velocity[pos+1]*velocity[pos+1]+velocity[pos+2]*velocity[pos+2]);

        NM_entry(subdiff_u, i, pos+1, -1.*scale_sub_diff*velocity[pos+1]/ub);
        NM_entry(subdiff_u, i, pos+2, -1.*scale_sub_diff*velocity[pos+2]/ub);

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


      if(arrow_r) { NM_free(arrow_r); arrow_r = NULL; }
      if(arrow_u) { NM_free(arrow_u); arrow_u = NULL; }
      if(subdiff_u) { NM_free(subdiff_u); subdiff_u = NULL; }



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
      cblas_dcopy(m + 2*nd + n, rhs, 1, rhs_2, 1);    // rhs_2 = old rhs

      // Solve
      double * rhs_tmp = (double*)calloc(m+2*nd+n,sizeof(double));
      cblas_dcopy(m+2*nd+n, rhs_2, 1, rhs_tmp, 1);
      for (int k=0; k<m+2*nd+n; sol[k] = 0., k++);  // reset sol

      int max_refine = 1;
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
      alpha_primal = getStepLength(velocity, d_velocity, nd, n, 1.);
      alpha_dual = getStepLength(reaction, d_reaction, nd, n, 1.);

      if (alpha_primal < alpha_dual) alpha_dual = alpha_primal;
      else alpha_primal = alpha_dual;

      gmm = gmmp1 + gmmp2 * alpha_primal;

      /* ----- Corrector step of Mehrotra ----- */
      cblas_dcopy(nd, velocity, 1, u_plus_du, 1);
      cblas_dcopy(nd, reaction, 1, r_plus_dr, 1);
      cblas_daxpy(nd, alpha_primal, d_velocity, 1, u_plus_du, 1);
      cblas_daxpy(nd, alpha_dual, d_reaction, 1, r_plus_dr, 1);

      barr_param_a = cblas_ddot(nd, u_plus_du, 1, r_plus_dr, 1) / n / 3;
      e = barr_param > sgmp1 ? fmax(1.0, sgmp2 * pow(alpha_primal,2)) : sgmp3;
      sigma = fmin(1.0, pow(barr_param_a / barr_param, e));


      // Compute rhs for 2nd linear system
      cblas_dcopy(m+2*nd+n, rhs_2, 1, rhs, 1);    // Get back value for rhs

      JA_prod(d_velocity, d_reaction, nd, n, dudr_jprod);
      cblas_daxpy(nd, -1.0, dudr_jprod, 1, rhs + m, 1);

      barr_param = cblas_ddot(nd, velocity, 1, reaction, 1) / n / 3;
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

      /* computing the affine step-length */
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

      primalResidual_s(velocity, H, globalVelocity, w, s, primalConstraint, &pinfeas, tol);
      dualResidual(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol);
      complem = complemResidualNorm(velocity, reaction, nd, n);
      udotr = cblas_ddot(nd, velocity, 1, reaction, 1)/n;

      JA_prod(velocity, reaction, nd, n, complemConstraint);
      max_uor_2mu = 0.0;
      for (int k = 0; k < nd; k+=d)
      {
        // complemConstraint[k] -= 2*barr_param;
        tmp_uor_2mu = cblas_dnrm2(3, complemConstraint+k, 1);

        if (tmp_uor_2mu > max_uor_2mu) max_uor_2mu = tmp_uor_2mu;
      }
      uor_mu = cblas_dnrm2(nd, complemConstraint, 1);


      diff_fixp = 0.;
      for (unsigned int i = 0; i<nd; i+=d)
      {
        nub = cblas_dnrm2(2, velocity+i+1, 1);
        diff_fixp += (s[i/d] - nub)*(s[i/d] - nub);
      }
      diff_fixp = sqrt(diff_fixp);

      break;
    }



    /*  This attemps is for TESTING 3.
     *  Its description is detailed in the solver printing routine.
     *  To reach this, need to search SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST3 until the place containing numerics_printf_verbose(-1,"TESTING: ...
     *
     *  Build the Jacobian matrix without reducing the linear system.
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
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST3:
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST4:
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST5:
    {
      // if (iteration == 0)
      //   barr_param = (cblas_ddot(nd, velocity, 1, reaction, 1) / n) * 0.3;

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
      scale_sub_diff = 1.;
      for(size_t i = 0; i < n; ++i)
      {
        pos = i * d;
        ub = sqrt(velocity[pos+1]*velocity[pos+1]+velocity[pos+2]*velocity[pos+2]);

        NM_entry(subdiff_u, i, pos+1, -1.*scale_sub_diff*velocity[pos+1]/ub);
        NM_entry(subdiff_u, i, pos+2, -1.*scale_sub_diff*velocity[pos+2]/ub);
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


      // SOLVE
      // NM_LU_solve(J, rhs, 1);


      double * rhs_tmp = (double*)calloc(m+2*nd+n,sizeof(double));
      cblas_dcopy(m+2*nd+n, rhs_2, 1, rhs_tmp, 1);
      for (int k=0; k<m+2*nd+n; sol[k] = 0., k++);  // reset sol

      int max_refine = 1;
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] == SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES) max_refine = 10;

      for (int itr = 0; itr < max_refine; itr++)
      {
        NM_LU_solve(J, rhs_tmp, 1);                       // rhs_tmp = d = solution of J*d = rhs
        cblas_daxpy(m+2*nd+n, 1.0, rhs_tmp, 1, sol, 1);   // sol = x_+ = x + d
        cblas_dcopy(m+2*nd+n, rhs_2, 1, rhs_tmp, 1);      // rhs_tmp = old rhs = b
        NM_gemv(-1.0, J, sol, 1.0, rhs_tmp);              // rhs_tmp = b - J*x_+
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

      /* computing the affine step-length */
      alpha_primal = getStepLength(velocity, d_velocity, nd, n, 0.95);
      alpha_dual = getStepLength(reaction, d_reaction, nd, n, 0.95);
      alpha_primal = alpha_dual = fmin(alpha_primal, alpha_dual);

      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST3 && iteration > 1)
      {
        alpha_complem = backtrack_linesearch(check_merite_function_complem, 1e-6, alpha_primal, 1e-4, 7, velocity, d_velocity, reaction, d_reaction, nd, n, 0.3);

        alpha_primal = fmin(alpha_primal, alpha_complem);

        alpha_diffixP = alpha_primal;
        double alpha_diffixP_tmp = alpha_primal;
        // for (size_t i = 0; i<n; i++)
        // {
        //   // printf("Cone %zu: ", i);
        //   alpha_diffixP_tmp = backtrack_linesearch(check_merite_function_diffixP, 1./10., alpha_diffixP_tmp, 0.01, 6, velocity+i*d, d_velocity+i*d, s+i, d_s+i, 3, 1);
        //   if (alpha_diffixP_tmp < alpha_diffixP) alpha_diffixP = alpha_diffixP_tmp;
        // }
        // printf("\n");
        alpha_diffixP = backtrack_linesearch(check_merite_function_diffixP, 1e-6, alpha_primal, 1e-4, 6, velocity, d_velocity, s, d_s, nd, n);
      }
      else if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST4 && iteration > 9)
      {
        alpha_primal = backtrack_linesearch(check_merite_function_Theta, 1e-6, alpha_primal, 1e-2,
                      20, globalVelocity, d_globalVelocity, velocity, d_velocity, reaction, d_reaction, s, d_s,
                      nd, n, 0.3, pinfeas, dinfeas, complem, diff_fixp, M, H, f, w, tol);
      }
      else if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST5 && iteration > 5)
      {
        alpha_primal = backtrack_linesearch(check_merite_function_Theta_nonMonotone, 1e-6, alpha_primal, 1e-2,
                      23, globalVelocity, d_globalVelocity, velocity, d_velocity, reaction, d_reaction, s, d_s,
                      nd, n, 0.3, pinfeas, dinfeas, complem, diff_fixp, iteration, arr_norm_theta, 3, M, H, f, w, tol);
      }


      cblas_daxpy(m, alpha_primal, d_globalVelocity, 1, globalVelocity, 1);
      cblas_daxpy(nd, alpha_primal, d_velocity, 1, velocity, 1);
      cblas_daxpy(nd, alpha_primal, d_reaction, 1, reaction, 1);
      cblas_daxpy(n, alpha_primal, d_s, 1, s, 1);

      if (NV_isnan(globalVelocity, m) | NV_isnan(velocity, nd) | NV_isnan(reaction, nd) | NV_isnan(s, n))
      {
        hasNotConverged = 2;
        break;
      }

      primalResidual_s(velocity, H, globalVelocity, w, s, primalConstraint, &pinfeas, tol);
      dualResidual(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol);
      complem = complemResidualNorm(velocity, reaction, nd, n);
      udotr = cblas_ddot(nd, velocity, 1, reaction, 1)/n;

      JA_prod(velocity, reaction, nd, n, complemConstraint);
      max_uor_2mu = 0.0;
      for (int k = 0; k < nd; k+=d)
      {
        if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] != SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST3)
          complemConstraint[k] -= 2*barr_param;

        tmp_uor_2mu = cblas_dnrm2(3, complemConstraint+k, 1);

        if (tmp_uor_2mu > max_uor_2mu) max_uor_2mu = tmp_uor_2mu;
      }
      uor_mu = cblas_dnrm2(nd, complemConstraint, 1);


      diff_fixp = 0.;
      for (unsigned int i = 0; i<nd; i+=d)
      {
        nub = cblas_dnrm2(2, velocity+i+1, 1);
        // diff_fixp += (s[i/d] - nub)*(s[i/d] - nub);
        diff_fixp += fabs(s[i/d] - nub);
      }
      // diff_fixp = sqrt(diff_fixp);
      diff_fixp /= n;

      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST4)
      {
        norm_theta = sqrt(pow(pinfeas,2)+pow(dinfeas,2)+pow(complem,2)+pow(diff_fixp,2));
      }

      else if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST5)
      {
        norm_theta = sqrt(pow(pinfeas,2)+pow(dinfeas,2)+pow(complem,2)+pow(diff_fixp,2));
        arr_norm_theta[iteration] = norm_theta;
      }

      break;
    }


    /*  Build the Jacobian matrix without reducing the linear system.
     *
     *
     *  In the case where the NT scaling is not used, the matrix to factorize is the following:
     *
     *         m            nd                  nd
     *      |  M            0                  -H^T      | m
     *      |                                            |
     *  J = |  0   Arw(r)+r(0 ub'/|ub|)   Arw(u)+|ub|I   | nd
     *      |                                            |
     *      | -H            I                   0        | nd
     *
     *
     * Case of a non-reduced system. Building the right-hand side related to the first system

       without NT scaling
       rhs = -
       [      M*v - H'*r - f         ]  m         dualConstraint
       [  (u + phi(u)) o r - 2 mu e  ]  nd        complemConstraint
       [       u - H*v - w           ]  nd        primalConstraint
    */
    case SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_TEST6:
    {
      double *u_plus_phiu = (double*)calloc(nd,sizeof(double));
      for (unsigned int i=0; i<nd; i++)
      {
        if (i%d==0)
        {
          nub = cblas_dnrm2(2, velocity+i+1, 1);
          u_plus_phiu[i] = velocity[i]+nub;
        }
        else
        {
          u_plus_phiu[i] = velocity[i];
        }
      }


      // First linear linear system
      J = NM_create(NM_SPARSE, m + 2*nd, m + 2*nd);
      J_nzmax = M_nzmax + 2*H_nzmax + 2*d*d*n + nd;
      NM_triplet_alloc(J, J_nzmax);
      J->matrix2->origin = NSM_TRIPLET;

      NumericsMatrix * R_plus_ru = NULL;
      NumericsMatrix * arrow_u_plus_phiu = Arrow_repr(u_plus_phiu, nd, n) ;

      /* Create R_plus_ru
       * Arw(r) + r(0 ub'/|ub|) = [ r0  rb'+r0*ub'/|ub|  ]
       *                          [ rb      rb*ub'/|ub|  ]
       */
      R_plus_ru = NM_create(NM_SPARSE, nd, nd);
      size_t R_plus_ru_nzmax = d*d*n;
      NM_triplet_alloc(R_plus_ru, R_plus_ru_nzmax);
      NM_fill(R_plus_ru, NM_SPARSE, nd, nd, R_plus_ru->matrix2);

      /* Matrix filling */
      size_t pos;
      scale_sub_diff = 1.; //0.95;
      for(size_t i = 0; i < n; i++)
      {
        pos = i * d;
        nub = cblas_dnrm2(2, velocity+pos+1, 1);

        NM_entry(R_plus_ru, pos  , pos, reaction[pos]);
        NM_entry(R_plus_ru, pos+1, pos, reaction[pos+1]);
        NM_entry(R_plus_ru, pos+2, pos, reaction[pos+2]);

        NM_entry(R_plus_ru, pos, pos+1, reaction[pos+1]+scale_sub_diff*reaction[pos]*velocity[pos+1]/nub);
        NM_entry(R_plus_ru, pos, pos+2, reaction[pos+2]+scale_sub_diff*reaction[pos]*velocity[pos+2]/nub);

        NM_entry(R_plus_ru, pos+1, pos+1, reaction[pos]+scale_sub_diff*reaction[pos+1]*velocity[pos+1]/nub);
        NM_entry(R_plus_ru, pos+1, pos+2, scale_sub_diff*reaction[pos+1]*velocity[pos+2]/nub);
        NM_entry(R_plus_ru, pos+2, pos+1, scale_sub_diff*reaction[pos+2]*velocity[pos+1]/nub);
        NM_entry(R_plus_ru, pos+2, pos+2, reaction[pos]+scale_sub_diff*reaction[pos+2]*velocity[pos+2]/nub);
      }

      // R_plus_ru = Arrow_repr(reaction, nd, n);
      // arrow_u_plus_phiu = Arrow_repr(velocity, nd, n) ;


      NM_insert(J, M, 0, 0);
      NM_insert(J, minus_H, m + nd, 0);

      NM_insert(J, R_plus_ru, m, m);
      NM_insert(J, eye_nd, m + nd, m);

      NM_insert(J, minus_Ht, 0, m + nd);
      NM_insert(J, arrow_u_plus_phiu, m, m + nd);


      jacobian_is_nan = NM_isnan(J);
      if (jacobian_is_nan)
      {
        numerics_printf_verbose(0, "The Jacobian matrix J contains NaN");
        break;
      }

      if (iteration == 0)
      {
        barr_param = (cblas_ddot(nd, u_plus_phiu, 1, reaction, 1)/n);
        dualResidual(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol);
        primalResidual(velocity, H, globalVelocity, w, primalConstraint, &pinfeas, tol);
        JA_prod(u_plus_phiu, reaction, nd, n, complemConstraint);
        // for (int k = 0; k < nd; complemConstraint[k] -= 2.*barr_param, k+=d);
      }

      cblas_dcopy(m, dualConstraint, 1, rhs, 1);
      cblas_dcopy(nd, complemConstraint, 1, rhs+m, 1);
      for (int k = 0; k < nd; rhs[m+k] -= 2.*barr_param, k+=d);
      cblas_dcopy(nd, primalConstraint, 1, rhs+m+nd, 1);

      cblas_dscal(m + 2*nd, -1.0, rhs, 1);
      cblas_dcopy(m + 2*nd, rhs, 1, rhs_2, 1);    // rhs_2 = old rhs


      // SOLVE
      // NM_LU_solve(J, rhs, 1);


      double * rhs_tmp = (double*)calloc(m+2*nd,sizeof(double));
      cblas_dcopy(m+2*nd, rhs_2, 1, rhs_tmp, 1);
      for (int k=0; k<m+2*nd; sol[k] = 0., k++);  // reset sol

      int max_refine = 1;
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] == SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES) max_refine = 10;

      for (int itr = 0; itr < max_refine; itr++)
      {
        NM_LU_solve(J, rhs_tmp, 1);                       // rhs_tmp = d = solution of J*d = rhs
        cblas_daxpy(m+2*nd, 1.0, rhs_tmp, 1, sol, 1);   // sol = x_+ = x + d
        cblas_dcopy(m+2*nd, rhs_2, 1, rhs_tmp, 1);      // rhs_tmp = old rhs = b
        NM_gemv(-1.0, J, sol, 1.0, rhs_tmp);              // rhs_tmp = b - J*x_+
        //printf("refinement iterations = %i %8.2e\n",itr, cblas_dnrm2(m+2*nd, rhs_tmp, 1));
        if (cblas_dnrm2(m+2*nd, rhs_tmp, 1) <= 1e-14)
        {
          // printf("\nrefinement iterations = %d %8.2e\n",itr+1, cblas_dnrm2(m+2*nd+n, rhs_tmp, 1));
          // if (rhs_tmp) {free(rhs_tmp); rhs_tmp = NULL;}
          break;
        }
      }
      if (rhs_tmp) {free(rhs_tmp); rhs_tmp = NULL;}


      NM_gemv(1.0, J, sol, -1.0, rhs_2);

      LS_norm_d = cblas_dnrm2(m, rhs_2, 1);
      LS_norm_c = cblas_dnrm2(nd, rhs_2+m, 1);
      LS_norm_p = cblas_dnrm2(nd, rhs_2+m+nd, 1);

      cblas_dcopy(m, sol, 1, d_globalVelocity, 1);
      cblas_dcopy(nd, sol+m, 1, d_velocity, 1);
      cblas_dcopy(nd, sol+m+nd, 1, d_reaction, 1);

      /* computing the affine step-length */
      alpha_dual = getStepLength(reaction, d_reaction, nd, n, 0.99);

      // Step-length for that u_k0 + alpha_primal*d_k0 >= (1-tau)u_k0 for all k = 1,...,n
      min_alpha_primal = 1.; alpha_primal = gmm = 1.;
      for (size_t i = 0; i<nd; i+=d)
      {
        if (d_velocity[i] < 0.)
        {
          // if (velocity[i] < tol) d_velocity[i] = -velocity[i];
          alpha_primal = gmm = fmin(gmm, -velocity[i] / d_velocity[i]);
          // alpha_primal = fmin(1., -0.99*velocity[i] / d_velocity[i]);
        }
        // else alpha_primal = 1.;
        // min_alpha_primal = (alpha_primal < min_alpha_primal) ? alpha_primal : min_alpha_primal;


        // if (alpha_primal < 1.)
        //   printf("Cone %zu: u0 = %e, du0 = %e, u/du = %e, alpha = %e\n", i/3, velocity[i], d_velocity[i], -velocity[i]/d_velocity[i], alpha_primal);
      }
      min_alpha_primal = gmm*alpha_primal;
      alpha_primal = fmin(min_alpha_primal, alpha_dual);

      // printf("\nu = "); NV_display(velocity, nd);
      // printf("\ndu = "); NV_display(d_velocity, nd); printf("\n");

      cblas_daxpy(m, alpha_primal, d_globalVelocity, 1, globalVelocity, 1);
      cblas_daxpy(nd, alpha_primal, d_velocity, 1, velocity, 1);
      cblas_daxpy(nd, alpha_primal, d_reaction, 1, reaction, 1);
      for (unsigned int i=0; i<nd; i++)
      {
        if (i%d==0)
        {
          nub = cblas_dnrm2(2, velocity+i+1, 1);
          u_plus_phiu[i] = velocity[i]+nub;
        }
        else
        {
          u_plus_phiu[i] = velocity[i];
        }
      }


      if (NV_isnan(globalVelocity, m) | NV_isnan(velocity, nd) | NV_isnan(reaction, nd) | NV_isnan(s, n))
      {
        hasNotConverged = 2;
        break;
      }

      primalResidual(velocity, H, globalVelocity, w, primalConstraint, &pinfeas, tol);
      dualResidual(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol);
      complem = complemResidualNorm(u_plus_phiu, reaction, nd, n);
      udotr = cblas_ddot(nd, u_plus_phiu, 1, reaction, 1)/n;
      if (iteration == 0)
        barr_param = udotr * 0.3;
      // printf("n = %d\n", n);
      // barr_param *= 0.3;

      JA_prod(u_plus_phiu, reaction, nd, n, complemConstraint);
      max_uor_2mu = 0.0;
      for (int k = 0; k < nd; k+=d)
      {
        // complemConstraint[k] -= 2.*barr_param;
        tmp_uor_2mu = cblas_dnrm2(3, complemConstraint+k, 1);

        if (tmp_uor_2mu > max_uor_2mu) max_uor_2mu = tmp_uor_2mu;
      }
      uor_mu = cblas_dnrm2(nd, complemConstraint, 1);


      diff_fixp = complem;

      if (u_plus_phiu) { free(u_plus_phiu); u_plus_phiu = NULL; }
      if(R_plus_ru) { NM_free(R_plus_ru); R_plus_ru = NULL; }
      if(arrow_u_plus_phiu) { NM_free(arrow_u_plus_phiu); arrow_u_plus_phiu = NULL; }
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



    // /* computing the affine step-length */
    // alpha_primal = getStepLength(velocity, d_velocity, nd, n, 0.95);
    // alpha_dual = getStepLength(reaction, d_reaction, nd, n, 0.95);
    // alpha_primal = alpha_dual = fmin(alpha_primal, alpha_dual);



    // /* ----- Update variables ----- */
    // // if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
    // //   printIteresProbMatlabFile(iteration, globalVelocity, velocity, reaction, s, d_globalVelocity, d_velocity, d_reaction, d_s, d, n, m, iterates);
    //   // printIteresProbMatlabFile(iteration, pinfeas, dinfeas, udotr, diff_fixp, d, n, m, iterates);



    // cblas_daxpy(m, alpha_primal, d_globalVelocity, 1, globalVelocity, 1);
    // cblas_daxpy(nd, alpha_primal, d_velocity, 1, velocity, 1);
    // cblas_daxpy(nd, alpha_primal, d_reaction, 1, reaction, 1);
    // cblas_daxpy(n, alpha_primal, d_s, 1, s, 1);

    if (NV_isnan(globalVelocity, m) | NV_isnan(velocity, nd) | NV_isnan(reaction, nd) | NV_isnan(s, n))
    {
      hasNotConverged = 2;
      break;
    }

    if (hasNotConverged == 0 || hasNotConverged == 2)
      break;

    // /* Computation of the values of
    //  - primal residual: u - H*v - w - phi(s)
    //  - dual residual: M*v - f - H'*r
    //  - duality gap: u'*r
    //  - true duality gap: (value_of_primal - value_of_dual)
    //  - complementarity: u o r
    //  - projection error: r - proj(r-u)
    // */
    // primalResidual_s(velocity, H, globalVelocity, w, s, primalConstraint, &pinfeas, tol);
    // dualResidual(M, globalVelocity, H, reaction, f, dualConstraint, &dinfeas, tol);
    // complem = complemResidualNorm(velocity, reaction, nd, n);
    // udotr = cblas_ddot(nd, velocity, 1, reaction, 1)/n;

    // JA_prod(velocity, reaction, nd, n, complemConstraint);
    // max_uor_2mu = 0.0;
    // for (int k = 0; k < nd; k+=d)
    // {
    //   // if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] != SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST2
    //   //  && options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] != SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST3)
    //     complemConstraint[k] -= 2*barr_param;
    //   tmp_uor_2mu = cblas_dnrm2(3, complemConstraint+k, 1);

    //   if (tmp_uor_2mu > max_uor_2mu) max_uor_2mu = tmp_uor_2mu;
    // }

    // uor_mu = cblas_dnrm2(nd, complemConstraint, 1);


    // diff_fixp = 0.;
    // for (unsigned int i = 0; i<nd; i+=d)
    // {

    //   if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST1)
    //     diff_fixp += 0.25 * pow(s[i/d]*s[i/d] - velocity[i+1]*velocity[i+1] - velocity[i+2]*velocity[i+2], 2.);

    //   else
    //   {
    //     nub = cblas_dnrm2(2, velocity+i+1, 1);
    //     diff_fixp += (s[i/d] - nub)*(s[i/d] - nub);
    //   }
    // }
    // diff_fixp = sqrt(diff_fixp);


    // if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST4)
    // {
    //   norm_theta = sqrt(pow(pinfeas,2)+pow(dinfeas,2)+pow(complem,2)+pow(diff_fixp,2));
    // }

    // else if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST5)
    // {
    //   norm_theta = sqrt(pow(pinfeas,2)+pow(dinfeas,2)+pow(complem,2)+pow(diff_fixp,2));
    //   arr_norm_theta[iteration] = norm_theta;
    // }


    // // compute Projection Error
    // NM_gemv(1.0, P_mu_inv, velocity, 0.0, data->tmp_point->t_velocity);
    // NM_gemv(1.0, P_mu, reaction, 0.0, data->tmp_point->t_reaction);
    // (*computeError)(problem,
    //                 data->tmp_point->t_reaction, data->tmp_point->t_velocity, globalVelocity,
    //                 tol, options, norm_q, norm_b, &projerr);




    // // totalresidual = fmax(fmax(fmax(pinfeas, dinfeas),udotr),diff_fixp);
    // totalresidual = fmax(fmax(fmax(fmax(pinfeas, dinfeas),diff_fixp),2.*max_uor_2mu),4.*barr_param);
    // totalresidual_mu = fmax(fmax(fmax(pinfeas, dinfeas),uor_mu),diff_fixp);

    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE])
    {
      if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST3)
      {
        fprintf(iterates,"%d %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e;\n",
              iteration, pinfeas, dinfeas, diff_fixp, uor_mu, udotr, projerr, barr_param, alpha_dual, alpha_primal, alpha_diffixP,
              fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
              fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
              fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]),
              fabs(d_s[cblas_idamax(n, d_s, 1)]),
              LS_norm_p, LS_norm_d, LS_norm_c, LS_norm_f);
      }
      else if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST4 ||
               options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST5)
      {
        fprintf(iterates,"%d %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e;\n",
              iteration, pinfeas, dinfeas, diff_fixp, uor_mu, norm_theta, udotr, projerr, barr_param, alpha_dual, alpha_primal,
              fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
              fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
              fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]),
              fabs(d_s[cblas_idamax(n, d_s, 1)]),
              LS_norm_p, LS_norm_d, LS_norm_c, LS_norm_f);
      }
      else if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_TEST6)
      {
        fprintf(iterates,"%d %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e;\n",
              iteration, pinfeas, dinfeas, complem, udotr, projerr, barr_param, alpha_dual, min_alpha_primal,
              fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
              fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
              fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]),
              LS_norm_p, LS_norm_d, LS_norm_c);
      }
      else
        fprintf(iterates,"%d %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e;\n",
              iteration, pinfeas, dinfeas, diff_fixp, uor_mu, udotr, projerr, barr_param, alpha_primal,
              fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
              fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
              fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]),
              fabs(d_s[cblas_idamax(n, d_s, 1)]),
              LS_norm_p, LS_norm_d, LS_norm_c, LS_norm_f);


      if (iterates_2)
      {
        // printIteresProbMatlabFile(iteration, globalVelocity, velocity, reaction, s, d_globalVelocity, d_velocity, d_reaction, d_s, d, n, m, iterates_2);

        // store s & u0 & ub in matlab file for inspection
        fprintf(iterates_2,"s(%3i,:) = [",iteration+1);
        for(int i = 0; i < n; i++)
        {
          fprintf(iterates_2, "%8.20e, ", s[i]);
        }
        fprintf(iterates_2,"];\n");

        fprintf(iterates_2,"u0(%3i,:) = [",iteration+1);
        for(int i = 0; i < nd; i+=d)
        {
          fprintf(iterates_2, "%8.20e, ", velocity[i]);
        }
        fprintf(iterates_2,"];\n");

        fprintf(iterates_2,"ub(%3i,:) = [",iteration+1);
        for(int i = 0; i < n; i++)
        {
          fprintf(iterates_2, "%8.20e, ", cblas_dnrm2(2, velocity+i*d+1, 1));
        }
        fprintf(iterates_2,"];\n");
      }
    }

    if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST3)
    {
      numerics_printf_verbose(-1, "| %3i%c| %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |",
                            iteration, fws, pinfeas, dinfeas, diff_fixp, uor_mu, udotr, projerr, barr_param, alpha_dual, alpha_primal, alpha_diffixP,
                            fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
                            fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
                            fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]),
                            fabs(d_s[cblas_idamax(n, d_s, 1)]),
             LS_norm_p, LS_norm_d, LS_norm_c, LS_norm_f);
    }
    else if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST4 ||
             options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST5)
    {
      numerics_printf_verbose(-1, "| %3i%c| %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |",
                            iteration, fws, pinfeas, dinfeas, diff_fixp, uor_mu, norm_theta, udotr, projerr, barr_param, alpha_dual, alpha_primal,
                            fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
                            fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
                            fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]),
                            fabs(d_s[cblas_idamax(n, d_s, 1)]),
             LS_norm_p, LS_norm_d, LS_norm_c, LS_norm_f);

      totalresidual = norm_theta;
    }
    else if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_TEST6)
    {
      numerics_printf_verbose(-1, "| %3i%c| %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |",
                            iteration, fws, pinfeas, dinfeas, complem, udotr, projerr, barr_param, alpha_dual, min_alpha_primal,
                            fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
                            fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
                            fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]),
                            LS_norm_p, LS_norm_d, LS_norm_c);

    }
    else
      numerics_printf_verbose(-1, "| %3i%c| %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e || %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e | %.1e |",
                            iteration, fws, pinfeas, dinfeas, diff_fixp, complem, udotr, projerr, barr_param, sigma, alpha_primal,
                            fabs(d_globalVelocity[cblas_idamax(m, d_globalVelocity, 1)]),
                            fabs(d_velocity[cblas_idamax(nd, d_velocity, 1)]),
                            fabs(d_reaction[cblas_idamax(nd, d_reaction, 1)]),
                            fabs(d_s[cblas_idamax(n, d_s, 1)]),
             LS_norm_p, LS_norm_d, LS_norm_c, LS_norm_f);


    // if ( totalresidual <= tol )
    // {
    //   if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_TEST6)
    //   {
    //     for (size_t i = 0; i<nd; i+=d)
    //     {
    //       nub = cblas_dnrm2(2, velocity+i+1, 1);
    //       velocity[i] += nub;
    //     }
    //   }

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
    //   break;
    // }


    // kappa_mu = 0.7;
    // // if (totalresidual_mu <= 1e-8)
    // if (totalresidual_mu <= 100.*barr_param)
    // // if (//totalresidual_mu <= 10*barr_param &&
    // //     options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] != SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST2
    // //  && options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] != SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST3
    // //  && options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] != SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST4)
    // // if (iteration >= 22 && totalresidual <= tol)
    // {


      if (alpha_primal < 1.)
        sigma = 0.49;
      else
        sigma = (0.1+sigma)/2.;
      // barr_param = (udotr / n)*kappa_mu;
    // //   // printf("abs(ub'd_ub - |ub|*|d_ub|) = %e\n", check_sub);
    //   printf("\n");
    //   if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST3)
    //   {
    //     numerics_printf_verbose(-1, "| it  | pinfeas | dinfeas |  |s-ub| |  |uor|  |  u'r/n  | prj err | barpram | alpha 1 | alpha 2 | alpha 3 |  |dv|   |  |du|   |  |dr|   |  |ds|   | ls prim | ls dual | ls comp | ls fixP |");
    //   }
    //   else if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST4)
    //   {
    //     numerics_printf_verbose(-1, "| it  | pinfeas | dinfeas |  |s-ub| |  |uor|  | |Theta| |  u'r/n  | prj err | barpram | alpha 1 | alpha 2 |  |dv|   |  |du|   |  |dr|   |  |ds|   | ls prim | ls dual | ls comp | ls fixP |");
    //   }
    //   else if (options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] == SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_TEST6)
    //   {
    //     numerics_printf_verbose(-1, "| it  | pinfeas | dinfeas |  |uor|  |  u'r/n  | prj err | barpram | alpha_r | alpha_u |  |dv|   |  |du|   |  |dr|   | ls prim | ls dual | ls comp |");
    //   }
    //   else
    //     numerics_printf_verbose(-1, "| it  | pinfeas | dinfeas |  |s-ub| | |uor-mu||2max|uor-mu||   4*mu  |  u'r/n  | prj err | barpram |  alpha  |  |dv|   |  |du|   |  |dr|   |  |ds|   | ls prim | ls dual | ls comp | ls fixP |");
    // }


    if(J) { J = NM_free(J); J = NULL;}
    if(J_dense) { J_dense = NM_free(J_dense); J_dense = NULL;}

    iteration++;

  } // while loop
}



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

  if (options->solverId == SICONOS_GLOBAL_FRICTION_3D_IPM_SNM_SEP)
  {
    options->solverData = (double *)malloc(sizeof(double));
    double *pinfeas_ptr = (double *)options->solverData;
    *pinfeas_ptr = pinfeas;
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
  if (complemConstraint_mu) free(complemConstraint_mu);

  *info = hasNotConverged;
}

void gfc3d_ipm_snm_set_default(SolverOptions* options)
{
  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_GET_PROBLEM_INFO] = SICONOS_FRICTION_3D_IPM_GET_PROBLEM_INFO_NO;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S] = SICONOS_FRICTION_3D_IPM_IPARAM_UPDATE_S_AT_EACH_ITE;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] = SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_NOSCAL;
  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] = SICONOS_FRICTION_3D_IPM_IPARAM_LS_4X4_TEST4;
  // options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_LS_FORM] = SICONOS_FRICTION_3D_IPM_IPARAM_LS_3X3_TEST6;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_MEHROTRA] = SICONOS_FRICTION_3D_IPM_IPARAM_MEHROTRA_NO;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_ITERATES_MATLAB_FILE] = 0;

  options->iparam[SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT] = SICONOS_FRICTION_3D_IPM_IPARAM_REFINEMENT_YES;

  options->iparam[SICONOS_IPARAM_MAX_ITER] = 200;
  options->dparam[SICONOS_DPARAM_TOL] = 1e-10;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_1] = 1e-10;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_2] = 3.;
  options->dparam[SICONOS_FRICTION_3D_IPM_SIGMA_PARAMETER_3] = 1.;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_1] = 0.9;
  options->dparam[SICONOS_FRICTION_3D_IPM_GAMMA_PARAMETER_2] = 0.09; //0.095

}
