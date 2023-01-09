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

typedef long double float_type;
/* typedef double float_type; */


/* Returns the square of 2-norm of a vector - uses long double - based on blas_dnrm2 */
static float_type dnrm2sqrl(const unsigned int n, const double * x)
{
  float_type norm, scale, ssq, absxi, quo;

  if (n < 1)
    norm = 0.0;
  else if (n == 1)
    norm = fabsl(x[0]);
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
      norm = scale * scale * ssq;
    }
  }
  return norm;
}



/* Returns the maximum step-length to the boundary reduced by a factor gamma. Uses long double. */
static double getStepLength(const double * const x, const double * const dx, const unsigned int vecSize,
                     const unsigned int varsCount, const double gamma)
{
  unsigned int dimension = (int)(vecSize / varsCount);
  unsigned int pos;
  float_type aL, bL, cL, dL, alphaL, nxb;
  double min_alpha;

  min_alpha = 1e20; //1.0;

  for(unsigned int i = 0; i < varsCount; ++i)
  {
    pos = i * dimension;
    aL = dnrm2l(dimension-1, dx+pos+1);
    aL = (dx[pos] - aL)*(dx[pos] + aL);
    bL = x[pos]*dx[pos];
    for (int k = 1; k < dimension; bL -= x[pos+k]*dx[pos+k], k++);
    nxb = dnrm2l(dimension-1, x+pos+1);
    // cL = (x[pos] - nxb)*(x[pos] + nxb);
    cL = x[pos] - nxb;
    if (cL <= 0.) cL = DBL_EPSILON*(x[pos] + nxb); else cL = (x[pos] - nxb)*(x[pos] + nxb); // to avoid negative number b/c of different data types
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



/* Rel gap = gapVal / (1 + abs(primal value) + abs(dual value)) */
static double relGap(NumericsMatrix * M, const double * f, const double * w, const double * globalVelocity, const double * reaction, const unsigned int nd, const unsigned int m, const double gapVal)
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


NumericsMatrix * compute_JQinv2Jt(const double *u1, const double *r1, const double *u2, const double *r2, const size_t vecSize, const size_t varsCount);

CS_INT cs_dupl_zeros (cs *A);
