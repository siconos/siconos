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

/* Returns the 2-norm of a vector - uses long double - based on blas_dnrm2 */
float_type dnrm2l(const unsigned int n, const double * x);


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
double getNewtonStepLength(const double * const x, const double * const dx,
                                  const unsigned int vecSize, const unsigned int varsCount, const double gamma);

/* Returns array of step-lengths to the boundary reduced by a factor gamma. Uses long double. */
double *array_getStepLength(const double * const x, const double * const dx, const unsigned int vecSize,
                     const unsigned int varsCount, const double gamma);


/* Returns the maximum step-length to the boundary reduced by a factor gamma. Uses long double. */
double getStepLength(const double * const x, const double * const dx, const unsigned int vecSize,
                            const unsigned int varsCount, const double gamma);


/**
 * Returns the primal constraint vector for global fricprob ( velocity - H @ globalVelocity - w )
 * \param velocity is the vector of relative velocities.
 * \param H is the constraint matrix.
 * \param globalVelocity is the vector of generalized velocities.
 * \param w is the constraint vector.
 * \param out is the result velocity - H x globalVelocity - w vector.
 * \param rnorm is the relative norm of out = |out|/max{|velocity|, |H x globalVelocity|, |w|}
 */
// void primalResidual(const double * velocity, NumericsMatrix * H, const double * globalVelocity, const double * w,
// 		    const double * s, double * out, double * rnorm, const double tol);
void primalResidual(const double * velocity, NumericsMatrix * H, const double * globalVelocity, const double * w,
                    double * out, double * rnorm, const double tol);


/**
 * Returns the dual constraint vector for global fricprob ( M @ globalVelocity + f - H @ reaction )
 * \param M is the mass matrix.
 * \param globalVelocity is the vector of generalized velocities.
 * \param H is the constraint matrix.
 * \param reaction is the vector of reaction forces at each contact point.
 * \param f is the constraint vector (vector of internal and external forces).
 * \param out os the result M x globalVelocity + f - H' x reaction vector.
 * \param rnorm is the relative 2-norm of out = |out| / max{|M x globalVelocity|, |f|, |H' x r|}
 */
void dualResidual(NumericsMatrix * M, const double * globalVelocity, NumericsMatrix * H, const double * reaction, const double * f,
		  double * out, double * rnorm, const double tol);


/**
 * Returns the Inf-norm of primal residual vecor ( velocity - H @ globalVelocity - w )
 * \param velocity is the vector of relative velocities.
 * \param H is the constraint matrix.
 * \param globalVelocity is the vector of generalized velocities.
 * \param w is the constraint vector.
 */
/* double primalResidualNorm(const double * velocity, NumericsMatrix * H, */
/*                                  const double * globalVelocity, const double * w); */


/**
 * Returns the Inf-norm of the dual residual vector ( M @ globalVelocity + f - H @ reaction )
 * \param M is the mass matrix.
 * \param globalVelocity is the vector of generalized velocities.
 * \param H is the constraint matrix.
 * \param reaction is the vector of reaction forces at each contact point.
 * \param f is the constraint vector (vector of internal and external forces).
 */
/* double dualResidualNorm(NumericsMatrix * M, const double * globalVelocity, */
/*                                NumericsMatrix * H, const double * reaction, const double * f); */


/**
 * Returns the Inf-norm of the cemplementarity residual vector ( velocity o reaction )
 * \param velocity is the vector of relative velocities at each contact point.
 * \param reaction is the vector of reaction forces at each contact point.
 * \param vecSize is the size of reaction and velocity vector.
 * \param varsCount is a count of variables concatenated into vectors reaction and velocity.
 */
double complemResidualNorm(const double * const velocity, const double * const reaction,
                                  const unsigned int vecSize, const unsigned int varsCount);


/* Returns the 2-norm of the complementarity residual vector = 2-norm of the Jordan product (Qp*velocity) o (Qp_inv * reaction)  */
double complemResidualNorm_p(const double * const velocity, const double * const reaction,
                                    const unsigned int vecSize, const unsigned int varsCount);


/* Returns the 2-norm of the complementarity residual vector = 2-norm of the Jordan product (Qp*velocity) o (Qp_inv * reaction)  */
/* This computation is done with the formula "F" */
double complemResidualNorm_p_F(NumericsMatrix * Qp, NumericsMatrix * Qpinv,
                                  const double * const velocity, const double * const reaction,
                                  const unsigned int vecSize, const unsigned int varsCount);


/* computation of the duality gap  */
double dualGap(NumericsMatrix * M, const double * f, const double * w, const double * globalVelocity, const double * reaction, const unsigned int nd, const unsigned int m);


/* computation of the relative gap  */
double relGap(NumericsMatrix * M, const double * f, const double * w, const double * globalVelocity, const double * reaction, const unsigned int nd, const unsigned int m, const double gapVal);


/* Establish an array of calculation errors  */
void setErrorArray(double * error, const double pinfeas, const double dinfeas,
		   const double dualgap, const double complem, const double complem_p, const double projerr);


/* Return the 2-norm of the difference between two vectors */
double norm2VecDiff (const double * vec1, const double * vec2, const unsigned int vecSize);

