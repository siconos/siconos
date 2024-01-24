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


typedef long double float_type;
/* typedef double float_type; */

#define MIN_RELATIVE_SCALING sqrt(DBL_EPSILON)
#define STANDARD 0
#define NORM_2 1
#define NORM_INF 2
#define NORM_2_INF 3


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
 * \param type is the norm type used: NORM_2 = L-2, NORM_INF = L-inf
 */
void primalResidual_s(const double * velocity, NumericsMatrix * H, const double * globalVelocity, const double * w,
		    const double * s, double * out, double * rnorm, const double tol);
void primalResidual(const double * velocity, NumericsMatrix * H, const double * globalVelocity, const double * w,
                    double * out, double * rnorm, const double tol);
void primalResidual_type(const double * velocity, NumericsMatrix * H, const double * globalVelocity, const double * w,
                    double * out, double * rnorm, const double tol, const int type);


/**
 * Returns the dual constraint vector for global fricprob ( M @ globalVelocity + f - H @ reaction )
 * \param M is the mass matrix.
 * \param globalVelocity is the vector of generalized velocities.
 * \param H is the constraint matrix.
 * \param reaction is the vector of reaction forces at each contact point.
 * \param f is the constraint vector (vector of internal and external forces).
 * \param out os the result M x globalVelocity + f - H' x reaction vector.
 * \param rnorm is the relative 2-norm of out = |out| / max{|M x globalVelocity|, |f|, |H' x r|}
 * \param type is the norm type used: NORM_2 = L-2, NORM_INF = L-inf
 */
void dualResidual(NumericsMatrix * M, const double * globalVelocity, NumericsMatrix * H, const double * reaction, const double * f,
		  double * out, double * rnorm, const double tol);
void dualResidual_type(NumericsMatrix * M, const double * globalVelocity, NumericsMatrix * H, const double * reaction, const double * f,
      double * out, double * rnorm, const double tol, const int type);

/**
 * Returns the scalar product of 2 vectors depending on:
 * type = STANDARD  : <x,y>_2     = sum_i ( xi * yi )
 * type = NORM_INF  : <x,y>_inf   = max_i ( xi * yi )
 * type = NORM_2_INF: <x,y>_2_inf = max_i norm_2 ( xi o yi )
 *
 */
double xdoty_type(const unsigned int varsCount, const unsigned int vecSize, const double * x, const double * y, const int type);

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

int gfc3d_compute_error_r(GlobalFrictionContactProblem* problem,
                        double*  reaction, double*  velocity,
                        double*  globalVelocity,
                        double tolerance,
                        SolverOptions * options,
                        double norm_q, double norm_b,
                        double* restrict error);

void NM_clear_cone_matrix_H(NumericsMatrix *H, unsigned int n_cones_to_clear, int *cones_to_clear);

NumericsMatrix * NM_extract(NumericsMatrix *A, int n_rows, int *target_rows, int n_cols, int *target_cols);

double projectionError(const double * velocity, const double * reaction, const unsigned int nc, const double tol);

/* Routine is to read matrix-block in hdf5 file
 * type = 0 : return number of blocks (blk_index is neglected)
 * = 1: contact indices, = 2: body indices (in this case, the 1st element is the size of vector)
 * = 3: rank of block-matrix Hc (blk_index is neglected)
 */
int *read_fricprob_block(const char* path, int type, int blk_index);


