/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

/*!\file GMPReduced.h
 * \brief Functions to solve GMP */

#ifndef NUMERICS_GMPREDUCED
#define NUMERICS_GMPREDUCED
#include "NumericsFwd.h"  // for GenericMechanicalProblem, SolverOptions

/* The equalities are eliminated.
 *
 * 0=(Me_1 Me_2)(Re Ri)' + Qe
 * Vi=(Mi_1 Mi_2)(Re Ri)' + Qi
 *
 * Re=-Me_1^{-1}(Me_2Ri+Qe)
 *
 * Vi=(Mi_2-Mi_1 Me_1^{-1} Me_2)Ri+Qi-Mi1 Me_1^{-1} Qe
 *
 */
void gmp_reduced_solve(GenericMechanicalProblem* problem, double* reaction, double* velocity,
                       int* info, SolverOptions* options);

/*  The equalities are assembled in an single block.
 *
 * 0=(Me_1 Me_2)(Re Ri)' + Qe
 * Vi=(Mi_1 Mi_2)(Re Ri)' + Qi
 *
 * and GS.
 */
void gmp_reduced_equality_solve(GenericMechanicalProblem* problem, double* reaction,
                                double* velocity, int* info, SolverOptions* options);

/* It converts the solution of the reduced problem to the initial problem.
 */
void gmp_reduced_convert_solution(GenericMechanicalProblem* problem, double* reaction,
                                  double* velocity, double* Re, double* Rreduced,
                                  double* Vreduced);

/* If the GMP is composed only of equalities and complementarities,
 * it is possible to used MLCP solvers.
 */
void gmp_as_mlcp(GenericMechanicalProblem* problem, double* reaction, double* velocity,
                 int* info, SolverOptions* options);

#endif
