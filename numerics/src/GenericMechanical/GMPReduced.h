/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

/*!\file GMPReduced.h
 * \brief Functions to solve GMP */

#ifndef NUMERICS_GMPREDUCED
#define NUMERICS_GMPREDUCED
#include <stdio.h>
/*
 * The equalities are eliminated.
 *
 *0=(Me_1 Me_2)(Re Ri)' + Qe
 *Vi=(Mi_1 Mi_2)(Re Ri)' + Qi
 *
 *Re=-Me_1^{-1}(Me_2Ri+Qe)
 *
 *Vi=(Mi_2-Mi_1 Me_1^{-1} Me_2)Ri+Qi-Mi1 Me_1^{-1} Qe
 *
 */
void GMPReducedSolve(GenericMechanicalProblem* pInProblem, double *reaction , double *velocity, int* info, SolverOptions* options);
/*
 * The equalities are assamblate in an unique block.
 *
 *0=(Me_1 Me_2)(Re Ri)' + Qe
 *Vi=(Mi_1 Mi_2)(Re Ri)' + Qi
 *
 *and GS.
 */
void GMPReducedEqualitySolve(GenericMechanicalProblem* pInProblem, double *reaction , double *velocity, int* info, SolverOptions* options);
/*
 * It converts the solution of the reduced problem to the initial problem.
 */
void GMPReducedSolToSol(GenericMechanicalProblem* pInProblem, double * reaction, double * velocity,
                        double * Re, double * Rreduced, double * Vreduced);

/*
 *If the GMP is composed only of equalities and complementarities, it is possible to used MLCP solvers.
 *
 */
void GMPasMLCP(GenericMechanicalProblem* pInProblem, double *reaction , double *velocity, int* info, SolverOptions* options);

/*
 * miscellaneous
 */
void printDenseMatrice(char* name, FILE * titi, double * m, int N, int M);


#endif
