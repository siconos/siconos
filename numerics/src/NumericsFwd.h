/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

/*!\file NumericsFwd.h
 * \brief Forward declaration of numerics structures*/

#ifndef NumericsFwd_h
#define NumericsFwd_h

#define TYPEDEF_STRUCT(X) \
  typedef struct X X; \

// Matrices storage
TYPEDEF_STRUCT(NumericsMatrix)
TYPEDEF_STRUCT(NumericsSparseMatrix)
TYPEDEF_STRUCT(NumericsSparseLinearSolverParams)
TYPEDEF_STRUCT(SparseBlockStructuredMatrix)
TYPEDEF_STRUCT(SparseBlockStructuredMatrixPred)
TYPEDEF_STRUCT(SparseBlockCoordinateMatrix)

// Nonsmooth solvers
TYPEDEF_STRUCT(SolverOptions)

// Nonsmooth problems 
TYPEDEF_STRUCT(SecondOrderConeLinearComplementarityProblem)
TYPEDEF_STRUCT(SecondOrderConeLinearComplementarityProblem_as_VI)
TYPEDEF_STRUCT(RelayProblem)
TYPEDEF_STRUCT(NonlinearComplementarityProblem)
TYPEDEF_STRUCT(MixedLinearComplementarityProblem)
TYPEDEF_STRUCT(MixedComplementarityProblem)
TYPEDEF_STRUCT(MixedComplementarityProblem2)
TYPEDEF_STRUCT(LinearComplementarityProblem)
TYPEDEF_STRUCT(GlobalFrictionContactProblem)
TYPEDEF_STRUCT(GenericMechanicalProblem)
TYPEDEF_STRUCT(listNumericsProblem)
TYPEDEF_STRUCT(FrictionContactProblem_as_VI)
TYPEDEF_STRUCT(FrictionContactProblem)
TYPEDEF_STRUCT(VariationalInequality)
TYPEDEF_STRUCT(AffineVariationalInequalities)



#endif
