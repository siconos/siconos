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

/*!\file FrictionalContact3D.cpp
  A very simple example to chow how to use SiconosNumerics to solve the time--discretized FrictionalContact3D problem with a dense storage.
  The problem: Find \f$(reaction,velocity)\f$ such that:\n

  \f$
  \left\lbrace
  \begin{array}{l}
  M \ reaction + q = velocity \\
  0 \le reaction_n \perp velocity_n \ge 0 \\
  -velocity_t \in \partial\psi_{D_(\mu reaction_n)}(reaction_t)\\
  D_(\mu reaction_n) = \{ reaction_t \mid  \|reaction_t\| \leq \mu reaction_n  \}
  \end{array}
  \right.
  \f$

  \f$ reaction, velocity, q\f$ are vectors of size n and \f$ M \f$ is a nXn matrix, with \f$ n = 2*nc or 3*nc \f$, nc being the number of contacts. \n
  \f$ reaction_n\f$ represents the normal part of the reaction while \f$ reaction_t\f$ is its tangential part.

  \f$ \mu \f$ is the friction coefficient (it may be different for each contact).

  Use the generic function frictionContact3D_driver() to call one the the specific solvers listed below:

  - frictionContact3D_nsgs() : non-smooth Gauss-Seidel solver

  (see the functions/solvers list in FrictionContact3D_Solvers.h)

  FrictionContact3D problems needs some specific parameters, given to the FrictionContact3D_driver() function thanks to a SolverOptions structure. \n
  They are:\n
     - the name of the solver (ex: NSGS), used to switch to the right solver function
     - iparam[0]: max. number of iterations allowed
     - iparam[1]:
     - dparam[0]: tolerance
     - isStorageSparse: 1 if a SparseBlockStructuredMatrix is used for M, else 0 (double* storage)


  \brief
*/

#include "SiconosNumerics.h"
#include "GlobalFrictionContactProblem.h"
#include "FrictionContactProblem.h"
#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"
#include "SparseBlockMatrix.h"
#include "SolverOptions.h"
#include "Friction_cst.h"
#include "fc3d_Solvers.h"

int main(int argc, char* argv[])
{


  // Problem Definition
  int NC = 3;//Number of contacts
  double M11[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; // Warning Fortran Storage
  double M22[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; // Warning Fortran Storage
  double M33[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; // Warning Fortran Storage

  double q[9] = { -1, 1, 3, -1, 1, 3, -1, 1, 3};
  double mu[3] = {0.1, 0.1, 0.1};

  int k;


  FrictionContactProblem NumericsProblem;
  NumericsProblem.numberOfContacts = NC;
  NumericsProblem.dimension = 3;
  NumericsProblem.mu = mu;
  NumericsProblem.q = q;

  NumericsMatrix *MM = (NumericsMatrix*)malloc(sizeof(*MM));
  MM->storageType = 1;
  MM->size0 = 3 * NC;
  MM->size1 = 3 * NC;
  NumericsProblem.M = MM;

  SparseBlockStructuredMatrix *MBlockMatrix = (SparseBlockStructuredMatrix*)malloc(sizeof(*MBlockMatrix));
  MM->matrix1 = MBlockMatrix;
  MBlockMatrix->nbblocks = 3;
  double * block[3] = {M11, M22, M33};
  MBlockMatrix->block = block;
  MBlockMatrix->blocknumber0 = 3;
  MBlockMatrix->blocknumber1 = 3;
  unsigned int blocksize[3] = {3, 6, 9} ;
  MBlockMatrix->blocksize0 = blocksize;
  MBlockMatrix->blocksize1 = blocksize;
  MBlockMatrix->filled1 = 4;
  MBlockMatrix->filled2 = 3;
  size_t index1_data[4] = {0, 1, 2, 3} ;
  size_t index2_data[3] = {0, 1, 2} ;
  MBlockMatrix->index1_data =  index1_data;
  MBlockMatrix->index2_data =  index2_data;



  // Unknown Declaration

  double *reaction = (double*)calloc(3 * NC, sizeof(double));
  double *velocity = (double*)calloc(3 * NC, sizeof(double));

  // Solver Options


  //Driver call
  FILE * ff =  fopen("FC3D_Example1_SBM.dat", "w");
  frictionContact_printInFile(&NumericsProblem, ff);
  fclose(ff);

  SolverOptions *numerics_solver_options  = (SolverOptions *)malloc(sizeof(SolverOptions));
  fc3d_setDefaultSolverOptions(numerics_solver_options, SICONOS_FRICTION_3D_NSGS);

  numerics_solver_options->dparam[0] = 1e-16;


  int info = fc3d_driver(&NumericsProblem,
			 reaction , velocity,
			 numerics_solver_options);


  solver_options_delete(numerics_solver_options);


  // Solver output
  printf("\n");
  for (k = 0 ; k < 3 * NC; k++) printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e \n ", k, velocity[k], k , reaction[k]);
  printf("\n");
  ff =  fopen("Example1_Fc3D_SBM-solution.dat", "w");
  for (k = 0 ; k < 3 * NC; k++) fprintf(ff, "%32.24e", velocity[k]);
  fprintf(ff, "\n");
  for (k = 0 ; k < 3 * NC; k++) fprintf(ff, "%32.24e", reaction[k]);
  fprintf(ff, "\n");
  fclose(ff);


  free(numerics_solver_options);
  free(reaction);
  free(velocity);
  free(MM);
  free(MBlockMatrix);
  return info;
}
