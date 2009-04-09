/* Siconos-sample version 3.0.0, Copyright INRIA 2005-2008.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
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

  FrictionContact3D problems needs some specific parameters, given to the FrictionContact3D_driver() function thanks to a Solver_Options structure. \n
  They are:\n
     - the name of the solver (ex: NSGS), used to switch to the right solver function
     - iparam[0]: max. number of iterations allowed
     - iparam[1]:
     - dparam[0]: tolerance
     - isStorageSparse: 1 if a SparseBlockStructuredMatrix is used for M, else 0 (double* storage)


  \brief
*/

#include "SiconosNumerics.h"

int main(int argc, char* argv[])
{


  // Problem Definition
  int NC = 3;//Number of contacts
  double M11[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; // Warning Fortran Storage
  double M22[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; // Warning Fortran Storage
  double M33[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; // Warning Fortran Storage

  double q[9] = { -1, 1, 3, -1, 1, 3, -1, 1, 3};
  double mu[3] = {0.1, 0.1, 0.1};

  int i, j, k;


  FrictionContact_Problem NumericsProblem;
  NumericsProblem.numberOfContacts = NC;
  NumericsProblem.isComplete = 0;
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
  MBlockMatrix->size = 3;
  int blocksize[3] = {3, 6, 9} ;
  MBlockMatrix->blocksize = blocksize;
  MBlockMatrix->filled1 = 4;
  MBlockMatrix->filled2 = 3;
  size_t index1_data[4] = {0, 1, 2, 3} ;
  size_t index2_data[3] = {0, 1, 2} ;
  MBlockMatrix->index1_data =  index1_data;
  MBlockMatrix->index2_data =  index2_data;



  // Unknown Declaration

  double *reaction = (double*)malloc(3 * NC * sizeof(double));
  double *velocity = (double*)malloc(3 * NC * sizeof(double));

  // Numerics and Solver Options

  Numerics_Options numerics_options;
  numerics_options.verboseMode = 1; // turn verbose mode to off by default


  Solver_Options numerics_solver_options;
  numerics_solver_options.filterOn = 0;
  numerics_solver_options.isSet = 1;

  strcpy(numerics_solver_options.solverName, "NSGS");

  numerics_solver_options.iSize = 5;
  numerics_solver_options.iparam = (int*)malloc(numerics_solver_options.iSize * sizeof(int));
  numerics_solver_options.dSize = 5;
  numerics_solver_options.dparam = (double*)malloc(numerics_solver_options.dSize * sizeof(double));

  int nmax = 10000; // Max number of iteration
  int localsolver = 0; // 0: projection on Cone, 1: Newton/AlartCurnier,  2: projection on Cone with local iteration, 2: projection on Disk  with diagonalization,
  double tolerance = 1e-10;
  double localtolerance = 1e-12;


  numerics_solver_options.iparam[0] = nmax ;
  numerics_solver_options.iparam[4] = localsolver ;
  numerics_solver_options.dparam[0] = tolerance ;
  numerics_solver_options.dparam[2] = localtolerance ;

  //Driver call
  frictionContact3D_driver(&NumericsProblem,
                           reaction , velocity,
                           &numerics_solver_options, &numerics_options);



  // Solver output
  printf("\n");
  for (k = 0 ; k < 3 * NC; k++) printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e \n ", k, velocity[k], k , reaction[k]);
  printf("\n");


  free(reaction);
  free(velocity);
  free(MM);
  free(MBlockMatrix);

}
