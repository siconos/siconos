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

#include "GlobalFrictionContactProblem.h"
#include "FrictionContactProblem.h"
#include "NumericsMatrix.h"
#include "SiconosNumerics.h"
#include "SparseBlockMatrix.h"
#include "NumericsSparseMatrix.h"
#include "SolverOptions.h"
#include "Friction_cst.h"
#include "gfc3d_Solvers.h"
#include "SiconosNumerics.h"

int main(int argc, char* argv[])
{


  // Problem Definition
  int info = -1;



  int NC = 1;//Number of contacts
  int m = 3;
  int n = 9;
  double M[81] = {1, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 1, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 1, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 1, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 1, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 1, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 1, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 1, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 1
                 };

  double H[27] = {1, 0, 0, 0, 0, 0, -1, 0, 0,
                  0, 1, 0, 0, 0, 0, 0, -1, 0,
                  0, 0, 1, 0, 0, 0, 0, 0, -1
                 };

  double q[9] = { -3, -3, -3, -1, 1, 3, -1, 1, 3};
  double b[3] = {0, 0, 0};
  double mu[1] = {0.1};

  /*     int NC = 3;//Number of contacts  */
  /*     int Ndof = 9;//Number of DOF  */
  /*     double M[81] = {1, 0, 0, 0, 0, 0, 0, 0, 0,  */
  /*        0, 1, 0, 0, 0, 0, 0, 0, 0,  */
  /*        0, 0, 1, 0, 0, 0, 0, 0, 0,  */
  /*        0, 0, 0, 1, 0, 0, 0, 0, 0,  */
  /*        0, 0, 0, 0, 1, 0, 0, 0, 0,  */
  /*        0, 0, 0, 0, 0, 1, 0, 0, 0,  */
  /*        0, 0, 0, 0, 0, 0, 1, 0, 0,  */
  /*        0, 0, 0, 0, 0, 0, 0, 1, 0,  */
  /*        0, 0, 0, 0, 0, 0, 0, 0, 1}; */
  /*     double H[81] = {1, 0, 0, 0, 0, 0, 0, 0, 0,  */
  /*        0, 1, 0, 0, 0, 0, 0, 0, 0,  */
  /*        0, 0, 1, 0, 0, 0, 0, 0, 0,  */
  /*        0, 0, 0, 1, 0, 0, 0, 0, 0,  */
  /*        0, 0, 0, 0, 1, 0, 0, 0, 0,  */
  /*        0, 0, 0, 0, 0, 1, 0, 0, 0,  */
  /*        0, 0, 0, 0, 0, 0, 1, 0, 0,  */
  /*        0, 0, 0, 0, 0, 0, 0, 1, 0,  */
  /*        0, 0, 0, 0, 0, 0, 0, 0, 1}; */


  /*     double q[9] = {-1, 1, 3, -1, 1, 3, -1, 1, 3}; */
  /*     double b[9] = {0, 0, 0,0, 0, 0,0, 0, 0 }; */
  /*     double mu[3] = {0.1,0.1,0.1};    */


  int i = 0, k = 0;

  GlobalFrictionContactProblem numericsProblem;
  globalFrictionContact_null(&numericsProblem);
  numericsProblem.numberOfContacts = NC;
  numericsProblem.dimension = 3;
  numericsProblem.mu = mu;
  numericsProblem.q = q;
  numericsProblem.b = b;


  numericsProblem.M = NM_new();
  NumericsMatrix *MM = numericsProblem.M ;
  MM->storageType = NM_DENSE;
  MM->matrix0 = M;
  MM->size0 = n;
  MM->size1 = n;


  numericsProblem.H  = NM_new();
  NumericsMatrix *HH = numericsProblem.H;
  HH->storageType = NM_DENSE;
  HH->matrix0 = H;
  HH->size0 = n;
  HH->size1 = m;


  /*     FILE * foutput = fopen("Example_GlobalFrictionContact.dat", "w"); */
  /*     globalFrictionContact_printInFile(&numericsProblem,  foutput ); */
  /*     fclose(foutput); */



  // Unknown Declaration

  double *reaction = (double*)calloc(m, sizeof(double));
  double *velocity = (double*)calloc(m, sizeof(double));
  double *globalVelocity = (double*)calloc(n, sizeof(double));

  // Solver Options

  SolverOptions * numerics_solver_options = (SolverOptions *)malloc(sizeof(SolverOptions)) ;
  //char solvername[10]= "NSGS";
  /*\warning Must be adpated  for future globalFrictionContact3D_setDefaultSolverOptions*/
  gfc3d_setDefaultSolverOptions(numerics_solver_options, SICONOS_GLOBAL_FRICTION_3D_NSGS);
  numerics_solver_options->dparam[0] = 1e-14;
  numerics_solver_options->iparam[0] = 100000;
  //Driver call
  i = 0;
  info = gfc3d_driver(&numericsProblem,
		      reaction , velocity, globalVelocity,
		      numerics_solver_options);

  solver_options_delete(numerics_solver_options);
  // Solver output
  printf("\n");
  for (k = 0 ; k < m; k++) printf("velocity[%i] = %12.8e \t \t reaction[%i] = %12.8e \n ", k, velocity[k], k , reaction[k]);
  for (k = 0 ; k < n; k++) printf("globalVelocity[%i] = %12.8e \t \n ", k, globalVelocity[k]);
  printf("\n");

  free(numerics_solver_options);
  free(reaction);
  free(velocity);
  free(globalVelocity);
  free(numericsProblem.M);
  free(numericsProblem.H);
  gfc3d_free_workspace(&numericsProblem);
  return info;


}
