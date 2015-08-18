/* Siconos-sample , Copyright INRIA 2005-2011.
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
  numericsProblem.numberOfContacts = NC;
  numericsProblem.dimension = 3;
  numericsProblem.mu = mu;
  numericsProblem.q = q;
  numericsProblem.b = b;


  numericsProblem.M = newNumericsMatrix();
  NumericsMatrix *MM = numericsProblem.M ;
  MM->storageType = 0;
  MM->matrix0 = M;
  MM->size0 = n;
  MM->size1 = n;


  numericsProblem.H  = newNumericsMatrix();
  NumericsMatrix *HH = numericsProblem.H;
  HH->storageType = 0;
  HH->matrix0 = H;
  HH->size0 = n;
  HH->size1 = m;


  /*     FILE * foutput = fopen("Example_GlobalFrictionContact.dat", "w"); */
  /*     globalFrictionContact_printInFile(&numericsProblem,  foutput ); */
  /*     fclose(foutput); */



  // Unknown Declaration

  double *reaction = (double*)malloc(m * sizeof(double));
  double *velocity = (double*)malloc(m * sizeof(double));
  double *globalVelocity = (double*)malloc(n * sizeof(double));
  for (k = 0 ; k < m; k++)
  {
    velocity[k] = 0.0;
    reaction[k] = 0.0;
  };
  for (k = 0 ; k < n; k++)
  {
    globalVelocity[k] = 0.0;
  };

  // Numerics and Solver Options

  NumericsOptions numerics_options;
  numerics_options.verboseMode = 1; // turn verbose mode to off by default


  SolverOptions * numerics_solver_options = (SolverOptions *)malloc(sizeof(SolverOptions)) ;
  //char solvername[10]= "NSGS";
  /*\warning Must be adpated  for future globalFrictionContact3D_setDefaultSolverOptions*/
  globalFrictionContact3D_setDefaultSolverOptions(numerics_solver_options, SICONOS_FRICTION_3D_GLOBAL_NSGS);
  numerics_solver_options->dparam[0] = 1e-14;
  numerics_solver_options->iparam[0] = 100000;
  //Driver call
  i = 0;
  info = globalFrictionContact3D_driver(&numericsProblem,
                                        reaction , velocity, globalVelocity,
                                        numerics_solver_options, &numerics_options);

  deleteSolverOptions(numerics_solver_options);
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
  return info;


}
