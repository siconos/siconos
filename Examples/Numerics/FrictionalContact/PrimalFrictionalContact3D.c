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
  int info = -1;



  int NC = 1;//Number of contacts
  int Ndof = 9;//Number of DOF
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


  int i = 0, j = 0, k = 0;

  PrimalFrictionContact_Problem NumericsProblem;
  NumericsProblem.numberOfContacts = NC;
  NumericsProblem.isComplete = 0;
  NumericsProblem.mu = mu;
  NumericsProblem.q = q;
  NumericsProblem.b = b;


  NumericsProblem.M = (NumericsMatrix*)malloc(sizeof(NumericsMatrix));
  NumericsMatrix *MM = NumericsProblem.M ;
  MM->storageType = 0;
  MM->matrix0 = M;
  MM->size0 = n;
  MM->size1 = n;


  NumericsProblem.H  = (NumericsMatrix*)malloc(sizeof(NumericsMatrix));
  NumericsMatrix *HH = NumericsProblem.H;
  HH->storageType = 0;
  HH->matrix0 = H;
  HH->size0 = m;
  HH->size1 = m;


  // Unknown Declaration

  double *reaction = (double*)malloc(m * sizeof(double));
  double *velocity = (double*)malloc(m * sizeof(double));
  double *globalVelocity = (double*)malloc(n * sizeof(double));

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

  for (i = 0; i < numerics_solver_options.iSize ; i++) numerics_solver_options.iparam[i] = 0;
  for (i = 0; i < numerics_solver_options.dSize ; i++) numerics_solver_options.dparam[i] = 0.0;


  numerics_solver_options.iparam[0] = nmax ;
  numerics_solver_options.iparam[4] = localsolver ;
  numerics_solver_options.dparam[0] = tolerance ;
  numerics_solver_options.dparam[2] = localtolerance ;

  //Driver call
  //i=0;
  //while (0==0){
  info = primalFrictionContact3D_driver(&NumericsProblem,
                                        reaction , velocity, globalVelocity,
                                        &numerics_solver_options, &numerics_options);
  //i++;
  //printf("i=%i\n", i);
  //}


  // Solver output
  printf("\n");
  for (k = 0 ; k < m; k++) printf("velocity[%i] = %12.8e \t \t reaction[%i] = %12.8e \n ", k, velocity[k], k , reaction[k]);
  for (k = 0 ; k < n; k++) printf("globalVelocity[%i] = %12.8e \t \n ", k, globalVelocity[k]);
  printf("\n");


  free(reaction);
  free(velocity);
  free(globalVelocity);
  free(NumericsProblem.M);
  free(NumericsProblem.H);
  free(numerics_solver_options.iparam);
  free(numerics_solver_options.dparam);


  return info;


}
