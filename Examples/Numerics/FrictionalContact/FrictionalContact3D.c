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
  A very simple example to chow how to use SiconosNumerics to solve the time--discretized FrictionalContact3D problem
  \brief
*/

#include "SiconosNumerics.h"

int main(int argc, char* argv[])
{


  // Problem Definition
  int NC = 3;//Number of contacts
  double W[9][9];// = { };
  double q[9];  // = { };
  double mu[3];   // = { };

  int i, j, k;
  for (j = 0 ; j < 9; j++)
  {
    for (k = 0 ; k < 9; k++)
      W[j][k] = 0.0;
  }

  for (i = 0 ; i < NC; i++)
  {
    mu[i] = 0.1;
    for (j = 0 ; j < 3; j++)
    {
      q[3 * i + j] = -1.0 + j * 2.0;
      W[3 * i + j][3 * i + j] = 1.0;
    }
  }


  FrictionContact_Problem NumericsProblem;
  NumericsProblem.numberOfContacts = NC;
  NumericsProblem.isComplete = 0;
  NumericsProblem.mu = mu;
  NumericsProblem.q = q;

  NumericsMatrix *MM = (NumericsMatrix*)malloc(sizeof(*MM));
  MM->storageType = 0;
  MM->matrix0 = W;
  MM->size0 = 3 * NC;
  MM->size1 = 3 * NC;

  NumericsProblem.M = MM;

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
  int localsolver = 1; // 0: projection, 1: Newton/AlartCurnier, 2: Newton/Fischer-Burmeister, 3: Path/Glocker

  double tolerance = 1e-6;
  double localtolerance = 1e-8;


  numerics_solver_options.iparam[0] = nmax ;
  numerics_solver_options.iparam[4] = localsolver ;
  numerics_solver_options.dparam[0] = tolerance ;
  numerics_solver_options.dparam[2] = localtolerance ;

  //Driver call
  frictionContact3D_driver(&NumericsProblem,
                           reaction , velocity,
                           &numerics_solver_options, &numerics_options);



  // Solver output
  for (k = 0 ; k < 3 * NC; k++) printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e \n ", k, velocity[k], k , reaction[k]);


  free(reaction);
  free(velocity);
  free(MM);


}
