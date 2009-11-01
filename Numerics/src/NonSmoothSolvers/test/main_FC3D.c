/* Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
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

/*
  Tests functions for NumericsMatrix structure

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "FrictionContact_Problem.h"
#include "LA.h"
#include <math.h>


int test_Series(FrictionContact_Problem* problem)
{
  int info = -1;


  return info;
}

int main(void)
{

  printf("========= Starts Numerics tests for FrictionalContact_problem ========= \n");



  FrictionContact_Problem * testproblem = NULL;
  testproblem = (FrictionContact_Problem*)malloc(sizeof(FrictionContact_Problem));

  int ntestfile = 1;
  char ** testfile = (char **)malloc(ntestfile * sizeof(char *));
  for (int i = 0 ; i < ntestfile; i++)
  {
    testfile[i] = (char *)malloc(20 * sizeof(char));
  }

  testfile[0] = "./DATA/Example1_Fc3D_SBM.dat" ;
  int info = -1;
  FILE * file;

  for (int i = 0 ; i < ntestfile; i++)
  {
    printf("Test file number %i -- %s\n", i, testfile[i]);
    file = fopen(testfile[i], "r");
    if (file)
    {
      frictionContact3D_newFromFile(testproblem, file);
    }
    fclose(file);


    file = fopen("test.dat", "w");
    frictionContact3D_printInFile(testproblem, file);
    fclose(file);



    info = test_Series(testproblem);


  }





  printf("========= End Numerics tests for FrictionalContact_problem ========= \n");
  return info;
}

