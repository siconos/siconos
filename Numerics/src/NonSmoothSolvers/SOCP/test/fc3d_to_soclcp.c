/* Siconos-Numerics, Copyright INRIA 2005-2012.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#include <stdio.h>
#include <stdlib.h>
#include "NonSmoothDrivers.h"
#include "soclcp_test_function.h"


#include "../FrictionContact/test-utils/frictionContact_test_function.h"



int main(void)
{
  int info = 0 ;
  printf("Test on ./data/Example1_Fc3D_SBM.dat\n");
  FILE * finput  =  fopen("./data/Example1_Fc3D_SBM.dat", "r");

  FrictionContactProblem* problem = (FrictionContactProblem *)malloc(sizeof(FrictionContactProblem));

  info = frictionContact_newFromFile(problem, finput);


  unsigned int * coneIndex = (unsigned int *) malloc((problem->numberOfContacts+1)*sizeof(unsigned int));

  for(int i = 0; i < problem->numberOfContacts+1; i++)
  {
    coneIndex[i]=i*3;
  }
  int n = coneIndex[problem->numberOfContacts];
  SecondOrderConeLinearComplementarityProblem* soclcp =  secondOrderConeLinearComplementarityProblem_new
      (n, problem->numberOfContacts, problem->M, problem->q, coneIndex, problem->mu);


  FILE * foutput  =  fopen("./data/Example1_SOCLCP_SBM.dat", "w");
  info = secondOrderConeLinearComplementarityProblem_printInFile(soclcp, foutput);







  fclose(finput);
  printf("\nEnd of test on ./data/Example1_Fc3D_SBM.dat\n");
  return info;
}
