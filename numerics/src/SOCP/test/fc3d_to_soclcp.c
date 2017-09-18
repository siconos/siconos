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
#include <stdio.h>
#include <stdlib.h>
#include "NonSmoothDrivers.h"
#include "soclcp_test_function.h"
#include "FrictionContactProblem.h"
#include "SecondOrderConeLinearComplementarityProblem.h"

#include "../FrictionContact/test-utils/frictionContact_test_function.h"



int main(void)
{
  int info = 0 ;
  printf("Test on ./data/FC3D_Example1_SBM.dat\n");
  FILE * finput  =  fopen("../../FrictionContact/test/data/Capsules-i122-1617.dat", "r");

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

//  secondOrderConeLinearComplementarityProblem_display(soclcp);

  FILE * foutput  =  fopen("./data/Capsules-i122-1617.dat", "w");
  info = secondOrderConeLinearComplementarityProblem_printInFile(soclcp, foutput);





  /* XXX should look for a better fix --xhub */
  soclcp->M = NULL;
  soclcp->q = NULL;
  soclcp->tau = NULL;
  freeSecondOrderConeLinearComplementarityProblem(soclcp);
  freeFrictionContactProblem(problem);


  fclose(finput);
  printf("\nEnd of test on ./data/FC3D_Example1_SBM.dat\n");
  return info;
}
