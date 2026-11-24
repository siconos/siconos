/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#include <stdio.h>   // for fopen, printf
#include <stdlib.h>  // for malloc

#include "FrictionContactProblem.h"                       // for FrictionCon...
#include "NumericsFwd.h"                                  // for FrictionCon...
#include "SecondOrderConeLinearComplementarityProblem.h"  // for SecondOrder...

int main(void) {
  int info = 0;
  printf("Test on ./data/FC3D_Example1_SBM.dat\n");
  const char* filename = "../../FrictionContact/test/data/Capsules-i122-1617.dat";

  FrictionContactProblem* problem = frictionContact_new_from_filename(filename);

  unsigned int* coneIndex =
      (unsigned int*)malloc((problem->numberOfContacts + 1) * sizeof(unsigned int));

  for (int i = 0; i < problem->numberOfContacts + 1; i++) {
    coneIndex[i] = i * 3;
  }
  int n = coneIndex[problem->numberOfContacts];
  SecondOrderConeLinearComplementarityProblem* soclcp =
      secondOrderConeLinearComplementarityProblem_new(n, problem->numberOfContacts, problem->M,
                                                      problem->q, coneIndex, problem->mu);

  //  secondOrderConeLinearComplementarityProblem_display(soclcp);

  FILE* foutput = fopen("./data/Capsules-i122-1617.dat", "w");
  info = secondOrderConeLinearComplementarityProblem_printInFile(soclcp, foutput);

  /* XXX should look for a better fix --xhub */
  soclcp->M = NULL;
  soclcp->q = NULL;
  soclcp->tau = NULL;
  freeSecondOrderConeLinearComplementarityProblem(soclcp);
  frictionContactProblem_free(problem);

  printf("\nEnd of test on ./data/FC3D_Example1_SBM.dat\n");
  return info;
}
