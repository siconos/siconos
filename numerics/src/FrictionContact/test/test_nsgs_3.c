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

#include "frictionContact_test_utils.h"

char *** test_collection(int n_data_1, char ** data_collection)
{
  int n_test=150;
  int n_entry = 50;
  char *** test_nsgs = (char ***)malloc(n_test*sizeof(char **));

  for (int n =0 ; n <n_test ; n++)
  {
    test_nsgs[n] = (char **)malloc(n_entry*sizeof(char *));
  }

  int n =0;
  int e=0;
  int d=0;/* "./data/FC3D_Example1_SBM.dat"; */
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "0";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGS);
  test_nsgs[n][e++] = "1e-16";
  test_nsgs[n][e++] = "10000";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone);
  test_nsgs[n][e++] = "0.0";
  test_nsgs[n][e++] = "0";

  test_nsgs[n][e++] = "---";
  n++;

  e=0;
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "0";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGS);
  test_nsgs[n][e++] = "1e-16";
  test_nsgs[n][e++] = "10000";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithDiagonalization);
  test_nsgs[n][e++] = "0.0";
  test_nsgs[n][e++] = "0";
  test_nsgs[n][e++] = "---";
  n++;

  e=0;
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "0";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGS);
  test_nsgs[n][e++] = "1e-16";
  test_nsgs[n][e++] = "10000";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);
  test_nsgs[n][e++] = "1e-3";
  test_nsgs[n][e++] = "10";
  test_nsgs[n][e++] = "---";
  n++;

  e=0;
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "0";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGS);
  test_nsgs[n][e++] = "1e-16";
  test_nsgs[n][e++] = "10000";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization);
  test_nsgs[n][e++] = "0.0";
  test_nsgs[n][e++] = "0.0";
  test_nsgs[n][e++] = "dparam";
  test_nsgs[n][e++] = "3";
  test_nsgs[n][e++] = "0.1";
  test_nsgs[n][e++] = "---";
  n++;

  /* e=0; */
  /* test_nsgs[n][e++] = data_collection[d]; */
  /* test_nsgs[n][e++] = "1"; */
  /* test_nsgs[n][e] = (char *)malloc(50*sizeof(char)); */
  /* sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGS); */
  /* test_nsgs[n][e++] = "1e-16"; */
  /* test_nsgs[n][e++] = "10000"; */
  /* test_nsgs[n][e] = (char *)malloc(50*sizeof(char)); */
  /* sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint); */
  /* test_nsgs[n][e++] = "0.0"; */
  /* test_nsgs[n][e++] = "10"; */
  /* test_nsgs[n][e++] = "---";  */
  /* n++; */

  /* e=0; */
  /* test_nsgs[n][e++] = data_collection[d]; */
  /* test_nsgs[n][e++] = "1"; */
  /* test_nsgs[n][e] = (char *)malloc(50*sizeof(char)); */
  /* sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGSV); */
  /* test_nsgs[n][e++] = "1e-05"; */
  /* test_nsgs[n][e++] = "10000"; */
  /* test_nsgs[n][e] = (char *)malloc(50*sizeof(char)); */
  /* sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone_velocity); */
  /* test_nsgs[n][e++] = "0.0"; */
  /* test_nsgs[n][e++] = "0.0"; */
  /* test_nsgs[n][e++] = "internal_iparam"; */
  /* test_nsgs[n][e++] = "0"; */
  /* test_nsgs[n][e++] = "0"; */
  /* test_nsgs[n][e++] = "internal_dparam"; */
  /* test_nsgs[n][e++] = "0"; */
  /* test_nsgs[n][e++] = "0"; */
  /* test_nsgs[n][e++] = "---"; */
  /* n++; */


  
  d=1;/* "./data/Capsules-i122-1617.dat"; */
  e=0;
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "1";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGS);
  test_nsgs[n][e++] = "1e-16";
  test_nsgs[n][e++] = "10000";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);
  test_nsgs[n][e++] = "1e-16";
  test_nsgs[n][e++] = "20";
  test_nsgs[n][e++] = "dparam";
  test_nsgs[n][e++] = "9";
  test_nsgs[n][e++] = "1.0";
  test_nsgs[n][e++] = "iparam";
  test_nsgs[n][e++] = "8";
  test_nsgs[n][e++] = "1";
  test_nsgs[n][e++] = "---";
  n++;

  d=2;/* "./data/Confeti-ex13-4contact-Fc3D-SBM.dat"; */
  e=0;
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "1";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGS);
  test_nsgs[n][e++] = "1e-05";
  test_nsgs[n][e++] = "10000";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone);
  test_nsgs[n][e++] = "0.0";
  test_nsgs[n][e++] = "0";
  test_nsgs[n][e++] = "dparam";
  test_nsgs[n][e++] = "9";
  test_nsgs[n][e++] = "1.0";
  test_nsgs[n][e++] = "iparam";
  test_nsgs[n][e++] = "8";
  test_nsgs[n][e++] = "1";
  test_nsgs[n][e++] = "---";
  n++;

  e=0;
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "0";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGS);
  test_nsgs[n][e++] = "1e-12";
  test_nsgs[n][e++] = "10000";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_NSN);
  test_nsgs[n][e++] = "1e-18";
  test_nsgs[n][e++] = "10";
  test_nsgs[n][e++] = "internal_iparam";
  test_nsgs[n][e++] = "1";
  test_nsgs[n][e++] = "1";
  test_nsgs[n][e++] = "---";

  n++;
  e=0;
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "0";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGS);
  test_nsgs[n][e++] = "1e-12";
  test_nsgs[n][e++] = "10000";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);
  test_nsgs[n][e++] = "1e-06";
  test_nsgs[n][e++] = "100";
  test_nsgs[n][e++] = "---";
  n++;


  e=0;
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "0";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGS);
  test_nsgs[n][e++] = "1e-12";
  test_nsgs[n][e++] = "10000";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization);
  test_nsgs[n][e++] = "0.0";
  test_nsgs[n][e++] = "0";
  test_nsgs[n][e++] = "---";
  n++;

  e=0;
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "0";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGS);
  test_nsgs[n][e++] = "1e-02";
  test_nsgs[n][e++] = "10000";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone);
  test_nsgs[n][e++] = "0.0";
  test_nsgs[n][e++] = "0";
  test_nsgs[n][e++] = "---";
  n++;

  e=0;
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "0";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGS);
  test_nsgs[n][e++] = "1e-05";
  test_nsgs[n][e++] = "1000";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_NSN);
  test_nsgs[n][e++] = "1e-16";
  test_nsgs[n][e++] = "10";
  test_nsgs[n][e++] = "---";
  n++;

   e=0;
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "0";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGS);
  test_nsgs[n][e++] = "1e-12";
  test_nsgs[n][e++] = "10000";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);
  test_nsgs[n][e++] = "1e-06";
  test_nsgs[n][e++] = "100";
  test_nsgs[n][e++] = "---";
  n++;
  
  e=0;
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "0";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGS);
  test_nsgs[n][e++] = "1e-12";
  test_nsgs[n][e++] = "10000";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);
  test_nsgs[n][e++] = "1e-16";
  test_nsgs[n][e++] = "100";
  test_nsgs[n][e++] = "---";
  n++;

  
  
  






  
  d=3 ; /* "./data/GFC3D_TwoRods1-condensed.dat"; */
  e=0;
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "0";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGS);
  test_nsgs[n][e++] = "1e-12";
  test_nsgs[n][e++] = "10000";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_NSN);
  test_nsgs[n][e++] = "1e-18";
  test_nsgs[n][e++] = "10";
  test_nsgs[n][e++] = "internal_iparam";
  test_nsgs[n][e++] = "1";
  test_nsgs[n][e++] = "1";
  test_nsgs[n][e++] = "---";


  d=4 ; /* "./data/FC3D_Example1.dat"; */
  e=0;
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "0";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGS);
  test_nsgs[n][e++] = "1e-12";
  test_nsgs[n][e++] = "10000";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_NSN);
  test_nsgs[n][e++] = "1e-18";
  test_nsgs[n][e++] = "10";




  d=5 ; /* "./data/Confeti-ex03-Fc3D-SBM.dat"; */


  e=0;
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "1";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGS);
  test_nsgs[n][e++] = "1e-05";
  test_nsgs[n][e++] = "10000";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone);
  test_nsgs[n][e++] = "0.0";
  test_nsgs[n][e++] = "0";
  test_nsgs[n][e++] = "---";
  n++;

  e=0;
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "1";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGS);
  test_nsgs[n][e++] = "1e-05";
  test_nsgs[n][e++] = "10000";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_NSN);
  test_nsgs[n][e++] = "1e-16";
  test_nsgs[n][e++] = "10";
  test_nsgs[n][e++] = "---";
  n++;

  e=0;
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "0";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGS);
  test_nsgs[n][e++] = "1e-05";
  test_nsgs[n][e++] = "10000";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);
  test_nsgs[n][e++] = "1e-12";
  test_nsgs[n][e++] = "10";
  test_nsgs[n][e++] = "---";
  n++;

  
  e=0;
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "1";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGS);
  test_nsgs[n][e++] = "1e-05";
  test_nsgs[n][e++] = "10000";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization);
  test_nsgs[n][e++] = "1e-08";
  test_nsgs[n][e++] = "10";
  test_nsgs[n][e++] = "---";
  n++;


  d=7;
  e=0;
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "1";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGS);
  test_nsgs[n][e++] = "1e-03";
  test_nsgs[n][e++] = "1000";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d",  SICONOS_FRICTION_3D_ONECONTACT_NSN_GP);
  test_nsgs[n][e++] = "1e-16";
  test_nsgs[n][e++] = "100";
  test_nsgs[n][e++] = "---";
  n++;
  
  e=0;
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "1";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGS);
  test_nsgs[n][e++] = "1e-03";
  test_nsgs[n][e++] = "1000";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d",  SICONOS_FRICTION_3D_ONECONTACT_NSN_GP);
  test_nsgs[n][e++] = "1e-16";
  test_nsgs[n][e++] = "1000";
  test_nsgs[n][e++] = "---";
  n++;
  
  e=0;
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "1";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGS);
  test_nsgs[n][e++] = "1e-03";
  test_nsgs[n][e++] = "2000";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d",  SICONOS_FRICTION_3D_ONECONTACT_NSN);
  test_nsgs[n][e++] = "1e-16";
  test_nsgs[n][e++] = "100";
  test_nsgs[n][e++] = "---";
  n++;
  
  e=0;
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "1";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGS);
  test_nsgs[n][e++] = "1e-03";
  test_nsgs[n][e++] = "2000";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d",  SICONOS_FRICTION_3D_ONECONTACT_NSN);
  test_nsgs[n][e++] = "1e-16";
  test_nsgs[n][e++] = "1000";
  test_nsgs[n][e++] = "---";
  n++;
  
  e=0;
  test_nsgs[n][e++] = data_collection[d];
  test_nsgs[n][e++] = "1";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_NSGS);
  test_nsgs[n][e++] = "1e-03";
  test_nsgs[n][e++] = "2000";
  test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_nsgs[n][e++], "%d",  SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);
  test_nsgs[n][e++] = "1e-06";
  test_nsgs[n][e++] = "100";
  test_nsgs[n][e++] = "---";
  n++;
  
  test_nsgs[n][0] ="---";
  return test_nsgs;

}
