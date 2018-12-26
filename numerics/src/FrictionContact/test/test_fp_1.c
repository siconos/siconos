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
#include <errno.h>
char *** test_collection(int, char **);

char *** test_collection(int n_data_1, char ** data_collection)
{
  int n_test=150;
  int n_entry = 50;
  char *** test_fp = (char ***)malloc(n_test*sizeof(char **));

  for (int n =0 ; n <n_test ; n++)
  {
    test_fp[n] = (char **)malloc(n_entry*sizeof(char *));
  }

  int n =0;
  int e=0;
  
  int d=0;/* "./data/FC3D_Example1_SBM.dat"; */
  e=0;
  test_fp[n][e++] = data_collection[d];
  test_fp[n][e++] = "0";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_FRICTION_3D_TFP);
  test_fp[n][e++] = "1e-16";
  test_fp[n][e++] = "100";
  test_fp[n][e++] = "---";
  n++;
  
  e=0;
  test_fp[n][e++] = data_collection[d];
  test_fp[n][e++] = "0";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_FRICTION_3D_PFP);
  test_fp[n][e++] = "1e-16";
  test_fp[n][e++] = "100";
  test_fp[n][e++] = "---";
  n++;

  e=0;
  test_fp[n][e++] = data_collection[d];
  test_fp[n][e++] = "0";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_FRICTION_3D_TFP);
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_FRICTION_3D_ConvexQP_PG_Cylinder);
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "---";
  n++;

  e=0;
  test_fp[n][e++] = data_collection[d];
  test_fp[n][e++] = "0";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_FRICTION_3D_TFP);
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_FRICTION_3D_ConvexQP_PG_Cylinder);
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "internal_iparam";
  test_fp[n][e++] = "2";
  test_fp[n][e++] = "20";
  test_fp[n][e++] = "internal_dparam";
  test_fp[n][e++] = "3";
  test_fp[n][e++] = "-1";
  test_fp[n][e++] = "internal_dparam";
  test_fp[n][e++] = "4";
  test_fp[n][e++] = "-1.e-6";
  test_fp[n][e++] = "---";
  n++;

  
  


  
  d=2;/* "./data/Confeti-ex13-4contact-Fc3D-SBM.dat"; */
  e=0;
  test_fp[n][e++] = data_collection[d];
  test_fp[n][e++] = "1";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_FRICTION_3D_TFP);
  test_fp[n][e++] = "1e-08";
  test_fp[n][e++] = "1000";
  test_fp[n][e++] = "---";
  n++;

  e=0;
  test_fp[n][e++] = data_collection[d];
  test_fp[n][e++] = "0";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_FRICTION_3D_PFP);
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "---";
  n++;
  
  e=0;
  test_fp[n][e++] = data_collection[d];
  test_fp[n][e++] = "0";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_FRICTION_3D_TFP);
  test_fp[n][e++] = "1e-4";
  test_fp[n][e++] = "100";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_FRICTION_3D_ConvexQP_PG_Cylinder);
  test_fp[n][e++] = "1e-6";
  test_fp[n][e++] = "200";
  test_fp[n][e++] = "---";
  n++;





  d=5;/* "./data/Confeti-ex03-Fc3D-SBM.dat";  */
  
  e=0;
  test_fp[n][e++] = data_collection[d];
  test_fp[n][e++] = "0";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_FRICTION_3D_TFP);
  test_fp[n][e++] = "1e-08";
  test_fp[n][e++] = "2000";
  test_fp[n][e++] = "---";
  n++;
  
  e=0;
  test_fp[n][e++] = data_collection[d];
  test_fp[n][e++] = "0";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_FRICTION_3D_PFP);
  test_fp[n][e++] = "1e-08";
  test_fp[n][e++] = "2000";
  test_fp[n][e++] = "---";
  n++;
  
  e=0;
  test_fp[n][e++] = data_collection[d];
  test_fp[n][e++] = "0";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_FRICTION_3D_TFP);
  test_fp[n][e++] = "1e-08";
  test_fp[n][e++] = "2000";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "iparam";
  test_fp[n][e++] = "1";
  test_fp[n][e++] = "1";
  test_fp[n][e++] = "---";
  n++;



  d=6;/* "./data/Confeti-ex03-Fc3D-SBM.dat";  */
  
  e=0;
  test_fp[n][e++] = data_collection[d];
  test_fp[n][e++] = "1";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_FRICTION_3D_TFP);
  test_fp[n][e++] = "1e-08";
  test_fp[n][e++] = "100";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "dparam";
  test_fp[n][e++] = "3";
  test_fp[n][e++] = "1e4";
  test_fp[n][e++] = "---";
  n++;


  
  d=9;
  e=0;
  test_fp[n][e++] = data_collection[d];
  test_fp[n][e++] = "1";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_FRICTION_3D_TFP);
  test_fp[n][e++] = "1e-08";
  test_fp[n][e++] = "10000";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "dparam";
  test_fp[n][e++] = "3";
  test_fp[n][e++] = "1e4";
  test_fp[n][e++] = "iparam";
  test_fp[n][e++] = "1";
  test_fp[n][e++] = "1";
  test_fp[n][e++] = "---";
  n++;

  d=6;
  e=0;
  test_fp[n][e++] = data_collection[d];
  test_fp[n][e++] = "1";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_FRICTION_3D_DSFP);
  test_fp[n][e++] = "1e-03";
  test_fp[n][e++] = "100000";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "dparam";
  test_fp[n][e++] = "3";
  test_fp[n][e++] = "8e4";
  test_fp[n][e++] = "---";
  n++;

  d=5;
  e=0;
  test_fp[n][e++] = data_collection[d];
  test_fp[n][e++] = "1";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_FRICTION_3D_DSFP);
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "dparam";
  test_fp[n][e++] = "3";
  test_fp[n][e++] = "1e2";
  test_fp[n][e++] = "---";
  n++;
  
  d=2;
  e=0;
  test_fp[n][e++] = data_collection[d];
  test_fp[n][e++] = "1";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_FRICTION_3D_DSFP);
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "dparam";
  test_fp[n][e++] = "3";
  test_fp[n][e++] = "5e3";
  test_fp[n][e++] = "---";
  n++;

  d=0;
  e=0;
  test_fp[n][e++] = data_collection[d];
  test_fp[n][e++] = "1";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_FRICTION_3D_DSFP);
  test_fp[n][e++] = "1e-08";
  test_fp[n][e++] = "100000";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "dparam";
  test_fp[n][e++] = "3";
  test_fp[n][e++] = "2.0";
  test_fp[n][e++] = "---";
  n++;


  d=6;
  e=0;
  test_fp[n][e++] = data_collection[d];
  test_fp[n][e++] = "1";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_FRICTION_3D_ACLMFP);
  test_fp[n][e++] = "1e-08";
  test_fp[n][e++] = "200";
  test_fp[n][e++] = "---";
  n++;

  d=5;
  e=0;
  test_fp[n][e++] = data_collection[d];
  test_fp[n][e++] = "1";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_FRICTION_3D_ACLMFP);
  test_fp[n][e++] = "1e-08";
  test_fp[n][e++] = "200";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "iparam";
  test_fp[n][e++] = "1";
  test_fp[n][e++] = "1";
  test_fp[n][e++] = "---";
  n++;
  
  d=2;
  e=0;
  test_fp[n][e++] = data_collection[d];
  test_fp[n][e++] = "1";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_FRICTION_3D_ACLMFP);
  test_fp[n][e++] = "1e-08";
  test_fp[n][e++] = "200";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_SOCLCP_VI_FPP);
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "---";
  n++;

  d=0;
  e=0;
  test_fp[n][e++] = data_collection[d];
  test_fp[n][e++] = "1";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_FRICTION_3D_ACLMFP);
  test_fp[n][e++] = "1e-08";
  test_fp[n][e++] = "200";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_SOCLCP_VI_EG);
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "---";
  n++;


  d=5;
  e=0;
  test_fp[n][e++] = data_collection[d];
  test_fp[n][e++] = "1";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_FRICTION_3D_SOCLCP);
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "0";
  test_fp[n][e++] = "---";
  n++;
  
  test_fp[n][0] ="---";
  return test_fp;

}
