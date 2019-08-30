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
#include "SOCLCP_cst.h"

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
  
  int d=6;
  e=0;
  test_fp[n][e++] = data_collection[d];
  test_fp[n][e++] = "1";
  test_fp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_fp[n][e++], "%d", SICONOS_FRICTION_3D_ACLMFP);
  test_fp[n][e++] = "1e-06";
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
