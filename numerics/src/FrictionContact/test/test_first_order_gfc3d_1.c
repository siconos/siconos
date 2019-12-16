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

char *** test_collection(int n_data_1, char ** data_collection_1)
{
  int n_test=150;
  int n_entry = 50;
  char *** test = (char ***)malloc(n_test*sizeof(char **));

  for (int n =0 ; n <n_test ; n++)
  {
    test[n] = (char **)malloc(n_entry*sizeof(char *));
  }

  int n =0;
  /* for ( int d =0; d <n_data_1; d++) */
  /* { */
  /*   int e=0; */
  /*   test[n][e++] = data_collection_1[d]; */
  /*   test[n][e++] = "0"; */
  /*   test[n][e] = (char *)malloc(50*sizeof(char)); */
  /*   sprintf(test[n][e++], "%d", SICONOS_GLOBAL_FRICTION_3D_NSGS); */
  /*   test[n][e++] = "0"; */
  /*   test[n][e++] = "10000"; */
  /*   test[n][e++] = "---"; */
  /*   n++; */
  /* } */
  /* for ( int d =0; d <n_data_1; d++) */
  /* { */
  /*   int e=0; */
  /*   test[n][e++] = data_collection_1[d]; */
  /*   test[n][e++] = "0"; */
  /*   test[n][e] = (char *)malloc(50*sizeof(char)); */
  /*   sprintf(test[n][e++], "%d", SICONOS_GLOBAL_FRICTION_3D_NSGS); */
  /*   test[n][e++] = "0"; */
  /*   test[n][e++] = "10000"; */
  /*   test[n][e] = (char *)malloc(50*sizeof(char)); */
  /*   sprintf(test[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone); */
  /*   test[n][e++] = "0.0"; */
  /*   test[n][e++] = "0"; */
  /*   test[n][e++] = "internal_iparam"; */
  /*   test[n][e++] = "0"; */
  /*   test[n][e++] = "0"; */
  /*   test[n][e++] = "internal_dparam"; */
  /*   test[n][e++] = "0"; */
  /*   test[n][e++] = "0"; */
  /*   test[n][e++] = "---"; */
  /*   n++; */
  /* } */
  /* for ( int d =0; d <n_data_1; d++) */
  /* { */
  /*   int e=0; */
  /*   test[n][e++] = data_collection_1[d]; */
  /*   test[n][e++] = "0"; */
  /*   test[n][e] = (char *)malloc(50*sizeof(char)); */
  /*   sprintf(test[n][e++], "%d", SICONOS_GLOBAL_FRICTION_3D_VI_EG); */
  /*   test[n][e++] = "0"; */
  /*   test[n][e++] = "0"; */
  /*   test[n][e++] = "---"; */
  /*   n++; */
  /* } */

  /* for ( int d =0; d <n_data_1; d++) */
  /* { */
  /*   int e=0; */
  /*   test[n][e++] = data_collection_1[d]; */
  /*   test[n][e++] = "0"; */
  /*   test[n][e] = (char *)malloc(50*sizeof(char)); */
  /*   sprintf(test[n][e++], "%d", SICONOS_GLOBAL_FRICTION_3D_VI_FPP); */
  /*   test[n][e++] = "40000"; */
  /*   test[n][e++] = "0"; */
  /*   test[n][e++] = "---"; */
  /*   n++; */
  /* } */

  /* for ( int d =0; d <n_data_1; d++) */
  /* { */
  /*   int e=0; */
  /*   test[n][e++] = data_collection_1[d]; */
  /*   test[n][e++] = "0"; */
  /*   test[n][e] = (char *)malloc(50*sizeof(char)); */
  /*   sprintf(test[n][e++], "%d", SICONOS_GLOBAL_FRICTION_3D_ACLMFP); */
  /*   test[n][e++] = "1e-5"; */
  /*   test[n][e++] = "1000"; */
  /*   test[n][e++] = "---"; */
  /*   n++; */
  /* } */
  /* for ( int d =0; d <n_data_1; d++) */
  /* { */
  /*   int e=0; */
  /*   test[n][e++] = data_collection_1[d]; */
  /*   test[n][e++] = "0"; */
  /*   test[n][e] = (char *)malloc(50*sizeof(char)); */
  /*   sprintf(test[n][e++], "%d", SICONOS_GLOBAL_FRICTION_3D_ADMM); */
  /*   test[n][e++] = "1e-12"; */
  /*   test[n][e++] = "100000"; */
  /*   test[n][e++] = "---"; */
  /*   n++; */
  /* } */
  /* for ( int d =0; d <n_data_1; d++) */
  /* { */
  /*   int e=0; */
  /*   test[n][e++] = data_collection_1[d]; */
  /*   test[n][e++] = "0"; */
  /*   test[n][e] = (char *)malloc(50*sizeof(char)); */
  /*   sprintf(test[n][e++], "%d", SICONOS_GLOBAL_FRICTION_3D_ADMM); */
  /*   test[n][e++] = "1e-12"; */
  /*   test[n][e++] = "100000"; */
  /*   test[n][e++] = "0"; */
  /*   test[n][e++] = "0"; */
  /*   test[n][e++] = "0"; */
  /*   test[n][e++] = "iparam"; */
  /*   test[n][e] = (char *)malloc(50*sizeof(char)); */
  /*   sprintf(test[n][e++], "%d", SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY ); */
  /*   test[n][e] = (char *)malloc(50*sizeof(char)); */
  /*   sprintf(test[n][e++], "%d", SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING); */
  /*   test[n][e++] = "---"; */
  /*   n++; */
  /* } */
  /* for ( int d =0; d <n_data_1; d++) */
  /* { */
  /*   int e=0; */
  /*   test[n][e++] = data_collection_1[d]; */
  /*   test[n][e++] = "0"; */
  /*   test[n][e] = (char *)malloc(50*sizeof(char)); */
  /*   sprintf(test[n][e++], "%d", SICONOS_GLOBAL_FRICTION_3D_ADMM); */
  /*   test[n][e++] = "1e-12"; */
  /*   test[n][e++] = "100000"; */
  /*   test[n][e++] = "0"; */
  /*   test[n][e++] = "0"; */
  /*   test[n][e++] = "0"; */
  /*   test[n][e++] = "iparam"; */
  /*   test[n][e] = (char *)malloc(50*sizeof(char)); */
  /*   sprintf(test[n][e++], "%d",SICONOS_FRICTION_3D_ADMM_IPARAM_UPDATE_S); */
  /*   test[n][e] = (char *)malloc(50*sizeof(char)); */
  /*   sprintf(test[n][e++], "%d",SICONOS_FRICTION_3D_ADMM_UPDATE_S_NO); */
  /*   test[n][e++] = "---"; */
  /*   n++; */
  /* } */
  /* for ( int d =0; d <n_data_1; d++) */
  /* { */
  /*   int e=0; */
  /*   test[n][e++] = data_collection_1[d]; */
  /*   test[n][e++] = "0"; */
  /*   test[n][e] = (char *)malloc(50*sizeof(char)); */
  /*   sprintf(test[n][e++], "%d", SICONOS_GLOBAL_FRICTION_3D_ADMM); */
  /*   test[n][e++] = "1e-12"; */
  /*   test[n][e++] = "1000"; */
  /*   test[n][e++] = "0"; */
  /*   test[n][e++] = "0"; */
  /*   test[n][e++] = "0"; */
  /*   test[n][e++] = "iparam"; */
  /*   test[n][e] = (char *)malloc(50*sizeof(char)); */
  /*   sprintf(test[n][e++], "%d", SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY); */
  /*   test[n][e] = (char *)malloc(50*sizeof(char)); */
  /*   sprintf(test[n][e++], "%d", SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING); */
  /*   test[n][e++] = "iparam"; */
  /*   test[n][e] = (char *)malloc(50*sizeof(char)); */
  /*   sprintf(test[n][e++], "%d", SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION); */
  /*   test[n][e] = (char *)malloc(50*sizeof(char)); */
  /*   sprintf(test[n][e++], "%d", SICONOS_FRICTION_3D_ADMM_NO_ACCELERATION); */
  /*   test[n][e++] = "---";  */
  /*   n++; */
  /* } */
    for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test[n][e++] = data_collection_1[d];
    test[n][e++] = "0";
    test[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test[n][e++], "%d", SICONOS_GLOBAL_FRICTION_3D_ADMM);
    test[n][e++] = "1e-08";
    test[n][e++] = "10000";
    test[n][e++] = "0";
    test[n][e++] = "0";
    test[n][e++] = "0";
    test[n][e++] = "iparam";
    test[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test[n][e++], "%d", SICONOS_FRICTION_3D_ADMM_IPARAM_UPDATE_S);
    test[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test[n][e++], "%d", SICONOS_FRICTION_3D_ADMM_UPDATE_S_NO);
    test[n][e++] = "iparam";
    test[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test[n][e++], "%d", SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY );
    test[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test[n][e++], "%d", SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING);
    test[n][e++] = "iparam";
    test[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test[n][e++], "%d", SICONOS_FRICTION_3D_IPARAM_RESCALING);
    test[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test[n][e++], "%d", SICONOS_FRICTION_3D_RESCALING_YES);
    test[n][e++] = "---";
    n++;
  }
    
    
  test[n][0] ="---";
  return test;

}
