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
  char *** test_admm = (char ***)malloc(n_test*sizeof(char **));

  for (int n =0 ; n <n_test ; n++)
  {
    test_admm[n] = (char **)malloc(n_entry*sizeof(char *));
  }

  int n =0;
  for ( int d =0; d <n_data_1; d++)
  {

    int e=0;
    test_admm[n][e++] = data_collection_1[d];
    test_admm[n][e++] = "0";
    test_admm[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_admm[n][e++], "%d", SICONOS_FRICTION_3D_ADMM);
    test_admm[n][e++] = "1e-5";
    test_admm[n][e++] = "10000";
    test_admm[n][e++] = "0";
    test_admm[n][e++] = "0";
    test_admm[n][e++] = "0";
    test_admm[n][e++] = "iparam";
    test_admm[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_admm[n][e++], "%d", SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY );
    test_admm[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_admm[n][e++], "%d", SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_NORM_INF);
    test_admm[n][e++] = "---"; 
    n++;
  }
  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_admm[n][e++] = data_collection_1[d];
    test_admm[n][e++] = "0";
    test_admm[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_admm[n][e++], "%d", SICONOS_FRICTION_3D_ADMM);
    test_admm[n][e++] = "1e-5";
    test_admm[n][e++] = "10000";
    test_admm[n][e++] = "0";
    test_admm[n][e++] = "0";
    test_admm[n][e++] = "0";
    test_admm[n][e++] = "iparam";
    test_admm[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_admm[n][e++], "%d", SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY );
    test_admm[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_admm[n][e++], "%d", SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING);
    test_admm[n][e++] = "---"; 
    n++;
  }
  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_admm[n][e++] = data_collection_1[d];
    test_admm[n][e++] = "0";
    test_admm[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_admm[n][e++], "%d", SICONOS_FRICTION_3D_ADMM);
    test_admm[n][e++] = "1e-5";
    test_admm[n][e++] = "10000";
    test_admm[n][e++] = "0";
    test_admm[n][e++] = "0";
    test_admm[n][e++] = "0";
    test_admm[n][e++] = "iparam";
    test_admm[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_admm[n][e++], "%d",  SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY );
    test_admm[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_admm[n][e++], "%d", SICONOS_FRICTION_3D_ADMM_FORCED_ASYMMETRY);
    test_admm[n][e++] = "---"; 
    n++;
  }
  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_admm[n][e++] = data_collection_1[d];
    test_admm[n][e++] = "0";
    test_admm[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_admm[n][e++], "%d", SICONOS_FRICTION_3D_ADMM);
    test_admm[n][e++] = "1e-5";
    test_admm[n][e++] = "10000";
    test_admm[n][e++] = "0";
    test_admm[n][e++] = "0";
    test_admm[n][e++] = "0";
    test_admm[n][e++] = "iparam";
    test_admm[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_admm[n][e++], "%d",  SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY );
    test_admm[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_admm[n][e++], "%d", SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING);
    test_admm[n][e++] = "iparam";
    test_admm[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_admm[n][e++], "%d",  SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY );
    test_admm[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_admm[n][e++], "%d", SICONOS_FRICTION_3D_ADMM_FORCED_ASYMMETRY);
    test_admm[n][e++] = "---"; 
    n++;
  }


  

  test_admm[n][0] ="---";
  return test_admm;

}
