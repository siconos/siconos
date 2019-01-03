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

char *** test_collection(int n_data_1, char ** data_collection_1)
{
  int n_test=150;
  int n_entry = 50;
  char *** test_nsgs = (char ***)malloc(n_test*sizeof(char **));

  for (int n =0 ; n <n_test ; n++)
  {
    test_nsgs[n] = (char **)malloc(n_entry*sizeof(char *));
  }

  int n =0;
  
  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_nsgs[n][e++] = data_collection_1[d];
    test_nsgs[n][e++] = "0";
    test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_nsgs[n][e++], "%d", SICONOS_GENERIC_MECHANICAL_NSGS);
    test_nsgs[n][e++] = "1e-5";
    test_nsgs[n][e++] = "10000";
    test_nsgs[n][e++] = "0";
    test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_QUARTIC);
    test_nsgs[n][e++] = "0.0";
    test_nsgs[n][e++] = "0";
    test_nsgs[n][e++] = "---";
    n++;
  }


#ifdef HAS_LAPACK_DGESVD
  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_nsgs[n][e++] = data_collection_1[d];
    test_nsgs[n][e++] = "0";
    test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_nsgs[n][e++], "%d", SICONOS_GENERIC_MECHANICAL_NSGS);
    test_nsgs[n][e++] = "1e-5";
    test_nsgs[n][e++] = "10000";
    test_nsgs[n][e++] = "1";
    test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_QUARTIC);
    test_nsgs[n][e++] = "0.0";
    test_nsgs[n][e++] = "0";
    test_nsgs[n][e++] = "---";
    n++;
  }

#endif
  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_nsgs[n][e++] = data_collection_1[d];
    test_nsgs[n][e++] = "0";
    test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_nsgs[n][e++], "%d", SICONOS_GENERIC_MECHANICAL_NSGS);
    test_nsgs[n][e++] = "1e-5";
    test_nsgs[n][e++] = "10000";
    test_nsgs[n][e++] = "2";
    test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_QUARTIC);
    test_nsgs[n][e++] = "0.0";
    test_nsgs[n][e++] = "0";
    test_nsgs[n][e++] = "---";
    n++;
  }
  for ( int d =0; d <n_data_1-3; d++)
  {
    int e=0;
    test_nsgs[n][e++] = data_collection_1[d];
    test_nsgs[n][e++] = "0";
    test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_nsgs[n][e++], "%d", SICONOS_GENERIC_MECHANICAL_NSGS);
    test_nsgs[n][e++] = "1e-5";
    test_nsgs[n][e++] = "10000";
    test_nsgs[n][e++] = "3";
    test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_QUARTIC);
    test_nsgs[n][e++] = "0.0";
    test_nsgs[n][e++] = "0";
    test_nsgs[n][e++] = "---";
    n++;
  }
  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_nsgs[n][e++] = data_collection_1[d];
    test_nsgs[n][e++] = "0";
    test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_nsgs[n][e++], "%d", SICONOS_GENERIC_MECHANICAL_NSGS);
    test_nsgs[n][e++] = "1e-5";
    test_nsgs[n][e++] = "10000";
    test_nsgs[n][e++] = "0";
    test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_NSN_GP);
    test_nsgs[n][e++] = "0.0";
    test_nsgs[n][e++] = "0";
    test_nsgs[n][e++] = "---";
    n++;
  }
  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_nsgs[n][e++] = data_collection_1[d];
    test_nsgs[n][e++] = "0";
    test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_nsgs[n][e++], "%d", SICONOS_GENERIC_MECHANICAL_NSGS);
    test_nsgs[n][e++] = "1e-5";
    test_nsgs[n][e++] = "10000";
    test_nsgs[n][e++] = "2";
    test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_NSN_GP);
    test_nsgs[n][e++] = "0.0";
    test_nsgs[n][e++] = "0";
    test_nsgs[n][e++] = "---";
    n++;
  }
  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_nsgs[n][e++] = data_collection_1[d];
    test_nsgs[n][e++] = "0";
    test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_nsgs[n][e++], "%d", SICONOS_GENERIC_MECHANICAL_NSGS);
    test_nsgs[n][e++] = "1e-5";
    test_nsgs[n][e++] = "10000";
    test_nsgs[n][e++] = "0";
    test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_NSN);
    test_nsgs[n][e++] = "0.0";
    test_nsgs[n][e++] = "0";
    test_nsgs[n][e++] = "---";
    n++;
  }
  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_nsgs[n][e++] = data_collection_1[d];
    test_nsgs[n][e++] = "0";
    test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_nsgs[n][e++], "%d", SICONOS_GENERIC_MECHANICAL_NSGS);
    test_nsgs[n][e++] = "1e-5";
    test_nsgs[n][e++] = "10000";
    test_nsgs[n][e++] = "2";
    test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_NSN);
    test_nsgs[n][e++] = "0.0";
    test_nsgs[n][e++] = "0";
    test_nsgs[n][e++] = "---";
    n++;
  }

  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_nsgs[n][e++] = data_collection_1[d];
    test_nsgs[n][e++] = "0";
    test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_nsgs[n][e++], "%d", SICONOS_GENERIC_MECHANICAL_NSGS);
    test_nsgs[n][e++] = "1e-5";
    test_nsgs[n][e++] = "10000";
    test_nsgs[n][e++] = "0";
    test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);
    test_nsgs[n][e++] = "0.0";
    test_nsgs[n][e++] = "0";
    test_nsgs[n][e++] = "---";
    n++;
  }
  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_nsgs[n][e++] = data_collection_1[d];
    test_nsgs[n][e++] = "0";
    test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_nsgs[n][e++], "%d", SICONOS_GENERIC_MECHANICAL_NSGS);
    test_nsgs[n][e++] = "1e-5";
    test_nsgs[n][e++] = "10000";
    test_nsgs[n][e++] = "2";
    test_nsgs[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_nsgs[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);
    test_nsgs[n][e++] = "0.0";
    test_nsgs[n][e++] = "0";
    test_nsgs[n][e++] = "---";
    n++;
  }
#ifdef HAS_LAPACK_DGESVD
  test_nsgs[5][1] = "1";
  test_nsgs[6][1] = "1";
  test_nsgs[19][1] = "1";
  test_nsgs[31][1] = "1";

  test_nsgs[45][1] = "1";
  test_nsgs[58][1] = "1";
  test_nsgs[59][1] = "1";
  test_nsgs[65][1] = "1";
#else
  test_nsgs[5][1] = "1";
  test_nsgs[6][1] = "1";
  test_nsgs[12][1] = "1";
  test_nsgs[24][1] = "1";
  test_nsgs[38][1] = "1";
  test_nsgs[51][1] = "1";
  test_nsgs[52][1] = "1";
  test_nsgs[58][1] = "1";
  
#endif





  
  test_nsgs[n][0] ="---";
  return test_nsgs;

}
