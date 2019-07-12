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
  char *** test_lcp = (char ***)malloc(n_test*sizeof(char **));

  for (int n =0 ; n <n_test ; n++)
  {
    test_lcp[n] = (char **)malloc(n_entry*sizeof(char *));
  }

  int n =0;
  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_lcp[n][e++] = data_collection_1[d];
    test_lcp[n][e++] = "0";
    test_lcp[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_lcp[n][e++], "%d", SICONOS_LCP_LEMKE);
    test_lcp[n][e++] = "---";
    n++;
  }

  /* for ( int d =0; d <n_data_1; d++) */
  /* { */
  /*   int e=0; */
  /*   test_lcp[n][e++] = data_collection_1[d]; */
  /*   test_lcp[n][e++] = "0"; */
  /*   test_lcp[n][e] = (char *)malloc(50*sizeof(char)); */
  /*   sprintf(test_lcp[n][e++], "%d", SICONOS_FRICTION_3D_NSGS); */
  /*   test_lcp[n][e++] = "1e-5"; */
  /*   test_lcp[n][e++] = "10000"; */
  /*   test_lcp[n][e] = (char *)malloc(50*sizeof(char)); */
  /*   sprintf(test_lcp[n][e++], "%d", SICONOS_FRICTION_3D_ONECONTACT_NSN_GP); */
  /*   test_lcp[n][e++] = "0.0"; */
  /*   test_lcp[n][e++] = "0"; */
  /*   test_lcp[n][e++] = "internal_iparam"; */
  /*   test_lcp[n][e++] = "10"; */
  /*   test_lcp[n][e++] = "1"; */
  /*   test_lcp[n][e++] = "---";  */
  /*   n++; */
  /* } */
  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_lcp[n][e++] = data_collection_1[d];
    test_lcp[n][e++] = "0";
    test_lcp[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_lcp[n][e++], "%d", SICONOS_LCP_CPG);
    test_lcp[n][e++] = "---";
    n++;
  }

  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_lcp[n][e++] = data_collection_1[d];
    test_lcp[n][e++] = "0";
    test_lcp[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_lcp[n][e++], "%d", SICONOS_LCP_PGS);
    test_lcp[n][e++] = "---";
    n++;
  }
  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_lcp[n][e++] = data_collection_1[d];
    test_lcp[n][e++] = "0";
    test_lcp[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_lcp[n][e++], "%d", SICONOS_LCP_RPGS);
    test_lcp[n][e++] = "---";
    n++;
  }

  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_lcp[n][e++] = data_collection_1[d];
    test_lcp[n][e++] = "0";
    test_lcp[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_lcp[n][e++], "%d", SICONOS_LCP_LATIN);
    test_lcp[n][e++] = "---";
    n++;
  }

  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_lcp[n][e++] = data_collection_1[d];
    test_lcp[n][e++] = "0";
    test_lcp[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_lcp[n][e++], "%d", SICONOS_LCP_LATIN_W);
    test_lcp[n][e++] = "---";
    n++;
  }
  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_lcp[n][e++] = data_collection_1[d];
    test_lcp[n][e++] = "0";
    test_lcp[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_lcp[n][e++], "%d", SICONOS_LCP_AVI_CAOFERRIS);
    test_lcp[n][e++] = "---";
    n++;
  }
  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_lcp[n][e++] = data_collection_1[d];
    test_lcp[n][e++] = "0";
    test_lcp[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_lcp[n][e++], "%d", SICONOS_LCP_NEWTONMIN);
    test_lcp[n][e++] = "---";
    n++;
  }
  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_lcp[n][e++] = data_collection_1[d];
    test_lcp[n][e++] = "0";
    test_lcp[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_lcp[n][e++], "%d", SICONOS_LCP_NEWTON_FB_FBLSA);
    test_lcp[n][e++] = "---";
    n++;
  }
  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_lcp[n][e++] = data_collection_1[d];
    test_lcp[n][e++] = "0";
    test_lcp[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_lcp[n][e++], "%d", SICONOS_LCP_NEWTON_MIN_FBLSA);
    test_lcp[n][e++] = "---";
    n++;
  }
  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_lcp[n][e++] = data_collection_1[d];
    test_lcp[n][e++] = "0";
    test_lcp[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_lcp[n][e++], "%d", SICONOS_LCP_BARD);
    test_lcp[n][e++] = "---";
    n++;
  }
  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_lcp[n][e++] = data_collection_1[d];
    test_lcp[n][e++] = "0";
    test_lcp[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_lcp[n][e++], "%d", SICONOS_LCP_MURTY);
    test_lcp[n][e++] = "---";
    n++;
  }
  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_lcp[n][e++] = data_collection_1[d];
    test_lcp[n][e++] = "0";
    test_lcp[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_lcp[n][e++], "%d", SICONOS_LCP_PIVOT);
    test_lcp[n][e++] = "---";
    n++;
  }
  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_lcp[n][e++] = data_collection_1[d];
    test_lcp[n][e++] = "0";
    test_lcp[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_lcp[n][e++], "%d", SICONOS_LCP_PIVOT_LUMOD);
    test_lcp[n][e++] = "---";
    n++;
  }
  for ( int d =0; d <n_data_1; d++)
  {
    int e=0;
    test_lcp[n][e++] = data_collection_1[d];
    test_lcp[n][e++] = "0";
    test_lcp[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_lcp[n][e++], "%d", SICONOS_LCP_PATHSEARCH);
    test_lcp[n][e++] = "---";
    n++;
  }

  test_lcp[n][0] ="---";
  return test_lcp;

}
