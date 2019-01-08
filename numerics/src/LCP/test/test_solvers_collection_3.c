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
  for ( int d =0; d <n_data_1-1; d++)
  {
    int e=0;
    test_lcp[n][e++] = data_collection_1[d];
    test_lcp[n][e++] = "0";
    test_lcp[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_lcp[n][e++], "%d", SICONOS_LCP_ENUM);
    test_lcp[n][e++] = "---"; 
    n++;
  }
  for ( int d =0; d <n_data_1-1; d++)
  {
    int e=0;
    test_lcp[n][e++] = data_collection_1[d];
    test_lcp[n][e++] = "0";
    test_lcp[n][e] = (char *)malloc(50*sizeof(char));
    sprintf(test_lcp[n][e++], "%d", SICONOS_LCP_LEMKE);
    test_lcp[n][e++] = "---"; 
    n++;
  }


  int d=3 ; /* tobenna */
  int e =0;
  test_lcp[n][e++] = data_collection_1[d];
  test_lcp[n][e++] = "0";
  test_lcp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_lcp[n][e++], "%d", SICONOS_LCP_PIVOT);
  test_lcp[n][e++] = "---"; 
  n++;


#ifdef HAVE_PATHFERRIS
  e =0;
  test_lcp[n][e++] = data_collection_1[d];
  test_lcp[n][e++] = "0";
  test_lcp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_lcp[n][e++], "%d", SICONOS_PATH);
  test_lcp[n][e++] = "---"; 
  n++;
#endif
#ifdef HAVE_GAMS_C_API
  e =0;
  test_lcp[n][e++] = data_collection_1[d];
  test_lcp[n][e++] = "0";
  test_lcp[n][e] = (char *)malloc(50*sizeof(char));
  sprintf(test_lcp[n][e++], "%d", SICONOS_GAMS);
  test_lcp[n][e++] = "---"; 
  n++;
#endif
  
  

  test_lcp[1][1]="1";
  test_lcp[4][1]="1";
  test_lcp[6][1]="1";

  test_lcp[n][0] ="---";
  return test_lcp;

}
