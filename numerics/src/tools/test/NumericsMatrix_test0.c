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

/*
  Tests functions for NumericsMatrix structure

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NumericsMatrix.h"
#include <math.h>
#include "numericsMatrixTestFunction.h"
int main(void)
{

  printf("========= Starts Numerics tests for NumericsMatrix ========= \n");

  int i, nmm = 4 ;
  NumericsMatrix ** NMM = (NumericsMatrix **)malloc(nmm * sizeof(NumericsMatrix *)) ;
  NumericsMatrix ** Mread = (NumericsMatrix **)malloc(nmm * sizeof(NumericsMatrix *)) ;


  for (i = 0 ; i < nmm; i++)
  {
    NMM[i] = NM_new();
    Mread[i] = NM_new();
  }


  int info = test_BuildNumericsMatrix(NMM);

  if (info != 0)
  {
    printf("Construction failed ...\n");
    return info;
  }
  printf("Construction ok ...\n");

  /* Test of various I/O functions */

  for (i = 0 ; i < nmm; i++)
  {

    printf("test on NMM[%i]\n", i);

    NM_display(NMM[i]);
    NM_display_row_by_row(NMM[i]);
    FILE * foutput = fopen("testprintInfile.dat", "w");
    NM_write_in_file(NMM[i], foutput);
    fclose(foutput);
    FILE * finput = fopen("testprintInfile.dat", "r");
    NM_read_in_file(NMM[i], finput);
    fclose(finput);
    FILE * finput2 = fopen("testprintInfile.dat", "r");
    NM_new_from_file(Mread[i], finput2);
    fclose(finput2);
    char  filename[50] = "testprintInfileName.dat";
    NM_write_in_filename(NMM[i], filename);
    NM_read_in_filename(NMM[i], filename);
    printf("end of test on NMM[%i]\n", i);

  }
  for (i = 0 ; i < nmm; i++, i++)
  {
    FILE * foutput2 = fopen("testprintInfileForScilab.dat", "w");
    NM_write_in_file_scilab(NMM[i], foutput2);
    fclose(foutput2);
  }



  /* free memory */

  for (i = 0 ; i < nmm; i++)
  {
    NM_free(NMM[i]);
    free(NMM[i]);
    NM_free(Mread[i]);
    free(Mread[i]);
  }

  free(NMM);
  free(Mread);



  printf("========= End Numerics tests for NumericsMatrix ========= \n");
  return info;
}

