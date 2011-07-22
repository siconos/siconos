/* Siconos-Numerics, Copyright INRIA 2005-2011.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/

/*
  Tests functions for NumericsMatrix structure

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NumericsMatrix.h"
#include "LA.h"
#include <math.h>
#include "numericsMatrixTestFunction.h"
int main(void)
{

  printf("========= Starts Numerics tests for NumericsMatrix ========= \n");

  int i, nmm = 4 ;
  NumericsMatrix ** NMM = malloc(nmm * sizeof(NumericsMatrix *)) ;
  NumericsMatrix ** Mread = malloc(nmm * sizeof(NumericsMatrix *)) ;


  for (i = 0 ; i < nmm; i++)
  {
    NMM[i] = malloc(sizeof(NumericsMatrix));
    Mread[i] = malloc(sizeof(NumericsMatrix));
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

    display(NMM[i]);
    displayRawbyRaw(NMM[i]);
    FILE * foutput = fopen("testprintInfile.dat", "w");
    printInFile(NMM[i], foutput);
    fclose(foutput);
    FILE * finput = fopen("testprintInfile.dat", "r");
    readInFile(NMM[i], finput);
    fclose(finput);
    FILE * finput2 = fopen("testprintInfile.dat", "r");
    newFromFile(Mread[i], finput2);
    fclose(finput2);
    char  filename[50] = "testprintInfileName.dat";
    printInFileName(NMM[i], filename);
    readInFileName(NMM[i], filename);
    printf("end of test on NMM[%i]\n", i);

  }
  for (i = 0 ; i < nmm; i++, i++)
  {
    FILE * foutput2 = fopen("testprintInfileForScilab.dat", "w");
    printInFileForScilab(NMM[i], foutput2);
    fclose(foutput2);
  }



  /* free memory */

  for (i = 0 ; i < nmm; i++)
  {
    freeNumericsMatrix(NMM[i]);
    free(NMM[i]);
    freeNumericsMatrix(Mread[i]);
    free(Mread[i]);
  }

  free(NMM);
  free(Mread);



  printf("========= End Numerics tests for NumericsMatrix ========= \n");
  return info;
}

