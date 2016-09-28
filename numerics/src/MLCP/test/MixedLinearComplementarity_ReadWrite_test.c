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
#include <string.h>
#include <assert.h>
#include "MixedLinearComplementarityProblem.h"
#include "SiconosCompat.h"

int write_newformat(char *filename);
int write_newformat(char *filename)
{
  printf("\n Start of test \n");
  printf("Test on %s\n", filename);
  int info = 0;
  size_t sizeoffilename = strlen(filename);
  printf("sizeoffilename " SN_SIZE_T_F "\n",  sizeoffilename);
  char  extension[4] = "ext";
  strncpy(extension, &filename[sizeoffilename - 3], 3);
  printf("extension %s\n",  extension);
  char * basename;

  if (strcmp(extension, "dat") == 0)
  {
    basename = (char *)malloc((sizeoffilename + 5) * sizeof(char));
    basename[sizeoffilename + 4] = 0;
    strcpy(basename, filename);
    strncpy(&basename[sizeoffilename - 4], ".dat.tmp", 8);
    printf("basename %s\n",  basename);

  }
  else
  {
    printf("Wrong file name extension %s\n", extension);
    return 0;
  }
  /* Remove file if it exists */
  FILE * foutput = fopen(basename, "w");
  fclose(foutput);

  FILE * f  =  fopen(filename, "r");
  MixedLinearComplementarityProblem * problem = (MixedLinearComplementarityProblem *)malloc(sizeof(MixedLinearComplementarityProblem));

  info = mixedLinearComplementarity_newFromFileOld(problem, f);
  mixedLinearComplementarity_display(problem);

  char  val[128];
  int t = 0;
  t++;
  t = fscanf(f , "%s" , val);
  int withSol = 0;
  int n = problem->n;
  int m = problem ->m;
  assert(n>0);
  assert(m>0);
  int i;
  double *sol  = (double*)malloc((n + m + m) * sizeof(double));

  if (!feof(f))
  {
    withSol = 1;


    sol[0] = atof(val);

    for (i = 1 ; i < n + m + m ; ++i)
    {
      t = fscanf(f , "%s" , val);
      sol[i] = atof(val);
    }
  }
  else
  {
    for (i = 0 ; i < (n + m + m) ; ++i) sol[i] = 0.0;
  }



  fclose(f);

  foutput = fopen(basename, "w");

  info = mixedLinearComplementarity_printInFile(problem, foutput);

  if (withSol)
  {
    for (i = 1 ; i < n + m + m ; ++i)
    {
      fprintf(foutput , "%lf " , sol[i]);
    }
  }

  free(sol);

  fclose(foutput);

  freeMixedLinearComplementarityProblem(problem);
  free(basename);
  printf("\n End of test \n");
  return info;


}


int main(int argc, char *argv[])
{
  int info = 0;
  printf("argc %i\n", argc);
  if (argc == 1)
  {
    char filename[50] = "./data/RLCD_mlcp.dat";
    info = write_newformat(filename);
  }
  else if (argc == 2)
  {
    // We assume argv[1] is a filename to open
    FILE *file = fopen(argv[1], "r");

    /* fopen returns 0, the NULL pointer, on failure */
    if (file == 0)
    {
      printf("Could not open file\n");
      fclose(file);
    }
    else
    {
      fclose(file);
      info = write_newformat(argv[1]);

    }

  }
  else
  {
    printf("usage: %s filename", argv[0]);
  }



  return info;
}



