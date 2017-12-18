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
#include "NonSmoothDrivers.h"
#include "fclib_interface.h"

static int write_test_fclib(char * filename)
{
  printf("\n Start of test \n");
  printf("Test on %s\n", filename);
  int info = 0;
  size_t sizeoffilename = strlen(filename);
  printf("sizeoffilename %ld\n",  sizeoffilename);
  char  extension[5]; // Don't forget the NUL byte ... --xhub
  strncpy(extension, &filename[sizeoffilename - sizeof(extension) + 1], sizeof(extension));
  char * basename;

  if (strcmp(extension, ".dat") == 0)
  {
    basename = (char *)malloc((sizeoffilename + 2) * sizeof(char *));
    strcpy(basename, filename);
    strncpy(&basename[sizeoffilename - 5], ".hdf5", 6);
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

  FrictionContactProblem* problem = (FrictionContactProblem *)malloc(sizeof(FrictionContactProblem));

  info = frictionContact_newFromFile(problem, f);
  fclose(f);

  int n = 100;
  char * title = (char *)malloc(n * sizeof(char *));
  strncpy(title, "Confeti-ex03-Fc3D-SBM", n);
  char * description = (char *)malloc(n * sizeof(char *));

  strncpy(description, "Rewriting Siconos Numerics test ", n);
  strncat(description, filename, n - strlen(description) - 1);
  strncat(description, " in FCLIB format", n - strlen(description) - 1);
  char * mathInfo = (char *)malloc(n * sizeof(char *));
  strncpy(mathInfo,  "unknown", n);

  frictionContact_fclib_write(problem,
                              title,
                              description,
                              mathInfo,
                              basename,0);

  /* read fclib problem */
  FrictionContactProblem* problem1 = frictionContact_fclib_read(basename);

  info += abs(problem->numberOfContacts - problem1->numberOfContacts);
  info += abs(problem->dimension - problem1->dimension);

  /*...*/

  /* an attempt : we write the problem in the file and read it again
     this does not work because SBM_from_csparse eliminates some zeros
  FrictionContactProblem* problem_from_file = frictionContact_fclib_read(basename);

  frictionContact_fclib_write(problem1,
                              title,
                              description,
                              mathInfo,
                              basename);


  SBM_print(problem1->M->matrix1);

  SBM_print(problem_from_file->M->matrix1);

  info += !(problem_from_file->M->matrix1->filled1 == problem1->M->matrix1->filled1);
  info += !(problem_from_file->M->matrix1->filled2 == problem1->M->matrix1->filled2);

  for(size_t i = 0; i<problem_from_file->M->matrix1->filled1; i++)
  {
    info += !(problem_from_file->M->matrix1->index1_data[i] == problem1->M->matrix1->index1_data[i]);
  }

  for(size_t i = 0; i<problem_from_file->M->matrix1->filled2; i++)
  {
    info += !(problem_from_file->M->matrix1->index2_data[i] == problem1->M->matrix1->index2_data[i]);
    }*/

  freeFrictionContactProblem(problem);
  freeFrictionContactProblem(problem1);
  free(basename);
  free(title);
  free(description);
  free(mathInfo);
  printf("\n End of test \n");
  return info;
}

int main(int argc, char *argv[])
{
  int info=-1;
  printf("argc %i\n", argc);
  if (argc == 1)
  {
    info = write_test_fclib("./data/Confeti-ex03-Fc3D-SBM.dat");
  }
  else if (argc == 2)
  {
    // We assume argv[1] is a filename to open
    FILE *file = fopen(argv[1], "r");

    /* fopen returns 0, the NULL pointer, on failure */
    if (file == 0)
    {
      printf("Could not open file\n");
    }
    else
    {
      info = write_test_fclib(argv[1]);
    }

  }
  else
  {
    printf("usage: %s filename", argv[0]);
  }

  return info;
}



