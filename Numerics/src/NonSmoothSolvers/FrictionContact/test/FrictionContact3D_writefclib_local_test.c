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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NonSmoothDrivers.h"
#include "fclib_interface.h"


int write_test_fclib(char * filename)
{
  printf("\n Start of test \n");
  printf("Test on %s\n", filename);
  int info = 0;
  int sizeoffilename = strlen(filename);
  printf("sizeoffilename %d\n",  sizeoffilename);
  char  extension[4];
  strncpy(extension, &filename[sizeoffilename - 4], 4);
  char * basename;

  if (strcmp(extension, ".dat") == 0)
  {
    basename = (char *)malloc((sizeoffilename + 1) * sizeof(char *));
    strcpy(basename, filename);
    strncpy(&basename[sizeoffilename - 4], ".hdf5", 5);
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
  strcpy(title, "Confeti-ex03-Fc3D-SBM");
  char * description = (char *)malloc(n * sizeof(char *));

  strcat(description, "Rewriting Siconos Numerics test ");
  strcat(description, filename);
  strcat(description, " in FCLIB format");
  char * math_info = (char *)malloc(n * sizeof(char *));
  strcpy(math_info,  "unknown");

  frictionContact_fclib_write(problem,
                              title,
                              description,
                              math_info,
                              basename);

  /* read fclib problem */
  FrictionContactProblem* problem1 = frictionContact_fclib_read(basename);

  info += abs(problem->numberOfContacts - problem1->numberOfContacts);
  info += abs(problem->dimension - problem1->dimension);

  /*...*/

  /* an attempt : we write the problem in the file and read it again
     this does not work because sparseToSBM eliminates some zeros
  FrictionContactProblem* problem_from_file = frictionContact_fclib_read(basename);

  frictionContact_fclib_write(problem1,
                              title,
                              description,
                              math_info,
                              basename);


  printSBM(problem1->M->matrix1);

  printSBM(problem_from_file->M->matrix1);

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
  free(basename);
  free(title);
  free(description);
  free(math_info);
  printf("\n End of test \n");
  return info;
}

int main(int argc, char *argv[])
{
  int info;
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



