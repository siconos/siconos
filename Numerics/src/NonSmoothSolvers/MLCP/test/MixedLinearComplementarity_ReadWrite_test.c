/* Siconos-Numerics, Copyright INRIA 2005-2012.
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
#include "MixedLinearComplementarityProblem.h"

int write_newformat(char *filename)
{
  printf("\n Start of test \n");
  printf("Test on %s\n", filename);
  int info = 0;
  int sizeoffilename = strlen(filename);
  printf("sizeoffilename %d\n",  sizeoffilename);
  char  extension[4];
  strncpy(extension, &filename[sizeoffilename - 4], 4);
  printf("extension %s\n",  extension);
  char * basename;

  if (strcmp(extension, ".dat") == 0)
  {
    basename = (char *)malloc((sizeoffilename + 1) * sizeof(char *));
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
  fscanf(f , "%s" , val);
  int withSol = 0;
  int n = problem->n;
  int m = problem ->m;
  int i;
  double *sol  = (double*)malloc((n + m + m) * sizeof(double));

  if (!feof(f))
  {
    withSol = 1;


    sol[0] = atof(val);

    for (i = 1 ; i < n + m + m ; ++i)
    {
      fscanf(f , "%s" , val);
      sol[i] = atof(val);
    }
  }
  else
  {
    for (i = 0 ; i < (n + m + m) ; ++i) sol[i] = 0.0;
  }

  free(sol);


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


  fclose(foutput);

  freeMixedLinearComplementarityProblem(problem);
  free(basename);
  printf("\n End of test \n");
  return info;


}


int main(int argc, char *argv[])
{
  int info;
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



