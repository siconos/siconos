#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "misc.h"

/* Default value for verbose mode: turned to off
Warning: global variable
*/
int verbose = 0;

void setNumericsVerbose(int newVerboseMode)
{
  verbose = newVerboseMode;
}

void numericsError(char * functionName, char* message)
{
  char output[200] = "Numerics error - ";
  strcat(output, functionName);
  strcat(output, message);
  strcat(output, ".\n");
  fprintf(stderr, "%s", output);
  exit(EXIT_FAILURE);
}

void numericsWarning(char * functionName, char* message)
{
  char output[200] = "Numerics warning - ";
  strcat(output, functionName);
  strcat(output, message);
  strcat(output, ".\n");
  fprintf(stderr, "%s", output);
  exit(EXIT_FAILURE);
}

void printm(unsigned int nl, unsigned int nc, double *m)
{
  for (unsigned int i = 0; i < nl; ++i)
  {
    for (unsigned int j = 0; j < nc; ++j)
    {
      printf("%32.24e ", *(m + j + nc * i));
    }
    printf("\n");
  }
}
