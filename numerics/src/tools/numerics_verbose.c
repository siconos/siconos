#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "numerics_verbose.h"

/* Default value for verbose mode: turned to off
Warning: global variable
*/
int verbose = 0;

void setNumericsVerbose(int newVerboseMode)
{
  verbose = newVerboseMode;
}

void numerics_set_verbose(int newVerboseMode)
{
  verbose = newVerboseMode;
}


void numerics_error(char * functionName, char* message)
{
  char output[300] = "Numerics error - ";
  strcat(output, functionName);
  strcat(output, " :\t");
  strcat(output, message);
  strcat(output, ".\n");
  fprintf(stderr, "%s", output);
  exit(EXIT_FAILURE);
}

void numerics_warning(char * functionName, char* message)
{
  char output[200] = "Numerics warning - ";
  strcat(output, functionName);
  strcat(output, message);
  strcat(output, ".\n");
  fprintf(stderr, "%s", output);
  exit(EXIT_FAILURE);
}
