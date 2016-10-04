#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "numerics_verbose.h"

#include <stdarg.h>
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

void numerics_warning(char * functionName, char* fmt, ...)
{
  char output[200] = "Numerics warning - ";
  strcat(output, functionName);
  /* strcat(output, message); */
  /* strcat(output, ".\n"); */
  /* fprintf(stderr, "%s", output); */
  va_list args;
  va_start(args, fmt);
  fputs(output, stderr);
  vfprintf(stderr, fmt, args);
  va_end(args);



}

/* the warning on vprintf is reported as a bug of clang ... --vacary */
#pragma clang diagnostic ignored "-Wformat-nonliteral"
void numerics_printf(const char * fmt, ...)
{
  if (verbose)
  {
    va_list args;
    va_start(args,fmt);
    printf("Numerics: ");
    vprintf(fmt,args);
    va_end(args);
  }
}

void numerics_printf_verbose(int verbose, const char * fmt, ...)
{
  if (verbose)
  {
    va_list args;
    va_start(args,fmt);
    printf("Numerics: ");
    vprintf(fmt,args);
    va_end(args);
  }
}
