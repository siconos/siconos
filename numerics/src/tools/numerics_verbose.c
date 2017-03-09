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
#include <string.h>
#include <stdlib.h>
#include "numerics_verbose.h"
#include "sn_error_handling.h"

#include <stdarg.h>
/* Default value for verbose mode: turned to off
Warning: global variable
*/
tlsvar int verbose = 0;

void setNumericsVerbose(int newVerboseMode)
{
  verbose = newVerboseMode;
}

void numerics_set_verbose(int newVerboseMode)
{
  verbose = newVerboseMode;
}


void numerics_error(const char* functionName, const char* message)
{
  char output[2048] = "Numerics error - ";
  strncat(output, functionName, strlen(functionName) - 1);
  strncat(output, ":\t", 3);
  strncat(output, message, strlen(message) - 1);
  strncat(output, ".\n", 3);
  sn_fatal_error(SN_UNKOWN_ERROR, output);
}

void numerics_error_nonfatal(const char* fn_name, const char* msg)
{
  size_t fn_name_len = strlen(fn_name);
  size_t msg_len = strlen(msg);
  char* output = (char*)malloc(fn_name_len + msg_len + 7);
  strncpy(output, fn_name, fn_name_len);
  strncat(output, ":\t", 3);
  strncat(output, msg, msg_len);
  strncat(output, ".\n", 3);
  fputs(output, stderr);
  free(output);

}

/* the warning on vprintf is reported as a bug of clang ... --vacary */
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wformat-nonliteral"
#endif
void numerics_warning(const char * functionName, char* fmt, ...)
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

void numerics_printf_verbose(int verbose_mode, const char * fmt, ...)
{
  if (verbose_mode)
  {
    va_list args;
    va_start(args,fmt);
    printf("Numerics: ");
    vprintf(fmt,args);
    va_end(args);
  }
}

#ifdef __clang__
#pragma clang diagnostic pop
#endif
