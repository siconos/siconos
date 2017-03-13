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
/* the warning on vprintf is reported as a bug of clang ... --vacary */
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wformat-nonliteral"
#endif

tlsvar int verbose = 0;
tlsvar FILE* logger_f = NULL;
tlsvar void* numerics_logger = NULL;
tlsvar enum numerics_loggers numerics_logger_type = NUMERICS_LOG_TO_SCREEN;

static void numerics_printf_internal(int level, const char* fmt, const char* extra_qual, va_list argp)
{
  switch (numerics_logger_type)
  {
  case NUMERICS_EXTERNAL_LOGGER:
  {
//    if (!numerics_logger)
//    {
      numerics_error("numerics_printf_internal", "unsupported custom logger");
//    }
    break;
  }
  case NUMERICS_LOG_TO_FILE:
  {
    if (!logger_f)
    {
      numerics_error("numerics_printf_internal", "no logger file opened!");
    }
    fputs("[Numerics]", logger_f);
//    fprintf(logger_f, "[%d]", level);
    if (extra_qual) fputs(extra_qual, logger_f);
    fputs(" ", logger_f);
    vfprintf(logger_f, fmt, argp);
    fputs("\n", logger_f);
    break;
  }
  case NUMERICS_LOG_TO_SCREEN:
  default:
  {
    puts("[Numerics]");
//    printf("[%d]", level);
    if (extra_qual) puts(extra_qual);
    puts(" ");
    vprintf(fmt, argp);
    puts("\n");
  }
  }
}

void numerics_set_verbose(int newVerboseMode)
{
  verbose = newVerboseMode;
}


void numerics_error(const char* fn_name, const char* msg, ...)
{
  char output[2048] = "[Numerics][fatal error] ";
  size_t fn_name_len = strlen(fn_name);

  strncat(output, fn_name, fn_name_len - 1);
  strncat(output, "::\t", 3);
  size_t cur_len = strlen(output);
  va_list args;
  va_start(args, msg);
  vsnprintf(&output[cur_len], 2048 - cur_len - 4, msg, args);
  va_end(args);
  strncat(output, ".\n", 3);
  sn_fatal_error(SN_UNKOWN_ERROR, output);
}

void numerics_error_nonfatal(const char* fn_name, const char* msg, ...)
{
  size_t fn_name_len = strlen(fn_name);
  const char* start_str = "[non-fatal error]. %s :: ";
  size_t total_len = fn_name_len + strlen(start_str);
  char* output = (char*)malloc(fn_name_len + strlen(start_str));
  snprintf(output, total_len, "[non-fatal error]. %s :: ", fn_name);

  va_list args;
  va_start(args, msg);

  numerics_printf_internal(0, msg, output, args);

  va_end(args);
  free(output);

}

void numerics_warning(const char * fn_name, char* msg, ...)
{
  size_t fn_name_len = strlen(fn_name);
  const char* start_str = "[warning]. %s :: ";
  size_t total_len = fn_name_len + strlen(start_str);
  char* output = (char*)malloc(fn_name_len + strlen(start_str));
  snprintf(output, total_len, "[warning]. %s :: ", fn_name);

  va_list args;
  va_start(args, msg);

  numerics_printf_internal(0, msg, output, args);

  va_end(args);
  free(output);

}


void numerics_printf(const char * fmt, ...)
{
  if (verbose)
  {
    va_list args;
    va_start(args,fmt);
    numerics_printf_internal(0, fmt, NULL, args);
    va_end(args);
  }
}

void numerics_printf_verbose(int verbose_mode, const char * fmt, ...)
{
  if (verbose_mode)
  {
    va_list args;
    va_start(args,fmt);
    numerics_printf_internal(0, fmt, NULL, args);
    va_end(args);
  }
}

#ifdef __clang__
#pragma clang diagnostic pop
#endif
