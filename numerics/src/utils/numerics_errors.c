/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

/*!\file numerics_errors.c
 * \brief Implementation of error logging functions
 */

#include "numerics_errors.h"
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

/* Static buffer for formatted messages */
static char msg_buffer[1024];

/* Current log level - can be set at runtime */
static int log_level = 2;  /* 0=errors only, 1=warnings, 2=info, 3=debug */

void numerics_log_error(const char* file, int line, const char* func,
                        const char* error_type, const char* fmt, ...) {
    /* Extract just the filename from the path */
    const char* filename = strrchr(file, '/');
    if (!filename) filename = strrchr(file, '\\');
    filename = filename ? filename + 1 : file;
    
    /* Format the message */
    va_list args;
    va_start(args, fmt);
    vsnprintf(msg_buffer, sizeof(msg_buffer), fmt, args);
    va_end(args);
    
    /* Print to stderr with location info */
    fprintf(stderr, "[ERROR] %s:%d (%s) - %s: %s\n",
            filename, line, func, error_type, msg_buffer);
    fflush(stderr);
}

void numerics_log_warning(const char* file, int line, const char* func,
                          const char* fmt, ...) {
    if (log_level < 1) return;
    
    const char* filename = strrchr(file, '/');
    if (!filename) filename = strrchr(file, '\\');
    filename = filename ? filename + 1 : file;
    
    va_list args;
    va_start(args, fmt);
    vsnprintf(msg_buffer, sizeof(msg_buffer), fmt, args);
    va_end(args);
    
    fprintf(stderr, "[WARN] %s:%d (%s) - %s\n",
            filename, line, func, msg_buffer);
}

void numerics_log_info(const char* file, int line, const char* func,
                       const char* fmt, ...) {
    if (log_level < 2) return;
    
    const char* filename = strrchr(file, '/');
    if (!filename) filename = strrchr(file, '\\');
    filename = filename ? filename + 1 : file;
    
    va_list args;
    va_start(args, fmt);
    vsnprintf(msg_buffer, sizeof(msg_buffer), fmt, args);
    va_end(args);
    
    fprintf(stdout, "[INFO] %s:%d (%s) - %s\n",
            filename, line, func, msg_buffer);
}
