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
#ifndef _NUMERICS_MISC_H_
#define _NUMERICS_MISC_H_

#include "NumericsConfig.h"
#include <errno.h>

#if defined(__cplusplus) && !defined (BUILD_AS_CPP)
extern "C"
{
#endif

/** print a dense matrix or vector
 */
void printm(unsigned int nl, unsigned int nc, double *m);

/** Check return code of an expression. */
#define CHECK_RETURN(EXPR)                                              \
  do                                                                    \
  {                                                                     \
    if (!EXPR)                                                          \
    {                                                                   \
      fprintf (stderr, "Siconos Numerics: Warning %s failed, %s:%d\n",  \
               #EXPR, __FILE__, __LINE__);                              \
    }                                                                   \
  } while (0)

/** check IO
 */
#define CHECK_IO(EXPR, ...)                                             \
  do                                                                    \
  {                                                                     \
    if (!EXPR)                                                          \
    {                                                                   \
      int * _arr_[] = {__VA_ARGS__};                                    \
      if (errno != 0)                                                   \
      {                                                                 \
        perror(#EXPR);                                                  \
        fprintf (stderr, "Siconos Numerics: Warning %s failed, %s:%d\n", #EXPR, __FILE__, __LINE__); \
        if (sizeof(_arr_)) { *_arr_[0] = errno; }                       \
      }                                                                 \
      else                                                              \
      {                                                                 \
        fprintf (stderr, "Siconos Numerics unknown error for %s, %s:%d\n", #EXPR, __FILE__, __LINE__); \
        if (sizeof(_arr_)) { *_arr_[0] = 1; }                           \
      }                                                                 \
    }                                                                   \
  } while (0)

/** ignore IO. Think carefully before using this ...
 */
#define IGNORE_IO(EXPR)                                                 \
  do                                                                    \
  {                                                                     \
    if(EXPR){};                                                         \
  } while (0)

#ifdef HAVE_MPI
#define CHECK_MPI(EXPR)                                                 \
  do                                                                    \
  {                                                                     \
    int error_code = EXPR;                                                  \
    MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);              \
    if (error_code != MPI_SUCCESS) {                                    \
      char error_string[1024];                                          \
      int length_of_error_string, error_class;                          \
      MPI_Error_class(error_code, &error_class);                        \
      MPI_Error_string(error_class, error_string, &length_of_error_string); \
      fprintf(stderr, "%3d: %s\n", 0, error_string);                    \
      MPI_Error_string(error_code, error_string, &length_of_error_string); \
      fprintf(stderr, "%3d: %s\n", 0, error_string);                    \
      MPI_Abort(MPI_COMM_WORLD, error_code);                            \
    };                                                                  \
  } while(0)
#endif

#ifdef __clang_analyzer__
#define NO_RETURN  __attribute__((analyzer_noreturn))
#else
#define NO_RETURN
#endif

#ifdef __GNUC__
#define MAYBE_UNUSED __attribute__((unused))
#else
#define MAYBE_UNUSED
#endif

#if defined(__cplusplus) && !defined (BUILD_AS_CPP)
}
#endif

#endif
