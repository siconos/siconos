/* Siconos-Kernel, Copyright INRIA 2005-2011.
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

/*! \file debug.h
  \brief Some debug facilities
*/

#ifndef DEBUG_H
#define DEBUG_H

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#ifdef DEBUG_MESSAGES
#ifdef DEBUG_WHERE_MESSAGES
#define DEBUG_WHERESTR  "%s:%d:\n"
#define DEBUG_WHEREARG  __FILE__, __LINE__
#else
#define DEBUG_WHERESTR  "%s"
#define DEBUG_WHEREARG  ""
#endif

#ifdef DEBUG_STDOUT

#ifdef DEBUG_NOCOLOR
#define DEBUG_INTERNAL_PRINTF(...)      printf(__VA_ARGS__);
#else
#define DEBUG_INTERNAL_PRINTF(...)      printf(ANSI_COLOR_RED); printf(__VA_ARGS__); printf(ANSI_COLOR_RESET);
#endif
#else

#ifdef DEBUG_NOCOLOR
#define DEBUG_INTERNAL_PRINTF(...)       fprintf(stderr, __VA_ARGS__);
#else
#define DEBUG_INTERNAL_PRINTF(...)       fprintf(stderr, ANSI_COLOR_RED); fprintf(stderr, __VA_ARGS__); fprintf(stderr, ANSI_COLOR_RESET);
#endif
#endif

#define DEBUG_BEGIN(M) DEBUG_INTERNAL_PRINTF("==== BEGIN === %s",M)
#define DEBUG_END(M)  DEBUG_INTERNAL_PRINTF("===   END  === %s",M)

#ifdef DEBUG_BEGIN_END_ONLY
#define DEBUG_PRINTF(_fmt, ...)
#define DEBUG_PRINT(M)
#define DEBUG_EXPR(E)
#define DEBUG_EXPR_WE(E)
#define DEBUG_GLOBAL_VAR_DECL(D)
#else
#define DEBUG_PRINTF(_fmt, ...)  DEBUG_INTERNAL_PRINTF(DEBUG_WHERESTR _fmt, DEBUG_WHEREARG, __VA_ARGS__)
#define DEBUG_PRINT(M)  DEBUG_PRINTF("%s",M)
#define DEBUG_EXPR(E) DEBUG_PRINTF("%s: ", #E) do { E ; } while(0)
#define DEBUG_EXPR_WE(E) do { E ; } while(0)
#define DEBUG_GLOBAL_VAR_DECL(D) D
#endif

#else
#define DEBUG_BEGIN(M)
#define DEBUG_END(M)
#define DEBUG_PRINTF(_fmt, ...)
#define DEBUG_PRINT(M)
#define DEBUG_EXPR(E)
#define DEBUG_EXPR_WE(E)
#define DEBUG_GLOBAL_VAR_DECL(D)
#endif

#define DEBUG_PRINT_MAT_STR(NAME, M, nrows, ncols) \
DEBUG_PRINT(#NAME " matrix\n"); \
DEBUG_EXPR_WE(for (unsigned i = 0; i < nrows; ++i) \
  { for(unsigned j = 0 ; j < ncols; ++j) \
  { DEBUG_PRINTF(ANSI_COLOR_BLUE " % 2.2e " ANSI_COLOR_RESET, M[i + j*nrows]) } \
   DEBUG_PRINT("\n")});

#define DEBUG_PRINT_MAT(M, nrows, ncols) DEBUG_PRINT_MAT_STR(#M, M, nrows, ncols)

#define DEBUG_PRINT_MAT_SMALL_STR(NAME, M, nrows, ncols, ...) \
DEBUG_PRINT(#NAME " matrix\n"); \
DEBUG_EXPR_WE(double* _arr_[] = {__VA_ARGS__}; for (unsigned i = 0; i < nrows; ++i) \
  { for (unsigned k = 0; k < sizeof(_arr_)/sizeof(double*); ++k) {DEBUG_PRINTF(ANSI_COLOR_YELLOW "% 1.0e " ANSI_COLOR_RESET, _arr_[k][i])} \
    for(unsigned j = 0 ; j < ncols; ++j) \
  { if (fabs(M[i + j*nrows]) > 2.2e-16) {DEBUG_PRINTF(ANSI_COLOR_YELLOW " % 2.f " ANSI_COLOR_RESET, M[i + j*nrows])} \
    else { DEBUG_PRINT(ANSI_COLOR_BLUE " . " ANSI_COLOR_RESET) } } \
   DEBUG_PRINT("\n")});

#define DEBUG_PRINT_SMALL_MAT(M, nrows, ncols) DEBUG_PRINT_MAT_SMALL_STR(#M, M, nrows, ncols)

#define DEBUG_PRINT_MAT_ROW_MAJOR_STR(NAME, M, nrows, ncols) \
DEBUG_PRINT(#NAME " matrix\n"); \
DEBUG_EXPR_WE(for (unsigned i = 0; i < nrows; ++i) \
  { for(unsigned j = 0 ; j < ncols; ++j) \
  { DEBUG_PRINTF(ANSI_COLOR_BLUE "% 2.2e " ANSI_COLOR_RESET, M[i*ncols + j]) } \
   DEBUG_PRINT("\n")});

#define DEBUG_PRINT_MAT_ROW_MAJOR(M, nrows, ncols) DEBUG_PRINT_MAT_ROW_MAJOR_STR(#M, M, nrows, ncols)

#define DEBUG_PRINT_MAT_ROW_MAJOR_NCOLS_STR(NAME, M, nrows, ncols, ncols_to_display) \
DEBUG_PRINT(#NAME " matrix\n"); \
DEBUG_EXPR_WE(for (unsigned i = 0; i < nrows; ++i) \
  { for(unsigned j = 0 ; j < ncols_to_display; ++j) \
  { DEBUG_PRINTF(ANSI_COLOR_BLUE "% 2.2e " ANSI_COLOR_RESET, M[i*ncols + j]) } \
   DEBUG_PRINT("\n")});

#define DEBUG_PRINT_MAT_ROW_MAJOR_NCOLS(M, nrows, ncols, ncols_to_display) DEBUG_PRINT_MAT_ROW_MAJOR_NCOLS_STR(#M, M, nrows, ncols, ncols_to_display)

#define DEBUG_PRINT_MAT_ROW_MAJOR_NCOLS_SMALL_STR(NAME, M, nrows, ncols, ncols_to_display) \
DEBUG_PRINT(#NAME " matrix\n"); \
DEBUG_EXPR_WE(for (unsigned i = 0; i < nrows; ++i) \
  { for(unsigned j = 0 ; j < ncols_to_display; ++j) \
  { if (fabs(M[i*ncols + j]) > 2.2e-16) {DEBUG_PRINTF(ANSI_COLOR_YELLOW "% 1.0e " ANSI_COLOR_RESET, M[i*ncols + j])} \
    else { DEBUG_PRINT(ANSI_COLOR_BLUE " 0     " ANSI_COLOR_RESET) } } \
   DEBUG_PRINT("\n")});

#define DEBUG_PRINT_MAT_ROW_MAJOR_NCOLS_SMALL2_STR(NAME, M, nrows, ncols, ncols_to_display, ...) \
DEBUG_PRINT(#NAME " matrix\n"); \
DEBUG_EXPR_WE( double* _arr_[] = {__VA_ARGS__}; for (unsigned i = 0; i < nrows; ++i) \
  { for (unsigned k = 0; k < sizeof(_arr_)/sizeof(double*); ++k) {DEBUG_PRINTF(ANSI_COLOR_YELLOW "% 1.0e " ANSI_COLOR_RESET, _arr_[k][i])} \
    for(unsigned j = 0 ; j < ncols_to_display; ++j) \
  { if (fabs(M[i*ncols + j]) > 2.2e-16) {DEBUG_PRINTF(ANSI_COLOR_YELLOW "% 1.0e " ANSI_COLOR_RESET, M[i*ncols + j])} \
    else { DEBUG_PRINT(ANSI_COLOR_BLUE " 0 " ANSI_COLOR_RESET) } } \
   DEBUG_PRINT("\n")});

#define DEBUG_PRINT_MAT_ROW_MAJOR_SMALL_NCOLS(M, nrows, ncols, ncols_to_display) DEBUG_PRINT_MAT_ROW_MAJOR_SMALL_NCOLS_STR(#M, M, nrows, ncols, ncols_to_display)

#define DEBUG_PRINT_VEC_STR(MSG, V, size) \
DEBUG_PRINT(MSG " vector\n"); \
DEBUG_EXPR_WE(for (unsigned i = 0; i < size; ++i) \
  { DEBUG_PRINTF(ANSI_COLOR_GREEN "% 2.2e " ANSI_COLOR_RESET, V[i]) }\
   DEBUG_PRINT("\n"));

#define DEBUG_PRINT_VEC(V, size) DEBUG_PRINT_VEC_STR(#V, V, size)

#define DEBUG_PRINT_VEC_INT_STR(MSG, V, size) \
DEBUG_PRINT(MSG " vector\n"); \
DEBUG_EXPR_WE(for (unsigned i = 0; i < size; ++i) \
  { DEBUG_PRINTF(ANSI_COLOR_CYAN "%d " ANSI_COLOR_RESET, V[i]) }\
   DEBUG_PRINT("\n"));

#define DEBUG_PRINT_VEC_INT(V, size) DEBUG_PRINT_VEC_INT_STR(#V, V, size)


#define DEBUG_OR_VERBOSE(X) if (verbose > 0) { X; } else { DEBUG_EXPR_WE(X); }

#endif
