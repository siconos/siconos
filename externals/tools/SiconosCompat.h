/* Siconos-Numerics, Copyright INRIA 2005-2015
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/
#ifndef _SICONOS_COMPAT
#define _SICONOS_COMPAT

#if defined(_MSC_VER) || defined(__MINGW32__) // check other mingw
    #define SN_SIZE_T_F    "%Iu"
    #define SN_SSIZE_T_F   "%Id"
    #define SN_PTRDIFF_T_F "%Id"

// Note FP : I think that this is obsolete, isn't it ?
// See https://docs.microsoft.com/en-us/cpp/c-runtime-library/format-specification-syntax-printf-and-wprintf-functions?view=vs-2019
// Could we remove this and the define ?

#else
    #define SN_SIZE_T_F    "%zu"
    #define SN_SSIZE_T_F   "%zd"
    #define SN_PTRDIFF_T_F "%zd"
#endif


#ifdef _WIN32
  #define DLLPRE ""
  #define DLLEXT ".dll"
#else
  #define DLLPRE "lib"
  #define DLLEXT ".so"
#endif

#define DLL_FROM_NAME(X) DLLPRE X  DLLEXT

#if defined(_MSC_VER)
// for M_PI
#define _USE_MATH_DEFINES

#endif /* defined(_MSC_VER) */

#if defined(__SUNPRO_CC)
#define INFINITY (DBL_MAX+DBL_MAX)
#define NAN (INFINITY-INFINITY)
#endif

#endif

