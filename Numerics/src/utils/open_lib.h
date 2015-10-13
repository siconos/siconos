/* Siconos-Numerics, Copyright INRIA 2005-2014
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

/** \file open_lib.h
 * \brief function to open library and find functions in then (useful for
 * plugin */

#ifndef OPEN_LIB_H
#define OPEN_LIB_H

#include "SiconosConfig.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** open a library and return an handle (casted as void*)
   * \param lib_name name of the library
   */
  void* open_library(const char* lib_name);

  /** get the address of a function in an already opened lib
   * \param plugin handle to the library
   * \param func name of function
   */
  void* get_function_address(void* plugin, const char* func);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
