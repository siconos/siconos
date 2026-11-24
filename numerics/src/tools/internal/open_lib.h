// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/** \file open_lib.h
 * \brief function to open library and find functions in then (useful for
 * plugin */

#ifndef OPEN_LIB_H
#define OPEN_LIB_H

#if defined(__cplusplus)
extern "C" {
#endif

/** open a library and return an handle (casted as void*)
 * \param lib_name name of the library
 * \param flags additional flags (for dlopen)
 */
void* open_library(const char* lib_name, const int flags);

/** get the address of a function in an already opened lib
 * \param plugin handle to the library
 * \param func name of function
 */
void* get_function_address(void* plugin, const char* func);

#if defined(__cplusplus)
}
#endif

#endif
