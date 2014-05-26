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



#include <stdio.h>
#include <stdlib.h>

#ifndef _WIN32
#include <dlfcn.h>
#endif

#ifdef _WIN32
#include <windows.h>
#define DLEXPORT __declspec(dllexport)
typedef HMODULE PluginHandle;
#else
#define  DLEXPORT
typedef void* PluginHandle;
#endif

void* open_library(const char* lib_name)
{
  void* HandleRes;
#ifdef _WIN32
  HandleRes = (void*) LoadLibrary(lib_name);
  if (!HandleRes)
  {
    int err = (int)GetLastError();
    printf("dlopen error number %d while trying to open %s\n", err, lib_name);
    exit(err);
  }
#else
  HandleRes = dlopen(lib_name, RTLD_LAZY);
  if (!HandleRes)
  {
    printf("dlopen error %s while trying to open %s\n", dlerror(), lib_name);
  }
#endif
  return HandleRes;
}

void* get_function_address(void* plugin, const char* func)
{
  void* ptr;
#ifdef _WIN32
  HMODULE pluginW = (HMODULE) plugin;
  ptr = (void*) GetProcAddress(pluginW, func);
  if (NULL == ptr)
  {
    DWORD err = GetLastError();
    printf("Error %d while trying to find procedure %s\n", func);
    exit(1);
  }
#else
  ptr = dlsym(plugin, func);
  if (ptr == NULL)
  {
    printf("Error %s while trying to find procedure %s\n", dlerror(), func);
    exit(EXIT_FAILURE);
  }
#endif
  return ptr;
}

