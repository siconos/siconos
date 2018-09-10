/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

#include "open_lib.h"

void* open_library(const char* lib_name, const int flags)
{
  void* HandleRes;
#ifdef _WIN32
  HandleRes = (void*) LoadLibrary(lib_name);
  if (!HandleRes)
  {
    int err = (int)GetLastError();
    printf("LoadLibrary error number %d while trying to open %s\n", err, lib_name);
    exit(err);
  }
#else
  /* -------------------------------------------------------------------------------------- *
   * For RTLD_DEEPBIND, see                                                                 *
   * https://stackoverflow.com/questions/34073051/when-we-are-supposed-to-use-rtld-deepbind *
   * We may want to change this behaviour                                                   *
   * -------------------------------------------------------------------------------------- */

#ifdef __APPLE__
  HandleRes = dlopen(lib_name, RTLD_LAZY | flags);
#else
  HandleRes = dlopen(lib_name, RTLD_LAZY | RTLD_DEEPBIND | flags);
#endif

  if (!HandleRes)
  {
    printf("dlopen error ``%s'' while trying to open library ``%s''\n", dlerror(), lib_name);
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
  if (!ptr)
  {
    DWORD err = GetLastError();
    printf("Error %d while trying to find procedure %s\n", func);
    exit(1);
  }
#else
  ptr = dlsym(plugin, func);
  if (!ptr)
  {
    printf("Error ``%s'' while trying to find procedure %s\n", dlerror(), func);
    exit(EXIT_FAILURE);
  }
#endif
  return ptr;
}

