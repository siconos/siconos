// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

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
#define DLEXPORT
typedef void* PluginHandle;
#endif

#include "numerics_errors.h"
#include "open_lib.h"

void* open_library(const char* lib_name, const int flags) {
  void* HandleRes;
#ifdef _WIN32
  HandleRes = (void*)LoadLibrary(lib_name);
  if (!HandleRes) {
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

  if (!HandleRes) {
    printf("dlopen error ``%s'' while trying to open library ``%s''\n", dlerror(), lib_name);
  }
#endif
  return HandleRes;
}

void* get_function_address(void* plugin, const char* func) {
  void* ptr;
#ifdef _WIN32
  HMODULE pluginW = (HMODULE)plugin;
  ptr = (void*)GetProcAddress(pluginW, func);
  if (!ptr) {
    DWORD err = GetLastError();
    printf("Error %d while trying to find procedure %s\n", func);
    exit(1);
  }
#else
  ptr = dlsym(plugin, func);
  if (!ptr) {
    printf("Error ``%s'' while trying to find procedure %s\n", dlerror(), func);
    exit(EXIT_FAILURE);
  }
#endif
  return ptr;
}
