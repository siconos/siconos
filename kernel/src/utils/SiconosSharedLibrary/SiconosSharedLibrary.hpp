/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

/*! \file SiconosSharedLibrary.hpp
*/

#ifndef SICONOSSHAREDLIBRARY_H
#define SICONOSSHAREDLIBRARY_H

#include <string>

#ifndef _WIN32
#define _SYS_UNX
#endif

#ifdef _WIN32
#include <windows.h>
#define DLEXPORT __declspec(dllexport)
typedef HMODULE PluginHandle;
#endif

#ifdef _SYS_UNX
#define  DLEXPORT
typedef void* PluginHandle;
#endif

/** Plug-in utilities
 *
 */
namespace SiconosSharedLibrary
{
  /** loads a plugin
   * \param pluginPath full path to plugin
   * \return PluginHandle the object to handle the plugin
   */
  PluginHandle loadPlugin(const std::string& pluginPath);

  /** Gets procedure address
   * \param plugin the plugin handle
   * \param procedure the procedure name
   * \return pointer on procedure
   */
  void * getProcAddress(PluginHandle plugin, const std::string& procedure);
  
  /**  Closes plugin
   * \param pluginFile the name of the plugin to close
   */
  void closePlugin(const std::string& pluginFile);
}

/** Alias for SiconosSharedLibrary */
namespace SSL = SiconosSharedLibrary;

#endif //SICONOSSHAREDLIBRARY_H
