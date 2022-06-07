/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#include "SiconosSharedLibrary.hpp"
#include "SiconosException.hpp"
#ifndef _WIN32
#include <dlfcn.h>                      // for dlerror, dlclose, dlopen, etc
#endif

#include <map>
#include <stddef.h>                     // for nullptr
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <utility>                      // for make_pair, pair
#include <cassert>

namespace SiconosSharedLibrary
{

std::multimap<const std::string, PluginHandle> openedPlugins;
typedef std::multimap<const std::string, PluginHandle>::iterator iter;

PluginHandle loadPlugin(const std::string& pluginPath)
{
  PluginHandle HandleRes;
#ifdef _WIN32
  HandleRes = LoadLibrary(pluginPath.c_str());
  if(!HandleRes)
  {
    DWORD err = GetLastError();
    std::cout << "SiconosSharedLibrary::loadPlugin Error returned : " << err << std::endl;
    std::cout << "Arguments: pluginPath = " << pluginPath << std::endl;
    THROW_EXCEPTION("can not open or find plugin");
  }
#endif
#ifdef _SYS_UNX
  /* -------------------------------------------------------------------------------------- *
   * For RTLD_DEEPBIND, see                                                                 *
   * https://stackoverflow.com/questions/34073051/when-we-are-supposed-to-use-rtld-deepbind *
   * We may want to change this behaviour                                                   *
   * -------------------------------------------------------------------------------------- */

#ifdef __APPLE__
  HandleRes = dlopen(pluginPath.c_str(), RTLD_LAZY);
#else
  HandleRes = dlopen(pluginPath.c_str(), RTLD_LAZY | RTLD_DEEPBIND);
#endif

  if(!HandleRes)
  {
    std::cout << "dlerror() :" << dlerror() <<std::endl;
    THROW_EXCEPTION("can not open or find plugin");
  }
#endif
  openedPlugins.insert(std::make_pair(pluginPath, HandleRes));
  return HandleRes;
}

void * getProcAddress(PluginHandle plugin, const std::string& procedure)
{
  void* ptr;
#ifdef _WIN32
  ptr = (void*) GetProcAddress(plugin, procedure.c_str());
  if(!ptr)
  {
    DWORD err = GetLastError();
    std::cout << "SiconosSharedLibrary::getProcAddress Error returned : " << err << std::endl;
    std::cout << "Arguments: procedure = " << procedure << std::endl;
    THROW_EXCEPTION("can not find procedure");
  }
#endif
#ifdef _SYS_UNX
  ptr = dlsym(plugin, procedure.c_str());
  if(!ptr)
  {
    std::cout << "SiconosSharedLibrary::getProcAddress Error returned : " << dlerror() << std::endl;
    std::cout << "Arguments: procedure = " << procedure << std::endl;
    THROW_EXCEPTION("can not find procedure procedure");
  }
#endif
  return ptr;
}

void closePlugin(const std::string& pluginFile)
{
  iter it = openedPlugins.find(pluginFile);
  if(it == openedPlugins.end())
  {
    std::cout << "SiconosSharedLibrary::closePlugin - could not find an opened plugin named " << pluginFile << std::endl;
    std::cout << "Plugins in openedPlugins:" << std::endl;
    for(iter it2 = openedPlugins.begin(); it2 != openedPlugins.end(); ++it2) std::cout <<  it2->first << std::endl;
    THROW_EXCEPTION("could not find an opened plugin with this name");
  }
  PluginHandle plugin = it->second;
  assert(plugin);
#ifdef _WIN32
  FreeLibrary(plugin);
#endif
#ifdef _SYS_UNX
  dlclose(plugin);
#endif

  openedPlugins.erase(it);
}

}
