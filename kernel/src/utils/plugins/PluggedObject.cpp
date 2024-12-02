/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include "PluggedObject.hpp"

#include <assert.h>

#include "SiconosException.hpp"
#ifndef _WIN32
#include <dlfcn.h>  // for dlerror, dlclose, dlopen, etc
#endif

#include <stddef.h>  // for nullptr

#include <cassert>
#include <iostream>  // for operator<<, basic_ostream, etc
#include <map>
#include <utility>  // for make_pair, pair

std::multimap<const std::string, PluginHandle> openedPlugins;

PluginHandle siconos::plugins::internal::loadPlugin(const std::string& pluginPath)
{
  PluginHandle HandleRes;
#ifdef _WIN32
  HandleRes = LoadLibrary(pluginPath.c_str());
  if (!HandleRes) {
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

  if (!HandleRes) {
    std::cout << "dlerror() :" << dlerror() << std::endl;
    THROW_EXCEPTION("can not open or find plugin");
  }
#endif
  openedPlugins.insert(std::make_pair(pluginPath, HandleRes));
  return HandleRes;
}

void* siconos::plugins::internal::getProcAddress(PluginHandle plugin,
                                                 const std::string& procedure)
{
  void* ptr;
#ifdef _WIN32
  ptr = (void*)GetProcAddress(plugin, procedure.c_str());
  if (!ptr) {
    DWORD err = GetLastError();
    std::cout << "SiconosSharedLibrary::getProcAddress Error returned : " << err << std::endl;
    std::cout << "Arguments: procedure = " << procedure << std::endl;
    THROW_EXCEPTION("can not find procedure");
  }
#endif
#ifdef _SYS_UNX
  ptr = dlsym(plugin, procedure.c_str());
  if (!ptr) {
    std::cout << "SiconosSharedLibrary::getProcAddress Error returned : " << dlerror()
              << std::endl;
    std::cout << "Arguments: procedure = " << procedure << std::endl;
    THROW_EXCEPTION("can not find procedure procedure");
  }
#endif
  return ptr;
}

void siconos::plugins::internal::closePlugin(const std::string& pluginFile)
{
  auto it = openedPlugins.find(pluginFile);
  if (it == openedPlugins.end()) {
    std::cout << "SiconosSharedLibrary::closePlugin - could not find an opened plugin named "
              << pluginFile << std::endl;
    std::cout << "Plugins in openedPlugins:" << std::endl;
    for (auto it2 : openedPlugins) std::cout << it2.first << std::endl;
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

const std::string siconos::plugins::getSharedLibraryExtension(void)
{
#ifdef _WIN32
  return ".dll";
#else
  return ".so";
#endif
}

const std::string siconos::plugins::getPluginName(const std::string& s)
{
  std::string res;

  if ((s.find("\n", 0) != std::string::npos) || (s.find("\t", 0) != std::string::npos) ||
      (s.find(" ", 0) != std::string::npos)) {
    THROW_EXCEPTION(
        "The 'string' which contains the plugin name contains '\\n' or '\\t' or ' '");
  }
  else if ((s.find(":", 0) == std::string::npos) && (s.rfind(":", 0) != s.rfind(":", 0))) {
    THROW_EXCEPTION(
        "The 'string' which contains the plugin name is not well formed. It must be like : "
        "plugin_name:plugin_function_name");
  }
  else {
    // return the plugin name
    int pos = s.find(":", 0);
    res = s.substr(0, pos);
    return res + getSharedLibraryExtension();
  }
}

const std::string siconos::plugins::getPluginFunctionName(const std::string& s)
{
  std::string res;

  if ((s.find("\n", 0) != std::string::npos) || (s.find("\t", 0) != std::string::npos) ||
      (s.find(" ", 0) != std::string::npos)) {
    THROW_EXCEPTION(
        "The 'string' which contains the plugin function name contains '\\n' or '\\t' or ' '");
  }
  else if ((s.find(":", 0) == std::string::npos) && (s.rfind(":", 0) != s.rfind(":", 0))) {
    THROW_EXCEPTION(
        "The 'string' which contains the plugin name is not well formed. It must be like : "
        "plugin_name:plugin_function_name");
  }
  else {
    // return the plugin function name
    int pos = s.find(":", 0);
    res = s.substr(pos + 1, s.length());
    return res;
  }
}

void siconos::plugins::setFunction(void* fPtr, const std::string& pluginPath,
                                   const std::string& fName)
{
  // load the library
  PluginHandle handle = siconos::plugins::internal::loadPlugin(pluginPath.c_str());
  // get the function pointer
  *(void**)(fPtr) = siconos::plugins::internal::getProcAddress(handle, fName.c_str());
}

void siconos::plugins::closePlugin(const std::string& pluginPath)
{
  siconos::plugins::internal::closePlugin(getPluginName(pluginPath));
}

siconos::plugins::PluggedObject::PluggedObject(const std::string& name) : _pluginName(name)
{
  setComputeFunction();
}

siconos::plugins::PluggedObject::PluggedObject(const PluggedObject& PO)
    : _pluginName(PO.pluginName())
{
  // we don't copy the fPtr since we need to increment the number of times we opened the plugin
  // file in the openedPlugins multimap
  if ((_pluginName.compare("unplugged") != 0) && (_pluginName.compare("Unknown") != 0))
    setComputeFunction();
}

siconos::plugins::PluggedObject::~PluggedObject() noexcept
{
  if ((_pluginName.compare("unplugged") != 0) && (_pluginName.compare("Unknown") != 0))
    siconos::plugins::closePlugin(_pluginName);
}

void siconos::plugins::PluggedObject::setComputeFunction(const std::string& pluginPath,
                                                         const std::string& functionName)
{
  std::string ext = siconos::plugins::getSharedLibraryExtension();
  if (ext.compare(pluginPath.substr(pluginPath.size() - ext.size())) == 0) {
    siconos::plugins::setFunction(&fPtr, pluginPath, functionName);
    _pluginName = pluginPath.substr(0, pluginPath.find_last_of(".")) + ":" + functionName;
  }
  else {
    siconos::plugins::setFunction(&fPtr, pluginPath + ext, functionName);
    _pluginName = pluginPath + ":" + functionName;
  }
}

void siconos::plugins::PluggedObject::setComputeFunction(const std::string& plugin)
{
  siconos::plugins::setFunction(&fPtr, siconos::plugins::getPluginName(plugin),
                                siconos::plugins::getPluginFunctionName(plugin));
  _pluginName = plugin;
}

void siconos::plugins::PluggedObject::setComputeFunction(void)
{
  assert(_pluginName != "unplugged" &&
         "PluggedObject::setComputeFunction error, try to plug an unnamed function.");
  siconos::plugins::setFunction(&fPtr, siconos::plugins::getPluginName(_pluginName),
                                siconos::plugins::getPluginFunctionName(_pluginName));
}
