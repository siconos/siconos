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

#include "PluggedObject.hpp"

#include <assert.h>

#include "SiconosSharedLibrary.hpp"

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

siconos::plugins::PluggedObject::~PluggedObject()
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
