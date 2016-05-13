/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

#include "SiconosSharedLibraryException.hpp"
#include "SSLH.hpp"
#include "SiconosSharedLibrary.hpp"

namespace SSLH
{
///////////////////////////////////////////////////////////////////////////
//
// new functions for the plugins
//

///////////////////////////////////////////////////////////////////////////
//
// getSharedLibraryExtension
//
const std::string getSharedLibraryExtension(void)
{
#ifdef _WIN32
  return ".dll";
#else
  return ".so";
#endif
}

const std::string getPluginName(const std::string& s)
{
  std::string res;

  if ((s.find("\n", 0) != std::string::npos) || (s.find("\t", 0) != std::string::npos) || (s.find(" ", 0) != std::string::npos))
  {
    //raise an exception
    throw SiconosSharedLibraryException("% SharedLibrary managment - getPluginName - The 'string' which contains the plugin name contains '\\n' or '\\t' or ' '");
  }
  else if ((s.find(":", 0) == std::string::npos) && (s.rfind(":", 0) != s.rfind(":", 0)))
  {
    //raise an exception
    throw SiconosSharedLibraryException("% SharedLibrary managment - getPluginName - The 'string' which contains the plugin name is not well formed. It must be like : plugin_name:plugin_function_name");
  }
  else
  {
    // return the plugin name
    int pos = s.find(":", 0);
    res = s.substr(0, pos);
    return res + getSharedLibraryExtension();
  }
}

const std::string getPluginFunctionName(const std::string& s)
{
  std::string res;

  if ((s.find("\n", 0) != std::string::npos) || (s.find("\t", 0) != std::string::npos) || (s.find(" ", 0) != std::string::npos))
  {
    //raise an exception
    throw SiconosSharedLibraryException("% SharedLibrary managment - getPluginFunctionName - The 'string' which contains the plugin function name contains '\\n' or '\\t' or ' '");
  }
  else if ((s.find(":", 0) == std::string::npos) && (s.rfind(":", 0) != s.rfind(":", 0)))
  {
    //raise an exception
    throw SiconosSharedLibraryException("% SharedLibrary managment - getPluginFunctionName - The 'string' which contains the plugin name is not well formed. It must be like : plugin_name:plugin_function_name");
  }
  else
  {
    // return the plugin function name
    int pos = s.find(":", 0);
    res = s.substr(pos + 1, s.length());
    return res;
  }
}

void setFunction(void* fPtr, const std::string& pluginPath, const std::string& fName)
{
  // load the library
  PluginHandle handle = SiconosSharedLibrary::loadPlugin(pluginPath.c_str());
  // get the function pointer
  *(void **)(fPtr) = SiconosSharedLibrary::getProcAddress(handle, fName.c_str());
}

void closePlugin(const std::string& pluginPath)
{
  SiconosSharedLibrary::closePlugin(getPluginName(pluginPath));
}

}
