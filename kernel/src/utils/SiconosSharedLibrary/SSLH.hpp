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

#ifndef SSLH_H
#define SSLH_H

#include <string>

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
  const std::string getSharedLibraryExtension(void);

  const std::string getPluginName(const std::string& s);

  const std::string getPluginFunctionName(const std::string& s);

  void setFunction(void* fPtr, const std::string& pluginPath, const std::string& fName);

  void closePlugin(const std::string& pluginPath);
}

#endif
