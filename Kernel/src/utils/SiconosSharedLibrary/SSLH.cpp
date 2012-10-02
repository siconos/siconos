/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
    res = res + getSharedLibraryExtension();
    return res;
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

}
