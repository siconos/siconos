/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#include "SiconosSharedLibrary.h"
#include "SiconosSharedLibraryException.h"
using namespace std;

//*************************************************************************
//
// Public methods
//
//*************************************************************************

SiconosSharedLibrary::SiconosSharedLibrary(const unsigned int & n)
{
  if (n > 0)
    isPlugged.reserve(n);
}

SiconosSharedLibrary::~SiconosSharedLibrary()
{
  vector<PluginHandle>::iterator iter;
  for (iter = isPlugged.begin(); iter != isPlugged.end(); ++iter)
    dlclose(*iter);
  isPlugged.clear();
}

///////////////////////////////////////////////////////////////////////////
//
// loadPlugin
//
PluginHandle SiconosSharedLibrary::loadPlugin(const string& pluginPath)
{
  PluginHandle HandleRes;
  //  cout << "PluginPath   :  " << pluginPath << endl;
#ifdef _SYS_WNT
  HandleRes = LoadLibrary(pluginPath.c_str());
#endif
#ifdef _SYS_UNX
  HandleRes = dlopen(pluginPath.c_str(), RTLD_LAZY);
  //cout << "dlerror() :" <<dlerror()<< endl;
  if (HandleRes == NULL)
  {
    cout << "dlerror() :" << dlerror() << endl;
    SiconosSharedLibraryException::selfThrow("SiconosSharedLibrary::loadPlugin, can not open or found " + pluginPath);
  }
#endif
  isPlugged.push_back(HandleRes);
  return HandleRes;
}

///////////////////////////////////////////////////////////////////////////
//
// getProcAddress
//
void* SiconosSharedLibrary::getProcAddress(PluginHandle plugin, const string& procedure)
{
#ifdef _SYS_WNT
  return GetProcAddress(plugin, procedure.c_str());
#endif
#ifdef _SYS_UNX
  void* ptr = dlsym(plugin, procedure.c_str());
  if (ptr == NULL)
    throw SiconosSharedLibraryException(dlerror());
  return ptr;
#endif
}

///////////////////////////////////////////////////////////////////////////
//
// closePlugin
//
void  SiconosSharedLibrary::closePlugin(PluginHandle plugin)
{
#ifdef _SYS_UNX
  dlclose(plugin);
#endif
}

void  SiconosSharedLibrary::closeAllPlugins()
{
#ifdef _SYS_UNX
  vector<PluginHandle>::iterator iter;
  for (iter = isPlugged.begin(); iter != isPlugged.end(); ++iter)
    dlclose(*iter);
#endif
}

///////////////////////////////////////////////////////////////////////////
//
// getSharedLibraryExtension
//
const string SiconosSharedLibrary::getSharedLibraryExtension() const
{
#ifdef _SYS_WNT
  return ".dll";
#endif
#ifdef _SYS_UNX
  return ".so";
#endif
}

///////////////////////////////////////////////////////////////////////////
//
// link a function
//
void SiconosSharedLibrary::setFunction(void* fPtr, const string& pluginPath, const string& fName)
{
  // load the library
  PluginHandle handle = loadPlugin(pluginPath.c_str());
  // get the function pointer
  *(void **)(fPtr) = getProcAddress(handle, fName.c_str());
}

///////////////////////////////////////////////////////////////////////////
//
// new functions for the plugins
//
const string SiconosSharedLibrary::getPluginName(const string& s) const
{
  string res;

  if ((s.find("\n", 0) != string::npos) || (s.find("\t", 0) != string::npos) || (s.find(" ", 0) != string::npos))
  {
    //raise an exception
    throw SiconosSharedLibraryException("% SharedLibrary managment - getPluginName - The 'string' which contains the plugin name contains '\\n' or '\\t' or ' '");
  }
  else if ((s.find(":", 0) == string::npos) && (s.rfind(":", 0) != s.rfind(":", 0)))
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

const string SiconosSharedLibrary::getPluginFunctionName(const string& s) const
{
  string res;

  if ((s.find("\n", 0) != string::npos) || (s.find("\t", 0) != string::npos) || (s.find(" ", 0) != string::npos))
  {
    //raise an exception
    throw SiconosSharedLibraryException("% SharedLibrary managment - getPluginFunctionName - The 'string' which contains the plugin function name contains '\\n' or '\\t' or ' '");
  }
  else if ((s.find(":", 0) == string::npos) && (s.rfind(":", 0) != s.rfind(":", 0)))
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

/* eof --------------------------------------------------------------------*/
