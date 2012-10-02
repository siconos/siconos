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
#include "SiconosSharedLibrary.hpp"
#include "SiconosSharedLibraryException.hpp"

using namespace std;

namespace SiconosSharedLibrary
{

std::vector<PluginHandle> isPlugged;

PluginHandle loadPlugin(const string& pluginPath)
{
  PluginHandle HandleRes;
#ifdef _WIN32
  HandleRes = LoadLibrary(pluginPath.c_str());
  if (!HandleRes)
  {
    DWORD err = GetLastError();
    cout << "Error returned : " << err << endl;
    SiconosSharedLibraryException::selfThrow("SiconosSharedLibrary::loadPlugin, can not open or found " + pluginPath);
  }
#endif
#ifdef _SYS_UNX
  HandleRes = dlopen(pluginPath.c_str(), RTLD_LAZY);
  if (!HandleRes)
  {
    cout << "dlerror() :" << dlerror() << endl;
    SiconosSharedLibraryException::selfThrow("SiconosSharedLibrary::loadPlugin, can not open or found " + pluginPath);
  }
#endif
  isPlugged.push_back(HandleRes);
  return HandleRes;
}

void * getProcAddress(PluginHandle plugin, const string& procedure)
{
#ifdef _WIN32
  return (void*) GetProcAddress(plugin, procedure.c_str());
#endif
#ifdef _SYS_UNX
  void* ptr = dlsym(plugin, procedure.c_str());
  if (!ptr)
    throw SiconosSharedLibraryException(dlerror());
  return ptr;
#endif
}

void closePlugin(PluginHandle plugin)
{
#ifdef _WIN32
  FreeLibrary(plugin);
#endif
#ifdef _SYS_UNX
  dlclose(plugin);
#endif
}

void closeAllPlugins()
{
  vector<PluginHandle>::iterator iter;
  for (iter = isPlugged.begin(); iter != isPlugged.end(); ++iter)
  {
#ifdef _WIN32
    FreeLibrary(*iter);
#endif
#ifdef _SYS_UNX
    dlclose(*iter);
#endif
  }
}


}
