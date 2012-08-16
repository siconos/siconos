/* Siconos-Kernel, Copyright INRIA 2005-2011.
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

/*! \file SiconosSharedLibrary.hpp
*/

#ifndef SICONOSSHAREDLIBRARY_H
#define SICONOSSHAREDLIBRARY_H

#include <iostream>
#ifndef _WIN32
#include <dlfcn.h>
#endif
#include <vector>

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
 * \author SICONOS Development Team - copyright INRIA
 * \date (creation) 07/21/2006
 *  Matrices can be either block or Simple.
 *  See Derived classes for details.
 */
namespace SiconosSharedLibrary
{
/**  loads a plugin
 * \param pluginPath full plugin path name
 * \exception SiconosSharedLibraryException if plugin fail to open
 * \return PluginHandle : plugin handle
 */
PluginHandle loadPlugin(const std::string& pluginPath);

/** Gets procedure address
 * \param plugin the plugin handle
 * \param procedure the procedure name
 * \exception SiconosSharedLibraryException if procedure not found
 * \return pointer on procedure
 */
void * getProcAddress(PluginHandle plugin, const std::string& procedure);

/**  Closes plugin
 * \param plugin the plugin handle
 */
void closePlugin(PluginHandle plugin);

/** Closes all plugin set using the current object
 */
void closeAllPlugins();

}

/** Alias for SiconosSharedLibrary */
namespace SSL = SiconosSharedLibrary;

#endif //SICONOSSHAREDLIBRARY_H
