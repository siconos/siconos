/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2006.
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
#ifndef SICONOSSHAREDLIBRARY_H
#define SICONOSSHAREDLIBRARY_H

#include "SiconosSharedLibraryException.h"
#include <string>
#include <iostream>
#include <dlfcn.h>
#include <vector>

#define _SYS_UNX

#ifdef _SYS_WNT
#include <windows.h>
#define DLEXPORT __declspec(dllexport)
typedef HMODULE PluginHandle;
#endif

#ifdef _SYS_UNX
#define  DLEXPORT
typedef void* PluginHandle;
#endif

// --------------------------------------------------------------------------
class SiconosSharedLibrary
{
private:

  /* to manage plug-in object closing (for each call to dlopen one object is add into isPlugged)*/
  std::vector<PluginHandle> isPlugged;

public:

  /** \fn SiconosSharedLibrary(const unsigned int& n = 0)
   *  \brief constructor
   *    \param unsigned int n, the number of plug-in that will be "connected" through the current object (optional)
   */
  SiconosSharedLibrary(const unsigned int& = 0);

  /** \fn ~SiconosSharedLibrary(){}
   *  \brief destructor
   */
  ~SiconosSharedLibrary();

  /** \fn PluginHandle loadPlugin(const std::string& pluginPath)
   *  \brief load a plugin
   *  \param std::string pluginPath : full plugin path name
   *  \exception SiconosSharedLibraryException if plugin fail to open
   *  \return PluginHandle : plugin handle
   */
  PluginHandle loadPlugin(const std::string& pluginPath) ;

  /** \fn void* getProcAddress(PluginHandle plugin, const std::string& procedure)
   *  \brief get procedure address
   *  \param PluginHandle plugin : plugin handle
   *  \param std::string procedure : procedure name
   *  \exception SiconosSharedLibraryException if procedure not found
   *  \return pointer on procedure
   */
  void* getProcAddress(PluginHandle plugin, const std::string& procedure);

  /** \fn void closePlugin(PluginHandle plugin)
   *  \brief close plugin
   *  \param PluginHandle plugin : plugin handle
   */
  void closePlugin(PluginHandle plugin);

  /** \fn void closeAllPlugins()
   *  \brief close all plugin set using the current object
   */
  void closeAllPlugins();

  /** \fn std::string getSharedLibraryExtension()
   *  \brief get shared library extension
   *  \return library extension ("*.so" for UNIX or "*.dll" for WNT)
   */
  const std::string getSharedLibraryExtension() const ;

  /** \fn void setFunction(void* functionPtr, PluginHandle pluginHandle, const std::string pluginPath&, const std::string functionName&)
   *  \brief set a function pointer to a function in an external library. Don't use it directely
   *  \param functionPtr : pointer to the function (in-out)
   *  \param std::string : the complet path to the external plugin (in)
   *  \param std::string : the name of the function to reference (in)
   */
  void setFunction(void* functionPtr, const std::string& pluginPath, const std::string& functionName);

  /** \fn std::string getPluginName(std::string )
   *  \brief extract the plugin name from a std::string containing data to call a plugin function
   *  \param std::string : its form is : "pluginName:functionName"
   *  \return a std::string containing the plugin name
   */
  const std::string getPluginName(const std::string&) const ;

  /** \fn std::string getPluginFunctionName(std::string )
   *  \brief extract the function name from a std::string containing data to call a plugin function
   *  \param std::string : its form is : "pluginName:functionName"
   *  \return a std::string containing the function name
   */
  const std::string getPluginFunctionName(const std::string&) const ;

};

#endif //SICONOSSHAREDLIBRARY_H
