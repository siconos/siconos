/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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

/*! \file SiconosSharedLibrary.h
*/

#ifndef SICONOSSHAREDLIBRARY_H
#define SICONOSSHAREDLIBRARY_H

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

/** Plug-in utilities
 *
 * \author SICONOS Development Team - copyright INRIA
 * \date (creation) 07/21/2006
 *  Matrices can be either block or Simple.
 *  See Derived classes for details.
 */
class SiconosSharedLibrary
{
private:

  /* to manage plug-in object closing (for each call to dlopen one object is add into isPlugged)*/
  std::vector<PluginHandle> isPlugged;

public:

  /**   Constructor
   *    \param unsigned int n, the number of plug-in that will be "connected" through the current object (optional)
   */
  SiconosSharedLibrary(const unsigned int& = 0);

  /** Destructor
   */
  ~SiconosSharedLibrary();

  /**   loads a plugin
   *  \param std::string pluginPath : full plugin path name
   *  \exception SiconosSharedLibraryException if plugin fail to open
   *  \return PluginHandle : plugin handle
   */
  PluginHandle loadPlugin(const std::string& pluginPath) ;

  /**   Gets procedure address
   *  \param PluginHandle plugin : plugin handle
   *  \param std::string procedure : procedure name
   *  \exception SiconosSharedLibraryException if procedure not found
   *  \return pointer on procedure
   */
  void* getProcAddress(PluginHandle plugin, const std::string& procedure);

  /**  Closes plugin
   *  \param PluginHandle plugin : plugin handle
   */
  void closePlugin(PluginHandle plugin);

  /**   Closes all plugin set using the current object
   */
  void closeAllPlugins();

  /** Gets shared library extension
   *  \return library extension ("*.so" for UNIX or "*.dll" for WNT)
   */
  const std::string getSharedLibraryExtension() const ;

  /** set a function pointer to a function in an external library. Don't use it directely
  *  \param functionPtr : pointer to the function (in-out)
  *  \param std::string : the complet path to the external plugin (in)
  *  \param std::string : the name of the function to reference (in)
  */
  void setFunction(void* functionPtr, const std::string& pluginPath, const std::string& functionName);

  /** extract the plugin name from a std::string containing data to call a plugin function
  *  \param std::string : its form is : "pluginName:functionName"
  *  \return a std::string containing the plugin name
  */
  const std::string getPluginName(const std::string&) const ;

  /** extract the function name from a std::string containing data to call a plugin function
  *  \param std::string : its form is : "pluginName:functionName"
  *  \return a std::string containing the function name
  */
  const std::string getPluginFunctionName(const std::string&) const ;

};

#endif //SICONOSSHAREDLIBRARY_H
