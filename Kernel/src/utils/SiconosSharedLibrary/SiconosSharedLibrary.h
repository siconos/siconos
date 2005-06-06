#ifndef SICONOSSHAREDLIBRARY_H
#define SICONOSSHAREDLIBRARY_H

#include "SiconosSharedLibraryException.h"
#include <string>
#include <iostream>
#include <dlfcn.h>

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
public:

  /** \fn SiconosSharedLibrary(){}
   *  \brief constructor
   */
  SiconosSharedLibrary() {}

  /** \fn ~SiconosSharedLibrary(){}
   *  \brief destructor
   */
  ~SiconosSharedLibrary() {}

  /** \fn PluginHandle loadPlugin(const std::string& pluginPath)
   *  \brief load a plugin
   *  \param std::string pluginPath : full plugin path name
   *  \exception SiconosSharedLibraryException if plugin fail to open
   *  \return PluginHandle : plugin handle
   */
  PluginHandle loadPlugin(const std::string& pluginPath) const ;

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
