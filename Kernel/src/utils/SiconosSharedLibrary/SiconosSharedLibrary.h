#ifndef SICONOSSHAREDLIBRARY_H
#define SICONOSSHAREDLIBRARY_H

#include <string>
#include <iostream>
#include "SiconosSharedLibraryException.h"

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

using namespace std;

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

  /** \fn PluginHandle loadPlugin(const string& pluginPath)
   *  \brief load a plugin
   *  \param string pluginPath : full plugin path name
   *  \exception SiconosSharedLibraryException if plugin fail to open
   *  \return PluginHandle : plugin handle
   */
  PluginHandle loadPlugin(const string& pluginPath);

  /** \fn void* getProcAddress(PluginHandle plugin, const string& procedure)
   *  \brief get procedure address
   *  \param PluginHandle plugin : plugin handle
   *  \param string procedure : procedure name
   *  \exception SiconosSharedLibraryException if procedure not found
   *  \return pointer on procedure
   */
  void* getProcAddress(PluginHandle plugin, const string& procedure);

  /** \fn void closePlugin(PluginHandle plugin)
   *  \brief close plugin
   *  \param PluginHandle plugin : plugin handle
   */
  void closePlugin(PluginHandle plugin);

  /** \fn string getSharedLibraryExtension()
   *  \brief get shared library extension
   *  \return library extension ("*.so" for UNIX or "*.dll" for WNT)
   */
  string getSharedLibraryExtension();


  /** \fn void setFunction(void* functionPtr, PluginHandle pluginHandle, string pluginPath, string functionName)
   *  \brief set a function pointer to a function in an external library. Don't use it directely
   *  \param functionPtr : pointer to the function (in-out)
   *  \param string : the complet path to the external plugin (in)
   *  \param string : the name of the function to reference (in)
   */
  void setFunction(void* functionPtr, string pluginPath, string functionName);

  /** \fn string getPluginName(string )
   *  \brief extract the plugin name from a string containing data to call a plugin function
   *  \param string : its form is : "pluginName:functionName"
   *  \return a string containing the plugin name
   */
  string getPluginName(string);

  /** \fn string getPluginFunctionName(string )
   *  \brief extract the function name from a string containing data to call a plugin function
   *  \param string : its form is : "pluginName:functionName"
   *  \return a string containing the function name
   */
  string getPluginFunctionName(string);
};

#endif //SICONOSSHAREDLIBRARY_H
