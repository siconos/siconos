#include "SiconosSharedLibrary.h"
using namespace std;

//*************************************************************************
//
// Public methods
//
//*************************************************************************

///////////////////////////////////////////////////////////////////////////
//
// loadPlugin
//
PluginHandle SiconosSharedLibrary::loadPlugin(const string& pluginPath) const
{
  PluginHandle HandleRes;
  cout << "PluginPath   :  " << pluginPath << endl;
#ifdef _SYS_WNT
  HandleRes = LoadLibrary(pluginPath.c_str());
#endif
#ifdef _SYS_UNX
  HandleRes = dlopen(pluginPath.c_str(), RTLD_LAZY);
  //  cout << "dlerror() :" <<dlerror()<< endl;
  if (HandleRes == NULL)
  {
    SiconosSharedLibraryException::selfThrow("SiconosSharedLibrary::loadPlugin");
    //throw SiconosSharedLibraryException(dlerror());
  }
#endif
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
