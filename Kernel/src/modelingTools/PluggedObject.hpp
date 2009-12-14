#ifndef PluggedObject_H
#define PluggedObject_H

#include "SiconosPointers.hpp"
#include "SiconosSharedLibrary.hpp"
#include <string>

/*! \file PluggedObject.hpp
  \brief utilities for plugin definition to compute matrices or vectors from user-defined functions.
*/

/** Tool to connect a vector or a matrix to a plugin

    This object derived from SimpleMatrix or SimpleVector and handle a pointer to function
    connected to the function defined by pluginName.

    pluginName = "fileName:functionName", fileName being the name of the file, without extension, which
    contains the function functionName.

    The first template parameter defines the type of the function pointer, and the second the base class
    type (ie SimpleMatrix or SimpleVector).

*/
class PluggedObject
{
private:


  std::string pluginName;

public:

  /** plug-in */
  PluginHandle fPtr;

  /** Default Constructor
   */
  PluggedObject(): pluginName("unplugged")
  {
    fPtr = 0;
  }

  /** Constructor with the plugin name */
  PluggedObject(const std::string& name): pluginName(name)
  {
    fPtr = 0;
    setComputeFunction();
  }

  /** bool to checked if a function is properly connected to the current object
   */
  bool isPlugged() const
  {
    return fPtr;
  }

  /** destructor
   */
  ~PluggedObject() {};

  /** Connect a function to fPtr
   \param pluginPath: name of the file where the function is defined (WITH extension)
   \param functionName: name of the function to be connected
  */
  void setComputeFunction(const std::string& pluginPath, const std::string& functionName)
  {
    SSL::setFunction(&fPtr, pluginPath, functionName);
    pluginName = pluginPath.substr(0, pluginPath.length() - 3) + ":" + functionName;
  }
  void setComputeFunction(const std::string& plugin)
  {
    SSL::setFunction(&fPtr, SSL::getPluginName(plugin), SSL::getPluginFunctionName(plugin));
    pluginName = plugin;
  }

  /** Connect pluginName to fPtr => pluginName must have been set before !!
   */
  void setComputeFunction()
  {
    assert(pluginName != "unplugged" && "PluggedObject::setComputeFunction error, try to plug an unamed function.");
    SSL::setFunction(&fPtr, SSL::getPluginName(pluginName), SSL::getPluginFunctionName(pluginName));
  }

  /** Connect input function to fPtr
      \param a T functionPtr
   */
  void setComputeFunction(void* functionPtr)
  {
    fPtr = functionPtr;
    pluginName = "Unknown";
  }

  /** Return the name of the plugin used to compute fPtr
   */
  std::string getPluginName() const
  {
    return pluginName;
  }

  /** Set the name of the plugin function
      \param string of the form pluginFile:functionName, without extension for pluginFile
  */
  void setPluginName(std::string& name)
  {
    pluginName = name;
  }

};
TYPEDEF_SPTR(PluggedObject);
#endif
