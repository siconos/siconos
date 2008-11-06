#ifndef PluggedObject_H
#define PluggedObject_H

#include "SiconosPointers.hpp"
#include "SiconosSharedLibrary.h"
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
template <class T, class U> class PluggedObject : public U
{
private:

  bool plugged;

  std::string pluginName;

public:

  /** plug-in */
  T fPtr;

  /** Default Constructor
   */
  PluggedObject(): U(), plugged(false), pluginName("unplugged")
  {}

  /** Constructor with vector size */
  PluggedObject(unsigned int row): U(row), plugged(false), pluginName("unplugged")
  {}

  /** Constructor with matrix dim */
  PluggedObject(unsigned int row, unsigned int col): U(row, col), plugged(false), pluginName("unplugged")
  {}

  /** Copy-constructor from a U (base-class type) */
  PluggedObject(const U& v): U(v), plugged(false), pluginName("unplugged")
  {}

  /** Constructor with the plugin name */
  PluggedObject(const std::string& name): U(), plugged(false), pluginName(name)
  {
    setComputeFunction();
  }

  /** Copy from a U (base-class type) */
  PluggedObject operator=(const U& in)
  {
    U::operator=(in);
    plugged = false;
    pluginName = "unplugged";
    return *this;
  }

  /** bool to checked if a function is properly connected to the current object
   */
  bool isPlugged() const
  {
    return plugged;
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
    plugged = true;
  }

  /** Connect pluginName to fPtr => pluginName must have been set before !!
   */
  void setComputeFunction()
  {
    assert(pluginName != "unplugged" && "PluggedObject::setComputeFunction error, try to plug an unamed function.");
    SSL::setFunction(&fPtr, SSL::getPluginName(pluginName), SSL::getPluginFunctionName(pluginName));
    plugged = true;
  }

  /** Connect input function to fPtr
      \param a T functionPtr
   */
  void setComputeFunction(T functionPtr)
  {
    fPtr = functionPtr;
    plugged = true;
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

#endif
