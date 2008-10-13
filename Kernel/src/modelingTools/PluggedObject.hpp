#ifndef PluggedObject_H
#define PluggedObject_H

#include "SiconosPointers.h"
#include "SiconosSharedLibrary.h"
#include <string>

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

  PluggedObject(unsigned int row): U(row), plugged(false), pluginName("unplugged")
  {}

  PluggedObject(unsigned int row, unsigned int col): U(row, col), plugged(false), pluginName("unplugged")
  {}

  PluggedObject(const U& v): U(v), plugged(false), pluginName("unplugged")
  {}

  PluggedObject operator=(const U& in)
  {
    U::operator=(in);
    plugged = false;
    pluginName = "unplugged";
    return *this;
  }

  void setPlugged(bool in)
  {
    plugged = in;
  }

  bool isPlugged() const
  {
    return plugged;
  }

  /** destructor
   */
  ~PluggedObject() {};

  void setComputeFunction(const std::string& pluginPath, const std::string& functionName)
  {
    SSL::setFunction(&fPtr, pluginPath, functionName);
    pluginName = pluginPath.substr(0, pluginPath.length() - 3) + ":" + functionName;
    plugged = true;
  }

  void setComputeFunction(T functionPtr)
  {
    fPtr = functionPtr;
    plugged = true;
    pluginName = "Unknown";
  }

  std::string getPluginName() const
  {
    return pluginName;
  }

};

#endif
