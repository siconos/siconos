#ifndef PluggedVector_H
#define PluggedVector_H

#include "SiconosPointers.h"
#include "SiconosSharedLibrary.h"
#include "SimpleVector.h"
#include <string>


template <typename T> class PluggedVector : public SimpleVector
{
private:

  bool plugged;

  std::string pluginName;

public:

  /** plug-in */
  T fPtr;

  /** Default Constructor
   */
  PluggedVector(): SimpleVector(), plugged(false), pluginName("unplugged")
  {}

  PluggedVector(unsigned int size): SimpleVector(size), plugged(false), pluginName("unplugged")
  {}

  PluggedVector(const SimpleVector& v): SimpleVector(v), plugged(false), pluginName("unplugged")
  {}

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
  ~PluggedVector() {};

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
