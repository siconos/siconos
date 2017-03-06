/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#ifndef PluggedObject_H
#define PluggedObject_H

#include <string>
#include "SiconosSerialization.hpp"
#include "SiconosPointers.hpp"

/*! \file PluggedObject.hpp
  \brief utilities for plugin definition to compute matrices or vectors from user-defined functions.
*/


/** Class to deal with plugged functions 

A plugin is a C-function defined in some external file.

This object handles a function pointer to this C-function.
*/
class PluggedObject
{
private:

  /* serialization hooks */
  ACCEPT_SERIALIZATION(PluggedObject);

protected:

  /** Plugin name, should be of the form "fileName:functionName" */
  std::string _pluginName;

public:

  /** plug-in */
  void * fPtr;

  /** Default Constructor
   */
  PluggedObject();

  /** Constructor with the plugin name
   * \param name a std::string of the form "fileName:functionName", without an extension for pluginFile
   */
  PluggedObject(const std::string& name);

  /** Copy constructor
   * \param PO a PluggedObject we are going to copy
  */
  PluggedObject(const PluggedObject & PO);

  /** bool to checked if a function is properly connected to the current object
   * \return a boolean, true if fPtr is set
   */
  inline bool isPlugged() const
  {
    return (fPtr != NULL);
  };

  /** destructor
   */
  virtual ~PluggedObject();

  /** Connect a function to fPtr
   \param pluginPath name of the file where the function is defined (WITH extension)
   \param functionName name of the function to be connected
  */
  void setComputeFunction(const std::string& pluginPath, const std::string& functionName);

  /* Connect a function to fPtr
   * \param plugin a std::string of the form "fileName:functionName,  without an extension for pluginFile"
   */
  void setComputeFunction(const std::string& plugin);

  /** Connect _pluginName to fPtr
   * \warning _pluginName must have been set before !
   */
  void setComputeFunction(void);

  /** Connect input function to fPtr
      \param functionPtr a pointer to a C function
   */
  inline void setComputeFunction(void* functionPtr)
  {
    fPtr = functionPtr;
    _pluginName = "Unknown";
  };

  /** Return the name of the plugin used to compute fPtr
   * \return _pluginName (a std::string)
   */
  inline std::string pluginName(void) const
  {
    return _pluginName;
  };

};

#endif
