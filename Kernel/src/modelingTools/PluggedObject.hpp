/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

#ifndef PluggedObject_H
#define PluggedObject_H

#include <string>
#include "SiconosSerialization.hpp"
#include "SiconosPointers.hpp"

/*! \file PluggedObject.hpp
  \brief utilities for plugin definition to compute matrices or vectors from user-defined functions.
*/

class PluggedObject
{
private:

  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(PluggedObject);

protected:

  /** Plugin name of the form "fileName:functionName" */
  std::string _pluginName;

public:

  /** plug-in */
  void * fPtr;

  /** Default Constructor
   */
  PluggedObject();

  /** Constructor with the plugin name
   * \param name a string of the form "fileName:functionName", without an extension for pluginFile
   */
  PluggedObject(const std::string& name);

  /** Copy constructor
   * \param PO a PluggedObject we are going to copy
  */
  PluggedObject(const PluggedObject & PO):  _pluginName(PO.getPluginName()), fPtr(PO.fPtr) {};

  /** bool to checked if a function is properly connected to the current object
   * \return a boolean, true if fPtr is set
   */
  inline bool isPlugged() const
  {
    return (fPtr != 0);
  };

  /** destructor
   */
  ~PluggedObject() {};

  /** Connect a function to fPtr
   \param pluginPath name of the file where the function is defined (WITH extension)
   \param functionName name of the function to be connected
  */
  void setComputeFunction(const std::string& pluginPath, const std::string& functionName);

  /* Connect a function to fPtr
   * \param plugin a string of the form "fileName:functionName,  without an extension for pluginFile"
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
  inline std::string getPluginName(void) const
  {
    return _pluginName;
  };

  /** Set the name of the plugin function
      \param name a string of the form "pluginFile:functionName", without extension for pluginFile
  */
  inline void setPluginName(const std::string& name)
  {
    _pluginName = name;
  };
};

TYPEDEF_SPTR(PluggedObject);

#endif
