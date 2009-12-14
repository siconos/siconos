/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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

/*! \file DynamicalSystem.h
  \brief Abstract class - General interface for all Dynamical Systems.
*/

#ifndef PLUGIN_HPP
#define PLUGIN_HPP

#include "SiconosSharedLibrary.hpp"

class Plugin
{

public:

  static bool setFunction(void* f,  const std::string& pluginPath, const std::string& functionName, std::string& name)
  {
    SSL::setFunction(f, pluginPath, functionName);
    name = pluginPath.substr(0, pluginPath.length() - 3) + ":" + functionName;
    return true;
  }

  static bool setFunction(void* f,  const std::string& pluginPath, const std::string& functionName)
  {
    SSL::setFunction(f, pluginPath, functionName);
    return true;
  }
  static bool setFunction(void* f,  const std::string& Name)
  {
    SSL::setFunction(f, SSL::getPluginName(Name), SSL::getPluginFunctionName(Name));
    return true;
  }


};

#endif


