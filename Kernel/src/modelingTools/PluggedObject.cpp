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


#include "SSLH.hpp"

#include "PluggedObject.hpp"

PluggedObject::PluggedObject(): _pluginName("unplugged")
{
  fPtr = 0;
}

PluggedObject::PluggedObject(const std::string& name): _pluginName(name)
{
  fPtr = 0;
  setComputeFunction();
}

void PluggedObject::setComputeFunction(const std::string& pluginPath, const std::string& functionName)
{
  SSLH::setFunction(&fPtr, pluginPath, functionName);
  _pluginName = pluginPath.substr(0, pluginPath.find_last_of(".")) + ":" + functionName;
}

void PluggedObject::setComputeFunction(const std::string& plugin)
{
  SSLH::setFunction(&fPtr, SSLH::getPluginName(plugin), SSLH::getPluginFunctionName(plugin));
  _pluginName = plugin;
}

void PluggedObject::setComputeFunction(void)
{
  assert(_pluginName != "unplugged" && "PluggedObject::setComputeFunction error, try to plug an unamed function.");
  SSLH::setFunction(&fPtr, SSLH::getPluginName(_pluginName), SSLH::getPluginFunctionName(_pluginName));
}
