/* Siconos-Kernel, Copyright INRIA 2005-2015
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

/*! \file FirstOrderR_helpers.hpp
\brief Set of helpers for FirstOrderR
 */

#ifndef FirstOrderR_helpers_H
#define FirstOrderR_helpers_H

#include "FirstOrderR.hpp"

namespace FirstOrderRHelpers
{

static inline void JacglambdaSetter(FirstOrderR& rel, SP::SiconosMatrix B, std::string& pluginName)
{
  if (B)
  {
    rel.setBPtr(B);
  }
  else if (!pluginName.empty())
  {
    rel.setComputeJacglambdaFunction(SSLH::getPluginName(pluginName), SSLH::getPluginFunctionName(pluginName));
  }
  else
    RuntimeException::selfThrow("FirstOrderRHelpers::JacglambdaSetter no B or pluginJacglambda given");
}

static inline void JachxSetter(FirstOrderR& rel, SP::SiconosMatrix C, std::string& pluginName)
{
  if (C)
  {
    rel.setCPtr(C);
  }
  else if (!pluginName.empty())
  {
    rel.setComputeJachxFunction(SSLH::getPluginName(pluginName), SSLH::getPluginFunctionName(pluginName));
  }
  else
  {
    RuntimeException::selfThrow("FirstOrderRHelpers::JachxSetter no C or pluginJachx given");
  }
}

static inline void JachlambdaSetter(FirstOrderR& rel, SP::SiconosMatrix D, std::string& pluginName)
{
  if (D)
  {
    rel.setCPtr(D);
  }
  else if (!pluginName.empty())
  {
    rel.setComputeJachlambdaFunction(SSLH::getPluginName(pluginName), SSLH::getPluginFunctionName(pluginName));
  }
  else
  {
    RuntimeException::selfThrow("FirstOrderRHelpers::JachlambdaSetter no D or pluginJachlambda given");
  }
}
}

#endif
