/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

/*! \file FirstOrderR_helpers.hpp
\brief Set of helpers for FirstOrderR
 */

#ifndef FirstOrderR_helpers_H
#define FirstOrderR_helpers_H

#include "FirstOrderR.hpp"

namespace FirstOrderRHelpers
{

static inline void JacglambdaSetter(FirstOrderR& rel, SP::SimpleMatrix B, std::string& pluginName)
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
    THROW_EXCEPTION("FirstOrderRHelpers::JacglambdaSetter no B or pluginJacglambda given");
}

static inline void JachxSetter(FirstOrderR& rel, SP::SimpleMatrix C, std::string& pluginName)
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
    THROW_EXCEPTION("FirstOrderRHelpers::JachxSetter no C or pluginJachx given");
  }
}

static inline void JachlambdaSetter(FirstOrderR& rel, SP::SimpleMatrix D, std::string& pluginName)
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
    THROW_EXCEPTION("FirstOrderRHelpers::JachlambdaSetter no D or pluginJachlambda given");
  }
}
}

#endif
