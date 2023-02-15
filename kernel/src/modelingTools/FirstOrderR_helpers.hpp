/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#include "SiconosException.hpp"


namespace siconos::modeling::FirstOrderRHelpers {

static inline void JacglambdaSetter(siconos::modeling::FirstOrderR& rel,
                                    std::shared_ptr<siconos::algebra::SimpleMatrix> B,
                                    std::string& pluginName) {
  if (B) {
    rel.setBPtr(B);
  } else if (!pluginName.empty()) {
    rel.setComputeJacglambdaFunction(siconos::plugins::getPluginName(pluginName),
                                     siconos::plugins::getPluginFunctionName(pluginName));
  } else
    THROW_EXCEPTION("FirstOrderRHelpers::JacglambdaSetter no B or pluginJacglambda given");
}

static inline void JachxSetter(siconos::modeling::FirstOrderR& rel,
                               std::shared_ptr<siconos::algebra::SimpleMatrix> C,
                               std::string& pluginName) {
  if (C) {
    rel.setCPtr(C);
  } else if (!pluginName.empty()) {
    rel.setComputeJachxFunction(siconos::plugins::getPluginName(pluginName),
                                siconos::plugins::getPluginFunctionName(pluginName));
  } else {
    THROW_EXCEPTION("FirstOrderRHelpers::JachxSetter no C or pluginJachx given");
  }
}

static inline void JachlambdaSetter(siconos::modeling::FirstOrderR& rel,
                                    std::shared_ptr<siconos::algebra::SimpleMatrix> D,
                                    std::string& pluginName) {
  if (D) {
    rel.setCPtr(D);
  } else if (!pluginName.empty()) {
    rel.setComputeJachlambdaFunction(siconos::plugins::getPluginName(pluginName),
                                     siconos::plugins::getPluginFunctionName(pluginName));
  } else {
    THROW_EXCEPTION("FirstOrderRHelpers::JachlambdaSetter no D or pluginJachlambda given");
  }
}
}  // namespace siconos::modeling::FirstOrderRHelpers

#endif
