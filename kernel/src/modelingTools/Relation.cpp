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
#include "Relation.hpp"
#include "Interaction.hpp"
#include "PluggedObject.hpp"

#include <iostream>

// Default constructor
Relation::Relation(RELATION::TYPES type,
                   RELATION::SUBTYPES subtype):
  _relationType(type), _subType(subtype)
{
  _zeroPlugin();
}

void Relation::_zeroPlugin()
{
  _pluginh.reset(new PluggedObject());
  _pluginJachx.reset(new PluggedObject());
  _pluginJachz.reset(new PluggedObject());
  _pluginJachlambda.reset(new PluggedObject());
  _pluging.reset(new PluggedObject());
  _pluginJacgx.reset(new PluggedObject());
  _pluginJacglambda.reset(new PluggedObject());
  _pluginf.reset(new PluggedObject());
  _plugine.reset(new PluggedObject());
}

void Relation::display() const
{
  std::cout << "=====> Relation of type "
            << _relationType
            << " and subtype "
            << _subType <<std::endl;
}

void Relation::setComputeJachlambdaFunction(const std::string& pluginPath, const std::string& functionName)
{
  _pluginJachlambda->setComputeFunction(pluginPath, functionName);
}

void Relation::setComputeJachxFunction(const std::string& pluginPath, const std::string& functionName)
{
  _pluginJachx->setComputeFunction(pluginPath, functionName);
}

void Relation::setComputeJachzFunction(const std::string& pluginPath, const std::string& functionName)
{
  _pluginJachz->setComputeFunction(pluginPath, functionName);
}
void Relation::setComputegFunction(const std::string& pluginPath, const std::string& functionName)
{
  _pluging->setComputeFunction(pluginPath, functionName);
}
void Relation::setComputeFFunction(const std::string& pluginPath, const std::string& functionName)
{
  _pluginf->setComputeFunction(pluginPath, functionName);
}
void Relation::setComputeEFunction(const std::string& pluginPath, const std::string& functionName)
{
  _plugine->setComputeFunction(pluginPath, functionName);
}

void Relation::setComputeJacgxFunction(const std::string& pluginPath, const std::string& functionName)
{
  _pluginJacgx->setComputeFunction(pluginPath, functionName);
}

void Relation::setComputeJacglambdaFunction(const std::string& pluginPath, const std::string& functionName)
{
  _pluginJacglambda->setComputeFunction(pluginPath, functionName);
}

void Relation::setComputehFunction(const std::string& pluginPath, const std::string& functionName)
{
  _pluginh->setComputeFunction(pluginPath, functionName);
}
