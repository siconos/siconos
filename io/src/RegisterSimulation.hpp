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

#ifndef RegisterModel_hpp
#define RegisterModel_hpp

#include "SiconosConfig.h"
#ifdef WITH_SERIALIZATION

#include <fstream>
#include <SiconosFwd.hpp>

void RegisterSimulationOxml(std::ofstream& ofs, SP::Simulation&);
void RegisterSimulationObin(std::ofstream& ofs, SP::Simulation&);
void RegisterSimulationIxml(std::ifstream& ifs, SP::Simulation&);
void RegisterSimulationIbin(std::ifstream& ifs, SP::Simulation&);

#endif

#endif
