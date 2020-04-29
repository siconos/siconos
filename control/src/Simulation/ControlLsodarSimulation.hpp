/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

/*!\file ControlLsodarSimulation.hpp
 * \brief Simulation with Control using Lsodar
 */

#ifndef ControlLsodarSimulation_hpp
#define ControlLsodarSimulation_hpp

#include "ControlSimulation.hpp"

class ControlLsodarSimulation : public ControlSimulation {

public:
  ControlLsodarSimulation(double t0, double T, double h);

  virtual void run();

};
#endif
