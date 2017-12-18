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

/*! \file SiconosRestart.hpp
  \brief provides pre-compiled functions for a full Siconos Model
  serialization through libSiconosIO library */

#ifndef SICONOSRESTART_HPP
#define SICONOSRESTART_HPP

#include <SiconosFwd.hpp>
#include <string>

namespace Siconos
{

/** save a Siconos Model with the full simulation state into a
 *  file
 * \param model
 * \param filename with extension : .xml, .bin (binary archive)
 */
void save(SP::Model model, const std::string& filename);

/** load a Siconos Model with the full simulation state from file
 * \param filename
 * \return a SP::Model
 */
SP::Model load(const std::string& filename);

}

#endif
